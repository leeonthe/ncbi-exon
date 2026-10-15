"""Parse an NCBI Datasets v2 gene product_report payload into a flat row.

See README.md §3, §5 and `Approach_non_unique_transcripts.md` for the contract.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any

# GRCh38 primary chromosome accessions all start with NC_0000 — used to filter
# out alternate assemblies (T2T-CHM13 etc.) before exon overlap analysis.
GRCH38_PREFIX = "NC_0000"


@dataclass
class GeneRow:
    """One output row per input gene.

    Field order = CSV column order. See `Approach_non_unique_transcripts.md` §2
    for the meaning and provenance of every field.
    """

    input_symbol: str
    status: str                        # OK | NOT_FOUND | MULTI_MATCH | ERROR
    gene_id: str = ""
    resolved_symbol: str = ""
    description: str = ""
    tax_id: str = ""
    gene_type: str = ""
    is_exon: bool = False              # True iff nm_count == 1 (strict)
    nm_count: int = 0
    nr_count: int = 0
    xm_count: int = 0                  # predicted protein-coding (XM_)
    xr_count: int = 0                  # predicted non-coding (XR_)
    transcript_count: int = 0
    protein_count: int = 0
    protein_coding_count: int = 0
    non_coding_count: int = 0
    nm_accessions: str = ""            # ';'-joined
    nr_accessions: str = ""
    xm_accessions: str = ""
    xr_accessions: str = ""
    mane_select_accession: str = ""
    mane_plus_clinical_accessions: str = ""
    # Overlap analysis — Definition #1 (partial overlap, "3-rule"; share >=1 bp).
    # See Approach_non_unique_transcripts.md §3.
    has_non_unique_transcript: bool = False
    non_unique_transcript_count: int = 0
    non_unique_transcript_accessions: str = ""
    transcripts_compared_for_overlap: int = 0
    transcripts_excluded_from_overlap: int = 0
    # Overlap analysis — Definition #2 (containment, "2-rule", directional):
    # a candidate NM is non-unique vs T' iff EVERY exon of the candidate sits
    # fully INSIDE some exon of T' (partial overlap does NOT count). Same NM
    # candidate set and GRCh38 filtering as definition #1, so `compared` /
    # `excluded` above are shared. See Approach_non_unique_transcripts.md §3.1.
    has_non_unique_transcript_containment: bool = False
    non_unique_transcript_count_containment: int = 0
    non_unique_transcript_accessions_containment: str = ""
    error_detail: str = ""
    fetched_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds"))

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def merge_paged_payloads(payloads: list[dict[str, Any]]) -> dict[str, Any]:
    """Merge a list of paginated product_report responses into a single payload.

    Concatenates `reports[0].product.transcripts` across pages and recomputes
    `transcript_count`. Other product fields are taken from the first page.
    No-op for a single payload.
    """
    if not payloads:
        return {}
    if len(payloads) == 1:
        return payloads[0]
    base = json.loads(json.dumps(payloads[0]))  # deep copy
    base_reports = base.get("reports") or []
    if not base_reports:
        return base
    base_product = base_reports[0].setdefault("product", {})
    merged_transcripts = list(base_product.get("transcripts") or [])
    for extra in payloads[1:]:
        extra_reports = extra.get("reports") or []
        if not extra_reports:
            continue
        extra_product = extra_reports[0].get("product") or {}
        merged_transcripts.extend(extra_product.get("transcripts") or [])
    base_product["transcripts"] = merged_transcripts
    base_product["transcript_count"] = len(merged_transcripts)
    return base


def parse_product_report(
    input_symbol: str,
    payload: dict[str, Any],
    *,
    http_status: int = 200,
    error_detail: str = "",
) -> GeneRow:
    """Turn a raw product_report response into a GeneRow.

    Handles all four status outcomes (OK / NOT_FOUND / MULTI_MATCH / ERROR).
    """
    if error_detail or (http_status >= 500):
        return GeneRow(
            input_symbol=input_symbol,
            status="ERROR",
            error_detail=error_detail or f"HTTP {http_status}",
        )

    total = int(payload.get("total_count") or 0)
    reports = payload.get("reports") or []
    if http_status == 404 or total == 0 or not reports:
        return GeneRow(input_symbol=input_symbol, status="NOT_FOUND")

    # MULTI_MATCH disambiguation — see README §3.
    if total > 1 or len(reports) > 1:
        primary = [r for r in reports
                   if str((r.get("product") or {}).get("symbol", "")).upper() == input_symbol.upper()]
        all_ids = [str((r.get("product") or {}).get("gene_id", "")) for r in reports]
        if len(primary) == 1:
            product = primary[0].get("product") or {}
            row = _build_ok_row(input_symbol, product)
            other_ids = [g for g in all_ids if g and g != row.gene_id]
            if other_ids:
                row.error_detail = f"also matched as synonym of: {','.join(other_ids)}"
            return row
        return GeneRow(
            input_symbol=input_symbol,
            status="MULTI_MATCH",
            error_detail=f"matched gene_ids: {','.join(g for g in all_ids if g)}",
        )

    product = reports[0].get("product") or {}
    return _build_ok_row(input_symbol, product)


def _build_ok_row(input_symbol: str, product: dict[str, Any]) -> GeneRow:
    """Construct an OK GeneRow from a single `product` dict."""
    transcripts = product.get("transcripts") or []

    nm_accs: list[str] = []
    nr_accs: list[str] = []
    xm_accs: list[str] = []
    xr_accs: list[str] = []
    mane_select: str = ""
    mane_plus: list[str] = []

    for t in transcripts:
        acc = (t.get("accession_version") or "").strip()
        if not acc:
            continue
        if acc.startswith("NM_"):
            nm_accs.append(acc)
        elif acc.startswith("NR_"):
            nr_accs.append(acc)
        elif acc.startswith("XM_"):
            xm_accs.append(acc)
        elif acc.startswith("XR_"):
            xr_accs.append(acc)
        sel = t.get("select_category")
        if sel == "MANE_SELECT" and not mane_select:
            mane_select = acc
        elif sel == "MANE_PLUS_CLINICAL":
            mane_plus.append(acc)

    type_counts = {tc.get("type"): int(tc.get("count") or 0) for tc in (product.get("transcript_type_counts") or [])}

    nm_count = len(nm_accs)
    overlap = _analyze_overlap(transcripts)

    return GeneRow(
        input_symbol=input_symbol,
        status="OK",
        gene_id=str(product.get("gene_id", "")),
        resolved_symbol=str(product.get("symbol", "")),
        description=str(product.get("description", "")),
        tax_id=str(product.get("tax_id", "")),
        gene_type=str(product.get("type", "")),
        is_exon=(nm_count == 1),
        nm_count=nm_count,
        nr_count=len(nr_accs),
        xm_count=len(xm_accs),
        xr_count=len(xr_accs),
        transcript_count=int(product.get("transcript_count") or len(transcripts)),
        protein_count=int(product.get("protein_count") or 0),
        protein_coding_count=type_counts.get("PROTEIN_CODING", 0),
        non_coding_count=type_counts.get("NON_CODING", 0),
        nm_accessions=";".join(nm_accs),
        nr_accessions=";".join(nr_accs),
        xm_accessions=";".join(xm_accs),
        xr_accessions=";".join(xr_accs),
        mane_select_accession=mane_select,
        mane_plus_clinical_accessions=";".join(mane_plus),
        has_non_unique_transcript=overlap["has_non_unique"],
        non_unique_transcript_count=overlap["non_unique_count"],
        non_unique_transcript_accessions=";".join(overlap["non_unique_accs"]),
        transcripts_compared_for_overlap=overlap["compared"],
        transcripts_excluded_from_overlap=overlap["excluded"],
        has_non_unique_transcript_containment=overlap["has_non_unique_containment"],
        non_unique_transcript_count_containment=overlap["non_unique_count_containment"],
        non_unique_transcript_accessions_containment=";".join(overlap["non_unique_accs_containment"]),
    )


# ---------------------------------------------------------------- overlap

def _grch38_exons(transcript: dict[str, Any]) -> list[tuple[int, int]] | None:
    """Return the (begin, end) tuples for the GRCh38 primary location, or None
    if this transcript has no usable GRCh38 entry / no exons.
    """
    for loc in transcript.get("genomic_locations") or []:
        gav = loc.get("genomic_accession_version") or ""
        if not gav.startswith(GRCH38_PREFIX):
            continue
        exons = loc.get("exons") or []
        out: list[tuple[int, int]] = []
        for e in exons:
            try:
                b = int(e.get("begin"))
                t = int(e.get("end"))
            except (TypeError, ValueError):
                continue
            if b > t:
                b, t = t, b
            out.append((b, t))
        if out:
            return out
        return None
    return None


def _intervals_overlap(a: tuple[int, int], b: tuple[int, int]) -> bool:
    """Definition #1 predicate — 1-based inclusive interval overlap: share at
    least one base. (Partial overlap counts.)"""
    return a[0] <= b[1] and b[0] <= a[1]


def _interval_contained(inner: tuple[int, int], outer: tuple[int, int]) -> bool:
    """Definition #2 predicate — 1-based inclusive: True iff `inner` lies fully
    within `outer` (inner.begin >= outer.begin AND inner.end <= outer.end).
    Partial overlap and touching-but-not-sharing both return False."""
    return outer[0] <= inner[0] and inner[1] <= outer[1]


def _is_non_unique_against(t_exons: list[tuple[int, int]],
                           other_exons: list[tuple[int, int]]) -> bool:
    """Definition #1 (partial overlap, "3-rule"). Per
    Approach_non_unique_transcripts.md §3: t is non-unique vs `other` iff every
    exon of t overlaps at least one exon of other.
    """
    if not t_exons or not other_exons:
        return False
    for ea in t_exons:
        if not any(_intervals_overlap(ea, eb) for eb in other_exons):
            return False
    return True


def _is_non_unique_against_containment(t_exons: list[tuple[int, int]],
                                       other_exons: list[tuple[int, int]]) -> bool:
    """Definition #2 (containment, "2-rule", directional). Per
    Approach_non_unique_transcripts.md §3.1: t is non-unique vs `other` iff every
    exon of t is fully contained inside at least one exon of other.
    """
    if not t_exons or not other_exons:
        return False
    for ea in t_exons:
        if not any(_interval_contained(ea, eb) for eb in other_exons):
            return False
    return True


def _run_pairwise(usable: list[tuple[str, list[tuple[int, int]]]],
                  predicate) -> list[str]:
    """Return the accessions of transcripts that are non-unique under `predicate`.

    T is non-unique iff there exists at least one *other* T' such that
    `predicate(exons_of_T, exons_of_T')` is True. Pairwise, one-vs-one — never
    the union of all others (Approach_non_unique_transcripts.md §3, rule 2).
    """
    non_unique: list[str] = []
    for i, (acc_i, exons_i) in enumerate(usable):
        for j, (_acc_j, exons_j) in enumerate(usable):
            if i == j:
                continue
            if predicate(exons_i, exons_j):
                non_unique.append(acc_i)
                break
    return non_unique


def _analyze_overlap(transcripts: list[dict[str, Any]]) -> dict[str, Any]:
    """Run the pairwise non-unique check across the NM transcripts of a gene.

    Per Approach_non_unique_transcripts.md §3 (NM-only rule, updated):
    only curated protein-coding mRNAs (`NM_` prefix) participate in the overlap
    comparison. Predicted transcripts (XM/XR) are treated as "not yet
    confirmed" and curated non-coding RNAs (NR) are out of scope for the
    coding-isoform redundancy question. Non-NM transcripts are skipped at the
    candidate stage and do NOT count toward `transcripts_excluded_from_overlap`
    (that field only counts NM transcripts that had no usable GRCh38 data).
    """
    usable: list[tuple[str, list[tuple[int, int]]]] = []
    excluded = 0
    for t in transcripts:
        acc = (t.get("accession_version") or "").strip()
        if not acc.startswith("NM_"):
            continue  # NM-only rule — skip NR/XM/XR silently
        exons = _grch38_exons(t)
        if exons is None:
            excluded += 1
            continue
        usable.append((acc, exons))

    # Both definitions run over the same NM candidate set; only the per-pair
    # predicate differs (partial overlap vs full containment).
    partial = _run_pairwise(usable, _is_non_unique_against)
    contained = _run_pairwise(usable, _is_non_unique_against_containment)

    return {
        "compared": len(usable),
        "excluded": excluded,
        # Definition #1 — partial overlap ("3-rule").
        "non_unique_count": len(partial),
        "non_unique_accs": partial,
        "has_non_unique": bool(partial),
        # Definition #2 — containment ("2-rule", directional).
        "non_unique_count_containment": len(contained),
        "non_unique_accs_containment": contained,
        "has_non_unique_containment": bool(contained),
    }
