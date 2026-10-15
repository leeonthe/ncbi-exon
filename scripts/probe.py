"""Validation probe: hit a small curated set of genes and print results.

Goal: let a human eyeball whether the API + parser are producing the right data,
by cross-checking against https://www.ncbi.nlm.nih.gov/datasets/gene/<gene_id>/

Run from project root:
    NCBI_API_KEY=... python scripts/probe.py
    # or without a key (slower):
    python scripts/probe.py

Add more symbols at the bottom by editing PROBE_SYMBOLS or pass them as args:
    python scripts/probe.py TP53 BRCA1 EGFR
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Allow running as a script without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ncbi_exon import NCBIClient, NCBIError, parse_product_report  # noqa: E402
from ncbi_exon.env import load_dotenv  # noqa: E402

# Load NCBI_API_KEY (and anything else) from .env at project root, if present.
load_dotenv()

# Curated probe set:
#   - TP53, BRCA1, BRCA2, EGFR, MYC: well-characterized, many transcripts (is_exon should be False)
#   - GAPDH, ACTB: housekeeping; useful sanity check
#   - MIR21: microRNA — should have NR_ only, no NM_ (is_exon = False, nm_count = 0)
#   - HOTAIR: lncRNA — also NR_ only
#   - FAKE123: deliberately bogus, expect NOT_FOUND
PROBE_SYMBOLS = [
    "TP53",
    "BRCA1",
    "BRCA2",
    "EGFR",
    "MYC",
    "GAPDH",
    "ACTB",
    "MIR21",
    "HOTAIR",
    "FAKE123",
]


def main(argv: list[str]) -> int:
    symbols = argv[1:] or PROBE_SYMBOLS

    cache_dir = Path(__file__).resolve().parents[1] / "cache"
    client = NCBIClient(cache_dir=cache_dir)
    key_status = "with API key" if client.api_key else "no API key (5 rps)"
    print(f"# NCBI probe — {len(symbols)} symbols, {key_status}, cache={cache_dir}\n")

    header = f"{'SYMBOL':<10}{'STATUS':<12}{'GENE_ID':<10}{'NM':>4}{'NR':>4}{'TX':>5}{'is_exon':>9}  {'MANE_SELECT':<18}  DESCRIPTION"
    print(header)
    print("-" * len(header))

    rows = []
    for symbol in symbols:
        try:
            res = client.fetch_product_report(symbol)
        except NCBIError as e:
            row = parse_product_report(symbol, payload={}, error_detail=str(e))
        else:
            row = parse_product_report(
                symbol,
                payload=res.payload,
                http_status=res.status_code,
                error_detail="" if res.status_code in (200, 404) else f"HTTP {res.status_code}",
            )
            cache_tag = " (cache)" if res.from_cache else f" ({res.latency_ms}ms)"
            row.error_detail = (row.error_detail + cache_tag).strip()

        rows.append(row)
        print(
            f"{row.input_symbol:<10}"
            f"{row.status:<12}"
            f"{row.gene_id:<10}"
            f"{row.nm_count:>4}"
            f"{row.nr_count:>4}"
            f"{row.transcript_count:>5}"
            f"{str(row.is_exon):>9}  "
            f"{(row.mane_select_accession or '-'):<18}  "
            f"{row.description[:60]}"
        )

    print()
    print("# Sanity checks (printed for human review):")
    by_symbol = {r.input_symbol: r for r in rows}

    def check(label: str, ok: bool, detail: str = ""):
        mark = "PASS" if ok else "FAIL"
        print(f"  [{mark}] {label}" + (f" — {detail}" if detail else ""))

    if "TP53" in by_symbol:
        r = by_symbol["TP53"]
        check("TP53 resolves to gene_id 7157", r.gene_id == "7157", f"got {r.gene_id!r}")
        check("TP53 has multiple NM transcripts", r.nm_count > 1, f"nm_count={r.nm_count}")
        check("TP53 is_exon is False (nm_count != 1)", r.is_exon is False)
        check("TP53 MANE Select is NM_000546.x", r.mane_select_accession.startswith("NM_000546."),
              f"got {r.mane_select_accession!r}")

    if "FAKE123" in by_symbol:
        r = by_symbol["FAKE123"]
        check("FAKE123 is NOT_FOUND", r.status == "NOT_FOUND")

    if "MIR21" in by_symbol:
        r = by_symbol["MIR21"]
        check("MIR21 has nm_count == 0 (miRNA, no protein-coding mRNA)", r.nm_count == 0,
              f"nm_count={r.nm_count}")

    # Any gene in our probe set that happens to have exactly 1 NM should flip is_exon=True.
    one_nm = [r for r in rows if r.status == "OK" and r.nm_count == 1]
    if one_nm:
        print(f"  [INFO] {len(one_nm)} probe gene(s) had exactly 1 NM → is_exon=True:")
        for r in one_nm:
            print(f"         {r.input_symbol}: NM={r.nm_accessions}")
    else:
        print("  [INFO] No probe gene had exactly nm_count == 1 — pick one from the website and add it to verify the True branch.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
