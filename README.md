# NCBI Exon — Gene Transcript Annotation Pipeline

> **Read this entire document before writing or modifying any code in this project.**
> All design decisions, API choices, data contracts, and edge-case rules are defined here.
> If something in the code disagrees with this README, the README wins until it is explicitly updated.

---

## 1. Project Goal

Given a list of ~20,000 human gene symbols (CSV or XLSX input), for each gene:

1. Query NCBI for the gene's RefSeq products (transcripts).
2. Count how many transcript accessions begin with `NM_` and `NR_`.
3. Mark `is_exon = True` **iff the gene has exactly one `NM_` transcript**.
4. Store as much supporting data as practical (gene ID, total transcript count, MANE Select accession, etc.) so downstream analysis does not require re-querying.

The pipeline must be re-runnable, resumable, rate-limit-aware, and produce one output row per input gene (no silent drops).

---

## 1.5 Implementation Status (last updated 2026-06-02, end of day)

**Pipeline is complete (v2 schema), full 20,593-gene CSV is produced, and a thesis-grade visualization suite is generated.** Outstanding work is scientific extension, not engineering.

### ✅ Completed — engineering

| Item | Where | Notes |
|---|---|---|
| Sync NCBI Datasets v2 client | `ncbi_exon/client.py` | stdlib `urllib`, retries + `Retry-After`, file cache, rate-limit gate. Used by the probe. |
| Async NCBI client (httpx + semaphore) | `ncbi_exon/async_client.py` | async leaky-bucket rate limiter, concurrent semaphore, `next_page_token` pagination, shared file cache |
| Product-report parser + `GeneRow` schema (v1) | `ncbi_exon/parser.py` | NM/NR counting, `is_exon = (nm_count == 1)`, MANE Select/Plus, MULTI_MATCH disambiguation (NAT1 → gene_id 9) |
| **v2 schema: XM/XR counts + dual-definition overlap analyzer (NM-only rule, updated 2026-06-04)** | `ncbi_exon/parser.py` | Adds 12 columns: `xm_count`, `xr_count`, `xm_accessions`, `xr_accessions`, `has_non_unique_transcript`, `non_unique_transcript_count`, `non_unique_transcript_accessions`, `transcripts_compared_for_overlap`, `transcripts_excluded_from_overlap` (Definition #1, partial overlap), plus `has_non_unique_transcript_containment`, `non_unique_transcript_count_containment`, `non_unique_transcript_accessions_containment` (Definition #2, containment — every candidate exon fully inside a partner exon, directional). **Both definitions computed every run; only NM transcripts participate** (NR/XM/XR reported in their own count columns but not compared). GRCh38-only. See `Approach_non_unique_transcripts.md` §3 (Definition #1) and §3.1 (Definition #2). |
| Page-merge helper | `ncbi_exon/parser.merge_paged_payloads` | Concatenates transcripts across pages before parsing |
| CSV/XLSX input reader | `ncbi_exon/input_reader.py` | Auto-detects `Symbol`/`gene` column; filters `type_of_gene == 'protein-coding'`; dedupes; logs counts |
| Streaming CSV output writer | `ncbi_exon/output_writer.py` | Flushes per row → survives Ctrl-C |
| CLI entry point | `python -m ncbi_exon` | `--input`, `--output`, `--cache-dir`, `--limit`, `--refresh`, `--rps`, `--concurrency`, `--taxon`, `--sheet`, `--no-protein-coding-filter`, `--log` |
| In-order async output | `ncbi_exon/__main__.py` | Sliding-window buffer → rows written in input order regardless of completion order (idempotent reruns) |
| Progress + `run.log` + graceful Ctrl-C | `ncbi_exon/__main__.py` | `[i/N]` counter every 25, ts-prefixed log, SIGINT handler flushes partial output |
| `.env` loader + secrets protection | `ncbi_exon/env.py`, `.gitignore`, `.env.example` | `.env` gitignored; key verified to exist only in `.env` |
| Validation probe (10 curated genes) | `scripts/probe.py` | TP53 → 26 (matches screenshot); ACTB → `is_exon=True`; MIR21/HOTAIR/FAKE123 all behave |
| 50-gene v2 smoke test (NM-only) | `output/test50_v2.csv` | 50/50 OK, 11 `is_exon=True`, sanity checks (ACTB / MYC / GAPDH / BRCA2) vs Approach doc all pass |
| **Full 20,593-gene v2 run (NM-only rule)** | `output/all_protein_coding_v2.csv` | 20,583 OK, 9 NOT_FOUND (mitochondrial), 1 MULTI_MATCH (COX2). 7,805 `is_exon=True` (37.9 %), 10,443 `has_non_unique=True` (50.7 %). Cold 37 min (rate-limited at 10 rps); warm-cache regen 9 s. |
| Schema decision: drop `synonyms` / `ensembl_gene_id` | §5.1 | One API call per gene |

### ✅ Completed — analysis & visualization

| Item | Where | Notes |
|---|---|---|
| Approach doc for overlap rule | `Approach_non_unique_transcripts.md` | Logic in plain English; worked examples for ACTB / MYC / GAPDH using real exon coordinates from API |
| 6-figure thesis-grade visualization suite | `output/figures/*.{png,pdf}` | 300 dpi PNG + vector PDF. Reproducible via `scripts/make_figures.py`. |
| Visualization plan | `Visualization_plan.md` | Per-figure question / encoding / actual-data insight / argument; standards (viridis, log axes, Mann–Whitney); suggested thesis narrative order; 5 Phase-2 extensions identified |

### §8 Verification checklist — all green

- [x] Probe runs end-to-end on TP53, BRCA1, BRCA2, EGFR, MYC, GAPDH, ACTB, MIR21, HOTAIR, FAKE123
- [x] TP53 → `nm_count=25`, `is_exon=False`, MANE = NM_000546.6
- [x] ACTB → `nm_count=1`, `is_exon=True`, MANE = NM_001101.5, `has_non_unique=False` (single transcript)
- [x] MYC → 2 NMs both flagged `non_unique` (mutually overlap)
- [x] GAPDH → 4 of 5 NMs flagged `non_unique` (NM_001256799.3 has a unique 5′ start; NR_152150.2 excluded from comparison)
- [x] BRCA2 → 5 of 6 NMs flagged `non_unique` (NM_001432077.1 unique — different 5′ exon)
- [x] FAKE123 → `status=NOT_FOUND` with empty counts
- [x] Rerun on populated cache finishes in <15 s for full 20 k
- [x] Idempotent: two consecutive runs produce byte-identical CSVs (excluding `fetched_at`)
- [x] API key loaded from `NCBI_API_KEY` env / `.env`; verified to exist only in `.env`
- [x] v2 schema produces 32 columns (29 base + 3 Definition-#2 containment columns); new columns populated correctly for OK rows and empty for NOT_FOUND / MULTI_MATCH

### ✅ Final paper figure set (added 2026-08-27, per `PCG Transcripts Computational Analysis Research - Jun&Yoon.pdf`)

The reviewed figure suite was culled and redesigned per the Jun&Yoon revision PDF. Only **three figures** survive into the paper, all under **Definition #2 (containment)**, generated by `scripts/make_final_figures.py` into `output/figures_final/`:

| Final figure | Source | Requirements from the PDF |
|---|---|---|
| **Figure 1** (3 panels) | old Fig 1 + old Fig 2A | **1A**: RefSeq composition as a *horizontal bar chart* (one bar per prefix, in the style of the example mockup). **1B**: NM-per-gene histogram with *integer x ticks* (no 2.5); the 0-NM bar labeled as XM-predicted-only genes (n=962); the 1-NM bar highlighted and labeled "Mono-isoform" in plain words (no `is_exon`-style code speak). **1C**: old Fig 2A with y-axis "Fraction with structurally non-unique NM transcript(s)" and no overlapping text/labels. |
| **Figure 2** | old Fig 5, GAPDH panel only | GAPDH exon tracks only; legend wording "Structurally unique transcript" / "Structurally non-unique transcript". |
| **Figure 3** | old Fig 6, redesigned | Flowchart-style classification (like the PDF's example): all input → successfully processed → **0 NM vs 1 NM (mono-isoform) vs multi-NM** → only multi-NM further split into contains / doesn't contain structurally non-unique transcript, with quantitative (count + %) annotations. |

Cross-cutting rules from the PDF: figure **title and subtitle must be visible**; the term is "**structurally non-unique**" (not just "non-unique"); old Figs 3 & 4 are dropped; mono-isoform prose must count the 962 zero-NM (XM-only) genes as their own class, not lump them into "two or more NM".

```bash
# Final paper figures (Definition #2), from the v2 CSV
.venv/bin/python scripts/make_final_figures.py \
  --input output/all_protein_coding_v2.csv \
  --output-dir output/figures_final
```

The exploratory 6-figure suites (`output/figures/`, `output/figures_def2/`, `scripts/make_figures.py`) remain for reference.

### ⏳ Remaining — scientific extension (not engineering)

The pipeline is feature-complete. Open items are research directions, listed in `Visualization_plan.md` §2:

- [ ] **Figure 7** — chromosome × redundancy heatmap (requires parsing chromosome from cached payload)
- [ ] **Figure 8** — sensitivity analysis: compare the original mixed-prefix rule (NM + NR + XM + XR) against the current NM-only rule, quantifying how much of the previous 67.6 % non-unique rate collapsed to 50.7 % once predicted and non-coding transcripts were excluded from the comparison set
- [ ] **Figure 9** — gene length / exon-span vs `non_unique_transcript_count` correlation
- [ ] **Figure 10** — MANE Select coverage cross-tab
- [ ] **Figure 11** — junction-level disambiguation (stricter rule: exact exon-boundary matching) — methods comparison

None require additional API calls; all are derivable from the existing cache + CSV.

### Run commands

```bash
# Validation probe (10 curated genes, sub-second from cache)
.venv/bin/python scripts/probe.py

# 50-gene v2 smoke test (cached after first run)
.venv/bin/python -m ncbi_exon \
  --input "Protein Coding genes (human).xlsx" \
  --output output/test50_v2.csv \
  --cache-dir cache \
  --limit 50

# Full 20 k v2 run (12 s from cache; ~37 min cold)
.venv/bin/python -m ncbi_exon \
  --input "Protein Coding genes (human).xlsx" \
  --output output/all_protein_coding_v2.csv \
  --cache-dir cache

# Regenerate all 6 figures (≤30 s)
.venv/bin/python scripts/make_figures.py
```

### Headline result (from `output/all_protein_coding_v2.csv`, NM-only rule)

| Metric | Value |
|---|---|
| Genes processed | 20,593 unique symbols |
| OK | 20,583 (99.95 %) |
| `is_exon = True` (exactly 1 NM) | 7,805 (37.9 % of OK) |
| `has_non_unique_transcript = True` | 10,443 (50.7 % of OK) |
| **`is_exon = True` ⇒ `has_non_unique = False` automatically** (single-NM gene has no other NM to compare against) | the two flags are no longer independent |
| `is_exon = False` AND `has_non_unique = False` (multi-NM, all distinguishable on exons) | 2,335 (11.3 % of OK) |
| `is_exon = False` AND `has_non_unique = True` (multi-NM, redundant) | 10,443 (50.7 % of OK) |
| Genes with 0 NM transcripts (compared = 0; flag forced False) | 962 (4.7 % of OK) |
| **Multi-NM genes that are flagged non-unique** (denominator = 11,816 genes with `nm_count ≥ 2`) | **10,443 (88.4 %)** — the biologically meaningful "redundancy among comparable isoforms" rate |
| Curated mRNA transcripts (NM) total | 70,925 |
| Predicted mRNA transcripts (XM) total | 127,798 (1.8× the curated count, reported but NOT in overlap test) |
| Mann–Whitney U test: XM count by non-uniqueness | p ≈ 0 (XM load still associates with non-uniqueness even though XMs aren't in the rule itself) |

---

## 2. API Decision — Use NCBI Datasets v2 REST API

After comparing the available options, the chosen API is the **NCBI Datasets v2 REST API**.

### 2.1 Options considered

| API | Verdict | Why |
|---|---|---|
| **NCBI Datasets v2 REST** (`api.ncbi.nlm.nih.gov/datasets/v2/`) | ✅ **CHOSEN** | Same backend that powers the `ncbi.nlm.nih.gov/datasets/gene/<id>/#transcripts-and-proteins` page shown in the screenshots. Returns structured JSON with the exact `transcripts[].accession_version` field we need. Single endpoint takes gene symbol + taxon and returns the product report — no two-step lookup needed. |
| E-utilities (Entrez) — `esearch` + `efetch` (`eutils.ncbi.nlm.nih.gov`) | ❌ | XML/flat-file based, requires two calls, transcript list is harder to parse, schema is less stable. |
| `datasets` CLI tool | ⚠️ Backup only | Wraps the same REST API. Useful for manual spot-checks but a subprocess-per-gene approach for 20k genes is slower and less observable than direct HTTP. |
| BioPython `Bio.Entrez` | ❌ | Same E-utilities under the hood. |

### 2.2 Chosen endpoint

```
GET https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/{symbol}/taxon/{taxon}/product_report
```

- **`{symbol}`** — gene symbol (e.g. `TP53`). Case-sensitive.
- **`{taxon}`** — use `human` (or `9606` / `Homo sapiens`). For this project, **always `human`**.
- **Headers:**
  - `Accept: application/json`
  - `api-key: <NCBI_API_KEY>` — read from env (see §6).
- **Returns:** JSON with `reports[].product` containing `gene_id`, `transcripts[]`, `transcript_count`, `protein_count`, `transcript_type_counts[]`, and gene metadata.
- **Pagination:** the response has `next_page_token`; for a single gene symbol it almost always fits in one page, but the code MUST still handle pagination defensively.

### 2.3 Verified response shape (TP53 sample)

```json
{
  "reports": [
    {
      "product": {
        "gene_id": "7157",
        "symbol": "TP53",
        "description": "tumor protein p53",
        "tax_id": "9606",
        "taxname": "Homo sapiens",
        "type": "PROTEIN_CODING",
        "transcript_count": 26,
        "protein_count": 25,
        "transcript_type_counts": [
          {"type": "PROTEIN_CODING", "count": 25},
          {"type": "NON_CODING", "count": 1}
        ],
        "transcripts": [
          {
            "accession_version": "NM_000546.6",
            "name": "transcript variant 1",
            "length": 2512,
            "type": "PROTEIN_CODING",
            "select_category": "MANE_SELECT",
            "protein": { "accession_version": "NP_000537.3", "length": 393, ... },
            "ensembl_transcript": "ENST00000269305.9",
            "genomic_locations": [ { "exons": [...] } ]
          }
          // ... 25 more
        ]
      },
      "query": ["TP53"]
    }
  ],
  "total_count": 1
}
```

> ⚠️ `transcript_type_counts` groups by **biological type** (`PROTEIN_CODING` / `NON_CODING`), NOT by accession prefix. To count `NM` / `NR` we MUST iterate `transcripts[]` and inspect `accession_version`.

### 2.4 Rate limits & API key

| Mode | Limit | 20k genes ETA |
|---|---|---|
| No API key | 5 requests/sec | ~67 min |
| With API key (`api-key` header) | 10 requests/sec | ~33 min |

API key strategy for this project: **the user will provide a key via the `NCBI_API_KEY` environment variable**. The code MUST:
- Read `NCBI_API_KEY` from env on startup.
- If present, send `api-key: <key>` header and pace at 10 rps.
- If absent, fall back to 5 rps (log a warning).
- Never log the key.

Get a key from My NCBI → Account Settings → API Key Management.

---

## 3. Counting & Classification Rules

For each gene, after fetching the product report:

```python
nm_count = sum(1 for t in transcripts if t["accession_version"].startswith("NM_"))
nr_count = sum(1 for t in transcripts if t["accession_version"].startswith("NR_"))
is_exon  = (nm_count == 1)   # ← STRICT EQUALITY, per user spec
```

**Rules confirmed with the user (2026-06-01):**
- `is_exon = True` **only when `nm_count == 1` exactly**. A gene with 0 NMs, or with 2+ NMs (e.g. TP53 has many), is `is_exon = False`.
- Track **NM and NR** counts only. We are NOT counting XM / XR in this project. (`XM_` / `XR_` accessions are computationally predicted and are out of scope.)
- Every input gene gets exactly one output row, with a `status` column (see §5.2).

---

## 4. Input Specification

- **Format:** CSV or XLSX file containing gene symbols, path supplied by the user at runtime.
- **Expected column:** a column named `gene` or `symbol` (case-insensitive). If neither exists, fall back to the first column and log the assumption.
- **Volume:** up to ~20,000 rows.
- **Sanitization:** strip whitespace, uppercase the symbol (NCBI gene symbols are uppercase by convention), drop empty rows. Deduplicate (preserve first occurrence) and log the dedupe count.

---

## 5. Output Specification

### 5.1 Primary output — one row per input gene

CSV file (UTF-8). Mirror input order. Columns:

| Column | Type | Description |
|---|---|---|
| `input_symbol` | str | The symbol as it appeared in the input (after strip/upper). |
| `status` | str | `OK` / `NOT_FOUND` / `MULTI_MATCH` / `ERROR` (see §5.2). |
| `gene_id` | int | NCBI Gene ID (e.g. `7157`). Empty on non-OK. |
| `resolved_symbol` | str | Symbol returned by NCBI (may differ from input if the input was a synonym). |
| `description` | str | Gene description (e.g. "tumor protein p53"). |
| `tax_id` | int | Should be `9606`. |
| `gene_type` | str | e.g. `PROTEIN_CODING`. |
| `is_exon` | bool | `True` iff `nm_count == 1`. |
| `nm_count` | int | Count of transcripts whose `accession_version` starts with `NM_`. |
| `nr_count` | int | Count of transcripts whose `accession_version` starts with `NR_`. |
| `transcript_count` | int | Total transcripts (from `product.transcript_count`, sanity-checked vs `len(transcripts)`). |
| `protein_count` | int | From `product.protein_count`. |
| `protein_coding_count` | int | From `transcript_type_counts` where type is `PROTEIN_CODING`. |
| `non_coding_count` | int | From `transcript_type_counts` where type is `NON_CODING`. |
| `xm_count` | int | **(v2)** Count of transcripts whose `accession_version` starts with `XM_` (predicted protein-coding). |
| `xr_count` | int | **(v2)** Count of transcripts whose `accession_version` starts with `XR_` (predicted non-coding). |
| `nm_accessions` | str | `;`-joined NM accessions (versioned, e.g. `NM_000546.6;NM_001126112.3`). |
| `nr_accessions` | str | `;`-joined NR accessions. |
| `xm_accessions` | str | **(v2)** `;`-joined XM accessions. |
| `xr_accessions` | str | **(v2)** `;`-joined XR accessions. |
| `mane_select_accession` | str | The `accession_version` of the transcript where `select_category == "MANE_SELECT"`, if any. Empty otherwise. |
| `mane_plus_clinical_accessions` | str | `;`-joined accessions where `select_category == "MANE_PLUS_CLINICAL"`. |
| `has_non_unique_transcript` | bool | **(v2, NM-only, Definition #1 — partial overlap)** True iff at least one **NM transcript** of this gene is "non-unique" per the rule in `Approach_non_unique_transcripts.md` §3 (every exon of that NM **overlaps ≥1 bp** of an exon of some other NM of the same gene, on GRCh38). NR / XM / XR are reported in their own count columns but do not participate. |
| `non_unique_transcript_count` | int | **(v2)** How many of this gene's **NM** transcripts are non-unique under Definition #1. |
| `non_unique_transcript_accessions` | str | **(v2)** `;`-joined NM accessions of those transcripts (Definition #1). |
| `transcripts_compared_for_overlap` | int | **(v2)** How many **NM** transcripts had usable GRCh38 exon data and entered the pairwise comparison. Shared by both definitions. |
| `transcripts_excluded_from_overlap` | int | **(v2)** How many **NM** transcripts were dropped because they had no GRCh38 location or no exons. Does NOT include NR / XM / XR (those are skipped at the candidate stage). Shared by both definitions. |
| `has_non_unique_transcript_containment` | bool | **(v2, NM-only, Definition #2 — containment)** True iff at least one **NM transcript** is non-unique under the stricter **containment** rule in `Approach_non_unique_transcripts.md` §3.1: every exon of that NM lies **fully inside** an exon of some other NM (partial overlap does NOT count; directional). Always a subset of Definition #1. |
| `non_unique_transcript_count_containment` | int | **(v2)** How many of this gene's **NM** transcripts are non-unique under Definition #2 (containment). |
| `non_unique_transcript_accessions_containment` | str | **(v2)** `;`-joined NM accessions of those transcripts (Definition #2). |
| `error_detail` | str | Empty on OK. On `ERROR`, contains HTTP status / exception class / short message. On `MULTI_MATCH` disambiguated to a primary, contains `also matched as synonym of: <other gene_ids>`. Never contains the API key. |
| `fetched_at` | ISO-8601 timestamp (UTC) | When this row was fetched. Used for cache freshness. |

### 5.2 Status semantics

| `status` | When |
|---|---|
| `OK` | `total_count == 1`, transcripts parsed. |
| `NOT_FOUND` | `total_count == 0` or HTTP 404. |
| `MULTI_MATCH` | `total_count > 1` (rare with symbol+taxon; treat each match by gene_id; for v1 we record all gene_ids in `error_detail` and leave counts empty). |
| `ERROR` | Network failure, 5xx, parse error, timeout after retries. |

### 5.3 Secondary artifacts

- **Raw JSON cache** (one file per gene, e.g. `cache/<SYMBOL>.json`) so the pipeline is resumable and re-runnable without re-hitting the API. The runner MUST consult the cache first and only call the API on cache miss (or if `--refresh` is passed).
- **`run.log`** — one line per gene: `ts | symbol | status | latency_ms`. Plus any retries.

---

## 6. Runtime & Operational Requirements

The implementation MUST include all of the following. Do not skip any.

1. **Rate limiting** — token-bucket or async semaphore. Default 10 rps with API key, 5 rps without. Configurable via `--rps` flag.
2. **Retries** — on `429`, `500`, `502`, `503`, `504`, and transport errors. Exponential backoff (1s, 2s, 4s, 8s), max 4 retries. Respect `Retry-After` header if present.
3. **Timeouts** — 30s connect + read per request.
4. **Resumability** — file-based cache (see §5.3). Re-running the script with the same input + cache directory must skip already-fetched genes.
5. **Progress reporting** — print a counter every N genes (e.g. `tqdm` or a manual `[i/total] symbol — status`).
6. **Graceful interruption** — Ctrl-C should flush the in-memory output buffer to the CSV before exiting, so partial runs are not lost.
7. **No secrets in logs** — never print the API key.
8. **Idempotent** — re-running on a fully cached input must produce a byte-identical output CSV (modulo `fetched_at`).

### 6.1 Suggested tech stack (not mandatory, but recommended)

- Python 3.11+
- `httpx` (async) or `requests` (sync). Async is preferred at this volume.
- `pandas` for input/output (CSV + XLSX via `openpyxl`).
- `tenacity` for retries (or hand-rolled).
- `python-dotenv` for `.env` loading of `NCBI_API_KEY`.

### 6.2 Suggested CLI

```
python -m ncbi_exon \
  --input genes.csv \
  --output results.csv \
  --cache-dir ./cache \
  [--refresh] [--rps 10] [--taxon human]
```

---

## 7. Open Questions / Future Work

- Should we also export a per-transcript-level CSV (one row per transcript)? Out of scope for v1; raw JSON cache covers this.
- XM / XR support: explicitly OUT of scope per user decision on 2026-06-01.
- Symbol-as-synonym handling: NCBI auto-resolves common synonyms, but not all. If `NOT_FOUND` rate is high, consider an Entrez `esearch` fallback. Defer until we see real data.

---

## 8. Verification Checklist (before declaring v1 done)

- [ ] Runs end-to-end on a 10-gene sample including TP53, BRCA1, BRCA2, EGFR, MYC, GAPDH, ACTB, and one deliberately-bad symbol (e.g. `FAKE123`).
- [ ] TP53 row shows `nm_count >= 1`, `is_exon = False` (TP53 has multiple NMs).
- [ ] A gene known to have exactly 1 NM (pick one from the sample run and verify on the NCBI website) shows `is_exon = True`.
- [ ] `FAKE123` row exists with `status = NOT_FOUND` and empty counts.
- [ ] Re-running with cache populated finishes in <5 seconds and produces identical output.
- [ ] API key is loaded from env and never appears in logs.

---

## 9. References

### Project documents (read in this order)

| Doc | What it covers |
|---|---|
| `README.md` (this file) | Project goal, API decision, schema, run instructions, status |
| `Definition_of_non_unique_transcript.md` | User-supplied definition of the non-unique-transcript concept |
| `Approach_non_unique_transcripts.md` | Logic spec for the overlap rule, with worked examples using real API data |
| `Visualization_plan.md` | Per-figure rationale, statistical standards, suggested thesis narrative order, Phase-2 figure proposals |
| `How_to_generate_figures.md` | Copy-paste-ready code snippets (one Jupyter cell per figure) with overlap-free layouts; includes Plotly Sankey alternative for the classification flow |

### External references

- NCBI Datasets v2 REST API docs: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/
- OpenAPI spec (source of truth for endpoint shapes): https://github.com/ncbi/datasets/blob/master/datasets.openapi.yaml
- Gene product report schema: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/gene-product/
- API keys: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/
- Example UI page (TP53): https://www.ncbi.nlm.nih.gov/datasets/gene/7157/#transcripts-and-proteins
