# Visualization Plan — PhD-level Figures for the Project

Companion to `README.md`, `Approach_non_unique_transcripts.md`, and `output/all_protein_coding_v2.csv` (n = 20,593 unique human protein-coding gene symbols, of which 20,583 yielded an OK API response).

> **Note (updated 2026-06-02):** the overlap rule was tightened to **NM-only** — only curated protein-coding mRNAs participate in the pairwise comparison. Predicted (XM / XR) and curated non-coding (NR) transcripts are still reported in their own count columns but do not enter the overlap test. This changes the headline non-uniqueness rate from 67.6 % → 50.7 % and creates a hard implication: `is_exon = True ⇒ has_non_unique = False`. All numbers below reflect the NM-only rule.

This document proposes **six thesis-grade figures** that take the per-gene CSV from "a table" to "a story". For each figure I describe **the question it answers**, **what the visual encoding shows**, **the insight surfaced by the actual data**, and **what statistical / biological argument it supports**. All six are already generated under `output/figures/` (PNG @ 300 dpi + vector PDF, both committed to the same script for reproducibility).

---

## 0. The thesis-level question this set of figures answers

> *How structurally redundant is the annotated human protein-coding transcriptome, and how does redundancy distribute across (a) curated vs. computationally predicted RefSeq products and (b) the diversity gradient from single-isoform to highly-spliced genes?*

The CSV gives us three orthogonal axes per gene:

1. **Transcript inventory** — `nm_count`, `nr_count`, `xm_count`, `xr_count` (counts of the four RefSeq prefixes).
2. **Canonicality** — `is_exon` (exactly one curated mRNA → single-canonical class).
3. **Structural redundancy** — `has_non_unique_transcript`, `non_unique_transcript_count` (per the pairwise overlap rule in `Approach_non_unique_transcripts.md` §3).

Each figure isolates a slice of the joint distribution of those axes and tests a specific claim about the transcriptome.

---

## 1. Recommended figure suite

### Figure 1 — Transcriptome composition overview

| | |
|---|---|
| **Question** | What is the relative abundance of curated vs. predicted RefSeq transcripts, and how concentrated is the NM count? |
| **Panel A** | Stacked horizontal bar — total transcripts pooled across all OK genes, partitioned by prefix (NM / NR / XM / XR). |
| **Panel B** | Histogram (log-Y) of `nm_count` per gene, with a dashed line at `nm_count = 1` and annotated count of `is_exon = True` genes. |
| **Insight from actual data** | NM: 70,925 · NR: 7,878 · XM: **127,798** · XR: 8,758. **Predicted protein-coding transcripts (XM) outnumber curated mRNAs (NM) by 1.8×**, a counterintuitive but important framing for any downstream argument about "the human transcriptome". 7,805 genes (37.9 % of OK) have exactly one NM — the modal class. |
| **PhD-level argument** | The bulk of RefSeq's apparent transcript diversity comes from computational prediction, not manual curation. Any claim about "the human transcriptome" needs to specify which slice is being counted. |

### Figure 2 — Non-uniqueness landscape

| | |
|---|---|
| **Question** | How frequently does the non-unique-transcript condition fire, and how does it scale with the number of transcripts compared? |
| **Panel A** | Bar chart: fraction of OK genes with `has_non_unique_transcript = True`, binned by `transcripts_compared_for_overlap` (1, 2, 3, 4–5, 6–10, 11–20, 21–50, 51–100, 100+). |
| **Panel B** | Scatter (log-log) of `transcript_count` vs `non_unique_transcript_count`, with the y = x diagonal overlaid. Points on the diagonal = every transcript of the gene is non-unique. |
| **Panel C** | 2 × 2 heatmap of `is_exon` × `has_non_unique_transcript`, cells annotated with count and % of OK. |
| **Insight from actual data (NM-only)** | The fraction is 0 % at 1 NM transcript (single-NM genes can never be flagged by definition), rises sharply at 2 NMs, and saturates near 100 % for genes with 6+ NMs. Panel C is now a **triangular table** because `is_exon = True ⇒ has_non_unique = False` is automatic under NM-only: the (True, True) cell is empty by construction. The clean class (`is_exon = True`, `has_non_unique = False`) is **37.9 %** of OK; the redundant class (`is_exon = False`, `has_non_unique = True`) is **50.7 %**; the multi-NM-but-distinguishable class (`is_exon = False`, `has_non_unique = False`) is **11.3 %**. |
| **PhD-level argument** | Among multi-NM genes (those that *could* be flagged, i.e., `nm_count ≥ 2`, n = 11,816), **88.4 %** turn out to have at least one NM that is a structural subset of another. The redundancy is concentrated in the multi-isoform regime, while single-NM genes form a clean baseline that the rule cannot disturb. (If the denominator is broadened to *all* non-mono-isoform genes including the 962 zero-NM cases, the fraction drops to 81.7 % — but those genes are structurally exempt from the rule, so the 88.4 % figure is the biologically meaningful one.) |

### Figure 3 — Curated vs predicted contribution to redundancy

| | |
|---|---|
| **Question** | The rule itself is now NM-only, so XM/XR no longer *cause* the flag. But does the **presence** of predicted transcripts still correlate with non-uniqueness (suggesting some shared underlying property of highly-spliced loci)? |
| **Panel A** | Violin plot of `log(1 + xm_count)` partitioned by `has_non_unique_transcript`, with a Mann–Whitney U test reporting in the panel. |
| **Panel B** | Stacked bar showing the mean transcript-type composition (% NM / NR / XM / XR) within five NM-count buckets (0, 1 = is_exon, 2–3, 4–6, 7+). |
| **Insight from actual data (NM-only)** | Mann–Whitney **p ≈ 0** — non-unique genes carry meaningfully more XM transcripts (mean 8.98 vs 3.35) even though XMs are not part of the rule. This is now a **correlation**, not a confounder: highly-spliced loci tend to attract both more NM curation *and* more XM prediction, and both surface together. The composition stack shows the predicted share growing with NM count. |
| **PhD-level argument** | Under the original mixed-prefix rule the XM signal was partly causal (XMs literally added comparison candidates that could flag the gene). Under the NM-only rule the signal is purely associative — informative about the underlying gene biology, not a methodological artifact. |

### Figure 4 — Per-gene distribution of transcript counts

| | |
|---|---|
| **Question** | What is the per-gene distribution shape — Gaussian, exponential, or heavy-tailed? Does the non-uniqueness count follow the same scaling? |
| **Panel A** | Empirical CDF of three quantities on the same log-X axis: total transcripts, NM-only, non-unique count. |
| **Panel B** | Rank–frequency plot (log–log) of total transcripts per gene. A straight line on this axis is the signature of a power-law (Zipfian) distribution. |
| **Insight from actual data** | The rank–frequency plot is approximately linear on log–log between rank ~10 and ~10,000 (slope ≈ −0.6 to −0.8) before dropping at the tail. This is consistent with a Zipfian regime: a small number of "transcript-rich" loci dominate the total count. |
| **PhD-level argument** | Standard descriptive statistics (mean, sd) are misleading for this distribution; report the median, IQR, and the top-decile concentration. Power-law-like behavior also motivates log-scale axes throughout the rest of the figure suite. |

### Figure 5 — Exon-level worked examples (qualitative)

| | |
|---|---|
| **Question** | Why does the rule fire on a real gene? What does "non-unique" look like as a picture? |
| **Encoding** | Three stacked panels — ACTB, MYC, GAPDH — each showing every GRCh38 transcript as one horizontal track. Exons are filled blocks; introns are gray baselines. Non-unique transcripts are colored **red**; unique transcripts are **blue**. Y-axis labels are accession numbers. |
| **Insight from actual data** | ACTB (1 transcript): all blue (trivially clean). MYC (2 transcripts): both red — the only difference between them is the start of exon 1. GAPDH (6 transcripts): the topmost transcript (`NM_001256799.3`) is blue because its 5′-most exon starts ~750 bp downstream of the others' 5′ exons, giving it a unique exon vs every other transcript; the remaining 5 are red. |
| **PhD-level argument** | These three panels make the rule legible to a reviewer in five seconds. Best used as the very first figure shown when introducing the non-unique-transcript concept — it grounds the methodology in real coordinates before the aggregate statistics. |

### Figure 6 — Classification flow

| | |
|---|---|
| **Question** | What fraction of all input symbols ultimately lands in each terminal bucket? |
| **Encoding** | Four-stage bar diagram, each stage proportional to the previous: Stage 1 = all 20,593 input symbols; Stage 2 = API status (OK / NOT_FOUND / MULTI_MATCH / ERROR); Stage 3 = `is_exon` partition within OK; Stage 4 = `has_non_unique_transcript` within each `is_exon` bucket. |
| **Insight from actual data (NM-only)** | 99.95 % of inputs resolved (OK), 0.05 % dropped (all 9 NOT_FOUND are mitochondrial genes; 1 MULTI_MATCH is COX2 as legacy alias). Among OK: 37.9 % `is_exon = True` → 100 % of those land in `has_non_unique = False` (forced by the rule). Among genes with `nm_count ≥ 2` (the eligible set, n = 11,816): **88.4 % are `has_non_unique = True`**. The remaining 962 OK genes have 0 NMs and are structurally exempt. |
| **PhD-level argument** | A single figure that summarizes the entire pipeline outcome. The Sankey now has one fewer branch than under the old rule because the (`is_exon=True`, `has_non_unique=True`) branch is structurally empty — making the figure cleaner and the "clean single-canonical" class (= the entire `is_exon = True` column, 37.9 % of OK) immediately visible. |

---

## 2. Proposed extensions (Phase 2 figures for a full thesis)

These are not yet generated but follow naturally from the same CSV.

### Figure 7 — Genomic distribution (chromosome × redundancy)
- **Encoding**: heatmap of chromosomes (rows) × NM-count buckets (cols), coloring = fraction with `has_non_unique = True`.
- **Why it matters**: tests whether gene-dense chromosomes (1, 19) show different redundancy profiles from sparse ones (Y, MT). Requires chromosome parsing from the cached payload (`product.transcripts[].genomic_locations[].genomic_accession_version` → chromosome).

### Figure 8 — Sensitivity to curated-only analysis
- **Encoding**: side-by-side mosaic plots — left = current rule (NM + NR + XM + XR), right = re-run with XM/XR excluded from the overlap analysis.
- **Why it matters**: addresses the prediction-artifact concern from Figure 3. If `has_non_unique = True` drops substantially in the curated-only re-run, the headline figure should be qualified.

### Figure 9 — Gene length / exon-count correlates
- **Encoding**: scatter of total exon span (max end − min begin on GRCh38) vs `non_unique_transcript_count`, colored by `is_exon`.
- **Why it matters**: tests whether large loci are biologically more redundant or whether the redundancy reflects annotation density.

### Figure 10 — MANE Select coverage cross-tab
- **Encoding**: cross-tab of "has MANE_SELECT" vs `is_exon` and `has_non_unique`.
- **Why it matters**: MANE Select is the curated consensus canonical transcript. Genes with a MANE Select pick should have a defensible "primary" transcript regardless of redundancy.

### Figure 11 — Junction-level disambiguation (if pursued in thesis scope)
- **Encoding**: for each gene currently flagged `has_non_unique = True`, recompute the flag using exon-*boundary* matching instead of exon-*range* overlap.
- **Why it matters**: a stricter rule (exact boundary matching) better reflects biological isoform identity than the inclusive-overlap rule used in v2. This is a methods comparison rather than a new biological claim.

---

## 3. Standards used in the generated figures

| Aspect | Choice | Rationale |
|---|---|---|
| Color palette | `viridis` (perceptually uniform, color-blind friendly) | Avoids the rainbow / jet pitfalls; reproduces in B/W |
| Axes | Log-X / log-Y where distributions span >1 decade | The CSV has 0–368 NM per gene, 1–368 transcripts |
| Statistical annotation | Mann–Whitney U, n labels on each bar | Non-parametric — distributions are not Gaussian |
| Resolution | 300 dpi PNG + vector PDF | Thesis-print quality; PDFs scale losslessly |
| Reproducibility | `scripts/make_figures.py`, single run from CSV | Re-runs in <30 s, deterministic |
| Spine cleanup | Top + right removed | Reduces chart junk per Tufte |

---

## 4. Suggested thesis narrative order

The six figures slot naturally into a results section in this order:

1. **Figure 5** — introduce the rule with three concrete examples (visceral)
2. **Figure 6** — show where everything ended up (orientation)
3. **Figure 1** — describe what the transcriptome is made of (descriptive)
4. **Figure 4** — show how it's distributed per gene (statistical shape)
5. **Figure 2** — answer the central redundancy question (main result)
6. **Figure 3** — partition by curation status (mechanistic / methodological caveat)

This sequence walks the reader from concrete → orienting → descriptive → analytical → mechanistic, which is the conventional thesis Results structure.

---

## 5. What the figures collectively show

Reading the six together, the project's central finding (NM-only rule) can be stated as:

> Among the 20,583 human protein-coding genes with usable RefSeq product reports, **50.7 % carry at least one curated mRNA (NM) whose every exon overlaps an exon of another NM of the same gene** (the "non-unique transcript" condition). Because the rule operates only on NMs, single-NM genes are structurally exempt: the **37.9 %** of OK genes with `is_exon = True` automatically have `has_non_unique = False`, partitioning the transcriptome cleanly into three classes — single-NM (37.9 %), zero-NM (4.7 %, also exempt by construction), and multi-NM (57.4 %). Within the multi-NM class (n = 11,816, the only one the rule can actually flag), **88.4 %** are non-unique. Non-uniqueness is independently associated with the presence of predicted (XM) transcripts (Mann–Whitney p ≈ 0), suggesting that highly-spliced loci attract both more curated and more predicted annotation in parallel — a biological signal, not a methodological artifact of the rule itself.

That single paragraph, plus Figures 5 / 6 / 2 / 3, is enough material for one results subsection of a methods-focused chapter.

---

## 6. Files produced

```
output/
├── all_protein_coding_v2.csv     (3.7 MB, 20,593 rows, v2 schema)
├── all_protein_coding_v2.log     (run log, ts-prefixed)
├── test50_v2.csv                 (50-gene smoke test, v2)
└── figures/
    ├── fig1_composition.{png,pdf}
    ├── fig2_non_uniqueness.{png,pdf}
    ├── fig3_curated_vs_predicted.{png,pdf}
    ├── fig4_distributions.{png,pdf}
    ├── fig5_exon_tracks.{png,pdf}
    └── fig6_classification_flow.{png,pdf}

scripts/
└── make_figures.py    (single-file reproducible generator)
```

Re-run with: `.venv/bin/python scripts/make_figures.py`
