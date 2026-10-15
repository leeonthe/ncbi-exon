# How to Generate the Thesis Figures

Self-contained, copy-paste-ready snippets to regenerate every figure in `output/figures/`.
The snippets use `constrained_layout=True` and explicit positioning to fix the title/label overlap that affected the earlier outputs.

Companion to `Visualization_plan.md` (which explains *why* each figure exists) and `output/all_protein_coding_v2.csv` (the data source).

---

## 0. Best workflow: Jupyter notebook in the project venv

The recommended way to work with these figures for a thesis is a **Jupyter notebook** living next to the CSV. Reasons:

- **Inline preview** — you see each figure as you tweak it. No round trip through the file system.
- **Cell isolation** — re-run only the figure you're editing, not the whole pipeline.
- **Export options** — `File → Save and Export → PDF/HTML` ships a thesis-ready artifact in one click.
- **Reproducibility** — the notebook is the experiment record; you can commit it to git alongside the CSV.

### One-time setup

From the project root:

```bash
# 1. Install Jupyter into the project's existing venv
.venv/bin/pip install jupyterlab ipykernel

# 2. Register the venv as a kernel so notebooks can pick it
.venv/bin/python -m ipykernel install --user --name ncbi-exon --display-name "ncbi-exon (.venv)"

# 3. Launch JupyterLab from the project root
.venv/bin/jupyter lab
```

In the notebook UI:
- `File → New → Notebook` and pick the **"ncbi-exon (.venv)"** kernel.
- Save the notebook as `thesis_figures.ipynb` in the project root so relative paths work.

If you'd rather skip Jupyter, every code block below also works as a standalone `.py` script (just chain the **Setup cell** with the figure cell you want).

---

## 1. Setup cell — paste this first

Every figure cell below assumes you've run this cell first in the same notebook (or scripted these imports at the top of your `.py` file).

```python
# === Setup ===
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Thesis-quality plot defaults
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.right": False,
    "axes.spines.top": False,
})
sns.set_palette("viridis")

# Project paths — adjust if you're not running from the project root
PROJECT = Path.cwd()
CSV_PATH = PROJECT / "output" / "all_protein_coding_v2.csv"
FIG_DIR = PROJECT / "output" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Load the v2 CSV and keep the OK subset around — most figures use it
df = pd.read_csv(CSV_PATH)
ok = df[df["status"] == "OK"].copy()
print(f"loaded {len(df):,} rows ({len(ok):,} OK)")

# Helper — save both a 300dpi PNG and a vector PDF
def save_fig(fig, name):
    fig.savefig(FIG_DIR / f"{name}.png")
    fig.savefig(FIG_DIR / f"{name}.pdf")
    print(f"  saved: {name}.{{png,pdf}}")
```

---

## 2. Figure 1 — Transcript-type composition + NM distribution

**What this figure shows:** A pooled view of all RefSeq transcript types (NM / NR / XM / XR) across the human protein-coding transcriptome, alongside a per-gene histogram of NM count. The headline number — **predicted XM transcripts outnumber curated NM transcripts 1.8×** — is hidden in the API but visible at a glance here.

**What insight it gives you:** This figure forces you to be honest about which "transcriptome" you're talking about. If you describe "human transcript diversity" using a count that mixes NM and XM, you're conflating two epistemically very different things (manual curation vs computational prediction). The histogram with the `nm_count = 1` callout grounds the `is_exon` definition: the modal class is single-NM, but only 37.9 % of OK genes.

**Code (paste as one Jupyter cell):**

```python
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
fig.suptitle("Figure 1 — Transcript-type composition", fontsize=12, fontweight="bold")

# Panel A — stacked horizontal bar
totals = {
    "NM (curated mRNA)":    ok["nm_count"].sum(),
    "NR (curated ncRNA)":   ok["nr_count"].sum(),
    "XM (predicted mRNA)":  ok["xm_count"].sum(),
    "XR (predicted ncRNA)": ok["xr_count"].sum(),
}
colors = sns.color_palette("viridis", 4)
left = 0
for (label, val), c in zip(totals.items(), colors):
    axA.barh([0], val, left=left, color=c, label=f"{label}: {val:,}")
    left += val
axA.set_yticks([])
axA.set_xlabel("Total transcripts (all OK genes pooled)")
axA.set_title("A. RefSeq transcript composition")
axA.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False, fontsize=8)
for s in ("left", "right", "top"):
    axA.spines[s].set_visible(False)

# Panel B — log-Y histogram of NM count per gene
n_exon = int((ok["nm_count"] == 1).sum())
bins = np.arange(0, 22) - 0.5
axB.hist(ok["nm_count"], bins=bins, color=colors[0], edgecolor="white", linewidth=0.6)
axB.axvline(1, color="crimson", ls="--", lw=1, label="nm_count = 1 (is_exon)")
# Annotation positioned BELOW the suptitle to avoid overlap
axB.annotate(
    f"is_exon = True\nn = {n_exon:,}  ({100 * n_exon / len(ok):.1f}%)",
    xy=(1, n_exon),
    xytext=(6, n_exon * 0.35),   # ← lower-right, well clear of the title
    arrowprops=dict(arrowstyle="->", color="crimson", lw=0.8),
    fontsize=9, color="crimson",
)
axB.set_yscale("log")
axB.set_xlim(-0.5, 20.5)
axB.set_xlabel("NM transcripts per gene")
axB.set_ylabel("Number of genes (log)")
axB.set_title(f"B. NM distribution (n = {len(ok):,})")
axB.legend(loc="upper right", frameon=False)

save_fig(fig, "fig1_composition")
plt.show()
```

---

## 3. Figure 2 — Non-uniqueness landscape

**What this figure shows:** Three views of how often the non-unique-transcript rule fires. Panel A bins genes by transcript count and shows the fraction flagged. Panel B is a log-log scatter of total transcripts vs the number of non-unique transcripts in that gene; the dashed `y = x` line marks the "every transcript is non-unique" regime. Panel C is a 2 × 2 mosaic of `is_exon` × `has_non_unique_transcript`.

**What insight it gives you (NM-only rule):** Non-uniqueness is concentrated in the multi-NM regime — by construction, single-NM genes (`is_mono-isoform = True`, 37.9 % of OK) cannot be flagged because there is no other NM to compare against. Among the genuine *multi*-NM subset (`nm_count ≥ 2`, n = 11,816), **88.4 %** of genes are flagged non-unique. The 2 × 2 mosaic (Panel C) is now **triangular**: the (True, True) cell is structurally empty under the NM-only rule, and the "clean single-canonical" class (`is_mono-isoform = True`, `has_non_unique = False`) is the entire `is_mono-isoform = True` row — **37.9 %** of OK, not 7.2 %.

**Code:**

```python
fig = plt.figure(figsize=(13, 4.7), constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1, 1])
axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))
fig.suptitle("Figure 2 — Non-unique transcript landscape", fontsize=12, fontweight="bold")

# Panel A — non-uniqueness fraction by transcript-count bucket
bins   = [0, 1, 2, 3, 5, 10, 20, 50, 100, 1000]
labels = ["1", "2", "3", "4–5", "6–10", "11–20", "21–50", "51–100", "100+"]
ok["tx_bucket"] = pd.cut(ok["transcripts_compared_for_overlap"],
                          bins=bins, labels=labels, right=True, include_lowest=True)
agg = ok.groupby("tx_bucket", observed=True).agg(
    n=("input_symbol", "size"),
    n_nu=("has_non_unique_transcript", "sum"),
)
agg["frac"] = agg["n_nu"] / agg["n"]
bars = axA.bar(range(len(agg)), agg["frac"],
               color=sns.color_palette("viridis", len(agg)),
               edgecolor="white", linewidth=0.5)
axA.set_xticks(range(len(agg)))
axA.set_xticklabels(agg.index.astype(str), rotation=0)
axA.set_xlabel("Transcripts compared (binned)")
axA.set_ylabel("Fraction with non-unique transcript(s)")
axA.set_title(f"A. Fraction by transcript count")
axA.set_ylim(0, 1.05)
for bar, n in zip(bars, agg["n"]):
    axA.annotate(f"n={n:,}", xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                  xytext=(0, 3), textcoords="offset points",
                  ha="center", fontsize=7, color="dimgray")

# Panel B — scatter, log-log, with y=x diagonal
s = ok[ok["transcript_count"] > 0]
axB.scatter(s["transcript_count"], s["non_unique_transcript_count"] + 0.5,
            s=4, alpha=0.18, color=sns.color_palette("viridis", 5)[2], edgecolor="none")
axB.plot([1, s["transcript_count"].max()], [0.5, s["transcript_count"].max() + 0.5],
         ls="--", color="crimson", lw=0.8, label="y = x")
axB.set_xscale("log"); axB.set_yscale("log")
axB.set_xlabel("Transcripts per gene (log)")
axB.set_ylabel("Non-unique transcripts (+0.5)")
axB.set_title("B. Per-gene scaling")
axB.legend(loc="upper left", frameon=False)

# Panel C — joint classification heatmap
ct = pd.crosstab(ok["is_exon"], ok["has_non_unique_transcript"])
ct = ct.reindex(index=[True, False], columns=[True, False]).fillna(0).astype(int)
total_n = int(ct.values.sum())
pct = 100.0 * ct / max(total_n, 1)
annot = np.empty(ct.shape, dtype=object)
for i in range(ct.shape[0]):
    for j in range(ct.shape[1]):
        annot[i, j] = f"{ct.iat[i, j]:,}\n({pct.iat[i, j]:.1f}%)"
sns.heatmap(ct, annot=annot, fmt="", cmap="viridis", ax=axC,
            cbar_kws={"label": "gene count"}, linewidths=0.5, linecolor="white")
axC.set_xlabel("has_non_unique_transcript")
axC.set_ylabel("is_exon")
axC.set_title(f"C. Joint classification (n = {total_n:,})")

save_fig(fig, "fig2_non_uniqueness")
plt.show()
```

---

## 4. Figure 3 — Curated vs predicted contribution

**What this figure shows:** Panel A is a violin plot of `log(1 + xm_count)` partitioned by `has_non_unique_transcript`, with a Mann–Whitney U test reported in-panel. Panel B is a stacked composition of transcript types within five NM-count buckets.

**What insight it gives you (NM-only rule):** Even though XM transcripts are no longer in the overlap rule, **non-unique genes still carry significantly more XM** (mean 8.98 vs 3.35; Mann–Whitney p ≈ 0). This is now a pure correlation, not a confounder: highly-spliced loci attract more curated NM annotation *and* more XM prediction together, and both surface in the same gene class. Useful for arguing that the 50.7 % non-uniqueness rate reflects something biological about the underlying loci, not an annotation-pipeline artifact.

**Code:**

```python
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.7), constrained_layout=True)
fig.suptitle("Figure 3 — Curated (NM/NR) vs predicted (XM/XR) contributions",
             fontsize=12, fontweight="bold")

# Panel A — violinplot with M–W U test
pal = sns.color_palette("viridis", 2)
grp_t = ok.loc[ok["has_non_unique_transcript"], "xm_count"]
grp_f = ok.loc[~ok["has_non_unique_transcript"], "xm_count"]
parts = axA.violinplot([np.log1p(grp_f), np.log1p(grp_t)],
                       positions=[0, 1], widths=0.7,
                       showmeans=False, showmedians=True)
for pc, c in zip(parts["bodies"], pal):
    pc.set_facecolor(c); pc.set_alpha(0.7); pc.set_edgecolor("white")
axA.set_xticks([0, 1])
axA.set_xticklabels([f"has_non_unique = False\n(n = {len(grp_f):,})",
                     f"has_non_unique = True\n(n = {len(grp_t):,})"])
axA.set_ylabel("log(1 + XM transcripts per gene)")
axA.set_title("A. Predicted-transcript load by status")

from scipy.stats import mannwhitneyu
u, p = mannwhitneyu(grp_t, grp_f, alternative="two-sided")
axA.text(0.5, 0.96, f"Mann–Whitney U: p = {p:.2e}",
         transform=axA.transAxes, ha="center", va="top",
         fontsize=9, color="dimgray")

# Panel B — composition stack by NM bucket
bins_nm = [0, 1, 2, 4, 7, 100]
labels_nm = ["0 NM", "1 NM (is_exon)", "2–3 NM", "4–6 NM", "7+ NM"]
ok["nm_bucket"] = pd.cut(ok["nm_count"], bins=bins_nm, labels=labels_nm,
                         right=False, include_lowest=True)
grouped = ok.groupby("nm_bucket", observed=True)[["nm_count","nr_count","xm_count","xr_count"]].sum()
grouped_pct = grouped.div(grouped.sum(axis=1), axis=0) * 100
grouped_pct.plot(kind="bar", stacked=True, ax=axB,
                 color=sns.color_palette("viridis", 4),
                 edgecolor="white", width=0.85, legend=True)
axB.set_xlabel("Gene class (by NM count)")
axB.set_ylabel("Mean transcript-type share (%)")
axB.set_title("B. Composition by gene class")
axB.legend(["NM", "NR", "XM", "XR"], loc="lower right", frameon=False)
plt.setp(axB.get_xticklabels(), rotation=20, ha="right")

save_fig(fig, "fig3_curated_vs_predicted")
plt.show()
```

---

## 5. Figure 4 — Distributions (CDF + rank-frequency)

**What this figure shows:** Panel A is the empirical CDF of three per-gene quantities on a log-X axis. Panel B is a rank-frequency plot — a straight line on log-log is the Zipfian signature.

**What insight it gives you:** The per-gene transcript count is heavy-tailed. A small number of "transcript-rich" loci (BRCA1 with 368, several others over 100) dominate the total count. This motivates log axes throughout the rest of the figure suite and means standard statistics (mean, sd) are misleading — report medians and IQR.

**Code:**

```python
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
fig.suptitle("Figure 4 — Per-gene distributions", fontsize=12, fontweight="bold")
pal = sns.color_palette("viridis", 4)

# Panel A — empirical CDFs
for i, (col, label) in enumerate([
    ("transcript_count", "Total transcripts"),
    ("nm_count", "NM only"),
    ("non_unique_transcript_count", "Non-unique transcripts"),
]):
    x = np.sort(ok[col].to_numpy())
    y = np.arange(1, len(x) + 1) / len(x)
    axA.plot(np.clip(x, 0.5, None), y, color=pal[i], label=label, lw=1.8)
axA.set_xscale("log")
axA.set_xlabel("Count per gene (log)")
axA.set_ylabel("Cumulative fraction of OK genes")
axA.set_title("A. Cumulative distributions")
axA.legend(frameon=False, loc="lower right")
axA.set_ylim(0, 1.02)

# Panel B — rank-frequency
counts = np.sort(ok["transcript_count"].to_numpy())[::-1]
ranks = np.arange(1, len(counts) + 1)
axB.scatter(ranks, counts, s=2, alpha=0.4, color=pal[2], edgecolor="none")
axB.set_xscale("log"); axB.set_yscale("log")
axB.set_xlabel("Gene rank (log)")
axB.set_ylabel("Transcripts per gene (log)")
axB.set_title("B. Rank–frequency (Zipf-style)")

save_fig(fig, "fig4_distributions")
plt.show()
```

---

## 6. Figure 5 — Worked examples (exon-level tracks)

**What this figure shows:** The actual GRCh38 exon coordinates for three real genes — ACTB (clean), MYC (both non-unique), and GAPDH (5 of 6 non-unique). Exons are colored **red** if the transcript is flagged non-unique, **blue** otherwise.

**What insight it gives you:** This is the figure to put first when introducing the rule. It makes the abstract pairwise-overlap definition concrete: in MYC you can see at a glance that both transcripts hit every exon; in GAPDH you can see that `NM_001256799.3` (top, blue) is unique because its 5′ start is downstream of all the others.

**Code:**

```python
# Pulls exon data straight from the cached API responses
import sys
sys.path.insert(0, str(PROJECT))
from ncbi_exon.parser import _grch38_exons, parse_product_report

cache = PROJECT / "cache"
examples = [
    ("ACTB",  "Single-transcript gene (trivially clean)"),
    ("MYC",   "Two transcripts, both non-unique"),
    ("GAPDH", "Six transcripts (5 of 6 non-unique)"),
]
fig, axes = plt.subplots(len(examples), 1,
                          figsize=(11, 1.6 * len(examples) + 1.2),
                          constrained_layout=True)
fig.suptitle("Figure 5 — Worked examples: exon-level overlap (GRCh38)",
             fontsize=12, fontweight="bold")

for ax, (sym, subtitle) in zip(axes, examples):
    payload = json.loads((cache / f"{sym}__human.json").read_text())
    prod = payload["reports"][0]["product"]
    row = parse_product_report(sym, payload)
    nu = set(row.non_unique_transcript_accessions.split(";")) if row.non_unique_transcript_accessions else set()

    rows_data = []
    for t in prod.get("transcripts", []):
        ex = _grch38_exons(t)
        if ex:
            rows_data.append((t["accession_version"], ex))

    starts = [e[0] for _, exs in rows_data for e in exs]
    ends   = [e[1] for _, exs in rows_data for e in exs]
    x0, x1 = min(starts), max(ends)
    span = x1 - x0
    ax.set_xlim(x0 - 0.02 * span, x1 + 0.02 * span)
    ax.set_ylim(-0.5, len(rows_data) - 0.5)

    for i, (acc, exs) in enumerate(rows_data):
        y = len(rows_data) - 1 - i
        color = "crimson" if acc in nu else "#2b8cbe"
        ax.hlines(y, min(e[0] for e in exs), max(e[1] for e in exs),
                  color="lightgray", lw=0.8, zorder=1)
        for b, e in exs:
            ax.add_patch(plt.Rectangle((b, y - 0.3), e - b, 0.6,
                                        facecolor=color, edgecolor="none"))
        ax.text(x0 - 0.025 * span, y, acc, ha="right", va="center", fontsize=8,
                color=color, fontweight="bold" if acc in nu else "normal")

    ax.set_yticks([])
    chrom = prod["transcripts"][0]["genomic_locations"][0]["genomic_accession_version"]
    ax.set_xlabel(f"GRCh38 position on {chrom}", fontsize=9)
    ax.set_title(f"{sym} (gene_id {prod['gene_id']}) — {subtitle}", loc="left")

from matplotlib.patches import Patch
fig.legend(handles=[
    Patch(facecolor="#2b8cbe", label="transcript NOT flagged non-unique"),
    Patch(facecolor="crimson", label="transcript flagged non-unique"),
], loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.02))

save_fig(fig, "fig5_exon_tracks")
plt.show()
```

---

## 7. Figure 6 — Classification flow

This figure has two versions. **Use the Sankey** for the thesis itself — it's the cleanest way to show proportional flow and Plotly produces publication-quality PDF/SVG. The matplotlib version is included only as a no-extra-dependency fallback.

### 7a. Plotly Sankey (recommended)

**What this figure shows:** A flow diagram from "all input symbols" through API status → `is_exon` → `has_non_unique_transcript`, with each ribbon's thickness proportional to gene count.

**What insight it gives you:** The single picture that summarizes the entire pipeline outcome. Suitable as the thesis chapter's opening figure or for a 1-page methods summary.

Install once: `.venv/bin/pip install plotly kaleido` (kaleido is the headless renderer for static export).

```python
import plotly.graph_objects as go

total = len(df)
n_ok = (df["status"] == "OK").sum()
n_nf = (df["status"] == "NOT_FOUND").sum()
n_mm = (df["status"] == "MULTI_MATCH").sum()
n_isexon_t = int((ok["is_exon"]).sum())
n_isexon_f = int((~ok["is_exon"]).sum())
n_tu  = int((ok["is_exon"] & ~ok["has_non_unique_transcript"]).sum())
n_tnu = int((ok["is_exon"] & ok["has_non_unique_transcript"]).sum())
n_fu  = int((~ok["is_exon"] & ~ok["has_non_unique_transcript"]).sum())
n_fnu = int((~ok["is_exon"] & ok["has_non_unique_transcript"]).sum())

labels = [
    f"All input symbols<br>{total:,}",                                    # 0
    f"OK<br>{n_ok:,}",                                                    # 1
    f"NOT_FOUND<br>{n_nf}",                                               # 2
    f"MULTI_MATCH<br>{n_mm}",                                             # 3
    f"is_exon = True<br>{n_isexon_t:,}",                                  # 4
    f"is_exon = False<br>{n_isexon_f:,}",                                 # 5
    f"clean single-canonical<br>{n_tu:,} ({100*n_tu/total:.1f}%)",        # 6
    f"single NM, redundant<br>{n_tnu:,} ({100*n_tnu/total:.1f}%)",        # 7
    f"multi NM, unique exons<br>{n_fu:,} ({100*n_fu/total:.1f}%)",        # 8
    f"multi NM, redundant<br>{n_fnu:,} ({100*n_fnu/total:.1f}%)",         # 9
]
source = [0, 0, 0, 1, 1, 4, 4, 5, 5]
target = [1, 2, 3, 4, 5, 6, 7, 8, 9]
value  = [n_ok, n_nf, n_mm, n_isexon_t, n_isexon_f, n_tu, n_tnu, n_fu, n_fnu]

fig = go.Figure(go.Sankey(
    arrangement="snap",
    node=dict(
        label=labels, pad=18, thickness=22,
        color=["#440154", "#3b528b", "#999", "#999",
               "#21918c", "#5ec962",
               "#1f9e89", "#fde725", "#1f9e89", "#fde725"],
    ),
    link=dict(
        source=source, target=target, value=value,
        color=["rgba(59,82,139,0.4)"] * 3
              + ["rgba(33,145,140,0.4)", "rgba(94,201,98,0.4)"]
              + ["rgba(31,158,137,0.5)", "rgba(253,231,37,0.5)"] * 2,
    ),
))
fig.update_layout(
    title="Figure 6 — Classification flow (Sankey)",
    title_font_size=14,
    font=dict(size=11),
    width=1100, height=550,
    margin=dict(l=10, r=10, t=50, b=10),
)

# Save as PNG (raster), PDF (vector), and HTML (interactive)
fig.write_image(FIG_DIR / "fig6_sankey.png", scale=2)
fig.write_image(FIG_DIR / "fig6_sankey.pdf")
fig.write_html(FIG_DIR / "fig6_sankey.html")
fig.show()
```

### 7b. Matplotlib fallback (no Plotly required)

A four-stage proportional-bar diagram. Less elegant than the Sankey but works without any extra dependency.

```python
fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
ax.axis("off")
fig.suptitle("Figure 6 — Classification flow", fontsize=12, fontweight="bold")

total = len(df)
n_ok = (df["status"] == "OK").sum()
n_nf = (df["status"] == "NOT_FOUND").sum()
n_mm = (df["status"] == "MULTI_MATCH").sum()
n_te = int((ok["is_exon"]).sum())
n_fe = int((~ok["is_exon"]).sum())
n_tu  = int((ok["is_exon"] & ~ok["has_non_unique_transcript"]).sum())
n_tnu = int((ok["is_exon"] & ok["has_non_unique_transcript"]).sum())
n_fu  = int((~ok["is_exon"] & ~ok["has_non_unique_transcript"]).sum())
n_fnu = int((~ok["is_exon"] & ok["has_non_unique_transcript"]).sum())

pal = sns.color_palette("viridis", 6)
def bar(x, y0, h, w, color, label):
    ax.add_patch(plt.Rectangle((x, y0), w, h, facecolor=color, edgecolor="white"))
    ax.text(x + w / 2, y0 + h / 2, label, ha="center", va="center", color="white",
            fontsize=9, fontweight="bold")

bar(0.05, 0.1, 0.8, 0.12, pal[0], f"All input\n{total:,}\n100.0%")
yc = 0.1
for label, n, c in [("OK", n_ok, pal[1]), ("NOT_FOUND", n_nf, "lightgray"),
                    ("MULTI_MATCH", n_mm, "lightgray")]:
    if n == 0: continue
    h = 0.8 * n / total
    bar(0.27, yc, h, 0.12, c, f"{label}\n{n:,}\n{100 * n / total:.2f}%")
    yc += h
yc = 0.1
for label, n, c in [("is_exon=True", n_te, pal[2]), ("is_exon=False", n_fe, pal[3])]:
    h = 0.8 * (n / n_ok) * (n_ok / total)
    bar(0.49, yc, h, 0.12, c, f"{label}\n{n:,}\n{100*n/n_ok:.1f}% of OK")
    yc += h
yc = 0.1
h_te = 0.8 * (n_te / total)
bar(0.71, yc, 0.8 * n_tu / total, 0.12, pal[4],
    f"non_unique=F\n{n_tu:,}\n{100*n_tu/n_te:.1f}%")
bar(0.71, yc + 0.8 * n_tu / total, 0.8 * n_tnu / total, 0.12, pal[5],
    f"non_unique=T\n{n_tnu:,}\n{100*n_tnu/n_te:.1f}%")
yc += h_te
bar(0.71, yc, 0.8 * n_fu / total, 0.12, pal[4],
    f"non_unique=F\n{n_fu:,}\n{100*n_fu/n_fe:.1f}%")
bar(0.71, yc + 0.8 * n_fu / total, 0.8 * n_fnu / total, 0.12, pal[5],
    f"non_unique=T\n{n_fnu:,}\n{100*n_fnu/n_fe:.1f}%")

for x, label in [(0.05, "Stage 1\nInput"), (0.27, "Stage 2\nAPI status"),
                  (0.49, "Stage 3\nis_exon (NM=1)"),
                  (0.71, "Stage 4\nhas_non_unique_transcript")]:
    ax.text(x + 0.06, 0.95, label, ha="center", va="bottom",
             fontsize=10, fontweight="bold", color="#333")

ax.set_xlim(0, 0.92); ax.set_ylim(0, 1)
save_fig(fig, "fig6_classification_flow")
plt.show()
```

---

## 8. Best practices for thesis-quality figures

| | Choice | Why |
|---|---|---|
| Color | `viridis` (or `cividis` for colorblind audiences) | Perceptually uniform; survives B/W printing |
| Layout | `constrained_layout=True` on every `subplots()` | Auto-resolves title/label overlap (fixes the issue you saw) |
| Export | Save **both** PNG (300 dpi raster) **and** PDF (vector) | PDFs scale losslessly into a typeset thesis; PNGs embed in Word/Google Docs |
| Axes | Use log axes whenever the data span > 1 decade | Required for our power-law-like distributions |
| Fonts | Default sans-serif at 10 pt is correct for a 5.5" figure | Avoid the temptation to make labels small to fit more |
| Spines | Hide top & right (`axes.spines.top/right`) | Reduces chart junk (Tufte) |
| Stats annotations | Always include n and a test p-value where comparing groups | Reviewers will ask if you don't |

---

## 9. Why some of your earlier outputs overlapped

The original `scripts/make_figures.py` mixed `fig.suptitle(..., y=1.04)` with hard-coded annotation y-positions that landed in the same vertical strip — matplotlib draws annotations relative to the **axes**, but `suptitle` is in **figure** coordinates, so they collided whenever the data extended near the top of the panel (as in TP53's NM histogram).

The snippets above fix this by:

1. Using `constrained_layout=True` (matplotlib reserves real space for the suptitle, instead of overdrawing).
2. Removing the hard-coded `y=1.04` from `suptitle` — `constrained_layout` handles it.
3. Moving the Figure 1B annotation from `xytext=(4, n_is_exon * 0.9)` to `xytext=(6, n_is_exon * 0.35)` so it sits *below* the bar's top, not above it.

If you want to verify, run the regeneration script and compare against `output/figures/`:

```bash
.venv/bin/python scripts/make_figures.py
```

---

## 10. Recommended workflow summary

1. Run the **Setup cell** once per notebook session.
2. For each figure, paste the cell, run, inspect inline, tweak.
3. Each cell saves a 300 dpi PNG + a vector PDF into `output/figures/`.
4. For the thesis PDF, pull the PDFs — they scale losslessly. For slides, pull the PNGs.
5. For Figure 6 specifically, prefer the **Plotly Sankey** version — it embeds as an interactive HTML supplement to your thesis or as a static PDF in the print version.

Done.
