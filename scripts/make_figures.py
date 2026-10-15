"""Generate the thesis figure suite from a v2 CSV.

Produces 6 figures (PNG @ 300 dpi + vector PDF) plus an interactive Plotly
Sankey (PNG + PDF + HTML) into the chosen output directory.

Usage:
    .venv/bin/python scripts/make_figures.py                    # full 20k by default
    .venv/bin/python scripts/make_figures.py \\
        --input output/test50_v2.csv \\
        --output-dir output/figures_test50

Layouts use `constrained_layout=True` and explicit annotation positions so the
suptitle never collides with axis content (the issue in earlier output). See
`How_to_generate_figures.md` for the per-figure rationale.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

# Which overlap definition drives Figures 1–6. Definition #1 = partial overlap
# (share >=1 bp). Definition #2 = containment (every exon fully inside a
# partner exon). Set once in main() via set_definition(). The def-comparison
# figure always reads BOTH and ignores this global.
DEFINITION = 1
DEF_SUFFIX = ""


def set_definition(d: int) -> None:
    global DEFINITION, DEF_SUFFIX
    DEFINITION = d
    DEF_SUFFIX = "" if d == 1 else "   [Definition #2 — containment]"


def remap_definition(df: pd.DataFrame) -> pd.DataFrame:
    """Return a frame whose Definition-#1 column names carry the active
    definition's values, so the figure code can stay definition-agnostic.

    For Definition #2 the three `_containment` columns overwrite the base
    names. No-op for Definition #1.
    """
    if DEFINITION != 2:
        return df
    df = df.copy()
    df["has_non_unique_transcript"] = df["has_non_unique_transcript_containment"]
    df["non_unique_transcript_count"] = df["non_unique_transcript_count_containment"]
    df["non_unique_transcript_accessions"] = df["non_unique_transcript_accessions_containment"]
    return df


# ── style ───────────────────────────────────────────────────────────────
def apply_style() -> None:
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


def load_df(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    for c in ("is_exon", "has_non_unique_transcript", "has_non_unique_transcript_containment"):
        if c in df.columns and df[c].dtype != bool:
            df[c] = df[c].map({"True": True, "False": False, True: True, False: False}).fillna(False).astype(bool)
    return df


def save(fig: plt.Figure, name: str, out_dir: Path) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"{name}.{ext}")
    plt.close(fig)
    print(f"  saved: {name}.{{png,pdf}}")


# ── Figure 1 — composition + NM distribution ────────────────────────────
def fig1(df: pd.DataFrame, ok: pd.DataFrame, out_dir: Path) -> None:
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    fig.suptitle("Figure 1 — Transcript-type composition", fontsize=12, fontweight="bold")

    totals = {
        "NM (curated mRNA)":    ok["nm_count"].sum(),
        "NR (curated ncRNA)":   ok["nr_count"].sum(),
        "XM (predicted mRNA)":  ok["xm_count"].sum(),
        "XR (predicted ncRNA)": ok["xr_count"].sum(),
    }
    colors = sns.color_palette("viridis", 4)
    left = 0
    for (label, val), c in zip(totals.items(), colors):
        axA.barh([0], val, left=left, color=c, label=f"{label}: {int(val):,}")
        left += val
    axA.set_yticks([])
    axA.set_xlabel("Total transcripts (all OK genes pooled)")
    axA.set_title("A. RefSeq transcript composition")
    axA.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18),
               ncol=2, frameon=False, fontsize=8)
    for s in ("left", "right", "top"):
        axA.spines[s].set_visible(False)

    n_exon = int((ok["nm_count"] == 1).sum())
    max_nm = max(int(ok["nm_count"].max()) if len(ok) else 5, 5)
    bins = np.arange(0, min(max_nm + 2, 22)) - 0.5
    axB.hist(ok["nm_count"], bins=bins, color=colors[0],
             edgecolor="white", linewidth=0.6)
    axB.axvline(1, color="crimson", ls="--", lw=1, label="nm_count = 1 (is_mono-isoform)")
    if n_exon > 0:
        # IMPORTANT: keep xytext inside the visible xlim ([-0.5, 20.5]) — placing
        # it past the cap would cause `constrained_layout` to expand the canvas
        # to fit the off-axis annotation, producing an absurdly wide figure.
        axB.annotate(
            f"is_mono-isoform = True\nn = {n_exon:,}  ({100 * n_exon / max(len(ok), 1):.1f}%)",
            xy=(1, n_exon),
            xytext=(6, max(n_exon * 0.35, 1.5)),
            arrowprops=dict(arrowstyle="->", color="crimson", lw=0.8),
            fontsize=9, color="crimson",
        )
    axB.set_yscale("log")
    axB.set_xlim(-0.5, min(max_nm + 1, 20.5))
    axB.set_xlabel("NM transcripts per gene")
    axB.set_ylabel("Number of genes (log)")
    axB.set_title(f"B. NM distribution (n = {len(ok):,})")
    axB.legend(loc="upper right", frameon=False)

    save(fig, "fig1_composition", out_dir)


# ── Figure 2 — non-uniqueness landscape ─────────────────────────────────
def fig2(df: pd.DataFrame, ok: pd.DataFrame, out_dir: Path) -> None:
    fig = plt.figure(figsize=(13, 4.7), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1, 1])
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))
    fig.suptitle("Figure 2 — Non-unique transcript landscape" + DEF_SUFFIX, fontsize=12, fontweight="bold")

    # NM-only consistency: both panels use nm_count, not the mixed transcript_count.
    # Bins are constructed so each label maps to exactly the NM count(s) it claims.
    # With right=True the interval (a, b] contains b, so 0 is its own bin via [-1, 0].
    bins = [-1, 0, 1, 2, 3, 5, 10, 20, 50, 100, 1000]
    labels = ["0", "1", "2", "3", "4–5", "6–10", "11–20", "21–50", "51–100", "100+"]
    tmp = ok.copy()
    tmp["nm_bucket"] = pd.cut(tmp["nm_count"], bins=bins, labels=labels, right=True)
    agg = tmp.groupby("nm_bucket", observed=True).agg(
        n=("input_symbol", "size"),
        n_nu=("has_non_unique_transcript", "sum"),
    )
    agg["frac"] = agg["n_nu"] / agg["n"]
    bars = axA.bar(range(len(agg)), agg["frac"],
                   color=sns.color_palette("viridis", max(len(agg), 1)),
                   edgecolor="white", linewidth=0.5)
    axA.set_xticks(range(len(agg)))
    axA.set_xticklabels(agg.index.astype(str), rotation=0)
    axA.set_xlabel("NM transcripts per gene (binned)")
    axA.set_ylabel("Fraction with non-unique NM transcript(s)")
    axA.set_title("A. Fraction by NM count")
    axA.set_ylim(0, 1.08)
    for bar, n in zip(bars, agg["n"]):
        axA.annotate(f"n={int(n):,}",
                     xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                     xytext=(0, 3), textcoords="offset points",
                     ha="center", fontsize=7, color="dimgray")

    # Panel B — now uses nm_count on the x-axis for NM-only consistency.
    s = ok[ok["nm_count"] > 0]
    if len(s):
        axB.scatter(s["nm_count"], s["non_unique_transcript_count"] + 0.5,
                    s=8, alpha=0.35, color=sns.color_palette("viridis", 5)[2],
                    edgecolor="none")
        xmax = max(int(s["nm_count"].max()), 2)
        axB.plot([1, xmax], [0.5, xmax + 0.5], ls="--", color="crimson", lw=0.8,
                 label="y = x (every NM is non-unique)")
    axB.set_xscale("log"); axB.set_yscale("log")
    axB.set_xlabel("NM transcripts per gene (log)")
    axB.set_ylabel("Non-unique NM transcripts (+0.5)")
    axB.set_title("B. Per-gene scaling (NM only)")
    axB.legend(loc="upper left", frameon=False)

    ct = pd.crosstab(ok["is_exon"], ok["has_non_unique_transcript"])
    ct = ct.reindex(index=[True, False], columns=[True, False]).fillna(0).astype(int)
    total_n = max(int(ct.values.sum()), 1)
    pct = 100.0 * ct / total_n
    annot = np.empty(ct.shape, dtype=object)
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            annot[i, j] = f"{int(ct.iat[i, j]):,}\n({pct.iat[i, j]:.1f}%)"
    sns.heatmap(ct, annot=annot, fmt="", cmap="viridis", ax=axC,
                cbar_kws={"label": "gene count"},
                linewidths=0.5, linecolor="white")
    axC.set_xlabel("has_non_unique_transcript")
    axC.set_ylabel("is_mono-isoform")
    axC.set_title(f"C. Joint classification (n = {total_n:,})")

    save(fig, "fig2_non_uniqueness", out_dir)


# ── Figure 3 — NM-only redundancy structure ─────────────────────────────
def fig3(df: pd.DataFrame, ok: pd.DataFrame, out_dir: Path) -> None:
    """NM-only redesign: how does NM redundancy distribute among multi-NM genes?
    Panel A asks "do flagged genes carry more NMs than unflagged multi-NM genes?"
    Panel B asks "within a flagged gene, what fraction of its NMs are non-unique?"
    """
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.7), constrained_layout=True)
    fig.suptitle("Figure 3 — Redundancy structure among NM transcripts (NM-only)" + DEF_SUFFIX,
                 fontsize=12, fontweight="bold")

    # ─ Panel A — NM count distribution restricted to multi-NM (the eligible set) ─
    pal = sns.color_palette("viridis", 2)
    multi = ok[ok["nm_count"] >= 2]
    grp_f = multi.loc[~multi["has_non_unique_transcript"], "nm_count"]
    grp_t = multi.loc[multi["has_non_unique_transcript"], "nm_count"]

    if len(grp_t) >= 2 and len(grp_f) >= 2:
        parts = axA.violinplot([np.log1p(grp_f), np.log1p(grp_t)],
                                positions=[0, 1], widths=0.7,
                                showmeans=False, showmedians=True)
        for pc, c in zip(parts["bodies"], pal):
            pc.set_facecolor(c); pc.set_alpha(0.7); pc.set_edgecolor("white")
        try:
            from scipy.stats import mannwhitneyu
            u, p = mannwhitneyu(grp_t, grp_f, alternative="two-sided")
            axA.text(0.5, 0.96,
                     f"Mann–Whitney U: p = {p:.2e}\nmedian NM: "
                     f"F={int(grp_f.median())}, T={int(grp_t.median())}",
                     transform=axA.transAxes, ha="center", va="top",
                     fontsize=9, color="dimgray")
        except ImportError:
            pass
    else:
        axA.text(0.5, 0.5, "Not enough samples\nfor violin (need ≥2 per group)",
                 transform=axA.transAxes, ha="center", va="center", color="dimgray")
    axA.set_xticks([0, 1])
    axA.set_xticklabels([f"has_non_unique = False\n(n = {len(grp_f):,})",
                         f"has_non_unique = True\n(n = {len(grp_t):,})"])
    axA.set_ylabel("log(1 + NM transcripts per gene)")
    axA.set_title("A. NM load among multi-NM genes (nm_count ≥ 2)")

    # ─ Panel B — within-NM redundancy share by gene class ─
    # For each NM bucket: % of NMs that are non-unique vs unique inside that bucket.
    bins_nm   = [-1, 1, 2, 3, 5, 10, 20, 1000]
    labels_nm = ["1", "2", "3", "4–5", "6–10", "11–20", "21+"]
    tmp = ok[ok["nm_count"] >= 1].copy()
    tmp["nm_bucket"] = pd.cut(tmp["nm_count"], bins=bins_nm, labels=labels_nm, right=True)
    grouped = tmp.groupby("nm_bucket", observed=True).agg(
        total_nm=("nm_count", "sum"),
        total_nu=("non_unique_transcript_count", "sum"),
        n_genes=("nm_count", "size"),
    )
    grouped["unique_nm"] = grouped["total_nm"] - grouped["total_nu"]
    grouped["pct_non_unique"] = 100 * grouped["total_nu"] / grouped["total_nm"]
    grouped["pct_unique"] = 100 - grouped["pct_non_unique"]

    nu_color = sns.color_palette("viridis", 5)[4]   # bright yellow-green
    u_color  = sns.color_palette("viridis", 5)[0]   # deep purple
    x_pos = np.arange(len(grouped))
    axB.bar(x_pos, grouped["pct_unique"], color=u_color,
            edgecolor="white", linewidth=0.5, label="unique NMs")
    axB.bar(x_pos, grouped["pct_non_unique"], bottom=grouped["pct_unique"],
            color=nu_color, edgecolor="white", linewidth=0.5, label="non-unique NMs")
    axB.set_xticks(x_pos)
    axB.set_xticklabels(grouped.index.astype(str))
    axB.set_xlabel("Gene class (by NM count)")
    axB.set_ylabel("Share of NM transcripts in class (%)")
    axB.set_title("B. NM redundancy share by gene class")
    axB.set_ylim(0, 105)
    axB.legend(loc="lower right", frameon=False)
    for xi, n, total in zip(x_pos, grouped["n_genes"], grouped["total_nm"]):
        axB.annotate(f"n={int(n):,}\nNM={int(total):,}",
                     xy=(xi, 101), ha="center", va="bottom",
                     fontsize=7, color="dimgray")

    save(fig, "fig3_nm_redundancy_structure", out_dir)


# ── Figure 4 — distributions ────────────────────────────────────────────
def fig4(df: pd.DataFrame, ok: pd.DataFrame, out_dir: Path) -> None:
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    fig.suptitle("Figure 4 — Per-gene distributions" + DEF_SUFFIX, fontsize=12, fontweight="bold")
    pal = sns.color_palette("viridis", 4)

    for i, (col, label) in enumerate([
        ("transcript_count", "Total transcripts"),
        ("nm_count", "NM only"),
        ("non_unique_transcript_count", "Non-unique transcripts"),
    ]):
        x = np.sort(ok[col].to_numpy())
        if len(x) == 0:
            continue
        y = np.arange(1, len(x) + 1) / len(x)
        axA.plot(np.clip(x, 0.5, None), y, color=pal[i], label=label, lw=1.8)
    axA.set_xscale("log")
    axA.set_xlabel("Count per gene (log)")
    axA.set_ylabel("Cumulative fraction of OK genes")
    axA.set_title("A. Cumulative distributions")
    axA.legend(frameon=False, loc="lower right")
    axA.set_ylim(0, 1.02)

    counts = np.sort(ok["transcript_count"].to_numpy())[::-1]
    ranks = np.arange(1, len(counts) + 1)
    axB.scatter(ranks, counts, s=8, alpha=0.5, color=pal[2], edgecolor="none")
    axB.set_xscale("log"); axB.set_yscale("log")
    axB.set_xlabel("Gene rank (log)")
    axB.set_ylabel("Transcripts per gene (log)")
    axB.set_title("B. Rank–frequency (Zipf-style)")

    save(fig, "fig4_distributions", out_dir)


# ── Figure 5 — exon tracks for ACTB / MYC / GAPDH ───────────────────────
def fig5(df: pd.DataFrame, ok: pd.DataFrame, out_dir: Path) -> None:
    from ncbi_exon.parser import _grch38_exons, parse_product_report

    cache = PROJECT_ROOT / "cache"
    if DEFINITION == 2:
        examples = [
            ("ACTB",  "Single transcript (trivially clean)"),
            ("MYC",   "Two transcripts — neither contained → 0 non-unique under #2"),
            ("GAPDH", "Five NM — 3 of 5 non-unique under #2"),
        ]
    else:
        examples = [
            ("ACTB",  "Single-transcript gene (trivially clean)"),
            ("MYC",   "Two transcripts, both non-unique"),
            ("GAPDH", "Six transcripts (5 of 6 non-unique)"),
        ]
    # Filter to genes whose cache files actually exist
    examples = [(s, sub) for (s, sub) in examples if (cache / f"{s}__human.json").exists()]
    if not examples:
        print("  fig5 skipped — no cached examples")
        return

    fig, axes = plt.subplots(len(examples), 1,
                              figsize=(11, 1.6 * len(examples) + 1.2),
                              constrained_layout=True)
    if len(examples) == 1:
        axes = [axes]
    fig.suptitle("Figure 5 — Worked examples: exon-level overlap (GRCh38)" + DEF_SUFFIX,
                 fontsize=12, fontweight="bold")

    for ax, (sym, subtitle) in zip(axes, examples):
        payload = json.loads((cache / f"{sym}__human.json").read_text())
        prod = payload["reports"][0]["product"]
        row = parse_product_report(sym, payload)
        nu_accs = (row.non_unique_transcript_accessions_containment if DEFINITION == 2
                   else row.non_unique_transcript_accessions)
        nu = set(nu_accs.split(";")) if nu_accs else set()

        # NM-only: the overlap rule operates exclusively on curated mRNAs, so the
        # track plot mirrors the candidate set and omits NR/XM/XR transcripts.
        rows_data = []
        for t in prod.get("transcripts", []):
            acc = t.get("accession_version", "")
            if not acc.startswith("NM_"):
                continue
            ex = _grch38_exons(t)
            if ex:
                rows_data.append((acc, ex))

        if not rows_data:
            ax.text(0.5, 0.5, "no GRCh38 data", transform=ax.transAxes, ha="center")
            continue

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
            ax.text(x0 - 0.025 * span, y, acc, ha="right", va="center",
                     fontsize=8, color=color,
                     fontweight="bold" if acc in nu else "normal")

        ax.set_yticks([])
        chrom = prod["transcripts"][0]["genomic_locations"][0]["genomic_accession_version"]
        ax.set_xlabel(f"GRCh38 position on {chrom}", fontsize=9)
        ax.set_title(f"{sym} (gene_id {prod['gene_id']}) — {subtitle}", loc="left")

    from matplotlib.patches import Patch
    fig.legend(handles=[
        Patch(facecolor="#2b8cbe", label="transcript NOT flagged non-unique"),
        Patch(facecolor="crimson", label="transcript flagged non-unique"),
    ], loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.02))

    save(fig, "fig5_exon_tracks", out_dir)


# ── Figure 6 — classification flow (matplotlib fallback) ────────────────
def fig6_matplotlib(df: pd.DataFrame, ok: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    ax.axis("off")
    fig.suptitle("Figure 6 — Classification flow" + DEF_SUFFIX, fontsize=12, fontweight="bold")

    total = len(df)
    n_ok = int((df["status"] == "OK").sum())
    n_nf = int((df["status"] == "NOT_FOUND").sum())
    n_mm = int((df["status"] == "MULTI_MATCH").sum())
    n_te = int(ok["is_exon"].sum())
    n_fe = int((~ok["is_exon"]).sum())
    n_tu  = int((ok["is_exon"] & ~ok["has_non_unique_transcript"]).sum())
    n_tnu = int((ok["is_exon"] & ok["has_non_unique_transcript"]).sum())
    n_fu  = int((~ok["is_exon"] & ~ok["has_non_unique_transcript"]).sum())
    n_fnu = int((~ok["is_exon"] & ok["has_non_unique_transcript"]).sum())

    pal = sns.color_palette("viridis", 6)
    def bar(x, y0, h, w, color, label):
        ax.add_patch(plt.Rectangle((x, y0), w, h, facecolor=color, edgecolor="white"))
        ax.text(x + w / 2, y0 + h / 2, label, ha="center", va="center",
                 color="white", fontsize=9, fontweight="bold")

    bar(0.05, 0.1, 0.8, 0.12, pal[0], f"All input\n{total:,}\n100.0%")
    yc = 0.1
    for label, n, c in [("OK", n_ok, pal[1]),
                         ("NOT_FOUND", n_nf, "lightgray"),
                         ("MULTI_MATCH", n_mm, "lightgray")]:
        if n == 0: continue
        h = 0.8 * n / max(total, 1)
        bar(0.27, yc, h, 0.12, c, f"{label}\n{n:,}\n{100*n/total:.2f}%")
        yc += h
    yc = 0.1
    if n_ok > 0:
        for label, n, c in [("is_mono-isoform=True", n_te, pal[2]),
                             ("is_mono-isoform=False", n_fe, pal[3])]:
            h = 0.8 * (n / n_ok) * (n_ok / total)
            bar(0.49, yc, h, 0.12, c, f"{label}\n{n:,}\n{100*n/n_ok:.1f}% of OK")
            yc += h
    yc = 0.1
    h_te = 0.8 * (n_te / max(total, 1))
    if n_te > 0:
        bar(0.71, yc, 0.8 * n_tu / total, 0.12, pal[4],
            f"non_unique=F\n{n_tu:,}\n{100*n_tu/n_te:.1f}%")
        bar(0.71, yc + 0.8 * n_tu / total, 0.8 * n_tnu / total, 0.12, pal[5],
            f"non_unique=T\n{n_tnu:,}\n{100*n_tnu/n_te:.1f}%")
    yc += h_te
    if n_fe > 0:
        bar(0.71, yc, 0.8 * n_fu / total, 0.12, pal[4],
            f"non_unique=F\n{n_fu:,}\n{100*n_fu/n_fe:.1f}%")
        bar(0.71, yc + 0.8 * n_fu / total, 0.8 * n_fnu / total, 0.12, pal[5],
            f"non_unique=T\n{n_fnu:,}\n{100*n_fnu/n_fe:.1f}%")

    for x, label in [(0.05, "Stage 1\nInput"), (0.27, "Stage 2\nAPI status"),
                      (0.49, "Stage 3\nis_mono-isoform (NM=1)"),
                      (0.71, "Stage 4\nhas_non_unique_transcript")]:
        ax.text(x + 0.06, 0.95, label, ha="center", va="bottom",
                 fontsize=10, fontweight="bold", color="#333")

    ax.set_xlim(0, 0.92); ax.set_ylim(0, 1)
    save(fig, "fig6_classification_flow", out_dir)


# ── Figure 6 — Plotly Sankey (recommended) ──────────────────────────────
def fig6_sankey(df: pd.DataFrame, ok: pd.DataFrame, out_dir: Path) -> None:
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("  fig6 sankey skipped — install plotly + kaleido")
        return

    total = len(df)
    n_ok = int((df["status"] == "OK").sum())
    n_nf = int((df["status"] == "NOT_FOUND").sum())
    n_mm = int((df["status"] == "MULTI_MATCH").sum())
    n_te = int(ok["is_exon"].sum())
    n_fe = int((~ok["is_exon"]).sum())
    n_tu  = int((ok["is_exon"] & ~ok["has_non_unique_transcript"]).sum())
    n_tnu = int((ok["is_exon"] & ok["has_non_unique_transcript"]).sum())
    n_fu  = int((~ok["is_exon"] & ~ok["has_non_unique_transcript"]).sum())
    n_fnu = int((~ok["is_exon"] & ok["has_non_unique_transcript"]).sum())

    labels = [
        f"All input symbols<br>{total:,}",
        f"OK<br>{n_ok:,}",
        f"NOT_FOUND<br>{n_nf}",
        f"MULTI_MATCH<br>{n_mm}",
        f"is_mono-isoform = True<br>{n_te:,}",
        f"is_mono-isoform = False<br>{n_fe:,}",
        f"clean single-canonical<br>{n_tu:,} ({100*n_tu/max(total,1):.1f}%)",
        f"single NM, redundant<br>{n_tnu:,} ({100*n_tnu/max(total,1):.1f}%)",
        f"multi NM, unique exons<br>{n_fu:,} ({100*n_fu/max(total,1):.1f}%)",
        f"multi NM, redundant<br>{n_fnu:,} ({100*n_fnu/max(total,1):.1f}%)",
    ]
    source = [0, 0, 0, 1, 1, 4, 4, 5, 5]
    target = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    value  = [n_ok, n_nf, n_mm, n_te, n_fe, n_tu, n_tnu, n_fu, n_fnu]
    # Drop any zero links so Plotly doesn't render empty ghosts
    keep = [i for i, v in enumerate(value) if v > 0]
    source = [source[i] for i in keep]
    target = [target[i] for i in keep]
    value  = [value[i]  for i in keep]

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            label=labels, pad=18, thickness=22,
            color=["#440154", "#3b528b", "#999", "#999",
                   "#21918c", "#5ec962",
                   "#1f9e89", "#fde725", "#1f9e89", "#fde725"],
        ),
        link=dict(source=source, target=target, value=value,
                  color="rgba(80,80,80,0.35)"),
    ))
    fig.update_layout(
        title="Figure 6 — Classification flow (Sankey)" + DEF_SUFFIX,
        title_font_size=14,
        font=dict(size=11),
        width=1100, height=550,
        margin=dict(l=10, r=10, t=50, b=10),
    )
    try:
        fig.write_image(out_dir / "fig6_sankey.png", scale=2)
        fig.write_image(out_dir / "fig6_sankey.pdf")
    except Exception as e:
        print(f"  (static export of sankey skipped: {e}) — HTML still written")
    fig.write_html(out_dir / "fig6_sankey.html")
    print("  saved: fig6_sankey.{png,pdf,html}")


# ── Definition comparison — overlap rule #1 vs #2 ───────────────────────
def fig_def_comparison(df: pd.DataFrame, ok: pd.DataFrame, out_dir: Path) -> None:
    """Compare Definition #1 (partial overlap, "3-rule") against Definition #2
    (containment, "2-rule", directional). Containment is strictly stronger, so
    every gene/transcript flagged under #2 is also flagged under #1 — the figure
    quantifies how much the redundancy signal shrinks under the stricter rule.
    """
    col2 = "has_non_unique_transcript_containment"
    if col2 not in ok.columns:
        print("  def-comparison skipped — CSV has no definition #2 columns (re-run pipeline)")
        return

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.7), constrained_layout=True)
    fig.suptitle("Overlap rule comparison: Definition #1 (partial) vs #2 (containment)",
                 fontsize=12, fontweight="bold")
    pal = sns.color_palette("viridis", 5)
    c1, c2 = pal[1], pal[3]

    # ─ Panel A — gene-level flag counts, all OK vs multi-NM eligible set ─
    multi = ok[ok["nm_count"] >= 2]
    groups = [
        ("All OK genes", ok),
        ("Multi-NM (nm≥2)", multi),
    ]
    x = np.arange(len(groups))
    w = 0.38
    d1 = [int(g["has_non_unique_transcript"].sum()) for _, g in groups]
    d2 = [int(g[col2].sum()) for _, g in groups]
    bars1 = axA.bar(x - w / 2, d1, w, color=c1, edgecolor="white", label="Def #1 — partial overlap")
    bars2 = axA.bar(x + w / 2, d2, w, color=c2, edgecolor="white", label="Def #2 — containment")
    for bars, vals, base in ((bars1, d1, groups), (bars2, d2, groups)):
        for bar, v, (_, g) in zip(bars, vals, base):
            pct = 100 * v / max(len(g), 1)
            axA.annotate(f"{v:,}\n{pct:.1f}%",
                         xy=(bar.get_x() + bar.get_width() / 2, v),
                         xytext=(0, 3), textcoords="offset points",
                         ha="center", va="bottom", fontsize=8, color="dimgray")
    axA.set_xticks(x)
    axA.set_xticklabels([f"{name}\n(n={len(g):,})" for name, g in groups])
    axA.set_ylabel("Genes flagged has_non_unique = True")
    axA.set_title("A. Gene-level flag count by definition")
    axA.legend(frameon=False, loc="upper right")
    axA.set_ylim(0, max(d1 + [1]) * 1.25)

    # ─ Panel B — where def-#1 genes land under def-#2 (subset breakdown) ─
    flagged1 = ok[ok["has_non_unique_transcript"]]
    both = int(flagged1[col2].sum())
    only1 = int((~flagged1[col2]).sum())
    neither = int((~ok["has_non_unique_transcript"]).sum())
    # transcript-level totals
    t1 = int(ok["non_unique_transcript_count"].sum())
    t2 = int(ok["non_unique_transcript_count_containment"].sum())

    wedges = [
        (f"Flagged by BOTH\n#1 and #2\n{both:,}", both, c2),
        (f"Flagged by #1 ONLY\n(lost under #2)\n{only1:,}", only1, pal[4]),
        (f"Flagged by NEITHER\n{neither:,}", neither, "lightgray"),
    ]
    wedges = [w_ for w_ in wedges if w_[1] > 0]
    axB.pie([v for _, v, _ in wedges], labels=[lbl for lbl, _, _ in wedges],
            colors=[c for _, _, c in wedges], autopct="",
            startangle=90, counterclock=False,
            wedgeprops=dict(edgecolor="white", linewidth=1.2),
            textprops=dict(fontsize=8))
    axB.set_title("B. Gene fate under the stricter rule\n"
                  f"(non-unique transcripts: #1={t1:,} → #2={t2:,}, "
                  f"−{100*(t1-t2)/max(t1,1):.0f}%)")

    save(fig, "fig_def_comparison", out_dir)


# ── main ────────────────────────────────────────────────────────────────
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", "-i", type=Path,
                    default=PROJECT_ROOT / "output" / "all_protein_coding_v2.csv",
                    help="Input v2 CSV (default: output/all_protein_coding_v2.csv)")
    ap.add_argument("--output-dir", "-o", type=Path,
                    default=PROJECT_ROOT / "output" / "figures",
                    help="Output directory (default: output/figures)")
    ap.add_argument("--definition", "-d", type=int, choices=(1, 2), default=1,
                    help="Overlap rule driving Figures 1–6: 1 = partial overlap "
                         "(default), 2 = containment. The def-comparison figure "
                         "always shows both regardless of this flag.")
    args = ap.parse_args()

    if not args.input.exists():
        print(f"error: input CSV not found: {args.input}", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    apply_style()
    set_definition(args.definition)

    print(f"loading: {args.input}")
    df_raw = load_df(args.input)
    # Figures 1–6 operate on the active-definition view; the def-comparison
    # figure stays on the pristine frame (it needs both definitions side by side).
    df = remap_definition(df_raw)
    ok = df[df["status"] == "OK"].copy()
    ok_raw = df_raw[df_raw["status"] == "OK"].copy()
    rule = "partial overlap" if args.definition == 1 else "containment"
    print(f"  total rows: {len(df):,}  OK rows: {len(ok):,}")
    print(f"  definition #{args.definition} ({rule}) drives Figures 1–6")
    print(f"writing figures to: {args.output_dir}\n")

    print("Figure 1 — composition");                fig1(df, ok, args.output_dir)
    print("Figure 2 — non-uniqueness");             fig2(df, ok, args.output_dir)
    print("Figure 3 — curated vs predicted");       fig3(df, ok, args.output_dir)
    print("Figure 4 — distributions");              fig4(df, ok, args.output_dir)
    print("Figure 5 — exon tracks");                fig5(df, ok, args.output_dir)
    print("Figure 6 — classification flow");        fig6_matplotlib(df, ok, args.output_dir)
    print("Figure 6 — Sankey (Plotly)");            fig6_sankey(df, ok, args.output_dir)
    print("Definition comparison — #1 vs #2");      fig_def_comparison(df_raw, ok_raw, args.output_dir)
    print(f"\nDone. Output: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
