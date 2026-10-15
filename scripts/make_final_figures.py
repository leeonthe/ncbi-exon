"""Generate the FINAL paper figure set (3 figures) per the revision spec in
`PCG Transcripts Computational Analysis Research - Jun&Yoon.pdf` (see README
§ "Final paper figure set"). All figures use Definition #2 (containment).

    Figure 1 — Transcript landscape (A composition bars, B NM distribution,
               C fraction structurally non-unique by NM count)
    Figure 2 — GAPDH exon tracks (structurally unique vs non-unique)
    Figure 3 — Classification flowchart (0 NM / mono-isoform / multi-NM,
               multi-NM split by structural non-uniqueness)

Usage:
    .venv/bin/python scripts/make_final_figures.py \\
        --input output/all_protein_coding_v2.csv \\
        --output-dir output/figures_final
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

# Palette — calm, print-safe. Blue = unique/neutral, red = structurally
# non-unique, amber = mono-isoform highlight, gray = exempt/predicted-only.
C_BLUE = "#4878b0"
C_LIGHTBLUE = "#7fa8d9"
C_RED = "#c94040"
C_AMBER = "#e8a33d"
C_GRAY = "#9aa5b1"
C_DARK = "#2f3640"

SUBTITLE_KW = dict(fontsize=10, color="#555", style="italic")


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


def load_ok(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path, low_memory=False)
    for c in ("is_exon", "has_non_unique_transcript", "has_non_unique_transcript_containment"):
        if c in df.columns and df[c].dtype != bool:
            df[c] = df[c].map({"True": True, "False": False, True: True, False: False}).fillna(False).astype(bool)
    return df


def save(fig: plt.Figure, name: str, out_dir: Path) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"{name}.{ext}")
    plt.close(fig)
    print(f"  saved: {name}.{{png,pdf}}")


# ── Figure 1 — transcript landscape (A + B + C) ─────────────────────────
def fig1(ok: pd.DataFrame, out_dir: Path) -> None:
    fig = plt.figure(figsize=(15.5, 5.9), constrained_layout=True)
    # Reserve the top band for the title + subtitle so they never collide
    # with the panel titles.
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.89))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.15, 1.15])
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))

    fig.suptitle("Figure 1 — The curated and predicted transcript landscape of human protein-coding genes",
                 fontsize=13, fontweight="bold", y=0.985)
    fig.text(0.5, 0.925,
             f"RefSeq transcripts of the {len(ok):,} successfully processed protein-coding genes; "
             "structural non-uniqueness assessed on curated NM transcripts by exon-wise containment (GRCh38)",
             ha="center", **SUBTITLE_KW)

    # ─ Panel A — horizontal composition bars (per the example mockup) ─
    cats = [
        ("NM — curated mRNA",    int(ok["nm_count"].sum())),
        ("NR — curated ncRNA",   int(ok["nr_count"].sum())),
        ("XM — predicted mRNA",  int(ok["xm_count"].sum())),
        ("XR — predicted ncRNA", int(ok["xr_count"].sum())),
    ]
    ypos = np.arange(len(cats))[::-1]
    vals = [v for _, v in cats]
    axA.barh(ypos, vals, height=0.62, color=C_LIGHTBLUE, edgecolor="none", zorder=3)
    for y, v in zip(ypos, vals):
        axA.text(v + max(vals) * 0.015, y, f"{v:,}", va="center", ha="left",
                 fontsize=9, color=C_DARK)
    axA.set_yticks(ypos)
    axA.set_yticklabels([lbl for lbl, _ in cats])
    axA.set_xlim(0, max(vals) * 1.16)
    axA.set_xlabel("Total transcripts (all processed genes pooled)")
    axA.set_title("A. RefSeq transcript composition", loc="left")
    axA.xaxis.grid(True, ls=":", color="#ccc", zorder=0)
    axA.spines["left"].set_visible(False)
    axA.tick_params(axis="y", length=0)

    # ─ Panel B — NM distribution, integer ticks, 0/1 bars called out ─
    CAP = 20
    nm = ok["nm_count"].clip(upper=CAP)
    counts = nm.value_counts().reindex(range(CAP + 1), fill_value=0)
    colors = [C_GRAY if k == 0 else C_AMBER if k == 1 else C_BLUE for k in counts.index]
    axB.bar(counts.index, counts.values, color=colors, edgecolor="white",
            linewidth=0.6, zorder=3)
    axB.set_yscale("log")
    axB.set_xticks(range(0, CAP + 1))
    axB.set_xticklabels([str(k) if k < CAP else f"{CAP}+" for k in range(CAP + 1)],
                        fontsize=7.5)
    axB.set_xlim(-0.7, CAP + 0.7)
    axB.set_xlabel("NM transcripts per gene")
    axB.set_ylabel("Number of genes (log scale)")
    axB.set_title(f"B. NM transcripts per gene (n = {len(ok):,})", loc="left")

    n0 = int(counts[0]); n1 = int(counts[1])
    axB.set_ylim(top=counts.max() * 4)
    axB.text(0, n0 * 1.25, f"{n0:,}\n({100*n0/len(ok):.1f}%)", ha="center",
             va="bottom", fontsize=7.5, color="#5a6470")
    axB.text(1, n1 * 1.25, f"{n1:,}\n({100*n1/len(ok):.1f}%)", ha="center",
             va="bottom", fontsize=7.5, color="#a86f1c", fontweight="bold")
    axB.legend(handles=[
        Patch(facecolor=C_GRAY, label="No curated NM (XM-predicted only)"),
        Patch(facecolor=C_AMBER, label="Mono-isoform (1 NM)"),
        Patch(facecolor=C_BLUE, label="Multi-NM (≥2 NM)"),
    ], loc="upper right", frameon=False, fontsize=8)

    # ─ Panel C — fraction structurally non-unique by NM count (old Fig 2A) ─
    bins = [-1, 0, 1, 2, 3, 5, 10, 20, 50, 1000]
    labels = ["0", "1", "2", "3", "4–5", "6–10", "11–20", "21–50", "50+"]
    tmp = ok.copy()
    tmp["nm_bucket"] = pd.cut(tmp["nm_count"], bins=bins, labels=labels, right=True)
    agg = tmp.groupby("nm_bucket", observed=True).agg(
        n=("input_symbol", "size"),
        n_nu=("has_non_unique_transcript_containment", "sum"),
    )
    agg["frac"] = agg["n_nu"] / agg["n"]
    x = np.arange(len(agg))
    bar_colors = [C_GRAY if lbl == "0" else C_AMBER if lbl == "1" else C_RED
                  for lbl in agg.index.astype(str)]
    bars = axC.bar(x, agg["frac"], color=bar_colors, edgecolor="white",
                   linewidth=0.6, zorder=3)
    axC.set_xticks(x)
    axC.set_xticklabels(agg.index.astype(str), fontsize=8.5)
    axC.set_xlabel("NM transcripts per gene (binned)")
    axC.set_ylabel("Fraction with structurally non-unique\nNM transcript(s)")
    axC.set_title("C. Structural non-uniqueness by NM count", loc="left")
    axC.set_ylim(0, 1.14)
    axC.yaxis.grid(True, ls=":", color="#ddd", zorder=0)
    # n labels: above each bar; for the two zero-height exempt bins put the
    # label at the baseline so nothing overlaps.
    for xi, (bar, (lbl, row)) in enumerate(zip(bars, agg.iterrows())):
        h = bar.get_height()
        axC.annotate(f"n={int(row['n']):,}",
                     xy=(xi, h), xytext=(0, 4), textcoords="offset points",
                     ha="center", va="bottom", fontsize=7.5, color="dimgray")
    axC.text(0.5, 0.075, "0–1 NM: exempt —\nno sibling NM to compare",
             transform=axC.transData, ha="center", va="bottom",
             fontsize=8, color="#5a6470",
             bbox=dict(boxstyle="round,pad=0.3", fc="#f2f4f7", ec="none"))

    save(fig, "figure1_transcript_landscape", out_dir)


# ── Figure 2 — GAPDH exon tracks ────────────────────────────────────────
def fig2(ok: pd.DataFrame, out_dir: Path) -> None:
    from ncbi_exon.parser import _grch38_exons, parse_product_report

    cache_file = PROJECT_ROOT / "cache" / "GAPDH__human.json"
    if not cache_file.exists():
        print("  figure 2 skipped — GAPDH cache missing")
        return
    payload = json.loads(cache_file.read_text())
    prod = payload["reports"][0]["product"]
    row = parse_product_report("GAPDH", payload)
    nu_accs = row.non_unique_transcript_accessions_containment
    nu = set(nu_accs.split(";")) if nu_accs else set()

    rows_data = []
    for t in prod.get("transcripts", []):
        acc = t.get("accession_version", "")
        if not acc.startswith("NM_"):
            continue
        ex = _grch38_exons(t)
        if ex:
            rows_data.append((acc, ex))

    fig, ax = plt.subplots(figsize=(11.5, 4.4), constrained_layout=True)
    # Top band for title + subtitle, bottom band for the legend.
    fig.get_layout_engine().set(rect=(0, 0.09, 1, 0.77))
    fig.suptitle("Figure 2 — Exon-level transcript structure of GAPDH (GRCh38)",
                 fontsize=13, fontweight="bold", y=0.975)
    fig.text(0.5, 0.895,
             "A structurally non-unique transcript has every exon fully contained "
             "within an exon of another curated (NM) transcript of the same gene",
             ha="center", **SUBTITLE_KW)

    starts = [e[0] for _, exs in rows_data for e in exs]
    ends = [e[1] for _, exs in rows_data for e in exs]
    x0, x1 = min(starts), max(ends)
    span = x1 - x0
    ax.set_xlim(x0 - 0.03 * span, x1 + 0.02 * span)
    ax.set_ylim(-0.6, len(rows_data) - 0.4)

    for i, (acc, exs) in enumerate(rows_data):
        y = len(rows_data) - 1 - i
        color = C_RED if acc in nu else C_BLUE
        ax.hlines(y, min(e[0] for e in exs), max(e[1] for e in exs),
                  color="#c9ced6", lw=1.0, zorder=1)
        for b, e in exs:
            ax.add_patch(plt.Rectangle((b, y - 0.28), e - b, 0.56,
                                       facecolor=color, edgecolor="none", zorder=3))
        ax.text(x0 - 0.035 * span, y, acc, ha="right", va="center",
                fontsize=9, color=color,
                fontweight="bold" if acc in nu else "normal")

    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    chrom = prod["transcripts"][0]["genomic_locations"][0]["genomic_accession_version"]
    ax.set_xlabel(f"GRCh38 genomic position on chromosome 12 ({chrom})")
    ax.set_title(f"GAPDH (gene_id {prod['gene_id']}) — {len(rows_data)} curated NM transcripts, "
                 f"{len(nu)} structurally non-unique", loc="left", fontsize=10)

    fig.legend(handles=[
        Patch(facecolor=C_BLUE, label="Structurally unique transcript"),
        Patch(facecolor=C_RED, label="Structurally non-unique transcript"),
    ], loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0))

    save(fig, "figure2_gapdh_exon_tracks", out_dir)


# ── Figure 3 — classification flowchart ─────────────────────────────────
def fig3(df: pd.DataFrame, ok: pd.DataFrame, out_dir: Path) -> None:
    total = len(df)
    n_ok = len(ok)
    n_drop = total - n_ok
    n0 = int((ok["nm_count"] == 0).sum())
    n1 = int((ok["nm_count"] == 1).sum())
    multi = ok[ok["nm_count"] >= 2]
    nm_ = len(multi)
    n_nu = int(multi["has_non_unique_transcript_containment"].sum())
    n_u = nm_ - n_nu

    fig, ax = plt.subplots(figsize=(12, 8.2))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.axis("off")
    fig.suptitle("Figure 3 — Classification of human protein-coding genes by curated transcript structure",
                 fontsize=13, fontweight="bold", y=0.97)
    fig.text(0.5, 0.925,
             "Genes are partitioned by curated (NM) transcript count; multi-NM genes are then screened for "
             "structurally non-unique transcripts (exon-wise containment, GRCh38)",
             ha="center", **SUBTITLE_KW)

    def box(cx, cy, w, h, fc, title, lines, n_title_lines=1, title_color="white",
            text_color="white", title_fs=10.5, text_fs=9):
        ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                    boxstyle="round,pad=0.6,rounding_size=1.2",
                                    facecolor=fc, edgecolor="none", zorder=3))
        ax.text(cx, cy + h / 2 - 1.8, title, ha="center", va="top",
                fontsize=title_fs, fontweight="bold", color=title_color, zorder=4)
        # Body text is anchored BELOW the title band so the two never collide.
        y_body = cy + h / 2 - 1.8 - 3.1 * n_title_lines - 0.8
        ax.text(cx, y_body, "\n".join(lines), ha="center", va="top",
                fontsize=text_fs, color=text_color, zorder=4, linespacing=1.4)

    def arrow(x0, y0, x1, y1, label=None, lx_off=1.5, color="#8a94a1"):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1),
                                     arrowstyle="-|>", mutation_scale=14,
                                     lw=1.4, color=color, zorder=2))
        if label:
            ax.text((x0 + x1) / 2 + lx_off, (y0 + y1) / 2, label,
                    ha="left", va="center", fontsize=8.5, color="#5a6470")

    # Stage 1 — input
    box(50, 89.5, 42, 9, C_DARK, "Human protein-coding genes",
        [f"n = {total:,}"])
    arrow(50, 84.5, 50, 77)
    ax.text(52, 80.8, "NCBI Datasets retrieval & filtering\n"
            f"({n_drop} genes excluded: 9 not found — mitochondrial; 1 ambiguous symbol match)",
            ha="left", va="center", fontsize=8.3, color="#5a6470", linespacing=1.4)

    # Stage 2 — processed
    box(50, 71.5, 42, 9, "#4a5b74", "Successfully processed",
        [f"n = {n_ok:,}  ({100*n_ok/total:.2f}%)"])

    # Stage 3 — split by NM count
    y3 = 47
    h3 = 15
    xs = {"zero": 16.5, "mono": 45.5, "multi": 79}
    arrow(50, 66.5, xs["zero"], y3 + h3 / 2 + 0.5)
    arrow(50, 66.5, xs["mono"], y3 + h3 / 2 + 0.5)
    arrow(50, 66.5, xs["multi"], y3 + h3 / 2 + 0.5)
    ax.text(50, 60.5, "Curated (NM) transcript count", ha="center",
            fontsize=9, color="#5a6470", style="italic",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none"), zorder=5)

    box(xs["zero"], y3, 24, h3, C_GRAY, "No curated NM",
        [f"n = {n0:,}  ({100*n0/n_ok:.1f}%)", "XM-predicted transcripts only", "— exempt from comparison —"],
        text_fs=8.3)
    box(xs["mono"], y3, 26, h3, C_AMBER, "Mono-isoform",
        [f"n = {n1:,}  ({100*n1/n_ok:.1f}%)", "single NM transcript", "— no sibling to compare —"],
        text_fs=8.3)
    box(xs["multi"], y3, 30, h3, C_BLUE, "Multi-NM",
        [f"n = {nm_:,}  ({100*nm_/n_ok:.1f}%)", "≥ 2 NM transcripts", "→ pairwise exon comparison"],
        text_fs=8.3)

    # Stage 4 — multi-NM split
    y4 = 16.5
    h4 = 17.5
    arrow(xs["multi"], y3 - h3 / 2 - 0.5, 61, y4 + h4 / 2 + 0.5)
    arrow(xs["multi"], y3 - h3 / 2 - 0.5, 91, y4 + h4 / 2 + 0.5)
    ax.text(76, 32, "Exon-wise containment analysis\n(Definition #2)", ha="center",
            fontsize=9, color="#5a6470", style="italic",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none"), zorder=5)

    box(61, y4, 29, h4, C_RED, "≥ 1 structurally\nnon-unique transcript",
        [f"n = {n_nu:,}", f"{100*n_nu/nm_:.1f}% of multi-NM", f"{100*n_nu/n_ok:.1f}% of all processed"],
        n_title_lines=2, title_fs=10, text_fs=8.6)
    box(91, y4, 25.5, h4, "#6b8f6e", "All transcripts\nstructurally unique",
        [f"n = {n_u:,}", f"{100*n_u/nm_:.1f}% of multi-NM", f"{100*n_u/n_ok:.1f}% of all processed"],
        n_title_lines=2, title_fs=10, text_fs=8.6)

    # Proportional strip along the bottom — the quantitative summary
    y_strip, h_strip = 3.2, 4.2
    segs = [
        (n0, C_GRAY, "No NM"),
        (n1, C_AMBER, "Mono-isoform"),
        (n_nu, C_RED, "Multi-NM, non-unique"),
        (n_u, "#6b8f6e", "Multi-NM, all unique"),
    ]
    left = 8
    width_total = 84
    for n, c, lbl in segs:
        w = width_total * n / n_ok
        ax.add_patch(plt.Rectangle((left, y_strip), w, h_strip, facecolor=c,
                                   edgecolor="white", linewidth=1.2, zorder=3))
        pct = 100 * n / n_ok
        if w > 7:
            ax.text(left + w / 2, y_strip + h_strip / 2, f"{pct:.1f}%",
                    ha="center", va="center", fontsize=8.3, color="white",
                    fontweight="bold", zorder=4)
        else:
            ax.text(left + w / 2, y_strip - 0.6, f"{pct:.1f}%",
                    ha="center", va="top", fontsize=7.5, color="#5a6470", zorder=4)
        left += w
    ax.text(8, y_strip + h_strip + 1.6, f"All processed genes (n = {n_ok:,}), to scale:",
            ha="left", va="bottom", fontsize=8.5, color="#5a6470")

    save(fig, "figure3_classification_flow", out_dir)


# ── main ────────────────────────────────────────────────────────────────
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", "-i", type=Path,
                    default=PROJECT_ROOT / "output" / "all_protein_coding_v2.csv")
    ap.add_argument("--output-dir", "-o", type=Path,
                    default=PROJECT_ROOT / "output" / "figures_final")
    args = ap.parse_args()

    if not args.input.exists():
        print(f"error: input CSV not found: {args.input}", file=sys.stderr)
        return 1
    args.output_dir.mkdir(parents=True, exist_ok=True)
    apply_style()

    df = load_ok(args.input)
    ok = df[df["status"] == "OK"].copy()
    print(f"loaded {len(df):,} rows ({len(ok):,} OK) — Definition #2 (containment)")
    print(f"writing figures to: {args.output_dir}\n")

    print("Figure 1 — transcript landscape"); fig1(ok, args.output_dir)
    print("Figure 2 — GAPDH exon tracks");    fig2(ok, args.output_dir)
    print("Figure 3 — classification flow");  fig3(df, ok, args.output_dir)
    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
