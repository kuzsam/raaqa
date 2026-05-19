import argparse
import math
import os
import re
import sys

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

try:
    from importlib.metadata import version
    VERSION = version("raaqa")
except Exception:
    VERSION = "unknown"

# ============================================================
# Shared styling
# ============================================================

BG = "#ffffff"       # figure background
AX_BG = "#f6f8fa"    # axes background
GRID = "#d0d7de"     # grid lines
TEXT = "#1f2328"     # titles and primary labels
TEXT_DIM = "#57606a" # axis labels and tick labels
REF_LINE = "#cf222e" # mean reference axhlines
MED_LINE = "#8250df" # median reference axhlines
FONT_MAIN = "DejaVu Sans"

MEAN_STYLE = dict(color=REF_LINE, linewidth=1.4, linestyle="-", zorder=5)
MEDIAN_STYLE = dict(color=MED_LINE, linewidth=1.4, linestyle="-", zorder=5)
MEAN_STYLE_DASH = dict(color=REF_LINE, linewidth=1, linestyle="--", zorder=5)
MEDIAN_STYLE_DASH = dict(color=MED_LINE, linewidth=1, linestyle="--", zorder=5)


# ============================================================
# Shared utilities
# ============================================================

def _check_columns(df, required, filename):
    """Exit with an error if any required columns are missing."""
    missing = required - set(df.columns)
    if missing:
        sys.exit(
            f"[ERROR] {filename} is missing required columns: {', '.join(sorted(missing))}"
        )


def apply_style(fig, axes):
    """Apply shared figure and axes styling: background, grid, tick and label colours."""
    fig.patch.set_facecolor(BG)
    for ax in axes if hasattr(axes, "__iter__") else [axes]:
        ax.set_facecolor(AX_BG)
        ax.tick_params(colors=TEXT_DIM, labelsize=8)
        ax.xaxis.label.set_color(TEXT_DIM)
        ax.yaxis.label.set_color(TEXT_DIM)
        ax.title.set_color(TEXT)
        for spine in ax.spines.values():
            spine.set_edgecolor(GRID)
        ax.grid(True, color=GRID, linewidth=0.5, alpha=1.0)
        ax.set_axisbelow(True)


def smart_bp_formatter(contig_length_bp):
    """Return (divisor, unit) for scaling a base-pair value to bp, kb, or Mb."""
    if contig_length_bp >= 1_000_000:
        return 1_000_000, "Mb"
    elif contig_length_bp >= 1_000:
        return 1_000, "kb"
    else:
        return 1, "bp"


def _sanitise_filename(name):
    """Replace characters not safe in filenames with underscores."""
    return re.sub(r"[^\w.\-]", "_", name)


def _safe_nanmax(values, default):
    """Return the max of finite values in values, or default if none are finite."""
    arr = np.asarray(values, dtype=float)
    finite = arr[np.isfinite(arr)]
    return float(finite.max()) if len(finite) > 0 else default


def make_figures_folder(input_folder):
    """Create and return the figures/ subfolder inside the input directory."""
    figures_folder = os.path.join(input_folder, "figures")
    os.makedirs(figures_folder, exist_ok=True)
    return figures_folder


# ============================================================
# Argument parsing and input validation
# ============================================================

def parse_arguments():
    """Parse and return command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate QC figures from mapq_softclip or hese module outputs."
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"RAAQA: {VERSION}"
    )
    parser.add_argument(
        "-m",
        "--module",
        required=True,
        choices=["mapq_softclip", "hese"],
        help="Which module's output to visualise",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="DIR",
        help="Path to the output folder from the respective module run",
    )
    parser.add_argument(
        "-r",
        "--rolling",
        type=float,
        default=1000,
        metavar="FLOAT",
        help="Genomic distance in kb to span for rolling mean in rolling figures, mapq_softclip module only (default: 1000)",
    )
    return parser.parse_args()


def check_input(input_folder, module):
    """Validate that the input folder and required output files exist for the given module."""
    if not os.path.isdir(input_folder):
        sys.exit(f"[ERROR] Input folder not found: {input_folder}")
    if module == "mapq_softclip":
        for f in ["window_stats.csv", "summary_stats.csv"]:
            if not os.path.isfile(os.path.join(input_folder, f)):
                sys.exit(f"[ERROR] Required file not found: {f}")
    elif module == "hese":
        for f in ["haplotig_labels.csv", "haplotype_summary.csv", "unitig_labels.csv"]:
            if not os.path.isfile(os.path.join(input_folder, f)):
                sys.exit(f"[ERROR] Required file not found: {f}")


# ============================================================
# mapq_softclip: styling
# ============================================================

COV_A = "#0969da"    # coverage fill, mean coverage dots/lines in genome-wide
COV_B = "#54aeff"    # median coverage dots/lines in genome-wide
COV_LINE = "#0550ae" # coverage overflow triangle markers
MAPQ_A = "#79a80a"   # median MAPQ fill, median MAPQ dots/lines in genome-wide
MAPQ_B = "#0e9488"   # mean MAPQ fill, mean MAPQ dots/lines in genome-wide
SC_A = "#e36209"     # softclip fill, softclip dots/lines in genome-wide
SC_LINE = "#953800"  # softclip reference axhlines and softclip overflow markers

SOFTCLIP_STYLE = dict(color=SC_LINE, linewidth=1, linestyle="-", zorder=5)

REQUIRED_WIN_COLS = {
    "Chromosome",
    "Start",
    "End",
    "Total_Bases",
    "Median_MAPQ",
    "Mean_MAPQ",
    "Softclip_%"
}
REQUIRED_SUM_COLS = {"Chromosome", "Mean_MAPQ", "Median_MAPQ", "Softclip_%"}


# ============================================================
# mapq_softclip: helpers
# ============================================================

def _build_genome_xaxis(df, chroms):
    """Assign each contig an integer x-axis index for equal-spaced genome-wide plots."""
    dividers = []
    chrom_ticks = []
    chrom_labels = []
    chrom_x = {}
    idx = 0

    for chrom in chroms:
        cdf = df[df["Chromosome"] == chrom]
        if cdf.empty:
            continue
        chrom_x[chrom] = (cdf, idx)
        chrom_ticks.append(idx)
        chrom_labels.append(chrom)
        dividers.append(idx + 0.5)
        idx += 1

    return idx, dividers, chrom_ticks, chrom_labels, chrom_x


def _apply_genome_xticks(ax, chrom_ticks, chrom_labels):
    """Set contig name tick labels on ax, rotated 90 degrees for readability."""
    ax.set_xticks(chrom_ticks)
    ax.set_xticklabels(
        chrom_labels, rotation=90, fontsize=4.5, color=TEXT_DIM, fontfamily=FONT_MAIN
    )


def _draw_dividers(axes, dividers):
    """Draw faint vertical lines between contigs across all axes."""
    for d in dividers[:-1]:
        for ax in axes:
            ax.axvline(d, color=GRID, linewidth=0.5, alpha=0.6)


# ============================================================
# mapq_softclip: figures
# ============================================================

def fig_genome_three_panel(df_win, df_sum, figures_folder):
    """Save genome_summary.png with one data point per contig across coverage, MAPQ and softclip panels."""
    df = df_win[df_win["Chromosome"] != "GENOME"].copy()
    df["Coverage"] = df["Total_Bases"] / (df["End"] - df["Start"]).replace(0, np.nan)
    df_sum_ctg = df_sum[df_sum["Chromosome"] != "GENOME"].copy()

    # genome-wide reference values, pre-computed in the GENOME row of summary_stats.csv
    genome_row = df_sum[df_sum["Chromosome"] == "GENOME"]
    if genome_row.empty:
        print(
            "[WARN] No GENOME row found in summary_stats.csv -- genome-wide reference lines will be omitted."
        )
    genome_mean_mapq = (
        float(genome_row["Mean_MAPQ"].iloc[0]) if not genome_row.empty else np.nan
    )
    genome_med_mapq = (
        float(genome_row["Median_MAPQ"].iloc[0]) if not genome_row.empty else np.nan
    )
    genome_softclip = (
        float(genome_row["Softclip_%"].iloc[0]) if not genome_row.empty else np.nan
    )

    chroms = list(dict.fromkeys(df["Chromosome"].tolist()))

    fig, (ax_cov, ax_mapq, ax_sc) = plt.subplots(
        3, 1, figsize=(24, 10), sharex=True, constrained_layout=True
    )
    apply_style(fig, [ax_cov, ax_mapq, ax_sc])

    genome_pos, dividers, chrom_ticks, chrom_labels, chrom_x = _build_genome_xaxis(
        df, chroms
    )

    contig_mids = []
    contig_cov_means = []
    contig_cov_meds = []
    contig_means = []
    contig_meds = []
    contig_softclips = []

    for chrom in chroms:
        if chrom not in chrom_x:
            continue
        cdf, genome_offset = chrom_x[chrom]

        contig_mids.append(genome_offset)

        # coverage derived from window data, not read directly from CSV
        contig_cov_means.append(float(cdf["Coverage"].mean()))
        contig_cov_meds.append(float(cdf["Coverage"].median()))

        # MAPQ and softclip, pre-computed in summary_stats.csv
        row = df_sum_ctg[df_sum_ctg["Chromosome"] == chrom]
        contig_means.append(
            float(row["Mean_MAPQ"].iloc[0]) if not row.empty else np.nan
        )
        contig_meds.append(
            float(row["Median_MAPQ"].iloc[0]) if not row.empty else np.nan
        )
        contig_softclips.append(
            float(row["Softclip_%"].iloc[0]) if not row.empty else np.nan
        )

    _draw_dividers([ax_cov, ax_mapq, ax_sc], dividers)

    # genome-wide coverage reference lines
    mean_cov = df["Coverage"].mean()
    median_cov = df["Coverage"].median()
    ax_cov.axhline(
        mean_cov, **MEAN_STYLE_DASH, label=f"Genome mean coverage  {mean_cov:.1f}×"
    )
    ax_cov.axhline(
        median_cov,
        **MEDIAN_STYLE_DASH,
        label=f"Genome median coverage  {median_cov:.1f}×",
    )

    ax_cov.plot(
        contig_mids,
        contig_cov_means,
        color=COV_A,
        linewidth=1.2,
        linestyle="-",
        zorder=4,
        label="Mean coverage per contig",
    )
    ax_cov.plot(
        contig_mids,
        contig_cov_meds,
        color=COV_B,
        linewidth=1.2,
        linestyle="-",
        zorder=4,
        label="Median coverage per contig",
    )
    ax_cov.scatter(contig_mids, contig_cov_means, color=COV_A, s=12, zorder=5)
    ax_cov.scatter(contig_mids, contig_cov_meds, color=COV_B, s=12, zorder=5)

    ax_mapq.plot(
        contig_mids,
        contig_means,
        color=MAPQ_B,
        linewidth=1.2,
        linestyle="-",
        zorder=4,
        label="Mean MAPQ per contig",
    )
    ax_mapq.plot(
        contig_mids,
        contig_meds,
        color=MAPQ_A,
        linewidth=1.2,
        linestyle="-",
        zorder=4,
        label="Median MAPQ per contig",
    )
    ax_mapq.scatter(contig_mids, contig_means, color=MAPQ_B, s=12, zorder=5)
    ax_mapq.scatter(contig_mids, contig_meds, color=MAPQ_A, s=12, zorder=5)
    ax_mapq.axhline(
        genome_mean_mapq,
        **MEAN_STYLE_DASH,
        label=f"Genome mean MAPQ  {genome_mean_mapq:.2f}",
    )
    ax_mapq.axhline(
        genome_med_mapq,
        **MEDIAN_STYLE_DASH,
        label=f"Genome median MAPQ  {genome_med_mapq:.1f}",
    )

    ax_sc.plot(
        contig_mids,
        contig_softclips,
        color=SC_A,
        linewidth=1.2,
        linestyle="-",
        zorder=4,
        label="Softclip % per contig",
    )
    ax_sc.scatter(contig_mids, contig_softclips, color=SC_A, s=12, zorder=5)
    ax_sc.axhline(
        genome_softclip,
        **{**SOFTCLIP_STYLE, "linestyle": "--"},
        label=f"Genome softclip %  {genome_softclip:.2f}",
    )

    cov_vals = [
        v
        for v in contig_cov_means + contig_cov_meds + [mean_cov, median_cov]
        if np.isfinite(v)
    ]
    cov_ylim_top = max(cov_vals) * 1.3 if cov_vals else 1.0

    ax_cov.set_xlim(-0.5, genome_pos - 0.5)
    ax_cov.set_ylim(0, cov_ylim_top)
    ax_cov.set_ylabel("Coverage", color=TEXT_DIM, fontsize=10, fontfamily=FONT_MAIN)
    ax_cov.set_title(
        "Coverage, MAPQ & Softclip Summary — Genome-wide",
        color=TEXT,
        fontsize=13,
        fontweight="bold",
        fontfamily=FONT_MAIN,
    )
    ax_cov.legend(
        fontsize=9,
        facecolor=AX_BG,
        edgecolor=GRID,
        labelcolor=TEXT,
        framealpha=0.9,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        bbox_transform=ax_cov.transAxes,
    )

    mapq_data_max = _safe_nanmax(
        [*contig_means, *contig_meds, genome_mean_mapq, genome_med_mapq], default=np.nan
    )
    ylim_top_mapq = max(65, mapq_data_max * 1.1) if np.isfinite(mapq_data_max) else 65
    ax_mapq.set_ylim(0, ylim_top_mapq)
    ax_mapq.set_ylabel(
        "MAPQ per contig", color=TEXT_DIM, fontsize=10, fontfamily=FONT_MAIN
    )
    ax_mapq.legend(
        fontsize=9,
        facecolor=AX_BG,
        edgecolor=GRID,
        labelcolor=TEXT,
        framealpha=0.9,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        bbox_transform=ax_mapq.transAxes,
    )

    ax_sc.set_ylim(bottom=0)
    ax_sc.set_ylabel(
        "Softclip % per contig", color=TEXT_DIM, fontsize=10, fontfamily=FONT_MAIN
    )
    _apply_genome_xticks(ax_sc, chrom_ticks, chrom_labels)
    ax_sc.legend(
        fontsize=9,
        facecolor=AX_BG,
        edgecolor=GRID,
        labelcolor=TEXT,
        framealpha=0.9,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        bbox_transform=ax_sc.transAxes,
    )

    out = os.path.join(figures_folder, "genome_summary.png")
    fig.savefig(out, dpi=180, facecolor=BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out}")


def fig_per_contig_raw(df_win, df_sum, figures_folder):
    """Save per-contig 4-panel raw figures (coverage, median MAPQ, mean MAPQ, softclip) to figures/contigs/raw/."""
    df = df_win[df_win["Chromosome"] != "GENOME"].copy()
    df["Coverage"] = df["Total_Bases"] / (df["End"] - df["Start"]).replace(0, np.nan)
    df_sum_ctg = df_sum[df_sum["Chromosome"] != "GENOME"].copy()

    if df.empty:
        print(
            "[WARN] No per-contig data found in window_stats.csv -- skipping raw contig figures."
        )
        return

    contigs_folder = os.path.join(figures_folder, "contigs", "raw")
    try:
        os.makedirs(contigs_folder, exist_ok=True)
    except OSError as e:
        sys.exit(f"[ERROR] Cannot create output directory: {e}")

    chroms = list(dict.fromkeys(df["Chromosome"].tolist()))
    total = len(chroms)

    for i, chrom in enumerate(chroms, 1):
        cdf = df[df["Chromosome"] == chrom].copy()
        if cdf.empty:
            continue

        contig_len = int(cdf["End"].max())
        divisor, unit = smart_bp_formatter(contig_len)
        x = cdf["Start"].values / divisor

        win_bp = int(cdf["End"].iloc[0] - cdf["Start"].iloc[0])
        win_div, win_unit = smart_bp_formatter(win_bp)
        win_str = f"{win_bp / win_div:g} {win_unit}"
        if len(cdf) > 1:
            step_bp = int(cdf["Start"].iloc[1] - cdf["Start"].iloc[0])
            step_div, step_unit = smart_bp_formatter(step_bp)
            step_str = f"{step_bp / step_div:g} {step_unit}"
        else:
            step_str = "N/A (single window)"

        # coverage derived from window_stats.csv, not read directly
        cov = cdf["Coverage"].values
        mean_cov = cdf["Coverage"].mean()
        median_cov = cdf["Coverage"].median()

        # MAPQ and softclip per window, pre-computed in window_stats.csv
        med_mapq_win = cdf["Median_MAPQ"].values
        mean_mapq_win = cdf["Mean_MAPQ"].values
        sc_win = cdf["Softclip_%"].values

        # contig-level reference values, pre-computed in summary_stats.csv
        row = df_sum_ctg[df_sum_ctg["Chromosome"] == chrom]
        contig_med_mapq = float(row["Median_MAPQ"].iloc[0]) if not row.empty else np.nan
        contig_mean_mapq = float(row["Mean_MAPQ"].iloc[0]) if not row.empty else np.nan
        contig_softclip = float(row["Softclip_%"].iloc[0]) if not row.empty else np.nan

        fig, (ax_cov, ax_med, ax_mean, ax_sc) = plt.subplots(
            4, 1, figsize=(18, 11), sharex=True, constrained_layout=True
        )
        apply_style(fig, [ax_cov, ax_med, ax_mean, ax_sc])

        xlim = (cdf["Start"].min() / divisor, cdf["End"].max() / divisor)

        # coverage panel
        # Cap y-axis at 99th percentile
        p99_cov = np.nanpercentile(cov, 99)
        ylim_top_cov = (p99_cov * 1.2) if np.isfinite(p99_cov) else 1.0
        ylim_top_cov = max(ylim_top_cov, mean_cov * 2, 1.0)

        ax_cov.fill_between(x, cov, color=COV_A, alpha=0.95, linewidth=0)
        ax_cov.axhline(
            mean_cov, **MEAN_STYLE, label=f"Contig mean cov  {mean_cov:.1f}×"
        )
        ax_cov.axhline(
            median_cov, **MEDIAN_STYLE, label=f"Contig median cov  {median_cov:.1f}×"
        )

        overflow_cov = cov > ylim_top_cov
        if overflow_cov.any():
            cov_clip = cov[overflow_cov]
            ax_cov.scatter(
                x[overflow_cov],
                np.full(overflow_cov.sum(), ylim_top_cov * 0.96),
                marker="^",
                color=COV_LINE,
                s=20,
                zorder=6,
                alpha=0.35,
                label=f"Clipped above {ylim_top_cov:.0f}× (P99)\n{overflow_cov.sum()} windows: {cov_clip.min():.0f}×–{cov_clip.max():.0f}×",
            )

        ax_cov.set_ylabel("Coverage", color=TEXT_DIM, fontsize=9, fontfamily=FONT_MAIN)
        ax_cov.set_ylim(0, ylim_top_cov)
        ax_cov.set_xlim(*xlim)
        ax_cov.legend(
            fontsize=8,
            facecolor=AX_BG,
            edgecolor=GRID,
            labelcolor=TEXT,
            framealpha=0.9,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            bbox_transform=ax_cov.transAxes,
        )
        ax_cov.set_title(
            f"Coverage, MAPQ & Softclip — raw per window — {chrom}\n"
            f"window size: {win_str}  |  step size: {step_str}",
            color=TEXT,
            fontsize=12,
            fontweight="bold",
            fontfamily=FONT_MAIN,
        )

        mapq_data_max = _safe_nanmax(
            [*mean_mapq_win, *med_mapq_win, contig_mean_mapq, contig_med_mapq], default=np.nan
        )
        ylim_top_mapq = (
            max(65, mapq_data_max * 1.1) if np.isfinite(mapq_data_max) else 65
        )

        # median MAPQ panel
        ax_med.fill_between(x, med_mapq_win, color=MAPQ_A, alpha=0.95, linewidth=0)
        ax_med.axhline(
            contig_med_mapq,
            **MEDIAN_STYLE,
            label=f"Contig median MAPQ  {contig_med_mapq:.1f}",
        )
        ax_med.set_ylabel(
            "Median MAPQ", color=TEXT_DIM, fontsize=9, fontfamily=FONT_MAIN
        )
        ax_med.set_ylim(0, ylim_top_mapq)
        ax_med.legend(
            fontsize=8,
            facecolor=AX_BG,
            edgecolor=GRID,
            labelcolor=TEXT,
            framealpha=0.9,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            bbox_transform=ax_med.transAxes,
        )

        # mean MAPQ panel
        ax_mean.fill_between(x, mean_mapq_win, color=MAPQ_B, alpha=0.95, linewidth=0)
        ax_mean.axhline(
            contig_mean_mapq,
            **MEAN_STYLE,
            label=f"Contig mean MAPQ  {contig_mean_mapq:.2f}",
        )
        ax_mean.set_ylabel(
            "Mean MAPQ", color=TEXT_DIM, fontsize=9, fontfamily=FONT_MAIN
        )
        ax_mean.set_ylim(0, ylim_top_mapq)
        ax_mean.legend(
            fontsize=8,
            facecolor=AX_BG,
            edgecolor=GRID,
            labelcolor=TEXT,
            framealpha=0.9,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            bbox_transform=ax_mean.transAxes,
        )

        # softclip panel
        sc_max = _safe_nanmax(sc_win, default=0.1)
        ylim_top_sc = max(sc_max * 1.05, 0.1)

        ax_sc.fill_between(
            x, sc_win, color=SC_A, alpha=0.95, linewidth=0.5, edgecolor=SC_A
        )
        ax_sc.axhline(
            contig_softclip,
            **SOFTCLIP_STYLE,
            label=f"Contig softclip %  {contig_softclip:.2f}",
        )

        ax_sc.set_ylabel("Softclip %", color=TEXT_DIM, fontsize=9, fontfamily=FONT_MAIN)
        ax_sc.set_ylim(0, ylim_top_sc)
        ax_sc.legend(
            fontsize=8,
            facecolor=AX_BG,
            edgecolor=GRID,
            labelcolor=TEXT,
            framealpha=0.9,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            bbox_transform=ax_sc.transAxes,
        )
        ax_sc.set_xlabel(
            f"Position ({unit})", color=TEXT_DIM, fontsize=9, fontfamily=FONT_MAIN
        )
        ax_sc.xaxis.set_major_formatter(
            ticker.FuncFormatter(
                lambda v, _: f"{v:,.1f}" if divisor > 1 else f"{int(v):,}"
            )
        )

        safe = _sanitise_filename(chrom)
        out = os.path.join(contigs_folder, f"{safe}_per_window_raw.png")
        fig.savefig(out, dpi=150, facecolor=BG, bbox_inches="tight")
        plt.close(fig)

        if i % 50 == 0 or i == total:
            print(f"[INFO] Raw per-contig figures: {i}/{total}")

    print(f"[INFO] Raw per-contig figures saved to: {contigs_folder}")


def fig_per_contig_rolling(df_win, df_sum, figures_folder, rolling_target_bp):
    """Save per-contig 3-panel rolling mean figures to figures/contigs/rolling/.

    rolling_target_bp sets the genomic span in base pairs that the rolling window aims to cover.
    The number of windows used is derived from the step size and capped at contig length.
    """
    df = df_win[df_win["Chromosome"] != "GENOME"].copy()
    df["Coverage"] = df["Total_Bases"] / (df["End"] - df["Start"]).replace(0, np.nan)
    df_sum_ctg = df_sum[df_sum["Chromosome"] != "GENOME"].copy()

    if df.empty:
        print(
            "[WARN] No per-contig data found in window_stats.csv -- skipping rolling contig figures."
        )
        return

    rolling_folder = os.path.join(figures_folder, "contigs", "rolling")
    try:
        os.makedirs(rolling_folder, exist_ok=True)
    except OSError as e:
        sys.exit(f"[ERROR] Cannot create output directory: {e}")

    chroms = list(dict.fromkeys(df["Chromosome"].tolist()))
    total = len(chroms)
    rolling_capped = False

    for i, chrom in enumerate(chroms, 1):
        cdf = df[df["Chromosome"] == chrom].copy()
        if cdf.empty:
            continue

        contig_len = int(cdf["End"].max())
        divisor, unit = smart_bp_formatter(contig_len)
        x = cdf["Start"].values / divisor

        win_bp = int(cdf["End"].iloc[0] - cdf["Start"].iloc[0])
        win_div, win_unit = smart_bp_formatter(win_bp)
        win_str = f"{win_bp / win_div:g} {win_unit}"
        if len(cdf) > 1:
            step_bp = int(cdf["Start"].iloc[1] - cdf["Start"].iloc[0])
            step_div, step_unit = smart_bp_formatter(max(step_bp, 1))
            step_str = f"{step_bp / step_div:g} {step_unit}" if step_bp > 0 else "0 bp"
            if step_bp > 0:
                uncapped_n = max(3, round(rolling_target_bp / step_bp))
                rolling_n = min(uncapped_n, len(cdf))
                if uncapped_n > len(cdf):
                    rolling_capped = True
                rolling_span_bp = rolling_n * step_bp
            else:
                rolling_n = 1
                rolling_span_bp = win_bp
        else:
            step_str = "N/A (single window)"
            rolling_n = 1
            rolling_span_bp = win_bp

        rolling_span_div, rolling_span_unit = smart_bp_formatter(rolling_span_bp)
        rolling_span_str = f"{rolling_span_bp / rolling_span_div:g} {rolling_span_unit}"

        cov = cdf["Coverage"].values
        mean_cov = cdf["Coverage"].mean()
        median_cov = cdf["Coverage"].median()

        med_mapq_win = cdf["Median_MAPQ"].values
        mean_mapq_win = cdf["Mean_MAPQ"].values
        sc_win = cdf["Softclip_%"].values

        row = df_sum_ctg[df_sum_ctg["Chromosome"] == chrom]
        contig_med_mapq = float(row["Median_MAPQ"].iloc[0]) if not row.empty else np.nan
        contig_mean_mapq = float(row["Mean_MAPQ"].iloc[0]) if not row.empty else np.nan
        contig_softclip = float(row["Softclip_%"].iloc[0]) if not row.empty else np.nan

        cov_rolling = (
            pd.Series(cov).rolling(rolling_n, center=True, min_periods=1).mean().values
        )
        med_mapq_rolling = (
            pd.Series(med_mapq_win)
            .rolling(rolling_n, center=True, min_periods=1)
            .mean()
            .values
        )
        mean_mapq_rolling = (
            pd.Series(mean_mapq_win)
            .rolling(rolling_n, center=True, min_periods=1)
            .mean()
            .values
        )
        sc_rolling = (
            pd.Series(sc_win)
            .rolling(rolling_n, center=True, min_periods=1)
            .mean()
            .values
        )

        fig, (ax_cov, ax_mapq, ax_sc) = plt.subplots(
            3, 1, figsize=(18, 9), sharex=True, constrained_layout=True
        )
        apply_style(fig, [ax_cov, ax_mapq, ax_sc])

        xlim = (cdf["Start"].min() / divisor, cdf["End"].max() / divisor)

        # coverage panel
        cov_max = _safe_nanmax(cov_rolling, default=1.0)
        cov_ylim_top = max(cov_max * 1.05, 1.0)

        ax_cov.plot(
            x,
            cov_rolling,
            color=COV_A,
            linewidth=1.5,
            zorder=4,
            label="Rolling mean cov",
        )
        ax_cov.axhline(
            mean_cov, **MEAN_STYLE_DASH, label=f"Contig mean cov {mean_cov:.1f}×"
        )
        ax_cov.axhline(
            median_cov,
            **MEDIAN_STYLE_DASH,
            label=f"Contig median cov {median_cov:.1f}×",
        )
        ax_cov.set_ylabel("Coverage", color=TEXT_DIM, fontsize=9, fontfamily=FONT_MAIN)
        ax_cov.set_ylim(0, cov_ylim_top)
        ax_cov.set_xlim(*xlim)
        ax_cov.legend(
            fontsize=8,
            facecolor=AX_BG,
            edgecolor=GRID,
            labelcolor=TEXT,
            framealpha=0.9,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            bbox_transform=ax_cov.transAxes,
        )
        ax_cov.set_title(
            f"Coverage, MAPQ & Softclip — rolling mean — {chrom}\n"
            f"window size: {win_str}  |  step size: {step_str}  |  rolling mean: {rolling_n} windows ({rolling_span_str})",
            color=TEXT,
            fontsize=12,
            fontweight="bold",
            fontfamily=FONT_MAIN,
        )

        # MAPQ panel, mean and median combined with fill between
        ax_mapq.fill_between(
            x,
            mean_mapq_rolling,
            med_mapq_rolling,
            alpha=0.35,
            color=GRID,
            linewidth=0,
            zorder=3,
        )
        ax_mapq.plot(
            x,
            mean_mapq_rolling,
            color=MAPQ_B,
            linewidth=1.0,
            zorder=4,
            label="Rolling mean MAPQ",
        )
        ax_mapq.plot(
            x,
            med_mapq_rolling,
            color=MAPQ_A,
            linewidth=1.0,
            zorder=4,
            label="Rolling median MAPQ",
        )
        ax_mapq.axhline(
            contig_mean_mapq,
            **MEAN_STYLE_DASH,
            label=f"Contig mean MAPQ  {contig_mean_mapq:.2f}",
        )
        ax_mapq.axhline(
            contig_med_mapq,
            **MEDIAN_STYLE_DASH,
            label=f"Contig median MAPQ  {contig_med_mapq:.1f}",
        )
        mapq_data_max = _safe_nanmax(
            [*mean_mapq_rolling, *med_mapq_rolling, contig_mean_mapq, contig_med_mapq], default=np.nan
        )
        ylim_top_mapq = (
            max(65, mapq_data_max * 1.1) if np.isfinite(mapq_data_max) else 65
        )
        ax_mapq.set_ylabel("MAPQ", color=TEXT_DIM, fontsize=9, fontfamily=FONT_MAIN)
        ax_mapq.set_ylim(0, ylim_top_mapq)
        ax_mapq.legend(
            fontsize=8,
            facecolor=AX_BG,
            edgecolor=GRID,
            labelcolor=TEXT,
            framealpha=0.9,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            bbox_transform=ax_mapq.transAxes,
        )

        # softclip panel
        sc_max = _safe_nanmax(sc_rolling, default=0.1)
        sc_ylim_top = max(sc_max * 1.05, 0.1)

        ax_sc.plot(
            x,
            sc_rolling,
            color=SC_A,
            linewidth=1.2,
            zorder=4,
            label="Rolling mean softclip %",
        )
        ax_sc.axhline(
            contig_softclip,
            **{**SOFTCLIP_STYLE, "linestyle": "--"},
            label=f"Contig softclip %  {contig_softclip:.2f}",
        )
        ax_sc.set_ylabel("Softclip %", color=TEXT_DIM, fontsize=9, fontfamily=FONT_MAIN)
        ax_sc.set_ylim(0, sc_ylim_top)
        ax_sc.legend(
            fontsize=8,
            facecolor=AX_BG,
            edgecolor=GRID,
            labelcolor=TEXT,
            framealpha=0.9,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            bbox_transform=ax_sc.transAxes,
        )
        ax_sc.set_xlabel(
            f"Position ({unit})", color=TEXT_DIM, fontsize=9, fontfamily=FONT_MAIN
        )
        ax_sc.xaxis.set_major_formatter(
            ticker.FuncFormatter(
                lambda v, _: f"{v:,.1f}" if divisor > 1 else f"{int(v):,}"
            )
        )

        safe = _sanitise_filename(chrom)
        out = os.path.join(rolling_folder, f"{safe}_per_window_rolling.png")
        fig.savefig(out, dpi=150, facecolor=BG, bbox_inches="tight")
        plt.close(fig)

        if i % 50 == 0 or i == total:
            print(f"[INFO] Rolling per-contig figures: {i}/{total}")

    if rolling_capped:
        rolling_div, rolling_unit = smart_bp_formatter(int(rolling_target_bp))
        print(
            f"[WARN] -r/--rolling ({rolling_target_bp / rolling_div:g} {rolling_unit}) "
            f"exceeds the window count of one or more contigs. Rolling window was capped at contig length for those contigs."
        )
    print(f"[INFO] Rolling per-contig figures saved to: {rolling_folder}")


def run_mapq_softclip(input_folder, rolling_target_bp):
    """Load mapq_softclip output files and generate all genome-wide and per-contig figures."""
    window_file = os.path.join(input_folder, "window_stats.csv")
    summary_file = os.path.join(input_folder, "summary_stats.csv")
    try:
        figures_folder = make_figures_folder(input_folder)
    except OSError as e:
        sys.exit(f"[ERROR] Cannot create figures directory: {e}")

    print("[INFO] Loading window_stats.csv ...")
    try:
        df_win = pd.read_csv(window_file)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read window_stats.csv: {e}")
    df_win.columns = df_win.columns.str.strip()
    _check_columns(df_win, REQUIRED_WIN_COLS, "window_stats.csv")

    print("[INFO] Loading summary_stats.csv ...")
    try:
        df_sum = pd.read_csv(summary_file)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read summary_stats.csv: {e}")
    df_sum.columns = df_sum.columns.str.strip()
    _check_columns(df_sum, REQUIRED_SUM_COLS, "summary_stats.csv")

    print("[INFO] Generating genome-wide summary figure ...")
    fig_genome_three_panel(df_win, df_sum, figures_folder)

    print("[INFO] Generating raw per-contig figures ...")
    fig_per_contig_raw(df_win, df_sum, figures_folder)
    print("[INFO] Generating rolling per-contig figures ...")
    fig_per_contig_rolling(df_win, df_sum, figures_folder, rolling_target_bp)

    print(f"\n[INFO] All figures saved to: {figures_folder}")


# ============================================================
# hese: styling
# ============================================================

HAP_P1  = "#0969da"
HAP_P2  = "#e36209"
HAP_AMB = "#8c959f"

REQUIRED_HAP_COLS = {
    "hap_id", "hap_type", "n_unitigs", "n_p1", "n_p2", "n_amb", "label",
}
REQUIRED_UNITIG_COLS = {"unitig", "p1_norm", "p2_norm", "label"}
REQUIRED_SUM_COLS_HESE = {
    "hap", "hap_global_parent", "assigned_parent", "hamming_struct_%", "switch_struct_%",
}
REQUIRED_TRUTH_EVAL_COLS = {
    "overall_hamming_%", "overall_switch_%",
}

ERR_HAMMING        = "#cf222e"
ERR_SWITCH         = "#bf8700"
ERR_HAMMING_TRUTH  = "#1a7f37"
ERR_SWITCH_TRUTH   = "#0e7490"


# ============================================================
# hese: figures
# ============================================================

def fig_hese_label_balance(df_hap, figures_folder):
    """Save hese_label_balance.png, a horizontal stacked bar of P1, amb and P2 haplotig counts."""
    counts = df_hap["label"].value_counts()
    n_p1  = int(counts.get("P1",  0))
    n_p2  = int(counts.get("P2",  0))
    n_amb = int(counts.get("amb", 0))
    total = n_p1 + n_p2 + n_amb
    if total == 0:
        print("[WARN] No haplotigs found, skipping label balance figure.")
        return

    ratio_str = f"{n_p1/n_p2:.1f}:1" if n_p2 > 0 else "inf"

    fig, ax = plt.subplots(figsize=(5, 1.8), constrained_layout=True)
    apply_style(fig, ax)

    segments = [
        ("P1",  n_p1,  HAP_P1),
        ("amb", n_amb, HAP_AMB),
        ("P2",  n_p2,  HAP_P2),
    ]
    left = 0
    for label, count, color in segments:
        ax.barh(0, count, left=left, color=color, height=0.5, zorder=3)
        if count / total >= 0.04:
            ax.text(
                left + count / 2, 0,
                f"{label}\n{count}\n({count/total*100:.1f}%)",
                ha="center", va="center",
                fontsize=8, color="white", fontfamily=FONT_MAIN, fontweight="bold",
            )
        left += count

    ax.set_xlim(0, total)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel("Number of haplotigs", color=TEXT_DIM, fontsize=10, fontfamily=FONT_MAIN)
    ax.set_yticks([])
    ax.set_title(
        f"Haplotig label balance  P1:P2 = {ratio_str}  (n={total})",
        color=TEXT, fontsize=12, fontweight="bold", fontfamily=FONT_MAIN,
    )

    out = os.path.join(figures_folder, "hese_label_balance.png")
    fig.savefig(out, dpi=180, facecolor=BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out}")


def fig_hese_signal_penetration(df_hap, min_unitigs, figures_folder):
    """Save hese_signal_penetration.png, a histogram of informative unitig counts per haplotig with a threshold line at min_unitigs."""
    informative = (df_hap["n_p1"] + df_hap["n_p2"]).values
    total = len(informative)
    if total == 0:
        print("[WARN] No haplotigs found, skipping signal penetration figure.")
        return

    max_val = int(informative.max())
    bins = np.arange(0, max_val + 2) - 0.5
    n_below = int((informative < min_unitigs).sum())
    n_above = total - n_below

    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    apply_style(fig, ax)

    ax.hist(informative, bins=bins, color=HAP_AMB, zorder=3, linewidth=0)
    ax.axvline(
        min_unitigs - 0.5, color=REF_LINE, linewidth=1.4, linestyle="--", zorder=4,
        label=f"labeling threshold (min-unitigs = {min_unitigs})",
    )
    ax.text(
        0.98, 0.97,
        f"{n_below} haplotigs below threshold ({n_below/total*100:.1f}%)\n"
        f"{n_above} haplotigs above threshold ({n_above/total*100:.1f}%)",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=8, color=TEXT_DIM, fontfamily=FONT_MAIN,
    )

    ax.set_xlabel(
        "Informative (non-amb) unitigs per haplotig",
        color=TEXT_DIM, fontsize=10, fontfamily=FONT_MAIN,
    )
    ax.set_ylabel("Number of haplotigs", color=TEXT_DIM, fontsize=10, fontfamily=FONT_MAIN)
    ax.set_title(
        "Parental signal penetration per haplotig",
        color=TEXT, fontsize=12, fontweight="bold", fontfamily=FONT_MAIN,
    )
    ax.legend(fontsize=9, facecolor=AX_BG, edgecolor=GRID, labelcolor=TEXT, framealpha=0.9)
    ax.set_xlim(-0.5, max_val + 1)

    out = os.path.join(figures_folder, "hese_signal_penetration.png")
    fig.savefig(out, dpi=180, facecolor=BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out}")


def fig_hese_label_distributions(df_uni, df_hap, figures_folder):
    """Save hese_label_distributions.png, side-by-side bar charts of unitig and haplotig label counts."""
    label_order = ["P1", "P2", "amb"]
    colors = {"P1": HAP_P1, "P2": HAP_P2, "amb": HAP_AMB}

    uni_counts = df_uni["label"].value_counts()
    hap_counts = df_hap["label"].value_counts()

    fig, (ax_uni, ax_hap) = plt.subplots(1, 2, figsize=(7, 4), constrained_layout=True)
    apply_style(fig, [ax_uni, ax_hap])

    for ax, counts, total, title, ylabel in [
        (ax_uni, uni_counts, len(df_uni), "Unitig label distribution", "Number of unitigs"),
        (ax_hap, hap_counts, len(df_hap), "Haplotig label distribution", "Number of haplotigs"),
    ]:
        vals = [int(counts.get(lbl, 0)) for lbl in label_order]
        bars = ax.bar(
            label_order, vals,
            color=[colors[lbl] for lbl in label_order],
            zorder=3, width=0.5,
        )
        for bar, val in zip(bars, vals):
            if total > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + total * 0.01,
                    f"{val}\n({val/total*100:.1f}%)",
                    ha="center", va="bottom",
                    fontsize=8, color=TEXT_DIM, fontfamily=FONT_MAIN,
                )
        ax.set_ylabel(ylabel, color=TEXT_DIM, fontsize=10, fontfamily=FONT_MAIN)
        ax.set_title(title, color=TEXT, fontsize=12, fontweight="bold", fontfamily=FONT_MAIN)
        ax.set_ylim(0, max(vals) * 1.18 if max(vals) > 0 else 1)
        ax.tick_params(axis="x", labelsize=10)

    out = os.path.join(figures_folder, "hese_label_distributions.png")
    fig.savefig(out, dpi=180, facecolor=BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out}")


def fig_hese_signal_depth(df_hap, figures_folder):
    """Save hese_signal_depth.png, showing informative unitig counts per labeled haplotig as a jittered scatter plot."""
    labeled = df_hap[df_hap["label"].isin(["P1", "P2"])].copy()
    labeled["informative"] = labeled["n_p1"] + labeled["n_p2"]

    if labeled.empty:
        print("[WARN] No labeled (P1/P2) haplotigs found, skipping signal depth figure.")
        return

    fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
    apply_style(fig, ax)

    for lbl, color, xpos in [("P1", HAP_P1, 0), ("P2", HAP_P2, 1)]:
        sub = labeled[labeled["label"] == lbl]["informative"].values
        if len(sub) == 0:
            continue
        jitter = np.random.default_rng(42).uniform(-0.12, 0.12, size=len(sub))
        ax.scatter(
            np.full(len(sub), xpos) + jitter, sub,
            color=color, s=25, alpha=0.7, linewidths=0, zorder=3,
            label=f"{lbl}  (n={len(sub)})",
        )
        ax.plot(
            [xpos - 0.2, xpos + 0.2], [np.median(sub), np.median(sub)],
            color=color, linewidth=2, zorder=4,
        )

    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        ["haplotigs labeled P1", "haplotigs labeled P2"],
        fontsize=10, color=TEXT_DIM, fontfamily=FONT_MAIN,
    )
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel(
        "Informative (non-amb) unitigs inside the haplotig",
        color=TEXT_DIM, fontsize=10, fontfamily=FONT_MAIN,
    )
    ax.set_title(
        "Parental signal depth per labeled haplotig",
        color=TEXT, fontsize=12, fontweight="bold", fontfamily=FONT_MAIN,
    )
    ax.legend(fontsize=9, facecolor=AX_BG, edgecolor=GRID, labelcolor=TEXT, framealpha=0.9)

    out = os.path.join(figures_folder, "hese_signal_depth.png")
    fig.savefig(out, dpi=180, facecolor=BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out}")


def fig_hese_phasing_errors(df_sum, df_truth_eval, figures_folder):
    """Save hese_phasing_errors.png, Hamming and switch error bars per haplotype and parent assignment.

    If df_truth_eval is provided, overall truth reference lines are added to the plot.
    """
    groups = []
    global_parents = {}
    for hap in ["hap1", "hap2"]:
        hap_rows = df_sum[df_sum["hap"] == hap]
        if hap_rows.empty:
            continue
        gp = str(hap_rows["hap_global_parent"].iloc[0]).strip()
        global_parents[hap] = gp
        best = hap_rows[hap_rows["assigned_parent"] == gp]
        if best.empty:
            continue
        groups.append((hap, gp, best))

    bar_width = 0.3
    x = np.arange(len(groups))

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    apply_style(fig, ax)

    all_vals = []
    for i, (metric_col, color, legend_label) in enumerate([
        ("hamming_struct_%", ERR_HAMMING, "Hamming structural"),
        ("switch_struct_%",  ERR_SWITCH,  "Switch structural"),
    ]):
        offset = (i - 0.5) * bar_width
        vals = []
        for hap, assigned, row in groups:
            v = row[metric_col].iloc[0] if not row.empty else np.nan
            vals.append(float(v) if v != "NA" else np.nan)
        all_vals.extend(vals)
        bars = ax.bar(x + offset, vals, width=bar_width * 0.9, color=color, label=legend_label, zorder=3)
        for bar, val in zip(bars, vals):
            if np.isfinite(val) and val > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.3,
                    f"{val:.1f}%",
                    ha="center", va="bottom",
                    fontsize=8, color=TEXT_DIM, fontfamily=FONT_MAIN,
                )

    if df_truth_eval is not None and not df_truth_eval.empty:
        hamming_truth = float(df_truth_eval["overall_hamming_%"].iloc[0])
        switch_truth  = float(df_truth_eval["overall_switch_%"].iloc[0])
        ax.axhline(hamming_truth, color=ERR_HAMMING_TRUTH, linewidth=1.4, linestyle="--", zorder=4,
                   label=f"Hamming truth (overall)  {hamming_truth:.1f}%")
        ax.axhline(switch_truth,  color=ERR_SWITCH_TRUTH,  linewidth=1.4, linestyle="--", zorder=4,
                   label=f"Switch truth (overall)  {switch_truth:.1f}%")
        all_vals.extend([hamming_truth, switch_truth])

    xtick_labels = [f"{hap}\nassigned {assigned}" for hap, assigned, _ in groups]
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels, fontsize=9, color=TEXT_DIM, fontfamily=FONT_MAIN)
    ax.set_ylabel("Error %", color=TEXT_DIM, fontsize=10, fontfamily=FONT_MAIN)
    finite_vals = [v for v in all_vals if np.isfinite(v)]
    max_val = max(finite_vals) if finite_vals else 0
    ax.set_ylim(0, max(max_val * 1.35 + 0.5, 2.0))

    shared_parent = None
    if len(global_parents) == 2 and len(set(global_parents.values())) == 1:
        shared_parent = next(iter(global_parents.values()))

    if shared_parent is not None:
        title_text = f"Phasing errors per haplotype\nboth haplotypes claimed {shared_parent}"
        title_color = ERR_HAMMING
    else:
        title_text = "Phasing errors per haplotype"
        title_color = TEXT

    ax.set_title(
        title_text,
        color=title_color, fontsize=12, fontweight="bold", fontfamily=FONT_MAIN,
    )
    ax.legend(fontsize=9, facecolor=AX_BG, edgecolor=GRID, labelcolor=TEXT, framealpha=0.9)

    out = os.path.join(figures_folder, "hese_phasing_errors.png")
    fig.savefig(out, dpi=180, facecolor=BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out}")


def _read_min_unitigs(input_folder, default=3):
    """Parse the min-unitigs threshold used in the hese run from run_log.txt, default if not found."""
    log_path = os.path.join(input_folder, "run_log.txt")
    if not os.path.isfile(log_path):
        return default
    try:
        with open(log_path) as f:
            for line in f:
                if "Min unitigs" in line and "(-mu)" in line:
                    return int(line.split(":")[1].strip())
    except (OSError, ValueError, IndexError):
        pass
    return default


def run_hese(input_folder):
    """Load hese output files and generate all hese figures."""
    try:
        figures_folder = make_figures_folder(input_folder)
    except OSError as e:
        sys.exit(f"[ERROR] Cannot create figures directory: {e}")

    hap_file = os.path.join(input_folder, "haplotig_labels.csv")
    print("[INFO] Loading haplotig_labels.csv ...")
    try:
        df_hap = pd.read_csv(hap_file)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read haplotig_labels.csv: {e}")
    df_hap.columns = df_hap.columns.str.strip()
    _check_columns(df_hap, REQUIRED_HAP_COLS, "haplotig_labels.csv")

    uni_file = os.path.join(input_folder, "unitig_labels.csv")
    print("[INFO] Loading unitig_labels.csv ...")
    try:
        df_uni = pd.read_csv(uni_file)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read unitig_labels.csv: {e}")
    df_uni.columns = df_uni.columns.str.strip()
    _check_columns(df_uni, REQUIRED_UNITIG_COLS, "unitig_labels.csv")

    sum_file = os.path.join(input_folder, "haplotype_summary.csv")
    print("[INFO] Loading haplotype_summary.csv ...")
    try:
        df_sum = pd.read_csv(sum_file)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read haplotype_summary.csv: {e}")
    df_sum.columns = df_sum.columns.str.strip()
    _check_columns(df_sum, REQUIRED_SUM_COLS_HESE, "haplotype_summary.csv")

    truth_eval_file = os.path.join(input_folder, "truth_eval_summary.csv")
    df_truth_eval = None
    if os.path.isfile(truth_eval_file):
        print("[INFO] Loading truth_eval_summary.csv ...")
        try:
            df_truth_eval = pd.read_csv(truth_eval_file)
        except Exception as e:
            sys.exit(f"[ERROR] Failed to read truth_eval_summary.csv: {e}")
        df_truth_eval.columns = df_truth_eval.columns.str.strip()
        _check_columns(df_truth_eval, REQUIRED_TRUTH_EVAL_COLS, "truth_eval_summary.csv")

    min_unitigs = _read_min_unitigs(input_folder)
    print(f"[INFO] Using min-unitigs = {min_unitigs} (from run_log.txt)")

    print("[INFO] Generating label distributions figure ...")
    fig_hese_label_distributions(df_uni, df_hap, figures_folder)

    print("[INFO] Generating label balance figure ...")
    fig_hese_label_balance(df_hap, figures_folder)

    print("[INFO] Generating signal penetration figure ...")
    fig_hese_signal_penetration(df_hap, min_unitigs, figures_folder)

    print("[INFO] Generating signal depth figure ...")
    fig_hese_signal_depth(df_hap, figures_folder)

    print("[INFO] Generating phasing errors figure ...")
    fig_hese_phasing_errors(df_sum, df_truth_eval, figures_folder)

    print(f"\n[INFO] All figures saved to: {figures_folder}")


# ============================================================
# Entry point
# ============================================================

def main():
    """Entry point: parse arguments and dispatch to the appropriate module visualiser."""
    args = parse_arguments()
    check_input(args.input, args.module)

    if math.isnan(args.rolling) or math.isinf(args.rolling):
        sys.exit(f"[ERROR] -r/--rolling must be a finite number, got {args.rolling}")
    if args.rolling < 1:
        sys.exit("[ERROR] -r/--rolling must be at least 1 kb")

    try:
        if args.module == "mapq_softclip":
            run_mapq_softclip(args.input, args.rolling * 1000)
        elif args.module == "hese":
            run_hese(args.input)
    except KeyboardInterrupt:
        sys.exit("\n[INFO] Interrupted.")


if __name__ == "__main__":
    main()
