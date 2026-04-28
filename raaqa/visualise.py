import argparse
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

BG = "#ffffff"  # figure background
AX_BG = "#f6f8fa"  # axes background
GRID = "#d0d7de"  # grid lines
TEXT = "#1f2328"  # titles and primary labels
TEXT_DIM = "#57606a"  # axis labels and tick labels
COV_A = "#0969da"  # coverage fill; mean coverage dots/lines in genome-wide
COV_B = "#54aeff"  # median coverage dots/lines in genome-wide
COV_LINE = "#0550ae"  # coverage overflow triangle markers
MAPQ_A = "#79a80a"  # median MAPQ fill; median MAPQ dots/lines in genome-wide
MAPQ_B = "#0e9488"  # mean MAPQ fill; mean MAPQ dots/lines in genome-wide
SC_A = "#e36209"  # softclip fill; softclip dots/lines in genome-wide
SC_LINE = "#953800"  # softclip reference axhlines and softclip overflow markers
REF_LINE = "#cf222e"  # mean reference axhlines
MED_LINE = "#8250df"  # median reference axhlines
FONT_MAIN = "DejaVu Sans"

MEAN_STYLE = dict(color=REF_LINE, linewidth=1.4, linestyle="-", zorder=5)
MEDIAN_STYLE = dict(color=MED_LINE, linewidth=1.4, linestyle="-", zorder=5)
SOFTCLIP_STYLE = dict(color=SC_LINE, linewidth=1, linestyle="-", zorder=5)
MEAN_STYLE_DASH = dict(color=REF_LINE, linewidth=1, linestyle="--", zorder=5)
MEDIAN_STYLE_DASH = dict(color=MED_LINE, linewidth=1, linestyle="--", zorder=5)

REQUIRED_WIN_COLS = {
    "Chromosome",
    "Start",
    "End",
    "Total_Bases",
    "Median_MAPQ",
    "Mean_MAPQ",
    "Softclip_%",
}
REQUIRED_SUM_COLS = {"Chromosome", "Mean_MAPQ", "Median_MAPQ", "Softclip_%"}


def _check_columns(df, required, filename):
    missing = required - set(df.columns)
    if missing:
        sys.exit(
            f"[ERROR] {filename} is missing required columns: {', '.join(sorted(missing))}"
        )


def apply_style(fig, axes):
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
    if contig_length_bp >= 1_000_000:
        return 1_000_000, "Mb"
    elif contig_length_bp >= 1_000:
        return 1_000, "kb"
    else:
        return 1, "bp"


def _sanitise_filename(name):
    return re.sub(r"[^\w.\-]", "_", name)


def make_figures_folder(input_folder):
    figures_folder = os.path.join(input_folder, "figures")
    os.makedirs(figures_folder, exist_ok=True)
    return figures_folder


# Argument parsing and input validation
def parse_arguments():
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
        help="Path to the output folder from the respective module run",
    )
    parser.add_argument(
        "-f",
        "--figures",
        choices=["genome", "contigs", "all"],
        default="all",
        help="Which figures to generate: genome-wide summary only, per-contig only, or all (default: all)",
    )
    parser.add_argument(
        "-r",
        "--rolling",
        type=int,
        default=1000,
        help="Genomic distance in kb to span for rolling mean in rolling figures (default: 1000)",
    )
    return parser.parse_args()


def check_input(input_folder, module):
    if not os.path.isdir(input_folder):
        sys.exit(f"[ERROR] Input folder not found: {input_folder}")
    if module == "mapq_softclip":
        for f in ["window_stats.csv", "summary_stats.csv"]:
            if not os.path.exists(os.path.join(input_folder, f)):
                sys.exit(f"[ERROR] Required file not found: {f}")


# Shared helper — build genome x axis from window data
def _build_genome_xaxis(df, chroms):
    genome_pos = 0
    dividers = []
    chrom_ticks = []
    chrom_labels = []
    chrom_x = {}

    for chrom in chroms:
        cdf = df[df["Chromosome"] == chrom]
        if cdf.empty:
            continue
        chrom_x[chrom] = (cdf, genome_pos)
        chrom_ticks.append(genome_pos + int(cdf["End"].max()) / 2)
        chrom_labels.append(chrom)
        genome_pos += int(cdf["End"].max())
        dividers.append(genome_pos)

    return genome_pos, dividers, chrom_ticks, chrom_labels, chrom_x


def _apply_genome_xticks(ax, chrom_ticks, chrom_labels):
    ax.set_xticks(chrom_ticks)
    ax.set_xticklabels(
        chrom_labels, rotation=90, fontsize=4.5, color=TEXT_DIM, fontfamily=FONT_MAIN
    )


def _draw_dividers(axes, dividers):
    for d in dividers[:-1]:
        for ax in axes:
            ax.axvline(d, color=GRID, linewidth=0.5, alpha=0.6)


# Genome-wide 3-panel: coverage + per-contig MAPQ lines + per-contig softclip lines
def fig_genome_three_panel(df_win, df_sum, figures_folder):
    df = df_win[df_win["Chromosome"] != "GENOME"].copy()
    df["Coverage"] = df["Total_Bases"] / (df["End"] - df["Start"]).replace(0, np.nan)
    df_sum_ctg = df_sum[df_sum["Chromosome"] != "GENOME"].copy()

    # genome-wide reference values — pre-computed in GENOME row of summary_stats.csv
    genome_row = df_sum[df_sum["Chromosome"] == "GENOME"]
    if genome_row.empty:
        print(
            "[WARN] No GENOME row found in summary_stats.csv — genome-wide reference lines will be omitted."
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

        mid = genome_offset + int(cdf["End"].max()) / 2
        contig_mids.append(mid)

        # coverage — derived from window data
        contig_cov_means.append(float(cdf["Coverage"].mean()))
        contig_cov_meds.append(float(cdf["Coverage"].median()))

        # MAPQ and softclip — pre-computed in summary_stats.csv
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

    # genome-wide coverage reference lines — derived
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

    ax_cov.set_xlim(0, genome_pos)
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

    mapq_data_max = np.nanmax(
        [*contig_means, *contig_meds, genome_mean_mapq, genome_med_mapq]
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


# Per-contig 4-panel: coverage + median MAPQ + mean MAPQ + softclip %
def fig_per_contig_raw(df_win, df_sum, figures_folder):
    df = df_win[df_win["Chromosome"] != "GENOME"].copy()
    df["Coverage"] = df["Total_Bases"] / (df["End"] - df["Start"]).replace(0, np.nan)
    df_sum_ctg = df_sum[df_sum["Chromosome"] != "GENOME"].copy()

    if df.empty:
        print(
            "[WARN] No per-contig data found in window_stats.csv — skipping raw contig figures."
        )
        return

    contigs_folder = os.path.join(figures_folder, "contigs", "raw")
    os.makedirs(contigs_folder, exist_ok=True)

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

        # coverage — derived
        cov = cdf["Coverage"].values
        mean_cov = cdf["Coverage"].mean()
        median_cov = cdf["Coverage"].median()

        # MAPQ and softclip per window — pre-computed in window_stats.csv
        med_mapq_win = cdf["Median_MAPQ"].values
        mean_mapq_win = cdf["Mean_MAPQ"].values
        sc_win = cdf["Softclip_%"].values

        # contig-level reference values — pre-computed in summary_stats.csv
        row = df_sum_ctg[df_sum_ctg["Chromosome"] == chrom]
        contig_med_mapq = float(row["Median_MAPQ"].iloc[0]) if not row.empty else np.nan
        contig_mean_mapq = float(row["Mean_MAPQ"].iloc[0]) if not row.empty else np.nan
        contig_softclip = float(row["Softclip_%"].iloc[0]) if not row.empty else np.nan

        fig, (ax_cov, ax_med, ax_mean, ax_sc) = plt.subplots(
            4, 1, figsize=(18, 11), sharex=True, constrained_layout=True
        )
        apply_style(fig, [ax_cov, ax_med, ax_mean, ax_sc])

        xlim = (cdf["Start"].min() / divisor, cdf["End"].max() / divisor)

        # — Coverage panel —
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

        mapq_data_max = np.nanmax(
            [*mean_mapq_win, *med_mapq_win, contig_mean_mapq, contig_med_mapq]
        )
        ylim_top_mapq = (
            max(65, mapq_data_max * 1.1) if np.isfinite(mapq_data_max) else 65
        )

        # — Median MAPQ panel — pre-computed per window —
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

        # — Mean MAPQ panel — pre-computed per window —
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

        # — Softclip % panel — pre-computed per window —
        sc_max = np.nanmax(sc_win) if len(sc_win) else 0.1
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
    df = df_win[df_win["Chromosome"] != "GENOME"].copy()
    df["Coverage"] = df["Total_Bases"] / (df["End"] - df["Start"]).replace(0, np.nan)
    df_sum_ctg = df_sum[df_sum["Chromosome"] != "GENOME"].copy()

    if df.empty:
        print(
            "[WARN] No per-contig data found in window_stats.csv — skipping rolling contig figures."
        )
        return

    rolling_folder = os.path.join(figures_folder, "contigs", "rolling")
    os.makedirs(rolling_folder, exist_ok=True)

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
            rolling_n = min(max(3, round(rolling_target_bp / step_bp)), len(cdf))
            rolling_span_bp = rolling_n * step_bp
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

        # — Coverage panel —
        cov_max = np.nanmax(cov_rolling) if len(cov_rolling) else 1.0
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

        # — Merged MAPQ panel —
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
        mapq_data_max = np.nanmax(
            [*mean_mapq_rolling, *med_mapq_rolling, contig_mean_mapq, contig_med_mapq]
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

        # — Softclip panel —
        sc_max = np.nanmax(sc_rolling) if len(sc_rolling) else 0.1
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

    print(f"[INFO] Rolling per-contig figures saved to: {rolling_folder}")


# mapq_softclip
def run_mapq_softclip(input_folder, figures, rolling_target_bp):
    window_file = os.path.join(input_folder, "window_stats.csv")
    summary_file = os.path.join(input_folder, "summary_stats.csv")
    figures_folder = make_figures_folder(input_folder)

    print("[INFO] Loading window_stats.csv ...")
    df_win = pd.read_csv(window_file)
    _check_columns(df_win, REQUIRED_WIN_COLS, "window_stats.csv")

    print("[INFO] Loading summary_stats.csv ...")
    df_sum = pd.read_csv(summary_file)
    _check_columns(df_sum, REQUIRED_SUM_COLS, "summary_stats.csv")

    if figures in ("genome", "all"):
        print("[INFO] Generating genome-wide summary figure ...")
        fig_genome_three_panel(df_win, df_sum, figures_folder)

    if figures in ("contigs", "all"):
        print("[INFO] Generating raw per-contig figures ...")
        fig_per_contig_raw(df_win, df_sum, figures_folder)
        print("[INFO] Generating rolling per-contig figures ...")
        fig_per_contig_rolling(df_win, df_sum, figures_folder, rolling_target_bp)

    print(f"\n[INFO] All figures saved to: {figures_folder}")


# hese
def run_hese(input_folder):
    print("[INFO] HESE visualisation is not yet implemented.")


def main():
    args = parse_arguments()
    check_input(args.input, args.module)

    if args.rolling < 1:
        sys.exit("[ERROR] -r/--rolling must be a positive integer (kb)")

    if args.module == "mapq_softclip":
        run_mapq_softclip(args.input, args.figures, args.rolling * 1000)
    elif args.module == "hese":
        run_hese(args.input)


if __name__ == "__main__":
    main()
