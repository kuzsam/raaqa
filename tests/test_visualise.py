"""
Tests for the visualise module.

Organised as:
  UNIT TESTS
    mapq_softclip helpers  - smart_bp_formatter, _sanitise_filename, _check_columns,
                             make_figures_folder, check_input, _build_genome_xaxis
    hese utilities         - _to_float_or_nan, _safe_nanmax, _read_min_unitigs

  INTEGRATION TESTS
    mapq_softclip figures  - run_mapq_softclip called with hand-crafted CSVs or a
                             real BAM pipeline, asserting figure files are created
    hese figures           - figure functions called with hand-crafted DataFrames,
                             asserting file creation, warning on empty input, handling
                             NA metric values, and run_hese end-to-end behaviour
"""

import math
import os

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import pysam
import pytest

from raaqa.mapq_softclip import run_analysis
from raaqa.visualise import (
    smart_bp_formatter,
    _sanitise_filename,
    _safe_nanmax,
    _to_float_or_nan,
    _check_columns,
    make_figures_folder,
    check_input,
    _build_genome_xaxis,
    _read_min_unitigs,
    fig_hese_label_balance,
    fig_hese_signal_penetration,
    fig_hese_label_distributions,
    fig_hese_signal_depth,
    fig_hese_phasing_errors,
    run_hese,
    run_mapq_softclip,
)


# ══════════════════════════════════════════════════════════════════════════════
# UNIT TESTS
# ══════════════════════════════════════════════════════════════════════════════

# ── helpers ───────────────────────────────────────────────────────────────────

def _write_csv(path, content):
    with open(path, "w") as f:
        f.write(content)


def _make_minimal_csvs(folder, window_content=None, summary_content=None):
    """Write minimal valid window_stats.csv and summary_stats.csv to folder."""
    if window_content is None:
        window_content = (
            "Chromosome,Start,End,Mean_MAPQ,Median_MAPQ,Read_Count,"
            "Total_Bases,Softclip_Bases,Softclip_%,Flag\n"
            "chr1,0,1000,60.00,60,10,10000,0,0.00000,\n"
            "chr1,1000,2000,60.00,60,10,10000,0,0.00000,\n"
            "chr1,2000,3000,60.00,60,10,10000,0,0.00000,\n"
        )
    if summary_content is None:
        summary_content = (
            "Chromosome,Mean_MAPQ,Median_MAPQ,Reads_Seen,"
            "Total_Bases,Softclip_Bases,Softclip_%,Windows_Created\n"
            "chr1,60.00,60,30,30000,0,0.00000,3\n"
            "GENOME,60.00,60,30,30000,0,0.00000,3\n"
        )
    _write_csv(os.path.join(folder, "window_stats.csv"), window_content)
    _write_csv(os.path.join(folder, "summary_stats.csv"), summary_content)


# ── smart_bp_formatter ────────────────────────────────────────────────────────

class TestSmartBpFormatter:
    """Returns (divisor, unit) for bp/kb/Mb based on the largest coordinate value."""

    def test_below_1000_returns_bp(self):
        divisor, unit = smart_bp_formatter(999)
        assert divisor == 1
        assert unit == "bp"

    def test_exactly_1000_returns_kb(self):
        divisor, unit = smart_bp_formatter(1_000)
        assert divisor == 1_000
        assert unit == "kb"

    def test_between_1000_and_1m_returns_kb(self):
        divisor, unit = smart_bp_formatter(500_000)
        assert divisor == 1_000
        assert unit == "kb"

    def test_exactly_1m_returns_mb(self):
        divisor, unit = smart_bp_formatter(1_000_000)
        assert divisor == 1_000_000
        assert unit == "Mb"

    def test_above_1m_returns_mb(self):
        divisor, unit = smart_bp_formatter(250_000_000)
        assert divisor == 1_000_000
        assert unit == "Mb"

    def test_small_contig_stays_as_bp(self):
        divisor, unit = smart_bp_formatter(1)
        assert divisor == 1
        assert unit == "bp"


# ── _sanitise_filename ────────────────────────────────────────────────────────

class TestSanitiseFilename:
    """Filesystem-unsafe characters replaced with underscore; safe chars preserved."""

    def test_clean_name_unchanged(self):
        assert _sanitise_filename("chr1") == "chr1"

    def test_colon_replaced(self):
        assert _sanitise_filename("chr1:100-200") == "chr1_100-200"

    def test_spaces_replaced(self):
        assert _sanitise_filename("my contig") == "my_contig"

    def test_dot_and_hyphen_preserved(self):
        assert _sanitise_filename("scaffold.1-hap1") == "scaffold.1-hap1"

    def test_special_chars_replaced(self):
        assert _sanitise_filename("ctg|1|foo") == "ctg_1_foo"


# ── _check_columns ────────────────────────────────────────────────────────────

class TestCheckColumns:
    """Missing required columns exit with a message naming the missing column(s)."""

    def test_all_columns_present_no_exit(self):
        df = pd.DataFrame(columns=["A", "B", "C"])
        _check_columns(df, {"A", "B"}, "test.csv")   # must not raise

    def test_missing_one_column_exits(self):
        df = pd.DataFrame(columns=["A", "B"])
        with pytest.raises(SystemExit):
            _check_columns(df, {"A", "B", "C"}, "test.csv")

    def test_missing_column_name_in_error_message(self):
        df = pd.DataFrame(columns=["A"])
        with pytest.raises(SystemExit, match="missing_col"):
            _check_columns(df, {"A", "missing_col"}, "test.csv")

    def test_multiple_missing_columns_listed(self):
        df = pd.DataFrame(columns=["A"])
        try:
            _check_columns(df, {"A", "B", "C"}, "test.csv")
        except SystemExit as e:
            msg = str(e)
            assert "B" in msg or "C" in msg

    def test_extra_columns_in_df_are_ignored(self):
        # DataFrame with more columns than required must not exit
        df = pd.DataFrame(columns=["A", "B", "C", "D", "E"])
        _check_columns(df, {"A", "B"}, "test.csv")   # must not raise


# ── make_figures_folder ───────────────────────────────────────────────────────

class TestMakeFiguresFolder:
    """Creates <base>/figures/ and is idempotent on repeated calls."""

    def test_creates_figures_subfolder(self, tmp_path):
        result = make_figures_folder(str(tmp_path))
        assert os.path.isdir(result)

    def test_returns_correct_path(self, tmp_path):
        result = make_figures_folder(str(tmp_path))
        assert result == os.path.join(str(tmp_path), "figures")

    def test_idempotent_no_error_on_second_call(self, tmp_path):
        make_figures_folder(str(tmp_path))
        make_figures_folder(str(tmp_path))   # must not raise (exist_ok=True)


# ── check_input ───────────────────────────────────────────────────────────────

class TestCheckInput:
    """Exits if the input folder or required CSV files are missing."""

    def test_missing_folder_exits(self):
        with pytest.raises(SystemExit, match="Input folder not found"):
            check_input("/nonexistent/path/xyz", "mapq_softclip")

    def test_missing_window_stats_exits(self, tmp_path):
        # Only summary present, window_stats.csv absent
        (tmp_path / "summary_stats.csv").touch()
        with pytest.raises(SystemExit, match="window_stats.csv"):
            check_input(str(tmp_path), "mapq_softclip")

    def test_missing_summary_stats_exits(self, tmp_path):
        # Only window present, summary_stats.csv absent
        (tmp_path / "window_stats.csv").touch()
        with pytest.raises(SystemExit, match="summary_stats.csv"):
            check_input(str(tmp_path), "mapq_softclip")

    def test_both_files_present_no_exit(self, tmp_path):
        (tmp_path / "window_stats.csv").touch()
        (tmp_path / "summary_stats.csv").touch()
        check_input(str(tmp_path), "mapq_softclip")   # must not raise

    def test_hese_missing_file_exits(self, tmp_path):
        (tmp_path / "haplotype_summary.csv").touch()
        (tmp_path / "unitig_labels.csv").touch()
        with pytest.raises(SystemExit, match="haplotig_labels.csv"):
            check_input(str(tmp_path), "hese")

    def test_hese_all_files_present_no_exit(self, tmp_path):
        (tmp_path / "haplotig_labels.csv").touch()
        (tmp_path / "haplotype_summary.csv").touch()
        (tmp_path / "unitig_labels.csv").touch()
        check_input(str(tmp_path), "hese")   # must not raise


# ── _build_genome_xaxis ───────────────────────────────────────────────────────

class TestBuildGenomeXaxis:
    """Assigns integer x-axis indices for equal-spaced genome-wide plots."""

    def _win_df(self, chrom, ends):
        return pd.DataFrame({"Chromosome": chrom, "End": ends})

    def test_single_chrom_total_is_one(self):
        df = self._win_df("chr1", [1000, 2000, 3000])
        total, *_ = _build_genome_xaxis(df, ["chr1"])
        assert total == 1

    def test_single_chrom_tick_at_index_zero(self):
        df = self._win_df("chr1", [1000, 2000, 3000])
        _, _, ticks, _, _ = _build_genome_xaxis(df, ["chr1"])
        assert ticks == [0]

    def test_single_chrom_label(self):
        df = self._win_df("chr1", [1000, 2000, 3000])
        _, _, _, labels, _ = _build_genome_xaxis(df, ["chr1"])
        assert labels == ["chr1"]

    def test_single_chrom_index_is_zero(self):
        df = self._win_df("chr1", [1000, 2000, 3000])
        _, _, _, _, chrom_x = _build_genome_xaxis(df, ["chr1"])
        assert chrom_x["chr1"][1] == 0

    def test_two_chroms_total_is_two(self):
        df = pd.concat([
            self._win_df("chr1", [1000, 2000, 3000]),
            self._win_df("chr2", [1000, 2000]),
        ])
        total, *_ = _build_genome_xaxis(df, ["chr1", "chr2"])
        assert total == 2

    def test_two_chroms_second_index_is_one(self):
        df = pd.concat([
            self._win_df("chr1", [1000, 2000, 3000]),
            self._win_df("chr2", [1000, 2000]),
        ])
        _, _, _, _, chrom_x = _build_genome_xaxis(df, ["chr1", "chr2"])
        assert chrom_x["chr2"][1] == 1

    def test_two_chroms_ticks_and_labels_ordered(self):
        df = pd.concat([
            self._win_df("chr1", [1000, 2000, 3000]),
            self._win_df("chr2", [1000, 2000]),
        ])
        _, _, ticks, labels, _ = _build_genome_xaxis(df, ["chr1", "chr2"])
        assert labels == ["chr1", "chr2"]
        assert ticks[0] == 0
        assert ticks[1] == 1

    def test_chrom_absent_from_df_is_skipped(self):
        df = self._win_df("chr1", [1000, 2000, 3000])
        _, _, _, _, chrom_x = _build_genome_xaxis(df, ["chr1", "chrX"])
        assert "chrX" not in chrom_x
        assert "chr1" in chrom_x

    def test_dividers_count_matches_chrom_count(self):
        df = pd.concat([
            self._win_df("chr1", [1000, 2000, 3000]),
            self._win_df("chr2", [1000, 2000]),
        ])
        _, dividers, *_ = _build_genome_xaxis(df, ["chr1", "chr2"])
        assert len(dividers) == 2


# ══════════════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS
# ══════════════════════════════════════════════════════════════════════════════
#
# Tests verify that the expected output files are created at the correct paths.
# Pixel content is not tested - the behavioural contract is "file created, no crash".
# Two levels:
#   - hand-crafted CSVs - run_mapq_softclip (fast, no BAM needed)
#   - synthetic BAM - run_analysis - run_mapq_softclip (full pipeline)

@pytest.fixture(scope="module")
def csv_folder(tmp_path_factory):
    folder = str(tmp_path_factory.mktemp("vis_input"))
    _make_minimal_csvs(folder)
    return folder


class TestFigureFileCreation:
    """run_mapq_softclip always creates all expected figure files."""

    def test_genome_figure_created(self, tmp_path_factory):
        out = str(tmp_path_factory.mktemp("vis_genome"))
        _make_minimal_csvs(out)
        run_mapq_softclip(out, rolling_target_bp=1000)
        assert os.path.exists(os.path.join(out, "figures", "genome_summary.png"))

    def test_contigs_raw_figure_created(self, tmp_path_factory):
        out = str(tmp_path_factory.mktemp("vis_contigs"))
        _make_minimal_csvs(out)
        run_mapq_softclip(out, rolling_target_bp=1000)
        assert os.path.exists(
            os.path.join(out, "figures", "contigs", "raw", "chr1_per_window_raw.png")
        )

    def test_contigs_rolling_figure_created(self, tmp_path_factory):
        out = str(tmp_path_factory.mktemp("vis_rolling"))
        _make_minimal_csvs(out)
        run_mapq_softclip(out, rolling_target_bp=1000)
        assert os.path.exists(
            os.path.join(out, "figures", "contigs", "rolling", "chr1_per_window_rolling.png")
        )


class TestFigureEdgeCases:
    """Edge inputs that must not crash, and rolling window capping and warning behaviour."""

    def test_single_window_contig_does_not_crash(self, tmp_path_factory):
        # A contig with only one window cannot compute a step size from two rows.
        # The code falls back to "N/A (single window)" - it must not crash.
        out = str(tmp_path_factory.mktemp("vis_single_win"))
        _make_minimal_csvs(
            out,
            window_content=(
                "Chromosome,Start,End,Mean_MAPQ,Median_MAPQ,Read_Count,"
                "Total_Bases,Softclip_Bases,Softclip_%,Flag\n"
                "chr1,0,800,60.00,60,5,4000,0,0.00000,\n"
            ),
            summary_content=(
                "Chromosome,Mean_MAPQ,Median_MAPQ,Reads_Seen,"
                "Total_Bases,Softclip_Bases,Softclip_%,Windows_Created\n"
                "chr1,60.00,60,5,4000,0,0.00000,1\n"
                "GENOME,60.00,60,5,4000,0,0.00000,1\n"
            ),
        )
        run_mapq_softclip(out, rolling_target_bp=1000)
        assert os.path.exists(
            os.path.join(out, "figures", "contigs", "raw", "chr1_per_window_raw.png")
        )

    def test_float_rolling_value_does_not_crash(self, tmp_path_factory):
        # rolling_target_bp can be a float (e.g. 1.5 kb = 1500.5 bp)
        out = str(tmp_path_factory.mktemp("vis_float_rolling"))
        _make_minimal_csvs(out)
        run_mapq_softclip(out, rolling_target_bp=1500.5)
        assert os.path.exists(
            os.path.join(out, "figures", "contigs", "rolling", "chr1_per_window_rolling.png")
        )

    def test_oversized_rolling_capped_and_warns(self, tmp_path_factory, capsys):
        # rolling_target_bp larger than the contig must be silently capped and warn once
        out = str(tmp_path_factory.mktemp("vis_large_rolling"))
        _make_minimal_csvs(out)
        run_mapq_softclip(out, rolling_target_bp=999_999_999)
        assert os.path.exists(
            os.path.join(out, "figures", "contigs", "rolling", "chr1_per_window_rolling.png")
        )
        assert "[WARN]" in capsys.readouterr().out

    def test_no_warning_when_rolling_fits_all_contigs(self, tmp_path_factory, capsys):
        # rolling smaller than every contig, no [WARN] should appear
        out = str(tmp_path_factory.mktemp("vis_small_rolling"))
        _make_minimal_csvs(out)
        run_mapq_softclip(out, rolling_target_bp=500)   # step=1000bp, 3 windows, fits
        assert "[WARN]" not in capsys.readouterr().out

    def test_warning_message_contains_rolling_value(self, tmp_path_factory, capsys):
        # the [WARN] message must include the rolling value so the user knows what triggered it
        out = str(tmp_path_factory.mktemp("vis_warn_content"))
        _make_minimal_csvs(out)
        run_mapq_softclip(out, rolling_target_bp=999_999_999)
        out_text = capsys.readouterr().out
        assert "[WARN]" in out_text
        assert "--rolling" in out_text   # argument name appears so user knows what to adjust
        assert "capped" in out_text      # capping behaviour is mentioned explicitly

    def test_partial_capping_warns_when_only_some_contigs_exceeded(self, tmp_path_factory, capsys):
        # chr1 has 5 windows (step 1000bp), chr2 has 2 windows (step 1000bp)
        # rolling=3500bp: uncapped_n=4, fits chr1 (5 windows) but caps chr2 (2 windows)
        out = str(tmp_path_factory.mktemp("vis_partial_cap"))
        _make_minimal_csvs(
            out,
            window_content=(
                "Chromosome,Start,End,Mean_MAPQ,Median_MAPQ,Read_Count,"
                "Total_Bases,Softclip_Bases,Softclip_%,Flag\n"
                "chr1,0,1000,60.00,60,10,10000,0,0.00000,\n"
                "chr1,1000,2000,60.00,60,10,10000,0,0.00000,\n"
                "chr1,2000,3000,60.00,60,10,10000,0,0.00000,\n"
                "chr1,3000,4000,60.00,60,10,10000,0,0.00000,\n"
                "chr1,4000,5000,60.00,60,10,10000,0,0.00000,\n"
                "chr2,0,1000,60.00,60,10,10000,0,0.00000,\n"
                "chr2,1000,2000,60.00,60,10,10000,0,0.00000,\n"
            ),
            summary_content=(
                "Chromosome,Mean_MAPQ,Median_MAPQ,Reads_Seen,"
                "Total_Bases,Softclip_Bases,Softclip_%,Windows_Created\n"
                "chr1,60.00,60,50,50000,0,0.00000,5\n"
                "chr2,60.00,60,20,20000,0,0.00000,2\n"
                "GENOME,60.00,60,70,70000,0,0.00000,7\n"
            ),
        )
        run_mapq_softclip(out, rolling_target_bp=3500)
        assert "[WARN]" in capsys.readouterr().out

    def test_contig_name_with_special_chars_sanitised(self, tmp_path_factory):
        # Contig names like "chr1:1-5000" must be sanitised before use as filenames.
        # The expected filename is the sanitised version.
        out = str(tmp_path_factory.mktemp("vis_special_name"))
        _make_minimal_csvs(
            out,
            window_content=(
                "Chromosome,Start,End,Mean_MAPQ,Median_MAPQ,Read_Count,"
                "Total_Bases,Softclip_Bases,Softclip_%,Flag\n"
                "chr1:seg1,0,1000,60.00,60,10,10000,0,0.00000,\n"
                "chr1:seg1,1000,2000,60.00,60,10,10000,0,0.00000,\n"
            ),
            summary_content=(
                "Chromosome,Mean_MAPQ,Median_MAPQ,Reads_Seen,"
                "Total_Bases,Softclip_Bases,Softclip_%,Windows_Created\n"
                "chr1:seg1,60.00,60,20,20000,0,0.00000,2\n"
                "GENOME,60.00,60,20,20000,0,0.00000,2\n"
            ),
        )
        run_mapq_softclip(out, rolling_target_bp=1000)
        # colon must be replaced by underscore in the output filename
        assert os.path.exists(
            os.path.join(out, "figures", "contigs", "raw", "chr1_seg1_per_window_raw.png")
        )


# ── full pipeline: mapq_softclip to visualise ─────────────────────────────────
#
# These tests use real run_analysis output rather than hand-crafted CSVs.
# The fixture creates a two-contig BAM, runs mapq_softclip analysis, then
# runs visualise in "all" mode. Tests assert expected figure files exist.
#
# chr1 (3000 bp): 10 reads × 100M, MAPQ 60, one read per 100 bp
# chr2 (2000 bp):  8 reads × 100M, MAPQ 40
# windows: 1 kb / 1 kb step
# chr1: [0,1000) [1000,2000) [2000,3000)
# chr2: [0,1000) [1000,2000)

@pytest.fixture(scope="module")
def pipeline_output(tmp_path_factory):
    """BAM to run_analysis to run_mapq_softclip. Returns the output folder path."""
    base = tmp_path_factory.mktemp("vis_pipeline")
    bam_path = str(base / "test.bam")

    sq = [{"SN": "chr1", "LN": 3000}, {"SN": "chr2", "LN": 2000}]
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": sq,
    })

    reads = []
    for ref_id, n, mapq in [(0, 10, 60), (1, 8, 40)]:
        for i in range(n):
            r = pysam.AlignedSegment(header)
            r.query_name      = f"r{ref_id}_{i}"
            r.query_sequence  = "A" * 100
            r.flag            = 0
            r.reference_id    = ref_id
            r.reference_start = i * 100
            r.mapping_quality = mapq
            r.cigar           = [(0, 100)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            reads.append(r)

    tmp = bam_path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=header) as bam:
        for r in sorted(reads, key=lambda x: (x.reference_id, x.reference_start)):
            bam.write(r)
    pysam.sort("-o", bam_path, tmp)
    os.remove(tmp)
    pysam.index(bam_path)

    contigs_folder = str(base / "contigs")
    os.makedirs(contigs_folder)
    window_file  = str(base / "window_stats.csv")
    summary_file = str(base / "summary_stats.csv")
    run_analysis(bam_path, bam_path + ".bai", 1.0, 1.0, 1,
                 window_file, summary_file, contigs_folder)

    run_mapq_softclip(str(base), rolling_target_bp=1000)
    return str(base)


class TestVisualiseEndToEnd:
    """Figures generated from real mapq_softclip output (full BAM to analysis to visualise pipeline)."""

    def test_genome_figure_created(self, pipeline_output):
        assert os.path.exists(os.path.join(pipeline_output, "figures", "genome_summary.png"))

    def test_chr1_raw_figure_created(self, pipeline_output):
        assert os.path.exists(
            os.path.join(pipeline_output, "figures", "contigs", "raw", "chr1_per_window_raw.png")
        )

    def test_chr2_raw_figure_created(self, pipeline_output):
        # Second contig with different MAPQ must also produce a figure
        assert os.path.exists(
            os.path.join(pipeline_output, "figures", "contigs", "raw", "chr2_per_window_raw.png")
        )

    def test_chr1_rolling_figure_created(self, pipeline_output):
        assert os.path.exists(
            os.path.join(pipeline_output, "figures", "contigs", "rolling", "chr1_per_window_rolling.png")
        )

    def test_chr2_rolling_figure_created(self, pipeline_output):
        assert os.path.exists(
            os.path.join(pipeline_output, "figures", "contigs", "rolling", "chr2_per_window_rolling.png")
        )

    def test_no_crash_on_mixed_mapq_contigs(self, pipeline_output):
        # chr1 (MAPQ 60) and chr2 (MAPQ 40) in the same genome view must not crash
        # (verified implicitly by genome_summary.png existing, but made explicit here)
        assert os.path.exists(os.path.join(pipeline_output, "figures", "genome_summary.png"))


# ══════════════════════════════════════════════════════════════════════════════
# UNIT TESTS - hese visualise utilities
# ══════════════════════════════════════════════════════════════════════════════

# ── _to_float_or_nan ──────────────────────────────────────────────────────────

class TestToFloatOrNan:
    """Converts numeric values to float and returns nan for anything non-numeric."""

    def test_numeric_string(self):
        assert _to_float_or_nan("3.14") == pytest.approx(3.14)

    def test_float_passthrough(self):
        assert _to_float_or_nan(2.5) == pytest.approx(2.5)

    def test_integer_passthrough(self):
        assert _to_float_or_nan(5) == pytest.approx(5.0)

    def test_zero(self):
        assert _to_float_or_nan(0) == pytest.approx(0.0)

    def test_na_string_returns_nan(self):
        assert math.isnan(_to_float_or_nan("NA"))

    def test_none_returns_nan(self):
        assert math.isnan(_to_float_or_nan(None))

    def test_garbage_string_returns_nan(self):
        assert math.isnan(_to_float_or_nan("not_a_number"))

    def test_empty_string_returns_nan(self):
        assert math.isnan(_to_float_or_nan(""))


# ── _safe_nanmax ──────────────────────────────────────────────────────────────

class TestSafeNanmax:
    """Returns max of finite values, or the given default when none are finite."""

    def test_all_finite(self):
        assert _safe_nanmax([1.0, 2.0, 3.0], default=0.0) == pytest.approx(3.0)

    def test_with_nan_ignored(self):
        assert _safe_nanmax([1.0, float("nan"), 3.0], default=0.0) == pytest.approx(3.0)

    def test_all_nan_returns_default(self):
        assert _safe_nanmax([float("nan"), float("nan")], default=99.0) == pytest.approx(99.0)

    def test_empty_list_returns_default(self):
        assert _safe_nanmax([], default=42.0) == pytest.approx(42.0)

    def test_single_finite_value(self):
        assert _safe_nanmax([7.0], default=0.0) == pytest.approx(7.0)

    def test_inf_excluded(self):
        assert _safe_nanmax([1.0, float("inf"), 2.0], default=0.0) == pytest.approx(2.0)


# ── _read_min_unitigs ─────────────────────────────────────────────────────────

class TestReadMinUnitigs:
    """Parses min_unitigs from run_log.txt with a fallback to a default."""

    def test_reads_value_from_log(self, tmp_path):
        (tmp_path / "run_log.txt").write_text(
            "Min unitigs (-mu)      : 5\n"
        )
        assert _read_min_unitigs(str(tmp_path)) == 5

    def test_missing_file_returns_default(self, tmp_path):
        assert _read_min_unitigs(str(tmp_path), default=3) == 3

    def test_custom_default_returned_when_no_file(self, tmp_path):
        assert _read_min_unitigs(str(tmp_path), default=7) == 7

    def test_malformed_value_returns_default(self, tmp_path):
        (tmp_path / "run_log.txt").write_text(
            "Min unitigs (-mu)      : not_a_number\n"
        )
        assert _read_min_unitigs(str(tmp_path), default=3) == 3

    def test_ignores_unrelated_lines(self, tmp_path):
        (tmp_path / "run_log.txt").write_text(
            "Started : 2026-01-01 00:00:00\n"
            "Min unitigs (-mu)      : 4\n"
            "Hap frac (-hf)         : 0.60\n"
        )
        assert _read_min_unitigs(str(tmp_path)) == 4


# ══════════════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS - hese figures
# ══════════════════════════════════════════════════════════════════════════════
#
# Figure functions are called directly with hand-crafted DataFrames.
# The behavioural contract is "file created, no crash" - pixel content is not tested.
# Edge cases (empty inputs, NA values) verify graceful handling.

def _make_hese_dfs():
    """Return (df_uni, df_hap, df_sum) with minimal valid hese output data."""
    df_uni = pd.DataFrame({
        "unitig":  ["utg001", "utg002", "utg003"],
        "p1_norm": [1e-4, 2e-4, 1e-4],
        "p2_norm": [2e-4, 1e-4, 1e-4],
        "label":   ["P1", "P2", "amb"],
    })
    df_hap = pd.DataFrame({
        "hap_id":    ["h1tg000001l", "h1tg000002l", "h2tg000001l", "h2tg000002l"],
        "hap_type":  ["hap1", "hap1", "hap2", "hap2"],
        "n_unitigs": [10, 10, 10, 10],
        "n_p1":      [8, 7, 1, 2],
        "n_p2":      [1, 2, 8, 7],
        "n_amb":     [1, 1, 1, 1],
        "label":     ["P1", "P1", "P2", "P2"],
    })
    # haplotype_summary has four rows: hap1×P1, hap1×P2, hap2×P1, hap2×P2
    df_sum = pd.DataFrame({
        "hap":               ["hap1", "hap1", "hap2", "hap2"],
        "hap_global_parent": ["P1",   "P1",   "P2",   "P2"  ],
        "assigned_parent":   ["P1",   "P2",   "P1",   "P2"  ],
        "n_haplotigs_total": [2, 2, 2, 2],
        "n_haplotigs_p1":    [2, 2, 0, 0],
        "n_haplotigs_p2":    [0, 0, 2, 2],
        "n_haplotigs_amb":   [0, 0, 0, 0],
        "hamming_struct_%":  [0.0, 100.0, 100.0, 0.0],
        "switch_struct_%":   [0.0,   0.0,   0.0, 0.0],
    })
    return df_uni, df_hap, df_sum


def _write_hese_csvs(folder, df_uni, df_hap, df_sum, df_truth_eval=None):
    """Write hese output CSVs to folder for use with run_hese."""
    df_uni.to_csv(os.path.join(folder, "unitig_labels.csv"), index=False)
    df_hap.to_csv(os.path.join(folder, "haplotig_labels.csv"), index=False)
    df_sum.to_csv(os.path.join(folder, "haplotype_summary.csv"), index=False)
    if df_truth_eval is not None:
        df_truth_eval.to_csv(os.path.join(folder, "truth_eval_summary.csv"), index=False)


# ── figure file creation ──────────────────────────────────────────────────────

class TestFigHeseFileCreation:
    """Each hese figure function saves its expected PNG to the figures folder."""

    def test_label_distributions_created(self, tmp_path):
        df_uni, df_hap, _ = _make_hese_dfs()
        fig_hese_label_distributions(df_uni, df_hap, str(tmp_path))
        assert os.path.exists(os.path.join(str(tmp_path), "hese_label_distributions.png"))

    def test_label_balance_created(self, tmp_path):
        _, df_hap, _ = _make_hese_dfs()
        fig_hese_label_balance(df_hap, str(tmp_path))
        assert os.path.exists(os.path.join(str(tmp_path), "hese_label_balance.png"))

    def test_signal_penetration_created(self, tmp_path):
        _, df_hap, _ = _make_hese_dfs()
        fig_hese_signal_penetration(df_hap, min_unitigs=3, figures_folder=str(tmp_path))
        assert os.path.exists(os.path.join(str(tmp_path), "hese_signal_penetration.png"))

    def test_signal_depth_created(self, tmp_path):
        _, df_hap, _ = _make_hese_dfs()
        fig_hese_signal_depth(df_hap, str(tmp_path))
        assert os.path.exists(os.path.join(str(tmp_path), "hese_signal_depth.png"))

    def test_phasing_errors_no_truth_created(self, tmp_path):
        _, _, df_sum = _make_hese_dfs()
        fig_hese_phasing_errors(df_sum, df_truth_eval=None, figures_folder=str(tmp_path))
        assert os.path.exists(os.path.join(str(tmp_path), "hese_phasing_errors.png"))

    def test_phasing_errors_with_truth_created(self, tmp_path):
        _, _, df_sum = _make_hese_dfs()
        df_truth_eval = pd.DataFrame({
            "hap1_label": ["P1"], "hap2_label": ["P2"],
            "overall_wrong": [0], "overall_total": [4],
            "overall_hamming_%": [0.0],
            "overall_switches": [0], "overall_boundaries": [3],
            "overall_switch_%": [0.0],
        })
        fig_hese_phasing_errors(df_sum, df_truth_eval, str(tmp_path))
        assert os.path.exists(os.path.join(str(tmp_path), "hese_phasing_errors.png"))

    def test_phasing_errors_shared_parent_no_crash(self, tmp_path):
        # Both hap1 and hap2 claim P1, title turns red, must not crash
        df_sum = pd.DataFrame({
            "hap":               ["hap1", "hap1", "hap2", "hap2"],
            "hap_global_parent": ["P1",   "P1",   "P1",   "P1"  ],
            "assigned_parent":   ["P1",   "P2",   "P1",   "P2"  ],
            "n_haplotigs_total": [2, 2, 2, 2],
            "n_haplotigs_p1":    [2, 2, 2, 2],
            "n_haplotigs_p2":    [0, 0, 0, 0],
            "n_haplotigs_amb":   [0, 0, 0, 0],
            "hamming_struct_%":  [0.0, 100.0, 0.0, 100.0],
            "switch_struct_%":   [0.0,   0.0, 0.0,   0.0],
        })
        fig_hese_phasing_errors(df_sum, df_truth_eval=None, figures_folder=str(tmp_path))
        assert os.path.exists(os.path.join(str(tmp_path), "hese_phasing_errors.png"))


# ── edge cases ────────────────────────────────────────────────────────────────

class TestFigHeseEdgeCases:
    """Empty inputs skip with a warning; NA metric values do not crash."""

    def test_label_balance_empty_df_warns_and_skips(self, tmp_path, capsys):
        df_hap = pd.DataFrame(columns=["hap_id", "hap_type", "n_unitigs",
                                        "n_p1", "n_p2", "n_amb", "label"])
        fig_hese_label_balance(df_hap, str(tmp_path))
        assert "[WARN]" in capsys.readouterr().out
        assert not os.path.exists(os.path.join(str(tmp_path), "hese_label_balance.png"))

    def test_signal_penetration_empty_df_warns_and_skips(self, tmp_path, capsys):
        df_hap = pd.DataFrame(columns=["hap_id", "hap_type", "n_unitigs",
                                        "n_p1", "n_p2", "n_amb", "label"])
        fig_hese_signal_penetration(df_hap, min_unitigs=3, figures_folder=str(tmp_path))
        assert "[WARN]" in capsys.readouterr().out
        assert not os.path.exists(os.path.join(str(tmp_path), "hese_signal_penetration.png"))

    def test_signal_depth_all_amb_warns_and_skips(self, tmp_path, capsys):
        # No P1/P2 labeled haplotigs, signal depth must warn and skip
        df_hap = pd.DataFrame({
            "hap_id": ["h1", "h2"], "hap_type": ["hap1", "hap2"],
            "n_unitigs": [5, 5], "n_p1": [0, 0], "n_p2": [0, 0], "n_amb": [5, 5],
            "label": ["amb", "amb"],
        })
        fig_hese_signal_depth(df_hap, str(tmp_path))
        assert "[WARN]" in capsys.readouterr().out
        assert not os.path.exists(os.path.join(str(tmp_path), "hese_signal_depth.png"))

    def test_phasing_errors_na_metrics_no_crash(self, tmp_path):
        # hamming_struct_% and switch_struct_% are "NA", must not crash
        df_sum = pd.DataFrame({
            "hap":               ["hap1", "hap1", "hap2", "hap2"],
            "hap_global_parent": ["P1",   "P1",   "P2",   "P2"  ],
            "assigned_parent":   ["P1",   "P2",   "P1",   "P2"  ],
            "n_haplotigs_total": [1, 1, 1, 1],
            "n_haplotigs_p1":    [1, 1, 0, 0],
            "n_haplotigs_p2":    [0, 0, 1, 1],
            "n_haplotigs_amb":   [0, 0, 0, 0],
            "hamming_struct_%":  ["NA", "NA", "NA", "NA"],
            "switch_struct_%":   ["NA", "NA", "NA", "NA"],
        })
        fig_hese_phasing_errors(df_sum, df_truth_eval=None, figures_folder=str(tmp_path))
        assert os.path.exists(os.path.join(str(tmp_path), "hese_phasing_errors.png"))

    def test_phasing_errors_truth_na_switch_no_crash(self, tmp_path):
        # truth_eval_summary has overall_switch_% = "NA" (no boundaries), must not crash
        _, _, df_sum = _make_hese_dfs()
        df_truth_eval = pd.DataFrame({
            "hap1_label": ["P1"], "hap2_label": ["P2"],
            "overall_wrong": [0], "overall_total": [0],
            "overall_hamming_%": ["NA"],
            "overall_switches": [0], "overall_boundaries": [0],
            "overall_switch_%": ["NA"],
        })
        fig_hese_phasing_errors(df_sum, df_truth_eval, str(tmp_path))
        assert os.path.exists(os.path.join(str(tmp_path), "hese_phasing_errors.png"))

    def test_phasing_errors_amb_global_parent_no_crash(self, tmp_path):
        # hap_global_parent = "amb" means no best row matches, hap skipped, blank figure
        df_sum = pd.DataFrame({
            "hap":               ["hap1", "hap1", "hap2", "hap2"],
            "hap_global_parent": ["amb",  "amb",  "amb",  "amb" ],
            "assigned_parent":   ["P1",   "P2",   "P1",   "P2"  ],
            "n_haplotigs_total": [1, 1, 1, 1],
            "n_haplotigs_p1":    [0, 0, 0, 0],
            "n_haplotigs_p2":    [0, 0, 0, 0],
            "n_haplotigs_amb":   [1, 1, 1, 1],
            "hamming_struct_%":  ["NA", "NA", "NA", "NA"],
            "switch_struct_%":   ["NA", "NA", "NA", "NA"],
        })
        fig_hese_phasing_errors(df_sum, df_truth_eval=None, figures_folder=str(tmp_path))
        assert os.path.exists(os.path.join(str(tmp_path), "hese_phasing_errors.png"))


# ── run_hese end-to-end ───────────────────────────────────────────────────────

class TestRunHese:
    """run_hese generates all expected figures from a valid hese output folder."""

    def test_all_core_figures_created(self, tmp_path):
        df_uni, df_hap, df_sum = _make_hese_dfs()
        _write_hese_csvs(str(tmp_path), df_uni, df_hap, df_sum)
        run_hese(str(tmp_path))
        figs = tmp_path / "figures"
        assert (figs / "hese_label_distributions.png").exists()
        assert (figs / "hese_label_balance.png").exists()
        assert (figs / "hese_signal_penetration.png").exists()
        assert (figs / "hese_signal_depth.png").exists()
        assert (figs / "hese_phasing_errors.png").exists()

    def test_with_truth_eval_no_crash(self, tmp_path):
        df_uni, df_hap, df_sum = _make_hese_dfs()
        df_truth_eval = pd.DataFrame({
            "hap1_label": ["P1"], "hap2_label": ["P2"],
            "overall_wrong": [0], "overall_total": [4],
            "overall_hamming_%": [0.0],
            "overall_switches": [0], "overall_boundaries": [3],
            "overall_switch_%": [0.0],
        })
        _write_hese_csvs(str(tmp_path), df_uni, df_hap, df_sum, df_truth_eval)
        run_hese(str(tmp_path))
        assert (tmp_path / "figures" / "hese_phasing_errors.png").exists()

    def test_with_truth_na_values_no_crash(self, tmp_path):
        # truth_eval_summary.csv with NA switch/hamming must not crash run_hese
        df_uni, df_hap, df_sum = _make_hese_dfs()
        df_truth_eval = pd.DataFrame({
            "hap1_label": ["P1"], "hap2_label": ["P2"],
            "overall_wrong": [0], "overall_total": [0],
            "overall_hamming_%": ["NA"],
            "overall_switches": [0], "overall_boundaries": [0],
            "overall_switch_%": ["NA"],
        })
        _write_hese_csvs(str(tmp_path), df_uni, df_hap, df_sum, df_truth_eval)
        run_hese(str(tmp_path))
        assert (tmp_path / "figures" / "hese_phasing_errors.png").exists()

    def test_missing_required_file_exits(self, tmp_path):
        # haplotig_labels.csv missing, check_input should catch this before run_hese
        df_uni, _, df_sum = _make_hese_dfs()
        df_uni.to_csv(os.path.join(str(tmp_path), "unitig_labels.csv"), index=False)
        df_sum.to_csv(os.path.join(str(tmp_path), "haplotype_summary.csv"), index=False)
        with pytest.raises(SystemExit):
            check_input(str(tmp_path), "hese")