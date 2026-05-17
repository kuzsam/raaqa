"""
Tests for the visualise module.

Organised as:
  UNIT TESTS        — isolated function tests using direct inputs (no files)
  INTEGRATION TESTS — run_mapq_softclip called with hand-crafted CSVs or a real
                      BAM pipeline, verifying figure files are created correctly
"""

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
    _check_columns,
    make_figures_folder,
    check_input,
    _build_genome_xaxis,
    run_mapq_softclip
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
# Pixel content is not tested — the behavioural contract is "file created, no crash".
# Two levels:
#   - hand-crafted CSVs → run_mapq_softclip (fast, no BAM needed)
#   - synthetic BAM → run_analysis → run_mapq_softclip (full pipeline)

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
        # The code falls back to "N/A (single window)" — it must not crash.
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
        # rolling_target_bp can be a float (e.g. 1.5 kb → 1500.5 bp)
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
        # rolling smaller than every contig — no [WARN] should appear
        out = str(tmp_path_factory.mktemp("vis_small_rolling"))
        _make_minimal_csvs(out)
        run_mapq_softclip(out, rolling_target_bp=500)   # step=1000bp, 3 windows → fits
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
        # rolling=3500bp → uncapped_n=4, fits chr1 (5 windows) but caps chr2 (2 windows)
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


# ── full pipeline: mapq_softclip → visualise ──────────────────────────────────
#
# These tests use real run_analysis output rather than hand-crafted CSVs.
# The fixture creates a two-contig BAM, runs mapq_softclip analysis, then
# runs visualise in "all" mode. Tests assert expected figure files exist.
#
# chr1 (3000 bp): 10 reads × 100M, MAPQ 60, one read per 100 bp
# chr2 (2000 bp):  8 reads × 100M, MAPQ 40
# windows: 1 kb / 1 kb step
# → chr1: [0,1000) [1000,2000) [2000,3000)
# → chr2: [0,1000) [1000,2000)

@pytest.fixture(scope="module")
def pipeline_output(tmp_path_factory):
    """BAM → run_analysis → run_mapq_softclip('all'). Returns the output folder path."""
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
    """Figures generated from real mapq_softclip output (full BAM → analysis → visualise pipeline)."""

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