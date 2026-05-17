"""
Tests for the mapq_softclip module.

Organised as:
  UNIT TESTS        — isolated function tests using mocks and synthetic inputs
  INTEGRATION TESTS — full run_analysis pipeline against synthetic BAM files
"""

import argparse
import csv
import math
import os
from unittest.mock import MagicMock

import pysam
import pytest

from raaqa.mapq_softclip import (
    MAX_MAPQ,
    MAX_WINDOW_KB,
    LOW_COVERAGE_READ_THRESHOLD,
    LOW_DEPTH_THRESHOLD,
    compute_median_from_hist,
    compute_mean_from_hist,
    get_softclip_bases,
    _make_window,
    _create_windows,
    _flush_window_to_csv,
    _format_or_empty,
    _format_kb_value,
    _sanitise_filename,
    validate_inputs,
    _accumulate_read_into_windows,
    prepare_output_dirs,
    run_analysis,
)


# ══════════════════════════════════════════════════════════════════════════════
# UNIT TESTS
# ══════════════════════════════════════════════════════════════════════════════

# ── helpers ───────────────────────────────────────────────────────────────────


def _empty_hist():
    return [0] * (MAX_MAPQ + 1)


def _hist_at(mapq, count=1):
    h = _empty_hist()
    h[mapq] = count
    return h


def _hist_multi(pairs):
    h = _empty_hist()
    for mapq, count in pairs:
        h[mapq] = count
    return h


class _MockRead:
    def __init__(self, cigartuples, ref_start, ref_end=None, mapq=60):
        self.cigartuples = cigartuples
        self.reference_start = ref_start
        self.reference_end = ref_end
        self.mapping_quality = mapq


def _window_with(
    start, end, read_count=0, total_bases=0, softclip_bases=0, hist_pairs=None
):
    w = _make_window(start, end - start)
    w["read_count"] = read_count
    w["total_bases"] = total_bases
    w["softclip_bases"] = softclip_bases
    if hist_pairs:
        for mapq, count in hist_pairs:
            w["hist"][mapq] = count
    return w


def _flush_and_get_row(w, chrom="chr1"):
    writer = MagicMock()
    _flush_window_to_csv(w, chrom, writer)
    return writer.writerow.call_args[0][0]


def _make_args(bam, window=5.0, step=2.5, threads=1):
    return argparse.Namespace(bam=bam, window=window, step=step, threads=threads)


# ── compute_median_from_hist ──────────────────────────────────────────────────


class TestComputeMedianFromHist:
    """Histogram-based median: NaN on empty, correct value across boundary distributions."""

    def test_empty_returns_nan(self):
        assert math.isnan(compute_median_from_hist(_empty_hist()))

    def test_single_value_at_zero(self):
        assert compute_median_from_hist(_hist_at(0)) == 0

    def test_single_value_at_60(self):
        assert compute_median_from_hist(_hist_at(60)) == 60

    def test_single_value_at_max_mapq(self):
        assert compute_median_from_hist(_hist_at(MAX_MAPQ - 1)) == MAX_MAPQ - 1

    def test_equal_weight_returns_lower(self):
        # midpoint=1.0, cumulative at 40 reaches 1 >= 1.0 → 40
        h = _hist_multi([(40, 1), (60, 1)])
        assert compute_median_from_hist(h) == 40

    def test_skewed_weight_returns_higher(self):
        # total=4, midpoint=2.0 — cumulative at 40 is 1 < 2, at 60 is 4 >= 2 → 60
        h = _hist_multi([(40, 1), (60, 3)])
        assert compute_median_from_hist(h) == 60

    def test_heavy_lower_half(self):
        # total=4, midpoint=2.0 — cumulative at 20 is 3 >= 2 → 20
        h = _hist_multi([(20, 3), (80, 1)])
        assert compute_median_from_hist(h) == 20

    def test_large_equal_weight(self):
        h = _hist_multi([(30, 50), (70, 50)])
        assert compute_median_from_hist(h) == 30

    def test_multiple_bins_known_result(self):
        # total=4, midpoint=2.0 — at 10: cum=1<2, at 20: cum=3>=2 → 20
        h = _hist_multi([(10, 1), (20, 2), (30, 1)])
        assert compute_median_from_hist(h) == 20

    def test_many_reads_all_mapq_zero(self):
        # 100 bases all at MAPQ 0 — median must be 0, not NaN or any other value
        h = _hist_at(0, 100)
        assert compute_median_from_hist(h) == 0


# ── compute_mean_from_hist ────────────────────────────────────────────────────


class TestComputeMeanFromHist:
    """Histogram-based mean: NaN on empty, weighted average matches known results."""

    def test_empty_returns_nan(self):
        assert math.isnan(compute_mean_from_hist(_empty_hist()))

    def test_single_value(self):
        assert compute_mean_from_hist(_hist_at(40)) == pytest.approx(40.0)

    def test_equal_weight_two_values(self):
        h = _hist_multi([(20, 1), (60, 1)])
        assert compute_mean_from_hist(h) == pytest.approx(40.0)

    def test_known_weighted_mean(self):
        # (30*3 + 60*1) / 4 = 37.5
        h = _hist_multi([(30, 3), (60, 1)])
        assert compute_mean_from_hist(h) == pytest.approx(37.5)

    def test_all_at_zero(self):
        assert compute_mean_from_hist(_hist_at(0, 100)) == pytest.approx(0.0)

    def test_all_at_max_mapq(self):
        assert compute_mean_from_hist(_hist_at(MAX_MAPQ)) == pytest.approx(MAX_MAPQ)


# ── _format_or_empty ──────────────────────────────────────────────────────────


class TestFormatOrEmpty:
    """NaN → empty string (used for NO_COVERAGE fields); finite values formatted normally."""

    def test_nan_returns_empty(self):
        assert _format_or_empty(math.nan) == ""

    def test_nan_with_spec_returns_empty(self):
        assert _format_or_empty(math.nan, ".2f") == ""

    def test_zero_no_spec(self):
        assert _format_or_empty(0.0) == "0.0"

    def test_float_with_two_decimal_spec(self):
        assert _format_or_empty(42.1234, ".2f") == "42.12"

    def test_integer_value(self):
        assert _format_or_empty(60, "") == "60"

    def test_five_decimal_spec(self):
        assert _format_or_empty(1.23456789, ".5f") == "1.23457"


# ── _format_kb_value ──────────────────────────────────────────────────────────


class TestFormatKbValue:
    """Integer kb values use no decimal; fractional values use comma as separator."""

    def test_integer_float_no_decimal(self):
        assert _format_kb_value(5.0) == "5"

    def test_fractional_uses_comma(self):
        assert _format_kb_value(2.5) == "2,5"

    def test_large_integer(self):
        assert _format_kb_value(1000.0) == "1000"

    def test_small_fractional(self):
        assert _format_kb_value(0.1) == "0,1"

    def test_integer_10(self):
        assert _format_kb_value(10.0) == "10"


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
        assert _sanitise_filename("file.name-v1") == "file.name-v1"

    def test_special_chars_replaced(self):
        assert _sanitise_filename("contig!@#$") == "contig_"

    def test_consecutive_special_chars_collapsed(self):
        assert _sanitise_filename("contig!@#$name") == "contig_name"

    def test_alphanumeric_unchanged(self):
        assert _sanitise_filename("ABC123") == "ABC123"

    def test_underscore_preserved(self):
        assert _sanitise_filename("my_contig") == "my_contig"


# ── get_softclip_bases ────────────────────────────────────────────────────────


class TestGetSoftclipBases:
    """Counts only S (soft-clip, op 4) bases; H (hard-clip, op 5) must not be included."""

    def test_none_cigar_returns_zero(self):
        read = MagicMock()
        read.cigartuples = None
        assert get_softclip_bases(read) == 0

    def test_no_softclip_returns_zero(self):
        read = MagicMock()
        read.cigartuples = [(0, 100)]
        assert get_softclip_bases(read) == 0

    def test_left_softclip_only(self):
        read = MagicMock()
        read.cigartuples = [(4, 10), (0, 100)]
        assert get_softclip_bases(read) == 10

    def test_right_softclip_only(self):
        read = MagicMock()
        read.cigartuples = [(0, 100), (4, 15)]
        assert get_softclip_bases(read) == 15

    def test_both_softclips(self):
        read = MagicMock()
        read.cigartuples = [(4, 5), (0, 100), (4, 10)]
        assert get_softclip_bases(read) == 15

    def test_softclip_with_insertion(self):
        read = MagicMock()
        read.cigartuples = [(4, 5), (0, 50), (1, 3), (0, 47), (4, 10)]
        assert get_softclip_bases(read) == 15

    def test_hard_clip_not_counted(self):
        read = MagicMock()
        read.cigartuples = [(5, 10), (0, 100)]
        assert get_softclip_bases(read) == 0

    def test_hard_clip_combined_with_soft_clip(self):
        # 10H 5S 90M → only the 5 soft-clipped bases count; hard clip must not be added
        read = MagicMock()
        read.cigartuples = [(5, 10), (4, 5), (0, 90)]
        assert get_softclip_bases(read) == 5


# ── _make_window ──────────────────────────────────────────────────────────────


class TestMakeWindow:
    """Window dict is initialised with correct start/end and all counters at zero."""

    def test_start_and_end(self):
        w = _make_window(0, 5000)
        assert w["start"] == 0
        assert w["end"] == 5000

    def test_nonzero_start(self):
        w = _make_window(2500, 5000)
        assert w["start"] == 2500
        assert w["end"] == 7500

    def test_all_counters_zero(self):
        w = _make_window(0, 5000)
        assert w["read_count"] == 0
        assert w["total_bases"] == 0
        assert w["softclip_bases"] == 0

    def test_hist_length(self):
        w = _make_window(0, 5000)
        assert len(w["hist"]) == MAX_MAPQ + 1

    def test_hist_all_zeros(self):
        w = _make_window(0, 5000)
        assert all(v == 0 for v in w["hist"])


# ── _create_windows ───────────────────────────────────────────────────────────


class TestCreateWindows:
    """Window creation: stretching, clipping, step advancement, and boundary conditions."""

    def test_no_windows_when_start_exceeds_limit(self):
        windows, next_start, finished = _create_windows(
            next_window_start=6000,
            limit_start=5000,
            chrom_len=10000,
            window_bp=5000,
            step_bp=2500,
        )
        assert windows == []
        assert next_start == 6000
        assert finished is False

    def test_no_windows_when_start_at_chrom_end(self):
        windows, _, _ = _create_windows(
            next_window_start=10000,
            limit_start=10000,
            chrom_len=10000,
            window_bp=5000,
            step_bp=2500,
        )
        assert windows == []

    def test_single_window_basic(self):
        windows, next_start, finished = _create_windows(
            next_window_start=0,
            limit_start=0,
            chrom_len=10000,
            window_bp=5000,
            step_bp=2500,
        )
        assert len(windows) == 1
        assert windows[0]["start"] == 0
        assert windows[0]["end"] == 5000
        assert next_start == 2500
        assert finished is False

    def test_multiple_windows_within_limit(self):
        windows, next_start, finished = _create_windows(
            next_window_start=0,
            limit_start=4999,
            chrom_len=10000,
            window_bp=5000,
            step_bp=2500,
        )
        assert len(windows) == 2
        assert windows[0]["start"] == 0
        assert windows[1]["start"] == 2500
        assert next_start == 5000
        assert finished is False

    def test_stretching_near_chrom_end(self):
        # remainder = 5200-5000 = 200 <= max_stretch(1500) → stretch to chrom_len
        windows, _, finished = _create_windows(
            next_window_start=0,
            limit_start=4999,
            chrom_len=5200,
            window_bp=5000,
            step_bp=2500,
        )
        assert len(windows) == 1
        assert windows[0]["end"] == 5200
        assert finished is True

    def test_no_stretch_when_remainder_too_large(self):
        # remainder = 10000-5000 = 5000 > max_stretch(1500) → no stretch
        windows, _, finished = _create_windows(
            next_window_start=0,
            limit_start=0,
            chrom_len=10000,
            window_bp=5000,
            step_bp=2500,
        )
        assert windows[0]["end"] == 5000
        assert finished is False

    def test_window_clipped_to_chrom_end(self):
        # chrom_len=4800 < window_bp=5000, w_end overshoots → clip
        windows, _, finished = _create_windows(
            next_window_start=0,
            limit_start=4799,
            chrom_len=4800,
            window_bp=5000,
            step_bp=2500,
        )
        assert len(windows) == 1
        assert windows[0]["end"] == 4800
        assert finished is True

    def test_step_advances_next_window_start(self):
        _, next_start, _ = _create_windows(
            next_window_start=0,
            limit_start=2500,
            chrom_len=10000,
            window_bp=5000,
            step_bp=2500,
        )
        assert next_start == 5000


# ── _flush_window_to_csv — flags ──────────────────────────────────────────────


class TestFlushWindowFlags:
    """Flag logic: NO_COVERAGE, LOW_COVERAGE (<5 reads), LOW_DEPTH (<5x), and clean window."""

    def test_no_coverage_when_zero_reads(self):
        row = _flush_and_get_row(_window_with(0, 5000, read_count=0))
        assert row[-1] == "NO_COVERAGE"

    def test_low_coverage_only(self):
        # read_count=3 < 5, depth=30000/5000=6x >= 5
        row = _flush_and_get_row(_window_with(0, 5000, read_count=3, total_bases=30000))
        assert row[-1] == "LOW_COVERAGE"

    def test_low_depth_only(self):
        # read_count=10 >= 5, depth=10000/5000=2x < 5
        row = _flush_and_get_row(
            _window_with(0, 5000, read_count=10, total_bases=10000)
        )
        assert row[-1] == "LOW_DEPTH"

    def test_both_flags(self):
        # read_count=3 < 5, depth=5000/5000=1x < 5
        row = _flush_and_get_row(_window_with(0, 5000, read_count=3, total_bases=5000))
        assert row[-1] == "LOW_COVERAGE|LOW_DEPTH"

    def test_clean_window_empty_flag(self):
        # read_count=10 >= 5, depth=50000/5000=10x >= 5
        row = _flush_and_get_row(
            _window_with(0, 5000, read_count=10, total_bases=50000)
        )
        assert row[-1] == ""

    def test_exactly_at_thresholds_is_clean(self):
        # read_count == threshold (not <), depth == threshold (not <) → clean
        row = _flush_and_get_row(
            _window_with(
                0,
                5000,
                read_count=LOW_COVERAGE_READ_THRESHOLD,
                total_bases=int(LOW_DEPTH_THRESHOLD * 5000),
            )
        )
        assert row[-1] == ""

    def test_one_below_read_threshold(self):
        row = _flush_and_get_row(
            _window_with(
                0, 5000, read_count=LOW_COVERAGE_READ_THRESHOLD - 1, total_bases=50000
            )
        )
        assert "LOW_COVERAGE" in row[-1]

    def test_depth_just_below_threshold(self):
        # depth = 24999/5000 = 4.9998 < 5.0
        row = _flush_and_get_row(
            _window_with(0, 5000, read_count=10, total_bases=24999)
        )
        assert "LOW_DEPTH" in row[-1]


# ── _flush_window_to_csv — output values ──────────────────────────────────────


class TestFlushWindowOutput:
    """CSV row content: MAPQ/softclip fields empty for NO_COVERAGE, correct values otherwise."""

    def test_no_coverage_mapq_empty(self):
        row = _flush_and_get_row(_window_with(0, 5000, read_count=0))
        assert row[3] == ""  # Mean_MAPQ
        assert row[4] == ""  # Median_MAPQ

    def test_no_coverage_softclip_empty(self):
        row = _flush_and_get_row(_window_with(0, 5000, read_count=0))
        assert row[8] == ""  # Softclip_%

    def test_mapq_values_written_correctly(self):
        w = _window_with(
            0, 5000, read_count=10, total_bases=50000, hist_pairs=[(40, 100)]
        )
        row = _flush_and_get_row(w)
        assert row[3] == "40.00"  # Mean_MAPQ
        assert row[4] == "40"  # Median_MAPQ

    def test_softclip_pct_written_correctly(self):
        # 100 softclip / 1000 total = 10.00000%
        w = _window_with(
            0,
            5000,
            read_count=10,
            total_bases=1000,
            softclip_bases=100,
            hist_pairs=[(60, 1000)],
        )
        row = _flush_and_get_row(w)
        assert row[8] == "10.00000"

    def test_chrom_written_correctly(self):
        w = _window_with(0, 5000, read_count=10, total_bases=50000)
        row = _flush_and_get_row(w, chrom="chrX")
        assert row[0] == "chrX"

    def test_start_end_written(self):
        w = _window_with(2500, 7500, read_count=10, total_bases=50000)
        row = _flush_and_get_row(w)
        assert row[1] == 2500
        assert row[2] == 7500

    def test_row_has_ten_columns(self):
        w = _window_with(0, 5000, read_count=10, total_bases=50000)
        assert len(_flush_and_get_row(w)) == 10


# ── validate_inputs ───────────────────────────────────────────────────────────


class TestValidateInputs:
    """CLI argument validation: missing files and out-of-range parameters exit with a message."""

    def test_bam_not_found(self, tmp_path):
        bai = tmp_path / "test.bam.bai"
        bai.touch()
        with pytest.raises(SystemExit, match="BAM file not found"):
            validate_inputs(str(bai), _make_args(str(tmp_path / "missing.bam")))

    def test_bai_not_found(self, tmp_path):
        bam = tmp_path / "test.bam"
        bam.touch()
        with pytest.raises(SystemExit, match="BAM index not found"):
            validate_inputs(str(tmp_path / "missing.bai"), _make_args(str(bam)))

    def test_window_nan(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="must be a finite number"):
            validate_inputs(str(bai), _make_args(str(bam), window=float("nan")))

    def test_window_inf(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="must be a finite number"):
            validate_inputs(str(bai), _make_args(str(bam), window=float("inf")))

    def test_window_zero(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="must be greater than 0"):
            validate_inputs(str(bai), _make_args(str(bam), window=0.0))

    def test_window_negative(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="must be greater than 0"):
            validate_inputs(str(bai), _make_args(str(bam), window=-1.0))

    def test_window_below_min(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="must be at least"):
            validate_inputs(str(bai), _make_args(str(bam), window=0.5))

    def test_window_above_max(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="exceeds maximum"):
            validate_inputs(str(bai), _make_args(str(bam), window=MAX_WINDOW_KB + 1))

    def test_step_nan(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="must be a finite number"):
            validate_inputs(str(bai), _make_args(str(bam), step=float("nan")))

    def test_step_below_min(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="must be at least"):
            validate_inputs(str(bai), _make_args(str(bam), step=0.5))

    def test_step_above_max(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="exceeds maximum"):
            validate_inputs(str(bai), _make_args(str(bam), step=MAX_WINDOW_KB + 1))

    def test_threads_zero(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="must be greater than 0"):
            validate_inputs(str(bai), _make_args(str(bam), threads=0))

    def test_threads_negative(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        with pytest.raises(SystemExit, match="must be greater than 0"):
            validate_inputs(str(bai), _make_args(str(bam), threads=-1))

    def test_step_larger_than_window_warns(self, tmp_path, capsys):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        validate_inputs(str(bai), _make_args(str(bam), window=5.0, step=10.0))
        assert "[WARN]" in capsys.readouterr().out

    def test_valid_inputs_no_error(self, tmp_path):
        bam, bai = tmp_path / "t.bam", tmp_path / "t.bam.bai"
        bam.touch()
        bai.touch()
        validate_inputs(str(bai), _make_args(str(bam), window=5.0, step=2.5, threads=4))


# ── _accumulate_read_into_windows ─────────────────────────────────────────────


class TestAccumulateReadIntoWindows:
    """Base/softclip/MAPQ accumulation per CIGAR op; soft-clips anchored to read ends."""

    def test_empty_cigar_returns_early(self):
        read = _MockRead(cigartuples=[], ref_start=0, ref_end=100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["total_bases"] == 0

    def test_match_op_adds_bases(self):
        read = _MockRead([(0, 100)], ref_start=0, ref_end=100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["total_bases"] == 100
        assert windows[0]["hist"][60] == 100

    def test_match_op_spanning_two_windows(self):
        # 8000M: window [0,5000] gets 5000 bases, window [5000,10000] gets 3000
        read = _MockRead([(0, 8000)], ref_start=0, ref_end=8000)
        w1, w2 = _make_window(0, 5000), _make_window(5000, 5000)
        _accumulate_read_into_windows(read, [w1, w2], 60)
        assert w1["total_bases"] == 5000
        assert w2["total_bases"] == 3000

    def test_match_outside_window_not_counted(self):
        read = _MockRead([(0, 100)], ref_start=6000, ref_end=6100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["total_bases"] == 0

    def test_insertion_adds_bases_at_rpos(self):
        # 50M 10I 50M → 50+50 from M + 10 from I = 110
        read = _MockRead([(0, 50), (1, 10), (0, 50)], ref_start=0, ref_end=100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["total_bases"] == 110

    def test_deletion_advances_position_no_bases(self):
        # 50M 100D 50M → only the 2×50M bases counted
        read = _MockRead([(0, 50), (2, 100), (0, 50)], ref_start=0, ref_end=200)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["total_bases"] == 100

    def test_skip_op_advances_position_no_bases(self):
        # 50M 100N 50M (N = skipped/intron) → only 2×50M
        read = _MockRead([(0, 50), (3, 100), (0, 50)], ref_start=0, ref_end=200)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["total_bases"] == 100

    def test_mapq_none_skips_hist(self):
        read = _MockRead([(0, 100)], ref_start=0, ref_end=100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, None)
        assert windows[0]["total_bases"] == 100
        assert sum(windows[0]["hist"]) == 0

    def test_mapq_255_remapped_to_zero(self):
        # MAPQ 255 = "not available" in SAM spec; remapped to 0 before histogram
        read = _MockRead([(0, 100)], ref_start=0, ref_end=100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 255)
        assert windows[0]["hist"][0] == 100
        assert windows[0]["hist"][255] == 0

    def test_left_softclip_anchored_to_read_start(self):
        # 10S 100M, ref_start=50 → left clip anchored at 50, inside [0, 5000]
        read = _MockRead([(4, 10), (0, 100)], ref_start=50, ref_end=150)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["softclip_bases"] == 10

    def test_right_softclip_anchored_to_read_end(self):
        # 100M 15S, ref_start=0, ref_end=100 → right clip anchored at ref 99
        read = _MockRead([(0, 100), (4, 15)], ref_start=0, ref_end=100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["softclip_bases"] == 15

    def test_both_softclips(self):
        read = _MockRead([(4, 5), (0, 100), (4, 10)], ref_start=50, ref_end=150)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["softclip_bases"] == 15

    def test_hard_clip_before_softclip_left(self):
        # 10H 5S 100M → left_sc = 5
        read = _MockRead([(5, 10), (4, 5), (0, 100)], ref_start=0, ref_end=100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["softclip_bases"] == 5

    def test_sequence_match_op_7(self):
        # op=7 (=, sequence match) treated same as M
        read = _MockRead([(7, 100)], ref_start=0, ref_end=100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["total_bases"] == 100

    def test_sequence_mismatch_op_8(self):
        # op=8 (X, sequence mismatch) treated same as M
        read = _MockRead([(8, 100)], ref_start=0, ref_end=100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 60)
        assert windows[0]["total_bases"] == 100

    def test_mapq_zero_counted_at_hist_zero(self):
        # MAPQ 0 is a valid mapping quality (not a sentinel); must land in hist[0], not be dropped
        read = _MockRead([(0, 100)], ref_start=0, ref_end=100)
        windows = [_make_window(0, 5000)]
        _accumulate_read_into_windows(read, windows, 0)
        assert windows[0]["hist"][0] == 100
        assert sum(windows[0]["hist"]) == 100


# ══════════════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS
# ══════════════════════════════════════════════════════════════════════════════

# ── helpers ───────────────────────────────────────────────────────────────────


def _read_csv(path):
    with open(path, newline="") as f:
        return list(csv.reader(f))


def _make_reads(header, positions, mapq=60, seq_len=100, ref_id=0, cigar=None):
    """Return primary reads at the given reference positions.

    cigar defaults to [(0, seq_len)] (plain M).
    seq_len must equal the sum of all query-consuming CIGAR ops (M, I, S, etc.).
    """
    reads = []
    for i, pos in enumerate(positions):
        r = pysam.AlignedSegment(header)
        r.query_name = f"read_{ref_id}_{i}"
        r.query_sequence = "A" * seq_len
        r.flag = 0
        r.reference_id = ref_id
        r.reference_start = pos
        r.mapping_quality = mapq
        r.cigar = cigar if cigar is not None else [(0, seq_len)]
        r.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
        reads.append(r)
    return reads


def _write_sorted_bam(bam_path, sq_list, reads):
    """Write reads to a sorted, indexed BAM. sq_list: [{"SN": name, "LN": length}, ...]"""
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": sq_list,
        }
    )
    rebuilt = []
    for r in reads:
        seg = pysam.AlignedSegment(header)
        seg.query_name = r.query_name
        seg.query_sequence = r.query_sequence
        seg.flag = r.flag
        seg.reference_id = r.reference_id
        seg.reference_start = r.reference_start
        seg.mapping_quality = r.mapping_quality
        seg.cigar = r.cigar
        seg.query_qualities = r.query_qualities
        rebuilt.append(seg)

    tmp = bam_path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=header) as bam:
        for r in sorted(rebuilt, key=lambda x: (x.reference_id, x.reference_start)):
            bam.write(r)
    pysam.sort("-o", bam_path, tmp)
    os.remove(tmp)
    pysam.index(bam_path)


def _run_in(base, bam_path, window_kb=1.0, step_kb=1.0, threads=1):
    """Create output paths under base, run analysis, return output dict."""
    contigs_folder = str(base / "contigs")
    os.makedirs(contigs_folder, exist_ok=True)
    window_file = str(base / "window_stats.csv")
    summary_file = str(base / "summary_stats.csv")
    run_analysis(
        bam_path,
        bam_path + ".bai",
        window_kb,
        step_kb,
        threads,
        window_file,
        summary_file,
        contigs_folder,
    )
    return {
        "window_file": window_file,
        "summary_file": summary_file,
        "contigs_folder": contigs_folder,
    }


# ── output structure and CSV schema ──────────────────────────────────────────
#
# chr1 (5000 bp), two covered regions with a gap:
#   [0,   1000) — 10 reads × 100M, MAPQ 60
#   [1000,2000) — no reads → NO_COVERAGE
#   [2000,3000) — no reads → NO_COVERAGE
#   [3000,4000) —  5 reads × 100M, MAPQ 60
#   [4000,5000) — no reads reach here → window never opened
#
# Total: 15 reads, 4 windows, all MAPQ 60.


@pytest.fixture(scope="module")
def analysis_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("basic")
    bam_path = str(base / "test.bam")

    sq = [{"SN": "chr1", "LN": 5000}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )

    reads = _make_reads(header, [i * 100 for i in range(10)])  # [0, 1000)
    reads += _make_reads(header, [3000 + i * 100 for i in range(5)])  # [3000, 4000)
    _write_sorted_bam(bam_path, sq, reads)
    return _run_in(base, bam_path)


class TestPrepareOutputDirs:
    """prepare_output_dirs creates the expected folder tree and returns correct paths."""

    def test_creates_output_folder(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        folder, _, _, _ = prepare_output_dirs("sample.bam", 5.0, 2.5)
        assert os.path.isdir(folder)

    def test_creates_contigs_subfolder(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        folder, contigs_folder, _, _ = prepare_output_dirs("sample.bam", 5.0, 2.5)
        assert os.path.isdir(contigs_folder)
        assert contigs_folder == os.path.join(folder, "contigs")

    def test_window_file_path_inside_folder(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        folder, _, window_file, _ = prepare_output_dirs("sample.bam", 5.0, 2.5)
        assert window_file == os.path.join(folder, "window_stats.csv")

    def test_summary_file_path_inside_folder(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        folder, _, _, summary_file = prepare_output_dirs("sample.bam", 5.0, 2.5)
        assert summary_file == os.path.join(folder, "summary_stats.csv")

    def test_folder_name_contains_bam_prefix(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        folder, _, _, _ = prepare_output_dirs("sample.bam", 5.0, 2.5)
        assert "sample" in folder

    def test_folder_name_contains_window_and_step(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        folder, _, _, _ = prepare_output_dirs("sample.bam", 5.0, 2.5)
        assert "5kb" in folder
        assert "2,5kb" in folder


class TestOutputStructure:
    """All expected output files are created at their documented paths."""

    def test_window_stats_csv_exists(self, analysis_output):
        assert os.path.exists(analysis_output["window_file"])

    def test_summary_stats_csv_exists(self, analysis_output):
        assert os.path.exists(analysis_output["summary_file"])

    def test_contig_windows_csv_exists(self, analysis_output):
        assert os.path.exists(
            os.path.join(analysis_output["contigs_folder"], "chr1.windows.csv")
        )

    def test_contig_summary_csv_exists(self, analysis_output):
        assert os.path.exists(
            os.path.join(analysis_output["contigs_folder"], "chr1.summary.csv")
        )


class TestCsvHeaders:
    """CSV column names match the documented schema; contig and global headers are identical."""

    def test_window_stats_header(self, analysis_output):
        rows = _read_csv(analysis_output["window_file"])
        assert rows[0] == [
            "Chromosome",
            "Start",
            "End",
            "Mean_MAPQ",
            "Median_MAPQ",
            "Read_Count",
            "Total_Bases",
            "Softclip_Bases",
            "Softclip_%",
            "Flag",
        ]

    def test_summary_stats_header(self, analysis_output):
        rows = _read_csv(analysis_output["summary_file"])
        assert rows[0] == [
            "Chromosome",
            "Mean_MAPQ",
            "Median_MAPQ",
            "Reads_Seen",
            "Total_Bases",
            "Softclip_Bases",
            "Softclip_%",
            "Windows_Created",
        ]

    def test_contig_windows_header_matches_global(self, analysis_output):
        global_rows = _read_csv(analysis_output["window_file"])
        contig_rows = _read_csv(
            os.path.join(analysis_output["contigs_folder"], "chr1.windows.csv")
        )
        assert contig_rows[0] == global_rows[0]


class TestSummaryStats:
    """summary_stats.csv contains correct per-chrom and GENOME aggregate rows."""

    def _summary(self, analysis_output):
        return _read_csv(analysis_output["summary_file"])

    def test_genome_row_exists(self, analysis_output):
        assert any(r[0] == "GENOME" for r in self._summary(analysis_output)[1:])

    def test_chr1_row_exists(self, analysis_output):
        assert any(r[0] == "chr1" for r in self._summary(analysis_output)[1:])

    def test_genome_total_reads(self, analysis_output):
        genome_row = next(r for r in self._summary(analysis_output) if r[0] == "GENOME")
        assert int(genome_row[3]) == 15  # 10 + 5 reads

    def test_genome_windows_created(self, analysis_output):
        genome_row = next(r for r in self._summary(analysis_output) if r[0] == "GENOME")
        assert int(genome_row[7]) == 4  # [0,1000) [1000,2000) [2000,3000) [3000,4000)

    def test_genome_mapq_is_60(self, analysis_output):
        genome_row = next(r for r in self._summary(analysis_output) if r[0] == "GENOME")
        assert genome_row[1] == "60.00"  # Mean_MAPQ
        assert genome_row[2] == "60"  # Median_MAPQ


class TestWindowFlags:
    """Flag column: gap windows → NO_COVERAGE; covered windows are not."""

    def _rows(self, analysis_output):
        return _read_csv(analysis_output["window_file"])[1:]

    def test_no_coverage_windows_exist(self, analysis_output):
        assert "NO_COVERAGE" in [r[-1] for r in self._rows(analysis_output)]

    def test_gap_windows_are_no_coverage(self, analysis_output):
        gap_rows = [r for r in self._rows(analysis_output) if r[1] in ("1000", "2000")]
        assert len(gap_rows) == 2
        for row in gap_rows:
            assert row[-1] == "NO_COVERAGE"

    def test_covered_windows_not_no_coverage(self, analysis_output):
        covered = [r for r in self._rows(analysis_output) if r[1] in ("0", "3000")]
        assert len(covered) == 2
        for row in covered:
            assert row[-1] != "NO_COVERAGE"

    def test_total_window_count(self, analysis_output):
        assert len(self._rows(analysis_output)) == 4


# ── window stretching ─────────────────────────────────────────────────────────
#
# chrom_len=5200, window=1kb, step=1kb
# remainder after last full window = 5200-5000 = 200 ≤ max_stretch(300) → stretched
# expected windows: [0,1000) [1000,2000) [2000,3000) [3000,4000) [4000,5200)


@pytest.fixture(scope="module")
def stretch_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("stretch")
    bam_path = str(base / "stretch.bam")
    sq = [{"SN": "chr1", "LN": 5200}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )
    reads = _make_reads(header, [0, 1000, 2000, 3000, 4000])
    _write_sorted_bam(bam_path, sq, reads)
    return _run_in(base, bam_path)


class TestWindowStretching:
    """Remainder ≤ 30% of window_bp → last window stretched to chrom_len, not truncated."""

    def _rows(self, stretch_output):
        return _read_csv(stretch_output["window_file"])[1:]

    def test_five_windows_created(self, stretch_output):
        assert len(self._rows(stretch_output)) == 5

    def test_last_window_stretched_to_chrom_end(self, stretch_output):
        last = self._rows(stretch_output)[-1]
        assert last[1] == "4000"  # Start
        assert last[2] == "5200"  # End — stretched, not 5000

    def test_no_window_ends_at_5000(self, stretch_output):
        assert "5000" not in [r[2] for r in self._rows(stretch_output)]


# ── window clipping ───────────────────────────────────────────────────────────
#
# chrom_len=4800, window=1kb, step=1kb
# last window overshoots: w_end=5000 > 4800 → remainder=800 > max_stretch(300) → clipped
# expected windows: [0,1000) [1000,2000) [2000,3000) [3000,4000) [4000,4800)


@pytest.fixture(scope="module")
def clip_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("clip")
    bam_path = str(base / "clip.bam")
    sq = [{"SN": "chr1", "LN": 4800}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )
    reads = _make_reads(header, [0, 1000, 2000, 3000, 4000])
    _write_sorted_bam(bam_path, sq, reads)
    return _run_in(base, bam_path)


class TestWindowClipping:
    """Remainder > 30% of window_bp → overshooting window clipped to chrom_len, not stretched."""

    def _rows(self, clip_output):
        return _read_csv(clip_output["window_file"])[1:]

    def test_five_windows_created(self, clip_output):
        assert len(self._rows(clip_output)) == 5

    def test_last_window_clipped_to_chrom_end(self, clip_output):
        last = self._rows(clip_output)[-1]
        assert last[1] == "4000"  # Start
        assert last[2] == "4800"  # End — clipped, not 5000

    def test_no_window_ends_at_5000(self, clip_output):
        assert "5000" not in [r[2] for r in self._rows(clip_output)]


# ── read spanning window boundary ─────────────────────────────────────────────
#
# One 3000M read at position 0 spans all three 1kb windows.
# Each window receives exactly 1000 bases and 1 read.
# All windows flagged LOW_COVERAGE|LOW_DEPTH (1 read < 5, depth=1x < 5x).


@pytest.fixture(scope="module")
def span_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("span")
    bam_path = str(base / "span.bam")
    sq = [{"SN": "chr1", "LN": 3000}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )

    r = pysam.AlignedSegment(header)
    r.query_name = "long_read"
    r.query_sequence = "A" * 3000
    r.flag = 0
    r.reference_id = 0
    r.reference_start = 0
    r.mapping_quality = 60
    r.cigar = [(0, 3000)]
    r.query_qualities = pysam.qualitystring_to_array("I" * 3000)

    _write_sorted_bam(bam_path, sq, [r])
    return _run_in(base, bam_path)


class TestReadSpanningBoundary:
    """A single read spanning multiple windows contributes bases to each window it overlaps."""

    def _rows(self, span_output):
        return _read_csv(span_output["window_file"])[1:]

    def test_three_windows_created(self, span_output):
        assert len(self._rows(span_output)) == 3

    def test_each_window_has_one_read(self, span_output):
        for row in self._rows(span_output):
            assert row[5] == "1"  # Read_Count

    def test_each_window_gets_1000_bases(self, span_output):
        for row in self._rows(span_output):
            assert row[6] == "1000"  # Total_Bases

    def test_all_windows_flagged_low_coverage_and_low_depth(self, span_output):
        for row in self._rows(span_output):
            assert row[-1] == "LOW_COVERAGE|LOW_DEPTH"


# ── multiple contigs ──────────────────────────────────────────────────────────
#
# chr1 (5000 bp): 10 reads in [0,1000)
# chr2 (3500 bp):  5 reads in [0,1000) + 1 read at 3000
#   → remainder after [2000,3000) = 500 > 300 → no stretch
#   → last window [3000,4000) clipped to 3500
# total reads: 16


@pytest.fixture(scope="module")
def multi_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("multi")
    bam_path = str(base / "multi.bam")
    sq = [{"SN": "chr1", "LN": 5000}, {"SN": "chr2", "LN": 3500}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )

    reads = _make_reads(header, [i * 100 for i in range(10)], ref_id=0)
    reads += _make_reads(header, [i * 100 for i in range(5)], ref_id=1)
    reads += _make_reads(header, [3000], ref_id=1)
    _write_sorted_bam(bam_path, sq, reads)
    return _run_in(base, bam_path)


class TestMultipleContigs:
    """Multiple contigs produce per-contig files and correct GENOME aggregate totals."""

    def test_chr1_contig_file_exists(self, multi_output):
        assert os.path.exists(
            os.path.join(multi_output["contigs_folder"], "chr1.windows.csv")
        )

    def test_chr2_contig_file_exists(self, multi_output):
        assert os.path.exists(
            os.path.join(multi_output["contigs_folder"], "chr2.windows.csv")
        )

    def test_both_chroms_in_global_window_stats(self, multi_output):
        chroms = {r[0] for r in _read_csv(multi_output["window_file"])[1:]}
        assert "chr1" in chroms
        assert "chr2" in chroms

    def test_genome_total_reads(self, multi_output):
        genome_row = next(
            r for r in _read_csv(multi_output["summary_file"]) if r[0] == "GENOME"
        )
        assert int(genome_row[3]) == 16  # 10 chr1 + 5 + 1 chr2

    def test_chr2_last_window_clipped(self, multi_output):
        rows = _read_csv(
            os.path.join(multi_output["contigs_folder"], "chr2.windows.csv")
        )[1:]
        last = rows[-1]
        assert last[1] == "3000"  # Start
        assert last[2] == "3500"  # End — clipped to chr2 length


# ── contig shorter than one window ───────────────────────────────────────────
#
# chr1 (800 bp), window=1kb — the single window overshoots (w_end=1000 > 800)
# and is clipped to [0, 800). 5 reads × 100bp: read_count=5 (no LOW_COVERAGE),
# depth=500/800=0.625x → LOW_DEPTH fires.


@pytest.fixture(scope="module")
def short_contig_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("short_contig")
    bam_path = str(base / "short.bam")
    sq = [{"SN": "chr1", "LN": 800}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )
    reads = _make_reads(header, [0, 100, 200, 300, 400])
    _write_sorted_bam(bam_path, sq, reads)
    return _run_in(base, bam_path)


class TestShortContig:
    """A contig shorter than one window produces exactly one clipped window."""

    def _rows(self, short_contig_output):
        return _read_csv(short_contig_output["window_file"])[1:]

    def test_exactly_one_window(self, short_contig_output):
        assert len(self._rows(short_contig_output)) == 1

    def test_window_starts_at_zero(self, short_contig_output):
        assert self._rows(short_contig_output)[0][1] == "0"

    def test_window_end_clipped_to_contig_length(self, short_contig_output):
        # w_end=0+1000=1000 overshoots 800 → must be clipped, not left as 1000
        assert self._rows(short_contig_output)[0][2] == "800"

    def test_reads_are_counted(self, short_contig_output):
        assert int(self._rows(short_contig_output)[0][5]) == 5  # Read_Count

    def test_flag_is_not_no_coverage(self, short_contig_output):
        assert self._rows(short_contig_output)[0][-1] != "NO_COVERAGE"


# ── LOW_COVERAGE read-count boundary ─────────────────────────────────────────
#
# chr1 (2000 bp), window=1kb:
#   [0,   1000): 4 reads → read_count=4 < 5  → LOW_COVERAGE fires
#   [1000,2000): 5 reads → read_count=5 ≥ 5  → LOW_COVERAGE must NOT fire
# Both windows are below 5x depth → LOW_DEPTH fires in both.


@pytest.fixture(scope="module")
def cov_boundary_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("cov_boundary")
    bam_path = str(base / "boundary.bam")
    sq = [{"SN": "chr1", "LN": 2000}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )
    reads = _make_reads(header, [0, 100, 200, 300])  # 4 reads
    reads += _make_reads(header, [1000, 1100, 1200, 1300, 1400])  # 5 reads
    _write_sorted_bam(bam_path, sq, reads)
    return _run_in(base, bam_path)


class TestLowCoverageBoundary:
    """LOW_COVERAGE fires at read_count < 5 (strict); exactly 5 reads must not trigger it."""

    def _by_start(self, cov_boundary_output):
        return {r[1]: r for r in _read_csv(cov_boundary_output["window_file"])[1:]}

    def test_four_read_window_triggers_low_coverage(self, cov_boundary_output):
        assert "LOW_COVERAGE" in self._by_start(cov_boundary_output)["0"][-1]

    def test_five_read_window_does_not_trigger_low_coverage(self, cov_boundary_output):
        # condition is < 5, so exactly 5 must not fire
        assert "LOW_COVERAGE" not in self._by_start(cov_boundary_output)["1000"][-1]

    def test_four_read_window_read_count(self, cov_boundary_output):
        assert int(self._by_start(cov_boundary_output)["0"][5]) == 4

    def test_five_read_window_read_count(self, cov_boundary_output):
        assert int(self._by_start(cov_boundary_output)["1000"][5]) == 5


# ── softclip base accounting ──────────────────────────────────────────────────
#
# chr1 (1000 bp), 10 reads with CIGAR 10S90M.
# S does not consume the reference; each read spans 90 bp on the reference.
#
# In window [0, 1000):
#   total_bases    = 10 × (90M + 10S) = 1000   (aligned + soft-clip both counted)
#   softclip_bases = 10 × 10          =  100
#   Softclip_%     = 100/1000 × 100   = 10.00000%
#
# total_bases includes soft-clip because the code anchors S bases to the read
# boundary and adds them alongside M bases (see _accumulate_read_into_windows).


@pytest.fixture(scope="module")
def softclip_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("softclip")
    bam_path = str(base / "softclip.bam")
    sq = [{"SN": "chr1", "LN": 1000}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )
    # 10S90M — query length=100, reference length=90
    reads = _make_reads(
        header, [i * 90 for i in range(10)], seq_len=100, cigar=[(4, 10), (0, 90)]
    )
    _write_sorted_bam(bam_path, sq, reads)
    return _run_in(base, bam_path)


class TestSoftclipContent:
    """Soft-clip bases are counted in both Softclip_Bases and Total_Bases; percentage is correct."""

    def _rows(self, softclip_output):
        return _read_csv(softclip_output["window_file"])[1:]

    def test_softclip_bases_exact(self, softclip_output):
        # 10 reads × 10 soft-clipped bases each
        assert int(self._rows(softclip_output)[0][7]) == 100

    def test_total_bases_includes_softclip(self, softclip_output):
        # aligned (90M) + soft-clip (10S) both counted → 10 × 100 = 1000
        assert int(self._rows(softclip_output)[0][6]) == 1000

    def test_softclip_pct_exact(self, softclip_output):
        # 100 / 1000 × 100 = 10.00000% (5 decimal places)
        assert self._rows(softclip_output)[0][8] == "10.00000"

    def test_softclip_pct_nonzero(self, softclip_output):
        assert self._rows(softclip_output)[0][8] != "0.00000"


# ── MAPQ 0 reads ──────────────────────────────────────────────────────────────
#
# MAPQ 0 is a valid mapping quality (not the unmapped sentinel 255).
# Reads must be counted and produce numeric MAPQ fields — not empty strings,
# which are reserved for NO_COVERAGE windows.


@pytest.fixture(scope="module")
def mapq_zero_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("mapq_zero")
    bam_path = str(base / "mapq0.bam")
    sq = [{"SN": "chr1", "LN": 1000}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )
    reads = _make_reads(header, [i * 90 for i in range(10)], mapq=0)
    _write_sorted_bam(bam_path, sq, reads)
    return _run_in(base, bam_path)


class TestMapqZeroReads:
    """MAPQ 0 is a valid mapping quality; reads must be counted and produce numeric output fields."""

    def _row(self, mapq_zero_output):
        return next(
            r for r in _read_csv(mapq_zero_output["window_file"])[1:] if r[5] != "0"
        )

    def test_reads_are_counted(self, mapq_zero_output):
        assert int(self._row(mapq_zero_output)[5]) > 0

    def test_mean_mapq_is_zero_not_empty(self, mapq_zero_output):
        # empty string is reserved for NO_COVERAGE; a covered window must have a value
        assert self._row(mapq_zero_output)[3] == "0.00"

    def test_median_mapq_is_zero_not_empty(self, mapq_zero_output):
        assert self._row(mapq_zero_output)[4] == "0"

    def test_flag_is_not_no_coverage(self, mapq_zero_output):
        assert self._row(mapq_zero_output)[-1] != "NO_COVERAGE"


# ── MAPQ 255 remapping ────────────────────────────────────────────────────────
#
# MAPQ 255 = "not available" in the SAM spec. The tool remaps it to 0 so it is
# not treated as a perfect mapping. The unit-level remap is tested above in
# TestAccumulateReadIntoWindows. This test verifies the remap survives the full
# pipeline and appears correctly in the CSV output.


@pytest.fixture(scope="module")
def mapq255_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("mapq255")
    bam_path = str(base / "mapq255.bam")
    sq = [{"SN": "chr1", "LN": 1000}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )
    reads = _make_reads(header, [i * 90 for i in range(10)], mapq=255)
    _write_sorted_bam(bam_path, sq, reads)
    return _run_in(base, bam_path)


class TestMapq255Remapping:
    """MAPQ 255 reads are counted but their quality is remapped to 0 in CSV output."""

    def _row(self, mapq255_output):
        return next(
            r for r in _read_csv(mapq255_output["window_file"])[1:] if r[5] != "0"
        )

    def test_reads_counted(self, mapq255_output):
        assert int(self._row(mapq255_output)[5]) == 10

    def test_mean_mapq_remapped_to_zero(self, mapq255_output):
        # MAPQ 255 → remapped to 0; mean must be 0.00, not 255.00
        assert self._row(mapq255_output)[3] == "0.00"

    def test_median_mapq_remapped_to_zero(self, mapq255_output):
        assert self._row(mapq255_output)[4] == "0"


# ── secondary and supplementary read filtering ────────────────────────────────
#
# Reads flagged secondary (0x100) or supplementary (0x800) must be skipped.
# They must not contribute to read_count, total_bases, or MAPQ.
#
# chr1 (1000 bp):
#   5 primary reads × 100M, MAPQ 60   → counted
#   5 secondary reads at same pos     → skipped
#   5 supplementary reads at same pos → skipped
# Expected: Read_Count=5, Total_Bases=500


@pytest.fixture(scope="module")
def filtered_reads_output(tmp_path_factory):
    base = tmp_path_factory.mktemp("filtered_reads")
    bam_path = str(base / "filtered.bam")
    sq = [{"SN": "chr1", "LN": 1000}]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    )

    positions = [i * 100 for i in range(5)]
    primary_reads = _make_reads(header, positions, mapq=60)

    extra = []
    for flag, label in [(256, "sec"), (2048, "sup")]:
        for i, pos in enumerate(positions):
            r = pysam.AlignedSegment(header)
            r.query_name = f"{label}_{i}"
            r.query_sequence = "A" * 100
            r.flag = flag
            r.reference_id = 0
            r.reference_start = pos
            r.mapping_quality = 60
            r.cigar = [(0, 100)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            extra.append(r)

    _write_sorted_bam(bam_path, sq, primary_reads + extra)
    return _run_in(base, bam_path)


class TestFilteredReads:
    """Secondary and supplementary reads are silently skipped; only primary reads are counted."""

    def _row(self, filtered_reads_output):
        return _read_csv(filtered_reads_output["window_file"])[1:][0]

    def test_only_primary_reads_counted(self, filtered_reads_output):
        # 5 primary + 5 secondary + 5 supplementary → only 5 must appear
        assert int(self._row(filtered_reads_output)[5]) == 5  # Read_Count

    def test_bases_from_primary_reads_only(self, filtered_reads_output):
        # 5 × 100M = 500; secondary/supplementary bases must not be included
        assert int(self._row(filtered_reads_output)[6]) == 500  # Total_Bases

    def test_flag_is_not_no_coverage(self, filtered_reads_output):
        assert self._row(filtered_reads_output)[-1] != "NO_COVERAGE"
