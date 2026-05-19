"""
Tests for the hese module.

Organised as:
  UNIT TESTS        - isolated function tests using synthetic in-memory inputs
  INTEGRATION TESTS - full pipeline calls against synthetic text files written to tmp_path
"""

import argparse
import os

import pytest

from raaqa.hese import (
    _TRUTH_SUFFIX_PAIRS,
    _scan_paf_tnames,
    build_chrom_groups,
    build_unitig_table,
    carry_forward_amb,
    choose_best_orientation,
    chrom_metrics,
    compute_unitig_label,
    detect_truth_suffixes,
    haplotype_struct_metrics,
    parse_path_paf,
    parse_truth_assignments,
    read_idxstats,
    resolve_truth_suffixes,
    score_orientation,
    summarize_haplotigs,
    validate_inputs,
    write_haplotig_table,
    write_haplotype_summary,
    write_truth_assignments,
    write_truth_chrom_metrics,
    write_truth_eval_summary,
    write_unitig_table,
)


# ══════════════════════════════════════════════════════════════════════════════
# UNIT TESTS
# ══════════════════════════════════════════════════════════════════════════════

# ── helpers ───────────────────────────────────────────────────────────────────


def _write_file(path, content):
    with open(path, "w") as f:
        f.write(content)


def _make_args(**kwargs):
    defaults = dict(
        idx_p1="p1.idxstats",
        idx_p2="p2.idxstats",
        hap1_paf="hap1.paf",
        hap2_paf="hap2.paf",
        truth_hap1_paf=None,
        truth_hap2_paf=None,
        hap1_suffix=None,
        hap2_suffix=None,
        truth_min_span=50000,
        truth_min_best_frac=0.60,
        min_unitig_length=2000,
        min_total_reads=50,
        unitig_frac_threshold=0.60,
        min_unitigs_per_haplotig=3,
        hap_frac_threshold=0.60,
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


# ── read_idxstats ─────────────────────────────────────────────────────────────


class TestReadIdxstats:
    """Parses samtools idxstats: skips * row, accumulates total_mapped, handles edge cases."""

    def test_normal_lines_parsed(self, tmp_path):
        p = tmp_path / "p1.idxstats"
        _write_file(p, "utg1\t10000\t120\t3\nutg2\t5000\t80\t1\n")
        stats, total = read_idxstats(str(p))
        assert stats["utg1"] == {"length": 10000, "mapped": 120, "unmapped": 3}
        assert stats["utg2"] == {"length": 5000, "mapped": 80, "unmapped": 1}

    def test_total_mapped_accumulated(self, tmp_path):
        p = tmp_path / "p1.idxstats"
        _write_file(p, "utg1\t10000\t120\t3\nutg2\t5000\t80\t1\n")
        _, total = read_idxstats(str(p))
        assert total == 200

    def test_star_row_skipped(self, tmp_path):
        p = tmp_path / "p1.idxstats"
        _write_file(p, "utg1\t10000\t50\t2\n*\t0\t0\t999\n")
        stats, total = read_idxstats(str(p))
        assert "*" not in stats
        assert total == 50

    def test_empty_file_returns_empty(self, tmp_path):
        p = tmp_path / "empty.idxstats"
        _write_file(p, "")
        stats, total = read_idxstats(str(p))
        assert stats == {}
        assert total == 0

    def test_blank_lines_skipped(self, tmp_path):
        p = tmp_path / "p1.idxstats"
        _write_file(p, "\nutg1\t10000\t10\t0\n\n")
        stats, _ = read_idxstats(str(p))
        assert list(stats.keys()) == ["utg1"]

    def test_short_lines_skipped(self, tmp_path):
        p = tmp_path / "p1.idxstats"
        _write_file(p, "utg1\t10000\t10\nutg2\t5000\t20\t1\n")
        stats, _ = read_idxstats(str(p))
        assert "utg1" not in stats
        assert "utg2" in stats

    def test_malformed_int_exits(self, tmp_path):
        p = tmp_path / "p1.idxstats"
        _write_file(p, "utg1\t10000\tXXX\t0\n")
        with pytest.raises(SystemExit, match="malformed idxstats"):
            read_idxstats(str(p))

    def test_malformed_error_includes_line_number(self, tmp_path):
        p = tmp_path / "p1.idxstats"
        _write_file(p, "utg1\t10000\t50\t0\nutg2\t5000\tBAD\t0\n")
        with pytest.raises(SystemExit, match="2"):
            read_idxstats(str(p))


# ── parse_path_paf ────────────────────────────────────────────────────────────


class TestParsePathPaf:
    """Parses hifiasm path PAF: groups by hap_id, sorts by qstart, skips short lines."""

    def test_normal_entries_parsed(self, tmp_path):
        p = tmp_path / "hap1.paf"
        _write_file(p, "h1\t5000\t0\t5000\t+\tutg1\t10000\t0\t5000\t5000\t5000\t60\n")
        paths = parse_path_paf(str(p))
        assert paths == {"h1": ["utg1"]}

    def test_multiple_entries_sorted_by_qstart(self, tmp_path):
        p = tmp_path / "hap1.paf"
        _write_file(p,
            "h1\t9000\t5000\t9000\t+\tutg2\t5000\t0\t4000\t4000\t4000\t60\n"
            "h1\t9000\t0\t5000\t+\tutg1\t5000\t0\t5000\t5000\t5000\t60\n"
        )
        paths = parse_path_paf(str(p))
        assert paths["h1"] == ["utg1", "utg2"]

    def test_multiple_hap_ids_grouped(self, tmp_path):
        p = tmp_path / "hap1.paf"
        _write_file(p,
            "h1\t5000\t0\t5000\t+\tutg1\t5000\t0\t5000\t5000\t5000\t60\n"
            "h2\t4000\t0\t4000\t+\tutg2\t4000\t0\t4000\t4000\t4000\t60\n"
        )
        paths = parse_path_paf(str(p))
        assert set(paths.keys()) == {"h1", "h2"}

    def test_short_lines_skipped(self, tmp_path):
        p = tmp_path / "hap1.paf"
        _write_file(p, "h1\t5000\t0\t5000\t+\n")
        paths = parse_path_paf(str(p))
        assert paths == {}

    def test_blank_lines_skipped(self, tmp_path):
        p = tmp_path / "hap1.paf"
        _write_file(p, "\nh1\t5000\t0\t5000\t+\tutg1\t5000\t0\t5000\t5000\t5000\t60\n\n")
        paths = parse_path_paf(str(p))
        assert "h1" in paths

    def test_empty_file_returns_empty(self, tmp_path):
        p = tmp_path / "hap1.paf"
        _write_file(p, "")
        assert parse_path_paf(str(p)) == {}

    def test_malformed_qstart_exits(self, tmp_path):
        p = tmp_path / "hap1.paf"
        _write_file(p, "h1\t5000\tBAD\t5000\t+\tutg1\t5000\t0\t5000\t5000\t5000\t60\n")
        with pytest.raises(SystemExit, match="malformed PAF"):
            parse_path_paf(str(p))


# ── compute_unitig_label ──────────────────────────────────────────────────────


class TestComputeUnitigLabel:
    """P1/P2/amb labelling: length/read thresholds, library-size normalisation, boundary values."""

    _DEFAULTS = dict(
        total_p1_reads=1000,
        total_p2_reads=1000,
        min_unitig_length=2000,
        min_total_reads=50,
        frac_threshold=0.60,
    )

    def _call(self, p1, p2, length, **overrides):
        kw = {**self._DEFAULTS, **overrides}
        return compute_unitig_label(p1, p2, length, **kw)

    def test_below_length_threshold_is_amb(self):
        label, _, _ = self._call(500, 10, length=1999)
        assert label == "amb"

    def test_below_read_threshold_is_amb(self):
        label, _, _ = self._call(30, 10, length=5000)
        assert label == "amb"

    def test_strong_p1_fraction_is_p1(self):
        label, _, _ = self._call(800, 100, length=5000)
        assert label == "P1"

    def test_strong_p2_fraction_is_p2(self):
        label, _, _ = self._call(100, 800, length=5000)
        assert label == "P2"

    def test_balanced_fracs_below_threshold_is_amb(self):
        # equal reads: frac_p1 = frac_p2 = 0.5, both below 0.60
        label, _, _ = self._call(500, 500, length=5000)
        assert label == "amb"

    def test_zero_total_p1_library_p2_wins(self):
        # p1_norm forced to 0.0: frac_p2 = 1.0, P2 labelled
        label, _, _ = self._call(100, 500, length=5000, total_p1_reads=0)
        assert label == "P2"

    def test_zero_total_p2_library_p1_wins(self):
        # p2_norm forced to 0.0: frac_p1 = 1.0, P1 labelled
        label, _, _ = self._call(500, 100, length=5000, total_p2_reads=0)
        assert label == "P1"

    def test_zero_norm_denom_is_amb(self):
        # both libraries zero: p1_norm + p2_norm == 0
        label, _, _ = self._call(100, 100, length=5000, total_p1_reads=0, total_p2_reads=0)
        assert label == "amb"

    def test_frac_exactly_at_threshold_is_labelled(self):
        # frac_p1 must be exactly 0.60 and > frac_p2, labels P1
        # p1_norm/p2_norm = 0.60/0.40: frac_p1 = 0.6/1.0 = 0.60
        label, _, _ = self._call(600, 400, length=5000, total_p1_reads=1000, total_p2_reads=1000)
        assert label == "P1"

    def test_norm_values_returned(self):
        # p1_norm = 600/1000 = 0.6, p2_norm = 100/1000 = 0.1
        _, p1n, p2n = self._call(600, 100, length=5000)
        assert p1n == pytest.approx(0.6)
        assert p2n == pytest.approx(0.1)


# ── build_unitig_table ────────────────────────────────────────────────────────


class TestBuildUnitigTable:
    """Merges p1/p2 index dicts: handles unitigs present in one or both, delegates labelling."""

    _SHARED = dict(
        total_p1_reads=1000,
        total_p2_reads=1000,
        min_unitig_length=2000,
        min_total_reads=10,
        frac_threshold=0.60,
    )

    def _build(self, p1_idx, p2_idx, **overrides):
        kw = {**self._SHARED, **overrides}
        return build_unitig_table(p1_idx, kw.pop("total_p1_reads"), p2_idx, kw.pop("total_p2_reads"), **kw)

    def test_unitig_only_in_p1(self):
        p1 = {"utg1": {"length": 5000, "mapped": 500, "unmapped": 0}}
        result = self._build(p1, {})
        assert result["utg1"]["p2_reads"] == 0
        assert result["utg1"]["length"] == 5000

    def test_unitig_only_in_p2(self):
        p2 = {"utg1": {"length": 4000, "mapped": 400, "unmapped": 0}}
        result = self._build({}, p2)
        assert result["utg1"]["p1_reads"] == 0
        assert result["utg1"]["length"] == 4000

    def test_unitig_in_both(self):
        p1 = {"utg1": {"length": 5000, "mapped": 600, "unmapped": 0}}
        p2 = {"utg1": {"length": 5000, "mapped": 100, "unmapped": 0}}
        result = self._build(p1, p2)
        assert result["utg1"]["p1_reads"] == 600
        assert result["utg1"]["p2_reads"] == 100

    def test_all_unitigs_included(self):
        p1 = {"utg1": {"length": 5000, "mapped": 600, "unmapped": 0}}
        p2 = {"utg2": {"length": 3000, "mapped": 400, "unmapped": 0}}
        result = self._build(p1, p2)
        assert set(result.keys()) == {"utg1", "utg2"}

    def test_total_reads_field(self):
        p1 = {"utg1": {"length": 5000, "mapped": 300, "unmapped": 0}}
        p2 = {"utg1": {"length": 5000, "mapped": 200, "unmapped": 0}}
        result = self._build(p1, p2)
        assert result["utg1"]["total_reads"] == 500

    def test_label_assigned(self):
        p1 = {"utg1": {"length": 5000, "mapped": 800, "unmapped": 0}}
        p2 = {"utg1": {"length": 5000, "mapped": 100, "unmapped": 0}}
        result = self._build(p1, p2)
        assert result["utg1"]["label"] in ("P1", "P2", "amb")


# ── summarize_haplotigs ───────────────────────────────────────────────────────


class TestSummarizeHaplotigs:
    """Labels haplotigs P1/P2/amb based on their unitig composition."""

    def _unitigs(self, labels):
        return {f"utg{i}": {"label": lbl} for i, lbl in enumerate(labels)}

    def test_majority_p1_labelled_p1(self):
        unitigs = self._unitigs(["P1", "P1", "P1", "P2"])
        paths = {"h1": [f"utg{i}" for i in range(4)]}
        result = summarize_haplotigs(paths, unitigs, min_unitigs_per_haplotig=3, hap_frac_threshold=0.60)
        assert result["h1"]["label"] == "P1"

    def test_majority_p2_labelled_p2(self):
        unitigs = self._unitigs(["P2", "P2", "P2", "P1"])
        paths = {"h1": [f"utg{i}" for i in range(4)]}
        result = summarize_haplotigs(paths, unitigs, min_unitigs_per_haplotig=3, hap_frac_threshold=0.60)
        assert result["h1"]["label"] == "P2"

    def test_too_few_informative_is_amb(self):
        # 2 informative (P1+P2) < min_unitigs_per_haplotig=3, label is amb
        unitigs = self._unitigs(["P1", "P2", "amb", "amb"])
        paths = {"h1": [f"utg{i}" for i in range(4)]}
        result = summarize_haplotigs(paths, unitigs, min_unitigs_per_haplotig=3, hap_frac_threshold=0.60)
        assert result["h1"]["label"] == "amb"

    def test_balanced_p1_p2_is_amb(self):
        # frac_p1 = frac_p2 = 0.5 < 0.60, label is amb
        unitigs = self._unitigs(["P1", "P1", "P1", "P2", "P2", "P2"])
        paths = {"h1": [f"utg{i}" for i in range(6)]}
        result = summarize_haplotigs(paths, unitigs, min_unitigs_per_haplotig=3, hap_frac_threshold=0.60)
        assert result["h1"]["label"] == "amb"

    def test_unknown_unitig_treated_as_amb(self):
        # utg_missing not in unitigs dict: .get returns None, counted as amb
        unitigs = {"utg0": {"label": "P1"}, "utg1": {"label": "P1"}, "utg2": {"label": "P1"}}
        paths = {"h1": ["utg0", "utg1", "utg2", "utg_missing"]}
        result = summarize_haplotigs(paths, unitigs, min_unitigs_per_haplotig=3, hap_frac_threshold=0.60)
        assert result["h1"]["n_amb"] == 1
        assert result["h1"]["label"] == "P1"

    def test_all_amb_unitigs_is_amb(self):
        # n_inf = 0, below any min_unitigs threshold, label is amb
        unitigs = self._unitigs(["amb", "amb", "amb", "amb"])
        paths = {"h1": [f"utg{i}" for i in range(4)]}
        result = summarize_haplotigs(paths, unitigs, min_unitigs_per_haplotig=3, hap_frac_threshold=0.60)
        assert result["h1"]["label"] == "amb"

    def test_counts_recorded(self):
        unitigs = self._unitigs(["P1", "P1", "P2", "amb"])
        paths = {"h1": [f"utg{i}" for i in range(4)]}
        result = summarize_haplotigs(paths, unitigs, min_unitigs_per_haplotig=3, hap_frac_threshold=0.60)
        assert result["h1"]["n_p1"] == 2
        assert result["h1"]["n_p2"] == 1
        assert result["h1"]["n_amb"] == 1
        assert result["h1"]["n_unitigs"] == 4


# ── carry_forward_amb ─────────────────────────────────────────────────────────


class TestCarryForwardAmb:
    """Replaces amb with the preceding P1/P2; falls back to default_label at the start."""

    def test_no_amb_unchanged(self):
        assert carry_forward_amb(["P1", "P2", "P1"], "P1") == ["P1", "P2", "P1"]

    def test_amb_at_start_uses_default(self):
        assert carry_forward_amb(["amb", "P1"], "P2") == ["P2", "P1"]

    def test_amb_in_middle_uses_preceding(self):
        assert carry_forward_amb(["P1", "amb", "P2"], "P1") == ["P1", "P1", "P2"]

    def test_consecutive_ambs_all_filled(self):
        assert carry_forward_amb(["P2", "amb", "amb", "P1"], "P1") == ["P2", "P2", "P2", "P1"]

    def test_all_amb_uses_default_throughout(self):
        assert carry_forward_amb(["amb", "amb", "amb"], "P1") == ["P1", "P1", "P1"]

    def test_empty_list_returns_empty(self):
        assert carry_forward_amb([], "P1") == []

    def test_trailing_amb_uses_last_known(self):
        assert carry_forward_amb(["P1", "amb"], "P2") == ["P1", "P1"]


# ── haplotype_struct_metrics ──────────────────────────────────────────────────


def _hap_entry(label):
    return {"label": label}


class TestHaplotypeStructMetrics:
    """Computes global parent, Hamming %, and switch % from a haplotig label dict."""

    def test_empty_dict_all_na(self):
        m = haplotype_struct_metrics({})
        assert m["hap_global_parent"] == "amb"
        assert m["hamming_struct_%"] == "NA"
        assert m["switch_struct_%"] == "NA"

    def test_single_haplotig_switch_is_na(self):
        m = haplotype_struct_metrics({"h1": _hap_entry("P1")})
        assert m["switch_struct_%"] == "NA"

    def test_all_p1_hamming_zero_switch_zero(self):
        haps = {f"h{i}": _hap_entry("P1") for i in range(4)}
        m = haplotype_struct_metrics(haps)
        assert m["hap_global_parent"] == "P1"
        assert m["hamming_struct_%"] == pytest.approx(0.0)
        assert m["switch_struct_%"] == pytest.approx(0.0)

    def test_all_p2_hamming_zero_switch_zero(self):
        haps = {f"h{i}": _hap_entry("P2") for i in range(4)}
        m = haplotype_struct_metrics(haps)
        assert m["hap_global_parent"] == "P2"
        assert m["hamming_struct_%"] == pytest.approx(0.0)

    def test_p1_majority_hamming_counts_p2_wrong(self):
        # 3 P1, 1 P2: global parent P1, hamming = 1/4 * 100 = 25%
        haps = {"h0": _hap_entry("P1"), "h1": _hap_entry("P1"),
                "h2": _hap_entry("P1"), "h3": _hap_entry("P2")}
        m = haplotype_struct_metrics(haps)
        assert m["hap_global_parent"] == "P1"
        assert m["hamming_struct_%"] == pytest.approx(25.0)

    def test_tied_p1_p2_is_amb_and_na(self):
        haps = {"h0": _hap_entry("P1"), "h1": _hap_entry("P2")}
        m = haplotype_struct_metrics(haps)
        assert m["hap_global_parent"] == "amb"
        assert m["hamming_struct_%"] == "NA"
        assert m["switch_struct_%"] == "NA"

    def test_switch_computed_correctly(self):
        # 1 P1, 3 P2: P2 majority, sorted: h0=P1 h1=P2 h2=P2 h3=P2, 1 switch / 3 boundaries
        haps = {"h0": _hap_entry("P1"), "h1": _hap_entry("P2"),
                "h2": _hap_entry("P2"), "h3": _hap_entry("P2")}
        m = haplotype_struct_metrics(haps)
        assert m["hap_global_parent"] == "P2"
        assert m["switch_struct_%"] == pytest.approx(100.0 * 1 / 3)

    def test_counts_correct(self):
        haps = {"h0": _hap_entry("P1"), "h1": _hap_entry("P2"),
                "h2": _hap_entry("amb"), "h3": _hap_entry("P1")}
        m = haplotype_struct_metrics(haps)
        assert m["n_haplotigs_total"] == 4
        assert m["n_haplotigs_p1"] == 2
        assert m["n_haplotigs_p2"] == 1
        assert m["n_haplotigs_amb"] == 1


# ── _scan_paf_tnames ──────────────────────────────────────────────────────────


class TestScanPafTnames:
    """Collects uppercased target names from column 6 of a PAF; skips comments, blanks, short lines."""

    def _paf_line(self, tname):
        return f"qname\t5000\t0\t5000\t+\t{tname}\t10000\t0\t5000\t5000\t5000\t60\n"

    def test_normal_line_tname_collected(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, self._paf_line("chr1_HAP1"))
        assert "CHR1_HAP1" in _scan_paf_tnames(str(p))

    def test_tname_uppercased(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, self._paf_line("chr1_hap1"))
        assert "CHR1_HAP1" in _scan_paf_tnames(str(p))

    def test_comment_lines_skipped(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, "# this is a comment\n" + self._paf_line("chr1_HAP1"))
        result = _scan_paf_tnames(str(p))
        assert len(result) == 1

    def test_blank_lines_skipped(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, "\n" + self._paf_line("chr1_HAP1") + "\n")
        result = _scan_paf_tnames(str(p))
        assert len(result) == 1

    def test_short_lines_skipped(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, "qname\t5000\t0\t5000\t+\n")
        assert _scan_paf_tnames(str(p)) == set()

    def test_empty_file_returns_empty_set(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, "")
        assert _scan_paf_tnames(str(p)) == set()

    def test_multiple_tnames_all_collected(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, self._paf_line("chr1_HAP1") + self._paf_line("chr2_HAP2"))
        result = _scan_paf_tnames(str(p))
        assert result == {"CHR1_HAP1", "CHR2_HAP2"}


# ── detect_truth_suffixes ─────────────────────────────────────────────────────


class TestDetectTruthSuffixes:
    """Returns the first known suffix pair where both suffixes appear in the PAF tnames."""

    def _paf_with_tnames(self, tmp_path, name, *tnames):
        p = tmp_path / name
        lines = ""
        for tname in tnames:
            lines += f"q\t5000\t0\t5000\t+\t{tname}\t10000\t0\t5000\t5000\t5000\t60\n"
        _write_file(p, lines)
        return str(p)

    def test_each_known_pair_detected(self, tmp_path):
        for s1, s2 in _TRUTH_SUFFIX_PAIRS:
            paf = self._paf_with_tnames(tmp_path, f"t{s1}.paf", f"chr1{s1}", f"chr1{s2}")
            r1, r2 = detect_truth_suffixes(paf)
            assert r1 == s1
            assert r2 == s2

    def test_only_one_suffix_returns_none(self, tmp_path):
        paf = self._paf_with_tnames(tmp_path, "one.paf", "chr1_HAP1")
        assert detect_truth_suffixes(paf) == (None, None)

    def test_unknown_naming_returns_none(self, tmp_path):
        paf = self._paf_with_tnames(tmp_path, "unk.paf", "chr1_CUSTOM1", "chr1_CUSTOM2")
        assert detect_truth_suffixes(paf) == (None, None)

    def test_empty_file_returns_none(self, tmp_path):
        p = tmp_path / "empty.paf"
        _write_file(p, "")
        assert detect_truth_suffixes(str(p)) == (None, None)

    def test_case_insensitive_detection(self, tmp_path):
        # _scan_paf_tnames uppercases, so lowercase suffixes in file still match
        paf = self._paf_with_tnames(tmp_path, "lower.paf", "chr1_hap1", "chr1_hap2")
        r1, r2 = detect_truth_suffixes(paf)
        assert r1 == "_HAP1"
        assert r2 == "_HAP2"


# ── parse_truth_assignments ───────────────────────────────────────────────────


def _truth_paf_line(qname, qstart, qend, tname, tstart, tend):
    qlen = qend
    tlen = 200000
    # 12-column MashMap PAF line
    return f"{qname}\t{qlen}\t{qstart}\t{qend}\t+\t{tname}\t{tlen}\t{tstart}\t{tend}\t{qend - qstart}\t{qend - qstart}\t60\n"


class TestParseTruthAssignments:
    """Filters by span and best-frac, assigns truth_bin from suffix, aggregates multi-segment."""

    def test_segment_below_min_span_filtered(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, _truth_paf_line("h1", 0, 10000, "chr1_HAP1", 0, 10000))
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert "h1" not in result

    def test_segment_above_min_span_included(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, _truth_paf_line("h1", 0, 60000, "chr1_HAP1", 0, 60000))
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert "h1" in result

    def test_hap1_suffix_assigned_hap1_bin(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, _truth_paf_line("h1", 0, 60000, "chr1_HAP1", 0, 60000))
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert result["h1"]["truth_bin"] == "HAP1"

    def test_hap2_suffix_assigned_hap2_bin(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, _truth_paf_line("h1", 0, 60000, "chr1_HAP2", 0, 60000))
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert result["h1"]["truth_bin"] == "HAP2"

    def test_chrom_stripped_of_suffix(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, _truth_paf_line("h1", 0, 60000, "chr5_HAP1", 0, 60000))
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert result["h1"]["chrom"] == "chr5"

    def test_neither_suffix_excluded(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, _truth_paf_line("h1", 0, 60000, "chr1_OTHER", 0, 60000))
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert "h1" not in result

    def test_suffix_match_is_case_insensitive(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, _truth_paf_line("h1", 0, 60000, "chr1_hap1", 0, 60000))
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert result["h1"]["truth_bin"] == "HAP1"

    def test_best_frac_below_threshold_excluded(self, tmp_path):
        p = tmp_path / "t.paf"
        # h1 splits equally across HAP1 and HAP2: best_frac = 0.5 < 0.60
        _write_file(p,
            _truth_paf_line("h1", 0, 60000, "chr1_HAP1", 0, 60000) +
            _truth_paf_line("h1", 60000, 120000, "chr1_HAP2", 0, 60000)
        )
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert "h1" not in result

    def test_multiple_segments_aggregated(self, tmp_path):
        p = tmp_path / "t.paf"
        # two segments both on HAP1: total bp = 120000, best_frac = 1.0
        _write_file(p,
            _truth_paf_line("h1", 0, 60000, "chr1_HAP1", 0, 60000) +
            _truth_paf_line("h1", 60000, 120000, "chr1_HAP1", 60000, 120000)
        )
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert result["h1"]["best_bp"] == 120000

    def test_tstart_tend_track_min_max(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p,
            _truth_paf_line("h1", 0, 60000, "chr1_HAP1", 100, 60100) +
            _truth_paf_line("h1", 60000, 120000, "chr1_HAP1", 55000, 115000)
        )
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert result["h1"]["tstart"] == 100
        assert result["h1"]["tend"] == 115000

    def test_short_lines_skipped(self, tmp_path):
        p = tmp_path / "t.paf"
        _write_file(p, "h1\t60000\t0\t60000\t+\tchr1_HAP1\t200000\t0\t60000\t60000\t60000\n")
        result = parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                         hap1_suffix="_HAP1", hap2_suffix="_HAP2")
        assert "h1" not in result

    def test_malformed_coord_exits(self, tmp_path):
        p = tmp_path / "t.paf"
        line = "h1\t60000\tBAD\t60000\t+\tchr1_HAP1\t200000\t0\t60000\t60000\t60000\t60\n"
        _write_file(p, line)
        with pytest.raises(SystemExit, match="malformed PAF"):
            parse_truth_assignments(str(p), min_span=50000, min_best_frac=0.60,
                                    hap1_suffix="_HAP1", hap2_suffix="_HAP2")


# ── build_chrom_groups ────────────────────────────────────────────────────────


class TestBuildChromGroups:
    """Groups truth haplotigs by (chrom, truth_bin), sorted by tstart within each group."""

    def _entry(self, chrom, truth_bin, tstart, tend=None):
        return {"chrom": chrom, "truth_bin": truth_bin, "tstart": tstart, "tend": tend or tstart + 1000,
                "track": f"{chrom}_{truth_bin}", "best_bp": 60000, "all_bp": 60000, "best_frac": 1.0}

    def test_single_entry_grouped(self):
        truth = {"h1": self._entry("chr1", "HAP1", 0)}
        groups = build_chrom_groups(truth)
        assert groups[("chr1", "HAP1")] == ["h1"]

    def test_separate_bins_separate_groups(self):
        truth = {
            "h1": self._entry("chr1", "HAP1", 0),
            "h2": self._entry("chr1", "HAP2", 0),
        }
        groups = build_chrom_groups(truth)
        assert "h1" in groups[("chr1", "HAP1")]
        assert "h2" in groups[("chr1", "HAP2")]

    def test_sorted_by_tstart(self):
        truth = {
            "h1": self._entry("chr1", "HAP1", tstart=50000),
            "h2": self._entry("chr1", "HAP1", tstart=10000),
            "h3": self._entry("chr1", "HAP1", tstart=30000),
        }
        groups = build_chrom_groups(truth)
        assert groups[("chr1", "HAP1")] == ["h2", "h3", "h1"]

    def test_separate_chroms_separate_groups(self):
        truth = {
            "h1": self._entry("chr1", "HAP1", 0),
            "h2": self._entry("chr2", "HAP1", 0),
        }
        groups = build_chrom_groups(truth)
        assert ("chr1", "HAP1") in groups
        assert ("chr2", "HAP1") in groups
        assert groups[("chr1", "HAP1")] == ["h1"]
        assert groups[("chr2", "HAP1")] == ["h2"]


# ── score_orientation ─────────────────────────────────────────────────────────


class TestScoreOrientation:
    """Counts wrong predictions for a given HAP1 to P1/P2 mapping. Skips absent and amb haplotigs."""

    _MAPPING = {"HAP1": "P1", "HAP2": "P2"}

    def _truth_entry(self, truth_bin):
        return {"truth_bin": truth_bin, "chrom": "chr1", "tstart": 0, "tend": 1000,
                "track": "chr1_HAP1", "best_bp": 60000, "all_bp": 60000, "best_frac": 1.0}

    def test_correct_predictions_no_wrong(self):
        truth = {"h1": self._truth_entry("HAP1")}
        pred = {"h1": {"label": "P1"}}
        rate, total, wrong = score_orientation(truth, pred, self._MAPPING)
        assert wrong == 0
        assert total == 1

    def test_wrong_prediction_counted(self):
        truth = {"h1": self._truth_entry("HAP1")}
        pred = {"h1": {"label": "P2"}}
        rate, total, wrong = score_orientation(truth, pred, self._MAPPING)
        assert wrong == 1

    def test_absent_haplotig_skipped(self):
        truth = {"h1": self._truth_entry("HAP1"), "h2": self._truth_entry("HAP1")}
        pred = {"h1": {"label": "P1"}}
        _, total, _ = score_orientation(truth, pred, self._MAPPING)
        assert total == 1

    def test_amb_label_not_informative(self):
        truth = {"h1": self._truth_entry("HAP1")}
        pred = {"h1": {"label": "amb"}}
        _, total, _ = score_orientation(truth, pred, self._MAPPING)
        assert total == 0

    def test_zero_informative_rate_is_one(self):
        truth = {"h1": self._truth_entry("HAP1")}
        pred = {"h1": {"label": "amb"}}
        rate, _, _ = score_orientation(truth, pred, self._MAPPING)
        assert rate == 1.0

    def test_wrong_rate_computed(self):
        truth = {
            "h1": self._truth_entry("HAP1"),
            "h2": self._truth_entry("HAP1"),
            "h3": self._truth_entry("HAP1"),
            "h4": self._truth_entry("HAP1"),
        }
        pred = {
            "h1": {"label": "P1"}, "h2": {"label": "P1"},
            "h3": {"label": "P2"}, "h4": {"label": "P1"},
        }
        rate, total, wrong = score_orientation(truth, pred, self._MAPPING)
        assert total == 4
        assert wrong == 1
        assert rate == pytest.approx(0.25)


# ── choose_best_orientation ───────────────────────────────────────────────────


class TestChooseBestOrientation:
    """Picks the HAP1/HAP2 to P1/P2 mapping with the lower wrong-rate, tie-breaks by total count."""

    def _truth_entry(self, truth_bin):
        return {"truth_bin": truth_bin, "chrom": "chr1", "tstart": 0, "tend": 1000,
                "track": "chr1_HAP1", "best_bp": 60000, "all_bp": 60000, "best_frac": 1.0}

    def test_hap1_p1_wins_when_better(self):
        # all HAP1 haplotigs labelled P1: HAP1 to P1 mapping has 0% error
        truth = {f"h{i}": self._truth_entry("HAP1") for i in range(4)}
        pred = {f"h{i}": {"label": "P1"} for i in range(4)}
        mapping, _ = choose_best_orientation(truth, pred)
        assert mapping == {"HAP1": "P1", "HAP2": "P2"}

    def test_hap1_p2_wins_when_better(self):
        # all HAP1 haplotigs labelled P2: HAP1 to P2 mapping has 0% error
        truth = {f"h{i}": self._truth_entry("HAP1") for i in range(4)}
        pred = {f"h{i}": {"label": "P2"} for i in range(4)}
        mapping, _ = choose_best_orientation(truth, pred)
        assert mapping == {"HAP1": "P2", "HAP2": "P1"}

    def test_tie_broken_by_higher_total(self):
        # equal rate, one candidate scores on more haplotigs, it wins
        truth = {
            "h1": self._truth_entry("HAP1"), "h2": self._truth_entry("HAP1"),
            "h3": self._truth_entry("HAP2"), "h4": self._truth_entry("HAP2"),
        }
        pred = {
            "h1": {"label": "P1"}, "h2": {"label": "P2"},
            "h3": {"label": "P1"}, "h4": {"label": "P2"},
        }
        # both mappings have 2 wrong / 4, rate 0.5, total equal, first candidate kept
        mapping, total = choose_best_orientation(truth, pred)
        assert total == 4

    def test_all_amb_returns_a_mapping(self):
        truth = {"h1": self._truth_entry("HAP1")}
        pred = {"h1": {"label": "amb"}}
        mapping, total = choose_best_orientation(truth, pred)
        assert mapping is not None
        assert total == 0


# ── chrom_metrics ─────────────────────────────────────────────────────────────


class TestChromMetrics:
    """Computes per-chrom Hamming and switch % from ordered haplotig predictions."""

    _MAPPING = {"HAP1": "P1", "HAP2": "P2"}

    def _truth_entry(self, chrom, truth_bin, tstart):
        return {"chrom": chrom, "truth_bin": truth_bin, "tstart": tstart, "tend": tstart + 1000,
                "track": f"{chrom}_{truth_bin}", "best_bp": 60000, "all_bp": 60000, "best_frac": 1.0}

    def test_all_correct_hamming_zero(self):
        truth = {"h1": self._truth_entry("chr1", "HAP1", 0),
                 "h2": self._truth_entry("chr1", "HAP1", 1000)}
        groups = build_chrom_groups(truth)
        pred = {"h1": {"label": "P1"}, "h2": {"label": "P1"}}
        rows, overall = chrom_metrics(groups, truth, pred, self._MAPPING)
        assert overall["hamming_%"] == pytest.approx(0.0)

    def test_all_correct_switch_zero(self):
        truth = {"h1": self._truth_entry("chr1", "HAP1", 0),
                 "h2": self._truth_entry("chr1", "HAP1", 1000),
                 "h3": self._truth_entry("chr1", "HAP1", 2000)}
        groups = build_chrom_groups(truth)
        pred = {"h1": {"label": "P1"}, "h2": {"label": "P1"}, "h3": {"label": "P1"}}
        _, overall = chrom_metrics(groups, truth, pred, self._MAPPING)
        assert overall["switch_%"] == pytest.approx(0.0)

    def test_wrong_haplotig_counted_in_hamming(self):
        truth = {"h1": self._truth_entry("chr1", "HAP1", 0),
                 "h2": self._truth_entry("chr1", "HAP1", 1000)}
        groups = build_chrom_groups(truth)
        pred = {"h1": {"label": "P1"}, "h2": {"label": "P2"}}
        rows, overall = chrom_metrics(groups, truth, pred, self._MAPPING)
        assert overall["wrong"] == 1
        assert overall["hamming_%"] == pytest.approx(50.0)

    def test_switch_counted_correctly(self):
        # P1 P2 P1: 2 switches / 2 boundaries = 100%
        truth = {"h1": self._truth_entry("chr1", "HAP1", 0),
                 "h2": self._truth_entry("chr1", "HAP1", 1000),
                 "h3": self._truth_entry("chr1", "HAP1", 2000)}
        groups = build_chrom_groups(truth)
        pred = {"h1": {"label": "P1"}, "h2": {"label": "P2"}, "h3": {"label": "P1"}}
        _, overall = chrom_metrics(groups, truth, pred, self._MAPPING)
        assert overall["switch_%"] == pytest.approx(100.0)

    def test_absent_haplotig_skipped(self):
        truth = {"h1": self._truth_entry("chr1", "HAP1", 0),
                 "h2": self._truth_entry("chr1", "HAP1", 1000)}
        groups = build_chrom_groups(truth)
        pred = {"h1": {"label": "P1"}}
        rows, overall = chrom_metrics(groups, truth, pred, self._MAPPING)
        assert overall["total"] == 1

    def test_no_haplotigs_in_chrom_na(self):
        truth = {"h1": self._truth_entry("chr1", "HAP1", 0)}
        groups = build_chrom_groups(truth)
        pred = {}
        rows, overall = chrom_metrics(groups, truth, pred, self._MAPPING)
        assert rows[0]["hamming_%"] == "NA"
        assert rows[0]["switch_%"] == "NA"

    def test_single_haplotig_switch_na(self):
        truth = {"h1": self._truth_entry("chr1", "HAP1", 0)}
        groups = build_chrom_groups(truth)
        pred = {"h1": {"label": "P1"}}
        rows, _ = chrom_metrics(groups, truth, pred, self._MAPPING)
        assert rows[0]["switch_%"] == "NA"

    def test_row_per_chrom_truth_bin(self):
        truth = {
            "h1": self._truth_entry("chr1", "HAP1", 0),
            "h2": self._truth_entry("chr1", "HAP2", 0),
        }
        groups = build_chrom_groups(truth)
        pred = {"h1": {"label": "P1"}, "h2": {"label": "P2"}}
        rows, _ = chrom_metrics(groups, truth, pred, self._MAPPING)
        keys = {(r["chrom"], r["truth_bin"]) for r in rows}
        assert ("chr1", "HAP1") in keys
        assert ("chr1", "HAP2") in keys


# ── validate_inputs ───────────────────────────────────────────────────────────


class TestValidateInputs:
    """Exits with a message for every invalid input combination; passes silently when valid."""

    def _touch(self, tmp_path, *names):
        paths = {}
        for name in names:
            p = tmp_path / name
            p.write_text("")
            paths[name] = str(p)
        return paths

    def test_valid_core_inputs_no_exit(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"])
        validate_inputs(args)  # must not raise

    def test_missing_idx_p1_exits(self, tmp_path):
        f = self._touch(tmp_path, "p2.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1="/nonexistent/p1.idx", idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"])
        with pytest.raises(SystemExit, match="not found"):
            validate_inputs(args)

    def test_missing_hap1_paf_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf="/nonexistent/h1.paf", hap2_paf=f["h2.paf"])
        with pytest.raises(SystemExit, match="not found"):
            validate_inputs(args)

    def test_only_truth_hap1_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf", "th1.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"],
                          truth_hap1_paf=f["th1.paf"])
        with pytest.raises(SystemExit, match="both"):
            validate_inputs(args)

    def test_only_hap1_suffix_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"],
                          hap1_suffix="_HAP1")
        with pytest.raises(SystemExit, match="both"):
            validate_inputs(args)

    def test_duplicate_paths_exit(self, tmp_path):
        f = self._touch(tmp_path, "same.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1=f["same.idx"], idx_p2=f["same.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"])
        with pytest.raises(SystemExit, match="same file"):
            validate_inputs(args)

    def test_truth_min_span_zero_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"],
                          truth_min_span=0)
        with pytest.raises(SystemExit, match="-ts"):
            validate_inputs(args)

    def test_truth_min_frac_zero_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"],
                          truth_min_best_frac=0.0)
        with pytest.raises(SystemExit, match="-tf"):
            validate_inputs(args)

    def test_min_unitig_length_zero_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"],
                          min_unitig_length=0)
        with pytest.raises(SystemExit, match="-ul"):
            validate_inputs(args)

    def test_min_total_reads_zero_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"],
                          min_total_reads=0)
        with pytest.raises(SystemExit, match="-mr"):
            validate_inputs(args)

    def test_unitig_frac_above_one_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"],
                          unitig_frac_threshold=1.1)
        with pytest.raises(SystemExit, match="-uf"):
            validate_inputs(args)

    def test_min_unitigs_zero_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"],
                          min_unitigs_per_haplotig=0)
        with pytest.raises(SystemExit, match="-mu"):
            validate_inputs(args)

    def test_hap_frac_zero_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"],
                          hap_frac_threshold=0.0)
        with pytest.raises(SystemExit, match="-hf"):
            validate_inputs(args)

    def test_missing_truth_paf_file_exits(self, tmp_path):
        f = self._touch(tmp_path, "p1.idx", "p2.idx", "h1.paf", "h2.paf", "th2.paf")
        args = _make_args(idx_p1=f["p1.idx"], idx_p2=f["p2.idx"],
                          hap1_paf=f["h1.paf"], hap2_paf=f["h2.paf"],
                          truth_hap1_paf="/nonexistent/th1.paf",
                          truth_hap2_paf=f["th2.paf"])
        with pytest.raises(SystemExit, match="not found"):
            validate_inputs(args)


# ── resolve_truth_suffixes ────────────────────────────────────────────────────


class TestResolveTruthSuffixes:
    """Returns suffix pair from user args, auto-detection, or exits on inconsistency/failure."""

    def _paf_with(self, tmp_path, name, *tnames):
        p = tmp_path / name
        lines = "".join(
            f"q\t5000\t0\t5000\t+\t{t}\t10000\t0\t5000\t5000\t5000\t60\n"
            for t in tnames
        )
        p.write_text(lines)
        return str(p)

    def test_no_truth_paf_returns_none_none(self):
        args = _make_args(truth_hap1_paf=None, truth_hap2_paf=None)
        assert resolve_truth_suffixes(args) == (None, None)

    def test_user_specified_suffixes_returned_directly(self, tmp_path):
        # user-specified suffixes bypass file scanning entirely
        h1 = self._paf_with(tmp_path, "h1.paf", "chr1_HAP1", "chr1_HAP2")
        h2 = self._paf_with(tmp_path, "h2.paf", "chr1_HAP1", "chr1_HAP2")
        args = _make_args(truth_hap1_paf=h1, truth_hap2_paf=h2,
                          hap1_suffix="_CUSTOM1", hap2_suffix="_CUSTOM2")
        s1, s2 = resolve_truth_suffixes(args)
        assert s1 == "_CUSTOM1"
        assert s2 == "_CUSTOM2"

    def test_auto_detected_suffixes_returned(self, tmp_path):
        h1 = self._paf_with(tmp_path, "h1.paf", "chr1_HAP1", "chr1_HAP2")
        h2 = self._paf_with(tmp_path, "h2.paf", "chr1_HAP1", "chr1_HAP2")
        args = _make_args(truth_hap1_paf=h1, truth_hap2_paf=h2)
        s1, s2 = resolve_truth_suffixes(args)
        assert s1 == "_HAP1"
        assert s2 == "_HAP2"

    def test_hap1_paf_detection_fails_exits(self, tmp_path):
        # unknown suffix: detect_truth_suffixes returns None, sys.exit
        h1 = self._paf_with(tmp_path, "h1.paf", "chr1_UNKNOWN")
        h2 = self._paf_with(tmp_path, "h2.paf", "chr1_HAP1", "chr1_HAP2")
        args = _make_args(truth_hap1_paf=h1, truth_hap2_paf=h2)
        with pytest.raises(SystemExit, match="hap1 PAF"):
            resolve_truth_suffixes(args)

    def test_hap2_paf_detection_fails_exits(self, tmp_path):
        h1 = self._paf_with(tmp_path, "h1.paf", "chr1_HAP1", "chr1_HAP2")
        h2 = self._paf_with(tmp_path, "h2.paf", "chr1_UNKNOWN")
        args = _make_args(truth_hap1_paf=h1, truth_hap2_paf=h2)
        with pytest.raises(SystemExit, match="hap2 PAF"):
            resolve_truth_suffixes(args)

    def test_inconsistent_suffix_pairs_exits(self, tmp_path):
        # hap1 PAF has _HAP1/_HAP2, hap2 PAF has _H1/_H2, mismatch, sys.exit
        h1 = self._paf_with(tmp_path, "h1.paf", "chr1_HAP1", "chr1_HAP2")
        h2 = self._paf_with(tmp_path, "h2.paf", "chr1_H1", "chr1_H2")
        args = _make_args(truth_hap1_paf=h1, truth_hap2_paf=h2)
        with pytest.raises(SystemExit, match="inconsistent"):
            resolve_truth_suffixes(args)


# ══════════════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS
# ══════════════════════════════════════════════════════════════════════════════

# ── synthetic data helpers ────────────────────────────────────────────────────
#
# Synthetic assembly:
#   utg1 (5000 bp): p1=800, p2=100, label P1
#   utg2 (5000 bp): p1=100, p2=800, label P2
#   utg3 (5000 bp): p1=200, p2=200, label amb (equal fracs)
#
# total_p1_reads = 1100, total_p2_reads = 1100
#
# Haplotigs:
#   hap1 / h1: [utg1, utg1, utg1, utg2]: 3 P1 + 1 P2, label P1
#   hap2 / h2: [utg2, utg2, utg2, utg1]: 3 P2 + 1 P1, label P2
#
# Truth PAFs (when used):
#   h1 maps to chr1_HAP1  (qspan 60 kb)
#   h2 maps to chr1_HAP2  (qspan 60 kb)


def _write_core_inputs(folder):
    p1_idx = os.path.join(folder, "p1.idxstats")
    p2_idx = os.path.join(folder, "p2.idxstats")
    hap1_paf = os.path.join(folder, "hap1.paf")
    hap2_paf = os.path.join(folder, "hap2.paf")

    _write_file(p1_idx,
        "utg1\t5000\t800\t10\n"
        "utg2\t5000\t100\t10\n"
        "utg3\t5000\t200\t10\n"
        "*\t0\t0\t5\n"
    )
    _write_file(p2_idx,
        "utg1\t5000\t100\t10\n"
        "utg2\t5000\t800\t10\n"
        "utg3\t5000\t200\t10\n"
        "*\t0\t0\t5\n"
    )

    def _paf_line(hap_id, qstart, utg_id):
        return f"{hap_id}\t20000\t{qstart}\t{qstart + 5000}\t+\t{utg_id}\t5000\t0\t5000\t5000\t5000\t60\n"

    _write_file(hap1_paf,
        _paf_line("h1", 0,     "utg1") +
        _paf_line("h1", 5000,  "utg1") +
        _paf_line("h1", 10000, "utg1") +
        _paf_line("h1", 15000, "utg2")
    )
    _write_file(hap2_paf,
        _paf_line("h2", 0,     "utg2") +
        _paf_line("h2", 5000,  "utg2") +
        _paf_line("h2", 10000, "utg2") +
        _paf_line("h2", 15000, "utg1")
    )
    return p1_idx, p2_idx, hap1_paf, hap2_paf


def _write_truth_pafs(folder):
    def _truth_line(qname, tname):
        return f"{qname}\t200000\t0\t200000\t+\t{tname}\t500000\t0\t200000\t200000\t200000\t60\n"

    th1 = os.path.join(folder, "truth_hap1.paf")
    th2 = os.path.join(folder, "truth_hap2.paf")
    _write_file(th1, _truth_line("h1", "chr1_HAP1"))
    _write_file(th2, _truth_line("h2", "chr1_HAP2"))
    return th1, th2


def _run_core_pipeline(folder, p1_idx, p2_idx, hap1_paf, hap2_paf,
                       min_unitig_length=2000, min_total_reads=50,
                       unitig_frac_threshold=0.60,
                       min_unitigs_per_haplotig=3, hap_frac_threshold=0.60):
    p1_idx_data, total_p1 = read_idxstats(p1_idx)
    p2_idx_data, total_p2 = read_idxstats(p2_idx)

    unitigs = build_unitig_table(
        p1_idx_data, total_p1, p2_idx_data, total_p2,
        min_unitig_length=min_unitig_length,
        min_total_reads=min_total_reads,
        frac_threshold=unitig_frac_threshold,
    )
    hap1_paths = parse_path_paf(hap1_paf)
    hap2_paths = parse_path_paf(hap2_paf)

    hap1_haplotigs = summarize_haplotigs(hap1_paths, unitigs, min_unitigs_per_haplotig, hap_frac_threshold)
    hap2_haplotigs = summarize_haplotigs(hap2_paths, unitigs, min_unitigs_per_haplotig, hap_frac_threshold)
    hap1_metrics = haplotype_struct_metrics(hap1_haplotigs)
    hap2_metrics = haplotype_struct_metrics(hap2_haplotigs)

    unitig_csv   = os.path.join(folder, "unitig_labels.csv")
    haplotig_csv = os.path.join(folder, "haplotig_labels.csv")
    summary_csv  = os.path.join(folder, "haplotype_summary.csv")

    write_unitig_table(unitigs, unitig_csv)
    write_haplotig_table(hap1_haplotigs, hap2_haplotigs, haplotig_csv)
    write_haplotype_summary([hap1_metrics], [hap2_metrics], summary_csv)

    return {
        "unitig_csv": unitig_csv,
        "haplotig_csv": haplotig_csv,
        "summary_csv": summary_csv,
        "unitigs": unitigs,
        "hap1_haplotigs": hap1_haplotigs,
        "hap2_haplotigs": hap2_haplotigs,
        "hap1_metrics": hap1_metrics,
        "hap2_metrics": hap2_metrics,
    }


def _read_csv_rows(path):
    import csv
    with open(path) as f:
        reader = csv.DictReader(f)
        return list(reader)


# ── core run fixture ──────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def core_output(tmp_path_factory):
    folder = str(tmp_path_factory.mktemp("core"))
    p1_idx, p2_idx, hap1_paf, hap2_paf = _write_core_inputs(folder)
    return _run_core_pipeline(folder, p1_idx, p2_idx, hap1_paf, hap2_paf)


# ── truth run fixture ─────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def truth_output(tmp_path_factory):
    folder = str(tmp_path_factory.mktemp("truth"))
    p1_idx, p2_idx, hap1_paf, hap2_paf = _write_core_inputs(folder)
    th1, th2 = _write_truth_pafs(folder)

    result = _run_core_pipeline(folder, p1_idx, p2_idx, hap1_paf, hap2_paf)

    hap1_suffix, hap2_suffix = "_HAP1", "_HAP2"
    truth_h1 = parse_truth_assignments(th1, min_span=50000, min_best_frac=0.60,
                                       hap1_suffix=hap1_suffix, hap2_suffix=hap2_suffix)
    truth_h2 = parse_truth_assignments(th2, min_span=50000, min_best_frac=0.60,
                                       hap1_suffix=hap1_suffix, hap2_suffix=hap2_suffix)
    truth = {**truth_h1, **truth_h2}

    groups = build_chrom_groups(truth)
    pred = {**result["hap1_haplotigs"], **result["hap2_haplotigs"]}
    mapping, _ = choose_best_orientation(truth, pred)
    chrom_rows, overall = chrom_metrics(groups, truth, pred, mapping)

    assignments_csv    = os.path.join(folder, "truth_assignments.csv")
    chrom_metrics_csv  = os.path.join(folder, "truth_chrom_metrics.csv")
    eval_summary_csv   = os.path.join(folder, "truth_eval_summary.csv")

    write_truth_assignments(truth, assignments_csv)
    write_truth_chrom_metrics(chrom_rows, chrom_metrics_csv)
    write_truth_eval_summary(mapping, overall, eval_summary_csv)

    return {
        **result,
        "truth": truth,
        "mapping": mapping,
        "chrom_rows": chrom_rows,
        "overall": overall,
        "assignments_csv": assignments_csv,
        "chrom_metrics_csv": chrom_metrics_csv,
        "eval_summary_csv": eval_summary_csv,
    }


# ── TestCoreRun ───────────────────────────────────────────────────────────────


class TestCoreRun:
    """Core pipeline (no truth): correct CSV files created with expected content."""

    def test_unitig_csv_created(self, core_output):
        assert os.path.exists(core_output["unitig_csv"])

    def test_haplotig_csv_created(self, core_output):
        assert os.path.exists(core_output["haplotig_csv"])

    def test_summary_csv_created(self, core_output):
        assert os.path.exists(core_output["summary_csv"])

    def test_unitig_csv_has_three_rows(self, core_output):
        rows = _read_csv_rows(core_output["unitig_csv"])
        assert len(rows) == 3

    def test_unitig_csv_columns(self, core_output):
        rows = _read_csv_rows(core_output["unitig_csv"])
        expected = {"unitig", "length", "p1_reads", "p2_reads", "total_reads", "p1_norm", "p2_norm", "label"}
        assert set(rows[0].keys()) == expected

    def test_utg1_labelled_p1(self, core_output):
        rows = {r["unitig"]: r for r in _read_csv_rows(core_output["unitig_csv"])}
        assert rows["utg1"]["label"] == "P1"

    def test_utg2_labelled_p2(self, core_output):
        rows = {r["unitig"]: r for r in _read_csv_rows(core_output["unitig_csv"])}
        assert rows["utg2"]["label"] == "P2"

    def test_utg3_labelled_amb(self, core_output):
        rows = {r["unitig"]: r for r in _read_csv_rows(core_output["unitig_csv"])}
        assert rows["utg3"]["label"] == "amb"

    def test_haplotig_csv_has_two_rows(self, core_output):
        rows = _read_csv_rows(core_output["haplotig_csv"])
        assert len(rows) == 2

    def test_haplotig_hap_type_values(self, core_output):
        rows = _read_csv_rows(core_output["haplotig_csv"])
        hap_types = {r["hap_type"] for r in rows}
        assert hap_types == {"hap1", "hap2"}

    def test_h1_labelled_p1(self, core_output):
        rows = {r["hap_id"]: r for r in _read_csv_rows(core_output["haplotig_csv"])}
        assert rows["h1"]["label"] == "P1"

    def test_h2_labelled_p2(self, core_output):
        rows = {r["hap_id"]: r for r in _read_csv_rows(core_output["haplotig_csv"])}
        assert rows["h2"]["label"] == "P2"

    def test_summary_csv_has_two_rows(self, core_output):
        rows = _read_csv_rows(core_output["summary_csv"])
        assert len(rows) == 2

    def test_summary_hap1_global_parent_p1(self, core_output):
        rows = {r["hap"]: r for r in _read_csv_rows(core_output["summary_csv"])}
        assert rows["hap1"]["hap_global_parent"] == "P1"

    def test_summary_hap2_global_parent_p2(self, core_output):
        rows = {r["hap"]: r for r in _read_csv_rows(core_output["summary_csv"])}
        assert rows["hap2"]["hap_global_parent"] == "P2"


# ── TestTruthRun ──────────────────────────────────────────────────────────────


class TestTruthRun:
    """Full pipeline with truth PAFs: correct truth CSV files and evaluation metrics."""

    def test_assignments_csv_created(self, truth_output):
        assert os.path.exists(truth_output["assignments_csv"])

    def test_chrom_metrics_csv_created(self, truth_output):
        assert os.path.exists(truth_output["chrom_metrics_csv"])

    def test_eval_summary_csv_created(self, truth_output):
        assert os.path.exists(truth_output["eval_summary_csv"])

    def test_both_haplotigs_assigned(self, truth_output):
        rows = _read_csv_rows(truth_output["assignments_csv"])
        assert len(rows) == 2

    def test_h1_assigned_hap1_bin(self, truth_output):
        rows = {r["hap_id"]: r for r in _read_csv_rows(truth_output["assignments_csv"])}
        assert rows["h1"]["truth_bin"] == "HAP1"

    def test_h2_assigned_hap2_bin(self, truth_output):
        rows = {r["hap_id"]: r for r in _read_csv_rows(truth_output["assignments_csv"])}
        assert rows["h2"]["truth_bin"] == "HAP2"

    def test_chrom_stripped_correctly(self, truth_output):
        rows = _read_csv_rows(truth_output["assignments_csv"])
        chroms = {r["chrom"] for r in rows}
        assert chroms == {"chr1"}

    def test_chrom_metrics_has_two_rows(self, truth_output):
        rows = _read_csv_rows(truth_output["chrom_metrics_csv"])
        assert len(rows) == 2

    def test_chrom_metrics_columns(self, truth_output):
        rows = _read_csv_rows(truth_output["chrom_metrics_csv"])
        expected = {"chrom", "truth_bin", "expected_pred", "n_haplotigs_used",
                    "wrong", "hamming_%", "switches", "boundaries", "switch_%"}
        assert set(rows[0].keys()) == expected

    def test_eval_summary_has_one_row(self, truth_output):
        rows = _read_csv_rows(truth_output["eval_summary_csv"])
        assert len(rows) == 1

    def test_orientation_mapping_consistent(self, truth_output):
        rows = _read_csv_rows(truth_output["eval_summary_csv"])
        hap1_label = rows[0]["hap1_label"]
        hap2_label = rows[0]["hap2_label"]
        assert {hap1_label, hap2_label} == {"P1", "P2"}

    def test_perfect_phasing_hamming_zero(self, truth_output):
        # h1-P1-HAP1, h2-P2-HAP2, correct orientation, 0 wrong
        rows = _read_csv_rows(truth_output["eval_summary_csv"])
        assert rows[0]["overall_wrong"] == "0"

    def test_perfect_phasing_switch_zero(self, truth_output):
        rows = _read_csv_rows(truth_output["eval_summary_csv"])
        # single haplotig per chrom per bin: no boundaries, switch = NA
        assert rows[0]["overall_switch_%"] in ("0.0000", "NA")


# ── TestEdgeCases ─────────────────────────────────────────────────────────────


class TestEdgeCases:
    """Edge inputs that must not crash: empty idxstats, empty path PAFs, all-amb unitigs."""

    def test_empty_idxstats_all_amb(self, tmp_path):
        folder = str(tmp_path)
        _write_file(os.path.join(folder, "p1.idxstats"), "")
        _write_file(os.path.join(folder, "p2.idxstats"), "")

        def _paf_line(hap_id, utg_id):
            return f"{hap_id}\t5000\t0\t5000\t+\t{utg_id}\t5000\t0\t5000\t5000\t5000\t60\n"

        _write_file(os.path.join(folder, "hap1.paf"), _paf_line("h1", "utg1"))
        _write_file(os.path.join(folder, "hap2.paf"), _paf_line("h2", "utg1"))

        result = _run_core_pipeline(
            folder,
            os.path.join(folder, "p1.idxstats"),
            os.path.join(folder, "p2.idxstats"),
            os.path.join(folder, "hap1.paf"),
            os.path.join(folder, "hap2.paf"),
        )
        assert os.path.exists(result["unitig_csv"])
        assert os.path.exists(result["haplotig_csv"])

    def test_empty_path_pafs_produce_empty_haplotig_table(self, tmp_path):
        folder = str(tmp_path)
        p1_idx, p2_idx, _, _ = _write_core_inputs(folder)
        _write_file(os.path.join(folder, "empty_hap1.paf"), "")
        _write_file(os.path.join(folder, "empty_hap2.paf"), "")

        result = _run_core_pipeline(
            folder, p1_idx, p2_idx,
            os.path.join(folder, "empty_hap1.paf"),
            os.path.join(folder, "empty_hap2.paf"),
        )
        rows = _read_csv_rows(result["haplotig_csv"])
        assert rows == []

    def test_all_unitigs_below_length_threshold_all_amb(self, tmp_path):
        folder = str(tmp_path)
        p1_idx, p2_idx, hap1_paf, hap2_paf = _write_core_inputs(folder)

        result = _run_core_pipeline(
            folder, p1_idx, p2_idx, hap1_paf, hap2_paf,
            min_unitig_length=999999,  # larger than any unitig
        )
        rows = _read_csv_rows(result["unitig_csv"])
        assert all(r["label"] == "amb" for r in rows)