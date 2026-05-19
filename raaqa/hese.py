import argparse
import csv
import datetime
import os
import sys
from collections import defaultdict

try:
    from importlib.metadata import version
    VERSION = version("raaqa")
except Exception:
    VERSION = "unknown"

_TRUTH_SUFFIX_PAIRS = [
    ("_MATERNAL", "_PATERNAL"),
    ("_MAT", "_PAT"),
    ("_HAPLOTYPE1", "_HAPLOTYPE2"),
    ("_HAP1", "_HAP2"),
    ("_H1", "_H2")
]


# ============================================================
# CLI
# ============================================================

def parse_args():
    """Parse and return command-line arguments."""
    p = argparse.ArgumentParser(
        description="Phasing quality assessment from parental idxstats, hap path PAFs, and optional truth PAFs."
    )
    p.add_argument("-v", "--version", action="version", version=f"RAAQA: {VERSION}")

    g = p.add_argument_group("core inputs")
    g.add_argument("-ip1", "--idx-p1", required=True, metavar="FILE",
                   help="samtools idxstats for Parent1 reads mapped to assembly unitigs")
    g.add_argument("-ip2", "--idx-p2", required=True, metavar="FILE",
                   help="samtools idxstats for Parent2 reads mapped to assembly unitigs")
    g.add_argument("-h1", "--hap1-paf", required=True, metavar="FILE",
                   help="hifiasm hap1 path PAF (haplotig -> unitig path)")
    g.add_argument("-h2", "--hap2-paf", required=True, metavar="FILE",
                   help="hifiasm hap2 path PAF (haplotig -> unitig path)")

    g = p.add_argument_group("truth inputs (optional - provide both or neither)")
    g.add_argument("-th1", "--truth-hap1-paf", default=None, metavar="FILE",
                   help="MashMap PAF: hap1 haplotigs mapped to truth reference")
    g.add_argument("-th2", "--truth-hap2-paf", default=None, metavar="FILE",
                   help="MashMap PAF: hap2 haplotigs mapped to truth reference")
    g.add_argument("-s1", "--hap1-suffix", default=None, metavar="STR",
                   help="target name suffix identifying hap1 in truth PAF (e.g. MATERNAL) [default: auto-detect]")
    g.add_argument("-s2", "--hap2-suffix", default=None, metavar="STR",
                   help="target name suffix identifying hap2 in truth PAF (e.g. PATERNAL) [default: auto-detect]")

    g = p.add_argument_group("truth thresholds")
    g.add_argument("-ts", "--truth-min-span", type=int, default=50000, metavar="N",
                   dest="truth_min_span",
                   help="min query span (bp) for a truth alignment segment [default: 50000]")
    g.add_argument("-tf", "--truth-min-frac", type=float, default=0.60, metavar="FLOAT",
                   dest="truth_min_best_frac",
                   help="min fraction of total truth bp supporting the best track [default: 0.60]")

    g = p.add_argument_group("unitig thresholds")
    g.add_argument("-ul", "--min-utg-len", type=int, default=2000, metavar="N",
                   dest="min_unitig_length",
                   help="minimum unitig length (bp) to consider [default: 2000]")
    g.add_argument("-mr", "--min-reads", type=int, default=50, metavar="N",
                   dest="min_total_reads",
                   help="minimum total reads (P1+P2) on a unitig to label it [default: 50]")
    g.add_argument("-uf", "--unitig-frac", type=float, default=0.60, metavar="FLOAT",
                   dest="unitig_frac_threshold",
                   help="min library-size-normalised fraction to call a unitig P1 or P2 [default: 0.60]")

    g = p.add_argument_group("haplotig thresholds")
    g.add_argument("-mu", "--min-unitigs", type=int, default=3, metavar="N",
                   dest="min_unitigs_per_haplotig",
                   help="min informative unitigs (P1+P2) required to label a haplotig [default: 3]")
    g.add_argument("-hf", "--hap-frac", type=float, default=0.60, metavar="FLOAT",
                   dest="hap_frac_threshold",
                   help="min fraction among informative unitigs to call a haplotig P1 or P2 [default: 0.60]")

    return p.parse_args()


def validate_inputs(args):
    """Exit with an error if any required file is missing or any threshold is out of range."""
    for label, path in [
        ("Parent1 idxstats", args.idx_p1),
        ("Parent2 idxstats", args.idx_p2),
        ("Hap1 path PAF", args.hap1_paf),
        ("Hap2 path PAF", args.hap2_paf)
    ]:
        if not os.path.isfile(path):
            sys.exit(f"Error: {label} file not found: {path}")

    if bool(args.truth_hap1_paf) != bool(args.truth_hap2_paf):
        sys.exit("Error: --truth-hap1-paf and --truth-hap2-paf must both be provided or both omitted.")
    if bool(args.hap1_suffix) != bool(args.hap2_suffix):
        sys.exit("Error: -s1/--hap1-suffix and -s2/--hap2-suffix must both be provided or both omitted.")
    if (args.hap1_suffix or args.hap2_suffix) and not args.truth_hap1_paf:
        print("[WARN] -s1/-s2 suffixes provided but no truth PAFs given -- suffixes will be ignored")
    for label, path in [
        ("Truth hap1 PAF", args.truth_hap1_paf),
        ("Truth hap2 PAF", args.truth_hap2_paf)
    ]:
        if path and not os.path.isfile(path):
            sys.exit(f"Error: {label} file not found: {path}")

    input_paths = [
        ("Parent1 idxstats", args.idx_p1),
        ("Parent2 idxstats", args.idx_p2),
        ("Hap1 path PAF", args.hap1_paf),
        ("Hap2 path PAF", args.hap2_paf)
    ]
    if args.truth_hap1_paf:
        input_paths += [
            ("Truth hap1 PAF", args.truth_hap1_paf),
            ("Truth hap2 PAF", args.truth_hap2_paf)
        ]
    seen_paths = {}
    for label, path in input_paths:
        abs_path = os.path.abspath(path)
        if abs_path in seen_paths:
            sys.exit(f"Error: {label} and {seen_paths[abs_path]} point to the same file: {path}")
        seen_paths[abs_path] = label

    if args.truth_min_span <= 0:
        sys.exit(f"Error: -ts/--truth-min-span must be > 0, got {args.truth_min_span}")
    if not 0.5 < args.truth_min_best_frac <= 1.0:
        sys.exit(f"Error: -tf/--truth-min-frac must be in (0.5, 1], got {args.truth_min_best_frac}")
    if args.min_unitig_length <= 0:
        sys.exit(f"Error: -ul/--min-utg-len must be > 0, got {args.min_unitig_length}")
    if args.min_total_reads <= 0:
        sys.exit(f"Error: -mr/--min-reads must be > 0, got {args.min_total_reads}")
    if not 0.5 < args.unitig_frac_threshold <= 1.0:
        sys.exit(f"Error: -uf/--unitig-frac must be in (0.5, 1], got {args.unitig_frac_threshold}")
    if args.min_unitigs_per_haplotig <= 0:
        sys.exit(f"Error: -mu/--min-unitigs must be > 0, got {args.min_unitigs_per_haplotig}")
    if not 0.5 < args.hap_frac_threshold <= 1.0:
        sys.exit(f"Error: -hf/--hap-frac must be in (0.5, 1], got {args.hap_frac_threshold}")


# ============================================================
# Output directory
# ============================================================

def prepare_output_dir():
    """Create and return a timestamped output folder (hese_YYYYMMDD_HHMMSS)."""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    folder = f"hese_{timestamp}"
    os.makedirs(folder, exist_ok=True)
    return folder


# ============================================================
# 1) Parsing: idxstats + path PAF
# ============================================================

def read_idxstats(path):
    """Parse a samtools idxstats file and return (stats_dict, total_mapped_reads).

    stats_dict maps unitig name to {length, mapped, unmapped}. The sentinel line (*) is skipped.
    """
    stats = {}
    total_mapped = 0
    try:
        f = open(path)
    except OSError as e:
        sys.exit(f"Error: cannot open idxstats file {path!r}: {e}")
    with f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            name, length_s, mapped_s, unmapped_s = parts[:4]
            if name == "*":
                continue
            try:
                length = int(length_s)
                mapped = int(mapped_s)
                unmapped = int(unmapped_s)
            except ValueError:
                sys.exit(f"Error: malformed idxstats line {lineno} in {path!r}: {line!r}")
            if length < 0 or mapped < 0 or unmapped < 0:
                sys.exit(f"Error: negative value on idxstats line {lineno} in {path!r}: {line!r}")
            stats[name] = {"length": length, "mapped": mapped, "unmapped": unmapped}
            total_mapped += mapped
    return stats, total_mapped


def parse_path_paf(path_paf):
    """Parse a hifiasm path PAF and return {haplotig_id: [unitig_id, ...]} ordered by query start."""
    tmp = defaultdict(list)

    try:
        f = open(path_paf)
    except OSError as e:
        sys.exit(f"Error: cannot open path PAF {path_paf!r}: {e}")
    with f:
        for lineno, line in enumerate(f, 1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            hap_id = parts[0]
            try:
                qstart = int(parts[2])
            except ValueError:
                sys.exit(f"Error: malformed PAF line {lineno} in {path_paf!r}: {line!r}")
            utg_id = parts[5]
            tmp[hap_id].append((qstart, utg_id))

    hap_paths = {}
    for hap_id, entries in tmp.items():
        entries.sort(key=lambda x: x[0])
        hap_paths[hap_id] = [utg for _, utg in entries]
    return hap_paths


# ============================================================
# 2) Unitig labeling (P1 / P2 / amb)
# ============================================================

def compute_unitig_label(p1_reads, p2_reads, length, total_p1_reads, total_p2_reads,
                         min_unitig_length, min_total_reads, frac_threshold):
    """Return (label, p1_norm, p2_norm) for one unitig.

    label is P1, P2, or amb. Reads are library-size-normalised before comparing fractions.
    Returns amb if the unitig is too short, has too few reads, or falls below frac_threshold.
    """
    total_reads = p1_reads + p2_reads
    if length < min_unitig_length or total_reads < min_total_reads:
        return "amb", 0.0, 0.0

    p1_norm = (p1_reads / total_p1_reads) if total_p1_reads > 0 else 0.0
    p2_norm = (p2_reads / total_p2_reads) if total_p2_reads > 0 else 0.0

    denom = p1_norm + p2_norm
    if denom == 0.0:
        return "amb", p1_norm, p2_norm

    frac_p1 = p1_norm / denom
    frac_p2 = p2_norm / denom

    if frac_p1 >= frac_threshold and frac_p1 > frac_p2:
        return "P1", p1_norm, p2_norm
    if frac_p2 >= frac_threshold and frac_p2 > frac_p1:
        return "P2", p1_norm, p2_norm

    return "amb", p1_norm, p2_norm


def build_unitig_table(p1_idx, total_p1_reads, p2_idx, total_p2_reads,
                       min_unitig_length, min_total_reads, frac_threshold):
    """Build and return the full unitig label dict over the union of P1 and P2 unitigs."""
    unitigs = {}
    all_utgs = set(p1_idx.keys()) | set(p2_idx.keys())

    for utg in all_utgs:
        if utg in p1_idx:
            length = p1_idx[utg]["length"]
            p1_reads = p1_idx[utg]["mapped"]
        else:
            length = p2_idx[utg]["length"]
            p1_reads = 0

        p2_reads = p2_idx[utg]["mapped"] if utg in p2_idx else 0

        label, p1_norm, p2_norm = compute_unitig_label(
            p1_reads=p1_reads,
            p2_reads=p2_reads,
            length=length,
            total_p1_reads=total_p1_reads,
            total_p2_reads=total_p2_reads,
            min_unitig_length=min_unitig_length,
            min_total_reads=min_total_reads,
            frac_threshold=frac_threshold
        )

        unitigs[utg] = {
            "length": length,
            "p1_reads": p1_reads,
            "p2_reads": p2_reads,
            "total_reads": p1_reads + p2_reads,
            "p1_norm": p1_norm,
            "p2_norm": p2_norm,
            "label": label
        }

    return unitigs


# ============================================================
# 3) Haplotig labeling (from unitig labels)
# ============================================================

def summarize_haplotigs(hap_paths, unitigs, min_unitigs_per_haplotig, hap_frac_threshold):
    """Label each haplotig P1, P2, or amb based on its constituent unitig labels.

    A haplotig needs at least min_unitigs_per_haplotig informative (P1 or P2) unitigs and
    a fraction >= hap_frac_threshold for the majority label to receive that label.
    """
    haplotigs = {}

    for hap_id, utg_list in hap_paths.items():
        n_p1 = n_p2 = n_amb = 0

        for utg in utg_list:
            info = unitigs.get(utg)
            label = info["label"] if info is not None else "amb"
            if label == "P1":
                n_p1 += 1
            elif label == "P2":
                n_p2 += 1
            else:
                n_amb += 1

        n_total = len(utg_list)
        n_inf = n_p1 + n_p2

        if n_inf < min_unitigs_per_haplotig:
            hap_label = "amb"
        else:
            frac_p1 = n_p1 / float(n_inf)
            frac_p2 = n_p2 / float(n_inf)
            if frac_p1 >= hap_frac_threshold and frac_p1 > frac_p2:
                hap_label = "P1"
            elif frac_p2 >= hap_frac_threshold and frac_p2 > frac_p1:
                hap_label = "P2"
            else:
                hap_label = "amb"

        haplotigs[hap_id] = {
            "n_unitigs": n_total,
            "n_p1": n_p1,
            "n_p2": n_p2,
            "n_amb": n_amb,
            "label": hap_label
        }

    return haplotigs


# ============================================================
# 4) Haplotype structural summary (no truth)
# ============================================================

def carry_forward_amb(labels, default_label):
    """Replace amb entries by propagating the last seen P1/P2 label, or default_label if none seen yet."""
    out = []
    last = None
    for x in labels:
        if x in ("P1", "P2"):
            last = x
            out.append(x)
        else:
            out.append(last if last is not None else default_label)
    return out


def haplotype_struct_metrics(haplotigs, forced_parent=None):
    """Compute structural Hamming and switch-error rates for one haplotype assembly.

    forced_parent pins the expected majority parent (P1 or P2). If None, it is inferred
    from the majority label. Returns a dict of counts and error rates.
    """
    n_total = len(haplotigs)
    n_p1 = sum(1 for v in haplotigs.values() if v["label"] == "P1")
    n_p2 = sum(1 for v in haplotigs.values() if v["label"] == "P2")
    n_amb = n_total - n_p1 - n_p2

    if n_p1 > n_p2:
        hap_global_parent = "P1"
    elif n_p2 > n_p1:
        hap_global_parent = "P2"
    else:
        hap_global_parent = "amb"

    assigned = forced_parent if forced_parent is not None else hap_global_parent

    if assigned in ("P1", "P2") and n_total > 0:
        wrong = n_p2 if assigned == "P1" else n_p1
        hamming = 100.0 * wrong / float(n_total)
    else:
        hamming = "NA"

    if assigned in ("P1", "P2") and n_total > 1:
        ordered_ids = sorted(haplotigs.keys())
        raw_labels = [haplotigs[hid]["label"] for hid in ordered_ids]
        eff_labels = carry_forward_amb(raw_labels, default_label=assigned)

        switches = sum(1 for i in range(1, len(eff_labels)) if eff_labels[i] != eff_labels[i - 1])
        switch_pct = 100.0 * switches / float(len(eff_labels) - 1)
    else:
        switch_pct = "NA"

    return {
        "hap_global_parent": hap_global_parent,
        "assigned_parent": assigned,
        "n_haplotigs_total": n_total,
        "n_haplotigs_p1": n_p1,
        "n_haplotigs_p2": n_p2,
        "n_haplotigs_amb": n_amb,
        "hamming_struct_%": hamming,
        "switch_struct_%": switch_pct
    }


def pick_best_assignment(metrics_p1, metrics_p2):
    """Return 'P1', 'P2', or 'amb' for the best parent assignment.

    Compares n_haplotigs_p1 vs n_haplotigs_p2 counts first. If tied, falls back to
    Hamming structural % then switch structural % as secondary tiebreakers.
    """
    n_p1 = metrics_p1["n_haplotigs_p1"]
    n_p2 = metrics_p1["n_haplotigs_p2"]

    if n_p1 > n_p2:
        return "P1"
    if n_p2 > n_p1:
        return "P2"

    h1 = metrics_p1["hamming_struct_%"]
    h2 = metrics_p2["hamming_struct_%"]
    s1 = metrics_p1["switch_struct_%"]
    s2 = metrics_p2["switch_struct_%"]

    h1_ok = isinstance(h1, float)
    h2_ok = isinstance(h2, float)
    s1_ok = isinstance(s1, float)
    s2_ok = isinstance(s2, float)

    score1 = (h1 if h1_ok else 0) + (s1 if s1_ok else 0)
    score2 = (h2 if h2_ok else 0) + (s2 if s2_ok else 0)
    has_any1 = h1_ok or s1_ok
    has_any2 = h2_ok or s2_ok

    if has_any1 and has_any2:
        if score1 < score2:
            return "P1"
        if score2 < score1:
            return "P2"
    elif has_any1:
        return "P1"
    elif has_any2:
        return "P2"

    return "amb"


# ============================================================
# 5) Truth parsing
# ============================================================

def _scan_paf_tnames(paf_path):
    """Return the set of upper-cased target names (column 6) from a PAF file."""
    tnames = set()
    try:
        f = open(paf_path)
    except OSError as e:
        sys.exit(f"Error: cannot open truth PAF {paf_path!r}: {e}")
    with f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 6:
                tnames.add(parts[5].upper())
    return tnames


def detect_truth_suffixes(paf_path):
    """Scan a truth PAF and return the first matching (hap1_suffix, hap2_suffix) pair from _TRUTH_SUFFIX_PAIRS.

    Returns (None, None) if no known pair is found.
    """
    tnames = _scan_paf_tnames(paf_path)
    for s1, s2 in _TRUTH_SUFFIX_PAIRS:
        if any(t.endswith(s1) for t in tnames) and any(t.endswith(s2) for t in tnames):
            return s1, s2
    return None, None


def resolve_truth_suffixes(args):
    """Return (hap1_suffix, hap2_suffix) from user flags or auto-detection. Exits on failure."""
    if not args.truth_hap1_paf:
        return None, None
    if args.hap1_suffix and args.hap2_suffix:
        return args.hap1_suffix, args.hap2_suffix
    known = ", ".join(f"{s1}/{s2}" for s1, s2 in _TRUTH_SUFFIX_PAIRS)
    hap1_suffix, hap2_suffix = detect_truth_suffixes(args.truth_hap1_paf)
    if hap1_suffix is None:
        sys.exit(
            f"Error: could not auto-detect haplotype suffix pairing from truth hap1 PAF.\n"
            f"  Scanned: {args.truth_hap1_paf}\n"
            f"  Known pairs: {known}\n"
            f"  Possible causes: file is empty/malformed, or suffixes don't match any known pair.\n"
            f"  Use -s1/-s2 to set suffixes explicitly."
        )
    hap1_suffix_h2, hap2_suffix_h2 = detect_truth_suffixes(args.truth_hap2_paf)
    if hap1_suffix_h2 is None:
        sys.exit(
            f"Error: could not auto-detect haplotype suffix pairing from truth hap2 PAF.\n"
            f"  Scanned: {args.truth_hap2_paf}\n"
            f"  Known pairs: {known}\n"
            f"  Possible causes: file is empty/malformed, or suffixes don't match any known pair.\n"
            f"  Use -s1/-s2 to set suffixes explicitly."
        )
    if (hap1_suffix_h2, hap2_suffix_h2) != (hap1_suffix, hap2_suffix):
        sys.exit(
            f"Error: truth hap1 and hap2 PAFs have inconsistent haplotype suffix pairings.\n"
            f"  Hap1 PAF detected: {hap1_suffix}/{hap2_suffix}  ({args.truth_hap1_paf})\n"
            f"  Hap2 PAF detected: {hap1_suffix_h2}/{hap2_suffix_h2}  ({args.truth_hap2_paf})\n"
            f"  Check that both truth PAFs were aligned to the same reference.\n"
            f"  Use -s1/-s2 to override if you are sure the files are correct."
        )
    return hap1_suffix, hap2_suffix


def parse_truth_assignments(paf_path, min_span, min_best_frac, hap1_suffix, hap2_suffix):
    """Parse a MashMap truth PAF and return {haplotig_id: truth_assignment_dict}.

    Alignments shorter than min_span are skipped. A haplotig is assigned only if its
    best-matching track accounts for >= min_best_frac of its total aligned bases and
    the track name ends with a known haplotype suffix.
    """
    agg = defaultdict(dict)

    try:
        f = open(paf_path)
    except OSError as e:
        sys.exit(f"Error: cannot open truth PAF {paf_path!r}: {e}")
    with f:
        for lineno, line in enumerate(f, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            qname = parts[0]
            try:
                qstart, qend = int(parts[2]), int(parts[3])
                tstart, tend = int(parts[7]), int(parts[8])
            except ValueError:
                sys.exit(f"Error: malformed PAF line {lineno} in {paf_path!r}: {line!r}")
            tname = parts[5]

            qspan = qend - qstart
            if qspan < min_span:
                continue

            if tname not in agg[qname]:
                agg[qname][tname] = [0, tstart, tend]

            agg[qname][tname][0] += qspan
            agg[qname][tname][1] = min(agg[qname][tname][1], tstart)
            agg[qname][tname][2] = max(agg[qname][tname][2], tend)

    truth = {}
    for hap_id, tracks in agg.items():
        all_bp = sum(v[0] for v in tracks.values())
        best_track, (best_bp, best_tstart, best_tend) = max(tracks.items(), key=lambda kv: kv[1][0])
        best_frac = best_bp / float(all_bp)
        if best_frac < min_best_frac:
            continue

        track_upper = best_track.upper()
        hap1_upper = hap1_suffix.upper()
        hap2_upper = hap2_suffix.upper()
        if track_upper.endswith(hap1_upper):
            chrom = best_track[:-len(hap1_suffix)]
            truth_bin = "HAP1"
        elif track_upper.endswith(hap2_upper):
            chrom = best_track[:-len(hap2_suffix)]
            truth_bin = "HAP2"
        else:
            continue

        truth[hap_id] = {
            "chrom": chrom,
            "truth_bin": truth_bin,
            "track": best_track,
            "tstart": best_tstart,
            "tend": best_tend,
            "best_bp": best_bp,
            "all_bp": all_bp,
            "best_frac": best_frac
        }

    return truth


def build_chrom_groups(truth):
    """Group truth-assigned haplotig IDs by (chrom, truth_bin), ordered by tstart."""
    groups = defaultdict(list)
    for hid, t in truth.items():
        key = (t["chrom"], t["truth_bin"])
        groups[key].append((t["tstart"], hid))

    return {
        key: [hid for _, hid in sorted(items)]
        for key, items in groups.items()
    }


# ============================================================
# 6) Truth evaluation: pick HAP1/HAP2 -> P1/P2 mapping + chrom metrics
# ============================================================

def score_orientation(truth, pred_haplotigs, mapping):
    """Return (wrong_rate, n_informative, n_wrong) for one HAP1/HAP2 -> P1/P2 mapping."""
    n_wrong = 0
    n_informative = 0
    for hid, t in truth.items():
        if hid not in pred_haplotigs:
            continue
        pred_label = pred_haplotigs[hid]["label"]
        if pred_label not in ("P1", "P2"):
            continue
        n_informative += 1
        if pred_label != mapping[t["truth_bin"]]:
            n_wrong += 1
    wrong_rate = (n_wrong / n_informative) if n_informative else 1.0
    return wrong_rate, n_informative, n_wrong


def choose_best_orientation(truth, pred_haplotigs):
    """Pick the HAP1/HAP2 -> P1/P2 orientation that minimises wrong-label rate. Returns (mapping, n_informative)."""
    candidates = [
        {"HAP1": "P1", "HAP2": "P2"},
        {"HAP1": "P2", "HAP2": "P1"}
    ]

    best_mapping = None
    best_rate = None
    best_total = None

    for mapping in candidates:
        rate, total, _ = score_orientation(truth, pred_haplotigs, mapping)
        if best_mapping is None or rate < best_rate or (rate == best_rate and total > best_total):
            best_mapping = mapping
            best_rate = rate
            best_total = total

    return best_mapping, best_total


def chrom_metrics(groups, truth, pred_haplotigs, mapping):
    """Compute per-chromosome Hamming and switch-error rates. Returns (rows, overall_dict)."""
    rows = []
    tot_wrong = tot_total = tot_switches = tot_boundaries = 0

    for (chrom, truth_bin), hap_ids in groups.items():
        expected = mapping[truth_bin]

        pred_seq = []
        n_total = 0
        n_wrong = 0

        for hid in hap_ids:
            if hid not in pred_haplotigs:
                continue
            pred_label = pred_haplotigs[hid]["label"]
            if pred_label not in ("P1", "P2", "amb"):
                continue
            pred_seq.append(pred_label)
            n_total += 1
            if pred_label in ("P1", "P2") and pred_label != expected:
                n_wrong += 1

        # count switches only between consecutive informative (P1/P2) states
        eff = carry_forward_amb(pred_seq, default_label=expected)
        switches = 0
        boundaries = 0
        prev = None
        for x in eff:
            if x not in ("P1", "P2"):
                continue
            if prev is None:
                prev = x
                continue
            boundaries += 1
            if x != prev:
                switches += 1
            prev = x

        ham_pct = (100.0 * n_wrong / n_total) if n_total else "NA"
        sw_pct = (100.0 * switches / boundaries) if boundaries else "NA"

        rows.append({
            "chrom": chrom,
            "truth_bin": truth_bin,
            "expected_pred": expected,
            "n_haplotigs_used": n_total,
            "wrong": n_wrong,
            "hamming_%": ham_pct,
            "switches": switches,
            "boundaries": boundaries,
            "switch_%": sw_pct
        })

        tot_wrong += n_wrong
        tot_total += n_total
        tot_switches += switches
        tot_boundaries += boundaries

    overall = {
        "wrong": tot_wrong,
        "total": tot_total,
        "hamming_%": (100.0 * tot_wrong / tot_total) if tot_total else "NA",
        "switches": tot_switches,
        "boundaries": tot_boundaries,
        "switch_%": (100.0 * tot_switches / tot_boundaries) if tot_boundaries else "NA"
    }
    return rows, overall


# ============================================================
# 7) Writers
# ============================================================

def _format_pct(value):
    """Format a float to 4 decimal places, or return the value unchanged if it is not a float (e.g. 'NA')."""
    return f"{value:.4f}" if isinstance(value, float) else value


def write_run_log(args, folder, timestamp, hap1_suffix=None, hap2_suffix=None):
    """Write a run_log.txt summarising all inputs and thresholds."""
    suffix_source = "user-specified" if args.hap1_suffix else "auto-detected"
    truth_section = (
        f"\n--- Truth thresholds ---\n"
        f"Truth min span (-ts)   : {args.truth_min_span} bp\n"
        f"Truth min frac (-tf)   : {args.truth_min_best_frac}\n"
        f"Hap1 suffix (-s1)      : {hap1_suffix}  ({suffix_source})\n"
        f"Hap2 suffix (-s2)      : {hap2_suffix}  ({suffix_source})\n"
    ) if args.truth_hap1_paf else ""

    content = (
        f"HESE run log\n"
        f"Started : {timestamp}\n"
        f"\n--- Inputs ---\n"
        f"Parent1 idxstats       : {args.idx_p1}\n"
        f"Parent2 idxstats       : {args.idx_p2}\n"
        f"Hap1 path PAF          : {args.hap1_paf}\n"
        f"Hap2 path PAF          : {args.hap2_paf}\n"
        f"Truth hap1 PAF         : {args.truth_hap1_paf or 'not provided'}\n"
        f"Truth hap2 PAF         : {args.truth_hap2_paf or 'not provided'}\n"
        f"\n--- Unitig thresholds ---\n"
        f"Min utg len (-ul)      : {args.min_unitig_length} bp\n"
        f"Min reads (-mr)        : {args.min_total_reads}\n"
        f"Unitig frac (-uf)      : {args.unitig_frac_threshold}\n"
        f"\n--- Haplotig thresholds ---\n"
        f"Min unitigs (-mu)      : {args.min_unitigs_per_haplotig}\n"
        f"Hap frac (-hf)         : {args.hap_frac_threshold}\n"
        f"{truth_section}"
    )

    with open(os.path.join(folder, "run_log.txt"), "w") as f:
        f.write(content)


def write_unitig_table(unitigs, out_path):
    """Write unitig labels and read counts to a CSV file."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["unitig", "length", "p1_reads", "p2_reads", "total_reads", "p1_norm", "p2_norm", "label"])
        for utg, info in sorted(unitigs.items()):
            w.writerow([
                utg,
                info["length"],
                info["p1_reads"],
                info["p2_reads"],
                info["total_reads"],
                f"{info['p1_norm']:.6e}",
                f"{info['p2_norm']:.6e}",
                info["label"]
            ])


def write_haplotig_table(hap1_haplotigs, hap2_haplotigs, out_path):
    """Write hap1 and hap2 haplotig labels and unitig counts to a CSV file."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hap_id", "hap_type", "n_unitigs", "n_p1", "n_p2", "n_amb", "label"])
        for hid, info in sorted(hap1_haplotigs.items()):
            w.writerow([
                hid,
                "hap1",
                info["n_unitigs"],
                info["n_p1"],
                info["n_p2"],
                info["n_amb"],
                info["label"]
            ])
        for hid, info in sorted(hap2_haplotigs.items()):
            w.writerow([
                hid,
                "hap2",
                info["n_unitigs"],
                info["n_p1"],
                info["n_p2"],
                info["n_amb"],
                info["label"]
            ])


def write_haplotype_summary(hap1_metrics_list, hap2_metrics_list, out_path):
    """Write structural Hamming and switch-error metrics for hap1 and hap2 to a CSV file."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hap", "hap_global_parent", "assigned_parent", "n_haplotigs_total",
                    "n_haplotigs_p1", "n_haplotigs_p2", "n_haplotigs_amb",
                    "hamming_struct_%", "switch_struct_%"])
        for hap_name, metrics_list in (("hap1", hap1_metrics_list), ("hap2", hap2_metrics_list)):
            for m in metrics_list:
                w.writerow([
                    hap_name,
                    m["hap_global_parent"],
                    m["assigned_parent"],
                    m["n_haplotigs_total"],
                    m["n_haplotigs_p1"],
                    m["n_haplotigs_p2"],
                    m["n_haplotigs_amb"],
                    m["hamming_struct_%"],
                    m["switch_struct_%"]
                ])


def write_truth_assignments(truth, out_path):
    """Write per-haplotig truth assignments (chrom, bin, track, coordinates, fractions) to a CSV file."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hap_id", "chrom", "truth_bin", "track", "tstart", "tend", "best_bp", "all_bp", "best_frac"])
        for hid, t in sorted(truth.items()):
            w.writerow([
                hid,
                t["chrom"],
                t["truth_bin"],
                t["track"],
                t["tstart"],
                t["tend"],
                t["best_bp"],
                t["all_bp"],
                f"{t['best_frac']:.6f}"
            ])


def write_truth_chrom_metrics(rows, out_path):
    """Write per-chromosome truth Hamming and switch-error rows to a CSV file."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["chrom", "truth_bin", "expected_pred", "n_haplotigs_used",
                    "wrong", "hamming_%", "switches", "boundaries", "switch_%"])
        for r in rows:
            w.writerow([
                r["chrom"],
                r["truth_bin"],
                r["expected_pred"],
                r["n_haplotigs_used"],
                r["wrong"],
                r["hamming_%"],
                r["switches"],
                r["boundaries"],
                r["switch_%"]
            ])


def write_truth_eval_summary(mapping, overall, out_path):
    """Write the orientation mapping and genome-wide truth error summary to a CSV file."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hap1_label", "hap2_label", "overall_wrong", "overall_total",
                    "overall_hamming_%", "overall_switches", "overall_boundaries", "overall_switch_%"])
        w.writerow([
            mapping["HAP1"],
            mapping["HAP2"],
            overall["wrong"],
            overall["total"],
            overall["hamming_%"],
            overall["switches"],
            overall["boundaries"],
            overall["switch_%"]
        ])


def run_analysis(args, folder, hap1_suffix, hap2_suffix):
    """Orchestrate all compute steps: idxstats loading, unitig/haplotig labelling, and optional truth evaluation."""
    # A) idxstats
    print("\n--- [1/5] Loading parental idxstats ---")

    p1_idx, total_p1_reads = read_idxstats(args.idx_p1)
    print(f"  P1: {len(p1_idx):,} unitigs, {total_p1_reads:,} mapped reads")
    if not p1_idx:
        print(f"  WARNING: no entries parsed from P1 idxstats - check file format: {args.idx_p1}")
    elif total_p1_reads == 0:
        print("  WARNING: P1 has zero mapped reads - all unitigs will be labelled amb or P2")

    p2_idx, total_p2_reads = read_idxstats(args.idx_p2)
    print(f"  P2: {len(p2_idx):,} unitigs, {total_p2_reads:,} mapped reads")
    if not p2_idx:
        print(f"  WARNING: no entries parsed from P2 idxstats - check file format: {args.idx_p2}")
    elif total_p2_reads == 0:
        print("  WARNING: P2 has zero mapped reads - all unitigs will be labelled amb or P1")

    # B) Unitig labels
    print("\n--- [2/5] Labelling unitigs ---")

    unitigs = build_unitig_table(
        p1_idx=p1_idx,
        total_p1_reads=total_p1_reads,
        p2_idx=p2_idx,
        total_p2_reads=total_p2_reads,
        min_unitig_length=args.min_unitig_length,
        min_total_reads=args.min_total_reads,
        frac_threshold=args.unitig_frac_threshold
    )

    label_counts = defaultdict(int)
    for v in unitigs.values():
        label_counts[v["label"]] += 1
    print(f"  Total unitigs : {len(unitigs):,}")
    print(f"  P1            : {label_counts['P1']:,}")
    print(f"  P2            : {label_counts['P2']:,}")
    print(f"  amb           : {label_counts['amb']:,}")
    if label_counts['P1'] == 0 and label_counts['P2'] == 0:
        print("  WARNING: no informative unitigs labelled - check -ul/--min-utg-len, -mr/--min-reads, -uf/--unitig-frac")

    # C) Path PAFs
    print("\n--- [3/5] Parsing hap path PAFs ---")

    hap1_paths = parse_path_paf(args.hap1_paf)
    print(f"  Hap1 haplotigs: {len(hap1_paths):,}")
    if not hap1_paths:
        print(f"  WARNING: no haplotigs parsed from hap1 path PAF - check file format: {args.hap1_paf}")

    hap2_paths = parse_path_paf(args.hap2_paf)
    print(f"  Hap2 haplotigs: {len(hap2_paths):,}")
    if not hap2_paths:
        print(f"  WARNING: no haplotigs parsed from hap2 path PAF - check file format: {args.hap2_paf}")

    # D) Haplotig labels
    print("\n--- [4/5] Labelling haplotigs ---")
    hap1_haplotigs = summarize_haplotigs(
        hap_paths=hap1_paths,
        unitigs=unitigs,
        min_unitigs_per_haplotig=args.min_unitigs_per_haplotig,
        hap_frac_threshold=args.hap_frac_threshold
    )
    hap2_haplotigs = summarize_haplotigs(
        hap_paths=hap2_paths,
        unitigs=unitigs,
        min_unitigs_per_haplotig=args.min_unitigs_per_haplotig,
        hap_frac_threshold=args.hap_frac_threshold
    )
    hap1_metrics_p1 = haplotype_struct_metrics(hap1_haplotigs, forced_parent="P1")
    hap1_metrics_p2 = haplotype_struct_metrics(hap1_haplotigs, forced_parent="P2")
    hap2_metrics_p1 = haplotype_struct_metrics(hap2_haplotigs, forced_parent="P1")
    hap2_metrics_p2 = haplotype_struct_metrics(hap2_haplotigs, forced_parent="P2")
    hap1_best = pick_best_assignment(hap1_metrics_p1, hap1_metrics_p2)
    hap2_best = pick_best_assignment(hap2_metrics_p1, hap2_metrics_p2)
    for m in (hap1_metrics_p1, hap1_metrics_p2):
        m["hap_global_parent"] = hap1_best
    for m in (hap2_metrics_p1, hap2_metrics_p2):
        m["hap_global_parent"] = hap2_best
    hap1_natural = hap1_metrics_p1 if hap1_best != "P2" else hap1_metrics_p2
    hap2_natural = hap2_metrics_p1 if hap2_best != "P2" else hap2_metrics_p2
    print(f"  Hap1 - global parent: {hap1_natural['hap_global_parent']}  "
          f"P1: {hap1_natural['n_haplotigs_p1']}  P2: {hap1_natural['n_haplotigs_p2']}  "
          f"amb: {hap1_natural['n_haplotigs_amb']}  "
          f"Hamming: {_format_pct(hap1_natural['hamming_struct_%'])}%  "
          f"Switch: {_format_pct(hap1_natural['switch_struct_%'])}%")
    if hap1_natural['n_haplotigs_p1'] == 0 and hap1_natural['n_haplotigs_p2'] == 0:
        print("  WARNING: no informative hap1 haplotigs labelled - check -mu/--min-unitigs, -hf/--hap-frac")
    print(f"  Hap2 - global parent: {hap2_natural['hap_global_parent']}  "
          f"P1: {hap2_natural['n_haplotigs_p1']}  P2: {hap2_natural['n_haplotigs_p2']}  "
          f"amb: {hap2_natural['n_haplotigs_amb']}  "
          f"Hamming: {_format_pct(hap2_natural['hamming_struct_%'])}%  "
          f"Switch: {_format_pct(hap2_natural['switch_struct_%'])}%")
    if hap2_natural['n_haplotigs_p1'] == 0 and hap2_natural['n_haplotigs_p2'] == 0:
        print("  WARNING: no informative hap2 haplotigs labelled - check -mu/--min-unitigs, -hf/--hap-frac")

    # E) Write core outputs
    print("\n--- [5/5] Writing outputs ---")
    try:
        write_unitig_table(unitigs, os.path.join(folder, "unitig_labels.csv"))
        write_haplotig_table(hap1_haplotigs, hap2_haplotigs, os.path.join(folder, "haplotig_labels.csv"))
        write_haplotype_summary(
            [hap1_metrics_p1, hap1_metrics_p2],
            [hap2_metrics_p1, hap2_metrics_p2],
            os.path.join(folder, "haplotype_summary.csv")
        )
    except OSError as e:
        sys.exit(f"[ERROR] Failed to write output files: {e}")
    print("  unitig_labels.csv")
    print("  haplotig_labels.csv")
    print("  haplotype_summary.csv")

    # F) Truth evaluation (optional)
    if args.truth_hap1_paf:
        print("\n--- [+] Truth evaluation ---")
        suffix_source = "user-specified" if args.hap1_suffix else "auto-detected"
        print(f"  Suffix pairing: HAP1={hap1_suffix}  HAP2={hap2_suffix}  ({suffix_source})")
        truth_h1 = parse_truth_assignments(args.truth_hap1_paf, args.truth_min_span, args.truth_min_best_frac, hap1_suffix, hap2_suffix)
        if not truth_h1:
            if args.hap1_suffix:
                print(f"  WARNING: no haplotigs assigned from truth hap1 PAF - possible causes: wrong -s1/-s2 suffixes, -ts/-tf thresholds too strict, or file malformed: {args.truth_hap1_paf}")
            else:
                print(f"  WARNING: no haplotigs assigned from truth hap1 PAF - possible causes: -ts/-tf thresholds too strict, or file partially malformed: {args.truth_hap1_paf}")
        truth_h2 = parse_truth_assignments(args.truth_hap2_paf, args.truth_min_span, args.truth_min_best_frac, hap1_suffix, hap2_suffix)
        if not truth_h2:
            if args.hap1_suffix:
                print(f"  WARNING: no haplotigs assigned from truth hap2 PAF - possible causes: wrong -s1/-s2 suffixes, -ts/-tf thresholds too strict, or file malformed: {args.truth_hap2_paf}")
            else:
                print(f"  WARNING: no haplotigs assigned from truth hap2 PAF - possible causes: -ts/-tf thresholds too strict, or file partially malformed: {args.truth_hap2_paf}")
        overlap = set(truth_h1) & set(truth_h2)
        if overlap:
            print(f"  WARNING: {len(overlap)} haplotig ID(s) appear in both truth PAFs -- hap2 assignments will overwrite hap1")
        truth = {**truth_h1, **truth_h2}
        print(f"  Truth-assigned haplotigs: {len(truth):,}  (hap1={len(truth_h1):,}, hap2={len(truth_h2):,})")

        groups = build_chrom_groups(truth)
        pred = {**hap1_haplotigs, **hap2_haplotigs}
        mapping, n_orientation_inf = choose_best_orientation(truth, pred)
        if n_orientation_inf == 0:
            print("  WARNING: no informative (P1/P2) haplotigs found among truth-assigned haplotigs - orientation is arbitrary")
        print(f"  Orientation mapping: HAP1={mapping['HAP1']}  HAP2={mapping['HAP2']}")

        chrom_rows, overall = chrom_metrics(groups, truth, pred, mapping)
        print(f"  Overall Hamming: {_format_pct(overall['hamming_%'])}%  "
              f"Switch: {_format_pct(overall['switch_%'])}%  "
              f"({overall['wrong']}/{overall['total']} wrong haplotigs)")
        n_skipped = len(truth) - overall['total']
        if n_skipped > 0:
            print(f"  Note: {n_skipped} haplotig(s) present in truth PAF but absent from path PAF - excluded from evaluation")

        try:
            write_truth_assignments(truth, os.path.join(folder, "truth_assignments.csv"))
            write_truth_chrom_metrics(chrom_rows, os.path.join(folder, "truth_chrom_metrics.csv"))
            write_truth_eval_summary(mapping, overall, os.path.join(folder, "truth_eval_summary.csv"))
        except OSError as e:
            sys.exit(f"[ERROR] Failed to write truth output files: {e}")
        print("  truth_assignments.csv")
        print("  truth_chrom_metrics.csv")
        print("  truth_eval_summary.csv")


def main():
    """Entry point: parse and validate arguments, create output directory, then run analysis."""
    args = parse_args()
    if args.hap1_suffix:
        args.hap1_suffix = args.hap1_suffix.strip() or None
    if args.hap1_suffix and not args.hap1_suffix.startswith("_"):
        args.hap1_suffix = f"_{args.hap1_suffix}"

    if args.hap2_suffix:
        args.hap2_suffix = args.hap2_suffix.strip() or None
    if args.hap2_suffix and not args.hap2_suffix.startswith("_"):
        args.hap2_suffix = f"_{args.hap2_suffix}"
    validate_inputs(args)

    hap1_suffix, hap2_suffix = resolve_truth_suffixes(args)

    start_time = datetime.datetime.now()
    try:
        folder = prepare_output_dir()
    except OSError as e:
        sys.exit(f"[ERROR] Cannot create output directory: {e}")
    try:
        write_run_log(args, folder, start_time.strftime("%Y-%m-%d %H:%M:%S"), hap1_suffix, hap2_suffix)
    except OSError as e:
        sys.exit(f"[ERROR] Failed to write run_log.txt: {e}")
    print(f"Output folder : {folder}")
    print(f"Started       : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    try:
        run_analysis(args, folder, hap1_suffix, hap2_suffix)
    except KeyboardInterrupt:
        sys.exit("\n[INFO] Interrupted.")

    end_time = datetime.datetime.now()
    print(f"\nDone. Finished in {str(end_time - start_time).split('.')[0]}")


if __name__ == "__main__":
    main()