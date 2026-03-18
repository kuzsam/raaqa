import pysam, datetime, csv, argparse, os, sys, math, re
from collections import deque
from concurrent.futures import ProcessPoolExecutor, as_completed

MAX_MAPQ = 255
HTSLIB_THREADS_PER_WORKER = 4
VERSION = "0.1.0"

# ------------------ Argument Parsing & Validation -------------------
def parse_arguments():
    parser = argparse.ArgumentParser(description="Compute MAPQ and softclip statistics in sliding windows from a BAM file.")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("-b", "--bam", required=True, help="BAM file")
    parser.add_argument("-i", "--bai", default=None, help="BAM index file (.bai)")
    parser.add_argument("-w", "--window", type=float, default=80.0, help="Window size in kb (default=80)")
    parser.add_argument("-s", "--step", type=float, default=20.0, help="Step size in kb (default=20)")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Max number of parallel contig workers/threads (default=1)")
    return parser.parse_args()

def validate_inputs(bai_file, args):
    if not os.path.exists(args.bam):
        sys.exit(f"Error: BAM file not found: {args.bam}")
    if not os.path.exists(bai_file):
        sys.exit(f"Error: BAM index not found: {bai_file}\nPlease run: samtools index {args.bam}")
    if args.window <= 0:
        sys.exit(f"Error: Window size must be greater than 0, got {args.window} kb")
    if args.step <= 0:
        sys.exit(f"Error: Step size must be greater than 0, got {args.step} kb")
    if args.threads <= 0:
        sys.exit(f"Error: Number of threads must be greater than 0, got {args.threads}")
    if args.step > args.window:
        print(f"[WARN] Step size ({args.step} kb) is larger than window size ({args.window} kb) — "
              f"regions between windows will not be covered and results will not be genome-wide complete")

# --------------------------- Output Setup ---------------------------
def _format_kb_value(value):
    if float(value).is_integer():
        return str(int(value))
    return str(value).replace('.', ',')

def _sanitise_filename(name):
    return re.sub(r'[^\w.\-]', '_', name)

def prepare_output_dirs(window, step):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    folder = f"mapq_softclip_{_format_kb_value(window)}kb_{_format_kb_value(step)}kb_{timestamp}"
    os.makedirs(folder, exist_ok=True)

    contigs_folder = os.path.join(folder, "contigs")
    os.makedirs(contigs_folder, exist_ok=True)
    window_file = os.path.join(folder, "window_stats.csv")
    summary_file = os.path.join(folder, "summary_stats.csv")

    return folder, contigs_folder, window_file, summary_file

# --------------------- Statistical Computation ----------------------
def compute_median_from_hist(hist):
    total = sum(hist)
    if total == 0:
        return math.nan
    midpoint = total / 2
    cumulative = 0
    for value, count in enumerate(hist):
        cumulative += count
        if cumulative >= midpoint:
            return value
    return math.nan

def compute_mean_from_hist(hist):
    total = sum(hist)
    if total == 0:
        return math.nan
    return sum(i * hist[i] for i in range(len(hist))) / total

def get_softclip_bases(read):
    if read.cigartuples is None:
        return 0
    return sum(length for op, length in read.cigartuples if op == 4)

# ------------------------ Window Management -------------------------
def _make_window(start, window_bp):
    return {
        "start": start,
        "end": start + window_bp,
        "read_count": 0,
        "hist": [0] * (MAX_MAPQ + 1),
        "total_bases": 0,
        "softclip_bases": 0
    }

def _create_windows(next_window_start, limit_start, chrom_len, window_bp, step_bp):
    windows = []
    max_stretch = int(window_bp * 0.3)
    chrom_end_reached = False

    while next_window_start <= limit_start and next_window_start < chrom_len:
        w_end = next_window_start + window_bp

        # Controlled stretching if last window would overshoot chromosome end
        remainder = chrom_len - w_end
        if remainder > 0 and remainder <= max_stretch:
            w_end = chrom_len
            chrom_end_reached = True

        # Prevent overshooting chromosome
        if w_end > chrom_len:
            w_end = chrom_len
            chrom_end_reached = True

        w = _make_window(next_window_start, w_end - next_window_start)
        w["end"] = w_end
        windows.append(w)

        next_window_start += step_bp

        # Stop creating more windows if chromosome end reached
        if chrom_end_reached or next_window_start >= chrom_len:
            break

    return windows, next_window_start, chrom_end_reached

def _flush_window_to_csv(w, chrom, window_writer):
    if sum(w["hist"]) == 0:
        mean_mapq = math.nan
        median_mapq = math.nan
    else:
        mean_mapq = compute_mean_from_hist(w["hist"])
        median_mapq = compute_median_from_hist(w["hist"])

    softclip_pct = (w["softclip_bases"] / w["total_bases"] * 100) if w["total_bases"] > 0 else math.nan

    window_writer.writerow([chrom, w["start"], w["end"], f"{mean_mapq:.2f}", 
                            median_mapq, w["read_count"], w["total_bases"], 
                            w["softclip_bases"], f"{softclip_pct:.5f}"])

# --------------------- Read & Contig Processing ---------------------
def _accumulate_read_into_windows(read, windows, mapq, chrom_hist):
    ct = read.cigartuples
    if not ct:
        return

    rpos = read.reference_start

    q = min(MAX_MAPQ, int(round(mapq)))

    # windows must be sorted by start index
    windows = sorted(windows, key=lambda w: w["start"])
    w_idx = 0
    n_windows = len(windows)

    for op, length in ct:
        if length <= 0:
            continue

        # M / 0 - alignment match (sequence match or mismatch)
        # = / 7 - sequence match
        # X / 8 - sequence mismatch
        # (consume query + ref)
        if op in (0, 7, 8):
            ref_start = rpos
            ref_end = rpos + length

            while w_idx < n_windows and windows[w_idx]["end"] <= ref_start:
                w_idx += 1

            w_tmp = w_idx
            while w_tmp < n_windows and windows[w_tmp]["start"] < ref_end:
                w = windows[w_tmp]

                ov_start = max(ref_start, w["start"])
                ov_end = min(ref_end, w["end"])

                if ov_start < ov_end:
                    ov_len = ov_end - ov_start

                    # MAPQ base-weighted
                    w["hist"][q] += ov_len
                    chrom_hist[q] += ov_len

                    w["total_bases"] += ov_len

                w_tmp += 1

            rpos += length

        # I / 1 - insertion to the reference 
        # (consume query only) -> assign to current rpos
        elif op == 1:
            for w in windows:
                if w["start"] <= rpos < w["end"]:
                    w["total_bases"] += length

        # D / 2 - deletion from the reference
        # N / 3 - skipped region from the reference
        # (consume ref only)
        elif op in (2, 3):
            rpos += length

        # S / 4 - soft clipping, clipped sequences present in SEQ
        # (consume query only) -> handled by anchoring below
        elif op == 4:
            pass

        # H / 5 - hard clipping, clipped sequences not present in SEQ (consume nothing)
        # P / 6 - padding, silent deletion from padded reference
        # ignore
        elif op in (5, 6):
            pass

    # softclip anchoring
    ct = read.cigartuples
 
    # left soft clip — check ct[0], if H check ct[1]
    if ct[0][0] == 4:
        left_sc = ct[0][1]
    elif ct[0][0] == 5 and len(ct) > 1 and ct[1][0] == 4:
        left_sc = ct[1][1]
    else:
        left_sc = 0
 
    # right soft clip — check ct[-1], if H check ct[-2]
    if ct[-1][0] == 4:
        right_sc = ct[-1][1]
    elif ct[-1][0] == 5 and len(ct) > 1 and ct[-2][0] == 4:
        right_sc = ct[-2][1]
    else:
        right_sc = 0
 
    if left_sc > 0:
        anchor = read.reference_start
        for w in windows:
            if w["start"] <= anchor < w["end"]:
                w["softclip_bases"] += left_sc
                w["total_bases"] += left_sc
 
    if right_sc > 0 and read.reference_end is not None \
            and read.reference_end > read.reference_start:
        anchor = read.reference_end - 1
        for w in windows:
            if w["start"] <= anchor < w["end"]:
                w["softclip_bases"] += right_sc
                w["total_bases"] += right_sc

def _compute_contig_window_stats(bam, chrom, chrom_len, window_bp, step_bp, window_writer):
    active = deque()
    next_window_start = 0
    total_windows_created = 0
    reads_seen = 0
    chrom_total_bases = 0
    chrom_softclip_bases = 0
    chrom_hist = [0] * (MAX_MAPQ + 1)
    chrom_finished = False

    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        reads_seen += 1
        mapq = read.mapping_quality

        # Whole-read totals
        softclip_bp_full = get_softclip_bases(read)
        read_len_full = read.query_length or 0
        chrom_total_bases += read_len_full
        chrom_softclip_bases += softclip_bp_full

        read_start = read.reference_start
        read_end = read.reference_end if read.reference_end is not None else (read_start + 1)

        # Flush windows that cannot overlap this read
        while active and active[0]["end"] <= read_start:
            finished = active.popleft()
            _flush_window_to_csv(finished, chrom, window_writer)

        # Create new windows if chromosome not finished
        if not chrom_finished:
            new_windows, next_window_start, chrom_finished = _create_windows(
                next_window_start, read_end - 1, chrom_len, window_bp, step_bp
            )
            total_windows_created += len(new_windows)
            active.extend(new_windows)

        # Collect overlapping windows
        overlapping = [
            w for w in active
            if read_start < w["end"] and read_end > w["start"]
        ]

        if overlapping:
            for w in overlapping:
                w["read_count"] += 1

            _accumulate_read_into_windows(read, overlapping, mapq, chrom_hist)

    # Flush remaining windows
    while active:
        finished = active.popleft()
        _flush_window_to_csv(finished, chrom, window_writer)

    return reads_seen, chrom_total_bases, chrom_softclip_bases, chrom_hist, total_windows_created

def _run_contig_worker(bam_file, bam_index, chrom, chrom_len, window_bp, step_bp, contigs_folder):
    safe_chrom = _sanitise_filename(chrom)
    win_path = os.path.join(contigs_folder, f"{safe_chrom}.windows.csv")
    sum_path = os.path.join(contigs_folder, f"{safe_chrom}.summary.csv")

    # Open BAM independently in the worker
    with pysam.AlignmentFile(bam_file, "rb", index_filename=bam_index, threads=HTSLIB_THREADS_PER_WORKER) as bam, \
         open(win_path, "w", newline="") as win_f, \
         open(sum_path, "w", newline="") as sum_f:

        window_writer = csv.writer(win_f)
        window_writer.writerow(["Chromosome", "Start", "End", "Mean_MAPQ", "Median_MAPQ",
                                "Read_Count", "Total_Bases", "Softclip_Bases", "Softclip_%"])

        summary_writer = csv.writer(sum_f)
        summary_writer.writerow(["Chromosome", "Mean_MAPQ", "Median_MAPQ", "Reads_Seen",
                                 "Total_Bases", "Softclip_Bases", "Softclip_%", "Windows_Created"])
        
        reads_seen, chrom_total_bases, chrom_softclip_bases, chrom_hist, windows_created = _compute_contig_window_stats(
            bam, chrom, chrom_len, window_bp, step_bp, window_writer 
        )

        chrom_median = compute_median_from_hist(chrom_hist)
        chrom_mean = compute_mean_from_hist(chrom_hist)
        softclip_pct = (chrom_softclip_bases / chrom_total_bases * 100) if chrom_total_bases > 0 else math.nan
        summary_writer.writerow([chrom, f"{chrom_mean:.2f}",
                                 chrom_median, reads_seen, chrom_total_bases,
                                 chrom_softclip_bases, f"{softclip_pct:.5f}", windows_created])

    return (chrom, reads_seen, chrom_mean, chrom_median, chrom_total_bases, chrom_softclip_bases, windows_created, chrom_hist)

# -------------------------- Orchestration ---------------------------
def run_analysis(bam_file, bam_index, window_kb, step_kb, requested_threads, window_file, summary_file, contigs_folder):
    window_bp = int(round(window_kb * 1000))
    step_bp   = int(round(step_kb * 1000))

    total_genome_reads = 0
    total_genome_bases = 0
    total_genome_softclip = 0
    total_windows_created_overall = 0
    global_hist = [0] * (MAX_MAPQ + 1)

    with pysam.AlignmentFile(bam_file, "rb", index_filename=bam_index) as bam:
        references = list(zip(bam.references, bam.lengths))
    
    if len(references) == 0:
        sys.exit("Error: BAM file contains no reference sequences — file may be unmapped or malformed")

    cpu_cap = max(1, os.cpu_count() or 1)
    max_workers = min(requested_threads, cpu_cap, len(references))
    if requested_threads > max_workers:
        print(f"Note: capping workers to {max_workers} (available CPUs: {cpu_cap}, contigs: {len(references)})")

    print(f"Parallel processing {len(references)} contigs with {max_workers} workers ...")

    futures = {}
    ex = ProcessPoolExecutor(max_workers=max_workers)
    try:
        for chrom, chrom_len in references:
            fut = ex.submit(_run_contig_worker, bam_file, bam_index, chrom, chrom_len,
                            window_bp, step_bp, contigs_folder)
            futures[fut] = chrom

        per_contig_stats = {}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                (chrom, reads_seen, chrom_mean, chrom_median, chrom_total_bases, chrom_softclip_bases,
                 windows_created, chrom_hist) = fut.result()
            except Exception as e:
                print(f"\n[ERROR] Worker failed on contig '{chrom}': {type(e).__name__}: {e}")
                print("[ERROR] Aborting — cancelling all remaining workers ...")
                ex.shutdown(wait=False, cancel_futures=True)
                sys.exit(1)

            per_contig_stats[chrom] = (reads_seen, chrom_mean, chrom_median, chrom_total_bases,
                                       chrom_softclip_bases, windows_created)

            total_genome_reads += reads_seen
            total_genome_bases += chrom_total_bases
            total_genome_softclip += chrom_softclip_bases
            total_windows_created_overall += windows_created

            for i in range(MAX_MAPQ + 1):
                global_hist[i] += chrom_hist[i]

            timestamp = datetime.datetime.now()
            print(f"[{timestamp.strftime('%Y-%m-%d %H:%M:%S')}] Finished {chrom} — Reads: {reads_seen:,}, Windows: {windows_created:,}")

    except KeyboardInterrupt:
        print("\n[INFO] Interrupted — cancelling pending workers and exiting ...")
        ex.shutdown(wait=False, cancel_futures=True)
        sys.exit(1)
    else:
        ex.shutdown(wait=True)

    with open(window_file, "w", newline="") as win_out:
        window_writer = csv.writer(win_out)
        window_writer.writerow(["Chromosome", "Start", "End", "Mean_MAPQ", "Median_MAPQ",
                                "Read_Count", "Total_Bases", "Softclip_Bases", "Softclip_%"])
        
        for chrom, _ in references:
            per_path = os.path.join(contigs_folder, f"{_sanitise_filename(chrom)}.windows.csv")
            if not os.path.exists(per_path):
                continue
            with open(per_path, "r", newline="") as f_in:
                r = csv.reader(f_in)
                next(r, None)  # skip header
                for row in r:
                    window_writer.writerow(row)

    genome_median_mapq = compute_median_from_hist(global_hist)
    genome_mean_mapq = compute_mean_from_hist(global_hist)
    genome_softclip_pct = (total_genome_softclip / total_genome_bases * 100) if total_genome_bases > 0 else math.nan

    with open(summary_file, "w", newline="") as sum_out:
        summary_writer = csv.writer(sum_out)
        summary_writer.writerow(["Chromosome", "Mean_MAPQ", "Median_MAPQ", "Reads_Seen",
                                 "Total_Bases", "Softclip_Bases", "Softclip_%", "Windows_Created"])
        for chrom, _ in references:
            if chrom not in per_contig_stats:
                continue
            reads_seen, chrom_mean, chrom_median, chrom_total_bases, chrom_softclip_bases, windows_created = per_contig_stats[chrom]
            softclip_pct = (chrom_softclip_bases / chrom_total_bases * 100) if chrom_total_bases > 0 else math.nan
            summary_writer.writerow([chrom, f"{chrom_mean:.2f}",
                                     chrom_median, reads_seen, chrom_total_bases,
                                     chrom_softclip_bases, f"{softclip_pct:.5f}", windows_created])

        summary_writer.writerow(["GENOME", f"{genome_mean_mapq:.2f}", genome_median_mapq, total_genome_reads,
                                 total_genome_bases, total_genome_softclip, f"{genome_softclip_pct:.5f}", 
                                 total_windows_created_overall])

    print("\nSummary completed.")
    print(f"Total reads seen: {total_genome_reads:,}")
    print(f"Total bases: {total_genome_bases:,}, Total softclip bases: {total_genome_softclip:,}, "
          f"Softclip %: {genome_softclip_pct:.5f}")
    print(f"Genome-wide mean MAPQ: {genome_mean_mapq:.2f}, median MAPQ: {genome_median_mapq}")
    print(f"Total windows created across genome: {total_windows_created_overall:,}")

def main():
    args = parse_arguments()
    bam_index = args.bai if args.bai else args.bam + ".bai"
    validate_inputs(bam_index, args)

    folder, contigs_folder, window_file, summary_file = prepare_output_dirs(args.window, args.step)
    print(f"Output folder: {folder}")
    print(f"Window size: {args.window} kb | Step size: {args.step} kb")

    start_time = datetime.datetime.now()
    print(f"[{start_time.strftime('%Y-%m-%d %H:%M:%S')}] Starting ...")

    run_analysis(args.bam, bam_index, args.window, args.step, args.threads, window_file, summary_file, contigs_folder)

    end_time = datetime.datetime.now()
    elapsed = end_time - start_time
    print(f"[{end_time.strftime('%Y-%m-%d %H:%M:%S')}] Finished in {str(elapsed).split('.')[0]}")

if __name__ == "__main__":
    main()