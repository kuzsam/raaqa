# RAAQA - Read Alignment-based Assembly Quality Assessment

RAAQA is a command-line tool for the quantitative assessment of genome assembly quality using read alignment-based metrics.

---

## Contents

- [Installation](#installation)
- [Modules](#modules)
  - [mapq\_softclip](#mapq_softclip)
  - [visualise](#visualise)
  - [hese](#hese)
- [Roadmap](#roadmap)
- [Citation](#citation)
- [License](#license)
 
---

## Installation

**From bioconda** (only v0.2.0 - `mapq_softclip` and `visualise` modules):

```bash
conda create -n raaqa -c conda-forge -c bioconda raaqa
conda activate raaqa
```

To install into an existing environment (may conflict with existing packages):

```bash
conda install -c conda-forge -c bioconda raaqa
```

> **Note:** The bioconda package includes only the raaqa analysis dependencies (`pysam`, `matplotlib`, `pandas`, `numpy`). The preprocessing tools required to generate raaqa inputs (`samtools`, `minimap2`, `mashmap`) are **not** included and must be installed separately.

**From source:**

```bash
git clone https://github.com/kuzsam/raaqa.git
cd raaqa
conda env create -f environment.yml
conda activate raaqa
pip install .
```

> **Note:** The `environment.yml` includes both raaqa's analysis dependencies and the preprocessing tools (`samtools`, `minimap2`, `mashmap`) needed to generate inputs for the analytical modules.

---


## Modules
 
| Module | Description | Status |
|---|---|---|
| `mapq_softclip` | Sliding-window MAPQ and soft-clip assessment | Available |
| `visualise` | Figures and plots from analytical module outputs | Available |
| `hese` | Phasing quality assessment for hifiasm assemblies | Available |
 
---


## `mapq_softclip`
 
Computes soft-clipped base percentage and mapping quality (MAPQ) in a
sliding-window framework across all contigs of a genome assembly. Produces
quality profiles at the window, contig, and genome-wide level.

### Options
 
```
usage: mapq_softclip [-h] [-v] -b FILE [-i FILE] [-w FLOAT] [-s FLOAT] [-t N]
 
options:
  -h, --help            show this help message and exit
  -v, --version         show program version and exit
  -b, --bam FILE        BAM file (required)
  -i, --bai FILE        BAM index file (.bai) — defaults to BAM path + .bai
  -w, --window FLOAT    window size in kb (default: 5)
  -s, --step FLOAT      step size in kb (default: 2.5)
  -t, --threads N       number of parallel contig workers (default: 1)
```
 
### Example
 
```bash
mapq_softclip -b file.bam -w 10 -s 2.5 -t 8
```
 
### Output
 
All output is written to a timestamped folder in the current working directory:
 
```
(bam name)_mapq_softc_(window)kb_(step)kb_(date)_(time)/
    window_stats.csv          # per-window statistics across all contigs
    summary_stats.csv         # per-contig and genome-wide summary
    contigs/
        (contig name).windows.csv      # per-window statistics for each contig
        (contig name).summary.csv      # per-contig summary for each contig
        ...
```
 
<details>
<summary><code>window_stats.csv</code> columns</summary>

| Column | Description |
|---|---|
| Chromosome | Contig name from BAM header |
| Start | Window start position (0-based) |
| End | Window end position |
| Mean_MAPQ | Base-weighted mean MAPQ |
| Median_MAPQ | Base-weighted median MAPQ |
| Read_Count | Number of primary reads |
| Total_Bases | Total aligned bases |
| Softclip_Bases | Total soft-clipped bases |
| Softclip_% | Soft-clipped bases as a percentage of total bases |
| Flag | Coverage quality flag(s) |

Flags are pipe-separated when multiple conditions apply (e.g. `LOW_COVERAGE|LOW_DEPTH`):

| Flag | Condition |
|---|---|
| `NO_COVERAGE` | Zero reads overlapping the window |
| `LOW_COVERAGE` | Read count > 0 but below 5 reads |
| `LOW_DEPTH` | Base-weighted depth (Total_Bases / window size) below 5x |

</details>

<details>
<summary><code>summary_stats.csv</code> columns</summary>

| Column | Description |
|---|---|
| Chromosome | Contig name or `GENOME`\* |
| Mean_MAPQ | Base-weighted mean MAPQ |
| Median_MAPQ | Base-weighted median MAPQ |
| Reads_Seen | Total primary reads processed |
| Total_Bases | Total aligned bases |
| Softclip_Bases | Total soft-clipped bases |
| Softclip_% | Soft-clipped bases as a percentage of total bases |
| Windows_Created | Number of windows created |

\* The final row reports genome-wide totals under the label `GENOME`.

</details>

---

## `visualise`

Generates figures and plots from the outputs of the analytical modules.

### Options

```
usage: visualise [-h] [-v] -m {mapq_softclip,hese} -i DIR [-r FLOAT]

options:
  -h, --help               show this help message and exit
  -v, --version            show program version and exit
  -m, --module             which module's output to visualise: mapq_softclip, hese
  -i, --input DIR          path to the output folder from the module run
  -r, --rolling FLOAT      rolling mean window in kb for rolling figures (mapq_softclip only) (default: 1000)
```

### Example

```bash
visualise -m mapq_softclip -i mapq_softc_output_folder/ -r 500
visualise -m hese -i hese_output_folder/
```

### Output

All figures are written to a `figures/` subfolder inside the input directory.

**mapq_softclip figures:**

```
figures/
    genome_summary.png                          # genome-wide 3-panel summary
    contigs/
        raw/
            (contig name)_per_window_raw.png    # raw per-window 4-panel per contig
            ...
        rolling/
            (contig name)_per_window_rolling.png  # rolling mean 3-panel per contig
            ...
```

<details>
<summary>mapq_softclip figure descriptions</summary>

**`genome_summary.png`** — three stacked panels, one data point per contig, equally spaced on the x-axis:
- **1** — mean and median coverage per contig
- **2** — mean and median MAPQ per contig
- **3** — softclip % per contig

**`(contig)_per_window_raw.png`** — four stacked panels showing raw per-window values across the contig:
- **1** — coverage depth per window, P99 clipping with overflow markers
- **2** — median MAPQ per window
- **3** — mean MAPQ per window
- **4** — softclip % per window

**`(contig)_per_window_rolling.png`** — three stacked panels showing rolling mean values across the contig:
- **1** — rolling mean coverage
- **2** — rolling mean and median MAPQ
- **3** — rolling mean softclip %

</details>

**hese figures:**

```
figures/
    hese_label_balance.png        # P1 / amb / P2 haplotig proportion bar
    hese_label_distributions.png  # unitig and haplotig label counts side by side
    hese_signal_penetration.png   # histogram of informative unitigs per haplotig
    hese_signal_depth.png         # informative unitig count per labeled haplotig (P1 vs P2)
    hese_phasing_errors.png       # Hamming and switch error % for all hap×parent combinations
```

<details>
<summary>hese figure descriptions</summary>

**`hese_label_balance.png`** — horizontal stacked bar showing the overall P1 / amb / P2 haplotig composition across both haplotypes, with the P1:P2 ratio.

**`hese_label_distributions.png`** — two side-by-side bar charts: unitig label counts and haplotig label counts. Annotated with count and percentage.

**`hese_signal_penetration.png`** — histogram of the number of informative (non-amb) unitigs per haplotig. A dashed threshold line marks the `--min-unitigs` cutoff below which haplotigs are labelled amb.

**`hese_signal_depth.png`** — scatter plot of informative unitig count for haplotigs labelled P1 (left) and P2 (right), with per-group medians overlaid. Shows how strongly parental signal is supported within labeled haplotigs.

**`hese_phasing_errors.png`** — grouped bar chart showing Hamming and switch error % for the best parent assignment of each haplotype. If truth PAFs were provided, dashed reference lines show truth-evaluated overall error rates. Title turns red with an additional line if both haplotypes claimed the same global parent (phasing separation failure).

</details>

---

## `hese`

Assesses haplotype separation quality for hifiasm assemblies using
structural Hamming and switch errors derived from parental read mappings.
Applicable specifically to hifiasm assemblies where parental sequencing data is available.

Each unitig is labelled P1 or P2 based on library-size-normalised parental read fractions. Haplotigs inherit a label from the majority vote of their constituent unitigs. Hamming and switch error rates are then computed both structurally (from hap path PAFs alone) and against a truth reference (from optional truth PAFs).

### Options

```
usage: hese [-h] [-v] -ip1 FILE -ip2 FILE -h1 FILE -h2 FILE
            [-th1 FILE] [-th2 FILE] [-s1 STR] [-s2 STR]
            [-ts N] [-tf FLOAT] [-ul N] [-mr N] [-uf FLOAT] [-mu N] [-hf FLOAT]

core inputs:
  -ip1, --idx-p1 FILE       samtools idxstats for Parent1 reads mapped to assembly unitigs
  -ip2, --idx-p2 FILE       samtools idxstats for Parent2 reads mapped to assembly unitigs
  -h1,  --hap1-paf FILE     hifiasm hap1 path PAF (haplotig -> unitig path)
  -h2,  --hap2-paf FILE     hifiasm hap2 path PAF (haplotig -> unitig path)

truth inputs (optional — provide both or neither):
  -th1, --truth-hap1-paf FILE   MashMap PAF: hap1 haplotigs mapped to truth reference
  -th2, --truth-hap2-paf FILE   MashMap PAF: hap2 haplotigs mapped to truth reference
  -s1,  --hap1-suffix STR       target name suffix identifying hap1 in truth PAF [default: auto-detect]
  -s2,  --hap2-suffix STR       target name suffix identifying hap2 in truth PAF [default: auto-detect]

truth thresholds:
  -ts,  --truth-min-span N      min query span (bp) for a truth alignment segment [default: 50000]
  -tf,  --truth-min-frac FLOAT  min fraction of total truth bp supporting the best track [default: 0.60]

unitig thresholds:
  -ul,  --min-utg-len N         minimum unitig length (bp) to consider [default: 2000]
  -mr,  --min-reads N           minimum total reads (P1+P2) on a unitig to label it [default: 50]
  -uf,  --unitig-frac FLOAT     min library-size-normalised fraction to call a unitig P1 or P2 [default: 0.60]

haplotig thresholds:
  -mu,  --min-unitigs N         min informative unitigs (P1+P2) required to label a haplotig [default: 3]
  -hf,  --hap-frac FLOAT        min fraction among informative unitigs to call a haplotig P1 or P2 [default: 0.60]
```

### Example

Without truth PAFs (structural assessment only):

```bash
hese -ip1 parent1.idxstats -ip2 parent2.idxstats -h1 hap1.p_ctg.path.paf -h2 hap2.p_ctg.path.paf
```

With truth PAFs (full evaluation):

```bash
hese -ip1 parent1.idxstats -ip2 parent2.idxstats \
     -h1 hap1.p_ctg.path.paf -h2 hap2.p_ctg.path.paf \
     -th1 hap1_vs_truth.paf -th2 hap2_vs_truth.paf
```

### Input preparation

The truth reference must have haplotype-specific contig/chromosome names with a consistent suffix pair, e.g. `chr1_MATERNAL` / `chr1_PATERNAL`. The suffix is auto-detected from a set of known pairs (`_MATERNAL/_PATERNAL`, `_MAT/_PAT`, `_HAPLOTYPE1/_HAPLOTYPE2`, `_HAP1/_HAP2`, `_H1/_H2`). Use `-s1`/`-s2` to override if your naming differs.

### Output

All output is written to a timestamped folder in the current working directory:

```
hese_{timestamp}/
     run_log.txt
     unitig_labels.csv
     haplotig_labels.csv
     haplotype_summary.csv
     truth_assignments.csv       (only with truth PAFs)
     truth_chrom_metrics.csv     (only with truth PAFs)
     truth_eval_summary.csv      (only with truth PAFs)
```

<details>
<summary><code>unitig_labels.csv</code> columns</summary>

| Column | Description |
|---|---|
| unitig | Unitig ID from idxstats |
| length | Unitig length in bp |
| p1_reads | Raw read count from Parent1 mapped to this unitig |
| p2_reads | Raw read count from Parent2 mapped to this unitig |
| total_reads | Total reads (p1_reads + p2_reads) |
| p1_norm | Library-size-normalised P1 read density (p1_reads / total_P1_reads) |
| p2_norm | Library-size-normalised P2 read density (p2_reads / total_P2_reads) |
| label | Unitig parental assignment: `P1`, `P2`, or `amb` |

A unitig is labelled `amb` if it is shorter than `--min-utg-len`, has fewer than `--min-reads` total reads, or neither parent exceeds the `--unitig-frac` threshold.

</details>

<details>
<summary><code>haplotig_labels.csv</code> columns</summary>

| Column | Description |
|---|---|
| hap_id | Haplotig ID from the path PAF |
| hap_type | `hap1` or `hap2` |
| n_unitigs | Total unitigs in this haplotig's path |
| n_p1 | Unitigs in path labelled P1 |
| n_p2 | Unitigs in path labelled P2 |
| n_amb | Unitigs in path labelled amb |
| label | Haplotig parental assignment: `P1`, `P2`, or `amb` |

A haplotig is labelled `amb` if the number of informative (P1+P2) unitigs is below `--min-unitigs`, or if neither parent exceeds the `--hap-frac` threshold among informative unitigs.

</details>

<details>
<summary><code>haplotype_summary.csv</code> columns</summary>

Contains 4 rows, one for each combination of haplotype (hap1/hap2) × assigned parent (P1/P2).

| Column | Description |
|---|---|
| hap | Haplotype: `hap1` or `hap2` |
| hap_global_parent | Majority-vote parental assignment for this haplotype (`P1`, `P2`, or `amb`), ties broken by combined Hamming and switch error |
| assigned_parent | Parent used for error calculation in this row (`P1` or `P2`) |
| n_haplotigs_total | Total haplotigs in this haplotype |
| n_haplotigs_p1 | Haplotigs labelled P1 |
| n_haplotigs_p2 | Haplotigs labelled P2 |
| n_haplotigs_amb | Haplotigs labelled amb |
| hamming_struct_% | Structural Hamming error rate: fraction of haplotigs labelled as the wrong parent relative to `assigned_parent` |
| switch_struct_% | Structural switch error rate: fraction of consecutive haplotig boundaries that change parent label |

The two rows where `assigned_parent = hap_global_parent` represent the natural assignment.
</details>

<details>
<summary><code>truth_assignments.csv</code> columns</summary>

One row per haplotig that was successfully assigned to a truth haplotype bin.

| Column | Description |
|---|---|
| hap_id | Haplotig ID |
| chrom | Chromosome name (from truth reference name with suffix stripped) |
| truth_bin | Truth haplotype: `HAP1` or `HAP2` |
| track | Full truth target name (e.g. `chr1_MATERNAL`) |
| tstart | Leftmost truth alignment start |
| tend | Rightmost truth alignment end |
| best_bp | Total aligned bp supporting the best-matching truth track |
| all_bp | Total aligned bp across all truth tracks for this haplotig |
| best_frac | Fraction of total bp on the best track |

Haplotigs where no single track exceeds `--truth-min-frac` are excluded (ambiguous truth assignment).

</details>

<details>
<summary><code>truth_chrom_metrics.csv</code> columns</summary>

One row per chromosome × truth haplotype bin combination.

| Column | Description |
|---|---|
| chrom | Chromosome name |
| truth_bin | Truth haplotype: `HAP1` or `HAP2` |
| expected_pred | Expected parental label for this truth bin |
| n_haplotigs_used | Haplotigs used for this chrom/bin (must appear in both truth and path PAFs) |
| wrong | Haplotigs labelled as the wrong parent |
| hamming_% | Hamming error rate for this chrom/bin |
| switches | Switch count among consecutive informative haplotigs |
| boundaries | Total consecutive informative haplotig boundaries |
| switch_% | Switch error rate for this chrom/bin |

</details>

<details>
<summary><code>truth_eval_summary.csv</code> columns</summary>

Single-row genome-wide summary of truth-based evaluation.

| Column | Description |
|---|---|
| hap1_label | How HAP1 in truth was mapped to parental label (`P1` or `P2`) |
| hap2_label | How HAP2 in truth was mapped to parental label (opposite of hap1_label) |
| overall_wrong | Total haplotigs labelled as wrong parent genome-wide |
| overall_total | Total informative haplotigs used genome-wide |
| overall_hamming_% | Genome-wide Hamming error rate |
| overall_switches | Total genome-wide switch count |
| overall_boundaries | Total genome-wide informative boundaries |
| overall_switch_% | Genome-wide switch error rate |

</details>

---
 
## Roadmap

---
 
## Citation

If you use RAAQA in your work, please cite it using the DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19952414.svg)](https://doi.org/10.5281/zenodo.19952414)

---

## License
 
MIT License - see [LICENSE](LICENSE) for details.
