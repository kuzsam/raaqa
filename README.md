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
- [License](#license)
 
---

## Installation
 
```bash
git clone https://github.com/kuzsam/raaqa.git
cd raaqa
conda env create -f environment.yml
conda activate raaqa
```
 
---


## Modules
 
| Module | Description | Status |
|---|---|---|
| `mapq_softclip.py` | Sliding-window MAPQ and soft-clip assessment | Available|
| `visualise.py` | Figures and plots from analytical module outputs | Coming soon |
| `hese.py` | Phasing quality assessment for hifiasm trio assemblies | Coming soon |
 
---


## `mapq_softclip`
 
Computes soft-clipped base percentage and mapping quality (MAPQ) in a
sliding-window framework across all contigs of a genome assembly. Produces
quality profiles at the window, contig, and genome-wide level.

### Options
 
```
usage: mapq_softclip.py [-h] [-v] -b BAM [-i BAI] [-w WINDOW] [-s STEP] [-t THREADS]
 
options:
  -h, --help            show this help message and exit
  -v, --version         show program version and exit
  -b, --bam BAM         BAM file (required)
  -i, --bai BAI         BAM index file (.bai) — defaults to BAM path + .bai
  -w, --window WINDOW   window size in kb (default: 80)
  -s, --step STEP       step size in kb (default: 20)
  -t, --threads THREADS number of parallel contig workers (default: 1)
```
 
### Example
 
```bash
python mapq_softclip.py -b file.bam -w 10 -s 2.5 -t 8
```
 
### Output
 
All output is written to a timestamped folder in the current working directory:
 
```
mapq_softclip_(window)kb_(step)kb_(date)_(time)/
    window_stats.csv          # per-window statistics across all contigs
    summary_stats.csv         # per-contig and genome-wide summary
    contigs/
        (contig name).windows.csv      # per-window statistics for each contig
        (contig name).summary.csv      # per-contig summary for each contig
        ...
```
 
#### `window_stats.csv` columns
 
| Column | Description |
|---|---|
| Chromosome | Contig name from BAM header |
| Start | Window start position (0-based) |
| End | Window end position |
| Mean_MAPQ | Base-weighted mean MAPQ across the window |
| Median_MAPQ | Base-weighted median MAPQ across the window |
| Read_Count | Number of primary reads overlapping the window |
| Total_Bases | Total aligned bases in the window |
| Softclip_Bases | Total soft-clipped bases anchored to the window |
| Softclip_% | Soft-clipped bases as a percentage of total bases |
 
#### `summary_stats.csv` columns
 
| Column | Description |
|---|---|
| Chromosome | Contig name, or `GENOME` for the genome-wide row |
| Mean_MAPQ | Base-weighted mean MAPQ across the contig or genome |
| Median_MAPQ | Base-weighted median MAPQ across the contig or genome |
| Reads_Seen | Total primary reads processed |
| Total_Bases | Total aligned bases |
| Softclip_Bases | Total soft-clipped bases |
| Softclip_% | Soft-clipped bases as a percentage of total bases |
| Windows_Created | Number of windows created for this contig or genome |

---

## `visualise` *(coming soon)*
 
Generates figures and plots from the outputs of the analytical modules.
 
---

## `hese` *(coming soon)*
 
Assesses haplotype separation quality for hifiasm trio assemblies using
Hamming and switch errors derived from parental read mappings. Applicable
specifically to hifiasm trio assemblies where parental sequencing data is
available.
 
---
 
## Roadmap
 
- **Visualisation module** — figures and plots from analytical module outputs
- **Phasing quality module** — Hamming and switch error assessment for hifiasm
  trio assemblies using parental read mappings
- Per-contig selection mode
- `--output-dir` argument
 
---
 
## License
 
MIT License — see [LICENSE](LICENSE) for details.
 