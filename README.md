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

> **Note:** bioconda packaging is in progress. Until available, use the install below.

```bash
git clone https://github.com/kuzsam/raaqa.git
cd raaqa
conda env create -f environment.yml
conda activate raaqa
pip install .
```

<details>
<summary>Dev install</summary>

Use `pip install -e .` instead to install in editable mode. Source changes take effect immediately without reinstalling.

</details>

---


## Modules
 
| Module | Description | Status |
|---|---|---|
| `mapq_softclip` | Sliding-window MAPQ and soft-clip assessment | Available |
| `visualise` | Figures and plots from analytical module outputs | Available (for mapq_softclip) |
| `hese` | Phasing quality assessment for hifiasm trio assemblies | Coming soon |
 
---


## `mapq_softclip`
 
Computes soft-clipped base percentage and mapping quality (MAPQ) in a
sliding-window framework across all contigs of a genome assembly. Produces
quality profiles at the window, contig, and genome-wide level.

### Options
 
```
usage: mapq_softclip [-h] [-v] -b BAM [-i BAI] [-w WINDOW] [-s STEP] [-t THREADS]
 
options:
  -h, --help            show this help message and exit
  -v, --version         show program version and exit
  -b, --bam BAM         BAM file (required)
  -i, --bai BAI         BAM index file (.bai) — defaults to BAM path + .bai
  -w, --window WINDOW   window size in kb (default: 5)
  -s, --step STEP       step size in kb (default: 2.5)
  -t, --threads THREADS number of parallel contig workers (default: 1)
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
| Mean_MAPQ | Base-weighted mean MAPQ across the window |
| Median_MAPQ | Base-weighted median MAPQ across the window |
| Read_Count | Number of primary reads overlapping the window |
| Total_Bases | Total aligned bases in the window |
| Softclip_Bases | Total soft-clipped bases anchored to the window |
| Softclip_% | Soft-clipped bases as a percentage of total bases |
| Flag | Quality flag(s) for the window — empty if no issues detected |

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
| Chromosome | Contig name, or `GENOME` for the genome-wide row |
| Mean_MAPQ | Base-weighted mean MAPQ across the contig or genome |
| Median_MAPQ | Base-weighted median MAPQ across the contig or genome |
| Reads_Seen | Total primary reads processed |
| Total_Bases | Total aligned bases |
| Softclip_Bases | Total soft-clipped bases |
| Softclip_% | Soft-clipped bases as a percentage of total bases |
| Windows_Created | Number of windows created for this contig or genome |

</details>

---

## `visualise`

Generates figures and plots from the outputs of the analytical modules. Currently supports `mapq_softclip` output.

### Options

```
usage: visualise [-h] [-v] -m MODULE -i INPUT [-f FIGURES] [-r ROLLING]

options:
  -h, --help               show this help message and exit
  -v, --version            show program version and exit
  -m, --module MODULE      which module's output to visualise: mapq_softclip
  -i, --input INPUT        path to the output folder from the module run
  -f, --figures FIGURES    figures to generate: genome | contigs | all (default: all)
  -r, --rolling ROLLING    rolling mean window in kb for rolling figures (default: 1000)
```

### Example

```bash
visualise -m mapq_softclip -i mapq_softc_output_folder/ -f all -r 500
```

### Output

All figures are written to a `figures/` subfolder inside the input directory:

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
<summary>Figure descriptions</summary>

**`genome_summary.png`** — three stacked panels, one data point per contig at the contig midpoint:
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

---

## `hese`

> **Note:** Coming soon

Assesses haplotype separation quality for hifiasm trio assemblies using
Hamming and switch errors derived from parental read mappings. Applicable
specifically to hifiasm trio assemblies where parental sequencing data is
available.
 
---
 
## Roadmap

- **Phasing quality module** - Hamming and switch error assessment for hifiasm
  trio assemblies using parental read mappings
- **Visualise hese support** - extend visualise module to support hese output
- **Bioconda package** - easy install through bioconda channel
 
---
 
## Citation

If you use RAAQA in your work, please cite it using the DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19952414.svg)](https://doi.org/10.5281/zenodo.19952414)

---

## License
 
MIT License - see [LICENSE](LICENSE) for details.
