# Changelog

All notable changes to this project will be documented in this file.

---

## [0.2.0] - unreleased

### Added
- `visualise` module — generates genome-wide summary and per-contig figures from `mapq_softclip` output
- Genome-wide 3-panel figure: coverage, MAPQ, and softclip % per contig
- Per-contig raw 4-panel figure: coverage, median MAPQ, mean MAPQ, and softclip % per window
- Per-contig rolling mean 3-panel figure: rolling mean coverage, MAPQ, and softclip %

### Fixed
- `mapq_softclip`: reads with no MAPQ value are now excluded from the MAPQ histogram
- `mapq_softclip`: NaN values in CSV outputs now written as empty fields instead of `nan` strings
- `mapq_softclip`: `ProcessPoolExecutor` prevents worker leak on unexpected exceptions
- `mapq_softclip`: malformed or truncated BAM files produce a clearer error message

### Changed
- `environment.yml` updated with `visualise` module dependencies: `matplotlib-base`, `pandas`, `numpy`

---

## [0.1.1] - unreleased

### Changed
- Output directory prefix is now derived from the BAM filename instead of a `mapq_softclip` prefix
- `environment.yml` now uses `conda-forge` and `bioconda` channels only with `nodefaults` — removes dependency conflicts caused by mixing with the `defaults` channel
- Default window size changed to 5 kb
- Default step size changed to 2.5 kb
- Validation error messages clarity
- MAPQ value of 255 is now remapped to 0 before accumulation

---

## [0.1.0] - unreleased

### Added
- `mapq_softclip` module — sliding-window MAPQ and soft-clip assessment across all contigs of a genome assembly
- Base-weighted mean and median MAPQ computed per window, per contig, and genome-wide
- Soft-clipped base percentage computed per window, per contig, and genome-wide
- Soft-clip anchoring — left clip anchored to read start, right clip anchored to read end
- Per-contig output files written to a `contigs/` subdirectory
- Genome-wide `window_stats.csv` and `summary_stats.csv` aggregated across all contigs
- Parallel contig processing via `-t/--threads`
- Controlled window stretching/clipping at chromosome ends to avoid narrow trailing windows
