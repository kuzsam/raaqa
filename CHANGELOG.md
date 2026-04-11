# Changelog

All notable changes to this project will be documented in this file.

---

## [0.1.1] - unreleased

### Changed
- Output directory prefix is now derived from the BAM filename instead of a `mapq_softclip` prefix
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
