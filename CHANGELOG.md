# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
this project adheres to [Semantic Versioning](https://semver.org/).

## [1.0.0] — unreleased

First public release. This is the version referenced by the accompanying
manuscript; the corresponding tag carries its own Zenodo DOI.

### Added

- `electrophy_gui.py` — point-and-click front-end wrapping the three pipelines
  as subprocesses, with an integrated helper for the group-mapping CSV format.
- Standalone Windows executable (PyInstaller, onedir), built by `build_exe.ps1`
  and distributed as a Release asset. No Python installation required.
- `--config` YAML support in `single_ap.py`, completing a uniform
  `--input` / `--output` / `--config` interface across all three pipelines.
- SHA-256 checksums of the group-mapping CSV and of the config file, recorded
  in the `Metadata` sheet, so a result can be tied to its exact inputs.
- Structured logging in `slow_depol.py`, aligning it with the other two
  pipelines: console output unchanged, timestamped trace written to
  `runinfo/pipeline.log`.
- `tests/` — pytest suite covering the AUC estimator against closed-form
  solutions, CLI/config precedence, the checksum helper, and a smoke test per
  pipeline.
- `examples/` — configuration files, a group-mapping template, and a minimal
  ABF dataset.

### Changed

- Synthetic AUC validation now sweeps **all six** signal archetypes across the
  four sampling rates (24 cases). Previously only three of the six were swept,
  yielding 12 cases plus 6 single-rate cases.
- Multi-sampling-rate figure now sizes its grid from the number of archetypes
  instead of a fixed three panels, which silently truncated the figure.

### Fixed

- Drift-tolerance threshold in the synthetic validation was applied
  inconsistently: the console label used 2 % while the Excel `pass` column used
  3 %, so the −10 mV / 120 s case (2.33 % error) printed `FAIL` while being
  recorded as `PASS`. Both now derive from a single `DRIFT_TOL_PCT` constant.

### Notes on the analysis itself

These carried over from the pre-release development and are documented here for
provenance:

- ADP detection window reduced from 200 ms to 80 ms, consistent with the
  canonical pyramidal-neuron ADP time constants; ≈ 70 % fewer false positives.
- Saturation-sweep criterion raised from 80 % to 90 % of maximum AP1 amplitude,
  because spike-amplitude decay begins earlier under 500 ms pulses than under
  the 800 ms pulses of the original protocol.
- Spike-frequency adaptation fit now uses a four-point multi-start strategy,
  preventing single-start fits from settling in unphysiological local minima.
- Slow-depolarization baseline replaced by a local sliding percentile, removing
  the drift dependence and the burst-merging artefact of the global-baseline
  version.
