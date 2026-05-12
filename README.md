# Haber_PP_MEA

MATLAB analysis pipeline for DOI and ketanserin multi-electrode array recordings.

This repository is focused on code for preprocessing, analysis, and visualization. It does not include raw recordings, caches, generated outputs, manuscript drafts, literature PDFs, or private lab documents.

## Layout

- `src/`: MATLAB source for configuration, preprocessing, analysis, figures, and utilities.
- `scripts/`: command-line entrypoints for smoke tests, preprocessing, and figure generation.
- `docs/`: selected technical notes for the analysis pipeline.

Python helper scripts live under `src/` and can be run with `uv run python <path>`.

See `src/README.md` for the detailed module map and workflow.
