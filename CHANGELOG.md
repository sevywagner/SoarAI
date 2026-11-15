# Changelog

## [1.1.0] - 2026-11-13
Manuscript cited version

### Added:
- Unified reagent models
- Removed LGBM and GB training and analysis
- Added script for generating peak fitting figure for manuscript (figure 4)
- Corrected old peak lists to dedup and correct ions and centroids values (x0) to account for charge
- Replaced web app peak fitting algo with a new AI-generated version

[1.1.0]: https://github.com/sevywagner/SoarAI/releases/tag/v1.1.0

## [1.0.1] - 2026-01-18
Docker patch

### Added:
- Docker container used to run the entirety of the reproduction pipeline with preinstalled dependencies

[1.0.1]: https://github.com/sevywagner/SoarAI/releases/tag/v1.0.1

## [1.0.0] - 2026-01-15
Official release of this project

### Added:
- Finished SoarAI pipeline
- Shell scripts for using the pipeline and reproducing manuscript results
- Scripts for plotting results and analysis
- Pipeline and determinism tests
- Training Data
- REST API for inference

[1.0.0]: https://github.com/sevywagner/SoarAI/releases/tag/v1.0.0

## [0.1.0] — 2025-09-15
Pre-release of this project.

### Added:
- Pre-processing and post-processing helpers
- A program for data preperation (peak list -> model ready data)
- A program for training an XGBoost binary classifier
- A program for inference
- A simple python program for model analysis

### Known limitations
- Sample peak list still not available
- Only tested on Gentoo Linux

[0.1.0]: https://github.com/sevywagner/SoarAI/releases/tag/v0.1.0
