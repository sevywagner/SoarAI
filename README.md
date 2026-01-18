# SoarAI — Ion Composition Identification for TOF Mass Spectra

**Status:** v1.0.1
**Platform:** Linux only (GCC)

SoarAI is an end-to-end machine-learning pipeline for identifying plausible ion compositions in time-of-flight (TOF) mass spectrometry peaks.

Given an **m/z** value, the pipeline enumerates elemental combinations within **±50 ppm** and returns each candidate formula along with a learned **probability** that it represents the peak.

This project is part of research conducted with the **Texas A&M Meteorology Department**, and it demonstrates an end-to-end pipeline using my C++ ML/numerical library <a href="https://github.com/sevywagner/CNum">**CNum**</a>.

---

## Table of contents
- [What this repo provides](#what-this-repo-provides)
- [Dependencies](#dependencies)
- [Supported platforms](#supported-platforms)
- [Repository layout](#repository-layout)
- [Quick start](#quick-start)
  - [Run the test pipeline](#run-the-test-pipeline)
  - [Reproduce manuscript results](#reproduce-manuscript-results)
- [Run the local web app](#run-the-local-web-app)
  - [API Endpoints](#api-endpoints)
- [Pipeline scripts](#pipeline-scripts)
- [Data preparation](#data-preparation)
  - [Peak list format](#peak-list-format)
- [Runtime configuration (YAML)](#runtime-configuration-yaml)
- [Security notes](#security-notes)
- [License](#license)

---

## What this repo provides
- **Data prep → training → evaluation → plots** for ion composition identification.
- Multiple model families (per reagent ion):
  - **CNum GBModel (XGBoost-style)**
  - **Scikit-learn GBDT**
  - **FLAML LightGBM**
  - **TensorFlow neural network**
- A local **REST API** for inference and spectrum processing.

---

## Dependencies
- <a href="https://github.com/sevywagner/CNum/releases/tag/v0.2.2">**CNum v0.2.2**</a> - required for CNum GBModel training + REST inference 
- <a href="https://github.com/jbeder/yaml-cpp">**yaml-cpp**</a> - required for parsing YAML configs in C++
- <a href="https://github.com/mikefarah/yq">**yq**</a> - required for parsing YAML configs in shell scripts

> Note: The scripts also create Python virtual environments as needed (see `prepare_python_venvs.sh`).

If you prefer not to install dependencies locally, the full pipeline can also be run in a **Docker container** with all dependencies preinstalled (see **Reproduce manuscript results** below).

---

## Supported platforms
- **Linux** (Gentoo, Arch, Ubuntu) using **GCC**

---

## Repository layout
- `tests/` — validate end-to-end pipeline + determinism on a small dataset
- `reproduce/` — reproduce all models, predictions, and figures from the manuscript
  - `reproduce/docker/` — Docker-based fully reproducible pipeline
- `web_app/` — run SoarAI locally via REST APIs
- `configs/` — YAML runtime configurations (paths, run IDs, seeds, hyperparameters)
- `data/` — peak lists output by Igor Pro used as inputs used by the pipeline

All provided bash scripts build required binaries and set up required Python environments for their respective tasks.

---

## Quick start

### Run the test pipeline
Recommended first step to confirm determinism and that everything builds/runs on your machine.

```bash
# From the SoarAI repo root:
cd tests
./shell_scripts/run_test_pipeline.sh
```

- Determinism test results print to the terminal.
- Artifacts are saved to `tests/test_results/` along with a copy of the runtime config used.

---

### Reproduce manuscript results

#### Option A — Run directly on your machine

To reproduce all models + figures using local dependencies:

```bash
# From the SoarAI repo root:
cd reproduce
./shell_scripts/reproduce_all.sh
```

Artifacts are written to the **run directory** (default: `./run/reproduced_results/`) with subdirectories:
- `models/` — trained model files
- `preds/` — model predictions on the test set
- `plots/` — reproduced manuscript figures
- `data/` — split peak lists and “combo files”

**Models produced (8 total):**  
For **NH4+** and **NO+** reagent ions (4 each):
- CNum XGBoost-style GBModel
- Scikit-learn GBDT
- FLAML LightGBM
- TensorFlow neural net

> If you modify the config, preserve **seeds**, **hyperparameters**, and **deterministic settings** if you want results to match the manuscript.

---

#### Option B — Run the entire pipeline in Docker (fully reproducible)

If you want a completely self-contained environment with **all dependencies preinstalled**, you can run the full pipeline inside a Docker container.

This option:
- builds CNum inside the container
- installs all system and Python dependencies
- runs the same reproduction scripts used in Option A

```bash
# From the SoarAI repo root:
cd reproduce
./docker/reproduce_all_docker.sh
```

The container executes the full reproduction pipeline automatically.  
All output artifacts are written to mounted directories so results persist outside the container.

**Important:** the Docker workflow ultimately invokes the same `reproduce_all.sh` pipeline used in Option A, inside the containerized environment. This ensures the code path and execution order are identical between native and containerized runs, differing only in how dependencies are provisioned.

This is the recommended option if:
- you want maximal reproducibility across machines
- you do not want to install CNum, compilers, or Python packages locally

---

## Run the local web app
Start the inference APIs locally:

```bash
# From the SoarAI repo root:
cd web_app
./shell_scripts/start_api.sh
```

This starts two REST APIs:
- **NH4+ model:** `localhost:18080`
- **NO+ model:** `localhost:18081`

### API Endpoints
- `predict/` — provided by the CNum `InferenceAPI` interface
- `process-graph/` — accepts a mass spectrum, performs naive peak fitting, and assigns formulas

---

## Pipeline scripts
If you want to run the pipeline step-by-step, these scripts are the “modules”:

- `data_prep.sh` — prepare data for both NO+ and NH4+ models
- `validate_config.sh` — validate directory formatting in YAML configs
- `build_cpp_bin.sh` — build all C++ modules
- `prepare_python_venvs.sh` — create Python venvs and install dependencies
- `train.sh` — train NH4+ and NO+ CNum GBModels
- `create_py_models.sh` — train Python models (LightGBM, GBDT, neural net) for both reagent ions
- `infer.sh` — run inference on a single input peak (m/z)
- `analysis.sh` — generate figures/plots for model results

---

## Data preparation
Example real MS experiment data is provided at:
- `data/peak_lists/master.txt`

This file can be used to run the full pipeline.

### Peak list format
Peak lists are currently exported from **Igor Pro** as **tab-separated values** in a `.txt` file.

Fields (as tagged by Igor Pro):

```
def	fit	ion	x0	tag	sumFormula	x_Lo	x_Hi	d2_Ctr	d2_Lo	d2_Hi	d1_Ctr	d1_Lo	d1_Hi	d0_Ctr	d0_Lo	d0_Hi	calFac	calUnit	charge	ionizFrac	fragOf	isotopeOf
```

---

## Runtime configuration (YAML)
Runtime configurations live in `configs/` and are used for:
- path setup
- experiment tracking via `run_id`
- setting seeds + hyperparameters

To run new experiments, change `run_id`. Artifacts (and a copy of the config) will be saved to a new run folder.

---

## Security notes
⚠️ **If you expose the API beyond your local machine or trusted network, add an additional security layer.**

- The CNum inference tools use the **Crow** C++ microframework, which has had vulnerabilities historically.
- The API uses an `exec` call to run a helper binary. It is handled safely in this project, but executing external programs always carries risk.

**Recommendation:** keep this version local-only, or place it behind:
- token-based authorization
- a tunnel/reverse proxy (e.g., Cloudflare Tunnel)
- network-level access controls

---

## License
Distributed under the **MIT License**.
