# SoarAI — Ion Composition Identification for TOF Mass Spectra

**Release:** v1.1.0
**Platform:** Linux (GCC, C++20)

SoarAI is a machine-learning pipeline for identifying plausible ion compositions in time-of-flight (TOF) mass spectra. Given a centroid **m/z** and a reagent ion, it generates candidate compositions within **±35 ppm** and ranks them using a trained classifier.

The supported reagent ions are **NH₄⁺, NO⁺, and H₃O⁺**.

This project supports atmospheric chemistry research at Texas A&M University and uses [CNum](https://github.com/sevywagner/CNum), a C++ numerical and machine-learning library developed by Severin Wagner.

## Immediate plans
- Improve memory efficiency with streams
- Improve combo file generation efficiency

## Table of contents

- [What this repository provides](#what-this-repository-provides)
- [Dependencies](#dependencies)
- [Repository layout](#repository-layout)
- [Quick start](#quick-start)
- [Docker: reproduce manuscript results](#docker-reproduction)
- [Docker: run tests](#docker-tests)
- [Run the pipeline step by step](#run-the-pipeline-step-by-step)
- [Command-line inference](#command-line-inference)
- [Run the local web app](#run-the-local-web-app)
- [Data preparation](#data-preparation)
- [Runtime configuration](#runtime-configuration)
- [Reproducibility](#reproducibility)
- [AI-generated code and documentation](#ai-generated-code-and-documentation)
- [License](#license)

## What this repository provides

- Preparation of labeled candidate formulas from reference peak lists.
- Training of a **CNum XGBoost-style gradient-boosted tree model** and a **TensorFlow/Keras neural network**.
- Evaluation of both models on the combined test set and separately for each reagent ion.
- Generation of the evaluation plots, mass-defect and Van Krevelen plots, and NH₄⁺ peak fitting visual used in the manuscript.
- Command-line inference and a web app using the trained CNum model. The REST API supports single-centroid input or paired lists of m/z and average-intensity points for end-to-end formula assignment.

The current reproduction pipeline trains **two models**, each using the combined training data. Reagent-specific predictions are evaluations of those same models on the corresponding subsets. The retained LightGBM and standard GBDT code is not part of the current reproduction workflow.

## Dependencies

- [**CNum v0.2.2**](https://github.com/sevywagner/CNum/releases/tag/v0.2.2) — C++ model training and REST inference; enable `DEPLOY_TOOLS` when building CNum for the API.
- [**yaml-cpp**](https://github.com/jbeder/yaml-cpp) — YAML configuration parsing in C++.
- [**yq**](https://github.com/kislyuk/yq) — YAML configuration parsing in shell scripts.
- [**wget**](https://www.gnu.org/software/wget/) — downloads the experimental `a.h5` dataset from Zenodo.
- [**Node.js**](https://nodejs.org/) (22 or newer) and npm — only needed for running web app frontend

The scripts create the required Python virtual environments and install their dependencies. The test and reproduction pipelines can also run in Docker, as shown below.

## Repository layout

| Path | Purpose |
| --- | --- |
| `src/`, `include/` | C++ candidate generation, preprocessing, training, inference, and API code |
| `configs/soar.yaml` | Full pipeline and API configuration |
| `configs/test_config.yaml` | Small test-pipeline configuration |
| `data/peak_lists/master_filtered.txt` | Cleaned reference peak list used by the full pipeline |
| `data/peak_lists/master.txt` | Original reference material |
| `data/test_data/test_peak_list.txt` | Small reference peak list used by the tests |
| `data/a.h5` | Experimental spectra used by the NH₄⁺ peak fitting visual |
| `py_models/` | Python model training and its virtual environment |
| `analysis/` | Evaluation and figure scripts, and their virtual environment |
| `reproduce/` | Full pipeline scripts and Docker workflow |
| `tests/` | Small pipeline and repeated-training comparison |
| `web_app/` | API startup scripts, deployed model, uploads, and peak-fitting utility |
| `web_app/soarai-frontend/` | React frontend for centroid input and spectrum uploads |

**Working directories matter.** Start each example from the repository root unless it says otherwise. Relative paths in the configuration are resolved from the process's working directory. Use `reproduce/` for the full pipeline, `tests/` for the test pipeline, and `web_app/` for the API.

## Quick start

### Run the test pipeline

After installing the native dependencies:

```bash
cd tests
./shell_scripts/run_test_pipeline.sh
```

This builds the C++ programs, prepares the Python environments, splits the small reference list, and trains each model twice. It compares the exported prediction-score files from the two runs and returns a nonzero status if a comparison fails. Results and a configuration copy are saved under **`tests/test_results/`** by default.

The supplied fixture contains 30 NH₄⁺ reference peaks. It checks the pipeline and repeated training on that fixture; it does not exercise all three reagent subsets or generate the manuscript figures.

#### Docker tests

To run the same test pipeline in Docker, start from the repository root:

```bash
cd tests
./docker/run_tests_docker.sh
```

The wrapper builds the image, runs the tests, and saves results to `tests/test_results/` on the host. Docker and host-side `yq` are required.

### Reproduce manuscript results

The full pipeline uses `data/peak_lists/master_filtered.txt` and `data/a.h5`.

The `a.h5` file contains experimental TOF-MS spectra used to generate the NH₄⁺ peak fitting visual in the SoarAI manuscript. Open the [Zenodo share link](https://zenodo.org/records/24214638?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImE0YTA5M2FlLWFhZDYtNDYyOC05ZmEwLTFkMjQ3NzMzYmMzZSIsImRhdGEiOnt9LCJyYW5kb20iOiI3YjEyYTkyMTM1ZmM0OGE3ZWE2MjkyYjk5NjM1MGY0YiJ9.jxDGs9qB7Q8LRj-2zTZsvETSNQV8AN5D-eVkBHIQ9IzxBvyc6By7K8Z94K1ouxL7wCobOyekJf3HI1_Ut2SSpQ) to access the dataset and download `a.h5`.

```bash
cd reproduce
./shell_scripts/reproduce_all.sh
```

The script prepares the Python environments, builds C++, prepares the data, trains both models, and runs all four analysis scripts. It also copies the HDF5 file into the run's data directory.

The default output directory is **`reproduce/manuscript_results/`**:

| Output | Contents |
| --- | --- |
| `data/peak_lists/` | `train.txt` and `test.txt` reference-peak splits |
| `data/combo_files/` | Combined and reagent-specific candidate/label CSV files |
| `data/combo_files/unidentified_compounds/` | Reference assignments not recovered during candidate generation |
| `data/a.h5` | Copy of the spectrum file used for the peak fitting visual |
| `models/` | `xgboost.cmod` and `nn.keras` |
| `preds/` | Predictions for the combined test set and each reagent subset |
| `preds/logits/` | Continuous prediction scores used for ROC curves and repeated-run comparisons |
| `plots/` | ROC curves, confusion matrices, mass-defect and Van Krevelen plots, and the NH₄⁺ peak fitting visual |
| `plots/demonstration_figure_data/` | Data for the peak fitting visual, instrument profile, fitted-peak table, and fit details |
| `soar.yaml` | Snapshot of the configuration used by the pipeline |

Combined prediction files are named `xgboost_preds.txt` and `nn_preds.txt`. Subset files add an `NH4_`, `NO_`, or `H3O_` prefix. Full-pipeline score dumps use the suffix `_logits_final.txt`.

The analysis prints precision, recall, ROC AUC, and accuracy. Figure filenames include `ROC_AUC.png` and its reagent-prefixed variants; `XGBoost.pdf`, `NN.pdf`, and their reagent-suffixed variants; `mass_defect.pdf`; `van_krevelen.pdf`; and `fig04.pdf`/`fig04.png` for the peak fitting visual.

#### Docker reproduction

To run the full pipeline in Docker, start from the repository root:

```bash
cd reproduce
./docker/reproduce_all_docker.sh
```

This builds the image, installs the dependencies, and runs the same reproduction pipeline inside the container. Results are saved to `reproduce/manuscript_results/` on the host. Docker and host-side `yq` are required.

## Run the pipeline step by step

Run these commands in order from the repository root:

```bash
cd reproduce
./shell_scripts/validate_config.sh ../configs/soar.yaml
./shell_scripts/download_h5.sh ../configs/soar.yaml
./shell_scripts/prepare_python_venvs.sh ../configs/soar.yaml
./shell_scripts/build_cpp_bin.sh
./shell_scripts/data_prep.sh ../configs/soar.yaml
./shell_scripts/train.sh ../configs/soar.yaml
./shell_scripts/create_py_models.sh ../configs/soar.yaml
./shell_scripts/analysis.sh ../configs/soar.yaml
```

| Script | Action |
| --- | --- |
| `validate_config.sh` | Reports missing trailing slashes in configured directory strings |
| `download_h5.sh` | Retrieves the HDF5 spectrum from Zenodo for the peak fitting visual |
| `prepare_python_venvs.sh` | Creates the model-training and analysis environments and installs their requirements |
| `build_cpp_bin.sh` | Configures and builds the C++ programs in `build/` |
| `data_prep.sh` | Splits the reference peaks, generates candidate datasets, and copies the configuration into the run directory |
| `train.sh` | Trains and evaluates the combined CNum model; writes score dumps tagged `final` |
| `create_py_models.sh` | Trains and evaluates the combined neural network; writes score dumps tagged `final` |
| `analysis.sh` | Runs all analysis/figure scripts |

These stage scripts rely on the outputs of the earlier stages. The full-pipeline wrappers use `configs/soar.yaml`; passing a different filename to `reproduce_all.sh` does not switch configurations.

### Rerun individual analysis scripts

After the pipeline has created the run outputs and copied the HDF5 file:

```bash
cd reproduce
source ../analysis/bin/activate
python3 ../analysis/analysis.py ../configs/soar.yaml final
python3 ../analysis/van_krevelen.py ../configs/soar.yaml
python3 ../analysis/mass_defect.py ../configs/soar.yaml
python3 ../analysis/plot_nh4_peak_fit.py ../configs/soar.yaml
```

The `final` argument selects the score-dump suffix used by the full pipeline. The other scripts do not take that argument.

### NH₄⁺ peak fitting visual

`analysis/plot_nh4_peak_fit.py` is run by `analysis.sh` and can also be invoked separately using the analysis commands above. It requires the configured HDF5 spectrum and reference peak list.

This script **fixes the component centers at the charged reference peak-list values** and fits amplitudes, widths, and a linear baseline. It produces `fig04.pdf`, `fig04.png`, and numerical/provenance files in the configured figure-data directory. It is a constrained illustration of the assigned peaks; it does not independently recover centroids or identify formulas.

## Command-line inference

After training the CNum model:

```bash
cd reproduce
./shell_scripts/infer.sh ../configs/soar.yaml
```

Enter the centroid m/z, then explicitly enter `NH4`, `NO`, or `H3O` at the reagent prompt. The program loads `xgboost.cmod` from the configured run directory and prints candidate formulas, classification scores on a 0–1 scale, and signed ppm errors, sorted by score.

## Run the local web app

The web app consists of the C++ REST API and a React frontend. The API uses the same combined CNum model for all three supported reagents. Place the frontend directory at `web_app/soarai-frontend/`.

### Start the API

Install CNum with deployment support, then run:

```bash
cd web_app
./shell_scripts/start_api.sh
```

Startup builds the C++ API, copies the trained `xgboost.cmod` from the configured reproduction run to `paths.api.model_path`, and builds the peak-fitting executable if it is missing. With the default settings, the model is copied from `reproduce/manuscript_results/models/xgboost.cmod` to `web_app/models/xgboost.cmod`.

If no model is available in the reproduction run, place a compatible trained CNum model at `web_app/models/xgboost.cmod` before starting. Startup stops if the destination model is still missing. An available model from the configured reproduction run replaces the destination copy at startup.

The server runs in the foreground on **port 18080** by default. Stop it with `Ctrl+C`. Restart after changing the model. The default configuration keeps 30 model instances in the API pool; adjust `api.n_model_instances` if needed for available memory and concurrency.

### Start the frontend

With the API running, open a second terminal and run from the repository root:

```bash
cd web_app/soarai-frontend
npm ci
npm run dev
```

Open [http://localhost:5173/mass-spec-web-ui/](http://localhost:5173/mass-spec-web-ui/), using the port printed by Vite if 5173 is unavailable. Enter a centroid or use the graph-upload page to submit the paired m/z and average-intensity CSV files. Select NH₄⁺, NO⁺, or H₃O⁺ for either workflow.

The frontend calls `http://localhost:18080`. If the API address changes, update the fetch URLs in `src/pages/ProcessMassValue.jsx` and `src/pages/Graphs.jsx` within the frontend directory.

The build is written to `dist/`. Append `/mass-spec-web-ui/` to the preview server address and keep the API running.

### Predict candidate formulas

Send a centroid and its reagent to `POST /predict`:

```bash
curl --fail-with-body -X POST http://127.0.0.1:18080/predict \
  -H 'Content-Type: application/json' \
  -d '{"mz_array":[174.040778],"reagent_ion":"NH4"}'
```

The JSON response contains parallel arrays:

| Field | Meaning |
| --- | --- |
| `compounds` | Simplified candidate formulas |
| `uCompounds` | Candidate formulas retaining the generated adduct notation |
| `scores` | Classification scores multiplied by 100 |
| `ppms` | Signed mass errors relative to the theoretical ion mass |
| `critereaEncodings` | Four chemical-constraint flags for each candidate, in the source-defined order |

Candidates are sorted by score. Multiple input m/z values are combined into one sorted list without an input-peak index in the response; submit one centroid per request when you need an unambiguous peak-to-candidate association. Scores are model outputs and should not be interpreted as calibrated identification certainty.

### Process a spectrum

`POST /process-graph` accepts a multipart upload with three fields:

| Field | Contents |
| --- | --- |
| `reagentIon` | `NH4`, `NO`, or `H3O` |
| `base` | One-column CSV with the header `mz_base`, followed by m/z values |
| `av` | One-column CSV with the header `mz_av`, followed by corresponding intensities |

Use equal-length numeric columns, a positive and strictly increasing m/z axis, finite values, and no blank data rows. Each uploaded file must be at most **7 MiB**. Replace the example filenames with your input paths:

```bash
curl --fail-with-body -X POST http://127.0.0.1:18080/process-graph \
  -F 'reagentIon=NH4' \
  -F 'base=@mz_base.csv' \
  -F 'av=@mz_av.csv'
```

The endpoint detects and fits peaks, classifies candidates for each fitted centroid, and returns plain-text tables with `Ion`, `PPM`, and `Confidence` columns. Here, `Confidence` is the model score on a **0–1 scale**; `/predict` expresses scores on a 0–100 scale.

Uploaded spectra and detected-centroid files are stored under `web_app/uploads/graphs/` and `web_app/uploads/detected_peaks/` by default. The response is returned after fitting and inference finish. Runtime depends on spectrum size and overlapping peak groups; dense spectra can take substantially longer than isolated-peak examples.

### API access

The supplied server has no authentication or TLS configuration, and the default allowed origin is `*`. Use it on a local machine or trusted network; configure access controls before making it publicly accessible.

## Data preparation

The full pipeline uses **`data/peak_lists/master_filtered.txt`**, containing 4,896 reference peaks across the three reagents. It splits reference-peak records into training and test sets before generating candidate rows. The current 90%/10% split gives 4,406 training peaks and 490 test peaks.

Candidate compositions matching the assigned elemental composition receive the positive label; other generated candidates receive the negative label. Consequently, reference-peak counts and reagent percentages differ from counts and percentages over all expanded candidate rows.

Each model input contains four chemical-constraint flags, 20 encoded composition entries, and a signed ppm error. Composition counts use the scale factor defined in `include/Chem.h`. Candidate limits, reagent masks, masses, and constraint logic are defined in `src/Helpers/Chem.cpp`; preprocessing and the ±35 ppm tolerance are implemented in `src/Helpers/Preprocess.cpp`.

Theoretical singly charged positive-ion masses include the electron-mass correction. Signed mass error is calculated as `(observed m/z - theoretical ion m/z) / theoretical ion m/z × 10^6`. Supply measured m/z coordinates directly; do not apply an additional electron correction to the observed spectrum.

### Peak list format

Use the supplied cleaned peak list as the template. Files are **tab-separated**, with a header and 29 columns. Preserve column order and empty fields. The parser uses positional columns, including:

| Column (1-based) | Header | Contents |
| --- | --- | --- |
| 3 | `ion` | Assigned ion formula, including the positive charge, e.g. `(NH4)CH4O+` |
| 4 | `x0` | Numeric centroid m/z |
| 29 | `Reagent` | `NH4`, `NO`, or `H3O` |

Keep the remaining Igor-export fields in their original positions. The `sumFormula` field is retained in the file format, but the assigned formula is read from `ion`. Include the `Reagent` column if it is known.

The supplied cleaned list has had equivalent duplicate records removed. For new datasets, deduplicate or group records that share elemental composition, numerical m/z, and reagent **before splitting**. The splitting function shuffles rows; it does not automatically deduplicate a new input list or enforce separation of experimental runs.

## Runtime configuration

Configurations are in `configs/`: use `soar.yaml` for reproduction and the web API, and `test_config.yaml` for the test pipeline. The YAML comments describe the settings.

To test a different configuration, change the relevant settings and choose a new `run.run_id` so its outputs are saved separately. The reproduction run defaults to `manuscript_results`. Each pipeline saves a copy of its configuration with the results.

## Reproducibility

Keep the release source, input data, configuration, seeds, dependency versions, and build settings with the resulting artifacts. Use the recorded release/tag or archived source when reproducing the manuscript.

The repeated-training test compares exported score files within the same environment. Exact numerical agreement across different compilers, hardware, Python packages, or CPU/GPU execution is not guaranteed. Docker helps reproduce the environment, but the supplied Dockerfiles fetch the current yaml-cpp source and operating-system packages; rebuilding later is not a fully frozen dependency snapshot.

## AI-generated code and documentation

AI-generated code in this repository is limited to these two peak-fitting scripts:

- `analysis/plot_nh4_peak_fit.py` — the NH₄⁺ peak fitting visual and its supporting outputs.
- `web_app/utils/fit_peaks.py` — the standalone peak-fitting utility packaged for the web API.

**This README was also AI-generated.** It is documentation; the core candidate-generation, model-training, and inference code was not AI-generated.

## License

Distributed under the [MIT License](LICENSE).