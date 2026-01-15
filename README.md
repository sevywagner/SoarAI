# SoarAI: Ion composition identification in time-of-flight mass spectra

> Status: v1.0.0. Currently only supported on Linux. 

This project is a part of the research I have been conducting with the Texas A&M meteorology department. This code demonstrates an end to end machine learning pipeline using my <a href="https://github.com/sevywagner/CNum">CNum library</a>. In this project I created a binary classifier for the identification of ion compositions in time-of-flight mass spectra. Given an m/z value my program will return a list of all the possible elemental combinations with a ppm within the range [-50, 50] along with the probability that they represent the peak at the m/z value.

## Dependencies
- <a href="https://github.com/sevywagner/CNum/releases/tag/v0.2.2">CNum v0.2.2</a> - Required for CNum GBModel training and using the REST API
- <a href="https://github.com/mikefarah/yq">yq</a> - Required for reading yaml configs in shell scripts

## Supported platforms
- Linux (Gentoo, Arch, Ubuntu) (GCC)

## Entry Points
- tests/ - Test the entire pipeline and determinism on a small dataset
- reproduce/ - Reproduce all of the results models and figures from the manuscript
- web_app/ - Locally host SoarAI

All provided bash scripts build all necessary binaries and prepare all necessary Python environments for their given tasks. 

## Tests
Testing the pipeline is recommended to ensure that the training is deterministic and all of the modules build and execute properly on your machine.

To run the tests navigate to the "tests" directory and run the provided shell script:
```bash
# From SoarAI source directory
cd tests
./shell_scripts/run_test_pipeline.sh
```

The results of the determinism tests are printed to the terminal. The artifacts produced by the test are saved to tests/test_results along with a copy of the test runtime configuration. 

## Reproduce manuscript results
To reproduce all models and figures in the manuscript navigate to the "reproduce" directory and run the provided shell script:
```bash
# From SoarAI source directory
cd reproduce
./shell_scripts/reproduce_all.sh
```

To reproduce the manuscript results use the provided yaml configuration, or at least preserve the seeds, hyperparameters, and deterministic setting from the default configuration if you must change the layout of the directory. 

All artifacts produced will be output to the run directory which is ./run/reproduced_results by default. Inside there will be "models", "preds", "plots", and "data" directories. 8 total models will be produced there is one CNum XGBoost, Sklearn GBDT, Flaml LightGBM, and TensorFlow neural network for both NH4+ and NO+ reagent ions. The models directory contains the trained model files. The "preds" directory contains the predictions of all of the models on the test dataset. The plots directory contains all of the reproduced figures from the manuscript, and the data directory contains the split peak lists, and combo files. Combo files are the files in which the preprocessed data passed directly to the model is stored. 

## Web App
You can also host SoarAI locally by navigating to the "web_app" directory and running the provided shell script
```bash
# From SoarAI source directory
cd web_app
./shell_scripts/start_api.sh
```

This will start two REST APIs, one for the NH4+ xgboost model and one for the NO+ xgboost model on ports 18080 and 18081, respectively.
### API endpoints
- predict/ - provided by the CNum InferenceAPI interface
- process-graph/ - takes in a mass spectrum, fits peaks (naively), and assigns formulas

### *Important*
The CNum inference API tools use the Crow C++ microframework which has shown vulnerabilites in the past. If you plan on hosting this and don't plan on it being only for your local network, an extra layer of security is highly recommended, for example token-based authorization and tunneling (i.e. via Cloudflare). The API also uses an exec function to execute a binary which is handled safely, but always has its inherent risks, so for this version only using the API locally is strongly recommended. 

## Using the pipeline piece by piece
If you prefer to use the pipeline piece by piece, here is a description of all the shell scripts:
- data_prep.sh - Prepare data for both NO and NH4 models
- validate_config.sh - Check to make sure the dirs in the yaml configs are properly formatted
- build_cpp_bin.sh - Build all C++ modules
- prepare_python_venvs.sh - Create Python venvs and install dependancies
- train.sh - Train NH4+ and NO+ CNum GBModels
- create_py_models.sh - Train all python models (LGBM, GBDT and Neural Net) (NH4+ and NO+)
- infer.sh - Make inference on a single peak (input m/z value at peak)
- analysis.sh - Make all of the figures for the model results

## Data preparation
An example of data from real MS expiriments has been provided in data/peak_lists/master.txt. This file can be used to run the entirety of pipeline. 

### Peak list format
The peak lists used are currently exported from Igor Pro as tab seperate values in a txt file. The structure of the samples is as follows (as tagged by Igor Pro):<br></br>
def	fit	ion	x0	tag	sumFormula	x_Lo	x_Hi	d2_Ctr	d2_Lo	d2_Hi	d1_Ctr	d1_Lo	d1_Hi	d0_Ctr	d0_Lo	d0_Hi	calFac	calUnit	charge	ionizFrac	fragOf	isotopeOf

### Yaml runtime configurations
The yaml runtime configurations in the "configs" directory can be used for setting up paths and expirement tracking. To run new expiriments simply change the run_id and the artifacts along with a copy of the config will be saved to a new folder. You can change the seeds and the hyperparameters of the models to get different results. 

### License
This project is distributed under the MIT license
