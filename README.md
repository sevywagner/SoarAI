# SoarAI: Ion composition identification in time-of-flight mass spectra

> Status: pre-alpha. Currently only tested on Linux (Gentoo). Data availability still pending

⚠️ License: Internal use only. Modification allowed; redistribution prohibited. See LICENSE.

This project is a part of the research I have been conducting with the Texas A&M meteorology department.
 This code demonstrates an end to end machine learning pipeline using my <a href="https://github.com/sevywagner/CNum">CNum library</a>. In this project I created a binary classifier for the identification of ion compositions in time-of-flight mass spectra. Given an m/z value my program will return a list of all the possible elemental combinations with a ppm within the range of positive or negative 50 along with the probability that they represent that peak.

## Build
```bash
mkdir build
cd build
cmake ..
make
```

## Data preparation
Sample peak lists will be provided in future releases.

### Peak list format
The peak lists used are currently exported from Igor Pro as tab seperate values in a txt file. The structure of the samples is as follows (as tagged by Igor Pro):<br></br>
def	fit	ion	x0	tag	sumFormula	x_Lo	x_Hi	d2_Ctr	d2_Lo	d2_Hi	d1_Ctr	d1_Lo	d1_Hi	d0_Ctr	d0_Lo	d0_Hi	calFac	calUnit	charge	ionizFrac	fragOf	isotopeOf

### Data prep from peak list
```bash
# Still inside build
./src/data_prep /path/to/peak/list.txt /path/to/combined/combo_file/output.csv /path/to/split/combo/output
```

The extension of the split combo output path will be added automatically i.e. if you use /home/soar the
generated files will be /home/soar_test_combos.csv and /home/soar_train_combos.csv

## Train model
```bash
# Still inside build
./src/train /path/to/train/combos.csv /path/to/test/combos.csv /path/to/model/output.cmod /path/to/pred/output.txt
```

A file containing the labeled values with the predicted values will be output by the train program.

## Analysis
```bash
# From source directory
cd analysis
python3 -m venv ./
source ./bin/activate
python3 -m pip install -r requirements.txt
python3 analysis.py /path/to/output/plot/directory/ /path/to/preds1.txt /path/to/preds2.txt ...
```

You can add as many prediction files as you'd like. This program will display binary classification metr
ics such as F1 and ROC AUC and generate a ROC plot and a confusion matrix for each pred file and output
them all to the output plot directory.

## Inference
```bash
# In build directory
./src/infer /path/to/model.cmod
```

Upon executing this program it will expect keyboard input. Input a m/z value and it will output its predictions. 
