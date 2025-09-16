#include "Chem.h"
#include "Preprocess.h"
#include <CNum.h>
#include <iostream>
#include <fstream>

using namespace CNum::Data;
using namespace CNum::Model;
using namespace CNum::Model::Tree;
using namespace CNum::Utils::ModelUtils;
using namespace CNum::DataStructs;

std::pair<int, int> check_bias(Matrix<double> &row_matrix);
void write_dataset_to_file(Matrix<double> &X, Matrix<double> &y, std::string output_path);

int main(int argc, char *argv[]) {
  if (argc != 5) {
    throw ::std::invalid_argument("Paths for the combo files, model output, and prediction output must be provided as arguments i.e. ./src/train /path/to/train/combos.csv /path/to/test/combos.csv /path/to/model/output.cmod /path/to/pred/output.txt");
  }
  
  auto train = CNum::Data::get_data(argv[1]);
  auto test = CNum::Data::get_data(argv[2]);
  
  GBModel<XGTreeBooster> xgboost("BCE",
				 200,
				 0.3,
				 1.0,
				 5,
				 3,
				 HIST,
				 "sigmoid",
				 0.0,
				 1.0,
				 0.0);


  xgboost.fit(train[0], train[1]);
  xgboost.save_model(argv[3]);
  
  auto preds = xgboost.predict(test[0]); 
  
  std::ofstream os(argv[4]);

  if (!os.is_open()) {
    std::cerr << "Error outputting to pred file" << std::endl;
    exit(1);
  }

  for (int i = 0; i < preds.get_rows(); i++) {
    double observed = test[1].get(i, 0);
    double predicted = (int) (preds.get(i, 0) >= 0.5);

    os << observed << " " << predicted << std::endl;
  }
  
  return 0;
}
