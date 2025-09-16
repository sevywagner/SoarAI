#include "Chem.h"
#include "Preprocess.h"
#include "Postprocess.h"
#include <CNum.h>
#include <iostream>

using namespace CNum::Data;
using namespace CNum::Model;
using namespace CNum::Model::Tree;
using namespace CNum::Utils::ModelUtils;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    throw ::std::invalid_argument("Must pass path to model as argument i.e. ./this_program /path/to/model.cmod");
  }
  
  double mz{ 0.0 };
  auto xgboost = GBModel<XGTreeBooster>::load_model(argv[1]);
  
  ::std::cin >> mz;
  
  auto data = Preprocess::mz_to_data(mz, Chem::PTR);
  data.encoded_compounds = Preprocess::simplify_compounds(data.encoded_compounds);
  auto preds = xgboost.predict(data.model_data);
  Postprocess::sort_preds(data.encoded_compounds, preds);

  auto dc = Preprocess::decode_compounds(data.encoded_compounds);

  for (int i = 0; i < dc.size(); i++)
    std::cout << dc[i].val << " " << preds.get(i, 0) << std::endl;
  
  return 0;
}
