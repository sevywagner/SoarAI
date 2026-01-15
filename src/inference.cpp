#include "Chem.h"
#include "Preprocess.h"
#include "Postprocess.h"
#include <CNum.h>
#include <iostream>
#include <yaml-cpp/yaml.h>

using namespace CNum::Data;
using namespace CNum::Model;
using namespace CNum::Model::Tree;
using namespace CNum::Utils::ModelUtils;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    throw ::std::invalid_argument("Invalid arguments. Usage: inference <path to yaml config> [NH4|NO]");
  }

  auto config = ::YAML::LoadFile(argv[1]);
  auto run_dir = config["paths"]["core"]["run_root"].as<::std::string>() + config["run"]["run_id"].as<::std::string>() + "/";
  auto model_output_dir = run_dir + config["paths"]["core"]["model_output_dir"].as<::std::string>();
  
  double mz{ 0.0 };
  auto xgboost = GBModel<XGTreeBooster>::load_model(model_output_dir + ::std::string(argv[2]) + "_xgboost.cmod");
  
  ::std::cin >> mz;
  
  auto data = Preprocess::mz_to_data(mz, { argv[2] });
  data.model_data.print_matrix();
  
  
  data.encoded_compounds = Preprocess::simplify_compounds(data.encoded_compounds);
  auto preds = xgboost.predict(data.model_data);
  Postprocess::sort_preds(data.encoded_compounds, preds, data.ppms);

  auto dc = Preprocess::decode_compounds(data.encoded_compounds);
  
  for (int i = 0; i < dc.size(); i++)
    std::cout << dc[i].val << " " << preds.get(i, 0) << " " << data.ppms.get(i, 0) << std::endl;
  
  return 0;
}
