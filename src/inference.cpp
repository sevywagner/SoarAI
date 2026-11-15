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
  if (argc != 2) {
    throw ::std::invalid_argument("Invalid arguments. Usage: inference <path to yaml config>");
  }

  auto config = ::YAML::LoadFile(argv[1]);
  auto run_dir = config["paths"]["core"]["run_root"].as<::std::string>() + config["run"]["run_id"].as<::std::string>() + "/";
  auto model_output_dir = run_dir + config["paths"]["core"]["model_output_dir"].as<::std::string>();
  
  double mz{ 0.0 };
  ::std::string mz_string;
  auto xgboost = GBModel<XGTreeBooster>::load_model(model_output_dir + "xgboost.cmod");

  ::Chem::unenc_compound reagent("");

  do {
    ::std::cout << "Enter a centroid m/z value: ";
    if (!getline(::std::cin, mz_string))
      exit(1);

    try {
      mz = ::std::stod(mz_string);
    } catch(::std::invalid_argument e) {
      ::std::cout << "Invalid centroid value." << ::std::endl;
      continue;
    }

    break;
  } while (true);
  
  do {
    ::std::cout << "Enter a reagent ion (NH4, NO, H3O, or leave blank): ";
    if (!getline(::std::cin, reagent.val))
      exit(1);
    
    if (::Chem::ChemMap::get_chem_map()->reagent_supported(reagent))
      break;

    ::std::cout << "Invalid reagent ion. Input must be NH4, NO, H3O, or leave blank." << ::std::endl;
  } while (true);
      
  auto data = Preprocess::mz_to_data(mz, reagent);
  
  data.encoded_compounds = Preprocess::simplify_compounds(data.encoded_compounds);
  auto preds = xgboost.predict(data.model_data);
  Postprocess::sort_preds(data.encoded_compounds, preds, data.ppms);

  auto dc = Preprocess::decode_compounds(data.encoded_compounds);
  
  for (int i = 0; i < dc.size(); i++)
    std::cout << dc[i].val << " " << preds.get(i, 0) << " " << data.ppms.get(i, 0) << std::endl;
  
  return 0;
}
