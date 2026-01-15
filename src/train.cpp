#include <yaml-cpp/yaml.h>
#include <CNum.h>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <filesystem>

#include "Chem.h"
#include "Preprocess.h"
#include "YamlHelpers.h"

using namespace CNum::Data;
using namespace CNum::Model;
using namespace CNum::Model::Tree;
using namespace CNum::Utils::ModelUtils;
using namespace CNum::DataStructs;

int main(int argc, char *argv[]) {
  if (argc != 3 && argc != 4) {
    throw ::std::invalid_argument("Invalid arguments. Usage: ./src/train <path to yaml config> <reagent ion (NH4|NO)> [ dump ]");
  }

  auto config = ::YAML::LoadFile(argv[1]);
  ::std::string reagent_ion = argv[2];
  
  ::std::vector<size_t> ones_indeces;
  ::std::unordered_set<size_t> ones_indeces_set;

  auto subsample = ::Preprocess::get_subsample_func(ones_indeces, ones_indeces_set);

  auto combo_dir = ::YamlHelpers::get_data_dir(config) + config["paths"]["core"]["data_combo_files_dir"].as<::std::string>();
  auto train_combos_path = combo_dir + config["paths"]["core"]["data_combo_" + reagent_ion + "_train"].as<::std::string>();
  auto test_combos_path = combo_dir + config["paths"]["core"]["data_combo_" + reagent_ion + "_test"].as<::std::string>();

  auto train = CNum::Data::get_data(train_combos_path);
  auto test = CNum::Data::get_data(test_combos_path);

  if (config["run"]["deterministic"].as<bool>()) {
    int seed = config["run"]["cnum_model_train_seed"].as<int>();
    CNum::Utils::Rand::RandomGenerator::set_global_seed(seed);
  }
    
  GBModel<XGTreeBooster> xgboost("BCE",
				 config["core"]["model_hyperparams"]["xgboost"][reagent_ion + "_reagent"]["n_learners"].as<int>(),
				 config["core"]["model_hyperparams"]["xgboost"][reagent_ion + "_reagent"]["learning_rate"].as<double>(),
				 config["core"]["model_hyperparams"]["xgboost"][reagent_ion + "_reagent"]["subsample"].as<double>(),
				 5,
				 3,
				 HIST,
				 "sigmoid",
				 0.0,
				 1.0,
				 0.0,
				 subsample);


  xgboost.fit(train[0], train[1], false);
  auto preds = xgboost.predict(test[0]);

  auto run_dir = ::YamlHelpers::get_run_dir(config);
  
  auto model_save_dir = run_dir + config["paths"]["core"]["model_output_dir"].as<::std::string>();
  ::std::filesystem::create_directory(model_save_dir);
  xgboost.save_model(model_save_dir + reagent_ion + "_xgboost.cmod");

  auto pred_output_dir = run_dir + config["paths"]["core"]["pred_output_dir"].as<::std::string>();
  ::std::filesystem::create_directory(pred_output_dir);
  std::ofstream os(pred_output_dir + reagent_ion + "_xgboost_preds.txt");

  if (!os.is_open()) {
    throw ::std::runtime_error("Error outputting to pred file in main function");
  }
  
  for (size_t i{}; i < preds.size(); i++) {
    double observed = test[1][i];
    double predicted = static_cast<int>(preds[i] >= 0.5);
    
    os << observed << " " << predicted << std::endl;
  }

  if (argc == 4) {
    auto logits_output_dir = pred_output_dir + config["paths"]["core"]["pred_logit_output_dir"].as<::std::string>();
    ::std::filesystem::create_directory(logits_output_dir);
    ::std::ofstream logits_of(logits_output_dir + reagent_ion + "_xgboost_logits_" + ::std::string(argv[3]) + ".txt");

    if (!logits_of.is_open())
      throw ::std::runtime_error("Logit output error in train -- Could not open logit path");

    for (size_t i{}; i < preds.size(); i++) {
      logits_of << preds[i] << ::std::endl;
    }
  }
  
  return 0;
}
