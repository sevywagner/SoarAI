#include <yaml-cpp/yaml.h>
#include <CNum.h>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <filesystem>
#include <stdexcept>

#include "Chem.h"
#include "Preprocess.h"
#include "YamlHelpers.h"

using namespace CNum::Data;
using namespace CNum::Model;
using namespace CNum::Model::Tree;
using namespace CNum::Utils::ModelUtils;
using namespace CNum::DataStructs;

int main(int argc, char *argv[]) {
  if (argc != 2 && argc != 3) {
    throw ::std::invalid_argument("Invalid arguments. Usage: ./src/train <path to yaml config> [ dump ]");
  }

  auto config = ::YAML::LoadFile(argv[1]);
  
  ::std::vector<size_t> ones_indeces;
  ::std::unordered_set<size_t> ones_indeces_set;

  auto subsample = ::Preprocess::get_subsample_func(ones_indeces, ones_indeces_set);

  auto combo_dir = ::YamlHelpers::get_data_dir(config) + config["paths"]["core"]["data_combo_files_dir"].as<::std::string>();
  auto train_combos_path = combo_dir + config["paths"]["core"]["data_combo_train"].as<::std::string>() + ".csv";
  auto test_combos_path = combo_dir + config["paths"]["core"]["data_combo_test"].as<::std::string>();

  auto train = CNum::Data::get_data(train_combos_path);
  
  if (config["run"]["deterministic"].as<bool>()) {
    int seed = config["run"]["cnum_model_train_seed"].as<int>();
    CNum::Utils::Rand::RandomGenerator::set_global_seed(seed);
  }
  
  GBModel<XGTreeBooster> xgboost("BCE",
				 config["core"]["model_hyperparams"]["xgboost"]["n_learners"].as<int>(),
				 config["core"]["model_hyperparams"]["xgboost"]["learning_rate"].as<double>(),
				 config["core"]["model_hyperparams"]["xgboost"]["subsample"].as<double>(),
				 5,
				 3,
				 HIST,
				 "sigmoid",
				 0.0,
				 1.0,
				 0.0,
				 subsample);


  xgboost.fit(train[0], train[1], true);
  
  auto run_dir = ::YamlHelpers::get_run_dir(config);
  
  auto model_save_dir = run_dir + config["paths"]["core"]["model_output_dir"].as<::std::string>();
  ::std::filesystem::create_directory(model_save_dir);  
  xgboost.save_model(model_save_dir + "xgboost.cmod");

  auto pred_output_dir = run_dir + config["paths"]["core"]["pred_output_dir"].as<::std::string>();
  ::std::filesystem::create_directory(pred_output_dir);

  auto *cm = ::Chem::ChemMap::get_chem_map();
  const auto &reagent_ions = cm->get_reagent_ions();
  ::std::vector<::std::string> combo_paths;
  combo_paths.reserve(reagent_ions.size() + 1);
  combo_paths.push_back("");
  for (const auto &ion: reagent_ions)
    combo_paths.push_back(ion.val);
    
  for (const auto &path: combo_paths) {
    std::ofstream os(pred_output_dir + (path != "" ? path + "_" : path) + "xgboost_preds.txt");

    if (!os.is_open()) {
      throw ::std::runtime_error("Error outputting to pred file in main function");
    }

    ::std::array< CNum::DataStructs::Matrix<double>, 2 > test;
    try {
      test = CNum::Data::get_data(test_combos_path + (path != "" ? "_" + path : path) + ".csv");
    } catch (::std::invalid_argument e) {
      continue; // file empty
    }

    auto preds = xgboost.predict(test[0]);
  
    for (size_t i{}; i < preds.size(); i++) {
      double observed = test[1][i];
      double predicted = static_cast<int>(preds[i] >= 0.5);
    
      os << observed << " " << predicted << std::endl;
    }

    if (argc == 3) {
      auto logits_output_dir = pred_output_dir + config["paths"]["core"]["pred_logit_output_dir"].as<::std::string>();
      ::std::filesystem::create_directory(logits_output_dir);
      ::std::ofstream logits_of(logits_output_dir + (path != "" ? path + "_" : path) + "xgboost_logits_" + argv[2] + ".txt");

      if (!logits_of.is_open())
	throw ::std::runtime_error("Logit output error in train -- Could not open logit path");

      for (size_t i{}; i < preds.size(); i++) {
	logits_of << preds[i] << ::std::endl;
      }
    }
  }
  
  return 0;
}
