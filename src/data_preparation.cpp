#include <CNum.h>
#include <stdexcept>
#include <yaml-cpp/yaml.h>
#include <filesystem>
#include <iostream>

#include "YamlHelpers.h"
#include "Preprocess.h"

int main(int argc, char *argv[]) {
  if (argc != 2)
    throw ::std::invalid_argument("Invalid arguments. Usage: ./src/data_prep <path to yaml config>");

  auto config = ::YAML::LoadFile(argv[1]);
  auto data_dir = ::YamlHelpers::get_data_dir(config);

  auto deterministic = config["run"]["deterministic"].as<bool>();

  if (deterministic)
    ::CNum::Utils::Rand::RandomGenerator::set_global_seed(config["run"]["data_prep_seed"].as<int>());
  
  auto peak_list_dir = data_dir + config["paths"]["core"]["data_out_peak_list_dir"].as<::std::string>();
  auto combo_files_dir = data_dir + config["paths"]["core"]["data_combo_files_dir"].as<::std::string>();
  auto unidentified_combos_dir = combo_files_dir + config["paths"]["core"]["data_combo_unidentified_dir"].as<::std::string>();
  auto in_peak_list_path = config["paths"]["core"]["in_peak_list_dir"].as<::std::string>() + config["paths"]["core"]["in_peak_list_filename"].as<::std::string>();
  ::std::filesystem::create_directories(peak_list_dir);
  ::std::filesystem::create_directory(combo_files_dir);
  ::std::filesystem::create_directory(unidentified_combos_dir);
  
  Preprocess::PrepareDataset::split_by_reagant_ion(in_peak_list_path,
						   peak_list_dir);

  auto *cm = Chem::ChemMap::get_chem_map();
  const auto &reagant_ions = cm->get_reagant_ions();

  int n_threads = deterministic ? 1 : ::std::thread::hardware_concurrency();
  ::std::array<::std::string, 2> test_train_ext({ "_test", "_train" });
  for (const auto &ion: reagant_ions) {
    auto peak_list_path = peak_list_dir + ion.val;
    auto combo_file_path = combo_files_dir + ion.val;
    auto unidentified_combos_path = unidentified_combos_dir + ion.val;

    Preprocess::PrepareDataset::peak_list_train_test_split(peak_list_path + ".txt", peak_list_path, .1);

    Preprocess::Bias bias;
    for (const auto &ext: test_train_ext) {
      auto peak_list_data = Preprocess::PrepareDataset::parse_peak_list(peak_list_path + ext + ".txt");
      auto encoded_unsimplified = Preprocess::encode_compounds(peak_list_data.compound_strings);
      auto encoded_simplified = Preprocess::simplify_compounds(encoded_unsimplified);

      bias = Preprocess::PrepareDataset::create_combo_file(combo_file_path + ext + ".csv",
							   unidentified_combos_path + ext + "_unidentified_compounds.csv",
							   ion,
							   peak_list_data,
							   encoded_unsimplified,
							   encoded_simplified,
							   n_threads);
      
      ::std::cout << ion.val << " Bias (" << ext.substr(1) << ")" << ": " << ::std::endl
		  << "Positive samples: " << bias.ones << ::std::endl
		  << "Negative samples: " << bias.zeros << ::std::endl;
    }
  }

  return 0;
}

