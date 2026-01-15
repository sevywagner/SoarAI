#include "YamlHelpers.h"

namespace YamlHelpers {
  ::std::string get_run_dir(const ::YAML::Node &config) {
    return config["paths"]["core"]["run_root"].as<::std::string>() + config["run"]["run_id"].as<::std::string>() + "/";
  }
  
  ::std::string get_data_dir(const ::YAML::Node &config) {
    return get_run_dir(config) + config["paths"]["core"]["data_dir"].as<::std::string>();
  }
}
