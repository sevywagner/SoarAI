#ifndef __YAML_HELPERS_H
#define __YAML_HELPERS_H

#include <yaml-cpp/yaml.h>
#include <string>

namespace YamlHelpers { 
  ::std::string get_run_dir(const ::YAML::Node &config);
  ::std::string get_data_dir(const ::YAML::Node &config);
}

#endif
