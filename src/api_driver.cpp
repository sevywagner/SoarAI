#include "InferenceAPI.h"
#include <CNum/Deploy.h>

#include <stdexcept>

using namespace ::CNum::Model::Tree;
using namespace ::InferenceAPI;

int main(int argc, char **argv) {
  if (argc != 3)
    throw ::std::invalid_argument("Incorrect arguments. Usage: ./src/api_driver <path to yaml config> <reagent ion <NH4|NO>>");

  auto config = ::YAML::LoadFile(argv[1]);
  ::std::string reagent_ion = argv[2];
  
  resolve_paths(config);
  
  CNum::Deploy::InferenceAPI< GBModel<XGTreeBooster>,
			      Storage > rest_api(config["paths"]["api"][reagent_ion + "_model_path"].as<::std::string>(),
						 preprocess_func,
						 postprocess_func,
						 config["api"]["allowed_origins"].as<::std::string>(),
						 config["api"]["n_model_instances"].as<int>(),
						 config["api"][reagent_ion + "_port"].as<unsigned short>());

  constexpr char url[::CNum::Deploy::MAX_URL_LEN] = "/process-graph"; // C-style string necessary here because ::std::string can't be constexpr until C++23
  constexpr ::CNum::Deploy::PathString url_path(url);
  rest_api.add_inference_route< url_path >(crow::HTTPMethod::Post, process_graph);
  rest_api.start();
  
  return 0;
}
