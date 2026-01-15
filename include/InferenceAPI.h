#ifndef __SOAR_INFERENCE_API_H
#define __SOAR_INFERENCE_API_H

#include <iostream>
#include <iomanip>
#include <ctime>
#include <CNum.h>
#include <crow.h>
#include <string.h>
#include <sys/stat.h>
#include <stdexcept>
#include <cmath>
#include <yaml-cpp/yaml.h>
#include <filesystem>

#include "Preprocess.h"
#include "Postprocess.h"
#include "SysUtils.h"
#include "Chem.h"

namespace InferenceAPI {
  constexpr int N_FILES = 2; // 2 files for mz_av and mz_base
  constexpr size_t MAX_FILE_SIZE = 7 * (1 << 20); // 7 MiB
  
  const ::std::array<size_t, 3> TABLE_WIDTHS = { 15, 11, 11 };
  const ::std::array<::std::string, 3> TABLE_HEADERS = { "Ion", "PPM", "Confidence" };
  constexpr size_t TABLE_MARGIN = 1;
  constexpr size_t TABLE_N_COLS = 3;

  struct Storage {
    ::CNum::DataStructs::Matrix<double> encoded_compounds;
    ::CNum::DataStructs::Matrix<double> ppms;
    ::CNum::DataStructs::Matrix<uint8_t> criterea_encodings;
  };

  extern char *python_executable_path; // to be used to c code hence the NULL over nulltpr
  extern ::std::string graph_upload_dir;
  extern ::std::string peak_output_dir;

  void resolve_paths(const ::YAML::Node &config);
  
  ::CNum::DataStructs::Matrix<double> preprocess_func(crow::json::rvalue &req_body,
						      crow::json::wvalue &res_body,
						      Storage &storage);
  void postprocess_func(::CNum::DataStructs::Matrix<double> &preds,
			crow::json::wvalue &res_body,
			Storage &storage);
  void process_graph(const crow::request &req,
		     crow::response &res,
		     ::CNum::Model::Tree::GBModel< ::CNum::Model::Tree::XGTreeBooster > *model);
}
 
#endif
