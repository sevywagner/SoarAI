#ifndef __PREPROCESS_H
#define __PREPROCESS_H

#include <CNum.h>
#include <future>
#include <string>
#include <sstream>
#include <fstream>
#include <mutex>
#include <malloc.h>

#include "Chem.h"

namespace Preprocess {
  struct PeakListData {
    ::std::vector<double> mz;
    ::std::vector< Chem::unenc_compound > compound_strings;
  };
  
  struct CompoundPermutations {
    ::CNum::DataStructs::Matrix<double> compounds;
    std::vector<double> masses;
  };
  
  struct MSData {
    ::CNum::DataStructs::Matrix<double> model_data;
    ::CNum::DataStructs::Matrix<double> encoded_compounds;
  };
  
  ::CNum::DataStructs::Matrix<double> encode_compounds(const std::vector< Chem::unenc_compound > &compound_strings);
  CompoundPermutations all_possible_elemental_combo(double mass, enum Chem::mode m);
  std::vector< Chem::unenc_compound > decode_compounds(const ::CNum::DataStructs::Matrix<double> &encoded_compounds);
  ::CNum::DataStructs::Matrix<double> simplify_compounds(const ::CNum::DataStructs::Matrix<double> &unsimplified_compounds);
  MSData mz_to_data(double mz, enum Chem::mode m, size_t n_features = 18);
  std::pair<int, int> check_bias(::CNum::DataStructs::Matrix<double> &row_matrix);

  namespace PrepareDataset {
    PeakListData parse_peak_list(::std::string path);
    void create_combo_file(::std::string output_path,
			   const PeakListData &peak_list_data,
			   const ::CNum::DataStructs::Matrix<double> &encoded_unsimplified,
			   const ::CNum::DataStructs::Matrix<double> &encoded_simplified);
    void negative_sample_reduction(::std::string path);
    void train_test_split(::std::string combo_file_path,
			  ::std::string output_path,
			  size_t n_test_pos = 548,
			  size_t n_test_neg = 10000);
  };
};

#endif
