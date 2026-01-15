#ifndef __PREPROCESS_H
#define __PREPROCESS_H

#include <CNum.h>
#include <future>
#include <string>
#include <sstream>
#include <fstream>
#include <mutex>
#include <malloc.h>
#include <span>
#include <unordered_set>

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
    ::CNum::DataStructs::Matrix<double> ppms;
  };

  struct Bias {
    size_t ones;
    size_t zeros;
  };

  ::CNum::Model::Tree::SubsampleFunction get_subsample_func(::std::vector<size_t> &ones_indeces,
							    ::std::unordered_set<size_t> &ones_indeces_set);
  
  ::CNum::DataStructs::Matrix<double> encode_compounds(const std::vector< Chem::unenc_compound > &compound_strings);
  CompoundPermutations all_possible_elemental_combo(double mass, Chem::unenc_compound reagant_ion);
  std::vector< Chem::unenc_compound > decode_compounds(const ::CNum::DataStructs::Matrix<double> &encoded_compounds);
  ::CNum::DataStructs::Matrix<double> simplify_compounds(const ::CNum::DataStructs::Matrix<double> &unsimplified_compounds);
  MSData mz_to_data(double mz, Chem::unenc_compound reagant_ion, size_t n_features = 18);
  Bias check_bias(::CNum::DataStructs::Matrix<double> &row_matrix);

  namespace PrepareDataset {
    void split_by_reagant_ion(::std::string combos_path, ::std::string output_path);
    PeakListData parse_peak_list(::std::string path);
    Bias create_combo_file(::std::string output_path,
			   ::std::string unidentified_combos_path,
			   ::Chem::unenc_compound reagant_ion,
			   const PeakListData &peak_list_data,
			   const ::CNum::DataStructs::Matrix<double> &encoded_unsimplified,
			   const ::CNum::DataStructs::Matrix<double> &encoded_simplified,
			   int n_threads = 10);
    void negative_sample_reduction(::std::string path);
    void train_test_split(::std::string combo_file_path,
			  ::std::string output_path,
			  size_t n_test_pos,
			  size_t n_test_neg);
    void peak_list_train_test_split(::std::string peak_list_path,
				    ::std::string output_path,
				    double test_percentage = 0.1);
  };
};

#endif
