#ifndef __CHEM_H
#define __CHEM_H

#include <CNum.h>
#include <map>
#include <string>
#include <algorithm>
#include <vector>

namespace Chem {
  constexpr uint8_t N_CRITEREA = 4;
  constexpr uint8_t TOTAL_CHEMS = 13;
  constexpr uint8_t PTR_END_IDX = 4; // ptr as in the mass spectrometry mode not pointer
  constexpr uint8_t POLYATOMIC_START_IDX = 9;
  constexpr uint8_t CHARGE = 8;
  constexpr double CHEM_SCALE_FACTOR = .01;

  struct unenc_compound {  
    std::string val;
    unenc_compound(std::string str) : val(str) {}
  };

  struct CritereaCheckRes {
    ::std::vector<bool> crit_mask;
    bool did_pass;
  };

  enum mode {
    PTR,
    NON_PTR
  };

  class ChemMap {
  private:
    static ChemMap *_chem_map;
    std::vector<double> _masses;
    std::vector<std::string> _chemicals;
    std::map<std::string, double> _chem_masses;
    std::map<std::string, double> _chem_idx;
    ::CNum::DataStructs::Mask<::CNum::DataStructs::IDX, uint32_t> _PTR_mask;

    ChemMap();

  public:
    static ChemMap *get_chem_map();

    ChemMap(const ChemMap &other) = delete;
    ChemMap &operator=(const ChemMap &other) = delete;
    ChemMap(ChemMap &&other) = delete;
    ChemMap &operator=(ChemMap &&other) = delete;

    double get_mass(std::string chemical);
    size_t get_idx(std::string chemical);

    std::vector<double> get_masses();
    std::vector<std::string> get_chems();
  };
  
  bool compounds_are_equal(const ::CNum::DataStructs::Matrix<double> &compound1, const ::CNum::DataStructs::Matrix<double> &compound2);
  bool is_PTR(const ::CNum::DataStructs::Matrix<double> &encoded_unsimplified);
  double get_ppm(double observed_mz, double theoretical_mz);
  double get_compound_mass(const ::CNum::DataStructs::Matrix<double> &compound);
  CritereaCheckRes check_criterea(const ::CNum::DataStructs::Matrix<double> &compound);
  ::CNum::DataStructs::Matrix<double> factor_polyatomics(const ::CNum::DataStructs::Matrix<double> &encoded_simplified);
};

#endif
