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
  constexpr uint8_t PTR_END_IDX = 4;
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

  struct ReagantIonMask {
    ::std::unique_ptr<uint8_t[]> indeces;
    uint8_t length;
  };

  class ChemMap {
  private:
    ::std::vector<double> _masses;
    ::std::vector<std::string> _chemicals;
    ::std::vector<unenc_compound> _reagant_ions;
    ::std::array<uint8_t, TOTAL_CHEMS> _proper_ordering;

    ::std::map<std::string, double> _chem_masses;
    ::std::map<std::string, double> _chem_idx;

    ::std::map<::std::string, ReagantIonMask> _ion_masks;

    ChemMap();

  public:
    static ChemMap *get_chem_map();

    ChemMap(const ChemMap &other) = delete;
    ChemMap &operator=(const ChemMap &other) = delete;
    ChemMap(ChemMap &&other) = delete;
    ChemMap &operator=(ChemMap &&other) = delete;

    double get_mass(std::string chemical);
    size_t get_idx(std::string chemical);
    const ReagantIonMask &get_reagant_ion_mask(unenc_compound reagant_ion);

    const ::std::vector<unenc_compound> &get_reagant_ions();
    const ::std::vector<double> &get_masses();
    const ::std::vector<std::string> &get_chems();
    const ::std::array<uint8_t, TOTAL_CHEMS> &get_proper_ordering();

    unenc_compound find_reagant_ion(::std::span<double> &encoded_compound_view);
  };
  
  bool compounds_are_equal(const ::std::span<double> &compound1,
			   const ::std::span<double> &compound2);
  bool uses_reagant_ion(const ::std::span<double> &encoded_unsimplified_view,
			unenc_compound reagant_ion);
  double get_ppm(double observed_mz, double theoretical_mz);
  double get_compound_mass(const ::std::span<double> &compound);
  CritereaCheckRes check_criterea(const ::std::span<double> &compound);
  CritereaCheckRes check_criterea(const ::CNum::DataStructs::Matrix<double> &compound);
  ::CNum::DataStructs::Matrix<double> factor_polyatomics(const ::CNum::DataStructs::Matrix<double> &encoded_simplified);
};

#endif
