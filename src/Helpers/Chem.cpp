#include "Chem.h"
#include "Preprocess.h"

using namespace CNum::DataStructs;

namespace Chem {
  // ----------------------
  // Chem Map Singleton
  // ----------------------

  // ---- Constructor ----
  ChemMap::ChemMap() {
    _chemicals.emplace_back("H");
    _chemicals.emplace_back("C");
    _chemicals.emplace_back("N");
    _chemicals.emplace_back("O");
    _chemicals.emplace_back("F");
    _chemicals.emplace_back("Si");
    _chemicals.emplace_back("S");
    _chemicals.emplace_back("Cl");
    _chemicals.emplace_back("+");
    _chemicals.emplace_back("NH4");
    _chemicals.emplace_back("NH3");
    _chemicals.emplace_back("H3O");
    _chemicals.emplace_back("H2O");

    _masses.push_back(1.00782503);
    _masses.push_back(12);
    _masses.push_back(14.003074);
    _masses.push_back(15.99491462);
    _masses.push_back(18.998403163);
    _masses.push_back(27.97692653465);
    _masses.push_back(31.97207117);
    _masses.push_back(34.968852682);
    _masses.push_back(-0.000548);
    _masses.push_back(18.03437412);
    _masses.push_back(17.02654909);
    _masses.push_back(19.01838971);
    _masses.push_back(18.01056468);

    for (int i{}; i < TOTAL_CHEMS; i++) {
      _chem_masses[_chemicals[i]] = _masses[i];
      _chem_idx[_chemicals[i]] = i;
    }

    constexpr uint8_t n_chems_detected_NH4_reagant = 4;
    constexpr uint8_t n_chems_detected_NO_reagant = 8;

    _reagant_ions = { {"NH4"}, {"NO"} };

    auto mask_ptr = std::make_unique<uint8_t[]>(n_chems_detected_NH4_reagant);
    ::std::iota(mask_ptr.get(), mask_ptr.get() + n_chems_detected_NH4_reagant, 0);
    _ion_masks["NH4"] = { ::std::move(mask_ptr), n_chems_detected_NH4_reagant };

    mask_ptr.reset();
    mask_ptr = ::std::make_unique<uint8_t[]>(n_chems_detected_NO_reagant);
    ::std::iota(mask_ptr.get(), mask_ptr.get() + n_chems_detected_NO_reagant, 0);
    _ion_masks["NO"] = { ::std::move(mask_ptr), n_chems_detected_NO_reagant };

    /* _chemicals is sorted by mass in ascending order for the elemental combo recursive function
       _proper_ordering can be used in Preprocess::decode_compounds to maintain standard chemical order
       in the decoded compound */
    _proper_ordering = { 9, 10, 11, 12, 1, 7, 4, 0, 2, 3, 6, 5, 8 };
  }

  // ---- Get singleton instance ----
  ChemMap *ChemMap::get_chem_map() {
    static auto *chem_map = new ChemMap(); // Leak by design to avoid shutdown order issues
    return chem_map;
  }

  // ---- Get mass of a chemical ----
  double ChemMap::get_mass(std::string chemical) {
    return _chem_masses[chemical];
  }

  // ---- Get the index of a chemical in the predefined order ----
  size_t ChemMap::get_idx(std::string chemical) {
    return _chem_idx[chemical];
  }

  const ReagantIonMask &ChemMap::get_reagant_ion_mask(unenc_compound reagant_ion) {
    return _ion_masks[reagant_ion.val];
  }
  
  // ---- Get masses vec ----
  const std::vector<double> &ChemMap::get_masses() {
    return _masses;
  }

  // ---- Get chem vec ----
  const std::vector<std::string> &ChemMap::get_chems() {
    return _chemicals;
  }

  const std::vector<unenc_compound> &ChemMap::get_reagant_ions() {
    return _reagant_ions;
  }

  const ::std::array<uint8_t, TOTAL_CHEMS> &ChemMap::get_proper_ordering() {
    return _proper_ordering;
  }

  unenc_compound ChemMap::find_reagant_ion(::std::span<double> &encoded_compound_view) {
    unenc_compound res{ "" };
    uint8_t ri_ctr{ 0 };
    
    while (ri_ctr < _reagant_ions.size() &&
	   !::Chem::uses_reagant_ion(encoded_compound_view, _reagant_ions[ri_ctr])) { ri_ctr++; }
    
    if (ri_ctr < _reagant_ions.size())
      res = _reagant_ions[ri_ctr];
    
    return res;
  }

  // ----------------
  // Chem utilities
  // ----------------

  // ---- Compare compounds ----
  bool compounds_are_equal(const ::std::span<double> &compound1, const ::std::span<double> &compound2) {
    if (compound1.size() != compound2.size()) {
      throw ::std::invalid_argument("Compound comparison err -- misaligned dims");
    }

    for (size_t i{}; i < TOTAL_CHEMS; i++) {
      int compound1_el_ct = ::std::round(compound1[i] / CHEM_SCALE_FACTOR);
      int compound2_el_ct = ::std::round(compound2[i] / CHEM_SCALE_FACTOR);
      if (compound1_el_ct != compound2_el_ct) {
	return false;
      }
    } 

    return true;
  }

  // ---- Get the ppm of a compound ----
  double get_ppm(double observed_mz, double theoretical_mz) {
    constexpr int ppm_normalization_factor = 1e6;
    return ((observed_mz - theoretical_mz) / theoretical_mz) * ppm_normalization_factor;
  }

  // ---- Get the mass of a compound ----
  double get_compound_mass(const ::std::span<double> &compound) {
    if (compound.size() != TOTAL_CHEMS) {
      throw ::std::invalid_argument("Malformed compound encoding in get_compound_mass -- size != TOTAL_CHEMS");
    }

    auto *cm = ChemMap::get_chem_map();
    auto masses = cm->get_masses();
  
    double total_mass{};
    for (size_t i{}; i < TOTAL_CHEMS; i++) {
      total_mass += (masses[i] * ::std::round(compound[i] / CHEM_SCALE_FACTOR));
    }

    return total_mass;
  }

  template <typename T> static CritereaCheckRes crit_check_logic(const T &compound) {
    std::vector<bool> criterea_mask;
    for (int i{}; i < 4; i++)
      criterea_mask.push_back(true);
  
    bool did_pass{true};
    auto *cm = ChemMap::get_chem_map();

    bool contains_nh4 = compound[cm->get_idx("NH4")] > 0;
    bool contains_nh3 = compound[cm->get_idx("NH3")] > 0;
    uint32_t n_hydrogen = ::std::round(compound[cm->get_idx("H")] / CHEM_SCALE_FACTOR);
    uint32_t n_carbon = ::std::round(compound[cm->get_idx("C")] / CHEM_SCALE_FACTOR);
    uint32_t n_nitrogen = ::std::round(compound[cm->get_idx("N")] / CHEM_SCALE_FACTOR);

    if (contains_nh4) {
      // has NH4 and the number of hydrogen is odd
      if (n_hydrogen % 2 == 1) { 
	criterea_mask[0] = false;
	did_pass = false;
      }

      // has NH4 and the number of hydrogen is greater than 2 + number of carbon
      if (n_hydrogen > (2 * n_carbon) + 2) { 
	criterea_mask[1] = false;
	did_pass = false;
      }
    }

    
    if (n_nitrogen == 0 && !contains_nh4 && !contains_nh3) {
      // no nitrogen and number of hydrogen is even
      if (n_hydrogen % 2 == 0) {
	criterea_mask[2] = false;
	did_pass = false;
      }

      // no nitrogen and number of hydrogen is greater than 2 times number of carbon plus 3
      if (n_hydrogen > (2 * n_carbon) + 3) {
	criterea_mask[3] = false;
	did_pass = false;
      }
    }

    return { std::move(criterea_mask), did_pass };
  }

  // ---- Check the 4 criterea of a compound ----
  CritereaCheckRes check_criterea(const ::std::span<double> &compound) {
    if (compound.size() != TOTAL_CHEMS) {
      throw ::std::invalid_argument("Malformed compound encoding in check_criterea -- size != TOTAL_CHEMS");
    }
    
    return crit_check_logic(compound);
  }

  CritereaCheckRes check_criterea(const ::CNum::DataStructs::Matrix<double> &compound) {
    if (compound.get_cols() > 1)
      throw ::std::invalid_argument("Check criterea error -- whole matrix provided, use matrix.get(ROW, i) to pass a single compound");
    
    return crit_check_logic(compound);
  }

  // ---- "Unsimplify" compoound by factoring polyatomics out and adjusting compound accordingly ----
  Matrix<double> factor_polyatomics(const Matrix<double> &encoded_simplified) {
    auto encoded_unsimplified = (Matrix<double>(std::cref(encoded_simplified))).move_ptr();

    constexpr int ammonium_idx = TOTAL_CHEMS - 1;
    auto *cm = ChemMap::get_chem_map();
    auto chems = cm->get_chems();
    std::vector<unenc_compound> polyatomics;

    for (int i = POLYATOMIC_START_IDX; i < POLYATOMIC_START_IDX + 1; i++) {
      polyatomics.emplace_back(chems[i]);
    }
  
    auto polyatomics_encoded = (::Preprocess::encode_compounds(polyatomics)).move_ptr();

    for (int i{}; i < encoded_simplified.get_rows(); i++) {
      for (int j = 0; j < polyatomics.size(); j++) {
	bool contains_polyatomic{ true };
      
	for(int k = 0; k < TOTAL_CHEMS; k++) {
	  int polyatomic_n_el = std::round(polyatomics_encoded[j * TOTAL_CHEMS + k] / CHEM_SCALE_FACTOR);
	  int compound_n_el = std::round(encoded_unsimplified[i * TOTAL_CHEMS + k] / CHEM_SCALE_FACTOR);
	  if (compound_n_el < polyatomic_n_el) {
	    contains_polyatomic = false;
	    break;
	  }
	}

	if (contains_polyatomic) {
	  int poly_idx = POLYATOMIC_START_IDX + j;
	  encoded_unsimplified[i * TOTAL_CHEMS + poly_idx] += CHEM_SCALE_FACTOR;
	
	  for(int k = 0; k < TOTAL_CHEMS; k++) {
	    encoded_unsimplified[i * TOTAL_CHEMS + k] -= polyatomics_encoded[j * TOTAL_CHEMS + k];
	  }
	}
      }
    }

    return Matrix<double>(encoded_simplified.get_rows(), encoded_simplified.get_cols(), std::move(encoded_unsimplified));
  }

  // ---- Check if a compound is for PTR mode ----
  bool uses_reagant_ion(const ::std::span<double> &encoded_unsimplified, unenc_compound reagant_ion) {
    auto *cm = ChemMap::get_chem_map();
    auto &mask = cm->get_reagant_ion_mask(reagant_ion);

    size_t ctr{};
    
    for (int i{}; i < POLYATOMIC_START_IDX; i++) {
      if (ctr >= mask.length || mask.indeces[ctr] != i) {
        if (encoded_unsimplified[i] > 0.0)
	  return false;
      } else {
	ctr++;
      }
    }
  
    return true;
  }
}
