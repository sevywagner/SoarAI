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
    _masses.push_back(18.998403);
    _masses.push_back(28.0855);
    _masses.push_back(31.97207117);
    _masses.push_back(35.4532);
    _masses.push_back(-0.000548);
    _masses.push_back(18.03437412);
    _masses.push_back(17.02654909);
    _masses.push_back(19.01838971);
    _masses.push_back(18.01056468);

    for (int i = 0; i < TOTAL_CHEMS; i++) {
      _chem_masses[_chemicals[i]] = _masses[i];
      _chem_idx[_chemicals[i]] = i;
    }

    auto mask_ptr = std::make_unique<uint32_t[]>(4);
    std::iota(mask_ptr.get(), mask_ptr.get() + 4, 0);
    _PTR_mask = Mask<IDX, uint32_t>(4, 0, std::move(mask_ptr));
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

  // ---- Get masses vec ----
  std::vector<double> ChemMap::get_masses() {
    return _masses;
  }

  // ---- Get chem vec ----
  std::vector<std::string> ChemMap::get_chems() {
    return _chemicals;
  }

  // ----------------
  // Chem utilities
  // ----------------

  // ---- Compare compounds ----
  bool compounds_are_equal(const Matrix<double> &compound1, const Matrix<double> &compound2) {
    if (compound1.get_rows() != compound2.get_rows() || compound1.get_cols() > 1 || compound1.get_cols() > 1) {
      std::cerr << "Compound comparison err -- misaligned dims" << std::endl;
      exit(1);
    }

    for (size_t i = 0; i < TOTAL_CHEMS; i++) {
      if (std::round(CHEM_SCALE_FACTOR / compound1.get(i, 0)) != std::round(CHEM_SCALE_FACTOR / compound2.get(i, 0))) {
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
  double get_compound_mass(const Matrix<double> &compound) {
    if (compound.get_rows() != TOTAL_CHEMS) {
      std::cerr << "Malformed compound encoding in get_compound_mass" << std::endl;
      exit(1);
    }

    auto *cm = ChemMap::get_chem_map();
    auto masses = cm->get_masses();
  
    double total_mass;
    for (size_t i = 0; i < TOTAL_CHEMS; i++) {
      total_mass += (masses[i] * (compound.get(i, 0) / CHEM_SCALE_FACTOR));
    }

    return total_mass;
  }

  // ---- Check the 4 criterea of a compound ----
  CritereaCheckRes check_criterea(const Matrix<double> &compound) {
    std::vector<bool> criterea_mask;
    for (int i = 0; i < 4; i++)
      criterea_mask.push_back(true);
  
    bool did_pass{true};
    auto *cm = ChemMap::get_chem_map();

    bool contains_nh4 = compound.get(cm->get_idx("NH4"), 0) > 0;
    bool contains_nh3 = compound.get(cm->get_idx("NH3"), 0) > 0;
    uint32_t n_hydrogen = static_cast<uint32_t>(compound.get(cm->get_idx("H"), 0) / CHEM_SCALE_FACTOR);
    uint32_t n_carbon = static_cast<uint32_t>(compound.get(cm->get_idx("C"), 0) / CHEM_SCALE_FACTOR);
    uint32_t n_nitrogen = static_cast<uint32_t>(compound.get(cm->get_idx("N"), 0) / CHEM_SCALE_FACTOR);

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

    for (int i = 0; i < encoded_simplified.get_rows(); i++) {
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
  bool is_PTR(const Matrix<double> &encoded_unsimplified) {
    auto *cm = ChemMap::get_chem_map();
    for (int i = PTR_END_IDX; i < POLYATOMIC_START_IDX; i++) {
      if (encoded_unsimplified.get(i, 0) > 0)
	return false;
    }
  
    return true;
  }
}
