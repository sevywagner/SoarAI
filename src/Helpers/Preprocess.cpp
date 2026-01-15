#include "Preprocess.h"

using namespace CNum::DataStructs;
using namespace CNum::Multithreading;
using namespace Chem;

// -----------------------
// Compound Manipulation
// -----------------------

namespace Preprocess {
  // ---- Turn compound string into an encoded compound ----
  Matrix<double> encode_compounds(const ::std::vector<unenc_compound> &compound_strings) {
    size_t ec_size = compound_strings.size() * TOTAL_CHEMS;
    auto encoded_compounds_ptr = ::std::make_unique<double[]>(ec_size);
    ::std::fill(encoded_compounds_ptr.get(), encoded_compounds_ptr.get() + ec_size, 0.0);

    auto *cm = ChemMap::get_chem_map();

    size_t comp_ctr = 0;
    for (auto compound: compound_strings) {
      auto open_per_it = compound.val.begin();
    
      while ((open_per_it = ::std::find(compound.val.begin(), compound.val.end(), '(')) != compound.val.end()) {
	auto end_per_it = ::std::find(open_per_it, compound.val.end(), ')');
	::std::string polyatomic = { open_per_it + 1, end_per_it };

	auto after_polyatomic = end_per_it + 1;
	::std::string num{""};
	while (after_polyatomic != compound.val.end() && ::std::isdigit(static_cast<unsigned char>(*after_polyatomic)) != 0) {
	  num += *after_polyatomic;
	  after_polyatomic++;
	}

	size_t idx = comp_ctr * TOTAL_CHEMS + cm->get_idx(polyatomic);
	encoded_compounds_ptr[idx] = (num == "" ? 1 : stoi(num)) *  CHEM_SCALE_FACTOR;
	compound.val = { end_per_it, compound.val.end() };
      }

      for (int i{}; i < compound.val.size(); i++) {
	::std::string chemical;
	chemical.push_back(compound.val[i]);
	if (i < compound.val.size() - 1 && ::std::islower(static_cast<unsigned char>(compound.val[i + 1])) != 0) {
	  chemical += compound.val[i + 1];
	}

	if (::std::isupper(static_cast<unsigned char>(chemical[0])) == 0) {
	  continue;
	}

	::std::string num{""};
	for (int j = i + chemical.size(); j < compound.val.size(); j++) {
	  if (::std::isdigit(static_cast<unsigned char>(compound.val[j])
			     ) == 0) {
	    break;
	  }

	  num += compound.val[j];
	}

	if (num == "") {
	  num = "1";
	}

	double val;
	try {
	  val = stof(num) * CHEM_SCALE_FACTOR;
	} catch(...) {
	  throw ::std::runtime_error("Encode compounds error -- converting " + num + " to string");
	}
	encoded_compounds_ptr[comp_ctr * TOTAL_CHEMS + cm->get_idx(chemical)] = val;
      }

      comp_ctr++;
    }
  
    return Matrix<double>(compound_strings.size(), TOTAL_CHEMS, ::std::move(encoded_compounds_ptr));
  }

  // ---- Recursive function for finding all possible elemental combos with a total mass close to the m/z ----
  void all_possible_el_recurse(::std::vector< ::std::vector<size_t> > &res,
			       ::std::vector<double> &theoretical_compound_masses,
			       const ::std::vector<double> &masses,
			       ::std::vector<size_t> &temp,
			       const ReagantIonMask &mask,
			       double mass,
			       double target_mass,
			       size_t idx = 0) {
    constexpr double pruning_tolerance = 1.5;
    constexpr double ppm_tolerance = 50;
  
    double theoretical_mass = target_mass - mass;

    if (::std::abs(Chem::get_ppm(target_mass, theoretical_mass)) <= ppm_tolerance) {    
      res.push_back(temp);
      theoretical_compound_masses.push_back(theoretical_mass);
      return;
    }
  
    for (int i = idx; i < mask.length; i++) {
      if (mass - masses[mask.indeces[i]] <= -pruning_tolerance) {
	continue;
      }
    
      temp.push_back(mask.indeces[i]);
      all_possible_el_recurse(res,
			      theoretical_compound_masses,
			      masses,
			      temp,
			      mask,
			      mass - masses[mask.indeces[i]],
			      target_mass,
			      i);
      temp.pop_back();
    }
  }

  // ---- Find all possible elemental combos with a total mass close to the m/z ----
  CompoundPermutations all_possible_elemental_combo(double mass, unenc_compound reagant_ion) {
    auto *cm = ChemMap::get_chem_map();
    const auto &masses = cm->get_masses();
    const auto &mask = cm->get_reagant_ion_mask(reagant_ion);
    
    ::std::vector<size_t> temp;
    ::std::vector<double> theoretical_compound_masses;
    ::std::vector< ::std::vector<size_t> > res;

    all_possible_el_recurse(res,
			    theoretical_compound_masses,
			    masses,
			    temp,
			    mask,
			    mass,
			    mass);
  
    auto encoded_compounds = ::std::make_unique<double[]>(TOTAL_CHEMS * res.size());
  
    for (int i{}; i < TOTAL_CHEMS; i++) {
      for (int j = 0; j < res.size(); j++) {
	encoded_compounds[j * TOTAL_CHEMS + i] = static_cast<double>(::std::count(res[j].begin(), res[j].end(), i))  * CHEM_SCALE_FACTOR;
      }
    }
  
    return { Matrix<double>(res.size(), TOTAL_CHEMS, ::std::move(encoded_compounds)),
	     ::std::move(theoretical_compound_masses) };
  }

  // ---- Take encoded compounds and decode them back into strings ----
  ::std::vector<unenc_compound> decode_compounds(const Matrix<double> &encoded_compounds) {
    ::std::vector<unenc_compound> decoded_compounds;
    decoded_compounds.reserve(encoded_compounds.get_rows());
    
    auto *cm = ChemMap::get_chem_map();
    const auto &chems = cm->get_chems();
    const auto &proper_order = cm->get_proper_ordering();

    for (int i{}; i < encoded_compounds.get_rows(); i++) {
      ::std::string compound_str{""};
      for (int j = 0; j < TOTAL_CHEMS; j++) {
	uint8_t idx = proper_order[j];
	double num_el_j = encoded_compounds.get(i, idx);
	if (num_el_j > 0) {
	  // Check if chem is a polyatomic (if its length is > 2)
	  if (chems[idx].size() > 2) {
	    compound_str += "(" + chems[idx] + ")";
	    continue;
	  }
	
	  compound_str += chems[idx];
	  if (num_el_j > CHEM_SCALE_FACTOR) {
	    compound_str += ::std::to_string((uint32_t) (num_el_j / CHEM_SCALE_FACTOR));
	  }
	}
      }

      decoded_compounds.emplace_back(compound_str);
    }
  
    return decoded_compounds;
  }

  // ---- Take out polyatomics from encoded compound and adjust the polyatomic chemicals accordingly ----
  Matrix<double> simplify_compounds(const Matrix<double> &unsimplified_compounds) {
    uint32_t n_compounds = unsimplified_compounds.get_rows();

    size_t simplified_size = n_compounds * TOTAL_CHEMS;
    auto simplified_ptr = ::std::make_unique<double[]>(simplified_size);
    ::std::fill(simplified_ptr.get(), simplified_ptr.get() + simplified_size, 0.0);

    auto *cm = ChemMap::get_chem_map();
    const auto &chems = cm->get_chems();

    ::std::vector<unenc_compound> polyatomics;

    for (size_t i = POLYATOMIC_START_IDX; i < TOTAL_CHEMS; i++) {
      polyatomics.emplace_back(chems[i]);
    }

    auto polyatomics_encoded = encode_compounds(polyatomics);

    for (size_t i{}; i < n_compounds; i++) {
      for (int j = 0; j < POLYATOMIC_START_IDX; j++) {
	simplified_ptr[i * TOTAL_CHEMS + j] = unsimplified_compounds.get(i, j);
      }
    
      for (size_t j = POLYATOMIC_START_IDX; j < TOTAL_CHEMS; j++) {
	double n_poly = unsimplified_compounds.get(i, j);
	if (n_poly > 0.0) {
	  auto poly_j = polyatomics_encoded.get(ROW, j - POLYATOMIC_START_IDX) * ::std::round(n_poly / CHEM_SCALE_FACTOR);
	  for (int k = 0; k < TOTAL_CHEMS; k++) {
	    simplified_ptr[i * TOTAL_CHEMS + k] += poly_j.get(k, 0);
	  }
	}
      }
    }

    return Matrix<double>(n_compounds, TOTAL_CHEMS, ::std::move(simplified_ptr));
  }

  // ---- Prepare data for training and inference ----
  MSData mz_to_data(double mz, unenc_compound ion, size_t n_features) {
    auto apc = all_possible_elemental_combo(mz, ion);
  
    apc.compounds = Chem::factor_polyatomics(apc.compounds);
    auto *all_possible_compounds = &apc.compounds;
    auto *theoretical_masses = &apc.masses;
    size_t total_possible_compounds = all_possible_compounds->get_rows();

    auto data = ::std::make_unique<double[]>(total_possible_compounds * n_features);
    auto ppms = ::std::make_unique<double[]>(total_possible_compounds);
      
    for (size_t j = 0; j < total_possible_compounds; j++) {
      auto permutation_matrix = all_possible_compounds->get(ROW, j);
      auto crit_check = Chem::check_criterea(permutation_matrix);

      auto permutation = permutation_matrix.move_ptr();
      auto *crit_mask = &crit_check.crit_mask;
      ::std::move(crit_mask->begin(), crit_mask->end(), data.get() + (j * n_features));
      ::std::move(permutation.get(), permutation.get() + TOTAL_CHEMS, data.get() + (j * n_features + N_CRITEREA));
      
      double ppm = Chem::get_ppm(mz, theoretical_masses->at(j));
      data[j * n_features + n_features - 1] = ppm;
      ppms[j] = ppm;
    }

    return { Matrix<double>(total_possible_compounds, n_features, ::std::move(data)),
	     ::std::move(apc.compounds),
	     Matrix<double>(total_possible_compounds, 1, ::std::move(ppms)) };
  }

  // ------------------
  // Dataset Creation
  // ------------------

  // ---- Go through peak list and collect compound strings and m/z values ----
  PeakListData PrepareDataset::parse_peak_list(::std::string path) {
    ::std::ifstream istream(path);

    if (!istream.is_open()) {
      throw ::std::runtime_error("Parse peak list error - opening file failed");
    }

    constexpr uint8_t ion_col = 3, x0_col = 4;
    ::std::vector<double> x0;
    ::std::vector<unenc_compound> assigned_formulas;

    ::std::string line;
    getline(istream, line, '\n');
    while (getline(istream, line, '\n')) {
      ::std::stringstream ss(line);
      ::std::string segment;
      for (int i{}; i < ion_col; i++) {
	getline(ss, segment, '\t');
      }

      assigned_formulas.emplace_back(segment);
      getline(ss, segment, '\t');
      try {
	x0.push_back(stof(segment));
      } catch(...) {
	throw ::std::runtime_error("Parse peak list error -- converting " + segment + " to double");
      }
    }

    return { ::std::move(x0), ::std::move(assigned_formulas) };
  }

  // ---- Prepare and output training/inference ready data for all peaks in a peak list to a file ----
  Bias PrepareDataset::create_combo_file(::std::string output_path,
					 ::std::string unidentified_combos_path,
					 ::Chem::unenc_compound reagant_ion,
					 const PeakListData &peak_list_data,
					 const Matrix<double> &encoded_unsimplified,
					 const Matrix<double> &encoded_simplified,
					 int n_threads) {
    if (n_threads == 0)
      throw ::std::invalid_argument("Combo file creation error -- n_threads cannot be 0");
    
    auto *x0 = &peak_list_data.mz;
    auto *assigned_formulas = &peak_list_data.compound_strings;

    size_t total_assigned = assigned_formulas->size();
    size_t assigned_per_thread = (total_assigned + n_threads - 1) / n_threads;
    
    if (total_assigned < n_threads)
      n_threads = static_cast<int>(total_assigned);
    
    ::std::ofstream ostream(output_path);
    
    if (!ostream.is_open()) {
      throw ::std::runtime_error("Combo file creation error -- error opening output file");
    }

    ::std::ofstream not_found_ostream(unidentified_combos_path);

    if (!not_found_ostream.is_open()) {
      throw ::std::runtime_error("Combo file creation error -- error opening unidentified compounds output file");
    }

    not_found_ostream << "Compound,PPM" << ::std::endl;

    ::std::vector< ::std::future<size_t> > workers;
    workers.reserve(n_threads);
    auto *tp = ThreadPool::get_thread_pool();

    ::std::mutex ostream_mtx;

    auto *cm = ChemMap::get_chem_map();

    for (uint8_t thread_num{}; thread_num < n_threads; thread_num++) {
      workers.push_back(tp->submit< size_t >([&, thread_num] (arena_t *arena) {
	size_t start = thread_num * assigned_per_thread;
	size_t end = ::std::min(start + assigned_per_thread, total_assigned);

	::std::ostringstream data("");
	::std::ostringstream unidentified("");

	size_t total_samples{};
	for (size_t i = start; i < end; i++) {
	  auto encoded_unsimplified_assigned_formula_view = encoded_unsimplified.get_row_view(i);
	  auto encoded_simplified_assigned_formula_view = encoded_simplified.get_row_view(i);

	  ::Chem::unenc_compound ion = reagant_ion;
	  
	  // change everything to work with views
	  if (reagant_ion.val.empty()) {
	    ion = cm->find_reagant_ion(encoded_unsimplified_assigned_formula_view);
	    if (ion.val.empty()) continue;
	  }
	
	  auto apc = all_possible_elemental_combo(x0->at(i), ion);
      
	  auto all_possible_compounds_unsimplified = Chem::factor_polyatomics(apc.compounds);
	  auto *all_possible_compounds_simplified = &apc.compounds;
	  auto *theoretical_masses = &apc.masses;

	  bool is_found{ false };
      
	  for (size_t j{}; j < all_possible_compounds_simplified->get_rows(); j++) {
	    auto permutation_simplified_view = all_possible_compounds_simplified->get_row_view(j);
	    auto permutation_unsimplified_view = all_possible_compounds_unsimplified.get_row_view(j);
	    bool are_same = Chem::compounds_are_equal(permutation_simplified_view, encoded_simplified_assigned_formula_view);

	    if (are_same) {
	      if (is_found) continue;
	      is_found = true;
	    }

	    auto crit_check = Chem::check_criterea(permutation_unsimplified_view);
	    auto *crit_mask = &crit_check.crit_mask;
	    for (int k{}; k < N_CRITEREA; k++) {
	      data << ::std::to_string((int) crit_mask->at(k)) + ",";
	    }
      
	    for (int k{}; k < TOTAL_CHEMS; k++) {
	      data << ::std::to_string(permutation_unsimplified_view[k]) + ",";
	    }
	    
	    data << ::std::to_string(Chem::get_ppm(x0->at(i), theoretical_masses->at(j)))
		 << "," + ::std::to_string(are_same) << ::std::endl;
	    
	    total_samples++;

	    // output and reset buffer if data buffer size > 1kb
	    if (data.tellp() > 1 << 10 || unidentified.tellp() > 1 << 10) {
	      {
		::std::lock_guard<::std::mutex> lg(ostream_mtx);
		ostream << data.str();
		not_found_ostream << unidentified.str();
	      }
	  
	      data.str("");
	      data.clear();
	      unidentified.str("");
	      unidentified.clear();
	    }
	  }

	  if (!is_found) {
	    double mass = Chem::get_compound_mass(encoded_simplified_assigned_formula_view);
	    double ppm = Chem::get_ppm(x0->at(i), mass);
	    unidentified << assigned_formulas->at(i).val << "," << ppm << ::std::endl;
	  }
	}

	{
	  ::std::lock_guard<::std::mutex> lg(ostream_mtx);
	  ostream << data.str();
	  not_found_ostream << unidentified.str();
	}
	
	data.str("");
	data.clear();
	unidentified.str("");
	unidentified.clear();
      
	return total_samples;
      }));
    }

    size_t total_samples{};
    for (auto &f: workers) {
      total_samples += f.get();
    }

    return { total_assigned, total_samples - total_assigned };
  }

  // ---- Downsample negative samples in a dataset *usually not a good idea + right now the ----
  // ---- dataset is too small but could possibly be used later for much larger ones*       ----
  void PrepareDataset::negative_sample_reduction(::std::string path) {
    constexpr int percent_reduced = 7;
    ::std::ifstream is(path);
    ::std::ofstream os("./data/reduced_negative_sample.csv");

    if (!is.is_open()) {
      throw ::std::runtime_error("Negative sample reduction error -- opening input file failed");
    }

    if (!os.is_open()) {
      throw ::std::runtime_error("Negative sample reduction error -- opening output file failed");
    }

    auto &rng = ::CNum::Utils::Rand::RandomGenerator::instance();
    ::std::uniform_int_distribution<int> dist(0, 100);

    ::std::string line;
    while (getline(is, line, '\n')) {
      if (line[line.size() - 1] == '0' && dist(rng) > percent_reduced) {
	continue;
      }
      
      os << line << ::std::endl;
    }
  }

  // ---- Seperate positive and negative samples ----
  void data_sample_seperation(::std::string path,
			      ::std::vector< ::std::vector<double> > &data_ones,
			      ::std::vector< ::std::vector<double> > &data_zeros) {
    ::std::ifstream is(path);

    if (!is.is_open()) {
      throw ::std::runtime_error("Sample seperation error -- opening input file failed");
    }

    ::std::string line;
    while (getline(is, line, '\n')) {
      ::std::vector<double> sample;
      ::std::istringstream iss(line);
      ::std::string segment;

      while (getline(iss, segment, ',')) {
	try {
	  sample.push_back(stod(segment));
	} catch (...) {
	  throw ::std::runtime_error("Train test split error -- error converting " + segment + " to double");
	}
      }

      if (sample.back() == 1.0) {
	data_ones.push_back(sample);
      } else {
	data_zeros.push_back(sample);
      }
    }
  }

  // ---- Seperate train and test datasets ----
  void PrepareDataset::train_test_split(::std::string combo_file_path,
					::std::string output_path,
					size_t n_test_pos,
					size_t n_test_neg) {
    ::std::vector< ::std::vector<double> > data_ones;
    ::std::vector< ::std::vector<double> > data_zeros;
    data_sample_seperation(combo_file_path, ::std::ref(data_ones), ::std::ref(data_zeros));
    assert(!data_ones.empty());
    assert(!data_zeros.empty());

    auto &rng = ::CNum::Utils::Rand::RandomGenerator::instance();
    
    ::std::array< ::std::vector< ::std::vector<double> >, 2 > data = {
      ::std::move(data_ones),
      ::std::move(data_zeros)
    };
    
    ::std::array< ::std::unordered_map<int, int>, 2 > used_idx{};
    ::std::array<size_t, 2> test_set_sizes = { n_test_pos, n_test_neg };
    
    std::string test_path = output_path + "_test_combos.csv";
    std::string train_path = output_path + "_train_combos.csv";
    ::std::ofstream os(test_path);

    if (!os.is_open()) {
      throw ::std::runtime_error("Train test split error -- problem opening test output file");
    }

    for (int x = 0; x < 2; x++) {
      ::std::uniform_int_distribution<int> dist(0, data[x].size() - 1);

      size_t i{};
      while (i < test_set_sizes[x]) {
	size_t idx = dist(rng);

	if (used_idx[x][idx] > 0) {
	  continue;
	}
	
	used_idx[x][idx]++;
    
	::std::string output{""};
	for (int j = 0; j < data[x][0].size(); j++) {
	  output += ::std::to_string(data[x][idx][j]) + ",";
	}
    
	output.pop_back();

	os << output << ::std::endl;
	i++;
      }
    }
    
    os.close();
    os = ::std::ofstream(train_path);

    if (!os.is_open()) {
      throw ::std::runtime_error("Train test split error -- problem opening train output file");
    }

    for (int x = 0; x < 2; x++) {
      for (size_t i{}; i < data[x].size(); i++) {
	if (used_idx[x][i] > 0) {
	  continue;
	}

	::std::string output{""};
	for (size_t j = 0; j < data[x][0].size(); j++) {
	  output += ::std::to_string(data[x][i][j]) + ",";
	}
    
	output.pop_back();

	os << output << ::std::endl;
      }
    }
  }

  // ---- Get the number of positive and number of negative samples ----
  Bias check_bias(Matrix<double> &row_matrix) {
    size_t ones{0};
    size_t zeros(0);

    for (size_t i{}; i < row_matrix.get_rows(); i++) {
      if (static_cast<uint8_t>(row_matrix.get(i, 0)) == 0) {
	zeros++;
      } else {
	ones++;
      }
    }

    return { ones, zeros };
  }

  void PrepareDataset::peak_list_train_test_split(::std::string peak_list_path,
						  ::std::string output_path,
						  double test_percentage) {
    ::std::ifstream istream(peak_list_path);
    if (!istream.is_open()) {
      throw ::std::runtime_error("Peak list train/test split error -- Could not open peak list");
    }

    ::std::ofstream train_out(output_path + "_train.txt");
    if (!train_out.is_open()) {
      throw ::std::runtime_error("Peak list train/test split error -- Could not open train output list");
    }

    ::std::ofstream test_out(output_path + "_test.txt");
    if (!test_out.is_open()) {
      throw ::std::runtime_error("Peak list train/test split error -- Could not open test output list");
    }
    
    ::std::vector<::std::string> lines;
    ::std::string line, headers;
    getline(istream, headers, '\n');
    
    
    while (getline(istream, line, '\n')) {
      lines.push_back(line);
    }

    auto &rng = ::CNum::Utils::Rand::RandomGenerator::instance();
    ::std::shuffle(lines.begin(), lines.end(), rng);
    
    for (size_t i{}; i < lines.size() * (1 - test_percentage); i++) {
      train_out << lines[i] << ::std::endl;
    }

    for (size_t i = lines.size() * (1 - test_percentage); i < lines.size(); i++) {
      test_out << lines[i] << ::std::endl;
    }

    istream.close();
    train_out.close();
    test_out.close();
  }

  // switch masks to spans, create individual decode and encode funcs
  
  void PrepareDataset::split_by_reagant_ion(::std::string peak_list_path, ::std::string output_path) {
    auto *cm = ChemMap::get_chem_map();
    const auto &reagant_ions = cm->get_reagant_ions();
    ::std::unordered_map<::std::string, ::std::ofstream> outputs;

    for (auto &ion: reagant_ions) {
      auto path = output_path + ion.val + ".txt";
      outputs[ion.val] = ::std::ofstream(path);

      if (!outputs[ion.val].is_open()) {
	throw ::std::runtime_error("Split by reagant ion error -- couldn't open " + ion.val + " output file");
      }
    }

    ::std::ifstream is(peak_list_path);
    if (!is.is_open()) {
      throw ::std::runtime_error("Split by reagant ion error -- couldn't open input file");
    }

    constexpr uint8_t ion_col = 3;
    
    ::std::string line;
    getline(is, line, '\n');
    for (auto &ion: reagant_ions)
      outputs[ion.val] << line << ::std::endl;
    
    while (getline(is, line, '\n')) {
      ::std::istringstream iss(line);
      unenc_compound c{ "" };
      for (int i{}; i < ion_col; i++) {
	getline(iss, c.val, ',');
      }

      ::std::vector< unenc_compound > dc = { c };
      auto ec = ::Preprocess::simplify_compounds(::Preprocess::encode_compounds(dc));
      auto ec_view = ec.get_row_view(0);

      auto ion = cm->find_reagant_ion(ec_view);
      if (ion.val.empty())
	continue;

      outputs[ion.val] << line << ::std::endl;
    }
  }

  ::CNum::Model::Tree::SubsampleFunction get_subsample_func(::std::vector<size_t> &ones_indeces,
							    ::std::unordered_set<size_t> &ones_indeces_set) {
    ::CNum::Model::Tree::SubsampleFunction subsample = [&ones_indeces, &ones_indeces_set] (size_t *pos_ptr,
											   size_t low,
											   size_t high,
											   size_t n_samples,
											   const Matrix<double> &y) -> void {
      auto &rng = ::CNum::Utils::Rand::RandomGenerator::instance();
      ::std::uniform_int_distribution<uint64_t> dist(low, high - 1);
      if (ones_indeces.empty()) {
	for (size_t i{}; i < y.get_rows(); i++) {
	  if (::std::round(y.get(i, 0)) == 1) {
	    ones_indeces.push_back(i);
	    ones_indeces_set.insert(i);
	  }
	}
      }

      ::std::unordered_set<size_t> used_indeces = ones_indeces_set;
      ::std::copy(ones_indeces.begin(), ones_indeces.end(), pos_ptr);

      size_t n = ones_indeces.size();
      while (n < n_samples) {
	size_t rand = dist(rng);
	if (used_indeces.contains(rand)) continue;

	pos_ptr[n] = rand;
	n++;
      }
    };

    return subsample;
  }
}
