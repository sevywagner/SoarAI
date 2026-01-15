#include "InferenceAPI.h"

using namespace CNum::DataStructs;
using namespace CNum::Model::Tree;

namespace InferenceAPI {
  char *python_executable_path = NULL;
  ::std::string graph_upload_dir = "";
  ::std::string peak_output_dir = "";
  
  // -----------------
  // File Validation
  // -----------------
  static void trim(::std::string &s) {
    while (!s.empty() && ::std::isspace(s.back()))
      s.pop_back();
  }
  
  static void validate_file_upload(crow::multipart::part p) {
    if (p.body.size() > MAX_FILE_SIZE)
      throw ::std::runtime_error("Uploaded CSV file is too large. The limit is 7mb. This file is " + ::std::to_string(p.body.size() / ::std::pow(2, 10)) + "mb");

    if (p.body.find('\0') != ::std::string::npos)
      throw ::std::runtime_error("Uploaded CSV has incorrect format. Null bytes found.");
    
    ::std::istringstream iss(p.body);
    ::std::string line;
    getline(iss, line, '\n');
    if (line.empty())
	throw ::std::runtime_error("Uploaded CSV has incorrect format. Empty line found.");
    trim(line);

    // check for UTF-8 BOM
    if (line.size() >= 3 && line[0] == 0xEF && line[1] == 0xBB && line[2] == 0xBF)
      line = line.substr(3);
    
    if (line != "mz_base" && line != "mz_av")
      throw ::std::runtime_error("Uploaded CSV has incorrect format. \"mz_base\" or \"mz_av\" required as 1st line with no extra whitespace");

    while (getline(iss, line, '\n')) {
      if (line.empty())
	throw ::std::runtime_error("Uploaded CSV has incorrect format. Empty line found.");
      trim(line);
      
      
      ::std::istringstream seg_iss(line);
      ::std::string segment;
      bool is_first_it{ true };
      while (getline(seg_iss, segment, ',')) {
	if (!is_first_it)
	  throw ::std::runtime_error("Uploaded CSV has incorrect format. Only one column is allowed. ");

	bool contains_digit{ false };
	for (unsigned char c: segment) {
	  if (::std::isdigit(c)) {
	    contains_digit = true;
	    break;
	  }
	}

	if (!contains_digit)
	  throw ::std::runtime_error("Uploaded CSV has incorrect format. Digitless line found");
	
	double val;
	try {
	  val = ::std::stod(segment);
	} catch (...) {
	  throw ::std::runtime_error("Bad float value in uploaded CSV. ");
	}

	if (::std::isnan(val) || ::std::isinf(val))
	  throw ::std::runtime_error("NaN/INF float value in uploaded CSV. ");

	is_first_it = false;
      }
    }
  }

  
  // ------------------
  // Data Processing
  // ------------------
  
  // ---- Preprocess data before using it as input to model ----
  Matrix<double> preprocess_func(crow::json::rvalue &req_body, crow::json::wvalue &res_body, Storage &storage) {
    auto mz_values = req_body["mz_array"];
    ::Chem::unenc_compound ion(req_body["reagant_ion"].s());

    // SoarAI currently only supports ammonium and nitrosyl reagant ions
    if (ion.val == "def" || !(ion.val == "NH4" || ion.val == "NO"))
      throw ::std::invalid_argument("Invalid reagent ion: " + ion.val);

    ::std::vector< Matrix<double> > model_data_matrices;
    ::std::vector< Matrix<double> > encoded_compounds_matrices;
    ::std::vector< Matrix<double> > ppm_matrices;
    model_data_matrices.reserve(mz_values.size());
    encoded_compounds_matrices.reserve(mz_values.size());
    ppm_matrices.reserve(mz_values.size());
  
    size_t total_rows{ 0 };
    res_body["critereaEncodings"] = crow::json::wvalue::list();

    for (size_t i = 0; i < mz_values.size(); i++) {
      auto data = Preprocess::mz_to_data(mz_values[i].d(), ion);

      size_t temp = total_rows;
      total_rows += data.model_data.get_rows();
    
      for (size_t j = temp; j < total_rows; j++) {
	res_body["critereaEncodings"][j] = crow::json::wvalue::list();
	for (int k = 0; k < 4; k++) {
	  res_body["critereaEncodings"][j][k] = data.model_data.get(j - temp, k);
	}
      }
    
      model_data_matrices.push_back(::std::move(data.model_data));
      encoded_compounds_matrices.push_back(::std::move(data.encoded_compounds));
      ppm_matrices.push_back(::std::move(data.ppms));
    }

    auto ec = Matrix<double>::combine_vertically(encoded_compounds_matrices, total_rows);
    auto ppms = Matrix<double>::combine_vertically(ppm_matrices, total_rows);
    auto res = Matrix<double>::combine_vertically(model_data_matrices, total_rows);
  
    storage.encoded_compounds = ::std::move(ec);
    storage.ppms = ::std::move(ppms);
  
    return res;
  }

  // ---- Postprocess data and save in the response ----
  void postprocess_func(Matrix<double> &preds, crow::json::wvalue &res_body, Storage &storage) {
    Postprocess::sort_preds(::std::ref(storage.encoded_compounds),
			    ::std::ref(preds),
			    ::std::ref(storage.ppms));

    auto decoded_unsimplified_compounds = Preprocess::decode_compounds(storage.encoded_compounds);
    auto decoded_simplified_compounds = Preprocess::decode_compounds(Preprocess::simplify_compounds(storage.encoded_compounds));
  
    res_body["scores"] = crow::json::wvalue::list();
    res_body["compounds"] = crow::json::wvalue::list();
    res_body["uCompounds"] = crow::json::wvalue::list();
    res_body["ppms"] = crow::json::wvalue::list();

    for (size_t i = 0; i < preds.get_rows(); i++) {
      res_body["scores"][i] = preds.get(i, 0) * 100; // convert to percentage
      res_body["compounds"][i] = decoded_simplified_compounds[i].val;
      res_body["uCompounds"][i] = decoded_unsimplified_compounds[i].val;
      res_body["ppms"][i] = storage.ppms.get(i, 0);
    }
  }

  // --------------------
  // Spectra Processing
  // --------------------

  // ---- Make the filename of uploads and stored files unique ----
  static ::std::string make_unique_filename(const ::std::string prefix, const ::std::string ext) {
    auto &rng = ::CNum::Utils::Rand::RandomGenerator::instance();
    ::std::string rand = ::std::to_string(rng());
    ::std::string time = ::std::to_string(static_cast<long long>(::std::time(nullptr)));
    return prefix + time + "_" + rand + ext;
  }

  // ---- Add the break in the table ----
  static void page_break(::std::ostringstream &oss) {
    oss.put('+').fill('-');
    for (auto w: TABLE_WIDTHS) {
      oss.width(w + TABLE_MARGIN * 2);
      oss << "-" << "+";
    }

    oss.put('\n').fill(' ');
  }

  // ---- Output a row of the table to a string stream ----
  static void print_row(::std::ostringstream &oss, const ::std::array<::std::string, 3> values) {
    oss << "|";
    for (int i = 0; i < TABLE_N_COLS; i++) {
      oss << " ";
      oss.width(TABLE_WIDTHS[i]);
      oss << values[i] << " |";
    }

    oss << ::std::endl;
  }

  // ---- Print the column headers of the table ----
  static void print_header(::std::ostringstream &oss) {
    page_break(oss);
    print_row(oss, TABLE_HEADERS);
    page_break(oss);
  }

  // ---- Take in a mass spectra, find and fit peaks, make predictions, and postprocess ----
  void process_graph(const crow::request &req, crow::response &res, GBModel<XGTreeBooster> *model) {
    const auto content_type = req.get_header_value("Content-Type");
    if (content_type.find("multipart/form-data") == ::std::string::npos) {
      res = crow::response(500, "Content-Type must be multipart/form-data");
      res.end();
      return;
    }
  
    crow::multipart::message msg(req);
    ::std::array<crow::multipart::part, N_FILES> parts;
    ::std::array<::std::string, N_FILES> filenames;
    
    ::Chem::unenc_compound reagentIon{ msg.get_part_by_name("reagentIon").body };
    parts[0] = msg.get_part_by_name("base");
    parts[1] = msg.get_part_by_name("av");

    for (int i = 0; i < N_FILES; i++) {
      const auto cd = parts[i].get_header_object("Content-Disposition");
      auto it = cd.params.find("filename");
      if (it != cd.params.end() && !it->second.empty())
	filenames[i] = make_unique_filename(::InferenceAPI::graph_upload_dir, ".csv");

      validate_file_upload(parts[i]);
      ::std::ofstream os(filenames[i], ::std::ios::binary);

      if (!os.is_open()) {
	::std::cerr << "Error in /process_graph - Failed to open file" << ::std::endl;
	res = crow::response(500, "Error in /process_graph - Failed to open " + filenames[i]);
	res.end();
	return;
      }
      
      os << parts[i].body;
      os.close();
    }

    auto peaks_file_path = make_unique_filename(peak_output_dir, ".txt");
    execute_peak_fit_bin(filenames[0], filenames[1], peaks_file_path);

    ::std::ifstream is(peaks_file_path, ::std::ios::binary);
    if (!is.is_open()) {
      res = crow::response(500, "Error in /process-graph - couldn't open " + peaks_file_path);
      res.end();
      return;
    }

    ::std::ostringstream oss(::std::ios::binary);
    ::std::string line;
    while (getline(is, line, '\n')) {
      double mz_value;
      mz_value = stod(line);

      auto data = Preprocess::mz_to_data(mz_value, reagentIon);
      if (data.model_data.get_rows() == 0) continue;
    
      auto preds = model->predict(data.model_data);

      if (data.model_data.get_rows() > 1)
	Postprocess::sort_preds(data.encoded_compounds, preds, data.ppms);
    
      auto decoded_compounds = Preprocess::decode_compounds(Preprocess::simplify_compounds(data.encoded_compounds));

      oss.setf(::std::ios::left, ::std::ios::adjustfield);
      oss << mz_value << ":" << ::std::endl;
      print_header(oss);
      for (size_t i = 0; i < decoded_compounds.size(); i++) {
	print_row(oss, { decoded_compounds[i].val, ::std::to_string(data.ppms.get(i, 0)), ::std::to_string(preds.get(i, 0)) });
      }
      page_break(oss);
    }

    is.close();
    if (!res.body.empty())
      res.body = "";
    res.write(oss.str());
    res.add_header("Content-Type", "text/plain");
    res.end();
  }

  void resolve_paths(const ::YAML::Node &config) {
    auto py_ex_path = config["paths"]["api"]["peak_fit_binary"].as<::std::string>();
    ::InferenceAPI::python_executable_path = (char *) malloc(sizeof(char) * (py_ex_path.size() + 1));
    ::std::atexit([] { free(python_executable_path); });
    strcpy(python_executable_path, py_ex_path.c_str());

    auto uploads_dir = config["paths"]["api"]["uploads_dir"].as<::std::string>();
    ::std::filesystem::create_directories(uploads_dir);
    
    ::InferenceAPI::graph_upload_dir = uploads_dir + config["paths"]["api"]["uploads_graphs_dir"].as<::std::string>();
    ::InferenceAPI::peak_output_dir = uploads_dir + config["paths"]["api"]["uploads_detected_peaks_dir"].as<::std::string>();

    ::std::filesystem::create_directories(::InferenceAPI::graph_upload_dir);
    ::std::filesystem::create_directories(::InferenceAPI::peak_output_dir);
  }
}
