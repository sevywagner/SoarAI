#include <CNum.h>
#include "Preprocess.h"

int main(int argc, char *argv[]) {
  if (argc != 4) {
    throw ::std::invalid_argument("Must pass 3 paths as arguments i.e. ./this_program /path/to/peak/list.txt /path/to/combined/combo/output.csv /path/to/train/test/split");
  }
  
  auto peak_list_data = Preprocess::PrepareDataset::parse_peak_list(argv[1]);
  auto encoded_unsimplified = Preprocess::encode_compounds(peak_list_data.compound_strings);
  auto encoded_simplified = Preprocess::simplify_compounds(encoded_unsimplified);

  Preprocess::PrepareDataset::create_combo_file(argv[2],
  						peak_list_data,
  						encoded_unsimplified,
  						encoded_simplified);
  
  Preprocess::PrepareDataset::train_test_split(argv[2], argv[3]);
}
