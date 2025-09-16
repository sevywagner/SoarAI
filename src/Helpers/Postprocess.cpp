#include "Postprocess.h"

using namespace CNum::DataStructs;

// ---- Sort the the prediction values and the compound they represent by the prediction values ----
void Postprocess::sort_preds(Matrix<double> &encoded_compounds, Matrix<double> &preds) {
  auto mask = preds.argsort();
  preds = preds[mask];
  encoded_compounds = encoded_compounds[mask];
}
