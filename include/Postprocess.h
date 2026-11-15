#ifndef __POSTPROCESS_H
#define __POSTPROCESS_H

#include <CNum.h>
#include <cstdint>

namespace Postprocess {
  void sort_preds(::CNum::DataStructs::Matrix<double> &encoded_compounds,
		  ::CNum::DataStructs::Matrix<double> &preds,
		  ::CNum::DataStructs::Matrix<double> &ppms);

  void sort_preds(::CNum::DataStructs::Matrix<double> &encoded_compounds,
		  ::CNum::DataStructs::Matrix<double> &preds,
		  ::CNum::DataStructs::Matrix<double> &ppms,
		  ::CNum::DataStructs::Matrix<uint8_t> &criterea_encodings);
};

#endif
