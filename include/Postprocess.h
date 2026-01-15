#ifndef __POSTPROCESS_H
#define __POSTPROCESS_H

#include <CNum.h>

namespace Postprocess {
  void sort_preds(::CNum::DataStructs::Matrix<double> &encoded_compounds,
		  ::CNum::DataStructs::Matrix<double> &preds,
		  ::CNum::DataStructs::Matrix<double> &ppms);
};

#endif
