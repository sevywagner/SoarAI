#ifndef __SYS_UTILS_H
#define __SYS_UTILS_H

#include <string>
#include <stdexcept>
#include "InferenceAPI.h"

void execute_peak_fit_bin(::std::string mz_b,
			  ::std::string mz_av,
			  ::std::string output_path);
  
#endif
