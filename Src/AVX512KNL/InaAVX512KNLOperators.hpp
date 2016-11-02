///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAAVX512KNLOPERATORS_HPP
#define INAAVX512KNLOPERATORS_HPP

#include "InastempConfig.h"

#ifndef INASTEMP_USE_AVX512KNL
#error InaAVX512KNLOperators is included but AVX512KNL is not enable in the configuration
#endif

#include <immintrin.h>
#include <cmath>

#ifdef INASTEMP_USE_AVX512KNL_OPERATORS

#include "AVX512COMMON/InaVecAVX512COMMONOperators.hpp"

#endif

#endif
