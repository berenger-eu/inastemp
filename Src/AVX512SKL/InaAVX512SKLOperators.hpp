///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAAVX512SKLOPERATORS_HPP
#define INAAVX512SKLOPERATORS_HPP

#include "InastempConfig.h"

#ifndef INASTEMP_USE_AVX512SKL
#error InaAVX512SKLOperators is included but AVX512SKL is not enable in the configuration
#endif

#include <immintrin.h>
#include <cmath>

#ifdef INASTEMP_USE_AVX512SKL_OPERATORS

#include "AVX512COMMON/InaVecAVX512COMMONOperators.hpp"

#endif

#endif
