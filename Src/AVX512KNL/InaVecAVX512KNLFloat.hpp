///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVX512KNLFLOAT_HPP
#define INAVECAVX512KNLFLOAT_HPP

#include "InastempConfig.h"
#include "InaAVX512KNLOperators.hpp"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"
#include "AVX512COMMON/InaVecAVX512COMMONFloat.hpp"

#ifndef INASTEMP_USE_AVX512KNL
#error InaVecAVX512KNL512KNL<float> is included but AVX512KNL is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <immintrin.h>

#include <cmath>

// Forward declarations
template <class RealType>
using InaVecMaskAVX512KNL = InaVecMaskAVX512COMMON<RealType>;

template <class RealType>
using InaVecAVX512KNL = InaVecAVX512COMMON<RealType>;



#endif
