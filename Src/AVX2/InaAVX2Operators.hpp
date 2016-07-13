///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAAVX2OPERATORS_HPP
#define INAAVX2OPERATORS_HPP

#include "InastempConfig.h"

#ifndef INASTEMP_USE_AVX2
#error InaAVX2Operators is included but AVX2 is not enable in the configuration
#endif


#ifdef INASTEMP_USE_AVX2_OPERATORS

#include "AVX/InaAVXOperators.hpp"

#endif

#endif
