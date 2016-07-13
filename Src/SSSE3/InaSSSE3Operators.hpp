///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INASSSE3OPERATORS_HPP
#define INASSSE3OPERATORS_HPP

#include "InastempConfig.h"

#ifndef INASTEMP_USE_SSSE3
#error InaSSSE3Operators is included but SSSE3 is not enable in the configuration
#endif


#ifdef INASTEMP_USE_SSSE3_OPERATORS

#include "SSE3/InaSSE3Operators.hpp"

#endif

#endif
