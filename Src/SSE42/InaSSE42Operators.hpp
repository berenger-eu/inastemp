///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INASSE41OPERATORS_HPP
#define INASSE41OPERATORS_HPP

#include "InastempConfig.h"

#ifndef INASTEMP_USE_SSE41
#error InaSSE41Operators is included but SSE41 is not enable in the configuration
#endif


#ifdef INASTEMP_USE_SSE41_OPERATORS

#include "SSE3/InaSSE3Operators.hpp"

#endif

#endif
