///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAALTIVECOPERATORS_HPP
#define INAALTIVECOPERATORS_HPP

#include "InastempConfig.h"

#ifndef INASTEMP_USE_ALTIVEC
#error InaALTIVECOperators is included but ALTIVEC is not enable in the configuration
#endif

#include <altivec.h>
#undef bool
#undef vector
#undef pixel

#ifdef INASTEMP_USE_ALTIVEC_OPERATORS

#error Operators for ALTIVEC must exist, we cannot overload this type operators.

#endif

#endif
