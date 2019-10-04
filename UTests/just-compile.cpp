///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

#include "InastempGlobal.h"


#include "@TYPE@/InaVec@TYPE@Double.hpp"
#include "@TYPE@/InaVec@TYPE@Float.hpp"


int main() {
    // clang-format off
    {
        InaVec@TYPE@<double> floatVec;
        InaVec@TYPE@<double> vecOne = 1;
        typename InaVec@TYPE@<double>::MaskType msk = (1 < vecOne);
        typename InaVec@TYPE@<double>::MaskType msk2 = msk;
        (void)msk2;
    }
    {
        InaVec@TYPE@<float> floatVec;
        InaVec@TYPE@<float> vecOne = 1;
        typename InaVec@TYPE@<float>::MaskType msk = (1 < vecOne);
        typename InaVec@TYPE@<float>::MaskType msk2 = msk;
        (void)msk2;
    }
    // clang-format on

    return 0;
}
