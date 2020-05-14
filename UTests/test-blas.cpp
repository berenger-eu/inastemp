///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#include "InastempGlobal.h"

#include "@TYPE@/InaVec@TYPE@Double.hpp"
#include "@TYPE@/InaVec@TYPE@Float.hpp"

#include "core-test-blas.hpp"
#include "FLOPS/InaVecFLOPS.hpp"

int main() {
    // clang-format off
    TestBlas<InaVec@TYPE@<double>> testerDouble;
    TestBlas<InaVec@TYPE@<float>> testerSingle;

    // clang-format on
    return testerDouble.Run() + testerSingle.Run();
}
