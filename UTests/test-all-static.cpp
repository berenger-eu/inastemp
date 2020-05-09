///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#include "InastempStaticConfig.h"

#include "core-test-all-static.hpp"
#include "FLOPS/InaVecFLOPS.hpp"

int main() {
    // clang-format off
    //TestAll<InaVec@TYPE@<float>> testerSingle;
    //return testerSingle.Run();
    // TODO
    TestAll<InaVec@TYPE@<double>> testerDouble;
    return testerDouble.Run();

    TestAll<InaVecFLOPS<InaVec@TYPE@<double>>> testerDoubleFlops;
    TestAll<InaVecFLOPS<InaVec@TYPE@<float>>> testerSingleFlops;
    // clang-format on
    return testerDouble.Run() + testerSingle.Run()
            + testerDoubleFlops.Run() + testerSingleFlops.Run();
}
