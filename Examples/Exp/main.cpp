///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

// In this example we time the duration to compute a given number of Exponential
// We compare the scalar version and the BestType

#include "InastempConfig.h"
#include "SCALAR/InaVecSCALARDouble.hpp"
#include "Common/InaTimer.hpp"

#include <memory>
#include <iostream>


int main(int /*argc*/, char* /*argv*/ []) {
    using RealType = double;

    const size_t NbOverLoop = 5;
    const size_t NbExp      = 1000000;

    /////////////////////////////////////////////////////////////

    std::unique_ptr< RealType[] > resScalar(new RealType[NbExp]);
    {
        InaTimer timer;

        InaVecSCALAR<double> vectorizer;

        for (size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for (size_t idx = 0; idx < NbExp; ++idx) {
                resScalar[idx] = static_cast<double>(InaVecSCALAR<double>(double(idx % 200)).exp());
            }
        }

        timer.stop();
        std::cout << "Scalar for " << NbExp * NbOverLoop << " exp took " << timer.getElapsed() << "s\n";
    }

    /////////////////////////////////////////////////////////////

    // Note : we increase the length of the vector to avoid checking the loop size
    std::unique_ptr< RealType[] > resIna(new RealType[NbExp + InaVecBestType<double>::VecLength]);
    {
        InaTimer timer;

        alignas(64) double bufferX[InaVecBestType<double>::VecLength];


        for (size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for (size_t idx = 0; idx < NbExp; idx += InaVecBestType<double>::VecLength) {
                // Copy value into a buffer since we do it on the fly
                for (size_t idxX = 0; idxX < InaVecBestType<double>::VecLength; ++idxX) {
                    bufferX[idxX] = static_cast<double>((idx + idxX) % 200);
                }
                InaVecBestType<double>(bufferX).exp().storeInArray(&resIna[idx]);
            }
        }

        timer.stop();
        std::cout << "Vector " << InaVecBestType<double>::GetName() << " for " << NbExp * NbOverLoop << " exp took " << timer.getElapsed() << "s\n";
    }

    return 0;
}
