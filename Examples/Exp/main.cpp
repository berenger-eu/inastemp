///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

// In this example we time the duration to compute a given number of Exponential
// We compare the scalar version and the BestType

#include "InastempConfig.h"
#include "SCALAR/InaVecSCALARDouble.hpp"
#include "SCALAR/InaVecSCALARFloat.hpp"
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

        InaVecSCALAR<RealType> vectorizer;

        for (size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for (size_t idx = 0; idx < NbExp; ++idx) {
                resScalar[idx] = static_cast<RealType>(InaVecSCALAR<RealType>(RealType(idx % 200)).exp());
            }
        }

        timer.stop();
        std::cout << "Scalar for " << NbExp * NbOverLoop 
                  << " exp took " << timer.getElapsed() << "s (" << timer.getElapsed()/double(NbExp * NbOverLoop) << "s per exp)\n";
    }

    /////////////////////////////////////////////////////////////

    // Note : we increase the length of the vector to avoid checking the loop size
    std::unique_ptr< RealType[] > resIna(new RealType[NbExp + InaVecBestType<RealType>::VecLength]);
    {
        InaTimer timer;

        alignas(64) RealType bufferX[InaVecBestType<RealType>::VecLength];


        for (size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for (size_t idx = 0; idx < NbExp; idx += InaVecBestType<RealType>::VecLength) {
                // Copy value into a buffer since we do it on the fly
                for (size_t idxX = 0; idxX < InaVecBestType<RealType>::VecLength; ++idxX) {
                    bufferX[idxX] = static_cast<RealType>((idx + idxX) % 200);
                }
                InaVecBestType<RealType>(bufferX).exp().storeInArray(&resIna[idx]);
            }
        }

        timer.stop();
        std::cout << "Vector " << InaVecBestType<RealType>::GetName() << " for " << NbExp * NbOverLoop
                  << " exp took " << timer.getElapsed() << "s (" << timer.getElapsed()/double(NbExp * NbOverLoop) << "s per exp)\n";
    }

    /////////////////////////////////////////////////////////////

    // Note : we increase the length of the vector to avoid checking the loop size
    std::unique_ptr< RealType[] > resInaLowAcc(new RealType[NbExp + InaVecBestType<RealType>::VecLength]);
    {
        InaTimer timer;

        alignas(64) RealType bufferX[InaVecBestType<RealType>::VecLength];


        for (size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for (size_t idx = 0; idx < NbExp; idx += InaVecBestType<RealType>::VecLength) {
                // Copy value into a buffer since we do it on the fly
                for (size_t idxX = 0; idxX < InaVecBestType<RealType>::VecLength; ++idxX) {
                    bufferX[idxX] = static_cast<RealType>((idx + idxX) % 200);
                }
                InaVecBestType<RealType>(bufferX).expLowAcc().storeInArray(&resInaLowAcc[idx]);
            }
        }

        timer.stop();
        std::cout << "Vector low acc " << InaVecBestType<RealType>::GetName() << " for " << NbExp * NbOverLoop 
                  << " exp took " << timer.getElapsed() << "s (" << timer.getElapsed()/double(NbExp * NbOverLoop) << "s per exp)\n";
    }

    return 0;
}
