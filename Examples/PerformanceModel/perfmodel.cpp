///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#include "InastempConfig.h"
#include "SCALAR/InaVecSCALARDouble.hpp"
#include "SCALAR/InaVecSCALARFloat.hpp"
#include "Common/InaTimer.hpp"

#include <cassert>
#include <iostream>
#include <memory>
#include <fstream>

#ifdef INASTEMP_USE_SSE3
#include "SSE3/InaVecSSE3Double.hpp"
#include "SSE3/InaVecSSE3Float.hpp"
#endif

#ifdef INASTEMP_USE_AVX
#include "AVX/InaVecAVXDouble.hpp"
#include "AVX/InaVecAVXFloat.hpp"
#include <immintrin.h>
#endif

#ifdef INASTEMP_USE_AVX512KNL
#include "AVX512KNL/InaVecAVX512KNLDouble.hpp"
#include "AVX512KNL/InaVecAVX512KNLFloat.hpp"
#endif


size_t cpt = 0;

template < class VecType >
VecType DummyPower(VecType inVal, const int inPow){
    VecType res = 1;
    for( int idx = 0 ; idx < inPow ; ++idx) {
        res *= inVal;
        cpt++;
    }
    return res;
}

template < class VecType, class MaskType,  class RealType >
double KernelAll(const RealType inVal, const MaskType msk, const int pow1, const int pow2,
                 const int pow3, const size_t nbLoops, VecType& res){
    InaTimer timer;

    VecType vec = inVal;
    for(size_t idxLoop = 0 ; idxLoop < nbLoops ; ++idxLoop){
        vec += DummyPower(inVal, pow1);
        vec += VecType::IfElse(msk, DummyPower(inVal, pow2), DummyPower(inVal, pow3));
    }
    res += vec;

    timer.stop();
    return timer.getElapsed();
}

template < class VecType, class MaskType, class RealType >
double KernelIfTrue(const RealType inVal, const MaskType msk, const int pow1, const int pow2,
                    const int pow3, const size_t nbLoops, VecType& res){
    InaTimer timer;

    VecType vec = inVal;
    for(size_t idxLoop = 0 ; idxLoop < nbLoops ; ++idxLoop){
        vec += DummyPower(inVal, pow1);
        if( msk.isAllTrue() ){
            vec += DummyPower(inVal, pow2);
        }
        else{
            vec += VecType::IfElse(msk, DummyPower(inVal, pow2), DummyPower(inVal, pow3));
        }
    }

    res += vec;

    timer.stop();
    return timer.getElapsed();
}

template < class VecType, class RealType >
int test(const size_t nbLoops, std::ofstream& myfile){
    std::cout << "Test " << VecType::GetName() << " in double" << std::endl;

    VecType res = 0;

    for(int idxSizeTrue = 0 ; idxSizeTrue <= VecType::VecLength ; ++idxSizeTrue){
        RealType mskValues[VecType::VecLength];
        for(int idx = 0 ; idx < VecType::VecLength ; ++idx){
            mskValues[idx] = (idx < idxSizeTrue ? 0. : 1.);
        }
        typename VecType::MaskType msk = VecType(mskValues).isZeroMask();

        for(int power2 = 1 ; power2 < (1 << 8 ) ; power2 *= 2){
            for(int power1p = 0 ; power1p <= 100 ; power1p += 25){
                const int power1 = (power2*power1p)/100;

                for(int power3p = 0 ; power3p <= 100 ; power3p += 25){
                    const int power3 = (power2*power3p)/100;

                    {
                        const double duration = KernelAll<VecType, typename VecType::MaskType, RealType>(
                                    RealType(1), msk, power1, power2, power3, nbLoops, res);

                        const size_t Flops = nbLoops * VecType::VecLength * (2 + power1 + power2 + power3);
                        const size_t EffFlops = nbLoops * (VecType::VecLength * (2 + power1) + power2*idxSizeTrue + power3*(VecType::VecLength-idxSizeTrue));

                        const double gflops = (double(Flops)/duration)/1E9;
                        const double effgflops = (double(EffFlops)/duration)/1E9;

                        std::cout << "[ALL] power 1 " << power1 << " power 2 " << power2 << " power 3 " << power3
                                  << " duration " << duration
                                  << " GFlops " << gflops
                                  << " Effective-GFlops " << effgflops << std::endl;

                        myfile << "\"all\","<< power1<<","<<power2<<","<<power3<<","<<duration<<","<<gflops<<","<<effgflops<<"\n";
                    }
                    {
                        const double duration = KernelIfTrue<VecType, typename VecType::MaskType, RealType>(
                                    RealType(1), msk, power1, power2, power3, nbLoops, res);

                        const size_t Flops = (msk.isAllTrue()?
                                                nbLoops * VecType::VecLength * (2 + power1 + power2) :
                                                nbLoops * VecType::VecLength * (2 + power1 + power2 + power3));
                        const size_t EffFlops = nbLoops * (VecType::VecLength * (2 + power1) + power2*idxSizeTrue + power3*(VecType::VecLength-idxSizeTrue));

                        const double gflops = (double(Flops)/duration)/1E9;
                        const double effgflops = (double(EffFlops)/duration)/1E9;

                        std::cout << "[IFT] power 1 " << power1 << " power 2 " << power2 << " power 3 " << power3
                                  << " duration " << duration
                                  << " GFlops " << gflops
                                  << " Effective-GFlops " << effgflops << std::endl;

                        myfile << "\"iftrue\","<< power1<<","<<power2<<","<<power3<<","<<duration<<","<<gflops<<","<<effgflops<<"\n";
                    }
                }
            }
        }
    }

    return int(res.horizontalSum());
}

int main() {
    std::ofstream myfile;
    myfile.open ("res.csv");

    myfile << "mode,power1,power2,power3,duration,gflops,effgflops\n";

    volatile size_t nbLoops = 100000000;

    int res = 0;

    res += test<InaVecSCALAR<double>, double>(nbLoops, myfile);

#ifdef INASTEMP_USE_SSE3
    res += test<InaVecSSE3<double>, double>(nbLoops, myfile);
#endif

#ifdef INASTEMP_USE_AVX
    res += test<InaVecAVX<double>, double>(nbLoops, myfile);
#endif

#ifdef INASTEMP_USE_AVX512KNL
    res += test<InaVecAVX512KNL<double>, double>(nbLoops, myfile);
#endif

    return res + int(cpt);
}
