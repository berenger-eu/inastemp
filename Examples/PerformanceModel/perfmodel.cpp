///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#include "InastempConfig.h"
#include "SCALAR/InaVecSCALARDouble.hpp"
#include "SCALAR/InaVecSCALARFloat.hpp"

#include <cassert>
#include <iostream>
#include <memory>

template < class VecType >
VecType DummyPower(VecType inVal, const int inPow){
    VecType res = 1;
    for( int idx = 0 ; idx < inPow ; ++idx) {
        res *= inVal;
    }
    return res;
}

template < class VecType, class MaskType,  class RealType >
void KernelAll(RealType inVal, MaskType msk, int pow1, int pow2, int pow3, size_t nbLoops){
    VecType vec = inVal;
    for(size_t idxLoop = 0 ; idxLoop < nbLoops ; ++idxLoop){
        vec += DummyPower(inVal, pow1);
        vec += VecType::IfElse(msk, DummyPower(inVal, pow2), DummyPower(inVal, pow3));
    }
}

template < class VecType, class MaskType, class RealType >
void KernelIfTrue(RealType inVal, MaskType msk, int pow1, int pow2, int pow3, size_t nbLoops){
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
}

template < class VecType, class RealType >
void test(const size_t nbLoops){
    std::cout << "The best vectorizer computes " << VecType::GetName() << " double values together." << std::endl;
    KernelAll<VecType, typename VecType::MaskType, RealType>(1, typename VecType::MaskType(), 1, 2, 3, nbLoops);
}

int main() {

    const size_t nbLoops = 10000;

    test<InaVecBestType<double>, double>(nbLoops);

    return 0;
}
