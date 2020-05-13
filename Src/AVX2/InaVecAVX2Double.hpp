///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVX2DOUBLE_HPP
#define INAVECAVX2DOUBLE_HPP

#include "InastempGlobal.h"
#include "AVX/InaVecAVXDouble.hpp"

#ifndef INASTEMP_USE_AVX2
#error InaVecAVX2<double> is included but AVX2 is not enable in the configuration
#endif

#include <tmmintrin.h>
#include <emmintrin.h>

template <class RealType>
class InaVecAVX2;

template <>
class alignas(32) InaVecAVX2<double> : public InaVecAVX<double> {
    using Parent = InaVecAVX<double>;

public:
    using Parent::GetVecLength;

    using InaVecAVX<double>::InaVecAVX;

    inline InaVecAVX2(){}

    inline InaVecAVX2(const InaVecAVX<double>& other)
        : Parent(other){}

    inline static const char* GetName() {
        return "InaVecAVX2<double>";
    }

#ifdef __FMA__
    static constexpr bool IsRealFma(){
        return true;
    }
#endif

    inline static InaIfElse< InaVecAVX2<double> >::ThenClass If(const typename Parent::MaskType& inTest) {
        return InaIfElse< InaVecAVX2<double> >::IfClass().If(inTest);
    }

#ifdef __FMA__
    inline static InaVecAVX2<double> Fma(const InaVecAVX2<double>& inValAdd, const InaVecAVX2<double>& inValMul1, const InaVecAVX2<double>& inValMul2){
        return _mm256_fmadd_pd(inValMul1.Parent::vec,inValMul2.Parent::vec, inValAdd.Parent::vec);
    }
#endif
};


#endif
