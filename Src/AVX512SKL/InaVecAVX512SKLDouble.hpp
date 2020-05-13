///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVX512SKLDOUBLE_HPP
#define INAVECAVX512SKLDOUBLE_HPP

#include "InastempGlobal.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"
#include "AVX512COMMON/InaVecAVX512COMMONDouble.hpp"

#ifndef INASTEMP_USE_AVX512SKL
#error InaVecAVX512SKL<double> is included but AVX512SKL is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <immintrin.h>

#include <cmath>

// Forward declarations
template <class RealType>
using InaVecMaskAVX512SKL = InaVecMaskAVX512COMMON<RealType>;

// Forward declarations
template <class RealType>
class InaVecAVX512SKL;

template <>
class alignas(InaVecAVX512COMMON<double>::Alignement) InaVecAVX512SKL<double> : public InaVecAVX512COMMON<double> {
    using Parent = InaVecAVX512COMMON<double>;

public:
    using Parent::GetVecLength;

    using InaVecAVX512COMMON<double>::InaVecAVX512COMMON;

    inline InaVecAVX512SKL(){}

    inline InaVecAVX512SKL(const InaVecAVX512COMMON<double>& other)
        : Parent(other){}

    inline static const char* GetName(){
        return "InaVecAVX512SKL<double>";
    }

    inline static InaIfElse< InaVecAVX512SKL<double> >::ThenClass If(const typename Parent::MaskType& inTest) {
        return InaIfElse< InaVecAVX512SKL<double> >::IfClass().If(inTest);
    }

    inline InaVecAVX512SKL exp() const {
         const __m512d COEFF_LOG2E = _mm512_set1_pd(double(InaFastExp::CoeffLog2E()));
         const __m512d COEFF_A     = _mm512_set1_pd(double(InaFastExp::CoeffA64()));
         const __m512d COEFF_B     = _mm512_set1_pd(double(InaFastExp::CoeffB64()));
         const __m512d COEFF_P5_X  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_8()));
         const __m512d COEFF_P5_Y  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_7()));
         const __m512d COEFF_P5_Z  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_6()));
         const __m512d COEFF_P5_A  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_5()));
         const __m512d COEFF_P5_B  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_4()));
         const __m512d COEFF_P5_C  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_3()));
         const __m512d COEFF_P5_D  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_2()));
         const __m512d COEFF_P5_E  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_1()));
         const __m512d COEFF_P5_F  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_0()));

         __m512d x = _mm512_mul_pd(Parent::vec, COEFF_LOG2E);

         const __m512d fractional_part = _mm512_sub_pd(x, InaVecAVX512SKL(x).floor().vec);

         __m512d factor = _mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(
                          _mm512_mul_pd(_mm512_add_pd( _mm512_mul_pd(_mm512_add_pd(
                          _mm512_mul_pd(_mm512_add_pd( _mm512_mul_pd(_mm512_add_pd(
                          _mm512_mul_pd(_mm512_add_pd( _mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(
                          COEFF_P5_X, fractional_part), COEFF_P5_Y), fractional_part),
                          COEFF_P5_Z),fractional_part), COEFF_P5_A), fractional_part),
                          COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                          COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part),
                          COEFF_P5_F);

         x = _mm512_sub_pd(x,factor);

         x = _mm512_add_pd(_mm512_mul_pd(COEFF_A, x), COEFF_B);

         alignas(64) double allvalreal[GetVecLength()];
         _mm512_store_pd(allvalreal, x);

         alignas(64) long long int allvalint[GetVecLength()] = { static_cast<long long int>(allvalreal[0]), static_cast<long long int>(allvalreal[1]),
                                                            static_cast<long long int>(allvalreal[2]), static_cast<long long int>(allvalreal[3]),
                                                            static_cast<long long int>(allvalreal[4]), static_cast<long long int>(allvalreal[5]),
                                                            static_cast<long long int>(allvalreal[6]), static_cast<long long int>(allvalreal[7]) };

         return _mm512_castsi512_pd(_mm512_load_epi64(reinterpret_cast<const __m512i*>(allvalint)));
    }

    inline InaVecAVX512SKL expLowAcc() const {
        const __m512d COEFF_LOG2E = _mm512_set1_pd(double(InaFastExp::CoeffLog2E()));
        const __m512d COEFF_A     = _mm512_set1_pd(double(InaFastExp::CoeffA64()));
        const __m512d COEFF_B     = _mm512_set1_pd(double(InaFastExp::CoeffB64()));
        const __m512d COEFF_P5_C  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_3()));
        const __m512d COEFF_P5_D  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_2()));
        const __m512d COEFF_P5_E  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_1()));
        const __m512d COEFF_P5_F  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_0()));

        __m512d x = _mm512_mul_pd(Parent::vec, COEFF_LOG2E);

        const __m512d fractional_part = _mm512_sub_pd(x, InaVecAVX512SKL(x).floor().vec);

        __m512d factor = _mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(
                         _mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(
                                         COEFF_P5_C, fractional_part),
                                         COEFF_P5_D), fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = _mm512_sub_pd(x,factor);

        x = _mm512_add_pd(_mm512_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) double allvalreal[GetVecLength()];
        _mm512_store_pd(allvalreal, x);

        alignas(64) long long int allvalint[GetVecLength()] = { static_cast<long long int>(allvalreal[0]), static_cast<long long int>(allvalreal[1]),
                                                           static_cast<long long int>(allvalreal[2]), static_cast<long long int>(allvalreal[3]),
                                                           static_cast<long long int>(allvalreal[4]), static_cast<long long int>(allvalreal[5]),
                                                           static_cast<long long int>(allvalreal[6]), static_cast<long long int>(allvalreal[7]) };

        return _mm512_castsi512_pd(_mm512_load_epi64(reinterpret_cast<const __m512i*>(allvalint)));
    }

    inline InaVecAVX512SKL floor() const {
        const __m512i vecConvLongInt = _mm512_cvt_roundpd_epi64(Parent::vec, (/*_MM_FROUND_TO_NEG_INF*/_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
        const __m512i valuesDec = _mm512_sub_epi64(vecConvLongInt, _mm512_set1_epi64(1));

        const __mmask8 maskPositive = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_LE_OQ);
        const __mmask8 maskNegative = __mmask8(~maskPositive);

        const __m512i valuesToConv = _mm512_or_epi64(_mm512_maskz_mov_epi64(maskPositive, vecConvLongInt),
                                                   _mm512_maskz_mov_epi64(maskNegative, valuesDec));

        return _mm512_cvtepi64_pd(valuesToConv);
    }
};

#endif
