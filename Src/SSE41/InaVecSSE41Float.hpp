///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSSE41FLOAT_HPP
#define INAVECSSE41FLOAT_HPP

#include "InastempConfig.h"
#include "SSSE3/InaVecSSSE3Float.hpp"
#include "InaSSE41Operators.hpp"

#ifndef INASTEMP_USE_SSE41
#error InaVecSSE41<float> is included but SSE41 is not enable in the configuration
#endif

#include <tmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

template <class RealType>
class InaVecSSE41;

template <>
class alignas(16) InaVecSSE41<float> : public InaVecSSSE3<float> {
    using Parent = InaVecSSSE3<float>;

public:
    using InaVecSSSE3<float>::InaVecSSSE3;

    inline InaVecSSE41(){}

    inline InaVecSSE41(const InaVecSSSE3<float>& other)
        : Parent(other){}

    // Re-put exp to benefit from Floor
    inline InaVecSSE41<float> exp() const {
#ifdef __INTEL_COMPILER
        return _mm_exp_ps(Parent::vec);
#else
        static const __m128 COEFF_LOG2E = _mm_set1_ps(float(InaFastExp::CoeffLog2E()));
        static const __m128 COEFF_A     = _mm_set1_ps(float(InaFastExp::CoeffA32()));
        static const __m128 COEFF_B     = _mm_set1_ps(float(InaFastExp::CoeffB32()));
        static const __m128 COEFF_P5_A  = _mm_set1_ps(float(InaFastExp::GetCoefficient6()[5]));
        static const __m128 COEFF_P5_B  = _mm_set1_ps(float(InaFastExp::GetCoefficient6()[4]));
        static const __m128 COEFF_P5_C  = _mm_set1_ps(float(InaFastExp::GetCoefficient6()[3]));
        static const __m128 COEFF_P5_D  = _mm_set1_ps(float(InaFastExp::GetCoefficient6()[2]));
        static const __m128 COEFF_P5_E  = _mm_set1_ps(float(InaFastExp::GetCoefficient6()[1]));
        static const __m128 COEFF_P5_F  = _mm_set1_ps(float(InaFastExp::GetCoefficient6()[0]));

        __m128 x = Parent::vec * COEFF_LOG2E;

        const __m128 fractional_part = x - InaVecSSE41(x).floor().vec;

        __m128 factor = COEFF_P5_A;
        factor        = (factor * fractional_part + COEFF_P5_B);
        factor        = (factor * fractional_part + COEFF_P5_C);
        factor        = (factor * fractional_part + COEFF_P5_D);
        factor        = (factor * fractional_part + COEFF_P5_E);
        factor        = (factor * fractional_part + COEFF_P5_F);

        x -= factor;

        __m128i castedInteger = _mm_cvtps_epi32(COEFF_A * x + COEFF_B);

        return _mm_castsi128_ps(castedInteger);
#endif
    }

    inline InaVecSSE41<float> ExpLowAcc() const {
        static const __m128 COEFF_LOG2E = _mm_set1_ps(float(InaFastExp::CoeffLog2E()));
        static const __m128 COEFF_A     = _mm_set1_ps(float(InaFastExp::CoeffA32()));
        static const __m128 COEFF_B     = _mm_set1_ps(float(InaFastExp::CoeffB32()));
        static const __m128 COEFF_P5_D  = _mm_set1_ps(float(InaFastExp::GetCoefficient3()[2]));
        static const __m128 COEFF_P5_E  = _mm_set1_ps(float(InaFastExp::GetCoefficient3()[1]));
        static const __m128 COEFF_P5_F  = _mm_set1_ps(float(InaFastExp::GetCoefficient3()[0]));

        __m128 x = Parent::vec * COEFF_LOG2E;

        const __m128 fractional_part = x - InaVecSSE41(x).floor().vec;

        __m128 factor = COEFF_P5_D;
        factor        = (factor * fractional_part + COEFF_P5_E);
        factor        = (factor * fractional_part + COEFF_P5_F);

        x -= factor;

        __m128i castedInteger = _mm_cvtps_epi32(COEFF_A * x + COEFF_B);

        return _mm_castsi128_ps(castedInteger);
    }

    inline InaVecSSE41<float> floor() const {
        return _mm_floor_ps(Parent::vec);
    }

    inline static const char* GetName() {
        return "InaVecSSE41<float>";
    }

    inline static InaIfElse< InaVecSSE41<float> >::ThenClass If(const typename Parent::MaskType inTest) {
        return InaIfElse< InaVecSSE41<float> >::IfClass().If(inTest);
    }
};

#endif
