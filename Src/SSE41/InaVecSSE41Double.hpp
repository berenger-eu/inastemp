///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSSE41DOUBLE_HPP
#define INAVECSSE41DOUBLE_HPP

#include "SSSE3/InaVecSSSE3Double.hpp"
#include "InaSSE41Operators.hpp"

#ifndef INASTEMP_USE_SSE41
#error InaVecSSE41<double> is included but SSE41 is not enable in the configuration
#endif

#include <tmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

template <class RealType>
class InaVecSSE41;

template <>
class alignas(16) InaVecSSE41<double> : public InaVecSSSE3<double> {
    using Parent = InaVecSSSE3<double>;

public:
    using InaVecSSSE3<double>::InaVecSSSE3;

    inline InaVecSSE41(){}

    inline InaVecSSE41(const InaVecSSSE3<double>& other)
        : Parent(other){}

    // Re-put exp to benefit from floor
    inline InaVecSSE41<double> exp() const {
#ifdef __INTEL_COMPILER
        return _mm_exp_pd(Parent::vec);
#else
        static const __m128d COEFF_LOG2E = _mm_set1_pd(double(InaFastExp::CoeffLog2E()));
        static const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
        static const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
        static const __m128d COEFF_P5_X  = _mm_set1_pd(double(InaFastExp::GetCoefficient9()[8]));
        static const __m128d COEFF_P5_Y  = _mm_set1_pd(double(InaFastExp::GetCoefficient9()[7]));
        static const __m128d COEFF_P5_Z  = _mm_set1_pd(double(InaFastExp::GetCoefficient9()[6]));
        static const __m128d COEFF_P5_A  = _mm_set1_pd(double(InaFastExp::GetCoefficient9()[5]));
        static const __m128d COEFF_P5_B  = _mm_set1_pd(double(InaFastExp::GetCoefficient9()[4]));
        static const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient9()[3]));
        static const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient9()[2]));
        static const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient9()[1]));
        static const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient9()[0]));

        __m128d x = Parent::vec * COEFF_LOG2E;

        const __m128d fractional_part = x - InaVecSSE41(x).floor().vec;

        __m128d factor = COEFF_P5_X;
        factor         = (factor * fractional_part + COEFF_P5_Y);
        factor         = (factor * fractional_part + COEFF_P5_Z);
        factor         = (factor * fractional_part + COEFF_P5_A);
        factor         = (factor * fractional_part + COEFF_P5_B);
        factor         = (factor * fractional_part + COEFF_P5_C);
        factor         = (factor * fractional_part + COEFF_P5_D);
        factor         = (factor * fractional_part + COEFF_P5_E);
        factor         = (factor * fractional_part + COEFF_P5_F);

        x -= factor;

        x = (COEFF_A * x + COEFF_B);

        alignas(64) long int allvalint[VecLength] = { _mm_cvtsd_si64(x),
                                                      _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

        return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
#endif
    }

    inline InaVecSSE41<double> expLowAcc() const {
        static const __m128d COEFF_LOG2E = _mm_set1_pd(double(InaFastExp::CoeffLog2E()));
        static const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
        static const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
        static const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient4()[3]));
        static const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient4()[2]));
        static const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient4()[1]));
        static const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient4()[0]));

        __m128d x = Parent::vec * COEFF_LOG2E;

        const __m128d fractional_part = x - InaVecSSE41(x).floor().vec;

        __m128d factor = COEFF_P5_C;
        factor         = (factor * fractional_part + COEFF_P5_D);
        factor         = (factor * fractional_part + COEFF_P5_E);
        factor         = (factor * fractional_part + COEFF_P5_F);

        x -= factor;

        x = (COEFF_A * x + COEFF_B);

        alignas(64) long int allvalint[VecLength] = { _mm_cvtsd_si64(x),
                                                      _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

        return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
    }

    inline InaVecSSE41<double> floor() const {
        return _mm_floor_pd(Parent::vec);
    }

    inline static const char* GetName(){
        return "InaVecSSE41<double>";
    }

    inline static InaIfElse< InaVecSSE41<double> >::ThenClass If(const typename Parent::MaskType& inTest) {
        return InaIfElse< InaVecSSE41<double> >::IfClass().If(inTest);
    }
};

#endif
