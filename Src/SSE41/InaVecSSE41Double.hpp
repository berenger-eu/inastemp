///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSSE41DOUBLE_HPP
#define INAVECSSE41DOUBLE_HPP

#include "SSSE3/InaVecSSSE3Double.hpp"

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
    using Parent::GetVecLength;

    using InaVecSSSE3<double>::InaVecSSSE3;

    inline InaVecSSE41(){}

    inline InaVecSSE41(const InaVecSSSE3<double>& other)
        : Parent(other){}

    // Re-put exp to benefit from floor
    inline InaVecSSE41<double> exp() const {
#ifdef __INTEL_COMPILER
        return _mm_exp_pd(Parent::vec);
#else
        const __m128d COEFF_LOG2E = _mm_set1_pd(double(InaFastExp::CoeffLog2E()));
        const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
        const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
        const __m128d COEFF_P5_X  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_8()));
        const __m128d COEFF_P5_Y  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_7()));
        const __m128d COEFF_P5_Z  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_6()));
        const __m128d COEFF_P5_A  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_5()));
        const __m128d COEFF_P5_B  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_4()));
        const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_3()));
        const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_2()));
        const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_1()));
        const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_0()));

        __m128d x = _mm_mul_pd(Parent::vec, COEFF_LOG2E);

        const __m128d fractional_part = _mm_sub_pd(x, InaVecSSE41(x).floor().vec);

        __m128d factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(_mm_mul_pd(
                         COEFF_P5_X, fractional_part), COEFF_P5_Y), fractional_part),
                         COEFF_P5_Z),fractional_part), COEFF_P5_A), fractional_part),
                         COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part),
                         COEFF_P5_F);

        x = _mm_sub_pd(x,factor);

        x = _mm_add_pd(_mm_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) long int allvalint[GetVecLength()] = { _mm_cvtsd_si64(x),
                                                      _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

        return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
#endif
    }

    inline InaVecSSE41<double> expLowAcc() const {
        const __m128d COEFF_LOG2E = _mm_set1_pd(double(InaFastExp::CoeffLog2E()));
        const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
        const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
        const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_3()));
        const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_2()));
        const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_1()));
        const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_0()));

        __m128d x = _mm_mul_pd(Parent::vec, COEFF_LOG2E);

        const __m128d fractional_part = _mm_sub_pd(x, InaVecSSE41(x).floor().vec);

        __m128d factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd(_mm_mul_pd(
                                         COEFF_P5_C, fractional_part),
                                         COEFF_P5_D), fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = _mm_sub_pd(x,factor);

        x = _mm_add_pd(_mm_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) long int allvalint[GetVecLength()] = { _mm_cvtsd_si64(x),
                                                      _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

        return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
    }

    // Re-put exp to benefit from floor
    inline InaVecSSE41<double> exp10() const {
#ifdef __INTEL_COMPILER
        return _mm_exp10_pd(Parent::vec);
#else
        const __m128d COEFF_LOG210 = _mm_set1_pd(double(InaFastExp::CoeffLog210()));
        const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
        const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
        const __m128d COEFF_P5_X  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_8()));
        const __m128d COEFF_P5_Y  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_7()));
        const __m128d COEFF_P5_Z  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_6()));
        const __m128d COEFF_P5_A  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_5()));
        const __m128d COEFF_P5_B  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_4()));
        const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_3()));
        const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_2()));
        const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_1()));
        const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_0()));

        __m128d x = _mm_mul_pd(Parent::vec, COEFF_LOG210);

        const __m128d fractional_part = _mm_sub_pd(x, InaVecSSE41(x).floor().vec);

        __m128d factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(_mm_mul_pd(
                         COEFF_P5_X, fractional_part), COEFF_P5_Y), fractional_part),
                         COEFF_P5_Z),fractional_part), COEFF_P5_A), fractional_part),
                         COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part),
                         COEFF_P5_F);

        x = _mm_sub_pd(x,factor);

        x = _mm_add_pd(_mm_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) long int allvalint[GetVecLength()] = { _mm_cvtsd_si64(x),
                                                      _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

        return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
#endif
    }

    inline InaVecSSE41<double> exp10LowAcc() const {
        const __m128d COEFF_LOG210 = _mm_set1_pd(double(InaFastExp::CoeffLog210()));
        const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
        const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
        const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_3()));
        const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_2()));
        const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_1()));
        const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_0()));

        __m128d x = _mm_mul_pd(Parent::vec, COEFF_LOG210);

        const __m128d fractional_part = _mm_sub_pd(x, InaVecSSE41(x).floor().vec);

        __m128d factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd(_mm_mul_pd(
                                         COEFF_P5_C, fractional_part),
                                         COEFF_P5_D), fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = _mm_sub_pd(x,factor);

        x = _mm_add_pd(_mm_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) long int allvalint[GetVecLength()] = { _mm_cvtsd_si64(x),
                                                      _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

        return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
    }


    inline InaVecSSE41<double> exp2() const {
#ifdef __INTEL_COMPILER
        return _mm_exp2_pd(Parent::vec);
#else
        const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
        const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
        const __m128d COEFF_P5_X  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_8()));
        const __m128d COEFF_P5_Y  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_7()));
        const __m128d COEFF_P5_Z  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_6()));
        const __m128d COEFF_P5_A  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_5()));
        const __m128d COEFF_P5_B  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_4()));
        const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_3()));
        const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_2()));
        const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_1()));
        const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_0()));

        __m128d x = Parent::vec;

        const __m128d fractional_part = _mm_sub_pd(x, InaVecSSE41(x).floor().vec);

        __m128d factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(_mm_mul_pd(
                         COEFF_P5_X, fractional_part), COEFF_P5_Y), fractional_part),
                         COEFF_P5_Z),fractional_part), COEFF_P5_A), fractional_part),
                         COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part),
                         COEFF_P5_F);

        x = _mm_sub_pd(x,factor);

        x = _mm_add_pd(_mm_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) long int allvalint[GetVecLength()] = { _mm_cvtsd_si64(x),
                                                      _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

        return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
#endif
    }

    inline InaVecSSE41<double> exp2LowAcc() const {
        const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
        const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
        const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_3()));
        const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_2()));
        const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_1()));
        const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_0()));

        __m128d x = Parent::vec;

        const __m128d fractional_part = _mm_sub_pd(x, InaVecSSE41(x).floor().vec);

        __m128d factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                         _mm_mul_pd(_mm_add_pd(_mm_mul_pd(
                                         COEFF_P5_C, fractional_part),
                                         COEFF_P5_D), fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = _mm_sub_pd(x,factor);

        x = _mm_add_pd(_mm_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) long int allvalint[GetVecLength()] = { _mm_cvtsd_si64(x),
                                                      _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

        return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
    }

    inline InaVecSSE41<double> log() const{
#ifdef __INTEL_COMPILER
        return _mm_log_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.log(*this);
#endif
    }

    inline InaVecSSE41<double> log2() const{
#ifdef __INTEL_COMPILER
        return _mm_log2_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.log2(*this);
#endif
    }

    inline InaVecSSE41<double> log10() const{
#ifdef __INTEL_COMPILER
        return _mm_log10_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.log10(*this);
#endif
    }

    inline InaVecSSE41<double> sin() const{
#ifdef __INTEL_COMPILER
        return _mm_sin_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.sin(*this);
#endif
    }
    inline InaVecSSE41<double> cos() const{
#ifdef __INTEL_COMPILER
        return _mm_cos_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.cos(*this);
#endif
    }
    inline InaVecSSE41<double> tan() const{
#ifdef __INTEL_COMPILER
        return _mm_tan_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.tan(*this);
#endif
    }
    inline InaVecSSE41<double> asin() const{
#ifdef __INTEL_COMPILER
        return _mm_asin_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.asin(*this);
#endif
    }
    inline InaVecSSE41<double> acos() const{
#ifdef __INTEL_COMPILER
        return _mm_acos_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.acos(*this);
#endif
    }
    inline InaVecSSE41<double> atan() const{
#ifdef __INTEL_COMPILER
        return _mm_atan_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.atan(*this);
#endif
    }
    
    inline InaVecSSE41<double> sinh() const{
#ifdef __INTEL_COMPILER
        return _mm_sinh_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.sinh(*this);
#endif
    }
    inline InaVecSSE41<double> cosh() const{
#ifdef __INTEL_COMPILER
        return _mm_cosh_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.cosh(*this);
#endif
    }
    inline InaVecSSE41<double> tanh() const{
#ifdef __INTEL_COMPILER
        return _mm_tanh_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.tanh(*this);
#endif
    }
    inline InaVecSSE41<double> asinh() const{
#ifdef __INTEL_COMPILER
        return _mm_asinh_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.asinh(*this);
#endif
    }
    inline InaVecSSE41<double> acosh() const{
#ifdef __INTEL_COMPILER
        return _mm_acosh_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.acosh(*this);
#endif
    }
    inline InaVecSSE41<double> atanh() const{
#ifdef __INTEL_COMPILER
        return _mm_atanh_pd(Parent::vec);
#else
        InaMath<InaVecSSE41<double>> a;
        return a.atanh(*this);
#endif
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
