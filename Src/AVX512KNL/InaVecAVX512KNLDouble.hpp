///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVX512KNLDOUBLE_HPP
#define INAVECAVX512KNLDOUBLE_HPP

#include "InastempConfig.h"
#include "InaAVX512KNLOperators.hpp"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_AVX512KNL
#error InaVecAVX512KNL<double> is included but AVX512KNL is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <immintrin.h>

#include <cmath>

// Forward declarations
template <class RealType>
class InaVecMaskAVX512KNL;

template <class RealType>
class InaVecAVX512KNL;

// Mask type
template <>
class alignas(64) InaVecMaskAVX512KNL<double> {
    __m512i mask;
public:
    // Classic constructors
    inline InaVecMaskAVX512KNL(){}

    inline InaVecMaskAVX512KNL(const InaVecMaskAVX512KNL&) = default;
    inline InaVecMaskAVX512KNL& operator=(const InaVecMaskAVX512KNL&) = default;

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskAVX512KNL(const __m512i inMask)
        : mask(inMask){}

    inline InaVecMaskAVX512KNL& operator=(const __m512i inMask){
        mask = inMask;
        return (*this);
    }

    inline explicit operator __m512i() const{
        return mask;
    }

    inline __m512i getMask() const{
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskAVX512KNL(const bool inBool){
        mask = (inBool? _mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)) : _mm512_setzero_si512());
    }

    inline InaVecMaskAVX512KNL& operator=(const bool inBool){
        mask = (inBool? _mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)) : _mm512_setzero_si512());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskAVX512KNL Not() const{
        return NotAnd(mask, _mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)));
    }

    inline bool isAllTrue() const{
        // true if all FF => !FF => 0 & FF => 0
        const __mmask8 testResult = _mm512_cmp_epu64_mask(mask, _mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)), _MM_CMPINT_EQ);
        return testResult == 0xFF;
    }

    inline bool isAllFalse() const{
        // true if all zero
        const __mmask8 testResult = _mm512_cmp_epu64_mask(mask, _mm512_setzero_si512(), _MM_CMPINT_EQ);
        return testResult == 0xFF;
    }

    // Double args methods
    inline static InaVecMaskAVX512KNL And(const InaVecMaskAVX512KNL& inMask1, const InaVecMaskAVX512KNL& inMask2){
        return InaVecMaskAVX512KNL(_mm512_and_si512(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskAVX512KNL NotAnd(const InaVecMaskAVX512KNL& inMask1, const InaVecMaskAVX512KNL& inMask2){
        return InaVecMaskAVX512KNL(_mm512_andnot_si512(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskAVX512KNL Or(const InaVecMaskAVX512KNL& inMask1, const InaVecMaskAVX512KNL& inMask2){
        return InaVecMaskAVX512KNL(_mm512_or_si512(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskAVX512KNL Xor(const InaVecMaskAVX512KNL& inMask1, const InaVecMaskAVX512KNL& inMask2){
        return InaVecMaskAVX512KNL(_mm512_xor_si512(inMask1.mask, inMask2.mask));
    }

    inline static bool IsEqual(const InaVecMaskAVX512KNL& inMask1, const InaVecMaskAVX512KNL& inMask2){
        const __mmask8 testResult = _mm512_cmp_epu64_mask(inMask1.mask, inMask2.mask, _MM_CMPINT_EQ);
        return testResult == 0xFF;
    }

    inline static bool IsNotEqual(const InaVecMaskAVX512KNL& inMask1, const InaVecMaskAVX512KNL& inMask2){
        const __mmask8 testResult = _mm512_cmp_epu64_mask(inMask1.mask, inMask2.mask, _MM_CMPINT_EQ);
        return testResult != 0xFF;
    }
};

// Mask must have operators
inline InaVecMaskAVX512KNL<double> operator&(const InaVecMaskAVX512KNL<double>& inMask1, const InaVecMaskAVX512KNL<double>& inMask2){
    return InaVecMaskAVX512KNL<double>::And(inMask1, inMask2);
}

inline InaVecMaskAVX512KNL<double> operator|(const InaVecMaskAVX512KNL<double>& inMask1, const InaVecMaskAVX512KNL<double>& inMask2){
    return InaVecMaskAVX512KNL<double>::Or(inMask1, inMask2);
}

inline InaVecMaskAVX512KNL<double> operator^(const InaVecMaskAVX512KNL<double>& inMask1, const InaVecMaskAVX512KNL<double>& inMask2){
    return InaVecMaskAVX512KNL<double>::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskAVX512KNL<double>& inMask1, const InaVecMaskAVX512KNL<double>& inMask2){
    return InaVecMaskAVX512KNL<double>::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskAVX512KNL<double>& inMask1, const InaVecMaskAVX512KNL<double>& inMask2){
    return InaVecMaskAVX512KNL<double>::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(64) InaVecAVX512KNL<double> {
protected:
    __m512d vec;

public:
    using VecRawType           = __m512d;
    using MaskType             = InaVecMaskAVX512KNL<double>;
    using RealType             = double;
    static const int VecLength = 8;
    static const int Alignement= 64;

    inline InaVecAVX512KNL(){}
    inline InaVecAVX512KNL(const InaVecAVX512KNL&) = default;
    inline InaVecAVX512KNL& operator = (const InaVecAVX512KNL&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecAVX512KNL(const __m512d inVec)
        : vec(inVec){
    }

    inline InaVecAVX512KNL& operator=(const __m512d inVec){
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __m512d inVec){
        vec = inVec;
    }

    inline explicit operator __m512d() const{
        return vec;
    }

    inline __m512d getVec() const{
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecAVX512KNL(const double val)
        : vec(_mm512_set1_pd(val)){
    }

    inline InaVecAVX512KNL& operator=(const double val){
        vec = _mm512_set1_pd(val);
        return *this;
    }

    inline void setFromScalar(const double val){
        vec = _mm512_set1_pd(val);
    }

    // Constructor from vec
    inline explicit InaVecAVX512KNL(const double ptr[])
        : vec(_mm512_loadu_pd(ptr)){
    }

    inline InaVecAVX512KNL& setFromArray(const double ptr[]){
        vec = _mm512_loadu_pd(ptr);
        return *this;
    }

    inline InaVecAVX512KNL& setFromAlignedArray(const double ptr[]){
        vec = _mm512_load_pd(ptr);
        return *this;
    }

    inline InaVecAVX512KNL& setFromIndirectArray(const double values[], const int inIndirection[]) {
        vec = _mm512_set_pd(
                    values[inIndirection[7]],
                    values[inIndirection[6]],
                    values[inIndirection[5]],
                    values[inIndirection[4]],
                    values[inIndirection[3]],
                    values[inIndirection[2]],
                    values[inIndirection[1]],
                    values[inIndirection[0]]);
        return *this;
    }

    inline InaVecAVX512KNL& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        vec = _mm512_set_pd(
                    inArray[inIndirection1[7] * inLeadingDimension + inIndirection2[7]],
                    inArray[inIndirection1[6] * inLeadingDimension + inIndirection2[6]],
                    inArray[inIndirection1[5] * inLeadingDimension + inIndirection2[5]],
                    inArray[inIndirection1[4] * inLeadingDimension + inIndirection2[4]],
                    inArray[inIndirection1[3] * inLeadingDimension + inIndirection2[3]],
                    inArray[inIndirection1[2] * inLeadingDimension + inIndirection2[2]],
                    inArray[inIndirection1[1] * inLeadingDimension + inIndirection2[1]],
                    inArray[inIndirection1[0] * inLeadingDimension + inIndirection2[0]]);
        return *this;
    }

    // Move back to array
    inline void storeInArray(double ptr[]) const {
        _mm512_storeu_pd(ptr, vec);
    }

    inline void storeInAlignedArray(double ptr[]) const {
        _mm512_store_pd(ptr, vec);
    }

    // Acce to individual values
    inline double at(const int index) const {
        alignas(Alignement) double allval[VecLength];
        _mm512_store_pd(allval, vec);
        return allval[index];
    }

    // Horizontal operation
    inline double horizontalSum() const {
#ifdef __INTEL_COMPILER
        return _mm512_reduce_add_pd(vec);
#else
        __m256d low  = _mm512_castpd512_pd256(vec);
        __m256d high = _mm512_extractf64x4_pd(vec, 1);
        __m256d val  = low + high;

        const __m128d valupper = _mm256_extractf128_pd(val, 1);
        const __m128d rest = _mm256_castpd256_pd128(val);
        // Not in 512 _mm256_zeroupper();
        const __m128d valval = _mm_add_pd(valupper, rest);
        const __m128d res    = _mm_add_pd(_mm_permute_pd(valval, 1), valval);
        return _mm_cvtsd_f64(res);
#endif
    }

    inline double horizontalMul() const {
#ifdef __INTEL_COMPILER
        return _mm512_reduce_mul_pd(vec);
#else
        __m256d low  = _mm512_castpd512_pd256(vec);
        __m256d high = _mm512_extractf64x4_pd(vec, 1);
        __m256d val  = low * high;

        const __m128d valupper = _mm256_extractf128_pd(val, 1);
        const __m128d rest = _mm256_castpd256_pd128(val);
        // Not in 512 _mm256_zeroupper();
        const __m128d valval = _mm_mul_pd(valupper, rest);
        const __m128d res    = _mm_mul_pd(_mm_permute_pd(valval, 1), valval);
        return _mm_cvtsd_f64(res);
#endif
    }

    inline InaVecAVX512KNL sqrt() const {
        return _mm512_sqrt_pd(vec);
    }

    inline InaVecAVX512KNL exp() const {
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

         __m512d x = vec * COEFF_LOG2E;

         const __m512d fractional_part = x - InaVecAVX512KNL(x).floor().vec;

         __m512d factor = ((((((((COEFF_P5_X * fractional_part + COEFF_P5_Y)
                                 * fractional_part + COEFF_P5_Z) * fractional_part + COEFF_P5_A)
                                 * fractional_part + COEFF_P5_B) * fractional_part + COEFF_P5_C)
                                 * fractional_part + COEFF_P5_D) * fractional_part + COEFF_P5_E)
                                 * fractional_part + COEFF_P5_F);

         x -= factor;

         x = (COEFF_A * x + COEFF_B);

         alignas(64) double allvalreal[VecLength];
         _mm512_store_pd(allvalreal, x);

         alignas(64) long long int allvalint[VecLength] = { static_cast<long long int>(allvalreal[0]), static_cast<long long int>(allvalreal[1]),
                                                            static_cast<long long int>(allvalreal[2]), static_cast<long long int>(allvalreal[3]),
                                                            static_cast<long long int>(allvalreal[4]), static_cast<long long int>(allvalreal[5]),
                                                            static_cast<long long int>(allvalreal[6]), static_cast<long long int>(allvalreal[7]) };

         return _mm512_castsi512_pd(_mm512_load_epi64(reinterpret_cast<const __m512i*>(allvalint)));
    }

    inline InaVecAVX512KNL expLowAcc() const {
        const __m512d COEFF_LOG2E = _mm512_set1_pd(double(InaFastExp::CoeffLog2E()));
        const __m512d COEFF_A     = _mm512_set1_pd(double(InaFastExp::CoeffA64()));
        const __m512d COEFF_B     = _mm512_set1_pd(double(InaFastExp::CoeffB64()));
        const __m512d COEFF_P5_C  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_3()));
        const __m512d COEFF_P5_D  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_2()));
        const __m512d COEFF_P5_E  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_1()));
        const __m512d COEFF_P5_F  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_0()));

        __m512d x = vec * COEFF_LOG2E;

        const __m512d fractional_part = x - InaVecAVX512KNL(x).floor().vec;

        __m512d factor = (((COEFF_P5_C * fractional_part + COEFF_P5_D)
                           * fractional_part + COEFF_P5_E)
                           * fractional_part + COEFF_P5_F);

        x -= factor;

        x = (COEFF_A * x + COEFF_B);

        alignas(64) double allvalreal[VecLength];
        _mm512_store_pd(allvalreal, x);

        alignas(64) long long int allvalint[VecLength] = { static_cast<long long int>(allvalreal[0]), static_cast<long long int>(allvalreal[1]),
                                                           static_cast<long long int>(allvalreal[2]), static_cast<long long int>(allvalreal[3]),
                                                           static_cast<long long int>(allvalreal[4]), static_cast<long long int>(allvalreal[5]),
                                                           static_cast<long long int>(allvalreal[6]), static_cast<long long int>(allvalreal[7]) };

        return _mm512_castsi512_pd(_mm512_load_epi64(reinterpret_cast<const __m512i*>(allvalint)));
    }

    inline InaVecAVX512KNL rsqrt() const {
        // _mm512_rsqrt28_pd(vec) => 1E-10 error
        return _mm512_set1_pd(1) / _mm512_sqrt_pd(vec);
    }

    inline InaVecAVX512KNL abs() const {
        const __m512d minus0 = _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0x8000000000000000L)));
        return _mm512_castsi512_pd(_mm512_andnot_epi32(_mm512_castpd_si512(minus0), _mm512_castpd_si512(vec)));
    }

    inline InaVecAVX512KNL floor() const {
        return _mm512_cvt_roundps_pd(
            _mm256_cvtepi32_ps(
                _mm512_cvt_roundpd_epi32(vec, (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC))),
            _MM_FROUND_NO_EXC);
    }

    inline InaVecAVX512KNL signOf() const {
        const __m512d minus0 = _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0x8000000000000000L)));
        const __m512d signs  = _mm512_castsi512_pd(_mm512_and_epi32(_mm512_castpd_si512(vec), _mm512_castpd_si512(minus0)));
        return _mm512_maskz_mov_pd(
            _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_NEQ_OQ),
            _mm512_castsi512_pd(_mm512_or_epi32(_mm512_castpd_si512(signs),
                                                _mm512_castpd_si512(_mm512_set1_pd(1)))));
    }

    inline InaVecAVX512KNL isPositive() const {
        const __mmask8 greater = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_LE_OQ);
        const __m512d ones       = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(greater, ones);
    }

    inline InaVecAVX512KNL isNegative() const {
        const __mmask8 less = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_GE_OQ);
        const __m512d ones    = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(less, ones);
    }

    inline InaVecAVX512KNL isPositiveStrict() const {
        const __mmask8 greater = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_LT_OQ);
        const __m512d ones       = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(greater, ones);
    }

    inline InaVecAVX512KNL isNegativeStrict() const {
        const __mmask8 less = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_GT_OQ);
        const __m512d ones    = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(less, ones);
    }

    inline InaVecAVX512KNL isZero() const {
        const __mmask8 equalZero = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_EQ_OQ);
        const __m512d ones         = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(equalZero, ones);
    }

    inline InaVecAVX512KNL isNotZero() const {
        const __mmask8 equalZero = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_NEQ_OQ);
        const __m512d ones         = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(equalZero, ones);
    }

    inline InaVecMaskAVX512KNL<double> isPositiveMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_LE_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline InaVecMaskAVX512KNL<double> isNegativeMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_GE_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline InaVecMaskAVX512KNL<double> isPositiveStrictMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_LT_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline InaVecMaskAVX512KNL<double> isNegativeStrictMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_GT_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline InaVecMaskAVX512KNL<double> isZeroMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_EQ_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline InaVecMaskAVX512KNL<double> isNotZeroMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_NEQ_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    // Static basic methods
    inline static InaVecAVX512KNL GetZero() {
        return InaVecAVX512KNL(_mm512_setzero_pd());
    }

    inline static InaVecAVX512KNL GetOne() {
        return InaVecAVX512KNL(_mm512_set1_pd(1));
    }

    inline static InaVecAVX512KNL Min(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_min_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX512KNL Max(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_max_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX512KNL IsLowerOrEqual(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_LE_OQ);
        const __m512d ones          = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecAVX512KNL IsLower(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_LT_OQ);
        const __m512d ones          = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecAVX512KNL IsGreaterOrEqual(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_GE_OQ);
        const __m512d ones          = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecAVX512KNL IsGreater(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_GT_OQ);
        const __m512d ones          = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecAVX512KNL IsEqual(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_EQ_OQ);
        const __m512d ones          = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecAVX512KNL IsNotEqual(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ);
        const __m512d ones          = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecMaskAVX512KNL<double> IsLowerOrEqualMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_LE_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecMaskAVX512KNL<double> IsLowerMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_LT_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecMaskAVX512KNL<double> IsGreaterOrEqualMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_GE_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecMaskAVX512KNL<double> IsGreaterMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_GT_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecMaskAVX512KNL<double> IsEqualMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_EQ_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecMaskAVX512KNL<double> IsNotEqualMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ),
                                   _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecAVX512KNL BitsAnd(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(inVec1.vec), _mm512_castpd_si512(inVec2.vec)));
    }

    inline static InaVecAVX512KNL BitsNotAnd(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castsi512_pd(_mm512_andnot_si512(_mm512_castpd_si512(inVec1.vec), _mm512_castpd_si512(inVec2.vec)));
    }

    inline static InaVecAVX512KNL BitsOr(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castsi512_pd(_mm512_or_si512(_mm512_castpd_si512(inVec1.vec), _mm512_castpd_si512(inVec2.vec)));
    }

    inline static InaVecAVX512KNL BitsXor(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(inVec1.vec), _mm512_castpd_si512(inVec2.vec)));
    }

    inline static  const char* GetName() {
        return "InaVecAVX512KNL<double>";
    }

    inline static  InaIfElse< InaVecAVX512KNL<double> >::ThenClass If(const InaVecMaskAVX512KNL<double>& inTest) {
        return InaIfElse< InaVecAVX512KNL<double> >::IfClass().If(inTest);
    }

    inline static InaVecAVX512KNL IfElse(const InaVecMaskAVX512KNL<double>& inMask, const InaVecAVX512KNL& inIfTrue, const InaVecAVX512KNL& inIfFalse) {
        return _mm512_castsi512_pd(_mm512_or_si512(_mm512_castpd_si512(IfTrue(inMask, inIfTrue.vec).vec),
                      _mm512_castpd_si512(IfFalse(inMask, inIfFalse.vec).vec)));
    }

    inline static InaVecAVX512KNL IfTrue(const InaVecMaskAVX512KNL<double>& inMask, const InaVecAVX512KNL& inIfTrue) {
        return _mm512_castsi512_pd(_mm512_and_si512(inMask.getMask(), _mm512_castpd_si512(inIfTrue.vec)));
    }

    inline static InaVecAVX512KNL IfFalse(const InaVecMaskAVX512KNL<double>& inMask, const InaVecAVX512KNL& inIfFalse) {
        return _mm512_castsi512_pd(_mm512_andnot_si512(inMask.getMask(), _mm512_castpd_si512(inIfFalse.vec)));
    }

    // Inner operators
    inline InaVecAVX512KNL<double>& operator+=(const InaVecAVX512KNL<double>& inVec){
        vec = _mm512_add_pd(vec,inVec.vec);
        return *this;
    }

    inline InaVecAVX512KNL<double>& operator-=(const InaVecAVX512KNL<double>& inVec){
        vec = _mm512_sub_pd(vec,inVec.vec);
        return *this;
    }

    inline InaVecAVX512KNL<double>& operator/=(const InaVecAVX512KNL<double>& inVec){
        vec = _mm512_div_pd(vec,inVec.vec);
        return *this;
    }

    inline InaVecAVX512KNL<double>& operator*=(const InaVecAVX512KNL<double>& inVec){
        vec = _mm512_mul_pd(vec,inVec.vec);
        return *this;
    }

    inline InaVecAVX512KNL<double> operator-() const {
        const __m512d minus0 = _mm512_castsi512_pd(_mm512_set1_epi64(static_cast<long long>(0x8000000000000000L)));
        return _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(vec), _mm512_castpd_si512(minus0)));
    }

    inline InaVecAVX512KNL<double> pow(size_t power) const{
        return InaUtils::FastPow<InaVecAVX512KNL<double>>(vec, power);
    }
};

// Bits operators
inline InaVecAVX512KNL<double> operator&(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return InaVecAVX512KNL<double>::BitsAnd(inVec1, inVec2);
}

inline InaVecAVX512KNL<double> operator|(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return InaVecAVX512KNL<double>::BitsOr(inVec1, inVec2);
}

inline InaVecAVX512KNL<double> operator^(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return InaVecAVX512KNL<double>::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecAVX512KNL<double> operator+(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return _mm512_add_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512KNL<double> operator-(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return _mm512_sub_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512KNL<double> operator/(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return _mm512_div_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512KNL<double> operator*(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return _mm512_mul_pd(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskAVX512KNL<double> operator<(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return InaVecAVX512KNL<double>::IsLowerMask(inVec1,inVec2);
}

inline InaVecMaskAVX512KNL<double> operator<=(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return InaVecAVX512KNL<double>::IsLowerOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskAVX512KNL<double> operator>(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return InaVecAVX512KNL<double>::IsGreaterMask(inVec1,inVec2);
}

inline InaVecMaskAVX512KNL<double> operator>=(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return InaVecAVX512KNL<double>::IsGreaterOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskAVX512KNL<double> operator==(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return InaVecAVX512KNL<double>::IsEqualMask(inVec1,inVec2);
}

inline InaVecMaskAVX512KNL<double> operator!=(const InaVecAVX512KNL<double>& inVec1, const InaVecAVX512KNL<double>& inVec2){
    return InaVecAVX512KNL<double>::IsNotEqualMask(inVec1,inVec2);
}


#endif
