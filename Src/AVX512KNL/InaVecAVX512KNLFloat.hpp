///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVX512KNLFLOAT_HPP
#define INAVECAVX512KNLFLOAT_HPP

#include "InastempConfig.h"
#include "InaAVX512KNLOperators.hpp"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_AVX512KNL
#error InaVecAVX512KNL512KNL<float> is included but AVX512KNL is not enable in the configuration
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
class alignas(64) InaVecMaskAVX512KNL<float> {
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
        mask = (inBool? _mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)) : _mm512_setzero_si512());
    }

    inline InaVecMaskAVX512KNL& operator=(const bool inBool){
        mask = (inBool? _mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)) : _mm512_setzero_si512());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskAVX512KNL Not() const{
        return NotAnd(mask, _mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)));
    }


    inline bool isAllTrue() const{
        // true if all FF => !FF => 0 & FF => 0
        // sde bug here with _mm512_cmp_epu32_mask return 0xFF instead of 0xFFFF so use 64 bits version
        const __mmask8 testResult = _mm512_cmp_epu64_mask(mask, _mm512_set1_epi64(0xFFFFFFFFFFFFFFFFUL), _MM_CMPINT_EQ);
        return testResult == 0xFF;
    }

    inline bool isAllFalse() const{
        // true if all zero
        // sde bug here with _mm512_cmp_epu32_mask return 0xFF instead of 0xFFFF so use 64 bits version
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
        // sde bug here with _mm512_cmp_epu32_mask return 0xFF instead of 0xFFFF so use 64 bits version
        const __mmask8 testResult = _mm512_cmp_epu64_mask(inMask1.mask, inMask2.mask, _MM_CMPINT_EQ);
        return testResult == 0xFF;
    }

    inline static bool IsNotEqual(const InaVecMaskAVX512KNL& inMask1, const InaVecMaskAVX512KNL& inMask2){
        // sde bug here with _mm512_cmp_epu32_mask return 0xFF instead of 0xFFFF so use 64 bits version
        const __mmask8 testResult = _mm512_cmp_epu64_mask(inMask1.mask, inMask2.mask, _MM_CMPINT_EQ);
        return testResult != 0xFF;
    }
};

// Mask must have operators
inline InaVecMaskAVX512KNL<float> operator&(const InaVecMaskAVX512KNL<float>& inMask1, const InaVecMaskAVX512KNL<float>& inMask2){
    return InaVecMaskAVX512KNL<float>::And(inMask1, inMask2);
}

inline InaVecMaskAVX512KNL<float> operator|(const InaVecMaskAVX512KNL<float>& inMask1, const InaVecMaskAVX512KNL<float>& inMask2){
    return InaVecMaskAVX512KNL<float>::Or(inMask1, inMask2);
}

inline InaVecMaskAVX512KNL<float> operator^(const InaVecMaskAVX512KNL<float>& inMask1, const InaVecMaskAVX512KNL<float>& inMask2){
    return InaVecMaskAVX512KNL<float>::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskAVX512KNL<float>& inMask1, const InaVecMaskAVX512KNL<float>& inMask2){
    return InaVecMaskAVX512KNL<float>::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskAVX512KNL<float>& inMask1, const InaVecMaskAVX512KNL<float>& inMask2){
    return InaVecMaskAVX512KNL<float>::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(64) InaVecAVX512KNL<float> {
protected:
    __m512 vec;

public:
    using VecRawType           = __m512;
    using MaskType             = InaVecMaskAVX512KNL<float>;
    using RealType             = float;
    static const int VecLength = 16;
    static const int Alignement= 64;

    inline InaVecAVX512KNL(){}
    inline InaVecAVX512KNL(const InaVecAVX512KNL&) = default;
    inline InaVecAVX512KNL& operator = (const InaVecAVX512KNL&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecAVX512KNL(const __m512 inVec)
        : vec(inVec){
    }

    inline InaVecAVX512KNL& operator=(const __m512 inVec){
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __m512 inVec){
        vec = inVec;
    }

    inline explicit operator __m512() const{
        return vec;
    }

    inline __m512 getVec() const{
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecAVX512KNL(const float val)
        : vec(_mm512_set1_ps(val)){
    }

    inline InaVecAVX512KNL& operator=(const float val){
        vec = _mm512_set1_ps(val);
        return *this;
    }

    inline void setFromScalar(const float val){
        vec = _mm512_set1_ps(val);
    }

    // Constructor from vec
    inline explicit InaVecAVX512KNL(const float ptr[])
        : vec(_mm512_loadu_ps(ptr)){
    }

    inline InaVecAVX512KNL& setFromArray(const float ptr[]){
        vec = _mm512_loadu_ps(ptr);
        return *this;
    }

    inline InaVecAVX512KNL& setFromAlignedArray(const float ptr[]){
        vec = _mm512_load_ps(ptr);
        return *this;
    }

    inline InaVecAVX512KNL& setFromIndirectArray(const float values[], const int inIndirection[]) {
        vec = _mm512_set_ps(
                    values[inIndirection[15]],
                    values[inIndirection[14]],
                    values[inIndirection[13]],
                    values[inIndirection[12]],
                    values[inIndirection[11]],
                    values[inIndirection[10]],
                    values[inIndirection[9]],
                    values[inIndirection[8]],
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

    inline InaVecAVX512KNL& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        vec = _mm512_set_ps(
                    inArray[inIndirection1[15] * inLeadingDimension + inIndirection2[15]],
                    inArray[inIndirection1[14] * inLeadingDimension + inIndirection2[14]],
                    inArray[inIndirection1[13] * inLeadingDimension + inIndirection2[13]],
                    inArray[inIndirection1[12] * inLeadingDimension + inIndirection2[12]],
                    inArray[inIndirection1[11] * inLeadingDimension + inIndirection2[11]],
                    inArray[inIndirection1[10] * inLeadingDimension + inIndirection2[10]],
                    inArray[inIndirection1[9] * inLeadingDimension + inIndirection2[9]],
                    inArray[inIndirection1[8] * inLeadingDimension + inIndirection2[8]],
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
    inline void storeInArray(float ptr[]) const {
        _mm512_storeu_ps(ptr, vec);
    }

    inline void storeInAlignedArray(float ptr[]) const {
        _mm512_store_ps(ptr, vec);
    }

    // Acce to individual values
    inline float at(const int index) const {
        alignas(Alignement) float allval[VecLength];
        _mm512_store_ps(allval, vec);
        return allval[index];
    }

    // Horizontal operation
    inline float horizontalSum() const {
#ifdef __INTEL_COMPILER
        return _mm512_reduce_add_ps(vec);
#else
        __m256 low  = _mm512_castps512_ps256(vec);
        __m256 high = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(vec), 1));
        __m256 val  = low + high;

        const __m128 valupper = _mm256_extractf128_ps(val, 1);
        _mm256_zeroupper(); // Could be moved after the _mm256_extractf128_ps
        const __m128 valval = _mm_add_ps(valupper,
                                         _mm256_extractf128_ps(val, 0));
        __m128 valsum = _mm_add_ps(_mm_permute_ps(valval, 0x1B), valval);
        __m128 res    = _mm_add_ps(_mm_permute_ps(valsum, 0xB1), valsum);
        return _mm_cvtss_f32(res);
#endif
    }

    inline float horizontalMul() const {
#ifdef __INTEL_COMPILER
        return _mm512_reduce_mul_ps(vec);
#else
        __m256 low  = _mm512_castps512_ps256(vec);
        __m256 high = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(vec), 1));
        __m256 val  = low * high;

        const __m128 valupper = _mm256_extractf128_ps(val, 1);
        _mm256_zeroupper(); // Could be moved after the _mm256_extractf128_ps
        const __m128 valval = _mm_mul_ps(valupper,
                                         _mm256_extractf128_ps(val, 0));
        __m128 valsum = _mm_mul_ps(_mm_permute_ps(valval, 0x1B), valval);
        __m128 res    = _mm_mul_ps(_mm_permute_ps(valsum, 0xB1), valsum);
        return _mm_cvtss_f32(res);
#endif
    }

    inline InaVecAVX512KNL sqrt() const {
        return _mm512_sqrt_ps(vec);
    }

    inline InaVecAVX512KNL exp() const {
         const __m512 COEFF_LOG2E = _mm512_set1_ps(float(InaFastExp::CoeffLog2E()));
         const __m512 COEFF_A     = _mm512_set1_ps(float(InaFastExp::CoeffA32()));
         const __m512 COEFF_B     = _mm512_set1_ps(float(InaFastExp::CoeffB32()));
         const __m512 COEFF_P5_A  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_5()));
         const __m512 COEFF_P5_B  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_4()));
         const __m512 COEFF_P5_C  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_3()));
         const __m512 COEFF_P5_D  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_2()));
         const __m512 COEFF_P5_E  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_1()));
         const __m512 COEFF_P5_F  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_0()));

         __m512 x = vec * COEFF_LOG2E;

         const __m512 fractional_part = x - InaVecAVX512KNL(x).floor().vec;

         __m512 factor = (((((COEFF_P5_A * fractional_part + COEFF_P5_B)
                             * fractional_part + COEFF_P5_C)
                             * fractional_part + COEFF_P5_D)
                             * fractional_part + COEFF_P5_E)
                             * fractional_part + COEFF_P5_F);

         x -= factor;

         __m512i castedInteger = _mm512_cvtps_epi32(COEFF_A * x + COEFF_B);

         return _mm512_castsi512_ps(castedInteger);
    }

    inline InaVecAVX512KNL expLowAcc() const {
        const __m512 COEFF_LOG2E = _mm512_set1_ps(float(InaFastExp::CoeffLog2E()));
        const __m512 COEFF_A     = _mm512_set1_ps(float(InaFastExp::CoeffA32()));
        const __m512 COEFF_B     = _mm512_set1_ps(float(InaFastExp::CoeffB32()));
        const __m512 COEFF_P5_D  = _mm512_set1_ps(float(InaFastExp::GetCoefficient3_2()));
        const __m512 COEFF_P5_E  = _mm512_set1_ps(float(InaFastExp::GetCoefficient3_1()));
        const __m512 COEFF_P5_F  = _mm512_set1_ps(float(InaFastExp::GetCoefficient3_0()));

        __m512 x = vec * COEFF_LOG2E;

        const __m512 fractional_part = x - InaVecAVX512KNL(x).floor().vec;

        __m512 factor = ((COEFF_P5_D * fractional_part + COEFF_P5_E) * fractional_part + COEFF_P5_F);

        x -= factor;

        __m512i castedInteger = _mm512_cvtps_epi32(COEFF_A * x + COEFF_B);

        return _mm512_castsi512_ps(castedInteger);
    }

    inline InaVecAVX512KNL rsqrt() const {
        return _mm512_rsqrt28_ps(vec);
    }

    inline InaVecAVX512KNL abs() const {
        const __m512 minus0 = _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0x80000000)));
        return _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(minus0), _mm512_castps_si512(vec)));
    }

    inline InaVecAVX512KNL floor() const {
        return _mm512_cvt_roundepi32_ps(
            _mm512_cvt_roundps_epi32(vec, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)),
            (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    }

    inline InaVecAVX512KNL signOf() const {
        const __m512 minus0 = _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0x80000000)));
        const __m512 signs  = _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(vec), _mm512_castps_si512(minus0)));
        return _mm512_maskz_mov_ps(
            _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_NEQ_OQ),
            _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(signs),
                                                _mm512_castps_si512(_mm512_set1_ps(1)))));
    }

    inline InaVecAVX512KNL isPositive() const {
        const __mmask16 greater = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_LE_OQ);
        const __m512 ones       = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(greater, ones);
    }

    inline InaVecAVX512KNL isNegative() const {
        const __mmask16 less = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_GE_OQ);
        const __m512 ones    = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(less, ones);
    }

    inline InaVecAVX512KNL isPositiveStrict() const {
        const __mmask16 greater = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_LT_OQ);
        const __m512 ones       = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(greater, ones);
    }

    inline InaVecAVX512KNL isNegativeStrict() const {
        const __mmask16 less = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_GT_OQ);
        const __m512 ones    = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(less, ones);
    }

    inline InaVecAVX512KNL isZero() const {
        const __mmask16 equalZero = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_EQ_OQ);
        const __m512 ones         = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(equalZero, ones);
    }

    inline InaVecAVX512KNL isNotZero() const {
        const __mmask16 equalZero = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_NEQ_OQ);
        const __m512 ones         = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(equalZero, ones);
    }

    inline InaVecMaskAVX512KNL<float> isPositiveMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_LE_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline InaVecMaskAVX512KNL<float> isNegativeMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_GE_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline InaVecMaskAVX512KNL<float> isPositiveStrictMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_LT_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline InaVecMaskAVX512KNL<float> isNegativeStrictMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_GT_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline InaVecMaskAVX512KNL<float> isZeroMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_EQ_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline InaVecMaskAVX512KNL<float> isNotZeroMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_NEQ_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    // Static basic methods
    inline static InaVecAVX512KNL GetZero() {
        return InaVecAVX512KNL(_mm512_setzero_ps());
    }

    inline static InaVecAVX512KNL GetOne() {
        return InaVecAVX512KNL(_mm512_set1_ps(1));
    }

    inline static InaVecAVX512KNL Min(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_min_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX512KNL Max(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_max_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX512KNL IsLowerOrEqual(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_LE_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecAVX512KNL IsLower(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_LT_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecAVX512KNL IsGreaterOrEqual(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_GE_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecAVX512KNL IsGreater(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_GT_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecAVX512KNL IsEqual(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_EQ_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecAVX512KNL IsNotEqual(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecMaskAVX512KNL<float> IsLowerOrEqualMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_LE_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline static InaVecMaskAVX512KNL<float> IsLowerMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_LT_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline static InaVecMaskAVX512KNL<float> IsGreaterOrEqualMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_GE_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline static InaVecMaskAVX512KNL<float> IsGreaterMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_GT_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline static InaVecMaskAVX512KNL<float> IsEqualMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_EQ_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline static InaVecMaskAVX512KNL<float> IsNotEqualMask(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ),
                                   _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0xFFFFFFFF)))));
    }

    inline static InaVecAVX512KNL BitsAnd(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(inVec1.vec), _mm512_castps_si512(inVec2.vec)));
    }

    inline static InaVecAVX512KNL BitsNotAnd(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castsi512_ps(_mm512_andnot_si512(_mm512_castps_si512(inVec1.vec), _mm512_castps_si512(inVec2.vec)));
    }

    inline static InaVecAVX512KNL BitsOr(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castsi512_ps(_mm512_or_si512(_mm512_castps_si512(inVec1.vec), _mm512_castps_si512(inVec2.vec)));
    }

    inline static InaVecAVX512KNL BitsXor(const InaVecAVX512KNL& inVec1, const InaVecAVX512KNL& inVec2) {
        return _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(inVec1.vec), _mm512_castps_si512(inVec2.vec)));
    }

    inline static  const char* GetName() {
        return "InaVecAVX512KNL<float>";
    }

    inline static  InaIfElse< InaVecAVX512KNL<float> >::ThenClass If(const InaVecMaskAVX512KNL<float>& inTest) {
        return InaIfElse< InaVecAVX512KNL<float> >::IfClass().If(inTest);
    }

    inline static InaVecAVX512KNL IfElse(const InaVecMaskAVX512KNL<float>& inMask, const InaVecAVX512KNL& inIfTrue, const InaVecAVX512KNL& inIfFalse) {
        return _mm512_castsi512_ps(_mm512_or_si512(_mm512_castps_si512(IfTrue(inMask, inIfTrue.vec).vec),
                      _mm512_castps_si512(IfFalse(inMask, inIfFalse.vec).vec)));
    }

    inline static InaVecAVX512KNL IfTrue(const InaVecMaskAVX512KNL<float>& inMask, const InaVecAVX512KNL& inIfTrue) {
        return _mm512_castsi512_ps(_mm512_and_si512(inMask.getMask(), _mm512_castps_si512(inIfTrue.vec)));
    }

    inline static InaVecAVX512KNL IfFalse(const InaVecMaskAVX512KNL<float>& inMask, const InaVecAVX512KNL& inIfFalse) {
        return _mm512_castsi512_ps(_mm512_andnot_si512(inMask.getMask(), _mm512_castps_si512(inIfFalse.vec)));
    }

    // Inner operators
    inline InaVecAVX512KNL<float>& operator+=(const InaVecAVX512KNL<float>& inVec){
        vec = _mm512_add_ps(vec,inVec.vec);
        return *this;
    }

    inline InaVecAVX512KNL<float>& operator-=(const InaVecAVX512KNL<float>& inVec){
        vec = _mm512_sub_ps(vec,inVec.vec);
        return *this;
    }

    inline InaVecAVX512KNL<float>& operator/=(const InaVecAVX512KNL<float>& inVec){
        vec = _mm512_div_ps(vec,inVec.vec);
        return *this;
    }

    inline InaVecAVX512KNL<float>& operator*=(const InaVecAVX512KNL<float>& inVec){
        vec = _mm512_mul_ps(vec,inVec.vec);
        return *this;
    }

    inline InaVecAVX512KNL<float> operator-() const {
        const __m512 minus0 = _mm512_castsi512_ps(_mm512_set1_epi32(static_cast<int>(0x80000000)));
        return _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(vec), _mm512_castps_si512(minus0)));
    }

    inline InaVecAVX512KNL<float> pow(size_t power) const{
        return InaUtils::FastPow<InaVecAVX512KNL<float>>(vec, power);
    }
};

// Bits operators
inline InaVecAVX512KNL<float> operator&(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return InaVecAVX512KNL<float>::BitsAnd(inVec1, inVec2);
}

inline InaVecAVX512KNL<float> operator|(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return InaVecAVX512KNL<float>::BitsOr(inVec1, inVec2);
}

inline InaVecAVX512KNL<float> operator^(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return InaVecAVX512KNL<float>::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecAVX512KNL<float> operator+(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return _mm512_add_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512KNL<float> operator-(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return _mm512_sub_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512KNL<float> operator/(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return _mm512_div_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512KNL<float> operator*(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return _mm512_mul_ps(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskAVX512KNL<float> operator<(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return InaVecAVX512KNL<float>::IsLowerMask(inVec1,inVec2);
}

inline InaVecMaskAVX512KNL<float> operator<=(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return InaVecAVX512KNL<float>::IsLowerOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskAVX512KNL<float> operator>(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return InaVecAVX512KNL<float>::IsGreaterMask(inVec1,inVec2);
}

inline InaVecMaskAVX512KNL<float> operator>=(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return InaVecAVX512KNL<float>::IsGreaterOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskAVX512KNL<float> operator==(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return InaVecAVX512KNL<float>::IsEqualMask(inVec1,inVec2);
}

inline InaVecMaskAVX512KNL<float> operator!=(const InaVecAVX512KNL<float>& inVec1, const InaVecAVX512KNL<float>& inVec2){
    return InaVecAVX512KNL<float>::IsNotEqualMask(inVec1,inVec2);
}


#endif
