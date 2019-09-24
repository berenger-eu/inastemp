///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSVEFLOAT_HPP
#define INAVECSVEFLOAT_HPP

#include "InastempGlobal.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_SVE
#error InaVecSVE<float> is included but SVE is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <arm_sve.h>
#include <cmath>
#include <initializer_list>

// Forward declarations
template <class RealType>
class InaVecMaskSVE;

template <class RealType>
class InaVecSVE;

// Mask type
template <>
class alignas(32) InaVecMaskSVE<float> {
  svbool_t mask;
public:
  // Classic constructors
  inline InaVecMaskSVE(){}

  inline InaVecMaskSVE(const InaVecMaskSVE&) = default;
  inline InaVecMaskSVE& operator=(const InaVecMaskSVE&) = default;

  // Native data type compatibility
  inline /*not explicit*/ InaVecMaskSVE(const svbool_t inMask)
      : mask(inMask){}

  inline InaVecMaskSVE& operator=(const svbool_t inMask){
      mask = inMask;
      return (*this);
  }

  inline explicit operator svbool_t() const{
      return mask;
  }

  inline svbool_t getMask() const{
      return mask;
  }

  // Bool data type compatibility
  inline explicit InaVecMaskSVE(const bool inBool){
      mask = (inBool? svptrue_b32() : svpfalse_b32());
  }

  inline InaVecMaskSVE& operator=(const bool inBool){
      mask = (inBool? svptrue_b32() : svpfalse_b32());
      return (*this);
  }

  // Binary methods
  inline InaVecMaskSVE Not() const{
      return svnot_b32_z(mask);
  }

  inline bool isAllTrue() const{
      // true if all FF => !FF => 0 & FF => 0
      return svcnt_b32_z(mask) == svcntw();
  }

  inline bool isAllFalse() const{
      // true if all zero
      return svcnt_b32_z(mask) == 0;
  }

  // Double args methods
  inline static InaVecMaskSVE And(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
      return svand_b32(inMask1.mask,inMask2.mask);
  }

  inline static InaVecMaskSVE NotAnd(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
      return svbic_b32_z(inMask1.mask,inMask2.mask);
  }

  inline static InaVecMaskSVE Or(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
      return svorr_b32_z(inMask1.mask,inMask2.mask);
  }

  inline static InaVecMaskSVE Xor(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
      return sveor_b32_z(inMask1.mask,inMask2.mask);
  }

  inline static bool IsEqual(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
      return svcmpeq_b32_z(inMask1.mask,inMask2.mask);
  }

  inline static bool IsNotEqual(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
      return svcmpne_s32(inMask1.mask,inMask2.mask);
  }
};

// Mask must have operators
inline InaVecMaskSVE<float> operator&(const InaVecMaskSVE<float>& inMask1, const InaVecMaskSVE<float>& inMask2){
    return InaVecMaskSVE<float>::And(inMask1, inMask2);
}

inline InaVecMaskSVE<float> operator|(const InaVecMaskSVE<float>& inMask1, const InaVecMaskSVE<float>& inMask2){
    return InaVecMaskSVE<float>::Or(inMask1, inMask2);
}

inline InaVecMaskSVE<float> operator^(const InaVecMaskSVE<float>& inMask1, const InaVecMaskSVE<float>& inMask2){
    return InaVecMaskSVE<float>::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskSVE<float>& inMask1, const InaVecMaskSVE<float>& inMask2){
    return InaVecMaskSVE<float>::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskSVE<float>& inMask1, const InaVecMaskSVE<float>& inMask2){
    return InaVecMaskSVE<float>::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(32) InaVecSVE<float> {
protected:
    svfloat32_t vec;

public:
    using VecRawType           = svfloat32_t;
    using MaskType             = InaVecMaskSVE<float>;
    using RealType             = float;
    static const int VecLength = svcntw();
    static const int Alignement= 32;

    inline InaVecSVE(){}
    inline InaVecSVE(const InaVecSVE&) = default;
    inline InaVecSVE& operator = (const InaVecSVE&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecSVE(const svfloat32_t inVec)
        : vec(inVec){
    }

    inline InaVecSVE& operator=(const svfloat32_t inVec){
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const svfloat32_t inVec){
        vec = inVec;
    }

    inline explicit operator svfloat32_t() const{
        return vec;
    }

    inline svfloat32_t getVec() const{
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecSVE(const float val)
        : vec(svdup_f32(val)){
    }

    inline InaVecSVE& operator=(const float val){
        vec = svdup_f32(val);
        return *this;
    }

    inline void setFromScalar(const float val){
        vec = svdup_f32(val);
    }

    // Constructor from vec
    inline InaVecSVE(const std::initializer_list<float> lst)
        : InaVecSVE(lst.begin()){
    }

    inline explicit InaVecSVE(const float ptr[])
        : vec(svld1_f32_z(svptrue_b32(),ptr)){
    }

    inline InaVecSVE& setFromArray(const float ptr[]){
        vec = svld1_f32_z(svptrue_b32(),ptr);
        return *this;
    }

    inline InaVecSVE& setFromAlignedArray(const float ptr[]){
        vec = svld1_f32_z(svptrue_b32(),ptr);
        return *this;
    }

    inline InaVecSVE& setFromIndirectArray(const float values[], const int inIndirection[]) {
        vec = _mm256_set_ps(
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

    inline InaVecSVE& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        vec = _mm256_set_ps(
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
        svst1_f32(svptrue_b32(),ptr, vec);
    }

    inline void storeInAlignedArray(float ptr[]) const {
        svst1_f32(svptrue_b32(),ptr, vec);
    }

    // Acce to individual values
    inline float at(const int index) const {
        alignas(Alignement) float allval[VecLength];
        svst1_f32(svptrue_b32(),allval, vec);
        return allval[index];
    }

    // Horizontal operation
    inline float horizontalSum() const {
        return svadda_f32(svptrue_b32(),res);
    }

    inline float horizontalMul() const {
      return sfmla_f32(svptrue_b32(),res);
    }

    inline InaVecSVE sqrt() const {
        return svsqrt_f32_z(svptrue_b32(),vec);
    }

    inline InaVecSVE exp() const {
#ifdef __INTEL_COMPILER
        return _mm256_exp_ps(vec);
#else
        const svfloat32_t COEFF_LOG2E = svdup_f32(float(InaFastExp::CoeffLog2E()));
        const svfloat32_t COEFF_A     = svdup_f32(float(InaFastExp::CoeffA32()));
        const svfloat32_t COEFF_B     = svdup_f32(float(InaFastExp::CoeffB32()));
        const svfloat32_t COEFF_P5_A  = svdup_f32(float(InaFastExp::GetCoefficient6_5()));
        const svfloat32_t COEFF_P5_B  = svdup_f32(float(InaFastExp::GetCoefficient6_4()));
        const svfloat32_t COEFF_P5_C  = svdup_f32(float(InaFastExp::GetCoefficient6_3()));
        const svfloat32_t COEFF_P5_D  = svdup_f32(float(InaFastExp::GetCoefficient6_2()));
        const svfloat32_t COEFF_P5_E  = svdup_f32(float(InaFastExp::GetCoefficient6_1()));
        const svfloat32_t COEFF_P5_F  = svdup_f32(float(InaFastExp::GetCoefficient6_0()));

        svfloat32_t x = svmul_f32_z(vec, COEFF_LOG2E);

        const svfloat32_t fractional_part = svsub_f32_z(x, InaVecSVE(x).floor().vec);

        svfloat32_t factor = svadd_f32_z(svmul_f32_z(svadd_f32_z( svmul_f32_z(svadd_f32_z(
                         svmul_f32_z(svadd_f32_z( svmul_f32_z(svadd_f32_z(svmul_f32_z(
                         COEFF_P5_A, fractional_part), COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part), COEFF_P5_F);

        x = svsub_f32_z(x,factor);

        svbool_t castedInteger = _mm256_cvtps_epi32(svadd_f32_z(svmul_f32_z(COEFF_A, x), COEFF_B));

        return _mm256_castsi256_ps(castedInteger);
#endif
    }

    inline InaVecSVE expLowAcc() const {
        const svfloat32_t COEFF_LOG2E = svdup_f32(float(InaFastExp::CoeffLog2E()));
        const svfloat32_t COEFF_A     = svdup_f32(float(InaFastExp::CoeffA32()));
        const svfloat32_t COEFF_B     = svdup_f32(float(InaFastExp::CoeffB32()));
        const svfloat32_t COEFF_P5_D  = svdup_f32(float(InaFastExp::GetCoefficient3_2()));
        const svfloat32_t COEFF_P5_E  = svdup_f32(float(InaFastExp::GetCoefficient3_1()));
        const svfloat32_t COEFF_P5_F  = svdup_f32(float(InaFastExp::GetCoefficient3_0()));

        svfloat32_t x = svmul_f32_z(vec, COEFF_LOG2E);

        const svfloat32_t fractional_part = svsub_f32_z(x, InaVecSVE(x).floor().vec);

        svfloat32_t factor = svadd_f32_z(svmul_f32_z(
                         svadd_f32_z(svmul_f32_z(
                                         COEFF_P5_D, fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = svsub_f32_z(x,factor);

        svbool_t castedInteger = _mm256_cvtps_epi32(svadd_f32_z(svmul_f32_z(COEFF_A, x), COEFF_B));

        return _mm256_castsi256_ps(castedInteger);
    }

    inline InaVecSVE rsqrt() const {
        return svdup_f32(1) / svsqrt_f32_z(svptrue_b32(),vec); // _mm256_rsqrt_ps(val); not accurate enough
    }

    inline InaVecSVE abs() const {
        const svfloat32_t minus0 = _mm256_castsi256_ps(_mm256_set1_epi32(static_cast<int>(0x80000000)));
        return _mm256_andnot_ps(minus0, vec);
    }

    inline InaVecSVE floor() const {
        return _mm256_floor_ps(vec);
    }

    inline InaVecSVE signOf() const {
        const svfloat32_t minus0 = _mm256_castsi256_ps(_mm256_set1_epi32(static_cast<int>(0x80000000)));
        const svfloat32_t signs  = _mm256_and_ps(vec, minus0);
        return _mm256_and_ps(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_NEQ_OQ),
                             _mm256_or_ps(signs, svdup_f32(1)));
    }

    inline InaVecSVE isPositive() const {
        const svfloat32_t greater = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_LE_OQ);
        const svfloat32_t ones    = svdup_f32(1);
        return _mm256_and_ps(greater, ones);
    }

    inline InaVecSVE isNegative() const {
        const svfloat32_t less = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_GE_OQ);
        const svfloat32_t ones = svdup_f32(1);
        return _mm256_and_ps(less, ones);
    }

    inline InaVecSVE isPositiveStrict() const {
        const svfloat32_t greater = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_LT_OQ);
        const svfloat32_t ones    = svdup_f32(1);
        return _mm256_and_ps(greater, ones);
    }

    inline InaVecSVE isNegativeStrict() const {
        const svfloat32_t less = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_GT_OQ);
        const svfloat32_t ones = svdup_f32(1);
        return _mm256_and_ps(less, ones);
    }

    inline InaVecSVE isZero() const {
        const svfloat32_t equalZero = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_EQ_OQ);
        const svfloat32_t ones      = svdup_f32(1);
        return _mm256_and_ps(equalZero, ones);
    }

    inline InaVecSVE isNotZero() const {
        const svfloat32_t equalZero = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_NEQ_OQ);
        const svfloat32_t ones      = svdup_f32(1);
        return _mm256_and_ps(equalZero, ones);
    }

    inline InaVecMaskSVE<float> isPositiveMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_LE_OQ));
    }

    inline InaVecMaskSVE<float> isNegativeMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_GE_OQ));
    }

    inline InaVecMaskSVE<float> isPositiveStrictMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_LT_OQ));
    }

    inline InaVecMaskSVE<float> isNegativeStrictMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_GT_OQ));
    }

    inline InaVecMaskSVE<float> isZeroMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_EQ_OQ));
    }

    inline InaVecMaskSVE<float> isNotZeroMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_NEQ_OQ));
    }

    // Static basic methods
    inline static InaVecSVE GetZero() {
        return InaVecSVE(_mm256_setzero_ps());
    }

    inline static InaVecSVE GetOne() {
        return InaVecSVE(svdup_f32(1));
    }

    inline static InaVecSVE Min(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_min_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE Max(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_max_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE IsLowerOrEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat32_t testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_LE_OQ);
        const svfloat32_t ones       = svdup_f32(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecSVE IsLower(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat32_t testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_LT_OQ);
        const svfloat32_t ones       = svdup_f32(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecSVE IsGreaterOrEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat32_t testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_GE_OQ);
        const svfloat32_t ones       = svdup_f32(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecSVE IsGreater(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat32_t testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_GT_OQ);
        const svfloat32_t ones       = svdup_f32(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecSVE IsEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat32_t testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_EQ_OQ);
        const svfloat32_t ones       = svdup_f32(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecSVE IsNotEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat32_t testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ);
        const svfloat32_t ones       = svdup_f32(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecMaskSVE<float> IsLowerOrEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_LE_OQ));
    }

    inline static InaVecMaskSVE<float> IsLowerMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_LT_OQ));
    }

    inline static InaVecMaskSVE<float> IsGreaterOrEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_GE_OQ));
    }

    inline static InaVecMaskSVE<float> IsGreaterMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_GT_OQ));
    }

    inline static InaVecMaskSVE<float> IsEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_EQ_OQ));
    }

    inline static InaVecMaskSVE<float> IsNotEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ));
    }

    inline static InaVecSVE BitsAnd(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_and_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE BitsNotAnd(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_andnot_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE BitsOr(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_or_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE BitsXor(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_xor_ps(inVec1.vec, inVec2.vec);
    }

    inline static  const char* GetName() {
        return "InaVecSVE<float>";
    }

    inline static  InaIfElse< InaVecSVE<float> >::ThenClass If(const InaVecMaskSVE<float>& inTest) {
        return InaIfElse< InaVecSVE<float> >::IfClass().If(inTest);
    }

    inline static InaVecSVE IfElse(const InaVecMaskSVE<float>& inMask, const InaVecSVE& inIfTrue, const InaVecSVE& inIfFalse) {
        return _mm256_or_ps(IfTrue(inMask, inIfTrue.vec).vec,
                      IfFalse(inMask, inIfFalse.vec).vec);
    }

    inline static InaVecSVE IfTrue(const InaVecMaskSVE<float>& inMask, const InaVecSVE& inIfTrue) {
        return _mm256_and_ps(_mm256_castsi256_ps(inMask.getMask()), inIfTrue.vec);
    }

    inline static InaVecSVE IfFalse(const InaVecMaskSVE<float>& inMask, const InaVecSVE& inIfFalse) {
        return _mm256_andnot_ps(_mm256_castsi256_ps(inMask.getMask()), inIfFalse.vec);
    }

    // Inner operators
    inline InaVecSVE<float>& operator+=(const InaVecSVE<float>& inVec){
        vec = svadd_f32_z(vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<float>& operator-=(const InaVecSVE<float>& inVec){
        vec = svsub_f32_z(vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<float>& operator/=(const InaVecSVE<float>& inVec){
        vec = svdiv_f32_z(vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<float>& operator*=(const InaVecSVE<float>& inVec){
        vec = svmul_f32_z(vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<float> operator-() const {
        const svfloat32_t minus0 = _mm256_castsi256_ps(_mm256_set1_epi32(static_cast<int>(0x80000000)));
        return _mm256_xor_ps(vec, minus0);
    }

    inline InaVecSVE<float> pow(std::size_t power) const{
        return InaUtils::FastPow<InaVecSVE<float>>(*this, power);
    }

    // Multiple sum
    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecSVE<float>& inVec1,
                                          const InaVecSVE<float>& inVec2, const InaVecSVE<float>& inVec3,
                                          const InaVecSVE<float>& inVec4, const InaVecSVE<float>& inVec5,
                                          const InaVecSVE<float>& inVec6, const InaVecSVE<float>& inVec7,
                                          const InaVecSVE<float>& inVec8, Args ...args){
       const svfloat32_t val_a01_a23_b01_b23_a45_a67_b45_b67 = _mm256_hadd_ps(inVec1.vec, inVec2.vec);
        const svfloat32_t val_c01_c23_d01_d23_c45_c67_d45_d67 = _mm256_hadd_ps(inVec3.vec, inVec4.vec);

        const svfloat32_t val_e01_e23_f01_f23_e45_e67_f45_f67 = _mm256_hadd_ps(inVec5.vec, inVec6.vec);
        const svfloat32_t val_g01_g23_h01_h23_g45_g67_h45_h67 = _mm256_hadd_ps(inVec7.vec, inVec8.vec);

        const svfloat32_t val_a0123_b01b23_c0123_d01b23_a4567_b4567_c4567_d4567 = _mm256_hadd_ps(val_a01_a23_b01_b23_a45_a67_b45_b67,
                                                           val_c01_c23_d01_d23_c45_c67_d45_d67);

        const svfloat32_t val_e0123_f01b23_g0123_h01b23_e4567_f4567_g4567_h4567 = _mm256_hadd_ps(val_e01_e23_f01_f23_e45_e67_f45_f67,
                                                           val_g01_g23_h01_h23_g45_g67_h45_h67);

        const svfloat32_t val_a0123_b01b23_c0123_d01b23_e0123_f01b23_g0123_h01b23 =
                                            _mm256_permute2f128_ps(val_a0123_b01b23_c0123_d01b23_a4567_b4567_c4567_d4567,
                                                                   val_e0123_f01b23_g0123_h01b23_e4567_f4567_g4567_h4567, 0x20);// 010.0000
        const svfloat32_t val_a4567_b4567_c4567_d4567_e4567_f4567_g4567_h4567 =
                                             _mm256_permute2f128_ps(val_a0123_b01b23_c0123_d01b23_a4567_b4567_c4567_d4567,
                                                                   val_e0123_f01b23_g0123_h01b23_e4567_f4567_g4567_h4567, 0x31);// 000.0001

        const svfloat32_t sum_a_b_c_d_e_f_g_h = val_a0123_b01b23_c0123_d01b23_e0123_f01b23_g0123_h01b23 + val_a4567_b4567_c4567_d4567_e4567_f4567_g4567_h4567;

        svfloat32_t vecBuffer = svld1_f32_z(svptrue_b32(),sumRes);
        vecBuffer += sum_a_b_c_d_e_f_g_h;
        svst1_f32(svptrue_b32(),sumRes, vecBuffer);

        MultiHorizontalSum(&sumRes[8], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecSVE<float>& inVec1,
                                          const InaVecSVE<float>& inVec2, const InaVecSVE<float>& inVec3,
                                          const InaVecSVE<float>& inVec4, Args ...args){
        const svfloat32_t val_a01_a23_b01_b23_a45_a67_b45_b67 = _mm256_hadd_ps(inVec1.vec, inVec2.vec);
        const svfloat32_t val_c01_c23_d01_d23_c45_c67_d45_d67 = _mm256_hadd_ps(inVec3.vec, inVec4.vec);

        const svfloat32_t val_a0123_b01b23_c0123_d01b23_a4567_b4567_c4567_d4567 = _mm256_hadd_ps(val_a01_a23_b01_b23_a45_a67_b45_b67,
                                                           val_c01_c23_d01_d23_c45_c67_d45_d67);

        __m128 valupper = _mm256_extractf128_ps(val_a0123_b01b23_c0123_d01b23_a4567_b4567_c4567_d4567, 1);
        __m128 vallower = _mm256_castps256_ps128(val_a0123_b01b23_c0123_d01b23_a4567_b4567_c4567_d4567);

        __m128 vecBuffer = _mm_loadu_ps(sumRes);
        vecBuffer += valupper + vallower;
        _mm_storeu_ps(sumRes, vecBuffer);

        MultiHorizontalSum(&sumRes[4], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecSVE<float>& inVec1,
                                          const InaVecSVE<float>& inVec2, Args ...args){

        const svfloat32_t val_a01_a23_b01_b23_a45_a67_b45_b67 = _mm256_hadd_ps(inVec1.vec, inVec2.vec);

        const __m128 valupper = _mm256_extractf128_ps(val_a01_a23_b01_b23_a45_a67_b45_b67, 1);
        const __m128 vallower = _mm256_castps256_ps128(val_a01_a23_b01_b23_a45_a67_b45_b67);

        const __m128 val_a0123_b0123_a4567_b4567 = _mm_hadd_ps(valupper, vallower);

        const __m128 val_a4567_b4567_a0123_b0123 = _mm_shuffle_ps(val_a0123_b0123_a4567_b4567, val_a0123_b0123_a4567_b4567, 0x9E);// 10.01.11.10

        const __m128 val_suma_x_sumb_x = _mm_add_ps(val_a0123_b0123_a4567_b4567, val_a4567_b4567_a0123_b0123);

        alignas(Alignement) float buffer[VecLength] = {0};
        buffer[0] = sumRes[0];
        buffer[1] = sumRes[1];
        __m128 vecBuffer = _mm_load_ps(buffer);
        vecBuffer += val_suma_x_sumb_x;
        _mm_store_ps(buffer, vecBuffer);
        sumRes[0] = buffer[0];
        sumRes[1] = buffer[1];

        MultiHorizontalSum(&sumRes[2], args... );
    }

    inline static void MultiHorizontalSum(float sumRes[], const InaVecSVE<float>& inVec){
        sumRes[0] += inVec.horizontalSum();
    }

    inline static void MultiHorizontalSum(float /*sumRes*/[]){
    }
};

// Bits operators
inline InaVecSVE<float> operator&(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return InaVecSVE<float>::BitsAnd(inVec1, inVec2);
}

inline InaVecSVE<float> operator|(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return InaVecSVE<float>::BitsOr(inVec1, inVec2);
}

inline InaVecSVE<float> operator^(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return InaVecSVE<float>::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecSVE<float> operator+(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return svadd_f32_z(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSVE<float> operator-(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return svsub_f32_z(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSVE<float> operator/(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return svdiv_f32_z(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSVE<float> operator*(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return svmul_f32_z(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskSVE<float> operator<(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return InaVecSVE<float>::IsLowerMask(inVec1,inVec2);
}

inline InaVecMaskSVE<float> operator<=(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return InaVecSVE<float>::IsLowerOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskSVE<float> operator>(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return InaVecSVE<float>::IsGreaterMask(inVec1,inVec2);
}

inline InaVecMaskSVE<float> operator>=(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return InaVecSVE<float>::IsGreaterOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskSVE<float> operator==(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return InaVecSVE<float>::IsEqualMask(inVec1,inVec2);
}

inline InaVecMaskSVE<float> operator!=(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return InaVecSVE<float>::IsNotEqualMask(inVec1,inVec2);
}


#endif
