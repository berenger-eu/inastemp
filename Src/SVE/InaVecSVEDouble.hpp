///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSVEDOUBLE_HPP
#define INAVECSVEDOUBLE_HPP

#include "InastempGlobal.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_SVE
#error InaVecSVE<double> is included but SVE is not enable in the configuration
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
class alignas(32) InaVecMaskSVE<double> {
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
        mask = (inBool? svptrue_b64() : svpfalse_b64());
    }

    inline InaVecMaskSVE& operator=(const bool inBool){
        mask = (inBool? svptrue_b64() : svpfalse_b64());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskSVE Not() const{
        return svnot_b64_z(mask);
    }

    inline bool isAllTrue() const{
        // true if all FF => !FF => 0 & FF => 0
        return svcnt_b64_z(mask) == svcntd();
    }

    inline bool isAllFalse() const{
        // true if all zero
        return svcnt_b64_z(mask) == 0;
    }

    // Double args methods
    inline static InaVecMaskSVE And(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return svand_b64(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSVE NotAnd(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return svbic_b64_z(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSVE Or(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return svorr_b64_z(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSVE Xor(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return sveor_b64_z(inMask1.mask,inMask2.mask);
    }

    inline static bool IsEqual(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return svcmpeq_b64_z(inMask1.mask,inMask2.mask);
    }

    inline static bool IsNotEqual(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return svcmpne_s64(inMask1.mask,inMask2.mask);
    }
};

// Mask must have operators
inline InaVecMaskSVE<double> operator&(const InaVecMaskSVE<double>& inMask1, const InaVecMaskSVE<double>& inMask2){
    return InaVecMaskSVE<double>::And(inMask1, inMask2);
}

inline InaVecMaskSVE<double> operator|(const InaVecMaskSVE<double>& inMask1, const InaVecMaskSVE<double>& inMask2){
    return InaVecMaskSVE<double>::Or(inMask1, inMask2);
}

inline InaVecMaskSVE<double> operator^(const InaVecMaskSVE<double>& inMask1, const InaVecMaskSVE<double>& inMask2){
    return InaVecMaskSVE<double>::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskSVE<double>& inMask1, const InaVecMaskSVE<double>& inMask2){
    return InaVecMaskSVE<double>::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskSVE<double>& inMask1, const InaVecMaskSVE<double>& inMask2){
    return InaVecMaskSVE<double>::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(32) InaVecSVE<double> {
protected:
    svfloat64_t vec;

public:
    using VecRawType           = svfloat64_t;
    using MaskType             = InaVecMaskSVE<double>;
    using RealType             = double;
    static const int VecLength = svcntd();
    static const int Alignement= 32;

    inline InaVecSVE(){}
    inline InaVecSVE(const InaVecSVE&) = default;
    inline InaVecSVE& operator = (const InaVecSVE&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecSVE(const svfloat64_t inVec)
        : vec(inVec){
    }

    inline InaVecSVE& operator=(const svfloat64_t inVec){
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const svfloat64_t inVec){
        vec = inVec;
    }

    inline explicit operator svfloat64_t() const{
        return vec;
    }

    inline svfloat64_t getVec() const{
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecSVE(const double val)
        : vec(svdup_f64(val)){
    }

    inline InaVecSVE& operator=(const double val){
        vec = svdup_f64(val);
        return *this;
    }

    inline void setFromScalar(const double val){
        vec = svdup_f64(val);
    }

    // Constructor from vec
    inline InaVecSVE(const std::initializer_list<double> lst)
        : InaVecSVE(lst.begin()){
    }

    inline explicit InaVecSVE(const double ptr[])
        : vec(svld1_f64_z(svptrue_b64(),ptr)){
    }

    inline InaVecSVE& setFromArray(const double ptr[]){
        vec = svld1_f64_z(svptrue_b64(),ptr);
        return *this;
    }

    inline InaVecSVE& setFromAlignedArray(const double ptr[]){
        vec = svld1_f64_z(svptrue_b64(),ptr);
        return *this;
    }

    inline InaVecSVE& setFromIndirectArray(const double values[], const int inIndirection[]) {
        vec = _mm256_set_pd(
                    values[inIndirection[3]],
                    values[inIndirection[2]],
                    values[inIndirection[1]],
                    values[inIndirection[0]]);
        return *this;
    }

    inline InaVecSVE& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        vec = _mm256_set_pd(
                    inArray[inIndirection1[3] * inLeadingDimension + inIndirection2[3]],
                    inArray[inIndirection1[2] * inLeadingDimension + inIndirection2[2]],
                    inArray[inIndirection1[1] * inLeadingDimension + inIndirection2[1]],
                    inArray[inIndirection1[0] * inLeadingDimension + inIndirection2[0]]);
        return *this;
    }

    // Move back to array
    inline void storeInArray(double ptr[]) const {
        svst1_f64(svptrue_b64(),ptr, vec);
    }

    inline void storeInAlignedArray(double ptr[]) const {
        svst1_f64(svptrue_b64(),ptr, vec);
    }

    // Acce to individual values
    inline double at(const int index) const {
        alignas(Alignement) double allval[VecLength];
        svst1_f64(svptrue_b64(),allval, vec);
        return allval[index];
    }

    // Horizontal operation
    inline double horizontalSum() const {
      return svadda_f64(svptrue_b64(),res);
    }

    inline double horizontalMul() const {
      return sfmla_f64(svptrue_b64(),res);
    }

    inline InaVecSVE sqrt() const {
        return svsqrt_f64_z(svptrue_b64(),vec);
    }

    inline InaVecSVE exp() const {
#ifdef __INTEL_COMPILER
        return _mm256_exp_pd(vec);
#else
        const svfloat64_t COEFF_LOG2E = svdup_f64(double(InaFastExp::CoeffLog2E()));
        const svfloat64_t COEFF_A     = svdup_f64(double(InaFastExp::CoeffA64()));
        const svfloat64_t COEFF_B     = svdup_f64(double(InaFastExp::CoeffB64()));
        const svfloat64_t COEFF_P5_X  = svdup_f64(double(InaFastExp::GetCoefficient9_8()));
        const svfloat64_t COEFF_P5_Y  = svdup_f64(double(InaFastExp::GetCoefficient9_7()));
        const svfloat64_t COEFF_P5_Z  = svdup_f64(double(InaFastExp::GetCoefficient9_6()));
        const svfloat64_t COEFF_P5_A  = svdup_f64(double(InaFastExp::GetCoefficient9_5()));
        const svfloat64_t COEFF_P5_B  = svdup_f64(double(InaFastExp::GetCoefficient9_4()));
        const svfloat64_t COEFF_P5_C  = svdup_f64(double(InaFastExp::GetCoefficient9_3()));
        const svfloat64_t COEFF_P5_D  = svdup_f64(double(InaFastExp::GetCoefficient9_2()));
        const svfloat64_t COEFF_P5_E  = svdup_f64(double(InaFastExp::GetCoefficient9_1()));
        const svfloat64_t COEFF_P5_F  = svdup_f64(double(InaFastExp::GetCoefficient9_0()));

        svfloat64_t x = svmul_f64_z(vec, COEFF_LOG2E);

        const svfloat64_t fractional_part = svsub_f64_z(x, InaVecSVE(x).floor().vec);

        svfloat64_t factor = svadd_f64_z(svmul_f64_z(svadd_f64_z(
                         svmul_f64_z(svadd_f64_z( svmul_f64_z(svadd_f64_z(
                         svmul_f64_z(svadd_f64_z( svmul_f64_z(svadd_f64_z(
                         svmul_f64_z(svadd_f64_z( svmul_f64_z(svadd_f64_z(svmul_f64_z(
                         COEFF_P5_X, fractional_part), COEFF_P5_Y), fractional_part),
                         COEFF_P5_Z),fractional_part), COEFF_P5_A), fractional_part),
                         COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part),
                         COEFF_P5_F);

        x = svsub_f64_z(x,factor);

        x = svadd_f64_z(svmul_f64_z(COEFF_A, x), COEFF_B);

        __m128d valupper = _mm256_extractf128_pd(x, 1);
        __m128d vallower = _mm256_castpd256_pd128(x);

        alignas(64) long long int allvalint[VecLength] = { _mm_cvtsd_si64(vallower),
                                                           _mm_cvtsd_si64(_mm_shuffle_pd(vallower, vallower, 1)),
                                                           _mm_cvtsd_si64(valupper),
                                                           _mm_cvtsd_si64(_mm_shuffle_pd(valupper, valupper, 1)) };

        return _mm256_castsi256_pd(_mm256_load_si256(reinterpret_cast<const svbool_t*>(allvalint)));
#endif
    }

    inline InaVecSVE expLowAcc() const {
        const svfloat64_t COEFF_LOG2E = svdup_f64(double(InaFastExp::CoeffLog2E()));
        const svfloat64_t COEFF_A     = svdup_f64(double(InaFastExp::CoeffA64()));
        const svfloat64_t COEFF_B     = svdup_f64(double(InaFastExp::CoeffB64()));
        const svfloat64_t COEFF_P5_C  = svdup_f64(double(InaFastExp::GetCoefficient4_3()));
        const svfloat64_t COEFF_P5_D  = svdup_f64(double(InaFastExp::GetCoefficient4_2()));
        const svfloat64_t COEFF_P5_E  = svdup_f64(double(InaFastExp::GetCoefficient4_1()));
        const svfloat64_t COEFF_P5_F  = svdup_f64(double(InaFastExp::GetCoefficient4_0()));

        svfloat64_t x = svmul_f64_z(vec, COEFF_LOG2E);

        const svfloat64_t fractional_part = svsub_f64_z(x, InaVecSVE(x).floor().vec);

        svfloat64_t factor = svadd_f64_z(svmul_f64_z(svadd_f64_z(
                         svmul_f64_z(svadd_f64_z(svmul_f64_z(
                                         COEFF_P5_C, fractional_part),
                                         COEFF_P5_D), fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = svsub_f64_z(x,factor);

        x = svadd_f64_z(svmul_f64_z(COEFF_A, x), COEFF_B);

        __m128d valupper = _mm256_extractf128_pd(x, 1);
        __m128d vallower = _mm256_castpd256_pd128(x);

        alignas(64) long long int allvalint[VecLength] = { _mm_cvtsd_si64(vallower),
                                                           _mm_cvtsd_si64(_mm_shuffle_pd(vallower, vallower, 1)),
                                                           _mm_cvtsd_si64(valupper),
                                                           _mm_cvtsd_si64(_mm_shuffle_pd(valupper, valupper, 1)) };

        return _mm256_castsi256_pd(_mm256_load_si256(reinterpret_cast<const svbool_t*>(allvalint)));
    }

    inline InaVecSVE rsqrt() const {
        return svdup_f64(1) / svsqrt_f64_z(svptrue_b64(),vec);
    }

    inline InaVecSVE abs() const {
        const svfloat64_t minus0 = _mm256_castsi256_pd(_mm256_set1_epi64x(static_cast<long long>(0x8000000000000000L)));
        return _mm256_andnot_pd(minus0, vec);
    }

    inline InaVecSVE floor() const {
        return _mm256_floor_pd(vec);
    }

    inline InaVecSVE signOf() const {
        const svfloat64_t minus0 = _mm256_castsi256_pd(_mm256_set1_epi64x(static_cast<long long>(0x8000000000000000L)));
        const svfloat64_t signs  = _mm256_and_pd(vec, minus0);
        return _mm256_and_pd(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_NEQ_OQ),
                             _mm256_or_pd(signs, svdup_f64(1)));
    }

    inline InaVecSVE isPositive() const {
        const svfloat64_t greater = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_LE_OQ);
        const svfloat64_t ones    = svdup_f64(1);
        return _mm256_and_pd(greater, ones);
    }

    inline InaVecSVE isNegative() const {
        const svfloat64_t less = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_GE_OQ);
        const svfloat64_t ones = svdup_f64(1);
        return _mm256_and_pd(less, ones);
    }

    inline InaVecSVE isPositiveStrict() const {
        const svfloat64_t greater = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_LT_OQ);
        const svfloat64_t ones    = svdup_f64(1);
        return _mm256_and_pd(greater, ones);
    }

    inline InaVecSVE isNegativeStrict() const {
        const svfloat64_t less = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_GT_OQ);
        const svfloat64_t ones = svdup_f64(1);
        return _mm256_and_pd(less, ones);
    }

    inline InaVecSVE isZero() const {
        const svfloat64_t equalZero = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_EQ_OQ);
        const svfloat64_t ones      = svdup_f64(1);
        return _mm256_and_pd(equalZero, ones);
    }

    inline InaVecSVE isNotZero() const {
        const svfloat64_t equalZero = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_NEQ_OQ);
        const svfloat64_t ones      = svdup_f64(1);
        return _mm256_and_pd(equalZero, ones);
    }

    inline InaVecMaskSVE<double> isPositiveMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_LE_OQ));
    }

    inline InaVecMaskSVE<double> isNegativeMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_GE_OQ));
    }

    inline InaVecMaskSVE<double> isPositiveStrictMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_LT_OQ));
    }

    inline InaVecMaskSVE<double> isNegativeStrictMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_GT_OQ));
    }

    inline InaVecMaskSVE<double> isZeroMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_EQ_OQ));
    }

    inline InaVecMaskSVE<double> isNotZeroMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_NEQ_OQ));
    }

    // Static basic methods
    inline static InaVecSVE GetZero() {
        return InaVecSVE(_mm256_setzero_pd());
    }

    inline static InaVecSVE GetOne() {
        return InaVecSVE(svdup_f64(1));
    }

    inline static InaVecSVE Min(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_min_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE Max(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_max_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE IsLowerOrEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat64_t testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_LE_OQ);
        const svfloat64_t ones       = svdup_f64(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecSVE IsLower(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat64_t testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_LT_OQ);
        const svfloat64_t ones       = svdup_f64(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecSVE IsGreaterOrEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat64_t testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_GE_OQ);
        const svfloat64_t ones       = svdup_f64(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecSVE IsGreater(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat64_t testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_GT_OQ);
        const svfloat64_t ones       = svdup_f64(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecSVE IsEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat64_t testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_EQ_OQ);
        const svfloat64_t ones       = svdup_f64(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecSVE IsNotEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        const svfloat64_t testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ);
        const svfloat64_t ones       = svdup_f64(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecMaskSVE<double> IsLowerOrEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_LE_OQ));
    }

    inline static InaVecMaskSVE<double> IsLowerMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_LT_OQ));
    }

    inline static InaVecMaskSVE<double> IsGreaterOrEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_GE_OQ));
    }

    inline static InaVecMaskSVE<double> IsGreaterMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_GT_OQ));
    }

    inline static InaVecMaskSVE<double> IsEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_EQ_OQ));
    }

    inline static InaVecMaskSVE<double> IsNotEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ));
    }

    inline static InaVecSVE BitsAnd(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_and_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE BitsNotAnd(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_andnot_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE BitsOr(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_or_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE BitsXor(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return _mm256_xor_pd(inVec1.vec, inVec2.vec);
    }

    inline static  const char* GetName() {
        return "InaVecSVE<double>";
    }

    inline static  InaIfElse< InaVecSVE<double> >::ThenClass If(const InaVecMaskSVE<double>& inTest) {
        return InaIfElse< InaVecSVE<double> >::IfClass().If(inTest);
    }

    inline static InaVecSVE IfElse(const InaVecMaskSVE<double>& inMask, const InaVecSVE& inIfTrue, const InaVecSVE& inIfFalse) {
        return _mm256_or_pd(IfTrue(inMask, inIfTrue.vec).vec,
                      IfFalse(inMask, inIfFalse.vec).vec);
    }

    inline static InaVecSVE IfTrue(const InaVecMaskSVE<double>& inMask, const InaVecSVE& inIfTrue) {
        return _mm256_and_pd(_mm256_castsi256_pd(inMask.getMask()), inIfTrue.vec);
    }

    inline static InaVecSVE IfFalse(const InaVecMaskSVE<double>& inMask, const InaVecSVE& inIfFalse) {
        return _mm256_andnot_pd(_mm256_castsi256_pd(inMask.getMask()), inIfFalse.vec);
    }

    // Inner operators
    inline InaVecSVE<double>& operator+=(const InaVecSVE<double>& inVec){
        vec = svadd_f64_z(vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<double>& operator-=(const InaVecSVE<double>& inVec){
        vec = svsub_f64_z(vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<double>& operator/=(const InaVecSVE<double>& inVec){
        vec = svdiv_f64_z(vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<double>& operator*=(const InaVecSVE<double>& inVec){
        vec = svmul_f64_z(vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<double> operator-() const {
        const svfloat64_t minus0 = _mm256_castsi256_pd(_mm256_set1_epi64x(static_cast<long long>(0x8000000000000000L)));
        return _mm256_xor_pd(vec, minus0);
    }

    inline InaVecSVE<double> pow(std::size_t power) const{
        return InaUtils::FastPow<InaVecSVE<double>>(*this, power);
    }

    // Multiple sum
    template <class ... Args>
    inline static void MultiHorizontalSum(double sumRes[], const InaVecSVE<double>& inVec1,
                                          const InaVecSVE<double>& inVec2, const InaVecSVE<double>& inVec3,
                                          const InaVecSVE<double>& inVec4, Args ...args){
        const svfloat64_t val_a01_b01_a23_b23 = _mm256_hadd_pd(inVec1.vec, inVec2.vec);
        const svfloat64_t val_c01_d01_c23_d23 = _mm256_hadd_pd(inVec3.vec, inVec4.vec);

        const svfloat64_t val_a01_b01_c23_d23 = _mm256_permute2f128_pd(val_a01_b01_a23_b23, val_c01_d01_c23_d23, 0x30);// 011.0000
        const svfloat64_t val_a23_b23_d01_d23 = _mm256_permute2f128_pd(val_a01_b01_a23_b23, val_c01_d01_c23_d23, 0x21);// 010.0001

        const svfloat64_t val_suma_sumb_sumc_sumd = val_a01_b01_c23_d23 + val_a23_b23_d01_d23;

        svfloat64_t vecBuffer = svld1_f64_z(svptrue_b64(),sumRes);
        vecBuffer += val_suma_sumb_sumc_sumd;
        svst1_f64(svptrue_b64(),sumRes, vecBuffer);

        MultiHorizontalSum(&sumRes[4], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(double sumRes[], const InaVecSVE<double>& inVec1,
                                          const InaVecSVE<double>& inVec2, Args ...args){
        const svfloat64_t val_a01_b01_a23_b23 = _mm256_hadd_pd(inVec1.vec, inVec2.vec);

        __m128d valupper = _mm256_extractf128_pd(val_a01_b01_a23_b23, 1);
        __m128d vallower = _mm256_castpd256_pd128(val_a01_b01_a23_b23);

        __m128d vecBuffer = _mm_loadu_pd(sumRes);
        vecBuffer += valupper+vallower;
        _mm_storeu_pd(sumRes, vecBuffer);

        MultiHorizontalSum(&sumRes[2], args... );
    }

    inline static void MultiHorizontalSum(double sumRes[], const InaVecSVE<double>& inVec){
        sumRes[0] += inVec.horizontalSum();
    }

    inline static void MultiHorizontalSum(double /*sumRes*/[]){
    }
};

// Bits operators
inline InaVecSVE<double> operator&(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return InaVecSVE<double>::BitsAnd(inVec1, inVec2);
}

inline InaVecSVE<double> operator|(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return InaVecSVE<double>::BitsOr(inVec1, inVec2);
}

inline InaVecSVE<double> operator^(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return InaVecSVE<double>::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecSVE<double> operator+(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return svadd_f64_z(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSVE<double> operator-(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return svsub_f64_z(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSVE<double> operator/(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return svdiv_f64_z(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSVE<double> operator*(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return svmul_f64_z(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskSVE<double> operator<(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return InaVecSVE<double>::IsLowerMask(inVec1,inVec2);
}

inline InaVecMaskSVE<double> operator<=(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return InaVecSVE<double>::IsLowerOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskSVE<double> operator>(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return InaVecSVE<double>::IsGreaterMask(inVec1,inVec2);
}

inline InaVecMaskSVE<double> operator>=(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return InaVecSVE<double>::IsGreaterOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskSVE<double> operator==(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return InaVecSVE<double>::IsEqualMask(inVec1,inVec2);
}

inline InaVecMaskSVE<double> operator!=(const InaVecSVE<double>& inVec1, const InaVecSVE<double>& inVec2){
    return InaVecSVE<double>::IsNotEqualMask(inVec1,inVec2);
}


#endif
