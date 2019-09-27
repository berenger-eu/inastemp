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
        vec = svld1_gather_f64_index_z(svptrue_b64(), values, svld1_s64_u64(svptrue_b64(),inIndirection));
        return *this;
    }

    inline InaVecSVE& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        vec = svld1_gather_f64_index_z(svptrue_b64(), inArray,
                                       svld1_s64_u64(svptrue_b64(),inIndirection1)*svdup_s64(inLeadingDimension)
                                       +svld1_s64_u64(svptrue_b64(),inIndirection2));
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

        svint32_t castedInteger = svcvt_f64_s64_z(x);

        return svreinterpret_s64_f64(castedInteger);
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

        svint32_t castedInteger = svcvt_f64_s64_z(x);

        return svreinterpret_s64_f64(castedInteger);
    }

    inline InaVecSVE rsqrt() const {
        return svrsqrte_f64(svptrue_b64(),vec);
    }

    inline InaVecSVE abs() const {
      return svabs_f64_z(vec, svptrue_b64());
    }

    inline InaVecSVE floor() const {
        svbool_t mask = (InaVecSVE(double(LLONG_MIN)) < vec && vec < InaVecSVE(double(LLONG_MAX)));
        svint64_t n = svcvt_s64_f64_z(svptrue_b64(), vec);
        svint64_t d = svcvt_f64_s64_z(svptrue_b64(), n);
        svbool_t lower = svacgt_f64(svptrue_b64(), svdup_s64(0), d);
        return svsel_f64(mask, svsel_f64(Lower, vec - svdup_f64(1), d), vec);
    }

    inline InaVecSVE signOf() const {
        return svsel_f64(svacgt_f64(svptrue_b64(), vec, svdup_f64(0)),
                  svdup_f64(1), svsel_f64(svacgt_f64(svptrue_b64(), svdup_f64(0), vec),
                                          svdup_f64(-1), svdup_f64(0)));
    }

    inline InaVecSVE isPositive() const {
        return svsel_f64(svacge_f64(svptrue_b64(), vec, svdup_f64(0)),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecSVE isNegative() const {
        return svsel_f64(svacgt_f64(svptrue_b64(), svdup_f64(0), vec),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecSVE isPositiveStrict() const {
        return svsel_f64(svacgt_f64(svptrue_b64(), vec, svdup_f64(0)),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecSVE isNegativeStrict() const {
        return svsel_f64(svacge_f64(svptrue_b64(), svdup_f64(0), vec),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecSVE isZero() const {
        return svsel_f64(svcmpeq_f64(svptrue_b64(), svdup_f64(0), vec),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecSVE isNotZero() const {
        return svsel_f64(svcmpne_f64(svptrue_b64(), svdup_f64(0), vec),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecMaskSVE<double> isPositiveMask() const {
        return svacge_f64(svptrue_b64(), vec, svdup_f64(0));
    }

    inline InaVecMaskSVE<double> isNegativeMask() const {
        return svacgt_f64(svptrue_b64(), svdup_f64(0), vec);
    }

    inline InaVecMaskSVE<double> isPositiveStrictMask() const {
        return svacgt_f64(svptrue_b64(), vec, svdup_f64(0));
    }

    inline InaVecMaskSVE<double> isNegativeStrictMask() const {
        return svacge_f64(svptrue_b64(), svdup_f64(0), vec);
    }

    inline InaVecMaskSVE<double> isZeroMask() const {
        return svcmpeq_f64(svptrue_b64(), svdup_f64(0), vec);
    }

    inline InaVecMaskSVE<double> isNotZeroMask() const {
        return svcmpne_f64(svptrue_b64(), svdup_f64(0), vec);
    }

    // Static basic methods
    inline static InaVecSVE GetZero() {
        return InaVecSVE(svdup_f64(0));
    }

    inline static InaVecSVE GetOne() {
        return InaVecSVE(svdup_f64(1));
    }

    inline static InaVecSVE Min(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svminp_f64_z(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE Max(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svmaxp_f64_z(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE IsLowerOrEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f64_z(svacgt_f64(svptrue_b64(), inVec2.vec, inVec1.vec), svdup_f64(1));
    }

    inline static InaVecSVE IsLower(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f64_z(svacge_f64(svptrue_b64(), inVec2.vec, inVec1.vec), svdup_f64(1));
    }

    inline static InaVecSVE IsGreaterOrEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f64_z(svacge_f64(svptrue_b64(), inVec1.vec, inVec2.vec), svdup_f64(1));
    }

    inline static InaVecSVE IsGreater(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f64_z(svacgt_f64(svptrue_b64(), inVec1.vec, inVec2.vec), svdup_f64(1));
    }

    inline static InaVecSVE IsEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f64_z(svcmpeq_f64(svptrue_b64(), inVec1.vec, inVec2.vec), svdup_f64(1));
    }

    inline static InaVecSVE IsNotEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f64_z(svcmpne_f64(svptrue_b64(), inVec1.vec, inVec2.vec), svdup_f64(1));
    }

    inline static InaVecMaskSVE<double> IsLowerOrEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svacgt_f64(svptrue_b64(), inVec2.vec, inVec1.vec);
    }

    inline static InaVecMaskSVE<double> IsLowerMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svacge_f64(svptrue_b64(), inVec2.vec, inVec1.vec);
    }

    inline static InaVecMaskSVE<double> IsGreaterOrEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svacge_f64(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSVE<double> IsGreaterMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svacgt_f64(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSVE<double> IsEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svcmpeq_f64(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSVE<double> IsNotEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svcmpne_f64(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE BitsAnd(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svreinterpret_f64_s64(svand_b64(svptrue_b64(), svreinterpret_s64_f64(inVec1.vec), svreinterpret_s64_f64(inVec2.vec)));
    }

    inline static InaVecSVE BitsNotAnd(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svreinterpret_f64_s64(svbic_b64(svptrue_b64(), svreinterpret_s64_f64(inVec1.vec), svreinterpret_s64_f64(inVec2.vec)));
    }

    inline static InaVecSVE BitsOr(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svreinterpret_f64_s64(svorr_b64(svptrue_b64(), svreinterpret_s64_f64(inVec1.vec), svreinterpret_s64_f64(inVec2.vec)));
    }

    inline static InaVecSVE BitsXor(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svreinterpret_f64_s64(sveor_b64(svptrue_b64(), svreinterpret_s64_f64(inVec1.vec), svreinterpret_s64_f64(inVec2.vec)));
    }

    inline static  const char* GetName() {
        return "InaVecSVE<double>";
    }

    inline static  InaIfElse< InaVecSVE<double> >::ThenClass If(const InaVecMaskSVE<double>& inTest) {
        return InaIfElse< InaVecSVE<double> >::IfClass().If(inTest);
    }

    inline static InaVecSVE IfElse(const InaVecMaskSVE<double>& inMask, const InaVecSVE& inIfTrue, const InaVecSVE& inIfFalse) {
        return svsel_f64(inMask, inIfTrue.vec, inIfFalse.vec);
    }

    inline static InaVecSVE IfTrue(const InaVecMaskSVE<double>& inMask, const InaVecSVE& inIfTrue) {
        return svsel_f64(inMask, inIfTrue.vec, svdup_f64(0));
    }

    inline static InaVecSVE IfFalse(const InaVecMaskSVE<double>& inMask, const InaVecSVE& inIfFalse) {
        return svsel_f64(inMask, svdup_f64(0), inIfFalse.vec);
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
        return svneg_f64(svptrue_b64(), vec);
    }

    inline InaVecSVE<double> pow(std::size_t power) const{
        return InaUtils::FastPow<InaVecSVE<double>>(*this, power);
    }

    // Multiple sum
    template <class ... Args>
    inline static void MultiHorizontalSum(double sumRes[], const InaVecSVE<double>& inVec1,
                                          const InaVecSVE<double>& inVec2, const InaVecSVE<double>& inVec3,
                                          const InaVecSVE<double>& inVec4, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1, inVec2 );
        MultiHorizontalSum(&sumRes[2], inVec3, inVec4 );

        MultiHorizontalSum(&sumRes[4], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(double sumRes[], const InaVecSVE<double>& inVec1,
                                          const InaVecSVE<double>& inVec2, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1);
        MultiHorizontalSum(&sumRes[1], inVec2);

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
