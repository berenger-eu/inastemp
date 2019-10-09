///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSXAFLOAT_HPP
#define INAVECSXAFLOAT_HPP

#include "InastempGlobal.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_SXA
#error InaVecSXA<float> is included but SXA is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <velintrin.h>
#include <cmath>
#include <initializer_list>
#include <limits>

// Forward declarations
template <class RealType>
class InaVecMaskSXA;

template <class RealType>
class InaVecSXA;


// Mask type
template <>
class InaVecMaskSXA<float> {
    __vr mask;

public:
    // Classic constructors
    inline InaVecMaskSXA(){ _ve_lvl(256); }

    inline InaVecMaskSXA(const InaVecMaskSXA& inMask){
        _ve_lvl(256);
        mask = inMask.mask;
    }

    inline InaVecMaskSXA& operator=(const InaVecMaskSXA& inMask){
        mask = inMask.mask;
        return *this;
    }

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskSXA(const __vr inMask)
        : InaVecMaskSXA() {
        mask = (inMask);
    }

    inline InaVecMaskSXA& operator=(const __vr inMask){
        mask = inMask;
        return (*this);
    }

    inline explicit operator __vr() const{
        return mask;
    }

    inline __vr getMask() const{
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskSXA(const bool inBool) : InaVecMaskSXA() {
        mask = (inBool? svptrue_b32() : svpfalse());
    }

    inline InaVecMaskSXA& operator=(const bool inBool){
        mask = (inBool? svptrue_b32() : svpfalse());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskSXA Not() const{
        return svnot_z(svptrue_b32(), mask);
    }

    inline bool isAllTrue() const{
        // Could with svptest_any(svptrue_b32(), pg)
        return svcntp_b32(svptrue_b32(), mask) == svcntd();
    }

    inline bool isAllFalse() const{
        // true if all zero
        return svcntp_b32(svptrue_b32(), mask) == 0;
    }

    // Double args methods
    inline static InaVecMaskSXA And(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return svand_z(svptrue_b32(), inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSXA NotAnd(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return svbic_z(svptrue_b32(), inMask2.mask,inMask1.mask);
    }

    inline static InaVecMaskSXA Or(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return svorr_z(svptrue_b32(), inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSXA Xor(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return sveor_z(svptrue_b32(), inMask1.mask,inMask2.mask);
    }

    inline static bool IsEqual(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return svcntp_b32(svptrue_b32(), sveor_z(svptrue_b32(), inMask1.mask,inMask2.mask)) == 0;
    }

    inline static bool IsNotEqual(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return svcntp_b32(svptrue_b32(), sveor_z(svptrue_b32(), inMask1.mask,inMask2.mask)) != 0;
    }
};

// Mask must have operators
inline InaVecMaskSXA<float> operator&(const InaVecMaskSXA<float>& inMask1, const InaVecMaskSXA<float>& inMask2){
    return InaVecMaskSXA<float>::And(inMask1, inMask2);
}

inline InaVecMaskSXA<float> operator|(const InaVecMaskSXA<float>& inMask1, const InaVecMaskSXA<float>& inMask2){
    return InaVecMaskSXA<float>::Or(inMask1, inMask2);
}

inline InaVecMaskSXA<float> operator^(const InaVecMaskSXA<float>& inMask1, const InaVecMaskSXA<float>& inMask2){
    return InaVecMaskSXA<float>::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskSXA<float>& inMask1, const InaVecMaskSXA<float>& inMask2){
    return InaVecMaskSXA<float>::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskSXA<float>& inMask1, const InaVecMaskSXA<float>& inMask2){
    return InaVecMaskSXA<float>::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class InaVecSXA<float> {
protected:
    __vr vec;

public:
    using VecRawType           = __vr;
    using MaskType             = InaVecMaskSXA<float>;
    using RealType             = float;
    static const int Alignement= 1;
    static const bool IsOfFixedSize = false;

    static constexpr int GetVecLength(){
        return 256;
    }

    inline InaVecSXA() { _ve_lvl(256);  }
    inline InaVecSXA(const InaVecSXA& inVec){
        _ve_lvl(256);
        vec = inVec.vec;
    }

    inline InaVecSXA& operator=(const InaVecSXA& inVec){
        vec = inVec.vec;
        return *this;
    }

    // Constructor from raw type
    inline /*not explicit*/ InaVecSXA(const __vr inVec)
        : InaVecSXA() {
        vec = (inVec);
    }

    inline InaVecSXA& operator=(const __vr inVec){
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __vr inVec){
        vec = inVec;
    }

    inline explicit operator __vr() const{
        return vec;
    }

    inline __vr getVec() const{
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecSXA(const float val)
        : InaVecSXA() {
        vec = (svdup_f32(val));
    }

    inline InaVecSXA& operator=(const float val){
        vec = svdup_f32(val);
        return *this;
    }

    inline void setFromScalar(const float val){
        vec = svdup_f32(val);
    }

    // Constructor from vec
    inline InaVecSXA(const std::initializer_list<float> lst)
        : InaVecSXA(lst.begin()){
    }

    inline explicit InaVecSXA(const float ptr[])
        : InaVecSXA() {
        vec = (svld1_f32(svptrue_b32(),ptr));
    }

    inline InaVecSXA& setFromArray(const float ptr[]){
        vec = svld1_f32(svptrue_b32(),ptr);
        return *this;
    }

    inline InaVecSXA& setFromAlignedArray(const float ptr[]){
        vec = svld1_f32(svptrue_b32(),ptr);
        return *this;
    }

    inline InaVecSXA& setFromIndirectArray(const float values[], const int inIndirection[]) {
        vec = svld1_gather_s32index_f32(svptrue_b32(), values, svld1_s32(svptrue_b32(),inIndirection));
        return *this;
    }

    inline InaVecSXA& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        vec = svld1_gather_s32index_f32(svptrue_b32(), inArray,
                                       svadd_s32_z(svptrue_b32(),svmul_s32_z(svptrue_b32(),svld1_s32(svptrue_b32(),inIndirection1),svdup_s32(inLeadingDimension)),
                                       svld1_s32(svptrue_b32(),inIndirection2)));
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
        return svlasta_f32(svwhilelt_b32(0, index),vec);
    }

    // Horizontal operation
    inline float horizontalSum() const {
        return svadda_f32(svptrue_b32(), 0, vec);
    }

    inline float horizontalMul() const {
        float sum = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            sum *= at(idx);
        }
        return sum;
    }

    inline InaVecSXA sqrt() const {
        return svsqrt_f32_z(svptrue_b32(),vec);
    }

    inline InaVecSXA exp() const {
        const __vr COEFF_LOG2E = svdup_f32(float(InaFastExp::CoeffLog2E()));
        const __vr COEFF_A     = svdup_f32(float(InaFastExp::CoeffA32()));
        const __vr COEFF_B     = svdup_f32(float(InaFastExp::CoeffB32()));
        const __vr COEFF_P5_A  = svdup_f32(float(InaFastExp::GetCoefficient6_5()));
        const __vr COEFF_P5_B  = svdup_f32(float(InaFastExp::GetCoefficient6_4()));
        const __vr COEFF_P5_C  = svdup_f32(float(InaFastExp::GetCoefficient6_3()));
        const __vr COEFF_P5_D  = svdup_f32(float(InaFastExp::GetCoefficient6_2()));
        const __vr COEFF_P5_E  = svdup_f32(float(InaFastExp::GetCoefficient6_1()));
        const __vr COEFF_P5_F  = svdup_f32(float(InaFastExp::GetCoefficient6_0()));

        __vr x = svmul_f32_z(svptrue_b32(), vec, COEFF_LOG2E);

        const __vr fractional_part = svsub_f32_z(svptrue_b32(), x, InaVecSXA(x).floor().vec);

        __vr factor = svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(),
                         svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(),svmul_f32_z(svptrue_b32(),
                         COEFF_P5_A, fractional_part), COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part), COEFF_P5_F);

        x = svsub_f32_z(svptrue_b32(), x,factor);

        svint32_t castedInteger = svcvt_s32_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), COEFF_A, x), COEFF_B));

        return svreinterpret_f32_s32(castedInteger);
    }

    inline InaVecSXA expLowAcc() const {
        const __vr COEFF_LOG2E = svdup_f32(float(InaFastExp::CoeffLog2E()));
        const __vr COEFF_A     = svdup_f32(float(InaFastExp::CoeffA32()));
        const __vr COEFF_B     = svdup_f32(float(InaFastExp::CoeffB32()));
        const __vr COEFF_P5_D  = svdup_f32(float(InaFastExp::GetCoefficient3_2()));
        const __vr COEFF_P5_E  = svdup_f32(float(InaFastExp::GetCoefficient3_1()));
        const __vr COEFF_P5_F  = svdup_f32(float(InaFastExp::GetCoefficient3_0()));

        __vr x = svmul_f32_z(svptrue_b32(), vec, COEFF_LOG2E);

        const __vr fractional_part = svsub_f32_z(svptrue_b32(), x, InaVecSXA(x).floor().vec);

        __vr factor = svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),
                         svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),
                                         COEFF_P5_D, fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = svsub_f32_z(svptrue_b32(), x,factor);

        svint32_t castedInteger = svcvt_s32_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), COEFF_A, x), COEFF_B));

        return svreinterpret_f32_s32(castedInteger);
    }

    inline InaVecSXA rsqrt() const {
        // too low acc svrsqrte_f32(vec);
        return  svdiv_f32_z(svptrue_b32(), svdup_f32(1), svsqrt_f32_z(svptrue_b32(),vec));
    }

    inline InaVecSXA abs() const {
        return svabs_f32_z(svptrue_b32(), vec);
    }

    inline InaVecSXA floor() const {
        __vr maskInLongInt = svand_z(svptrue_b32(),
                                svcmple_f32(svptrue_b32(), svdup_f32(float(std::numeric_limits<int>::min())), vec),
                                svcmple_f32(svptrue_b32(), vec, svdup_f32(float(std::numeric_limits<int>::max()))));
        svint32_t vecConvLongInt = svcvt_s32_f32_z(maskInLongInt, vec);
        __vr vecConvLongIntDouble = svcvt_f32_s32_z(maskInLongInt, vecConvLongInt);
        __vr maskTooLarge = svcmpgt_f32(svptrue_b32(), vecConvLongIntDouble, vec);
        return svsel_f32(maskInLongInt, svsel_f32(maskTooLarge,  svsub_f32_z(svptrue_b32(), vecConvLongIntDouble, svdup_f32(1)), vecConvLongIntDouble), vec);
    }

    inline InaVecSXA signOf() const {
        return svsel_f32(svcmplt_f32(svptrue_b32(), svdup_f32(0), vec),
                  svdup_f32(1), svsel_f32(svcmpgt_f32(svptrue_b32(), svdup_f32(0), vec),
                                          svdup_f32(-1), svdup_f32(0)));
    }

    inline InaVecSXA isPositive() const {
        return svsel_f32(svcmple_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecSXA isNegative() const {
        return svsel_f32(svcmpge_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecSXA isPositiveStrict() const {
        return svsel_f32(svcmplt_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecSXA isNegativeStrict() const {
        return svsel_f32(svcmpgt_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecSXA isZero() const {
        return svsel_f32(svcmpeq_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecSXA isNotZero() const {
        return svsel_f32(svcmpne_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecMaskSXA<float> isPositiveMask() const {
        return svorr_z(svptrue_b32(), svcmpeq_f32(svptrue_b32(), svdup_f32(0), vec),
                       svcmple_f32(svptrue_b32(), svdup_f32(0), vec));
    }

    inline InaVecMaskSXA<float> isNegativeMask() const {
        return svorr_z(svptrue_b32(), svcmpeq_f32(svptrue_b32(), svdup_f32(0), vec),
                       svcmpge_f32(svptrue_b32(), svdup_f32(0), vec));
    }

    inline InaVecMaskSXA<float> isPositiveStrictMask() const {
        return svcmplt_f32(svptrue_b32(), svdup_f32(0), vec);
    }

    inline InaVecMaskSXA<float> isNegativeStrictMask() const {
        return svcmpgt_f32(svptrue_b32(), svdup_f32(0), vec);
    }

    inline InaVecMaskSXA<float> isZeroMask() const {
        return svcmpeq_f32(svptrue_b32(), svdup_f32(0), vec);
    }

    inline InaVecMaskSXA<float> isNotZeroMask() const {
        return svcmpne_f32(svptrue_b32(), svdup_f32(0), vec);
    }

    // Static basic methods
    inline static InaVecSXA GetZero() {
        return InaVecSXA(svdup_f32(0));
    }

    inline static InaVecSXA GetOne() {
        return InaVecSXA(svdup_f32(1));
    }

    inline static InaVecSXA Min(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svmin_f32_z(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSXA Max(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svmax_f32_z(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSXA IsLowerOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f32_z(svcmple_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsLower(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f32_z(svcmplt_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsGreaterOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f32_z(svcmpge_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsGreater(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f32_z(svcmpgt_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f32_z(svcmpeq_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsNotEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f32_z(svcmpne_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecMaskSXA<float> IsLowerOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmple_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<float> IsLowerMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmplt_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<float> IsGreaterOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpge_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<float> IsGreaterMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpgt_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<float> IsEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpeq_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<float> IsNotEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpne_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSXA BitsAnd(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f32_s32(svand_s32_z(svptrue_b32(), svreinterpret_s32_f32(inVec1.vec), svreinterpret_s32_f32(inVec2.vec)));
    }

    inline static InaVecSXA BitsNotAnd(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f32_s32(svbic_s32_z(svptrue_b32(), svreinterpret_s32_f32(inVec2.vec), svreinterpret_s32_f32(inVec1.vec)));
    }

    inline static InaVecSXA BitsOr(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f32_s32(svorr_s32_z(svptrue_b32(), svreinterpret_s32_f32(inVec1.vec), svreinterpret_s32_f32(inVec2.vec)));
    }

    inline static InaVecSXA BitsXor(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f32_s32(sveor_s32_z(svptrue_b32(), svreinterpret_s32_f32(inVec1.vec), svreinterpret_s32_f32(inVec2.vec)));
    }

    inline static  const char* GetName() {
        return "InaVecSXA<float>";
    }

    inline static  InaIfElse< InaVecSXA<float> >::ThenClass If(const InaVecMaskSXA<float>& inTest) {
        return InaIfElse< InaVecSXA<float> >::IfClass().If(inTest);
    }

    inline static InaVecSXA IfElse(const InaVecMaskSXA<float>& inMask, const InaVecSXA& inIfTrue, const InaVecSXA& inIfFalse) {
        return svsel_f32(__vr(inMask), inIfTrue.vec, inIfFalse.vec);
    }

    inline static InaVecSXA IfTrue(const InaVecMaskSXA<float>& inMask, const InaVecSXA& inIfTrue) {
        return svsel_f32(__vr(inMask), inIfTrue.vec, svdup_f32(0));
    }

    inline static InaVecSXA IfFalse(const InaVecMaskSXA<float>& inMask, const InaVecSXA& inIfFalse) {
        return svsel_f32(__vr(inMask), svdup_f32(0), inIfFalse.vec);
    }

    // Inner operators
    inline InaVecSXA<float>& operator+=(const InaVecSXA<float>& inVec){
        vec = svadd_f32_z(svptrue_b32(), vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<float>& operator-=(const InaVecSXA<float>& inVec){
        vec = svsub_f32_z(svptrue_b32(), vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<float>& operator/=(const InaVecSXA<float>& inVec){
        vec = svdiv_f32_z(svptrue_b32(), vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<float>& operator*=(const InaVecSXA<float>& inVec){
        vec = svmul_f32_z(svptrue_b32(), vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<float> operator-() const {
        return svneg_f32_z(svptrue_b32(), vec);
    }

    inline InaVecSXA<float> pow(std::size_t power) const{
        return InaUtils::FastPow<InaVecSXA<float>>(*this, power);
    }

    // Multiple sum
    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecSXA<float>& inVec1,
                                          const InaVecSXA<float>& inVec2, const InaVecSXA<float>& inVec3,
                                          const InaVecSXA<float>& inVec4, const InaVecSXA<float>& inVec5,
                                          const InaVecSXA<float>& inVec6, const InaVecSXA<float>& inVec7,
                                          const InaVecSXA<float>& inVec8, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1, inVec2, inVec3, inVec4 );
        MultiHorizontalSum(&sumRes[4], inVec5, inVec6, inVec7, inVec8 );

        MultiHorizontalSum(&sumRes[8], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecSXA<float>& inVec1,
                                          const InaVecSXA<float>& inVec2, const InaVecSXA<float>& inVec3,
                                          const InaVecSXA<float>& inVec4, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1, inVec2 );
        MultiHorizontalSum(&sumRes[2], inVec3, inVec4 );

        MultiHorizontalSum(&sumRes[4], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecSXA<float>& inVec1,
                                          const InaVecSXA<float>& inVec2, Args ...args){

        MultiHorizontalSum(&sumRes[0], inVec1);
        MultiHorizontalSum(&sumRes[1], inVec2);

        MultiHorizontalSum(&sumRes[2], args... );
    }

    inline static void MultiHorizontalSum(float sumRes[], const InaVecSXA<float>& inVec){
        sumRes[0] += inVec.horizontalSum();
    }

    inline static void MultiHorizontalSum(float /*sumRes*/[]){
    }
};

// Bits operators
inline InaVecSXA<float> operator&(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return InaVecSXA<float>::BitsAnd(inVec1, inVec2);
}

inline InaVecSXA<float> operator|(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return InaVecSXA<float>::BitsOr(inVec1, inVec2);
}

inline InaVecSXA<float> operator^(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return InaVecSXA<float>::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecSXA<float> operator+(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return svadd_f32_z(svptrue_b32(), inVec1.getVec(), inVec2.getVec());
}

inline InaVecSXA<float> operator-(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return svsub_f32_z(svptrue_b32(), inVec1.getVec(), inVec2.getVec());
}

inline InaVecSXA<float> operator/(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return svdiv_f32_z(svptrue_b32(), inVec1.getVec(), inVec2.getVec());
}

inline InaVecSXA<float> operator*(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return svmul_f32_z(svptrue_b32(), inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskSXA<float> operator<(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return InaVecSXA<float>::IsLowerMask(inVec1,inVec2);
}

inline InaVecMaskSXA<float> operator<=(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return InaVecSXA<float>::IsLowerOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskSXA<float> operator>(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return InaVecSXA<float>::IsGreaterMask(inVec1,inVec2);
}

inline InaVecMaskSXA<float> operator>=(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return InaVecSXA<float>::IsGreaterOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskSXA<float> operator==(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return InaVecSXA<float>::IsEqualMask(inVec1,inVec2);
}

inline InaVecMaskSXA<float> operator!=(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return InaVecSXA<float>::IsNotEqualMask(inVec1,inVec2);
}


#endif
