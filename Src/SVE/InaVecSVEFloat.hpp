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
#include <limits>

// Forward declarations
template <class RealType>
class InaVecMaskSVE;

template <class RealType>
class InaVecSVE;


// Mask type
template <>
class InaVecMaskSVE<float> {
    // Ugly trick
    // Should simply be svbool_t mask;
    unsigned char maskData[2048/sizeof(unsigned char)];
    svbool_t& mask;

public:
    // Classic constructors
    inline InaVecMaskSVE() : mask(*reinterpret_cast<svbool_t*>(maskData)){}

    inline InaVecMaskSVE(const InaVecMaskSVE&) = default;

    inline InaVecMaskSVE& operator=(const InaVecMaskSVE& inMask){
        mask = inMask.mask;
        return *this;
    }

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskSVE(const svbool_t inMask)
        : InaVecMaskSVE() {
        mask = (inMask);
    }

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
    inline explicit InaVecMaskSVE(const bool inBool) : InaVecMaskSVE() {
        mask = (inBool? svptrue_b32() : svpfalse());
    }

    inline InaVecMaskSVE& operator=(const bool inBool){
        mask = (inBool? svptrue_b32() : svpfalse());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskSVE Not() const{
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
    inline static InaVecMaskSVE And(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return svand_z(svptrue_b32(), inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSVE NotAnd(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return svbic_z(svptrue_b32(), inMask2.mask,inMask1.mask);
    }

    inline static InaVecMaskSVE Or(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return svorr_z(svptrue_b32(), inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSVE Xor(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return sveor_z(svptrue_b32(), inMask1.mask,inMask2.mask);
    }

    inline static bool IsEqual(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return svcntp_b32(svptrue_b32(), sveor_z(svptrue_b32(), inMask1.mask,inMask2.mask)) == 0;
    }

    inline static bool IsNotEqual(const InaVecMaskSVE& inMask1, const InaVecMaskSVE& inMask2){
        return svcntp_b32(svptrue_b32(), sveor_z(svptrue_b32(), inMask1.mask,inMask2.mask)) != 0;
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
class InaVecSVE<float> {
protected:
    // Ugly trick
    // Should simply be svfloat32_t vec;
    unsigned char vecData[2048/sizeof(unsigned char)];
    svfloat32_t& vec;

public:
    using VecRawType           = svfloat32_t;
    using MaskType             = InaVecMaskSVE<float>;
    using RealType             = float;
    static const int Alignement= 1;
    static const bool IsOfFixedSize = false;

    static int GetVecLength(){
        return int(svcntw());
    }

    static constexpr bool IsRealFma(){
        return true;
    }

    inline InaVecSVE() : vec(*reinterpret_cast<svfloat32_t*>(vecData)) {}
    inline InaVecSVE(const InaVecSVE&) = default;

    inline InaVecSVE& operator=(const InaVecSVE& inVec){
        vec = inVec.vec;
        return *this;
    }

    // Constructor from raw type
    inline /*not explicit*/ InaVecSVE(const svfloat32_t inVec)
        : InaVecSVE() {
        vec = (inVec);
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
        : InaVecSVE() {
        vec = (svdup_f32(val));
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
        : InaVecSVE() {
        vec = (svld1_f32(svptrue_b32(),ptr));
    }

    inline InaVecSVE& setFromArray(const float ptr[]){
        vec = svld1_f32(svptrue_b32(),ptr);
        return *this;
    }

    inline InaVecSVE& setFromAlignedArray(const float ptr[]){
        vec = svld1_f32(svptrue_b32(),ptr);
        return *this;
    }

    inline InaVecSVE& setFromIndirectArray(const float values[], const int inIndirection[]) {
        vec = svld1_gather_s32index_f32(svptrue_b32(), values, svld1_s32(svptrue_b32(),inIndirection));
        return *this;
    }

    inline InaVecSVE& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
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

    inline float minInVec() const {
        return svminv_f32(svptrue_b64(), vec);
    }

    inline float maxInVec() const {
        return svmaxv_f32(svptrue_b64(), vec);
    }

    inline float horizontalMul() const {
        float sum = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            sum *= at(idx);
        }
        return sum;
    }

    inline InaVecSVE sqrt() const {
        return svsqrt_f32_z(svptrue_b32(),vec);
    }

    inline InaVecSVE exp() const {
        const svfloat32_t COEFF_LOG2E = svdup_f32(float(InaFastExp::CoeffLog2E()));
        const svfloat32_t COEFF_A     = svdup_f32(float(InaFastExp::CoeffA32()));
        const svfloat32_t COEFF_B     = svdup_f32(float(InaFastExp::CoeffB32()));
        const svfloat32_t COEFF_P5_A  = svdup_f32(float(InaFastExp::GetCoefficient6_5()));
        const svfloat32_t COEFF_P5_B  = svdup_f32(float(InaFastExp::GetCoefficient6_4()));
        const svfloat32_t COEFF_P5_C  = svdup_f32(float(InaFastExp::GetCoefficient6_3()));
        const svfloat32_t COEFF_P5_D  = svdup_f32(float(InaFastExp::GetCoefficient6_2()));
        const svfloat32_t COEFF_P5_E  = svdup_f32(float(InaFastExp::GetCoefficient6_1()));
        const svfloat32_t COEFF_P5_F  = svdup_f32(float(InaFastExp::GetCoefficient6_0()));

        svfloat32_t x = svmul_f32_z(svptrue_b32(), vec, COEFF_LOG2E);

        const svfloat32_t fractional_part = svsub_f32_z(svptrue_b32(), x, InaVecSVE(x).floor().vec);

        svfloat32_t factor = svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(),
                         svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(),svmul_f32_z(svptrue_b32(),
                         COEFF_P5_A, fractional_part), COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part), COEFF_P5_F);

        x = svsub_f32_z(svptrue_b32(), x,factor);

        svint32_t castedInteger = svcvt_s32_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), COEFF_A, x), COEFF_B));

        return svreinterpret_f32_s32(castedInteger);
    }

    inline InaVecSVE expLowAcc() const {
        const svfloat32_t COEFF_LOG2E = svdup_f32(float(InaFastExp::CoeffLog2E()));
        const svfloat32_t COEFF_A     = svdup_f32(float(InaFastExp::CoeffA32()));
        const svfloat32_t COEFF_B     = svdup_f32(float(InaFastExp::CoeffB32()));
        const svfloat32_t COEFF_P5_D  = svdup_f32(float(InaFastExp::GetCoefficient3_2()));
        const svfloat32_t COEFF_P5_E  = svdup_f32(float(InaFastExp::GetCoefficient3_1()));
        const svfloat32_t COEFF_P5_F  = svdup_f32(float(InaFastExp::GetCoefficient3_0()));

        svfloat32_t x = svmul_f32_z(svptrue_b32(), vec, COEFF_LOG2E);

        const svfloat32_t fractional_part = svsub_f32_z(svptrue_b32(), x, InaVecSVE(x).floor().vec);

        svfloat32_t factor = svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),
                         svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),
                                         COEFF_P5_D, fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = svsub_f32_z(svptrue_b32(), x,factor);

        svint32_t castedInteger = svcvt_s32_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), COEFF_A, x), COEFF_B));

        return svreinterpret_f32_s32(castedInteger);
    }

    inline InaVecSVE exp10() const {
        const svfloat32_t COEFF_LOG210 = svdup_f32(float(InaFastExp::CoeffLog210()));
        const svfloat32_t COEFF_A     = svdup_f32(float(InaFastExp::CoeffA32()));
        const svfloat32_t COEFF_B     = svdup_f32(float(InaFastExp::CoeffB32()));
        const svfloat32_t COEFF_P5_A  = svdup_f32(float(InaFastExp::GetCoefficient6_5()));
        const svfloat32_t COEFF_P5_B  = svdup_f32(float(InaFastExp::GetCoefficient6_4()));
        const svfloat32_t COEFF_P5_C  = svdup_f32(float(InaFastExp::GetCoefficient6_3()));
        const svfloat32_t COEFF_P5_D  = svdup_f32(float(InaFastExp::GetCoefficient6_2()));
        const svfloat32_t COEFF_P5_E  = svdup_f32(float(InaFastExp::GetCoefficient6_1()));
        const svfloat32_t COEFF_P5_F  = svdup_f32(float(InaFastExp::GetCoefficient6_0()));

        svfloat32_t x = svmul_f32_z(svptrue_b32(), vec, COEFF_LOG210);

        const svfloat32_t fractional_part = svsub_f32_z(svptrue_b32(), x, InaVecSVE(x).floor().vec);

        svfloat32_t factor = svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(),
                         svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(),svmul_f32_z(svptrue_b32(),
                         COEFF_P5_A, fractional_part), COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part), COEFF_P5_F);

        x = svsub_f32_z(svptrue_b32(), x,factor);

        svint32_t castedInteger = svcvt_s32_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), COEFF_A, x), COEFF_B));

        return svreinterpret_f32_s32(castedInteger);
    }

    inline InaVecSVE exp10LowAcc() const {
        const svfloat32_t COEFF_LOG210 = svdup_f32(float(InaFastExp::CoeffLog210()));
        const svfloat32_t COEFF_A     = svdup_f32(float(InaFastExp::CoeffA32()));
        const svfloat32_t COEFF_B     = svdup_f32(float(InaFastExp::CoeffB32()));
        const svfloat32_t COEFF_P5_D  = svdup_f32(float(InaFastExp::GetCoefficient3_2()));
        const svfloat32_t COEFF_P5_E  = svdup_f32(float(InaFastExp::GetCoefficient3_1()));
        const svfloat32_t COEFF_P5_F  = svdup_f32(float(InaFastExp::GetCoefficient3_0()));

        svfloat32_t x = svmul_f32_z(svptrue_b32(), vec, COEFF_LOG210);

        const svfloat32_t fractional_part = svsub_f32_z(svptrue_b32(), x, InaVecSVE(x).floor().vec);

        svfloat32_t factor = svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),
                         svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),
                                         COEFF_P5_D, fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = svsub_f32_z(svptrue_b32(), x,factor);

        svint32_t castedInteger = svcvt_s32_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), COEFF_A, x), COEFF_B));

        return svreinterpret_f32_s32(castedInteger);
    }

    inline InaVecSVE exp2() const {
        const svfloat32_t COEFF_A     = svdup_f32(float(InaFastExp::CoeffA32()));
        const svfloat32_t COEFF_B     = svdup_f32(float(InaFastExp::CoeffB32()));
        const svfloat32_t COEFF_P5_A  = svdup_f32(float(InaFastExp::GetCoefficient6_5()));
        const svfloat32_t COEFF_P5_B  = svdup_f32(float(InaFastExp::GetCoefficient6_4()));
        const svfloat32_t COEFF_P5_C  = svdup_f32(float(InaFastExp::GetCoefficient6_3()));
        const svfloat32_t COEFF_P5_D  = svdup_f32(float(InaFastExp::GetCoefficient6_2()));
        const svfloat32_t COEFF_P5_E  = svdup_f32(float(InaFastExp::GetCoefficient6_1()));
        const svfloat32_t COEFF_P5_F  = svdup_f32(float(InaFastExp::GetCoefficient6_0()));

        svfloat32_t x = vec;

        const svfloat32_t fractional_part = svsub_f32_z(svptrue_b32(), x, InaVecSVE(x).floor().vec);

        svfloat32_t factor = svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(),
                         svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),svadd_f32_z(svptrue_b32(),svmul_f32_z(svptrue_b32(),
                         COEFF_P5_A, fractional_part), COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part), COEFF_P5_F);

        x = svsub_f32_z(svptrue_b32(), x,factor);

        svint32_t castedInteger = svcvt_s32_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), COEFF_A, x), COEFF_B));

        return svreinterpret_f32_s32(castedInteger);
    }

    inline InaVecSVE exp2LowAcc() const {
        const svfloat32_t COEFF_A     = svdup_f32(float(InaFastExp::CoeffA32()));
        const svfloat32_t COEFF_B     = svdup_f32(float(InaFastExp::CoeffB32()));
        const svfloat32_t COEFF_P5_D  = svdup_f32(float(InaFastExp::GetCoefficient3_2()));
        const svfloat32_t COEFF_P5_E  = svdup_f32(float(InaFastExp::GetCoefficient3_1()));
        const svfloat32_t COEFF_P5_F  = svdup_f32(float(InaFastExp::GetCoefficient3_0()));

        svfloat32_t x = vec;

        const svfloat32_t fractional_part = svsub_f32_z(svptrue_b32(), x, InaVecSVE(x).floor().vec);

        svfloat32_t factor = svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),
                         svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(),
                                         COEFF_P5_D, fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = svsub_f32_z(svptrue_b32(), x,factor);

        svint32_t castedInteger = svcvt_s32_f32_z(svptrue_b32(), svadd_f32_z(svptrue_b32(), svmul_f32_z(svptrue_b32(), COEFF_A, x), COEFF_B));

        return svreinterpret_f32_s32(castedInteger);
    }

    inline InaVecSVE rsqrt() const {
        // too low acc svrsqrte_f32(vec);
        return  svdiv_f32_z(svptrue_b32(), svdup_f32(1), svsqrt_f32_z(svptrue_b32(),vec));
    }

    inline InaVecSVE abs() const {
        return svabs_f32_z(svptrue_b32(), vec);
    }

    inline InaVecSVE floor() const {
        svbool_t maskInLongInt = svand_z(svptrue_b32(),
                                svcmple_f32(svptrue_b32(), svdup_f32(float(std::numeric_limits<int>::min())), vec),
                                svcmple_f32(svptrue_b32(), vec, svdup_f32(float(std::numeric_limits<int>::max()))));
        svint32_t vecConvLongInt = svcvt_s32_f32_z(maskInLongInt, vec);
        svfloat32_t vecConvLongIntDouble = svcvt_f32_s32_z(maskInLongInt, vecConvLongInt);
        svbool_t maskTooLarge = svcmpgt_f32(svptrue_b32(), vecConvLongIntDouble, vec);
        return svsel_f32(maskInLongInt, svsel_f32(maskTooLarge,  svsub_f32_z(svptrue_b32(), vecConvLongIntDouble, svdup_f32(1)), vecConvLongIntDouble), vec);
    }

    inline InaVecSVE signOf() const {
        return svsel_f32(svcmplt_f32(svptrue_b32(), svdup_f32(0), vec),
                  svdup_f32(1), svsel_f32(svcmpgt_f32(svptrue_b32(), svdup_f32(0), vec),
                                          svdup_f32(-1), svdup_f32(0)));
    }

    inline InaVecSVE isPositive() const {
        return svsel_f32(svcmple_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecSVE isNegative() const {
        return svsel_f32(svcmpge_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecSVE isPositiveStrict() const {
        return svsel_f32(svcmplt_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecSVE isNegativeStrict() const {
        return svsel_f32(svcmpgt_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecSVE isZero() const {
        return svsel_f32(svcmpeq_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecSVE isNotZero() const {
        return svsel_f32(svcmpne_f32(svptrue_b32(), svdup_f32(0), vec),
                         svdup_f32(1), svdup_f32(0));
    }

    inline InaVecMaskSVE<float> isPositiveMask() const {
        return svorr_z(svptrue_b32(), svcmpeq_f32(svptrue_b32(), svdup_f32(0), vec),
                       svcmple_f32(svptrue_b32(), svdup_f32(0), vec));
    }

    inline InaVecMaskSVE<float> isNegativeMask() const {
        return svorr_z(svptrue_b32(), svcmpeq_f32(svptrue_b32(), svdup_f32(0), vec),
                       svcmpge_f32(svptrue_b32(), svdup_f32(0), vec));
    }

    inline InaVecMaskSVE<float> isPositiveStrictMask() const {
        return svcmplt_f32(svptrue_b32(), svdup_f32(0), vec);
    }

    inline InaVecMaskSVE<float> isNegativeStrictMask() const {
        return svcmpgt_f32(svptrue_b32(), svdup_f32(0), vec);
    }

    inline InaVecMaskSVE<float> isZeroMask() const {
        return svcmpeq_f32(svptrue_b32(), svdup_f32(0), vec);
    }

    inline InaVecMaskSVE<float> isNotZeroMask() const {
        return svcmpne_f32(svptrue_b32(), svdup_f32(0), vec);
    }

    // Static basic methods
    inline static InaVecSVE GetZero() {
        return InaVecSVE(svdup_f32(0));
    }

    inline static InaVecSVE GetOne() {
        return InaVecSVE(svdup_f32(1));
    }

    inline static InaVecSVE Min(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svmin_f32_z(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE Max(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svmax_f32_z(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE IsLowerOrEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f32_z(svcmple_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSVE IsLower(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f32_z(svcmplt_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSVE IsGreaterOrEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f32_z(svcmpge_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSVE IsGreater(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f32_z(svcmpgt_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSVE IsEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f32_z(svcmpeq_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSVE IsNotEqual(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svdup_f32_z(svcmpne_f32(svptrue_b32(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecMaskSVE<float> IsLowerOrEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svcmple_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSVE<float> IsLowerMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svcmplt_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSVE<float> IsGreaterOrEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svcmpge_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSVE<float> IsGreaterMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svcmpgt_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSVE<float> IsEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svcmpeq_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSVE<float> IsNotEqualMask(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svcmpne_f32(svptrue_b32(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSVE BitsAnd(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svreinterpret_f32_s32(svand_s32_z(svptrue_b32(), svreinterpret_s32_f32(inVec1.vec), svreinterpret_s32_f32(inVec2.vec)));
    }

    inline static InaVecSVE BitsNotAnd(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svreinterpret_f32_s32(svbic_s32_z(svptrue_b32(), svreinterpret_s32_f32(inVec2.vec), svreinterpret_s32_f32(inVec1.vec)));
    }

    inline static InaVecSVE BitsOr(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svreinterpret_f32_s32(svorr_s32_z(svptrue_b32(), svreinterpret_s32_f32(inVec1.vec), svreinterpret_s32_f32(inVec2.vec)));
    }

    inline static InaVecSVE BitsXor(const InaVecSVE& inVec1, const InaVecSVE& inVec2) {
        return svreinterpret_f32_s32(sveor_s32_z(svptrue_b32(), svreinterpret_s32_f32(inVec1.vec), svreinterpret_s32_f32(inVec2.vec)));
    }

    inline static  const char* GetName() {
        return "InaVecSVE<float>";
    }

    inline static  InaIfElse< InaVecSVE<float> >::ThenClass If(const InaVecMaskSVE<float>& inTest) {
        return InaIfElse< InaVecSVE<float> >::IfClass().If(inTest);
    }

    inline static InaVecSVE IfElse(const InaVecMaskSVE<float>& inMask, const InaVecSVE& inIfTrue, const InaVecSVE& inIfFalse) {
        return svsel_f32(svbool_t(inMask), inIfTrue.vec, inIfFalse.vec);
    }

    inline static InaVecSVE IfTrue(const InaVecMaskSVE<float>& inMask, const InaVecSVE& inIfTrue) {
        return svsel_f32(svbool_t(inMask), inIfTrue.vec, svdup_f32(0));
    }

    inline static InaVecSVE IfFalse(const InaVecMaskSVE<float>& inMask, const InaVecSVE& inIfFalse) {
        return svsel_f32(svbool_t(inMask), svdup_f32(0), inIfFalse.vec);
    }

    // Inner operators
    inline InaVecSVE<float>& operator+=(const InaVecSVE<float>& inVec){
        vec = svadd_f32_z(svptrue_b32(), vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<float>& operator-=(const InaVecSVE<float>& inVec){
        vec = svsub_f32_z(svptrue_b32(), vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<float>& operator/=(const InaVecSVE<float>& inVec){
        vec = svdiv_f32_z(svptrue_b32(), vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<float>& operator*=(const InaVecSVE<float>& inVec){
        vec = svmul_f32_z(svptrue_b32(), vec,inVec.vec);
        return *this;
    }

    inline InaVecSVE<float> operator-() const {
        return svneg_f32_z(svptrue_b32(), vec);
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
        MultiHorizontalSum(&sumRes[0], inVec1, inVec2, inVec3, inVec4 );
        MultiHorizontalSum(&sumRes[4], inVec5, inVec6, inVec7, inVec8 );

        MultiHorizontalSum(&sumRes[8], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecSVE<float>& inVec1,
                                          const InaVecSVE<float>& inVec2, const InaVecSVE<float>& inVec3,
                                          const InaVecSVE<float>& inVec4, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1, inVec2 );
        MultiHorizontalSum(&sumRes[2], inVec3, inVec4 );

        MultiHorizontalSum(&sumRes[4], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecSVE<float>& inVec1,
                                          const InaVecSVE<float>& inVec2, Args ...args){

        MultiHorizontalSum(&sumRes[0], inVec1);
        MultiHorizontalSum(&sumRes[1], inVec2);

        MultiHorizontalSum(&sumRes[2], args... );
    }

    inline static void MultiHorizontalSum(float sumRes[], const InaVecSVE<float>& inVec){
        sumRes[0] += inVec.horizontalSum();
    }

    inline static void MultiHorizontalSum(float /*sumRes*/[]){
    }

    inline static InaVecSVE<float> Fma(const InaVecSVE<float>& inValAdd, const InaVecSVE<float>& inValMul1, const InaVecSVE<float>& inValMul2){
        return svmad_f32_z(inValMul1.vec, inValMul2.vec, inValAdd.vec, svptrue_b32());
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
    return svadd_f32_z(svptrue_b32(), inVec1.getVec(), inVec2.getVec());
}

inline InaVecSVE<float> operator-(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return svsub_f32_z(svptrue_b32(), inVec1.getVec(), inVec2.getVec());
}

inline InaVecSVE<float> operator/(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return svdiv_f32_z(svptrue_b32(), inVec1.getVec(), inVec2.getVec());
}

inline InaVecSVE<float> operator*(const InaVecSVE<float>& inVec1, const InaVecSVE<float>& inVec2){
    return svmul_f32_z(svptrue_b32(), inVec1.getVec(), inVec2.getVec());
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
