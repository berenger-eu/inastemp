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
class alignas(2048) InaVecMaskSXA<float> {
    __vm256 mask;

public:
    // Classic constructors
    inline InaVecMaskSXA() {
        mask = _vel_xorm_mmm(mask,mask);
    }

    inline InaVecMaskSXA(const InaVecMaskSXA& inMask){
        mask = inMask.mask;
    }

    inline InaVecMaskSXA& operator=(const InaVecMaskSXA& inMask){
        mask = inMask.mask;
        return *this;
    }

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskSXA(const __vm256 inMask)
        : InaVecMaskSXA() {
        mask = (inMask);
    }

    inline InaVecMaskSXA& operator=(const __vm256 inMask){
        mask = inMask;
        return (*this);
    }

    inline explicit operator __vm256() const{
        return mask;
    }

    inline __vm256 getMask() const{
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskSXA(const bool inBool) : InaVecMaskSXA() {
        mask = (inBool? _vel_negm_mm(_vel_xorm_mmm(mask, mask)) : _vel_xorm_mmm(mask, mask));
    }

    inline InaVecMaskSXA& operator=(const bool inBool){
        mask = (inBool? _vel_negm_mm(_vel_xorm_mmm(mask, mask)) : _vel_xorm_mmm(mask, mask));
        return (*this);
    }

    // Binary methods
    inline InaVecMaskSXA Not() const{
        return _vel_negm_mm(mask);
    }

    inline bool isAllTrue() const{
        return _vel_pcvm_sml(mask, 256) == 256;
    }

    inline bool isAllFalse() const{
        // true if all zero
        return _vel_pcvm_sml(mask, 256) == 0;
    }

    // Float args methods
    inline static InaVecMaskSXA And(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_andm_mmm(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSXA NotAnd(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        __vm256 one = _vel_vfmklat_ml(0); // dumbe init
        one = _vel_negm_mm(_vel_xorm_mmm(one, one)); // set to zero than negate
        return _vel_andm_mmm(_vel_xorm_mmm(inMask1.mask, one),inMask2.mask);
    }

    inline static InaVecMaskSXA Or(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_orm_mmm(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSXA Xor(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_xorm_mmm(inMask1.mask,inMask2.mask);
    }

    inline static bool IsEqual(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_pcvm_sml(_vel_xorm_mmm(inMask1.mask,inMask2.mask), 256) == 0;
    }

    inline static bool IsNotEqual(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_pcvm_sml(_vel_xorm_mmm(inMask1.mask,inMask2.mask), 256) != 0;
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
class alignas(2048) InaVecSXA<float> {
protected:
    __vr vec;

    static __vm256 _vel_vfmklgt_mvl_256(__vr vz, int vl){
        return _vel_vfmklgt_mvl(vz, vl); // __vm256();
    }

public:
    using VecRawType           = __vr;
    using MaskType             = InaVecMaskSXA<float>;
    using RealType             = float;
    static const int Alignement= 1;
    static const bool IsOfFixedSize = true;

    static constexpr int GetVecLength(){
        return 256;
    }

    static constexpr bool IsRealFma(){
        return true;
    }

    inline InaVecSXA() {
        vec = _vel_vbrds_vsl(0,256);
    }
    inline InaVecSXA(const InaVecSXA& inVec){
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
        vec = _vel_vbrds_vsl(val,256);
    }

    inline InaVecSXA& operator=(const float val){
        vec = _vel_vbrds_vsl(val,256);
        return *this;
    }

    inline void setFromScalar(const float val){
        vec = _vel_vbrds_vsl(val,256);
    }

    // Constructor from vec
    inline InaVecSXA(const std::initializer_list<float> lst)
        : InaVecSXA(lst.begin()){
    }

    inline explicit InaVecSXA(const float ptr[])
        : InaVecSXA() {
        vec = _vel_vldu_vssl(4, ptr, 256);
    }

    inline InaVecSXA& setFromArray(const float ptr[]){
        vec = _vel_vldu_vssl(4, ptr, 256);
        return *this;
    }

    inline InaVecSXA& setFromAlignedArray(const float ptr[]){
        vec = _vel_vldu_vssl(4, ptr, 256);
        return *this;
    }

    inline InaVecSXA& setFromIndirectArray(const float values[], const unsigned long int inIndirection[]) {
        __vr offset = _vel_vld_vssl(8, inIndirection, 256);
        __vr address = _vel_vsfa_vvssl(offset, 2, reinterpret_cast<unsigned long>(values), 256);
        vec = _vel_vgtu_vvssl(address, 0, 0, 256);

        return *this;
    }

    inline InaVecSXA& setFromIndirectArray(const float values[], const int inIndirection[]) {
        __vr offset = _vel_vldu_vssl(4, inIndirection, 256);
        offset = _vel_vsrl_vvsl(offset, 32, 256);
        __vr address = _vel_vsfa_vvssl(offset, 2, reinterpret_cast<unsigned long>(values), 256);
        vec = _vel_vgtu_vvssl(address, 0, 0, 256);

        return *this;
    }

    inline InaVecSXA& setFromIndirect2DArray(const float inArray[], const long int inIndirection1[],
                                 const int inLeadingDimension, const long int inIndirection2[]){
        __vr offset = _vel_vaddsl_vvvl(_vel_vld_vssl(8, inIndirection2, 256),
                     _vel_vmulul_vvvl(_vel_vbrdl_vsl(inLeadingDimension, 256),
                                      _vel_vld_vssl(8, inIndirection1, 256),
                                      256),256);
        __vr address = _vel_vsfa_vvssl(offset, 2, reinterpret_cast<unsigned long>(inArray), 256);
        vec = _vel_vgtu_vvssl(address, 0, 0, 256);
        return *this;
    }

    inline InaVecSXA& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        __vr offset1 = _vel_vldu_vssl(4, inIndirection1, 256);
        offset1 = _vel_vsrl_vvsl(offset1, 32, 256);

        __vr offset2 = _vel_vldu_vssl(4, inIndirection2, 256);
        offset2 = _vel_vsrl_vvsl(offset2, 32, 256);

        __vr offset = _vel_vaddsl_vvvl(offset2,
                     _vel_vmulul_vvvl(_vel_vbrdl_vsl(inLeadingDimension, 256),
                                      offset1,
                                      256),256);
        __vr address = _vel_vsfa_vvssl(offset, 2, reinterpret_cast<unsigned long>(inArray), 256);
        vec = _vel_vgtu_vvssl(address, 0, 0, 256);

        return *this;
    }

    // Move back to array
    inline void storeInArray(float ptr[]) const {
        _vel_vstu_vssl(vec, 4, ptr, 256);
    }

    inline void storeInAlignedArray(float ptr[]) const {
        _vel_vstu_vssl(vec, 4, ptr, 256);
    }

    // Acce to individual values
    inline float at(const int index) const {
        return _vel_lvss_svs(vec, index);
    }

    // Horizontal operation
    inline float horizontalSum() const {
      return _vel_lvss_svs(_vel_vfsums_vvl(vec, 256),0);
    }

    inline float horizontalMul() const {
        // TODO use vfim when available
        float sum = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            sum *= at(idx);
        }
        return sum;
    }

    inline float minInVec() const {
        return _vel_lvss_svs(_vel_vfrminsfst_vvl(vec, 256),0);
    }

    inline float maxInVec() const {
        return _vel_lvss_svs(_vel_vfrmaxsfst_vvl(vec, 256),0);
    }

    inline InaVecSXA sqrt() const {
        return _vel_vfsqrts_vvl(vec, 256);
    }

    inline InaVecSXA exp() const {
        const __vr COEFF_LOG2E = _vel_vbrds_vsl(float(InaFastExp::CoeffLog2E()), 256);
        const __vr COEFF_A     = _vel_vbrds_vsl(float(InaFastExp::CoeffA32()), 256);
        const __vr COEFF_B     = _vel_vbrds_vsl(float(InaFastExp::CoeffB32()), 256);
        const __vr COEFF_P5_A  = _vel_vbrds_vsl(float(InaFastExp::GetCoefficient6_5()), 256);
        const __vr COEFF_P5_B  = _vel_vbrds_vsl(float(InaFastExp::GetCoefficient6_4()), 256);
        const __vr COEFF_P5_C  = _vel_vbrds_vsl(float(InaFastExp::GetCoefficient6_3()), 256);
        const __vr COEFF_P5_D  = _vel_vbrds_vsl(float(InaFastExp::GetCoefficient6_2()), 256);
        const __vr COEFF_P5_E  = _vel_vbrds_vsl(float(InaFastExp::GetCoefficient6_1()), 256);
        const __vr COEFF_P5_F  = _vel_vbrds_vsl(float(InaFastExp::GetCoefficient6_0()), 256);

        __vr x = _vel_vfmuls_vvvl(vec, COEFF_LOG2E, 256);

        const __vr fractional_part = _vel_vfsubs_vvvl(x, InaVecSXA(x).floor().vec, 256);

        __vr factor = _vel_vfadds_vvvl(
                         _vel_vfmuls_vvvl(_vel_vfadds_vvvl( _vel_vfmuls_vvvl(_vel_vfadds_vvvl(
                         _vel_vfmuls_vvvl(_vel_vfadds_vvvl( _vel_vfmuls_vvvl(_vel_vfadds_vvvl(
                         _vel_vfmuls_vvvl(COEFF_P5_A, fractional_part, 256),
                         COEFF_P5_B, 256), fractional_part, 256), COEFF_P5_C, 256),fractional_part, 256),
                         COEFF_P5_D, 256), fractional_part, 256), COEFF_P5_E, 256),fractional_part, 256),
                         COEFF_P5_F, 256);

        x = _vel_vfsubs_vvvl(x,factor, 256);

        x = _vel_vfadds_vvvl(_vel_vfmuls_vvvl(COEFF_A, x, 256), COEFF_B, 256);

        __vr castedInteger = _vel_vcvtwssxrz_vvl(x, 256);
        castedInteger = _vel_vsll_vvsl(castedInteger, 32, 256);

        return (castedInteger); // Automatically reinterpret not cast
    }

    inline InaVecSXA expLowAcc() const {
        const __vr COEFF_LOG2E = _vel_vbrds_vsl(float(InaFastExp::CoeffLog2E()), 256);
        const __vr COEFF_A     = _vel_vbrds_vsl(float(InaFastExp::CoeffA32()), 256);
        const __vr COEFF_B     = _vel_vbrds_vsl(float(InaFastExp::CoeffB32()), 256);
        const __vr COEFF_P5_D  = _vel_vbrds_vsl(float(InaFastExp::GetCoefficient3_2()), 256);
        const __vr COEFF_P5_E  = _vel_vbrds_vsl(float(InaFastExp::GetCoefficient3_1()), 256);
        const __vr COEFF_P5_F  = _vel_vbrds_vsl(float(InaFastExp::GetCoefficient3_0()), 256);

        __vr x = _vel_vfmuls_vvvl(vec, COEFF_LOG2E, 256);

        const __vr fractional_part = _vel_vfsubs_vvvl(x, InaVecSXA(x).floor().vec, 256);

        __vr factor = _vel_vfadds_vvvl(_vel_vfmuls_vvvl(_vel_vfadds_vvvl(
                         _vel_vfmuls_vvvl(fractional_part, COEFF_P5_D, 256),
                                         COEFF_P5_E, 256), fractional_part, 256),
                                         COEFF_P5_F, 256);

        x = _vel_vfsubs_vvvl(x,factor, 256);

        x = _vel_vfadds_vvvl(_vel_vfmuls_vvvl(COEFF_A, x, 256), COEFF_B, 256);

        __vr castedInteger = _vel_vcvtwssxrz_vvl(x, 256);
        castedInteger = _vel_vsll_vvsl(castedInteger, 32, 256);

        return (castedInteger); // Automatically reinterpret not cast
    }

    inline InaVecSXA rsqrt() const {
        // svrsqrte_f64(vec); seems low accurate
        return  _vel_vrsqrts_vvl(vec, 256);
        //const __vr one = _vel_vbrds_vsl(1.0, 256);
        //return  _vel_vfsqrts_vvl(_vel_vfdivs_vvvl(one, vec, 256), 256);
    }

    inline InaVecSXA abs() const {
        return _vel_vand_vvvl( vec, _vel_pvbrd_vsl(0x7FFFFFFF7FFFFFFFUL, 256), 256);
    }

    inline InaVecSXA floor() const {
        __vr valuesInIntervals = _vel_vfmins_vvvl(
                                    _vel_vfmaxs_vvvl( vec, _vel_vbrds_vsl(double(std::numeric_limits<int>::min()), 256), 256),
                                    _vel_vbrds_vsl(double(std::numeric_limits<int>::max()), 256), 256);
        __vr vecConvLongInt = _vel_vcvtwssxrz_vvl(valuesInIntervals, 256);
        __vr vecConvLongIntDouble = _vel_vcvtsw_vvl(vecConvLongInt, 256);
        __vm256 maskPositive = _vel_vfmklgt_mvl(_vel_vfcmps_vsvl( 0, vec, 256), 256);
        return _vel_vmrg_vvvml(vecConvLongIntDouble,
                               _vel_vfsubs_vvvl( vecConvLongIntDouble, _vel_vbrds_vsl(1, 256), 256),
                               maskPositive,
                               256);
    }

    inline InaVecSXA signOf() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        __vm256 maskPositive = _vel_vfmkllt_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);
        __vm256 maskNegative = _vel_vfmkllt_mvl(_vel_vfcmps_vvvl( vec, zero, 256), 256);

        return _vel_vmrg_vvvml(_vel_vmrg_vvvml(zero,
                                               _vel_vbrds_vsl(1, 256),
                                                maskPositive, 256),
                               _vel_vbrds_vsl(-1, 256),
                               maskNegative, 256);
    }

    inline InaVecSXA isPositive() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        __vm256 maskPositive = _vel_vfmklle_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);

        return _vel_vmrg_vvvml(zero,
                               _vel_vbrds_vsl(1, 256),
                               maskPositive, 256);
    }

    inline InaVecSXA isNegative() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        __vm256 maskNegative = _vel_vfmklle_mvl(_vel_vfcmps_vvvl( vec, zero, 256), 256);

        return _vel_vmrg_vvvml(zero,
                               _vel_vbrds_vsl(1, 256),
                               maskNegative, 256);
    }

    inline InaVecSXA isPositiveStrict() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        __vm256 maskPositive = _vel_vfmkllt_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);

        return _vel_vmrg_vvvml(zero,
                               _vel_vbrds_vsl(1, 256),
                               maskPositive, 256);
    }

    inline InaVecSXA isNegativeStrict() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        __vm256 maskNegative = _vel_vfmkllt_mvl(_vel_vfcmps_vvvl( vec, zero, 256), 256);

        return _vel_vmrg_vvvml( zero,
                                _vel_vbrds_vsl(1, 256),
                               maskNegative, 256);
    }

    inline InaVecSXA isZero() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        __vm256 maskEqual = _vel_vfmkleq_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);

        return _vel_vmrg_vvvml(zero,
                               _vel_vbrds_vsl(1, 256),
                               maskEqual, 256);
    }

    inline InaVecSXA isNotZero() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        __vm256 maskEqual = _vel_vfmkleq_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);

        return _vel_vmrg_vvvml(_vel_vbrds_vsl(1, 256),
                               zero,
                               maskEqual, 256);
    }

    inline InaVecMaskSXA<float> isPositiveMask() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        return _vel_vfmklle_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);
    }

    inline InaVecMaskSXA<float> isNegativeMask() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        return _vel_vfmklge_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);
    }

    inline InaVecMaskSXA<float> isPositiveStrictMask() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        return _vel_vfmkllt_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);
    }

    inline InaVecMaskSXA<float> isNegativeStrictMask() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        return _vel_vfmklgt_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);
    }

    inline InaVecMaskSXA<float> isZeroMask() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        return _vel_vfmkleq_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);
    }

    inline InaVecMaskSXA<float> isNotZeroMask() const {
        __vr zero = _vel_vbrds_vsl(0, 256);
        return _vel_vfmklne_mvl(_vel_vfcmps_vvvl( zero, vec, 256), 256);
    }

    // Static basic methods
    inline static InaVecSXA GetZero() {
        return InaVecSXA(_vel_vbrds_vsl(0, 256));
    }

    inline static InaVecSXA GetOne() {
        return InaVecSXA(_vel_vbrds_vsl(1, 256));
    }

    inline static InaVecSXA Min(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmins_vvvl( inVec1.vec, inVec2.vec, 256);
    }

    inline static InaVecSXA Max(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmaxs_vvvl( inVec1.vec, inVec2.vec, 256);
    }

    inline static InaVecSXA IsLowerOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmklle_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrds_vsl(0, 256),
                               _vel_vbrds_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecSXA IsLower(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmkllt_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrds_vsl(0, 256),
                               _vel_vbrds_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecSXA IsGreaterOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmklge_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrds_vsl(0, 256),
                               _vel_vbrds_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecSXA IsGreater(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmklgt_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrds_vsl(0, 256),
                               _vel_vbrds_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecSXA IsEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmkleq_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrds_vsl(0, 256),
                               _vel_vbrds_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecSXA IsNotEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmklne_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrds_vsl(0, 256),
                               _vel_vbrds_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecMaskSXA<float> IsLowerOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklle_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<float> IsLowerMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmkllt_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<float> IsGreaterOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklge_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<float> IsGreaterMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklgt_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<float> IsEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmkleq_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<float> IsNotEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklne_mvl(_vel_vfcmps_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecSXA BitsAnd(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vand_vvvl(inVec1.vec, inVec2.vec, 256);
    }

    inline static InaVecSXA BitsNotAnd(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vand_vvvl(_vel_vxor_vvvl(inVec1.vec, _vel_pvbrd_vsl(~0UL, 256), 256), inVec2.vec, 256);
    }

    inline static InaVecSXA BitsOr(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vor_vvvl(inVec1.vec, inVec2.vec, 256);
    }

    inline static InaVecSXA BitsXor(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vxor_vvvl(inVec1.vec, inVec2.vec, 256);
    }

    inline static  const char* GetName() {
        return "InaVecSXA<float>";
    }

    inline static  InaIfElse< InaVecSXA<float> >::ThenClass If(const InaVecMaskSXA<float>& inTest) {
       return InaIfElse< InaVecSXA<float> >::IfClass().If(inTest);
    }

    inline static InaVecSXA IfElse(const InaVecMaskSXA<float>& inMask, const InaVecSXA& inIfTrue, const InaVecSXA& inIfFalse) {
        return _vel_vmrg_vvvml(inIfFalse.vec, inIfTrue.vec, __vm256(inMask), 256);
    }

    inline static InaVecSXA IfTrue(const InaVecMaskSXA<float>& inMask, const InaVecSXA& inIfTrue) {
        return _vel_vmrg_vvvml(_vel_vbrds_vsl(0, 256), inIfTrue.vec, __vm256(inMask), 256);
    }

    inline static InaVecSXA IfFalse(const InaVecMaskSXA<float>& inMask, const InaVecSXA& inIfFalse) {
        return _vel_vmrg_vvvml(inIfFalse.vec, _vel_vbrds_vsl(0, 256), __vm256(inMask), 256);
    }

    // Inner operators
    inline InaVecSXA<float>& operator+=(const InaVecSXA<float>& inVec){
        vec = _vel_vfadds_vvvl(vec,inVec.vec, 256);
        return *this;
    }

    inline InaVecSXA<float>& operator-=(const InaVecSXA<float>& inVec){
        vec = _vel_vfsubs_vvvl(vec,inVec.vec, 256);
        return *this;
    }

    inline InaVecSXA<float>& operator/=(const InaVecSXA<float>& inVec){
        vec = _vel_vfdivs_vvvl(vec,inVec.vec, 256);
        return *this;
    }

    inline InaVecSXA<float>& operator*=(const InaVecSXA<float>& inVec){
        vec = _vel_vfmuls_vvvl(vec,inVec.vec, 256);
        return *this;
    }

    inline InaVecSXA<float> operator-() const {
        return _vel_vxor_vvvl(vec, _vel_pvbrd_vsl(0x8000000080000000UL, 256), 256);
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

    inline static InaVecSXA<float> Fma(const InaVecSXA<float>& inValAdd, const InaVecSXA<float>& inValMul1, const InaVecSXA<float>& inValMul2){
        return _vel_vfmads_vvvvl(inValAdd.vec, inValMul1.vec, inValMul2.vec);
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
    return _vel_vfadds_vvvl(inVec1.getVec(), inVec2.getVec(), 256);
}

inline InaVecSXA<float> operator-(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return _vel_vfsubs_vvvl(inVec1.getVec(), inVec2.getVec(), 256);
}

inline InaVecSXA<float> operator/(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return _vel_vfdivs_vvvl(inVec1.getVec(), inVec2.getVec(), 256);
}

inline InaVecSXA<float> operator*(const InaVecSXA<float>& inVec1, const InaVecSXA<float>& inVec2){
    return _vel_vfmuls_vvvl(inVec1.getVec(), inVec2.getVec(), 256);
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
