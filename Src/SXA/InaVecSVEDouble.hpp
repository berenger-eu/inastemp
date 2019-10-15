///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSXADOUBLE_HPP
#define INAVECSXADOUBLE_HPP

#include "InaIfElse.hpp"
#include "InaUtils.hpp"
#include "InaFastExp.hpp"

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
class InaVecMaskSXA<double> {
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
        mask = (inBool? _vel_vfmklat_ml(0xFFFFFFFFU) : _vel_vfmklat_ml(0));
    }

    inline InaVecMaskSXA& operator=(const bool inBool){
        mask = (inBool? _vel_vfmklat_ml(0xFFFFFFFFU) : _vel_vfmklat_ml(0));
        return (*this);
    }

    // Binary methods
    inline InaVecMaskSXA Not() const{
        return _vel_negm_mm(mask);
    }

    inline bool isAllTrue() const{
        return _vel_pcvm_sml(mask, 256) == 256*8*8;
    }

    inline bool isAllFalse() const{
        // true if all zero
        return _vel_pcvm_sml(mask, 256) == 0;
    }

    // Double args methods
    inline static InaVecMaskSXA And(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_andm_mmm(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSXA NotAnd(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_andm_mmm(_vel_xorm_mmm(inMask1.mask, _vel_vfmklat_ml(0xFFFFFFFFU)),inMask2.mask);
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
inline InaVecMaskSXA<double> operator&(const InaVecMaskSXA<double>& inMask1, const InaVecMaskSXA<double>& inMask2){
    return InaVecMaskSXA<double>::And(inMask1, inMask2);
}

inline InaVecMaskSXA<double> operator|(const InaVecMaskSXA<double>& inMask1, const InaVecMaskSXA<double>& inMask2){
    return InaVecMaskSXA<double>::Or(inMask1, inMask2);
}

inline InaVecMaskSXA<double> operator^(const InaVecMaskSXA<double>& inMask1, const InaVecMaskSXA<double>& inMask2){
    return InaVecMaskSXA<double>::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskSXA<double>& inMask1, const InaVecMaskSXA<double>& inMask2){
    return InaVecMaskSXA<double>::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskSXA<double>& inMask1, const InaVecMaskSXA<double>& inMask2){
    return InaVecMaskSXA<double>::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class InaVecSXA<double> {
protected:
    __vr vec;

public:
    using VecRawType           = __vr;
    using MaskType             = InaVecMaskSXA<double>;
    using RealType             = double;
    static const int Alignement= 1;
    static const bool IsOfFixedSize = true;

    static constexpr int GetVecLength(){
        return 256;
    }

    inline InaVecSXA() {
        vec = _vel_vbrdd_vsl(0,256);
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
    inline /*not explicit*/ InaVecSXA(const double val)
        : InaVecSXA() {
        vec = _vel_vbrdd_vsl(val,256);
    }

    inline InaVecSXA& operator=(const double val){
        vec = _vel_vbrdd_vsl(val,256);
        return *this;
    }

    inline void setFromScalar(const double val){
        vec = _vel_vbrdd_vsl(val,256);
    }

    // Constructor from vec
    inline InaVecSXA(const std::initializer_list<double> lst)
        : InaVecSXA(lst.begin()){
    }

    inline explicit InaVecSXA(const double ptr[])
        : InaVecSXA() {
        vec = _vel_vld_vssl(8, ptr, 256);
    }

    inline InaVecSXA& setFromArray(const double ptr[]){
        vec = _vel_vld_vssl(8, ptr, 256);
        return *this;
    }

    inline InaVecSXA& setFromAlignedArray(const double ptr[]){
        vec = _vel_vld_vssl(8, ptr, 256);
        return *this;
    }

    inline InaVecSXA& setFromIndirectArray(const double values[], const int inIndirection[]) {
        __vr offset = _vel_vld_vssl(4, inIndirection, 256);
        vec = _vel_vld_vssvl(0, values, offset, 256);
        return *this;
    }

    inline InaVecSXA& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        __vr offset = _vel_vmulul_vvvl(_vel_vld_vssl(4, inIndirection1, 256),
                                      _vel_vld_vssl(4, inIndirection2, 256),
                                      256);
        vec = _vel_vld_vssvl(0, inArray, offset, 256);
        return *this;
    }

    // Move back to array
    inline void storeInArray(double ptr[]) const {
        _vel_vst_vssl(vec, 8, ptr, 256);
    }

    inline void storeInAlignedArray(double ptr[]) const {
        _vel_vst_vssl(vec, 8, ptr, 256);
    }

    // Acce to individual values
    inline double at(const int index) const {
        return _vel_lvsl_svs(vec, index);
    }

    // Horizontal operation
    inline double horizontalSum() const {
      return _vel_lvsl_svs(_vel_vfsumd_vvl(vec, 256),0);
    }

    inline double horizontalMul() const {
        double sum = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            sum *= at(idx);
        }
        return sum;
    }

    inline InaVecSXA sqrt() const {
        return _vel_vfsqrtd_vvl(vec, 256);
    }

    inline InaVecSXA exp() const {
        const __vr COEFF_LOG2E = _vel_vbrdd_vsl(double(InaFastExp::CoeffLog2E()), 256);
        const __vr COEFF_A     = _vel_vbrdd_vsl(double(InaFastExp::CoeffA64()), 256);
        const __vr COEFF_B     = _vel_vbrdd_vsl(double(InaFastExp::CoeffB64()), 256);
        const __vr COEFF_P5_X  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient9_8()), 256);
        const __vr COEFF_P5_Y  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient9_7()), 256);
        const __vr COEFF_P5_Z  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient9_6()), 256);
        const __vr COEFF_P5_A  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient9_5()), 256);
        const __vr COEFF_P5_B  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient9_4()), 256);
        const __vr COEFF_P5_C  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient9_3()), 256);
        const __vr COEFF_P5_D  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient9_2()), 256);
        const __vr COEFF_P5_E  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient9_1()), 256);
        const __vr COEFF_P5_F  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient9_0()), 256);

        __vr x = _vel_vfmuld_vvvl(vec, COEFF_LOG2E, 256);

        const __vr fractional_part = _vel_vfsubd_vvvl(x, InaVecSXA(x).floor().vec, 256);

        __vr factor = _vel_vfaddd_vvvl(_vel_vfmuld_vvvl(_vel_vfaddd_vvvl(
                         _vel_vfmuld_vvvl(_vel_vfaddd_vvvl( _vel_vfmuld_vvvl(_vel_vfaddd_vvvl(
                         _vel_vfmuld_vvvl(_vel_vfaddd_vvvl( _vel_vfmuld_vvvl(_vel_vfaddd_vvvl(
                         _vel_vfmuld_vvvl(_vel_vfaddd_vvvl( _vel_vfmuld_vvvl(_vel_vfaddd_vvvl(_vel_vfmuld_vvvl(
                         COEFF_P5_X, fractional_part, 256), COEFF_P5_Y, 256), fractional_part, 256),
                         COEFF_P5_Z, 256),fractional_part, 256), COEFF_P5_A, 256), fractional_part, 256),
                         COEFF_P5_B, 256), fractional_part, 256), COEFF_P5_C, 256),fractional_part, 256),
                         COEFF_P5_D, 256), fractional_part, 256), COEFF_P5_E, 256),fractional_part, 256),
                         COEFF_P5_F, 256);

        x = _vel_vfsubd_vvvl(x,factor, 256);

        x = _vel_vfaddd_vvvl(_vel_vfmuld_vvvl(COEFF_A, x, 256), COEFF_B, 256);

        __vr castedInteger = _vel_vcvtldrz_vvl(x, 256);

        return (castedInteger); // Automatically reinterpret not cast
    }

    inline InaVecSXA expLowAcc() const {
        const __vr COEFF_LOG2E = _vel_vbrdd_vsl(double(InaFastExp::CoeffLog2E()), 256);
        const __vr COEFF_A     = _vel_vbrdd_vsl(double(InaFastExp::CoeffA64()), 256);
        const __vr COEFF_B     = _vel_vbrdd_vsl(double(InaFastExp::CoeffB64()), 256);
        const __vr COEFF_P5_C  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient4_3()), 256);
        const __vr COEFF_P5_D  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient4_2()), 256);
        const __vr COEFF_P5_E  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient4_1()), 256);
        const __vr COEFF_P5_F  = _vel_vbrdd_vsl(double(InaFastExp::GetCoefficient4_0()), 256);

        __vr x = _vel_vfmuld_vvvl(vec, COEFF_LOG2E, 256);

        const __vr fractional_part = _vel_vfsubd_vvvl(x, InaVecSXA(x).floor().vec, 256);

        __vr factor = _vel_vfaddd_vvvl(_vel_vfmuld_vvvl(_vel_vfaddd_vvvl(
                         _vel_vfmuld_vvvl(_vel_vfaddd_vvvl(_vel_vfmuld_vvvl(
                                         COEFF_P5_C, fractional_part, 256),
                                         COEFF_P5_D, 256), fractional_part, 256),
                                         COEFF_P5_E, 256), fractional_part, 256),
                                         COEFF_P5_F, 256);

        x = _vel_vfsubd_vvvl(x,factor, 256);

        x = _vel_vfaddd_vvvl(_vel_vfmuld_vvvl(COEFF_A, x, 256), COEFF_B, 256);

        __vr castedInteger = _vel_vcvtldrz_vvl(x, 256);

        return (castedInteger); // Automatically reinterpret not cast
    }

    inline InaVecSXA rsqrt() const {
        // svrsqrte_f64(vec); seems low accurate
        return  _vel_vrsqrtd_vvl(vec, 256);
    }

    inline InaVecSXA abs() const {
      return _vel_vand_vvvl( vec, _vel_pvbrd_vsl(0x7FFFFFFFFFFFFFFFUL, 256), 256);
    }

    inline InaVecSXA floor() const {
        __vm256 maskInLongInt = _vel_andm_mmm(
                                _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( _vel_vbrdd_vsl(double(std::numeric_limits<long int>::min()), 256), vec, 256), 256),
                                _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( vec, _vel_vbrdd_vsl(double(std::numeric_limits<long int>::max()), 256), 256), 256));
        __vr vecConvLongInt = _vel_vcvtldrz_vvl(vec, 256);
        __vr vecConvLongIntDouble = _vel_vcvtdl_vvl(vec, 256);
        __vm256 maskNegative = _vel_vfmklgt_mvl(_vel_vcmpsl_vsvl( 0, vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vmrg_vvvml(_vel_vfsubd_vvvl( vecConvLongIntDouble, _vel_vbrdd_vsl(1, 256), 256),
                                               vecConvLongIntDouble,
                                               maskNegative,
                                               256),
                               vec,
                               maskInLongInt,
                               256);
    }

    inline InaVecSXA signOf() const {
        return _vel_vfcmpd_vvvl(vec, _vel_vbrdd_vsl(0, 256), 256);
    }

    inline InaVecSXA isPositive() const {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(isPositiveMask()), 256);
    }

    inline InaVecSXA isNegative() const {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(isNegativeMask()), 256);
    }

    inline InaVecSXA isPositiveStrict() const {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(isPositiveStrictMask()), 256);
    }

    inline InaVecSXA isNegativeStrict() const {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(isNegativeStrictMask()), 256);
    }

    inline InaVecSXA isZero() const {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(isZeroMask()), 256);
    }

    inline InaVecSXA isNotZero() const {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(isNotZeroMask()), 256);
    }

    inline InaVecMaskSXA<double> isPositiveMask() const {
        return _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( _vel_vbrdd_vsl(0-std::numeric_limits<double>::epsilon(), 256), vec, 256), 256);
    }

    inline InaVecMaskSXA<double> isNegativeMask() const {
        return _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( vec, _vel_vbrdd_vsl(0+std::numeric_limits<double>::epsilon(), 256), 256), 256);
    }

    inline InaVecMaskSXA<double> isPositiveStrictMask() const {
        return _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( _vel_vbrdd_vsl(0-std::numeric_limits<double>::epsilon(), 256), vec, 256), 256);
    }

    inline InaVecMaskSXA<double> isNegativeStrictMask() const {
        return _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( vec, _vel_vbrdd_vsl(0+std::numeric_limits<double>::epsilon(), 256), 256), 256);
    }

    inline InaVecMaskSXA<double> isZeroMask() const {
        return _vel_negm_mm(_vel_vfmklgt_mvl(_vel_vand_vvvl( _vel_pvbrd_vsl(~0UL, 256), vec, 256), 256));
    }

    inline InaVecMaskSXA<double> isNotZeroMask() const {
        return _vel_vfmklgt_mvl(_vel_vand_vvvl( _vel_pvbrd_vsl(~0UL, 256), vec, 256), 256);
    }

    // Static basic methods
    inline static InaVecSXA GetZero() {
        return InaVecSXA(_vel_vbrdd_vsl(0, 256));
    }

    inline static InaVecSXA GetOne() {
        return InaVecSXA(_vel_vbrdd_vsl(1, 256));
    }

    inline static InaVecSXA Min(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vminswsx_vvvl( inVec1.vec, inVec2.vec, 256);
    }

    inline static InaVecSXA Max(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vmaxswsx_vvvl( inVec1.vec, inVec2.vec, 256);
    }

    inline static InaVecSXA IsLowerOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(IsLowerOrEqualMask(inVec1, inVec2)), 256);
    }

    inline static InaVecSXA IsLower(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(IsLowerMask(inVec1, inVec2)), 256);
    }

    inline static InaVecSXA IsGreaterOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(IsGreaterOrEqualMask(inVec1, inVec2)), 256);
    }

    inline static InaVecSXA IsGreater(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(IsGreaterMask(inVec1, inVec2)), 256);
    }

    inline static InaVecSXA IsEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(IsEqualMask(inVec1, inVec2)), 256);
    }

    inline static InaVecSXA IsNotEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256), _vel_vbrdd_vsl(0, 256), __vm256(IsNotEqualMask(inVec1, inVec2)), 256);
    }

    inline static InaVecMaskSXA<double> IsLowerOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( _vel_vfsubd_vvvl(inVec1.vec, _vel_vbrdd_vsl(std::numeric_limits<double>::epsilon(), 256), 256), inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<double> IsLowerMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<double> IsGreaterOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( inVec2.vec, inVec1.vec, 256), 256);
    }

    inline static InaVecMaskSXA<double> IsGreaterMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( _vel_vfsubd_vvvl(inVec2.vec, _vel_vbrdd_vsl(std::numeric_limits<double>::epsilon(), 256), 256), inVec1.vec, 256), 256);
    }

    inline static InaVecMaskSXA<double> IsEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return  InaVecMaskSXA<double>::And(IsGreaterOrEqualMask(inVec1, inVec2), IsLowerOrEqualMask(inVec1, inVec2));
    }

    inline static InaVecMaskSXA<double> IsNotEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return  InaVecMaskSXA<double>::Xor(IsGreaterOrEqualMask(inVec1, inVec2), IsLowerOrEqualMask(inVec1, inVec2));
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
        return "InaVecSXA<double>";
    }

    inline static  InaIfElse< InaVecSXA<double> >::ThenClass If(const InaVecMaskSXA<double>& inTest) {
       return InaIfElse< InaVecSXA<double> >::IfClass().If(inTest);
    }

    inline static InaVecSXA IfElse(const InaVecMaskSXA<double>& inMask, const InaVecSXA& inIfTrue, const InaVecSXA& inIfFalse) {
        return _vel_vmrg_vvvml(inIfTrue.vec, inIfFalse.vec, __vm256(inMask), 256);
    }

    inline static InaVecSXA IfTrue(const InaVecMaskSXA<double>& inMask, const InaVecSXA& inIfTrue) {
        return _vel_vmrg_vvvml(inIfTrue.vec, _vel_vbrdd_vsl(0, 256), __vm256(inMask), 256);
    }

    inline static InaVecSXA IfFalse(const InaVecMaskSXA<double>& inMask, const InaVecSXA& inIfFalse) {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(0, 256), inIfFalse.vec, __vm256(inMask), 256);
    }

    // Inner operators
    inline InaVecSXA<double>& operator+=(const InaVecSXA<double>& inVec){
        vec = _vel_vfaddd_vvvl(vec,inVec.vec, 256);
        return *this;
    }

    inline InaVecSXA<double>& operator-=(const InaVecSXA<double>& inVec){
        vec = _vel_vfsubd_vvvl(vec,inVec.vec, 256);
        return *this;
    }

    inline InaVecSXA<double>& operator/=(const InaVecSXA<double>& inVec){
        vec = _vel_vfdivd_vvvl(vec,inVec.vec, 256);
        return *this;
    }

    inline InaVecSXA<double>& operator*=(const InaVecSXA<double>& inVec){
        vec = _vel_vfmuld_vvvl(vec,inVec.vec, 256);
        return *this;
    }

    inline InaVecSXA<double> operator-() const {
        return _vel_vxor_vvvl(vec, _vel_pvbrd_vsl(0x8000000000000000UL, 256), 256);
    }

    inline InaVecSXA<double> pow(std::size_t power) const{
        return InaUtils::FastPow<InaVecSXA<double>>(*this, power);
    }

    // Multiple sum
    template <class ... Args>
    inline static void MultiHorizontalSum(double sumRes[], const InaVecSXA<double>& inVec1,
                                          const InaVecSXA<double>& inVec2, const InaVecSXA<double>& inVec3,
                                          const InaVecSXA<double>& inVec4, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1, inVec2 );
        MultiHorizontalSum(&sumRes[2], inVec3, inVec4 );

        MultiHorizontalSum(&sumRes[4], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(double sumRes[], const InaVecSXA<double>& inVec1,
                                          const InaVecSXA<double>& inVec2, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1);
        MultiHorizontalSum(&sumRes[1], inVec2);

        MultiHorizontalSum(&sumRes[2], args... );
    }

    inline static void MultiHorizontalSum(double sumRes[], const InaVecSXA<double>& inVec){
        sumRes[0] += inVec.horizontalSum();
    }

    inline static void MultiHorizontalSum(double /*sumRes*/[]){
    }
};

// Bits operators
inline InaVecSXA<double> operator&(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return InaVecSXA<double>::BitsAnd(inVec1, inVec2);
}

inline InaVecSXA<double> operator|(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return InaVecSXA<double>::BitsOr(inVec1, inVec2);
}

inline InaVecSXA<double> operator^(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return InaVecSXA<double>::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecSXA<double> operator+(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return _vel_vfaddd_vvvl(inVec1.getVec(), inVec2.getVec(), 256);
}

inline InaVecSXA<double> operator-(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return _vel_vfsubd_vvvl(inVec1.getVec(), inVec2.getVec(), 256);
}

inline InaVecSXA<double> operator/(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return _vel_vfdivd_vvvl(inVec1.getVec(), inVec2.getVec(), 256);
}

inline InaVecSXA<double> operator*(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return _vel_vfmuld_vvvl(inVec1.getVec(), inVec2.getVec(), 256);
}

// Tests and comparions
inline InaVecMaskSXA<double> operator<(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return InaVecSXA<double>::IsLowerMask(inVec1,inVec2);
}

inline InaVecMaskSXA<double> operator<=(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return InaVecSXA<double>::IsLowerOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskSXA<double> operator>(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return InaVecSXA<double>::IsGreaterMask(inVec1,inVec2);
}

inline InaVecMaskSXA<double> operator>=(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return InaVecSXA<double>::IsGreaterOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskSXA<double> operator==(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return InaVecSXA<double>::IsEqualMask(inVec1,inVec2);
}

inline InaVecMaskSXA<double> operator!=(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return InaVecSXA<double>::IsNotEqualMask(inVec1,inVec2);
}


#endif
