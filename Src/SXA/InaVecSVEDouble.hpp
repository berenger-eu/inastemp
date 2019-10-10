///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSXADOUBLE_HPP
#define INAVECSXADOUBLE_HPP

#include "InastempGlobal.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_SXA
#error InaVecSXA<double> is included but SXA is not enable in the configuration
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
class InaVecMaskSXA<double> {
    __vm mask;

public:
    // Classic constructors
    inline InaVecMaskSXA() {
        _ve_lvl(256);
        mask = _vel_vbrdw_vsl(0,256)
    }

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
        mask = (inBool? _vel_vbrdw_vsl(0xFFFFFFFFU,256) : _vel_vbrdw_vsl(0,256));
    }

    inline InaVecMaskSXA& operator=(const bool inBool){
        mask = (inBool? _vel_vbrdw_vsl(0xFFFFFFFFU,256) : _vel_vbrdw_vsl(0,256));
        return (*this);
    }

    // Binary methods
    inline InaVecMaskSXA Not() const{
        return _vel_vxor_vvvl(mask, _vel_vbrdw_vsl(0xFFFFFFFFU,256), 256);
    }

    inline bool isAllTrue() const{
        return _vel_lvsl_svs(_vel_vrand_vvl(mask, 256), 0) == 0xFFFFFFFFFFFFFFFFUL;
    }

    inline bool isAllFalse() const{
        // true if all zero
        return _vel_lvsl_svs(_vel_vror_vvl(mask, 256),0) == 0xFFFFFFFFFFFFFFFFUL == 0;
    }

    // Double args methods
    inline static InaVecMaskSXA And(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_vand_vvvl(inMask1.mask,inMask2.mask, 256);
    }

    inline static InaVecMaskSXA NotAnd(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_vand_vvvl(_vel_vxor_vvvl(inMask1.mask, _vel_vbrdw_vsl(0xFFFFFFFFU,256), 256),inMask2.mask, 256);
    }

    inline static InaVecMaskSXA Or(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_vor_vvvl(inMask1.mask,inMask2.mask, 256);
    }

    inline static InaVecMaskSXA Xor(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_vxor_vvvl(inMask1.mask,inMask2.mask, 256);
    }

    inline static bool IsEqual(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_lvsl_svs(_vel_vror_vvl(_vel_vcmpul_vvvl(inMask1.mask,inMask2.mask, 256),256), 0) == 0;
    }

    inline static bool IsNotEqual(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_lvsl_svs(_vel_vror_vvl(_vel_vcmpul_vvvl(inMask1.mask,inMask2.mask, 256),256), 0) != 0;
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
    static const bool IsOfFixedSize = false;
    
    static constexpr int GetVecLength(){
        return 256;
    }

    inline InaVecSXA() {
        _ve_lvl(256);
        mask = _vel_vbrdd_vsvl(0,256);
    }
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
    inline /*not explicit*/ InaVecSXA(const double val)
        : InaVecSXA() {
        vec = _vel_vbrdd_vsvl(val,256);
    }

    inline InaVecSXA& operator=(const double val){
        vec = _vel_vbrdd_vsvl(val,256);
        return *this;
    }

    inline void setFromScalar(const double val){
        vec = _vel_vbrdd_vsvl(val,256);
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
        _ve offset = _vel_vld_vssl(4, inIndirection, 256);
        vec = _vel_vld_vssvl(0, values, offset, 256);
        return *this;
    }

    inline InaVecSXA& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){                                       
        _ve offset = _vel_vmulul_vvvl(_vel_vld_vssl(4, inIndirection1, 256),
                                      _vel_vld_vssl(4, inIndirection2, 256),
                                      256);
        vec = _vel_vld_vssvl(0, values, offset, 256);
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
        return _vel_lvsl_svs(mask, index);
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
        const __vr COEFF_LOG2E = _vel_vbrdd_vsvl(double(InaFastExp::CoeffLog2E()), 256);
        const __vr COEFF_A     = _vel_vbrdd_vsvl(double(InaFastExp::CoeffA64()), 256);
        const __vr COEFF_B     = _vel_vbrdd_vsvl(double(InaFastExp::CoeffB64()), 256);
        const __vr COEFF_P5_X  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient9_8()), 256);
        const __vr COEFF_P5_Y  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient9_7()), 256);
        const __vr COEFF_P5_Z  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient9_6()), 256);
        const __vr COEFF_P5_A  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient9_5()), 256);
        const __vr COEFF_P5_B  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient9_4()), 256);
        const __vr COEFF_P5_C  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient9_3()), 256);
        const __vr COEFF_P5_D  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient9_2()), 256);
        const __vr COEFF_P5_E  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient9_1()), 256);
        const __vr COEFF_P5_F  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient9_0()), 256);

        __vr x = _vel_vfmuld_vvvl(vec, COEFF_LOG2E);

        const __vr fractional_part = _vel_vfsubd_vvvl(x, InaVecSXA(x).floor().vec);

        __vr factor = _vel_vfaddd_vvvl(_vel_vfmuld_vvvl(_vel_vfaddd_vvvl(
                         _vel_vfmuld_vvvl(_vel_vfaddd_vvvl( _vel_vfmuld_vvvl(_vel_vfaddd_vvvl(
                         _vel_vfmuld_vvvl(_vel_vfaddd_vvvl( _vel_vfmuld_vvvl(_vel_vfaddd_vvvl(
                         _vel_vfmuld_vvvl(_vel_vfaddd_vvvl( _vel_vfmuld_vvvl(_vel_vfaddd_vvvl(_vel_vfmuld_vvvl(
                         COEFF_P5_X, fractional_part), COEFF_P5_Y), fractional_part),
                         COEFF_P5_Z),fractional_part), COEFF_P5_A), fractional_part),
                         COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part),
                         COEFF_P5_F);

        x = _vel_vfsubd_vvvl(x,factor);

        x = _vel_vfaddd_vvvl(_vel_vfmuld_vvvl(COEFF_A, x), COEFF_B);

        __vr castedInteger = _vel_vcvtldrz_vvl(x, 256);

        return (castedInteger); // Automatically reinterpret not cast
    }

    inline InaVecSXA expLowAcc() const {
        const __vr COEFF_LOG2E = _vel_vbrdd_vsvl(double(InaFastExp::CoeffLog2E()));
        const __vr COEFF_A     = _vel_vbrdd_vsvl(double(InaFastExp::CoeffA64()));
        const __vr COEFF_B     = _vel_vbrdd_vsvl(double(InaFastExp::CoeffB64()));
        const __vr COEFF_P5_C  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient4_3()));
        const __vr COEFF_P5_D  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient4_2()));
        const __vr COEFF_P5_E  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient4_1()));
        const __vr COEFF_P5_F  = _vel_vbrdd_vsvl(double(InaFastExp::GetCoefficient4_0()));

        __vr x = _vel_vfmuld_vvvl(vec, COEFF_LOG2E);

        const __vr fractional_part = _vel_vfsubd_vvvl(x, InaVecSXA(x).floor().vec);

        __vr factor = _vel_vfaddd_vvvl(_vel_vfmuld_vvvl(_vel_vfaddd_vvvl(
                         _vel_vfmuld_vvvl(_vel_vfaddd_vvvl(_vel_vfmuld_vvvl(
                                         COEFF_P5_C, fractional_part),
                                         COEFF_P5_D), fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = _vel_vfsubd_vvvl(x,factor);

        x = _vel_vfaddd_vvvl(_vel_vfmuld_vvvl(COEFF_A, x), COEFF_B);

        __vr castedInteger = _vel_vcvtldrz_vvl(x, 256);

        return (castedInteger); // Automatically reinterpret not cast
    }

    inline InaVecSXA rsqrt() const {
        // svrsqrte_f64(vec); seems low accurate
        return  _vel_vrsqrtd_vvl(vec, 256);
    }

    inline InaVecSXA abs() const {
      return _vel_vand_vvvl( vec, _vel_vbrdw_vsl(0x7FFFFFFFFFFFFFFFL), 256);
    }

    inline InaVecSXA floor() const {
        __vr maskInLongInt = svand_z(
                                svcmple_f64( _vel_vbrdd_vsvl(double(std::numeric_limits<long int>::min())), vec),
                                svcmple_f64( vec, _vel_vbrdd_vsvl(double(std::numeric_limits<long int>::max()))));
        __vr vecConvLongInt = _vel_vcvtldrz_vvl(maskInLongInt, vec);
        __vr vecConvLongIntDouble = _vel_vcvtdl_vvl(maskInLongInt, vecConvLongInt);
        __vr maskNegative = _vel_vcmpsl_vsvl( 0, vec, 256);
        return svsel_f64(maskInLongInt, svsel_f64(maskNegative,  _vel_vfsubd_vvvl( vecConvLongIntDouble, _vel_vbrdd_vsvl(1)), vecConvLongIntDouble), vec);
    }

    inline InaVecSXA signOf() const {
        return _vel_vcmpul_vvvl(vec, _vel_vbrdd_vsvl(0), vec);
    }

    inline InaVecSXA isPositive() const {
        return svsel_f64(svcmple_f64( _vel_vbrdd_vsvl(0), vec),
                         _vel_vbrdd_vsvl(1), _vel_vbrdd_vsvl(0));
    }

    inline InaVecSXA isNegative() const {
        return svsel_f64(svcmpge_f64( _vel_vbrdd_vsvl(0), vec),
                         _vel_vbrdd_vsvl(1), _vel_vbrdd_vsvl(0));
    }

    inline InaVecSXA isPositiveStrict() const {
        return svsel_f64(svcmplt_f64( _vel_vbrdd_vsvl(0), vec),
                         _vel_vbrdd_vsvl(1), _vel_vbrdd_vsvl(0));
    }

    inline InaVecSXA isNegativeStrict() const {
        return svsel_f64(svcmpgt_f64( _vel_vbrdd_vsvl(0), vec),
                         _vel_vbrdd_vsvl(1), _vel_vbrdd_vsvl(0));
    }

    inline InaVecSXA isZero() const {
        return svsel_f64(svcmpeq_f64( _vel_vbrdd_vsvl(0), vec),
                         _vel_vbrdd_vsvl(1), _vel_vbrdd_vsvl(0));
    }

    inline InaVecSXA isNotZero() const {
        return svsel_f64(svcmpne_f64( _vel_vbrdd_vsvl(0), vec),
                         _vel_vbrdd_vsvl(1), _vel_vbrdd_vsvl(0));
    }

    inline InaVecMaskSXA<double> isPositiveMask() const {
        return svorr_z( svcmpeq_f64( _vel_vbrdd_vsvl(0), vec),
                       svcmple_f64( _vel_vbrdd_vsvl(0), vec));
    }

    inline InaVecMaskSXA<double> isNegativeMask() const {
        return svorr_z( svcmpeq_f64( _vel_vbrdd_vsvl(0), vec),
                       svcmpge_f64( _vel_vbrdd_vsvl(0), vec));
    }

    inline InaVecMaskSXA<double> isPositiveStrictMask() const {
        return svcmplt_f64( _vel_vbrdd_vsvl(0), vec);
    }

    inline InaVecMaskSXA<double> isNegativeStrictMask() const {
        return svcmpgt_f64( _vel_vbrdd_vsvl(0), vec);
    }

    inline InaVecMaskSXA<double> isZeroMask() const {
        return svcmpeq_f64( _vel_vbrdd_vsvl(0), vec);
    }

    inline InaVecMaskSXA<double> isNotZeroMask() const {
        return svcmpne_f64( _vel_vbrdd_vsvl(0), vec);
    }

    // Static basic methods
    inline static InaVecSXA GetZero() {
        return InaVecSXA(_vel_vbrdd_vsvl(0));
    }

    inline static InaVecSXA GetOne() {
        return InaVecSXA(_vel_vbrdd_vsvl(1));
    }

    inline static InaVecSXA Min(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svmin_f64_z( inVec1.vec, inVec2.vec);
    }

    inline static InaVecSXA Max(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svmax_f64_z( inVec1.vec, inVec2.vec);
    }

    inline static InaVecSXA IsLowerOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vbrdd_vsvl_z(svcmple_f64( inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsLower(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vbrdd_vsvl_z(svcmplt_f64( inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsGreaterOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vbrdd_vsvl_z(svcmpge_f64( inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsGreater(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vbrdd_vsvl_z(svcmpgt_f64( inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vbrdd_vsvl_z(svcmpeq_f64( inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsNotEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vbrdd_vsvl_z(svcmpne_f64( inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecMaskSXA<double> IsLowerOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmple_f64( inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<double> IsLowerMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmplt_f64( inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<double> IsGreaterOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpge_f64( inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<double> IsGreaterMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpgt_f64( inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<double> IsEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpeq_f64( inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<double> IsNotEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpne_f64( inVec1.vec, inVec2.vec);
    }

    inline static InaVecSXA BitsAnd(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f64_u64(svand_u64_z( svreinterpret_u64_f64(inVec1.vec), svreinterpret_u64_f64(inVec2.vec)));
    }

    inline static InaVecSXA BitsNotAnd(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f64_u64(svbic_u64_z( svreinterpret_u64_f64(inVec2.vec), svreinterpret_u64_f64(inVec1.vec)));
    }

    inline static InaVecSXA BitsOr(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f64_u64(svorr_u64_z( svreinterpret_u64_f64(inVec1.vec), svreinterpret_u64_f64(inVec2.vec)));
    }

    inline static InaVecSXA BitsXor(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f64_u64(sveor_u64_z( svreinterpret_u64_f64(inVec1.vec), svreinterpret_u64_f64(inVec2.vec)));
    }

    inline static  const char* GetName() {
        return "InaVecSXA<double>";
    }

    inline static  InaIfElse< InaVecSXA<double> >::ThenClass If(const InaVecMaskSXA<double>& inTest) {
        return InaIfElse< InaVecSXA<double> >::IfClass().If(inTest);
    }

    inline static InaVecSXA IfElse(const InaVecMaskSXA<double>& inMask, const InaVecSXA& inIfTrue, const InaVecSXA& inIfFalse) {
        return svsel_f64(__vr(inMask), inIfTrue.vec, inIfFalse.vec);
    }

    inline static InaVecSXA IfTrue(const InaVecMaskSXA<double>& inMask, const InaVecSXA& inIfTrue) {
        return svsel_f64(__vr(inMask), inIfTrue.vec, _vel_vbrdd_vsvl(0));
    }

    inline static InaVecSXA IfFalse(const InaVecMaskSXA<double>& inMask, const InaVecSXA& inIfFalse) {
        return svsel_f64(__vr(inMask), _vel_vbrdd_vsvl(0), inIfFalse.vec);
    }

    // Inner operators
    inline InaVecSXA<double>& operator+=(const InaVecSXA<double>& inVec){
        vec = _vel_vfaddd_vvvl(vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<double>& operator-=(const InaVecSXA<double>& inVec){
        vec = _vel_vfsubd_vvvl(vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<double>& operator/=(const InaVecSXA<double>& inVec){
        vec = _vel_vfdivd_vvvl(vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<double>& operator*=(const InaVecSXA<double>& inVec){
        vec = _vel_vfmuld_vvvl(vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<double> operator-() const {
        return svneg_f64_z( vec);
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
    return _vel_vfaddd_vvvl(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSXA<double> operator-(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return _vel_vfsubd_vvvl(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSXA<double> operator/(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return _vel_vfdivd_vvvl(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSXA<double> operator*(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return _vel_vfmuld_vvvl(inVec1.getVec(), inVec2.getVec());
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
