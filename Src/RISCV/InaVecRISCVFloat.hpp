///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECRISCVFLOAT_HPP
#define INAVECRISCVFLOAT_HPP

#include "InastempGlobal.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_RISCV
#error InaVecRISCV<float> is included but RISCV is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <riscv_vector.h>
#include <algorithm>
// #include <velintrin.h>
// #include <cmath>
// #include <initializer_list>
// #include <limits>

// Forward declarations
template <class RealType>
class InaVecMaskRISCV;

template <class RealType>
class InaVecRISCV;


// Mask type
template <>
class alignas(2048) InaVecMaskRISCV<float> {
    vbool4_t mask;

public:
    // Classic constructors
    inline InaVecMaskRISCV() {
        mask = vmxor_mm_b4(mask,mask);
    }

    inline InaVecMaskRISCV(const InaVecMaskRISCV& inMask){
        mask = inMask.mask;
    }

    inline InaVecMaskRISCV& operator=(const InaVecMaskRISCV& inMask){
        mask = inMask.mask;
        return *this;
    }

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskRISCV(const vbool4_t inMask)
        : InaVecMaskRISCV() {
        mask = (inMask);
    }

    inline InaVecMaskRISCV& operator=(const vbool4_t inMask){
        mask = inMask;
        return (*this);
    }

    inline explicit operator vbool4_t() const{
        return mask;
    }

    inline vbool4_t getMask() const{
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskRISCV(const bool inBool) : InaVecMaskRISCV() {
        mask = (inBool? vmxnor_mm_b4(mask, mask) : vmxor_mm_b4(mask, mask));
    }

    inline InaVecMaskRISCV& operator=(const bool inBool){
        mask = (inBool? vmxnor_mm_b4(mask, mask) : vmxor_mm_b4(mask, mask));
        return (*this);
    }

    // Binary methods
    inline InaVecMaskRISCV Not() const{
        return vmnot_mm_b4(mask);
    }

    inline bool isAllTrue() const{
        return vpopc_m_b4(mask) == 256;
    }

    inline bool isAllFalse() const{
        return vpopc_m_b4(mask) == 0;
    }
    
    // Float args methods
    inline static InaVecMaskRISCV And(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vmand_mm_b4(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskRISCV NotAnd(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vmnand_mm_b4(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskRISCV Or(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vmor_mm_b4(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskRISCV Xor(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vmxor_mm_b4(inMask1.mask,inMask2.mask);
    }

    inline static bool IsEqual(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vpopc_m_b4(vmxor_mm_b4(inMask1.mask,inMask2.mask)) == 0;
    }

    inline static bool IsNotEqual(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vpopc_m_b4(vmxor_mm_b4(inMask1.mask,inMask2.mask)) != 0;
    }
};

// Mask must have operators
inline InaVecMaskRISCV<float> operator&(const InaVecMaskRISCV<float>& inMask1, const InaVecMaskRISCV<float>& inMask2){
    return InaVecMaskRISCV<float>::And(inMask1, inMask2);
}

inline InaVecMaskRISCV<float> operator|(const InaVecMaskRISCV<float>& inMask1, const InaVecMaskRISCV<float>& inMask2){
    return InaVecMaskRISCV<float>::Or(inMask1, inMask2);
}

inline InaVecMaskRISCV<float> operator^(const InaVecMaskRISCV<float>& inMask1, const InaVecMaskRISCV<float>& inMask2){
    return InaVecMaskRISCV<float>::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskRISCV<float>& inMask1, const InaVecMaskRISCV<float>& inMask2){
    return InaVecMaskRISCV<float>::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskRISCV<float>& inMask1, const InaVecMaskRISCV<float>& inMask2){
    return InaVecMaskRISCV<float>::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(2048) InaVecRISCV<float> {
protected:
    vfloat32m8_t vec;

    // static __vm256 _vel_vfmklgt_mvl_256(__vr vz, int vl){
    //     return _vel_vfmklgt_mvl(vz, vl); // __vm256();
    // }

public:
    using VecRawType           = vfloat32m8_t;
    using MaskType             = InaVecMaskRISCV<float>;
    using RealType             = float;
    static const int Alignement= 1;
    static const bool IsOfFixedSize = true;

    static constexpr int GetVecLength(){
        return  vsetvlmax_e32m8();
    }

    static constexpr bool IsRealFma(){
        return true;
    }

    inline InaVecRISCV() {
        vec = vundefined_f32m8();
    }
    inline InaVecRISCV(const InaVecRISCV& inVec){
        vec = inVec.vec;
    }

    inline InaVecRISCV& operator=(const InaVecRISCV& inVec){
        vec = inVec.vec;
        return *this;
    }

    // Constructor from raw type
    inline /*not explicit*/ InaVecRISCV(const vfloat32m8_t inVec)
        : InaVecRISCV() {
        vec = (inVec);
    }

    inline InaVecRISCV& operator=(const vfloat32m8_t inVec){
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const vfloat32m8_t inVec){
        vec = inVec;
    }

    inline explicit operator vfloat32m8_t() const{
        return vec;
    }

    inline vfloat32m8_t getVec() const{
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecRISCV(const float val)
        : InaVecRISCV() {
        uint32_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vec = vlxei32_v_f32m8(&val,index);
    }

    inline InaVecRISCV& operator=(const float val){
        uint32_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vec = vlxei32_v_f32m8(&val,index);
        return *this;
    }

    inline void setFromScalar(const float val){
        uint32_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vec = vlxei32_v_f32m8(&val,index);
    }

    // Constructor from vec
    inline InaVecRISCV(const std::initializer_list<float> lst)
        : InaVecRISCV(lst.begin()){
    }

    inline explicit InaVecRISCV(const float ptr[])
        : InaVecRISCV() {
        vec = vle32_v_f32m8(ptr);
    }

    inline InaVecRISCV& setFromArray(const float ptr[]){
      vec = vle32_v_f32m8(ptr);
      return *this;
    }

    inline InaVecRISCV& setFromAlignedArray(const float ptr[]){
      vec = vle32_v_f32m8(ptr);
      return *this;
    }

    // inline InaVecSXA& setFromIndirectArray(const float values[], const unsigned long int inIndirection[]) {
    //     __vr offset = _vel_vld_vssl(8, inIndirection, 256);
    //     __vr address = _vel_vsfa_vvssl(offset, 2, reinterpret_cast<unsigned long>(values), 256);
    //     vec = _vel_vgtu_vvssl(address, 0, 0, 256);
    // 
    //     return *this;
    // }

    inline InaVecRISCV& setFromIndirectArray(const float values[], const int inIndirection[]) {
      vec = vlxei32_v_f32m8(values,Indirection);
      return *this;
    }

    // inline InaVecSXA& setFromIndirect2DArray(const float inArray[], const long int inIndirection1[],
    //                              const int inLeadingDimension, const long int inIndirection2[]){
    //     __vr offset = _vel_vaddsl_vvvl(_vel_vld_vssl(8, inIndirection2, 256),
    //                  _vel_vmulul_vvvl(_vel_vbrdl_vsl(inLeadingDimension, 256),
    //                                   _vel_vld_vssl(8, inIndirection1, 256),
    //                                   256),256);
    //     __vr address = _vel_vsfa_vvssl(offset, 2, reinterpret_cast<unsigned long>(inArray), 256);
    //     vec = _vel_vgtu_vvssl(address, 0, 0, 256);
    //     return *this;
    // 
    // }

    inline InaVecRISCV& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        vec = vlxei32_v_f32m8(inArray,vfadd_vv_f32m8(
          vfmul_vf_f32m8(vle32ff_v_f32m8(inIndirection1,32),inLeadingDimension),vle32ff_v_f32m8(inIndirection2,32)
        ));

        return *this;
    }

    // Move back to array
    inline void storeInArray(float ptr[]) const {
         vse32_v_f32m8(vec,ptr);
    }

    inline void storeInAlignedArray(float ptr[]) const {
        vse32_v_f32m8(vec,ptr);
    }

    // Acce to individual values
    inline float at(const int index) const {
        return vec[index];
    }

    // Horizontal operation
    inline float horizontalSum() const {
      float sum = at(0);
      for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
          sum+= at(idx);
      }
      return sum;
    }

    inline float horizontalMul() const {
      double sum = at(0);
      for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
          sum*= at(idx);
      }
      return sum;
    }

    inline float minInVec() const {
      float min = at(0);
      for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
          min = std::min(min,at(idx));
      }
      return min;
    }

    inline float maxInVec() const {
      float max = at(0);
      for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
          max = std::max(max,at(idx));
      }
      return max;
    }

    inline InaVecRISCV sqrt() const {
        return vfsqrt_v_f32m8(vec);
    }

    inline InaVecRISCV exp() const {

        const float32m8_t COEFF_LOG2E = InaFastExp::CoeffLog2E();
        const float32m8_t COEFF_A     = InaFastExp::CoeffA64();
        const float32m8_t COEFF_B     = InaFastExp::CoeffB64();
        const float32m8_t COEFF_P5_X  = InaFastExp::GetCoefficient9_8();
        const float32m8_t COEFF_P5_Y  = InaFastExp::GetCoefficient9_7();
        const float32m8_t COEFF_P5_Z  = InaFastExp::GetCoefficient9_6();
        const float32m8_t COEFF_P5_A  = InaFastExp::GetCoefficient9_5();
        const float32m8_t COEFF_P5_B  = InaFastExp::GetCoefficient9_4();
        const float32m8_t COEFF_P5_C  = InaFastExp::GetCoefficient9_3();
        const float32m8_t COEFF_P5_D  = InaFastExp::GetCoefficient9_2();
        const float32m8_t COEFF_P5_E  = InaFastExp::GetCoefficient9_1();
        const float32m8_t COEFF_P5_F  = InaFastExp::GetCoefficient9_0();

        vfloat32m8_t x = vfmul_vf_f32m8(vec, COEFF_LOG2E);

        const vfloat32m8_t fractional_part = vfsub_vv_f32m8(x, InaVecRISCV(x).floor().vec);

        vfloat32m8_t factor = vfadd_vf_f32m8( vfmul_vv_f32m8( vfadd_vf_f32m8(
                         vfmul_vv_f32m8( vfadd_vf_f32m8( vfmul_vv_f32m8( vfadd_vf_f32m8(
                         vfmul_vv_f32m8( vfadd_vf_f32m8( vfmul_vv_f32m8( vfadd_vf_f32m8(
                         vfmul_vv_f32m8( vfadd_vf_f32m8( vfmul_vv_f32m8( vfadd_vf_f32m8(
                         vfmul_vf_f32m8(fractional_part,COEFF_P5_X), COEFF_P5_Y),fractional_part),
                         COEFF_P5_Z),fractional_part), COEFF_P5_A), fractional_part),
                         COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part),
                         COEFF_P5_F);

        x = vfsub_vv_f32m8(x,factor);

        x = vfadd_vf_f32m8(vfmul_vf_f32m8(x, COEFF_A), COEFF_B);

        vint32m8_t castedInteger = vfcvt_rtz_x_f_v_i32m8(x);

        return vfcvt_f_x_v_f32m8(castedInteger);
    }

    inline InaVecRISCV expLowAcc() const {

        const float32m8_t COEFF_LOG2E = InaFastExp::CoeffLog2E();
        const float32m8_t COEFF_A     = InaFastExp::CoeffA64();
        const float32m8_t COEFF_B     = InaFastExp::CoeffB64();
        const float32m8_t COEFF_P5_C  = InaFastExp::GetCoefficient4_3();
        const float32m8_t COEFF_P5_D  = InaFastExp::GetCoefficient4_2();
        const float32m8_t COEFF_P5_E  = InaFastExp::GetCoefficient4_1();
        const float32m8_t COEFF_P5_F  = InaFastExp::GetCoefficient4_0();

        vfloat32m8_t x = vfmul_vf_f32m8(vec, COEFF_LOG2E);

        const vfloat32m8_t fractional_part = vfsub_vv_f32m8(x, InaVecRISCV(x).floor().vec);

        vfloat32m8_t factor = vfadd_vf_f32m8(vfmul_vv_f32m8(vfadd_vf_f32m8(
                         vfmul_vv_f32m8(vfadd_vf_f32m8(vfmul_vf_f32m8(
                                         fractional_part, COEFF_P5_C),
                                         COEFF_P5_D), fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = vfsub_vv_f32m8(x,factor);

        x = vfadd_vf_f32m8(vfmul_vf_f32m8(x, COEFF_A), COEFF_B);

        vint32m8_t castedInteger = vfcvt_rtz_x_f_v_i32m8(x);

        return  vfcvt_f_x_v_f32m8(castedInteger);
    }

    inline InaVecRISCV rsqrt() const {
      // We can use vfrsqrt7_v_f32m8_t
        uint32_t tabIndex [GetVecLength()];
        for (int i=0;i<int(GetVecLength());i++)
            tabIndex[i] = 0;
        const vuint32m8_t index = vle32_v_u32m8(tabIndex);

        const vfloat32m8_t one = vlxei32_v_f32m8(1.0,index);
        return  vfsqrt_v_f32m8(vfdiv_vv_f32m8(one, vec));
    }

    inline InaVecRISCV abs() const {
        return vfabs_v_f32m8(vec);
    }

    inline InaVecRISCV floor() const {
        vfloat32m8_t valuesInIntervals = vfmin_vf_f32m8(
                                    vfmax_vf_f32m8( vec,double(std::numeric_limits<long int>::min())),
                                    double(std::numeric_limits<long int>::max()));
        vint32m8_t vecConvLongInt = vfcvt_rtz_x_f_v_i32m8(valuesInIntervals);
        vfloat32m8_t vecConvLongIntDouble = vfcvt_f_x_v_f32m8(vecConvLongInt);
        vbool4_t maskPositive = vmflt_vf_f32m8_b4(vec,0);
        return vmerge_vvm_f32m8(maskPositive, vecConvLongIntDouble, vfsub_vv_f32m8(vecConvLongIntDouble,1.0));
    }

    inline InaVecRISCV signOf() const {
        vbool4_t maskPositive = vmfgt_vf_f32m8_b4(vec,0);
        vbool4_t maskNegative = vmflt_vf_f32m8_b4(vec,0);
        const float positif = 1.0;
        const float negatif = -1.0;

        uint32_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);

        vfloat32m8_t signPositive = vlxei32_v_f32m8(&positif,index);
        vfloat32m8_t signNegative = vlxei32_v_f32m8(&negatif,index);

        return vmerge_vvm_f32m8(maskNegative,signNegative,
            vfmerge_vfm_f32m8(maskPositive,signPositive,0));
    }

    inline InaVecRISCV isPositive() const {
        vbool4_t maskPositive = vmfge_vf_f32m8_b4(vec,0);
        const float positif = 1.0;

        uint32_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t signPositive = vlxei32_v_f32m8(1,index);

        return vfmerge_vfm_f32m8(maskPositive,signPositive,0);
    }

    inline InaVecRISCV isNegative() const {
        vbool4_t maskNegative = vmfle_vf_f32m8_b4(vec,0);
        const float negatif = -1.0;

        uint32_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t signNegative = vlxei32_v_f32m8(&negatif,index);

        return vmerge_vfm_f32m8(maskNegative,signNegative,0);
    }
// TODO continue 
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
        return _vel_vfmads_vvvvl(inValAdd.vec, inValMul1.vec, inValMul2.vec, 256);
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
