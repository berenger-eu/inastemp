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
#include <cmath>
#include <initializer_list>
#include <limits>

// Forward declarations
template <class RealType>
class InaVecMaskRISCV;

template <class RealType>
class InaVecRISCV;


// Mask type
template <>
class InaVecMaskRISCV<float> {
    unsigned char maskData[2048/sizeof(unsigned char)];
    vbool4_t& mask;
public:
    // Classic constructors
    // inline InaVecMaskRISCV() {
    //     // mask = vmxor_mm_b4(mask,mask);
    //     mask = vmset_m_b4();
    // }

    // inline InaVecMaskRISCV(const InaVecMaskRISCV& inMask){
    //     mask = inMask.mask;
    // }

    inline InaVecMaskRISCV() : mask(*reinterpret_cast<vbool4_t*>(maskData)){}

    inline InaVecMaskRISCV(const InaVecMaskRISCV&) = default;

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
        // return vmnot_mm_b4(mask);
        return mask;
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
class InaVecRISCV<float> {
protected:
    unsigned char vecData[2048/sizeof(unsigned char)];
    vfloat32m8_t& vec;

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
        return 64;
    }

    static constexpr bool IsRealFma(){
        return true;
    }

    // inline InaVecRISCV() {
    //     vec = vundefined_f32m8();
    // }
    // inline InaVecRISCV(const InaVecRISCV& inVec){
    //     vec = inVec.vec;
    // }

    inline InaVecRISCV() : vec(*reinterpret_cast<vfloat32m8_t*>(vecData)) {}

    inline InaVecRISCV(const InaVecRISCV&) = default;

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
        float32_t values [64];
        vsetvl_e32m8(64);
        for (int i=0;i<GetVecLength();i++){
          values[i] = val;
        }
        vec = vle32_v_f32m8(values);
    }

    inline InaVecRISCV& operator=(const float val){
        uint32_t tabIndex [64];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vec = vlxei32_v_f32m8(&val,index);
        return *this;
    }

    inline void setFromScalar(const float val){
        uint32_t tabIndex [64];
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

    inline InaVecRISCV& setFromIndirectArray(const float values[], const int inIndirection[]) {
      float32_t result[64];
      for(int i=0; i<64; i++){
          result[i] = values[inIndirection[i]];
      }
      vec = vle32_v_f32m8(result);
      return *this;
    }

    inline InaVecRISCV& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        float32_t result[64];
        for(int i=0; i<64; i++){
            result[i] = inArray[inIndirection1[i] * inLeadingDimension + inIndirection2[i]];
        }
        vec = vle32_v_f32m8(result);
        return *this;
    }

    // Move back to array
    inline void storeInArray(float ptr[]) const {
         vse32_v_f32m8(ptr,vec);
    }

    inline void storeInAlignedArray(float ptr[]) const {
        vse32_v_f32m8(ptr,vec);
    }

    // Acce to individual values
    inline float at(const int index) const {
        vsetvl_e32m8(64);
        float32_t tab[64];
        vse32_v_f32m8(tab,vec);
        return tab[index];
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
      float sum = at(0);
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

        const float32_t COEFF_LOG2E = float(InaFastExp::CoeffLog2E());
        const float32_t COEFF_A     = float(InaFastExp::CoeffA64());
        const float32_t COEFF_B     = float(InaFastExp::CoeffB64());
        const float32_t COEFF_P5_X  = float(InaFastExp::GetCoefficient9_8());
        const float32_t COEFF_P5_Y  = float(InaFastExp::GetCoefficient9_7());
        const float32_t COEFF_P5_Z  = float(InaFastExp::GetCoefficient9_6());
        const float32_t COEFF_P5_A  = float(InaFastExp::GetCoefficient9_5());
        const float32_t COEFF_P5_B  = float(InaFastExp::GetCoefficient9_4());
        const float32_t COEFF_P5_C  = float(InaFastExp::GetCoefficient9_3());
        const float32_t COEFF_P5_D  = float(InaFastExp::GetCoefficient9_2());
        const float32_t COEFF_P5_E  = float(InaFastExp::GetCoefficient9_1());
        const float32_t COEFF_P5_F  = float(InaFastExp::GetCoefficient9_0());

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

        const float32_t COEFF_LOG2E = float(InaFastExp::CoeffLog2E());
        const float32_t COEFF_A     = float(InaFastExp::CoeffA64());
        const float32_t COEFF_B     = float(InaFastExp::CoeffB64());
        const float32_t COEFF_P5_C  = float(InaFastExp::GetCoefficient4_3());
        const float32_t COEFF_P5_D  = float(InaFastExp::GetCoefficient4_2());
        const float32_t COEFF_P5_E  = float(InaFastExp::GetCoefficient4_1());
        const float32_t COEFF_P5_F  = float(InaFastExp::GetCoefficient4_0());

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
        uint32_t tabIndex [64];
        for (int i=0;i<int(GetVecLength());i++)
            tabIndex[i] = 0;
        const vuint32m8_t index = vle32_v_u32m8(tabIndex);
        const float positif = 1.0;
        const vfloat32m8_t one = vlxei32_v_f32m8(&positif,index);
        return  vfsqrt_v_f32m8(vfdiv_vv_f32m8(one, vec));
    }

    inline InaVecRISCV abs() const {
        // return vfabs(vec);
        float32_t values [64];
        vsetvl_e32m8(64);
        for (int i=0;i<GetVecLength();i++){
          values[i] = 0;
        }
        vfloat32m8_t vzero = vle32_v_f32m8(values);
        vfloat32m8_t min = vfsub_vv_f32m8(vzero,vec);
        vec = vfmax_vv_f32m8(vec,min);
        return vec;
    }

    inline InaVecRISCV floor() const {
        float32_t positif = 1.0;
        vfloat32m8_t valuesInIntervals = vfmin_vf_f32m8(
                                    vfmax_vf_f32m8( vec,double(std::numeric_limits<long int>::min())),
                                    double(std::numeric_limits<long int>::max()));
        vint32m8_t vecConvLongInt = vfcvt_rtz_x_f_v_i32m8(valuesInIntervals);
        vfloat32m8_t vecConvLongIntDouble = vfcvt_f_x_v_f32m8(vecConvLongInt);
        vbool4_t maskPositive = vmflt_vf_f32m8_b4(vec,0);
        return vmerge_vvm_f32m8(maskPositive, vecConvLongIntDouble, vfsub_vf_f32m8(vecConvLongIntDouble,positif));
    }

    inline InaVecRISCV signOf() const {
        float32_t values [64];
        vsetvl_e32m8(64);
        for (int i=0;i<GetVecLength();i++){
          values[i] = 0;
        }
        vfloat32m8_t vzero = vle32_v_f32m8(values);
        vbool4_t maskPositive = vmflt_vv_f32m8_b4(vzero,vec);
        vbool4_t maskNegative = vmfgt_vv_f32m8_b4(vzero,vec);

        float32_t tabPositif[64];
        float32_t tabNegatif[64];
        for (int i=0;i<GetVecLength();i++){
            tabPositif[i] = 1;
            tabNegatif[i] = -1;
        }
        vfloat32m8_t vpositif = vle32_v_f32m8(tabPositif);
        vfloat32m8_t vnegatif = vle32_v_f32m8(tabNegatif);

        return vmerge_vvm_f32m8(maskNegative,
            vmerge_vvm_f32m8(maskPositive,vzero,vpositif),vnegatif);
    }

    inline InaVecRISCV isPositive() const {
        float32_t values [64];
        vsetvl_e32m8(64);
        for (int i=0;i<GetVecLength();i++){
          values[i] = 0;
        }
        vfloat32m8_t vzero = vle32_v_f32m8(values);
        vbool4_t maskPositive = vmfge_vv_f32m8_b4(vec,vzero);

        float32_t tabPositif [64];
        for (int i=0;i<GetVecLength();i++)
            tabPositif[i] = 1;
        vfloat32m8_t vpositif = vle32_v_f32m8(tabPositif);
        return vmerge_vvm_f32m8(maskPositive,vzero,vpositif);
    }

    inline InaVecRISCV isNegative() const {
        float32_t values [64];
        vsetvl_e32m8(64);
        for (int i=0;i<GetVecLength();i++){
          values[i] = 0;
        }
        vfloat32m8_t vzero = vle32_v_f32m8(values);
        vbool4_t maskNegative = vmfle_vv_f32m8_b4(vec,vzero);

        float32_t tabNegatif[64];
        for (int i=0;i<GetVecLength();i++)
            tabNegatif[i] = 1;
        vfloat32m8_t vnegatif = vle32_v_f32m8(tabNegatif);

        return vmerge_vvm_f32m8(maskNegative,vzero,vnegatif);
    }

    inline InaVecRISCV isPositiveStrict() const {
        float32_t values [64];
        vsetvl_e32m8(64);
        for (int i=0;i<GetVecLength();i++){
          values[i] = 0;
        }
        vfloat32m8_t vzero = vle32_v_f32m8(values);
        vbool4_t maskPositive = vmfgt_vv_f32m8_b4(vec,vzero);

        float32_t tabPositif [64];
        for (int i=0;i<GetVecLength();i++)
            tabPositif[i] = 1;
        vfloat32m8_t vpositif = vle32_v_f32m8(tabPositif);
        return vmerge_vvm_f32m8(maskPositive,vzero,vpositif);
    }

    inline InaVecRISCV isNegativeStrict() const {
        float32_t values [64];
        vsetvl_e32m8(64);
        for (int i=0;i<GetVecLength();i++){
          values[i] = 0;
        }
        vfloat32m8_t vzero = vle32_v_f32m8(values);
        vbool4_t maskNegative = vmflt_vv_f32m8_b4(vec,vzero);

        float32_t tabNegatif[64];
        for (int i=0;i<GetVecLength();i++)
            tabNegatif[i] = 1;
        vfloat32m8_t vnegatif = vle32_v_f32m8(tabNegatif);

        return vmerge_vvm_f32m8(maskNegative,vzero,vnegatif);
    }

    inline InaVecRISCV isZero() const {
        float32_t values [64];
        vsetvl_e32m8(64);
        for (int i=0;i<GetVecLength();i++){
          values[i] = 0;
        }
        vfloat32m8_t vzero = vle32_v_f32m8(values);
        vbool4_t maskEqual = vmfeq_vv_f32m8_b4(vec,vzero);

        float32_t tabEqual[64];
        for (int i=0;i<GetVecLength();i++)
            tabEqual[i] = 1;
        vfloat32m8_t vequal = vle32_v_f32m8(tabEqual);

        return vmerge_vvm_f32m8(maskEqual,vzero,vequal);
    }

    inline InaVecRISCV isNotZero() const {
        float32_t values [64];
        vsetvl_e32m8(64);
        for (int i=0;i<GetVecLength();i++){
          values[i] = 0;
        }
        vfloat32m8_t vzero = vle32_v_f32m8(values);
        vbool4_t maskNotEqual = vmfne_vv_f32m8_b4(vec,vzero);

        float32_t tabNotEqual[64];
        for (int i=0;i<GetVecLength();i++)
            tabNotEqual[i] = 1;
        vfloat32m8_t vNotEqual = vle32_v_f32m8(tabNotEqual);

        return vmerge_vvm_f32m8(maskNotEqual,vzero,vNotEqual);
    }

    inline InaVecMaskRISCV<float> isPositiveMask() const {
        return vmfge_vf_f32m8_b4(vec, 0);
    }

    inline InaVecMaskRISCV<float> isNegativeMask() const {
        return vmfle_vf_f32m8_b4(vec,0);
    }

    inline InaVecMaskRISCV<float> isPositiveStrictMask() const {
        return vmfgt_vf_f32m8_b4(vec,0);
    }

    inline InaVecMaskRISCV<float> isNegativeStrictMask() const {
        return vmflt_vf_f32m8_b4(vec,0);
    }

    inline InaVecMaskRISCV<float> isZeroMask() const {
        return vmfeq_vf_f32m8_b4(vec,0);
    }

    inline InaVecMaskRISCV<float> isNotZeroMask() const {
        return vmfne_vf_f32m8_b4(vec,0);
    }

    // Static basic methods
    inline static InaVecRISCV GetZero() {
        const float zero = 0.0;
        uint32_t tabIndex [64];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t vecZero = vlxei32_v_f32m8(&zero,index);
        return InaVecRISCV(vecZero);
    }

    inline static InaVecRISCV GetOne() {
        const float one = 1.0;
        uint32_t tabIndex [64];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t vecOne = vlxei32_v_f32m8(&one,index);
        return InaVecRISCV(vecOne);
    }

    inline static InaVecRISCV Min(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfmin_vv_f32m8( inVec1.vec, inVec2.vec);
    }

    inline static InaVecRISCV Max(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfmax_vv_f32m8( inVec1.vec, inVec2.vec);
    }

    inline static InaVecRISCV IsLowerOrEqual(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool4_t mask = vmfle_vv_f32m8_b4(inVec1.vec,inVec2.vec);
        const float zero = 0.0;
        const float one = 1.0;

        uint32_t tabIndex [64];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t vecZero = vlxei32_v_f32m8(&zero,index);
        vfloat32m8_t vecOne = vlxei32_v_f32m8(&one,index);

        return vmerge_vvm_f32m8(mask,vecZero,vecOne);
    }

    inline static InaVecRISCV IsLower(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool4_t mask = vmflt_vv_f32m8_b4(inVec1.vec,inVec2.vec);
        const float zero = 0.0;
        const float one = 1.0;

        uint32_t tabIndex [64];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t vecZero = vlxei32_v_f32m8(&zero,index);
        vfloat32m8_t vecOne = vlxei32_v_f32m8(&one,index);

        return vmerge_vvm_f32m8(mask,vecZero,vecOne);
    }

    inline static InaVecRISCV IsGreaterOrEqual(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool4_t mask = vmfge_vv_f32m8_b4(inVec1.vec,inVec2.vec);
        const float zero = 0.0;
        const float one = 1.0;

        uint32_t tabIndex [64];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t vecZero = vlxei32_v_f32m8(&zero,index);
        vfloat32m8_t vecOne = vlxei32_v_f32m8(&one,index);

        return vmerge_vvm_f32m8(mask,vecZero,vecOne);
    }

    inline static InaVecRISCV IsGreater(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool4_t mask = vmfgt_vv_f32m8_b4(inVec1.vec,inVec2.vec);
        const float zero = 0.0;
        const float one = 1.0;

        uint32_t tabIndex [64];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t vecZero = vlxei32_v_f32m8(&zero,index);
        vfloat32m8_t vecOne = vlxei32_v_f32m8(&one,index);

        return vmerge_vvm_f32m8(mask,vecZero,vecOne);
    }

    inline static InaVecRISCV IsEqual(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool4_t mask = vmfeq_vv_f32m8_b4(inVec1.vec,inVec2.vec);
        const float zero = 0.0;
        const float one = 1.0;

        uint32_t tabIndex [64];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t vecZero = vlxei32_v_f32m8(&zero,index);
        vfloat32m8_t vecOne = vlxei32_v_f32m8(&one,index);

        return vmerge_vvm_f32m8(mask,vecZero,vecOne);
    }

    inline static InaVecRISCV IsNotEqual(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool4_t mask = vmfne_vv_f32m8_b4(inVec1.vec,inVec2.vec);
        const float zero = 0.0;
        const float one = 1.0;

        uint32_t tabIndex [64];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t vecZero = vlxei32_v_f32m8(&zero,index);
        vfloat32m8_t vecOne = vlxei32_v_f32m8(&one,index);

        return vmerge_vvm_f32m8(mask,vecZero,vecOne);
    }

    inline static InaVecMaskRISCV<float> IsLowerOrEqualMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return  vmfle_vv_f32m8_b4(inVec1.vec,inVec2.vec);
    }

    inline static InaVecMaskRISCV<float> IsLowerMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vmflt_vv_f32m8_b4(inVec1.vec,inVec2.vec);
    }

    inline static InaVecMaskRISCV<float> IsGreaterOrEqualMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vmfge_vv_f32m8_b4(inVec1.vec,inVec2.vec);
    }

    inline static InaVecMaskRISCV<float> IsGreaterMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vmfgt_vv_f32m8_b4(inVec1.vec,inVec2.vec);
    }

    inline static InaVecMaskRISCV<float> IsEqualMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vmfeq_vv_f32m8_b4(inVec1.vec,inVec2.vec);
    }

    inline static InaVecMaskRISCV<float> IsNotEqualMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vmfne_vv_f32m8_b4(inVec1.vec,inVec2.vec);
    }

    inline static InaVecRISCV BitsAnd(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfcvt_f_x_v_f32m8(vand_vv_i32m8(vfcvt_rtz_x_f_v_i32m8(inVec1.vec),vfcvt_rtz_x_f_v_i32m8(inVec2.vec)));
    }

    inline static InaVecRISCV BitsNotAnd(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfcvt_f_x_v_f32m8(vnot_v_i32m8(vand_vv_i32m8(vfcvt_rtz_x_f_v_i32m8(inVec1.vec),vfcvt_rtz_x_f_v_i32m8(inVec2.vec))));
    }

    inline static InaVecRISCV BitsOr(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfcvt_f_x_v_f32m8(vor_vv_i32m8(vfcvt_rtz_x_f_v_i32m8(inVec1.vec),vfcvt_rtz_x_f_v_i32m8(inVec2.vec)));
    }

    inline static InaVecRISCV BitsXor(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfcvt_f_x_v_f32m8(vxor_vv_i32m8(vfcvt_rtz_x_f_v_i32m8(inVec1.vec),vfcvt_rtz_x_f_v_i32m8(inVec2.vec)));
    }

    inline static  const char* GetName() {
        return "InaVecRISCV<float>";
    }

    inline static  InaIfElse< InaVecRISCV<float> >::ThenClass If(const InaVecMaskRISCV<float>& inTest) {
       return InaIfElse< InaVecRISCV<float> >::IfClass().If(inTest);
    }

    inline static InaVecRISCV IfElse(const InaVecMaskRISCV<float>& inMask, const InaVecRISCV& inIfTrue, const InaVecRISCV& inIfFalse) {
        return vmerge_vvm_f32m8(vbool4_t(inMask),inIfTrue.vec,inIfFalse.vec);
    }

    inline static InaVecRISCV IfTrue(const InaVecMaskRISCV<float>& inMask, const InaVecRISCV& inIfTrue) {
        return vfmerge_vfm_f32m8(vbool4_t(inMask),inIfTrue.vec,0);
    }

    inline static InaVecRISCV IfFalse(const InaVecMaskRISCV<float>& inMask, const InaVecRISCV& inIfFalse) {
        const float zero = 0.0;
        uint32_t tabIndex [64];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint32m8_t index = vle32_v_u32m8(tabIndex);
        vfloat32m8_t vecZero = vlxei32_v_f32m8(&zero,index);
        return vmerge_vvm_f32m8(vbool4_t(inMask),vecZero,inIfFalse.vec);
    }

    // Inner operators
    inline InaVecRISCV<float>& operator+=(const InaVecRISCV<float>& inVec){
        vec = vfadd_vv_f32m8(vec,inVec.vec);
        return *this;
    }

    inline InaVecRISCV<float>& operator-=(const InaVecRISCV<float>& inVec){
        vec = vfsub_vv_f32m8(vec,inVec.vec);
        return *this;
    }

    inline InaVecRISCV<float>& operator/=(const InaVecRISCV<float>& inVec){
        vec = vfdiv_vv_f32m8(vec,inVec.vec);
        return *this;
    }

    inline InaVecRISCV<float>& operator*=(const InaVecRISCV<float>& inVec){
        vec = vfmul_vv_f32m8(vec,inVec.vec);
        return *this;
    }

    inline InaVecRISCV<float> operator-() const {
        // return vfneg_v_f32m8(vec);
        return vec;
    }

    inline InaVecRISCV<float> pow(std::size_t power) const{
        return InaUtils::FastPow<InaVecRISCV<float>>(*this, power);
    }

    // Multiple sum
    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecRISCV<float>& inVec1,
                                          const InaVecRISCV<float>& inVec2, const InaVecRISCV<float>& inVec3,
                                          const InaVecRISCV<float>& inVec4, const InaVecRISCV<float>& inVec5,
                                          const InaVecRISCV<float>& inVec6, const InaVecRISCV<float>& inVec7,
                                          const InaVecRISCV<float>& inVec8, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1, inVec2, inVec3, inVec4 );
        MultiHorizontalSum(&sumRes[4], inVec5, inVec6, inVec7, inVec8 );

        MultiHorizontalSum(&sumRes[8], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecRISCV<float>& inVec1,
                                          const InaVecRISCV<float>& inVec2, const InaVecRISCV<float>& inVec3,
                                          const InaVecRISCV<float>& inVec4, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1, inVec2 );
        MultiHorizontalSum(&sumRes[2], inVec3, inVec4 );

        MultiHorizontalSum(&sumRes[4], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(float sumRes[], const InaVecRISCV<float>& inVec1,
                                          const InaVecRISCV<float>& inVec2, Args ...args){

        MultiHorizontalSum(&sumRes[0], inVec1);
        MultiHorizontalSum(&sumRes[1], inVec2);

        MultiHorizontalSum(&sumRes[2], args... );
    }

    inline static void MultiHorizontalSum(float sumRes[], const InaVecRISCV<float>& inVec){
        sumRes[0] += inVec.horizontalSum();
    }

    inline static void MultiHorizontalSum(float /*sumRes*/[]){
    }

    inline static InaVecRISCV<float> Fma(const InaVecRISCV<float>& inValAdd, const InaVecRISCV<float>& inValMul1, const InaVecRISCV<float>& inValMul2){
        return vfnmadd_vv_f32m8(inValAdd.vec, inValMul1.vec, inValMul2.vec);
    }
};

// Bits operators
inline InaVecRISCV<float> operator&(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return InaVecRISCV<float>::BitsAnd(inVec1, inVec2);
}

inline InaVecRISCV<float> operator|(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return InaVecRISCV<float>::BitsOr(inVec1, inVec2);
}

inline InaVecRISCV<float> operator^(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return InaVecRISCV<float>::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecRISCV<float> operator+(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return vfadd_vv_f32m8(inVec1.getVec(), inVec2.getVec());
}

inline InaVecRISCV<float> operator-(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return vfsub_vv_f32m8(inVec1.getVec(), inVec2.getVec());
}

inline InaVecRISCV<float> operator/(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return vfdiv_vv_f32m8(inVec1.getVec(), inVec2.getVec());
}

inline InaVecRISCV<float> operator*(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return vfmul_vv_f32m8(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskRISCV<float> operator<(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return InaVecRISCV<float>::IsLowerMask(inVec1,inVec2);
}

inline InaVecMaskRISCV<float> operator<=(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return InaVecRISCV<float>::IsLowerOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskRISCV<float> operator>(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return InaVecRISCV<float>::IsGreaterMask(inVec1,inVec2);
}

inline InaVecMaskRISCV<float> operator>=(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return InaVecRISCV<float>::IsGreaterOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskRISCV<float> operator==(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return InaVecRISCV<float>::IsEqualMask(inVec1,inVec2);
}

inline InaVecMaskRISCV<float> operator!=(const InaVecRISCV<float>& inVec1, const InaVecRISCV<float>& inVec2){
    return InaVecRISCV<float>::IsNotEqualMask(inVec1,inVec2);
}


#endif
