///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECRISCVDOUBLE_HPP
#define INAVECRISCVDOUBLE_HPP

#include "InastempGlobal.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_RISCV
#error InaVecRISCV<double> is included but RISCV is not enable in the configuration
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
class alignas(2048) InaVecMaskRISCV<double> {
    vbool8_t mask;

public:
    // Classic constructors
    inline InaVecMaskRISCV() {
        mask = vmxor_mm_b8(mask,mask);
    }

    inline InaVecMaskRISCV(const InaVecMaskRISCV& inMask){
        mask = inMask.mask;
    }

    inline InaVecMaskRISCV& operator=(const InaVecMaskRISCV& inMask){
        mask = inMask.mask;
        return *this;
    }

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskRISCV(const vbool8_t inMask)
        : InaVecMaskRISCV() {
        mask = (inMask);
    }

    inline InaVecMaskRISCV& operator=(const vbool8_t inMask){
        mask = inMask;
        return (*this);
    }

    inline explicit operator vbool8_t() const{
        return mask;
    }

    inline vbool8_t getMask() const{
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskRISCV(const bool inBool) : InaVecMaskRISCV() {
        mask = (inBool? vmxnor_mm_b8(mask, mask) : vmxor_mm_b8(mask, mask));
    }

    inline InaVecMaskRISCV& operator=(const bool inBool){
        mask = (inBool? vmxnor_mm_b8(mask, mask) : vmxor_mm_b8(mask, mask));
        return (*this);
    }

    // Binary methods
    inline InaVecMaskRISCV Not() const{
        return vmnot_mm_b8(mask);
    }

    inline bool isAllTrue() const{
        return vpopc_m_b8(mask) == 256;
    }

    inline bool isAllFalse() const{
        return vpopc_m_b8(mask) == 0;
    }

    // Double args methods
    inline static InaVecMaskRISCV And(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vmand_mm_b8(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskRISCV NotAnd(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vmnand_mm_b8(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskRISCV Or(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vmor_mm_b8(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskRISCV Xor(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vmxor_mm_b8(inMask1.mask,inMask2.mask);
    }

    inline static bool IsEqual(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vpopc_m_b8(vmxor_mm_b8(inMask1.mask,inMask2.mask)) == 0;
    }

    inline static bool IsNotEqual(const InaVecMaskRISCV& inMask1, const InaVecMaskRISCV& inMask2){
        return vpopc_m_b8(vmxor_mm_b8(inMask1.mask,inMask2.mask)) != 0;
    }
};

// Mask must have operators
inline InaVecMaskRISCV<double> operator&(const InaVecMaskRISCV<double>& inMask1, const InaVecMaskRISCV<double>& inMask2){
    return InaVecMaskRISCV<double>::And(inMask1, inMask2);
}

inline InaVecMaskRISCV<double> operator|(const InaVecMaskRISCV<double>& inMask1, const InaVecMaskRISCV<double>& inMask2){
    return InaVecMaskRISCV<double>::Or(inMask1, inMask2);
}

inline InaVecMaskRISCV<double> operator^(const InaVecMaskRISCV<double>& inMask1, const InaVecMaskRISCV<double>& inMask2){
    return InaVecMaskRISCV<double>::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskRISCV<double>& inMask1, const InaVecMaskRISCV<double>& inMask2){
    return InaVecMaskRISCV<double>::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskRISCV<double>& inMask1, const InaVecMaskRISCV<double>& inMask2){
    return InaVecMaskRISCV<double>::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(2048) InaVecRISCV<double> {
protected:
    vfloat64m8_t vec;

public:
    using VecRawType           = vfloat64m8_t;
    using MaskType             = InaVecMaskRISCV<double>;
    using RealType             = double;
    static const int Alignement= 1;
    static const bool IsOfFixedSize = true;

    static constexpr int GetVecLength(){
        return vsetvlmax_e64m8();
    }

    static constexpr bool IsRealFma(){
        return true;
    }

    inline InaVecRISCV() {
        vec = vundefined_f64m8();
    }
    inline InaVecRISCV(const InaVecRISCV& inVec){
        vec = inVec.vec;
    }

    inline InaVecRISCV& operator=(const InaVecRISCV& inVec){
        vec = inVec.vec;
        return *this;
    }

    // Constructor from raw type
    inline /*not explicit*/ InaVecRISCV(const vfloat64m8_t inVec)
        : InaVecRISCV() {
        vec = (inVec);
    }

    inline InaVecRISCV& operator=(const vfloat64m8_t inVec){
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const vfloat64m8_t inVec){
        vec = inVec;
    }

    inline explicit operator vfloat64m8_t() const{
        return vec;
    }

    inline vfloat64m8_t getVec() const{
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecRISCV(const double val)
        : InaVecRISCV() {
        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vec = vlxei64_v_f64m8(&val,index);
    }

    inline InaVecRISCV& operator=(const double val){
        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vec = vlxei64_v_f64m8(&val,index);
        return *this;
    }

    inline void setFromScalar(const double val){
        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vec = vlxei64_v_f64m8(&val,index);
    }

    // Constructor from vec
    inline InaVecRISCV(const std::initializer_list<double> lst)
        : InaVecRISCV(lst.begin()){
    }

    inline explicit InaVecRISCV(const double ptr[])
        : InaVecRISCV() {
        vec = vle64_v_f64m8(ptr);
    }

    inline InaVecRISCV& setFromArray(const double ptr[]){
        vec = vle64_v_f64m8(ptr);
        return *this;
    }

    inline InaVecRISCV& setFromAlignedArray(const double ptr[]){
        vec = vle64_v_f64m8(ptr);
        return *this;
    }

    // inline InaVecRISCV& setFromIndirectArray(const double values[], const unsigned long int inIndirection[]) {
    //     vec = vlxei64_v_f64m8(values,Indirection);
    //     return *this;
    // }

    inline InaVecRISCV& setFromIndirectArray(const double values[], const int inIndirection[]) {
        vec = vlxei64_v_f64m8(values,Indirection);
        return *this;
    }
// TODO a faire pour la suite

    // inline InaVecSXA& setFromIndirect2DArray(const double inArray[], const long int inIndirection1[],
    //                              const int inLeadingDimension, const long int inIndirection2[]){
    //     __vr offset = _vel_vaddsl_vvvl(_vel_vld_vssl(8, inIndirection2, 256),
    //                  _vel_vmulul_vvvl(_vel_vbrdl_vsl(inLeadingDimension, 256),
    //                                   _vel_vld_vssl(8, inIndirection1, 256),
    //                                   256),256);
    //     __vr address = _vel_vsfa_vvssl(offset, 3, reinterpret_cast<unsigned long>(inArray), 256);
    //     vec = _vel_vgt_vvssl(address, 0, 0, 256);
    //     return *this;
    // }

    inline InaVecSXA& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        vec = vlxei64_v_f64m8(inArray,vfadd_vv_f64m8(
          vfmul_vf_f64m8(vle64ff_v_f64m8(inIndirection1,32),inLeadingDimension),vle64ff_v_f64m8(inIndirection2,32)
        ));

        return *this;
    }

    // Move back to array
    inline void storeInArray(double ptr[]) const {
       vse64_v_f64m8(vec,ptr);
    }

    inline void storeInAlignedArray(double ptr[]) const {
        vse64_v_f64m8(vec,ptr);
    }

    // Acce to individual values
    inline double at(const int index) const {
        return vec[index];
    }

    // Horizontal operation
    inline double horizontalSum() const {
        double sum = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            sum+= at(idx);
        }
        return sum;
    }

    inline double horizontalMul() const {
        double sum = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            sum *= at(idx);
        }
        return sum;
    }

    inline double minInVec() const {
        double min = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            min = std::min(min,at(idx));
        }
        return min;
    }

    inline double maxInVec() const {
        double max = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            max = std::max(max,at(idx));
        }
        return max;
    }

    inline InaVecRISCV sqrt() const {
        return vfsqrt_v_f64m8(vec);
    }

    inline InaVecRISCV exp() const {

        const float64m8_t COEFF_LOG2E = InaFastExp::CoeffLog2E();
        const float64m8_t COEFF_A     = InaFastExp::CoeffA64();
        const float64m8_t COEFF_B     = InaFastExp::CoeffB64();
        const float64m8_t COEFF_P5_X  = InaFastExp::GetCoefficient9_8();
        const float64m8_t COEFF_P5_Y  = InaFastExp::GetCoefficient9_7();
        const float64m8_t COEFF_P5_Z  = InaFastExp::GetCoefficient9_6();
        const float64m8_t COEFF_P5_A  = InaFastExp::GetCoefficient9_5();
        const float64m8_t COEFF_P5_B  = InaFastExp::GetCoefficient9_4();
        const float64m8_t COEFF_P5_C  = InaFastExp::GetCoefficient9_3();
        const float64m8_t COEFF_P5_D  = InaFastExp::GetCoefficient9_2();
        const float64m8_t COEFF_P5_E  = InaFastExp::GetCoefficient9_1();
        const float64m8_t COEFF_P5_F  = InaFastExp::GetCoefficient9_0();

        vfloat64m8_t x = vfmul_vf_f64m8(vec, COEFF_LOG2E);

        const vfloat64m8_t fractional_part = vfsub_vv_f64m8(x, InaVecRISCV(x).floor().vec);

        vfloat64m8_t factor = vfadd_vf_f64m8( vfmul_vv_f64m8( vfadd_vf_f64m8(
                         vfmul_vv_f64m8( vfadd_vf_f64m8( vfmul_vv_f64m8( vfadd_vf_f64m8(
                         vfmul_vv_f64m8( vfadd_vf_f64m8( vfmul_vv_f64m8( vfadd_vf_f64m8(
                         vfmul_vv_f64m8( vfadd_vf_f64m8( vfmul_vv_f64m8( vfadd_vf_f64m8(
                         vfmul_vf_f64m8(fractional_part,COEFF_P5_X), COEFF_P5_Y),fractional_part),
                         COEFF_P5_Z),fractional_part), COEFF_P5_A), fractional_part),
                         COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part),
                         COEFF_P5_F);

        x = vfsub_vv_f64m8(x,factor);

        x = vfadd_vf_f64m8(vfmul_vf_f64m8(x, COEFF_A), COEFF_B);

        vint64m8_t castedInteger = vfcvt_rtz_x_f_v_i64m8(x);

        return vfcvt_f_x_v_f64m8(castedInteger);
    }

    inline InaVecRISCV expLowAcc() const {

        const float64m8_t COEFF_LOG2E = InaFastExp::CoeffLog2E();
        const float64m8_t COEFF_A     = InaFastExp::CoeffA64();
        const float64m8_t COEFF_B     = InaFastExp::CoeffB64();
        const float64m8_t COEFF_P5_C  = InaFastExp::GetCoefficient4_3();
        const float64m8_t COEFF_P5_D  = InaFastExp::GetCoefficient4_2();
        const float64m8_t COEFF_P5_E  = InaFastExp::GetCoefficient4_1();
        const float64m8_t COEFF_P5_F  = InaFastExp::GetCoefficient4_0();

        vfloat64m8_t x = vfmul_vf_f64m8(vec, COEFF_LOG2E);

        const vfloat64m8_t fractional_part = vfsub_vv_f64m8(x, InaVecRISCV(x).floor().vec);

        vfloat64m8_t factor = vfadd_vf_f64m8(vfmul_vv_f64m8(vfadd_vf_f64m8(
                         vfmul_vv_f64m8(vfadd_vf_f64m8(vfmul_vf_f64m8(
                                         fractional_part, COEFF_P5_C),
                                         COEFF_P5_D), fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = vfsub_vv_f64m8(x,factor);

        x = vfadd_vf_f64m8(vfmul_vf_f64m8(x, COEFF_A), COEFF_B);

        vint64m8_t castedInteger = vfcvt_rtz_x_f_v_i64m8(x);

        return  vfcvt_f_x_v_f64m8(castedInteger);
    }

    inline InaVecRISCV rsqrt() const {
      // We can use vfrsqrt7_v_f64m8_t
        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<int(GetVecLength());i++)
            tabIndex[i] = 0;
        const vuint64m8_t index = vle64_v_u64m8(tabIndex);

        const vfloat64m8_t one = vlxei64_v_f64m8(1.0,index);
        return  vfsqrt_v_f64m8(vfdiv_vv_f64m8(one, vec));
    }

    inline InaVecRISCV abs() const {
        return vfabs_v_f64m8(vec);
    }

    inline InaVecRISCV floor() const {
        vfloat64m8_t valuesInIntervals = vfmin_vf_f64m8(
                                    vfmax_vf_f64m8( vec,double(std::numeric_limits<long int>::min())),
                                    double(std::numeric_limits<long int>::max()));
        vint64m8_t vecConvLongInt = vfcvt_rtz_x_f_v_i64m8(valuesInIntervals);
        vfloat64m8_t vecConvLongIntDouble = vfcvt_f_x_v_f64m8(vecConvLongInt);
        vbool8_t maskPositive = vmflt_vf_f64m8_b8(vec,0);
        return vmerge_vvm_f64m8(maskPositive, vecConvLongIntDouble, vfsub_vv_f64m8(vecConvLongIntDouble,1.0));
    }

    inline InaVecRISCV signOf() const {
        vbool8_t maskPositive = vmfgt_vf_f64m8_b8(vec,0);
        vbool8_t maskNegative = vmflt_vf_f64m8_b8(vec,0);
        const double positif = 1.0;
        const double negatif = -1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);

        vfloat64m8_t signPositive = vlxei64_v_f64m8(&positif,index);
        vfloat64m8_t signNegative = vlxei64_v_f64m8(&negatif,index);

        return vmerge_vvm_f64m8(maskNegative,signNegative,
            vfmerge_vfm_f64m8(maskPositive,signPositive,0));
    }

    inline InaVecRISCV isPositive() const {
        vbool8_t maskPositive = vmfge_vf_f64m8_b8(vec,0);
        const double positif = 1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t signPositive = vlxei64_v_f64m8(1,index);

        return vfmerge_vfm_f64m8(maskPositive,signPositive,0);
    }

    inline InaVecRISCV isNegative() const {
        vbool8_t maskNegative = vmfle_vf_f64m8_b8(vec,0);
        const double negatif = -1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t signNegative = vlxei64_v_f64m8(&negatif,index);

        return vmerge_vfm_f64m8(maskNegative,signNegative,0);
    }

    inline InaVecRISCV isPositiveStrict() const {
        vbool8_t maskPositive = vmfgt_vf_f64m8_b8(vec,0);
        const double positif = 1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t signPositive = vlxei64_v_f64m8(1,index);

        return vfmerge_vfm_f64m8(maskPositive,signPositive,0);
    }

    inline InaVecRISCV isNegativeStrict() const {
        vbool8_t maskNegative = vmflt_vf_f64m8_b8(vec,0);
        const double negatif = -1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t signNegative = vlxei64_v_f64m8(&negatif,index);

        return vmerge_vfm_f64m8(maskNegative,signNegative,0);
    }

    inline InaVecRISCV isZero() const {
        vbool8_t maskEqual = vmfeq_vf_f64m8_b8(vec,0);
        const double Equal = 1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t signEqual = vlxei64_v_f64m8(1,index);

        return vfmerge_vfm_f64m8(maskEqual,signEqual,0);
    }

    inline InaVecRISCV isNotZero() const {
        vbool8_t maskNotEqual = vmfne_vf_f64m8_b8(vec,0);
        const double NotEqual = 1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t signNotEqual = vlxei64_v_f64m8(1,index);

        return vfmerge_vfm_f64m8(maskNotEqual,signNotEqual,0);
    }

    inline InaVecMaskRISCV<double> isPositiveMask() const {
        return vmfge_vf_f64m8_b8(vec, 0);
    }

    inline InaVecMaskRISCV<double> isNegativeMask() const {
        return vmfle_vf_f64m8_b8(vec,0);
    }

    inline InaVecMaskRISCV<double> isPositiveStrictMask() const {
        return vmfgt_vf_f64m8_b8(vec,0);
    }

    inline InaVecMaskRISCV<double> isNegativeStrictMask() const {
        return vmflt_vf_f64m8_b8(vec,0);
    }

    inline InaVecMaskRISCV<double> isZeroMask() const {
        return vmfeq_vf_f64m8_b8(vec,0);
    }

    inline InaVecMaskRISCV<double> isNotZeroMask() const {
        return vmfne_vf_f64m8_b8(vec,0);
    }

    // Static basic methods
    inline static InaVecRISCV GetZero() {
        const double zero = 0.0;
        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t vecZero = vlxei64_v_f64m8(&zero,index);
        return InaVecRISCV(vecZero);
    }

    inline static InaVecRISCV GetOne() {
        const double one = 1.0;
        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t vecOne = vlxei64_v_f64m8(&one,index);
        return InaVecRISCV(vecOne);
    }

    inline static InaVecRISCV Min(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfmin_vv_f64m8( inVec1.vec, inVec2.vec);
    }

    inline static InaVecRISCV Max(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfmax_vf_f64m8( inVec1.vec, inVec2.vec);
    }

    inline static InaVecRISCV IsLowerOrEqual(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool8_t mask = vmfle_vv_f64m8_b8(inVec1.vec,inVec2.vec);
        const double zero = 0.0;
        const double one = 1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t vecZero = vlxei64_v_f64m8(&zero,index);
        vfloat64m8_t vecOne = vlxei64_v_f64m8(&one,index);

        return vmerge_vvm_f64m8(mask,vecOne,vecZero);
    }

    inline static InaVecRISCV IsLower(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool8_t mask = vmflt_vv_f64m8_b8(inVec1.vec,inVec2.vec);
        const double zero = 0.0;
        const double one = 1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t vecZero = vlxei64_v_f64m8(&zero,index);
        vfloat64m8_t vecOne = vlxei64_v_f64m8(&one,index);

        return vmerge_vvm_f64m8(mask,vecOne,vecZero);
    }

    inline static InaVecRISCV IsGreaterOrEqual(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool8_t mask = vmfge_vv_f64m8_b8(inVec1.vec,inVec2.vec);
        const double zero = 0.0;
        const double one = 1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t vecZero = vlxei64_v_f64m8(&zero,index);
        vfloat64m8_t vecOne = vlxei64_v_f64m8(&one,index);

        return vmerge_vvm_f64m8(mask,vecOne,vecZero);
    }

    inline static InaVecRISCV IsGreater(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool8_t mask = vmfgt_vv_f64m8_b8(inVec1.vec,inVec2.vec);
        const double zero = 0.0;
        const double one = 1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t vecZero = vlxei64_v_f64m8(&zero,index);
        vfloat64m8_t vecOne = vlxei64_v_f64m8(&one,index);

        return vmerge_vvm_f64m8(mask,vecOne,vecZero);
    }

    inline static InaVecRISCV IsEqual(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool8_t mask = vmfeq_vv_f64m8_b8(inVec1.vec,inVec2.vec);
        const double zero = 0.0;
        const double one = 1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t vecZero = vlxei64_v_f64m8(&zero,index);
        vfloat64m8_t vecOne = vlxei64_v_f64m8(&one,index);

        return vmerge_vvm_f64m8(mask,vecOne,vecZero);
    }

    inline static InaVecRISCV IsNotEqual(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        vbool8_t mask = vmfne_vv_f64m8_b8(inVec1.vec,inVec2.vec);
        const double zero = 0.0;
        const double one = 1.0;

        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t vecZero = vlxei64_v_f64m8(&zero,index);
        vfloat64m8_t vecOne = vlxei64_v_f64m8(&one,index);

        return vmerge_vvm_f64m8(mask,vecOne,vecZero);
    }

    inline static InaVecMaskRISCV<double> IsLowerOrEqualMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return  vmfle_vv_f64m8_b8(inVec1.vec,inVec2.vec);
    }

    inline static InaVecMaskRISCV<double> IsLowerMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vmflt_vv_f64m8_b8(inVec1.vec,inVec2.vec);
    }

    inline static InaVecMaskRISCV<double> IsGreaterOrEqualMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vmfge_vv_f64m8_b8(inVec1.vec,inVec2.vec);
    }

    inline static InaVecMaskRISCV<double> IsGreaterMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vmfgt_vv_f64m8_b8(inVec1.vec,inVec2.vec);
    }

    inline static InaVecMaskRISCV<double> IsEqualMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vmfeq_vv_f64m8_b8(inVec1.vec,inVec2.vec);
    }

    inline static InaVecMaskRISCV<double> IsNotEqualMask(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vmfne_vv_f64m8_b8(inVec1.vec,inVec2.vec);
    }

    inline static InaVecRISCV BitsAnd(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfcvt_f_x_v_f64m8(vand_vv_i64m8(vfcvt_rtz_x_f_v_i64m8(inVec1.vec),vfcvt_rtz_x_f_v_i64m8(inVec2.vec)));
    }

    inline static InaVecRISCV BitsNotAnd(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfcvt_f_x_v_f64m8(vnot_v_i64m8(vand_vv_i64m8(vfcvt_rtz_x_f_v_i64m8(inVec1.vec),vfcvt_rtz_x_f_v_i64m8(inVec2.vec))));
    }

    inline static InaVecRISCV BitsOr(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfcvt_f_x_v_f64m8(vor_vv_i64m8(vfcvt_rtz_x_f_v_i64m8(inVec1.vec),vfcvt_rtz_x_f_v_i64m8(inVec2.vec)));
    }

    inline static InaVecRISCV BitsXor(const InaVecRISCV& inVec1, const InaVecRISCV& inVec2) {
        return vfcvt_f_x_v_f64m8(vxor_vv_i64m8(vfcvt_rtz_x_f_v_i64m8(inVec1.vec),vfcvt_rtz_x_f_v_i64m8(inVec2.vec)));
    }

    inline static const char* GetName() {
        return "InaVecRISCV<double>";
    }

    inline static InaIfElse< InaVecRISCV<double> >::ThenClass If(const InaVecMaskRIS<double>& inTest) {
       return InaIfElse< InaVecRISCV<double> >::IfClass().If(inTest);
    }

    inline static InaVecRISCV IfElse(const InaVecMaskRISCV<double>& inMask, const InaVecRISCV& inIfTrue, const InaVecRISCV& inIfFalse) {
        return vfmerge_vvm_f64m8(inMask,inIfTrue.vec,inIfFalse.vec);
    }

    inline static InaVecRISCV IfTrue(const InaVecMaskRISCV<double>& inMask, const InaVecRISCV& inIfTrue) {
        return vfmerge_vfm_f64m8(inMask,inIfTrue.vec,0);
    }

    inline static InaVecRISCV IfFalse(const InaVecMaskRISCV<double>& inMask, const InaVecRISCV& inIfFalse) {
        const double zero = 0.0;
        uint64_t tabIndex [GetVecLength()];
        for (int i=0;i<GetVecLength();i++)
            tabIndex[i] = 0;
        vuint64m8_t index = vle64_v_u64m8(tabIndex);
        vfloat64m8_t vecZero = vlxei64_v_f64m8(&zero,index);
        return vfmerge_vvm_f64m8(inMask,vecZero,inIfFalse.vec);
    }

    // Inner operators
    inline InaVecRISCV<double>& operator+=(const InaVecRISCV<double>& inVec){
        vec = vfadd_vv_f64m8(vec,inVec.vec);
        return *this;
    }

    inline InaVecRISCV<double>& operator-=(const InaVecRISCV<double>& inVec){
        vec = vfsub_vv_f64m8(vec,inVec.vec);
        return *this;
    }

    inline InaVecRISCV<double>& operator/=(const InaVecRISCV<double>& inVec){
        vec = vfdiv_vv_f64m8(vec,inVec.vec);
        return *this;
    }

    inline InaVecRISCV<double>& operator*=(const InaVecRISCV<double>& inVec){
        vec = vfmul_vv_f64m8(vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<double> operator-() const {
        return vfneg_v_f64m8(vec);
    }

    inline InaVecRISCV<double> pow(std::size_t power) const{
        return InaUtils::FastPow<InaVecRISCV<double>>(*this, power);
    }

    // Multiple sum
    template <class ... Args>
    inline static void MultiHorizontalSum(double sumRes[], const InaVecRISCV<double>& inVec1,
                                          const InaVecRISCV<double>& inVec2, const InaVecRISCV<double>& inVec3,
                                          const InaVecRISCV<double>& inVec4, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1, inVec2 );
        MultiHorizontalSum(&sumRes[2], inVec3, inVec4 );

        MultiHorizontalSum(&sumRes[4], args... );
    }

    template <class ... Args>
    inline static void MultiHorizontalSum(double sumRes[], const InaVecRISCV<double>& inVec1,
                                          const InaVecRISCV<double>& inVec2, Args ...args){
        MultiHorizontalSum(&sumRes[0], inVec1);
        MultiHorizontalSum(&sumRes[1], inVec2);

        MultiHorizontalSum(&sumRes[2], args... );
    }

    inline static void MultiHorizontalSum(double sumRes[], const InaVecRISCV<double>& inVec){
        sumRes[0] += inVec.horizontalSum();
    }

    inline static void MultiHorizontalSum(double /*sumRes*/[]){
    }

    inline static InaVecSXA<double> Fma(const InaVecSXA<double>& inValAdd, const InaVecSXA<double>& inValMul1, const InaVecSXA<double>& inValMul2){
        return vfnmadd_vv_f64m8(inValAdd.vec, inValMul1.vec, inValMul2.vec);
    }
};

// Bits operators
inline InaVecRISCV<double> operator&(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return InaVecRISCV<double>::BitsAnd(inVec1, inVec2);
}

inline InaVecRISCV<double> operator|(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return InaVecRISCV<double>::BitsOr(inVec1, inVec2);
}

inline InaVecRISCV<double> operator^(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return InaVecRISCV<double>::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecRISCV<double> operator+(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return vfadd_vv_f64m8(inVec1.getVec(), inVec2.getVec());
}

inline InaVecRISCV<double> operator-(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return vfsub_vv_f64m8(inVec1.getVec(), inVec2.getVec());
}

inline InaVecRISCV<double> operator/(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return vfdiv_vv_f64m8(inVec1.getVec(), inVec2.getVec());
}

inline InaVecRISCV<double> operator*(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return vfmul_vv_f64m8(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskRISCV<double> operator<(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return InaVecRISCV<double>::IsLowerMask(inVec1,inVec2);
}

inline InaVecMaskRISCV<double> operator<=(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return InaVecRISCV<double>::IsLowerOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskRISCV<double> operator>(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return InaVecRISCV<double>::IsGreaterMask(inVec1,inVec2);
}

inline InaVecMaskRISCV<double> operator>=(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return InaVecRISCV<double>::IsGreaterOrEqualMask(inVec1,inVec2);
}

inline InaVecMaskRISCV<double> operator==(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return InaVecRISCV<double>::IsEqualMask(inVec1,inVec2);
}

inline InaVecMaskRISCV<double> operator!=(const InaVecRISCV<double>& inVec1, const InaVecRISCV<double>& inVec2){
    return InaVecRISCV<double>::IsNotEqualMask(inVec1,inVec2);
}

#endif
