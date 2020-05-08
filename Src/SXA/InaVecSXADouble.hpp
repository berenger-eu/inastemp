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
class alignas(2048) InaVecMaskSXA<double> {
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

    // Double args methods
    inline static InaVecMaskSXA And(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return _vel_andm_mmm(inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSXA NotAnd(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        __vm256 one;
        one = _vel_negm_mm(_vel_xorm_mmm(one, one));
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
class alignas(2048) InaVecSXA<double> {
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

static void printVecInt(__vr vec){
    unsigned long int array[256];
    _vel_vst_vssl(vec, 8, array, 256);

    for(int idx = 0 ; idx < /*256*/10 ; ++idx){
        printf("[%d] = %lu (%lu)\n", idx, array[idx], _vel_lvsl_svs(vec, idx));
    }
}

static void printVec(__vr vec){
    double array[256];
    _vel_vst_vssl(vec, 8, array, 256);

    for(int idx = 0 ; idx < /*256*/10 ; ++idx){
        printf("[%d] = %e (%e)\n", idx, array[idx], _vel_lvsl_svs(vec, idx));
    }
}

    inline InaVecSXA& setFromIndirectArray(const double values[], const unsigned long int inIndirection[]) {
//        __vr offset = _vel_vld_vssl(8, inIndirection, 256);
//        __vr address = _vel_pvaddu_vvvl(_vel_pvbrd_vsl(reinterpret_cast<unsigned long int>(values), 256),
//                                            offset, 256);

//printf("values %p\n", values);
//printVecInt(address);
//unsigned long int vecRet = -1;
//__vr temp = _vel_vgt_vvssl(address, 0, vecRet, 256);
//printf("vecRet %p\n", vecRet);
//printVec(temp);

//        vec = _vel_vgt_vvssl(address, 0, vecRet, 256);
//printVec(vec);
//vec = temp;
//printVec(vec);
        __vr offset = _vel_vld_vssl(8, inIndirection, 256);
        __vr address = _vel_vsfa_vvssl(offset, 3, (unsigned long)values, 256);
        vec = _vel_vgt_vvssl(address, 0, 0, 256);

        return *this;
    }

    inline InaVecSXA& setFromIndirectArray(const double values[], const int inIndirection[]) {
        unsigned long int liIndirections[256];
        for(int idx = 0 ; idx < 256 ; ++idx){
            liIndirections[idx] = static_cast<unsigned long int>(inIndirection[idx]);
        }
        setFromIndirectArray(values, liIndirections);
        return *this;
    }

    inline InaVecSXA& setFromIndirect2DArray(const double inArray[], const long int inIndirection1[],
                                 const int inLeadingDimension, const long int inIndirection2[]){
        __vr offset = _vel_vaddsl_vvvl(_vel_vld_vssl(8, inIndirection2, 256),
                     _vel_vmulul_vvvl(_vel_vbrdl_vsl(inLeadingDimension, 256),
                                      _vel_vld_vssl(8, inIndirection1, 256),
                                      256),256);
        __vr address = _vel_vsfa_vvssl(offset, 3, (unsigned long)inArray, 256);
        vec = _vel_vgt_vvssl(address, 0, 0, 256);
        return *this;
    }

    inline InaVecSXA& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        long int liIndirections1[256];
        long int liIndirections2[256];
        for(int idx = 0 ; idx < 256 ; ++idx){
            liIndirections1[idx] = inIndirection1[idx];
            liIndirections2[idx] = inIndirection2[idx];
        }
        setFromIndirect2DArray(inArray, liIndirections1, inLeadingDimension, liIndirections2);
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
        return _vel_lvsd_svs(vec, index);
    }

    // Horizontal operation
    inline double horizontalSum() const {
      return _vel_lvsd_svs(_vel_vfsumd_vvl(vec, 256),0);
    }

    inline double horizontalMul() const {
        // TODO
        double sum = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            sum *= at(idx);
        }
        return sum;
    }

    inline double minInVec() const {
        // TODO
        double min = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            min = std::min(min, at(idx));
        }
        return min;
    }

    inline double maxInVec() const {
        // TODO
        double max = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            max = std::min(max, at(idx));
        }
        return max;
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
        const __vr one = _vel_vbrdd_vsl(1.0, 256);
        return  _vel_vfsqrtd_vvl(_vel_vfdivd_vvvl(one, vec, 256), 256);
    }

    inline InaVecSXA abs() const {
      return _vel_vand_vvvl( vec, _vel_pvbrd_vsl(0x7FFFFFFFFFFFFFFFUL, 256), 256);
    }

    inline InaVecSXA floor() const {
        __vr valuesInIntervals = _vel_vfmind_vvvl(
                                    _vel_vfmaxd_vvvl( vec, _vel_vbrdd_vsl(double(std::numeric_limits<long int>::min()), 256), 256),
                                    _vel_vbrdd_vsl(double(std::numeric_limits<long int>::max()), 256), 256);
        //__vm256 maskInLongInt = _vel_andm_mmm(
        //                        _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( , vec, 256), 256),
        //                        _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( vec, _vel_vbrdd_vsl(double(std::numeric_limits<long int>::max()), 256), 256), 256));
        __vr vecConvLongInt = _vel_vcvtldrz_vvl(valuesInIntervals, 256);
        __vr vecConvLongIntDouble = _vel_vcvtdl_vvl(vecConvLongInt, 256);

        __vm256 maskPositive = _vel_vfmklgt_mvl(_vel_vfcmpd_vsvl( 0, vec, 256), 256);

        return _vel_vmrg_vvvml(vecConvLongIntDouble,
                               _vel_vfsubd_vvvl( vecConvLongIntDouble, _vel_vbrdd_vsl(1, 256), 256),
                               maskPositive,
                               256);
//        return _vel_vmrg_vvvml(_vel_vmrg_vvvml(_vel_vfsubd_vvvl( vecConvLongIntDouble, _vel_vbrdd_vsl(1, 256), 256),
//                                               vecConvLongIntDouble,
//                                               maskNegative,
//                                               256),
//                               vec,
//                               maskInLongInt,
//                               256);
    }

    inline InaVecSXA signOf() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        __vm256 maskPositive = _vel_vfmkllt_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);
        __vm256 maskNegative = _vel_vfmkllt_mvl(_vel_vfcmpd_vvvl( vec, zero, 256), 256);

        return _vel_vmrg_vvvml(_vel_vmrg_vvvml(zero,
                                               _vel_vbrdd_vsl(1, 256),
                                                maskPositive, 256),
                               _vel_vbrdd_vsl(-1, 256),
                               maskNegative, 256);
    }

    inline InaVecSXA isPositive() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        __vm256 maskPositive = _vel_vfmklle_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);

        return _vel_vmrg_vvvml(zero,
                               _vel_vbrdd_vsl(1, 256),
                               maskPositive, 256);
    }

    inline InaVecSXA isNegative() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        __vm256 maskNegative = _vel_vfmklle_mvl(_vel_vfcmpd_vvvl( vec, zero, 256), 256);

        return _vel_vmrg_vvvml(zero,
                               _vel_vbrdd_vsl(1, 256),
                               maskNegative, 256);
    }

    inline InaVecSXA isPositiveStrict() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        __vm256 maskPositive = _vel_vfmkllt_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);

        return _vel_vmrg_vvvml(zero,
                               _vel_vbrdd_vsl(1, 256),
                               maskPositive, 256);
    }

    inline InaVecSXA isNegativeStrict() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        __vm256 maskNegative = _vel_vfmkllt_mvl(_vel_vfcmpd_vvvl( vec, zero, 256), 256);

        return _vel_vmrg_vvvml( zero,
                                _vel_vbrdd_vsl(1, 256),
                               maskNegative, 256);
    }

    inline InaVecSXA isZero() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        __vm256 maskEqual = _vel_vfmkleq_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);

        return _vel_vmrg_vvvml(zero,
                               _vel_vbrdd_vsl(1, 256),
                               maskEqual, 256);
    }

    inline InaVecSXA isNotZero() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        __vm256 maskEqual = _vel_vfmkleq_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);

        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(1, 256),
                               zero,
                               maskEqual, 256);
    }

    inline InaVecMaskSXA<double> isPositiveMask() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        printf("A %d\n", __LINE__);
        __vm256 mask = _vel_vfmklle_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);
        printf("B %d\n", __LINE__);
        InaVecMaskSXA<double> mask2(mask);
        printf("C %d\n", __LINE__);
        return mask;
    }

    inline InaVecMaskSXA<double> isNegativeMask() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        return _vel_vfmklge_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);
    }

    inline InaVecMaskSXA<double> isPositiveStrictMask() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        return _vel_vfmkllt_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);
    }

    inline InaVecMaskSXA<double> isNegativeStrictMask() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        return _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);
    }

    inline InaVecMaskSXA<double> isZeroMask() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        return _vel_vfmkleq_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);
    }

    inline InaVecMaskSXA<double> isNotZeroMask() const {
        __vr zero = _vel_vbrdd_vsl(0, 256);
        return _vel_vfmklne_mvl(_vel_vfcmpd_vvvl( zero, vec, 256), 256);
    }

    // Static basic methods
    inline static InaVecSXA GetZero() {
        return InaVecSXA(_vel_vbrdd_vsl(0, 256));
    }

    inline static InaVecSXA GetOne() {
        return InaVecSXA(_vel_vbrdd_vsl(1, 256));
    }

    inline static InaVecSXA Min(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmind_vvvl( inVec1.vec, inVec2.vec, 256);
    }

    inline static InaVecSXA Max(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmaxd_vvvl( inVec1.vec, inVec2.vec, 256);
    }

    inline static InaVecSXA IsLowerOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmklle_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(0, 256),
                               _vel_vbrdd_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecSXA IsLower(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmkllt_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(0, 256),
                               _vel_vbrdd_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecSXA IsGreaterOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmklge_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(0, 256),
                               _vel_vbrdd_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecSXA IsGreater(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(0, 256),
                               _vel_vbrdd_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecSXA IsEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmkleq_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(0, 256),
                               _vel_vbrdd_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecSXA IsNotEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        __vm256 mask = _vel_vfmklne_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(0, 256),
                               _vel_vbrdd_vsl(1, 256),
                               mask, 256);
    }

    inline static InaVecMaskSXA<double> IsLowerOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklle_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<double> IsLowerMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmkllt_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<double> IsGreaterOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklge_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<double> IsGreaterMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklgt_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<double> IsEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmkleq_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
    }

    inline static InaVecMaskSXA<double> IsNotEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return _vel_vfmklne_mvl(_vel_vfcmpd_vvvl( inVec1.vec, inVec2.vec, 256), 256);
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
        return _vel_vmrg_vvvml(inIfFalse.vec, inIfTrue.vec, __vm256(inMask), 256);
    }

    inline static InaVecSXA IfTrue(const InaVecMaskSXA<double>& inMask, const InaVecSXA& inIfTrue) {
        return _vel_vmrg_vvvml(_vel_vbrdd_vsl(0, 256), inIfTrue.vec, __vm256(inMask), 256);
    }

    inline static InaVecSXA IfFalse(const InaVecMaskSXA<double>& inMask, const InaVecSXA& inIfFalse) {
        return _vel_vmrg_vvvml(inIfFalse.vec, _vel_vbrdd_vsl(0, 256), __vm256(inMask), 256);
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
