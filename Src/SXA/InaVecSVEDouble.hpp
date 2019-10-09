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
    __vr mask;

public:
    // Classic constructors
    inline InaVecMaskSXA() { _ve_lvl(256);  }

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
        mask = (inBool? svptrue_b64() : svpfalse());
    }

    inline InaVecMaskSXA& operator=(const bool inBool){
        mask = (inBool? svptrue_b64() : svpfalse());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskSXA Not() const{
        return svnot_z(svptrue_b64(), mask);
    }

    inline bool isAllTrue() const{
        // Could with svptest_any(svptrue_b64(), pg)
        return svcntp_b64(svptrue_b64(), mask) == svcntd();
    }

    inline bool isAllFalse() const{
        // true if all zero
        return svcntp_b64(svptrue_b64(), mask) == 0;
    }

    // Double args methods
    inline static InaVecMaskSXA And(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return svand_z(svptrue_b64(), inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSXA NotAnd(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return svbic_z(svptrue_b64(), inMask2.mask,inMask1.mask);
    }

    inline static InaVecMaskSXA Or(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return svorr_z(svptrue_b64(), inMask1.mask,inMask2.mask);
    }

    inline static InaVecMaskSXA Xor(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return sveor_z(svptrue_b64(), inMask1.mask,inMask2.mask);
    }

    inline static bool IsEqual(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return svcntp_b64(svptrue_b64(), sveor_z(svptrue_b64(), inMask1.mask,inMask2.mask)) == 0;
    }

    inline static bool IsNotEqual(const InaVecMaskSXA& inMask1, const InaVecMaskSXA& inMask2){
        return svcntp_b64(svptrue_b64(), sveor_z(svptrue_b64(), inMask1.mask,inMask2.mask)) != 0;
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

    inline InaVecSXA() { _ve_lvl(256); }
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
        vec = (svdup_f64(val));
    }

    inline InaVecSXA& operator=(const double val){
        vec = svdup_f64(val);
        return *this;
    }

    inline void setFromScalar(const double val){
        vec = svdup_f64(val);
    }

    // Constructor from vec
    inline InaVecSXA(const std::initializer_list<double> lst)
        : InaVecSXA(lst.begin()){
    }

    inline explicit InaVecSXA(const double ptr[])
        : InaVecSXA() {
        vec = (svld1_f64(svptrue_b64(),ptr));
    }

    inline InaVecSXA& setFromArray(const double ptr[]){
        vec = svld1_f64(svptrue_b64(),ptr);
        return *this;
    }

    inline InaVecSXA& setFromAlignedArray(const double ptr[]){
        vec = svld1_f64(svptrue_b64(),ptr);
        return *this;
    }

    inline InaVecSXA& setFromIndirectArray(const double values[], const int inIndirection[]) {
        vec = svld1_gather_s64index_f64(svptrue_b64(), values, svldff1sw_s64(svptrue_b64(),inIndirection));
        return *this;
    }

    inline InaVecSXA& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){                                       
        vec = svld1_gather_s64index_f64(svptrue_b64(), inArray,
                                       svadd_s64_z(svptrue_b64(),svmul_s64_z(svptrue_b64(),svldff1sw_s64(svptrue_b64(),inIndirection1),svdup_s64(inLeadingDimension)),
                                       svldff1sw_s64(svptrue_b64(),inIndirection2)));
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
        return svlasta_f64(svwhilelt_b64(0, index),vec);
    }

    // Horizontal operation
    inline double horizontalSum() const {
      return svadda_f64(svptrue_b64(), 0, vec);
    }

    inline double horizontalMul() const {
        double sum = at(0);
        for(int idx = 1 ; idx < int(GetVecLength()) ; ++idx){
            sum *= at(idx);
        }
        return sum;
    }

    inline InaVecSXA sqrt() const {
        return svsqrt_f64_z(svptrue_b64(),vec);
    }

    inline InaVecSXA exp() const {
        const __vr COEFF_LOG2E = svdup_f64(double(InaFastExp::CoeffLog2E()));
        const __vr COEFF_A     = svdup_f64(double(InaFastExp::CoeffA64()));
        const __vr COEFF_B     = svdup_f64(double(InaFastExp::CoeffB64()));
        const __vr COEFF_P5_X  = svdup_f64(double(InaFastExp::GetCoefficient9_8()));
        const __vr COEFF_P5_Y  = svdup_f64(double(InaFastExp::GetCoefficient9_7()));
        const __vr COEFF_P5_Z  = svdup_f64(double(InaFastExp::GetCoefficient9_6()));
        const __vr COEFF_P5_A  = svdup_f64(double(InaFastExp::GetCoefficient9_5()));
        const __vr COEFF_P5_B  = svdup_f64(double(InaFastExp::GetCoefficient9_4()));
        const __vr COEFF_P5_C  = svdup_f64(double(InaFastExp::GetCoefficient9_3()));
        const __vr COEFF_P5_D  = svdup_f64(double(InaFastExp::GetCoefficient9_2()));
        const __vr COEFF_P5_E  = svdup_f64(double(InaFastExp::GetCoefficient9_1()));
        const __vr COEFF_P5_F  = svdup_f64(double(InaFastExp::GetCoefficient9_0()));

        __vr x = svmul_f64_z(svptrue_b64(),vec, COEFF_LOG2E);

        const __vr fractional_part = svsub_f64_z(svptrue_b64(),x, InaVecSXA(x).floor().vec);

        __vr factor = svadd_f64_z(svptrue_b64(),svmul_f64_z(svptrue_b64(),svadd_f64_z(svptrue_b64(),
                         svmul_f64_z(svptrue_b64(),svadd_f64_z(svptrue_b64(), svmul_f64_z(svptrue_b64(),svadd_f64_z(svptrue_b64(),
                         svmul_f64_z(svptrue_b64(),svadd_f64_z(svptrue_b64(), svmul_f64_z(svptrue_b64(),svadd_f64_z(svptrue_b64(),
                         svmul_f64_z(svptrue_b64(),svadd_f64_z(svptrue_b64(), svmul_f64_z(svptrue_b64(),svadd_f64_z(svptrue_b64(),svmul_f64_z(svptrue_b64(),
                         COEFF_P5_X, fractional_part), COEFF_P5_Y), fractional_part),
                         COEFF_P5_Z),fractional_part), COEFF_P5_A), fractional_part),
                         COEFF_P5_B), fractional_part), COEFF_P5_C),fractional_part),
                         COEFF_P5_D), fractional_part), COEFF_P5_E),fractional_part),
                         COEFF_P5_F);

        x = svsub_f64_z(svptrue_b64(),x,factor);

        x = svadd_f64_z(svptrue_b64(),svmul_f64_z(svptrue_b64(),COEFF_A, x), COEFF_B);

        svint64_t castedInteger = svcvt_s64_f64_z(svptrue_b64(),x);

        return svreinterpret_f64_s64(castedInteger);
    }

    inline InaVecSXA expLowAcc() const {
        const __vr COEFF_LOG2E = svdup_f64(double(InaFastExp::CoeffLog2E()));
        const __vr COEFF_A     = svdup_f64(double(InaFastExp::CoeffA64()));
        const __vr COEFF_B     = svdup_f64(double(InaFastExp::CoeffB64()));
        const __vr COEFF_P5_C  = svdup_f64(double(InaFastExp::GetCoefficient4_3()));
        const __vr COEFF_P5_D  = svdup_f64(double(InaFastExp::GetCoefficient4_2()));
        const __vr COEFF_P5_E  = svdup_f64(double(InaFastExp::GetCoefficient4_1()));
        const __vr COEFF_P5_F  = svdup_f64(double(InaFastExp::GetCoefficient4_0()));

        __vr x = svmul_f64_z(svptrue_b64(),vec, COEFF_LOG2E);

        const __vr fractional_part = svsub_f64_z(svptrue_b64(),x, InaVecSXA(x).floor().vec);

        __vr factor = svadd_f64_z(svptrue_b64(),svmul_f64_z(svptrue_b64(),svadd_f64_z(svptrue_b64(),
                         svmul_f64_z(svptrue_b64(),svadd_f64_z(svptrue_b64(),svmul_f64_z(svptrue_b64(),
                                         COEFF_P5_C, fractional_part),
                                         COEFF_P5_D), fractional_part),
                                         COEFF_P5_E), fractional_part),
                                         COEFF_P5_F);

        x = svsub_f64_z(svptrue_b64(),x,factor);

        x = svadd_f64_z(svptrue_b64(),svmul_f64_z(svptrue_b64(),COEFF_A, x), COEFF_B);

        svint64_t castedInteger = svcvt_s64_f64_z(svptrue_b64(),x);

        return svreinterpret_f64_s64(castedInteger);
    }

    inline InaVecSXA rsqrt() const {
        // svrsqrte_f64(vec); seems low accurate
        return  svdiv_f64_z(svptrue_b64(), svdup_f64(1), svsqrt_f64_z(svptrue_b64(),vec));
    }

    inline InaVecSXA abs() const {
      return svabs_f64_z(svptrue_b64(), vec);
    }

    inline InaVecSXA floor() const {
        __vr maskInLongInt = svand_z(svptrue_b64(),
                                svcmple_f64(svptrue_b64(), svdup_f64(double(std::numeric_limits<long int>::min())), vec),
                                svcmple_f64(svptrue_b64(), vec, svdup_f64(double(std::numeric_limits<long int>::max()))));
        svint64_t vecConvLongInt = svcvt_s64_f64_z(maskInLongInt, vec);
        __vr vecConvLongIntDouble = svcvt_f64_s64_z(maskInLongInt, vecConvLongInt);
        __vr maskNegative = svcmpgt_f64(svptrue_b64(), svdup_f64(0), vec);
        return svsel_f64(maskInLongInt, svsel_f64(maskNegative,  svsub_f64_z(svptrue_b64(), vecConvLongIntDouble, svdup_f64(1)), vecConvLongIntDouble), vec);
    }

    inline InaVecSXA signOf() const {
        return svsel_f64(svcmplt_f64(svptrue_b64(), svdup_f64(0), vec),
                  svdup_f64(1), svsel_f64(svcmpgt_f64(svptrue_b64(), svdup_f64(0), vec),
                                          svdup_f64(-1), svdup_f64(0)));
    }

    inline InaVecSXA isPositive() const {
        return svsel_f64(svcmple_f64(svptrue_b64(), svdup_f64(0), vec),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecSXA isNegative() const {
        return svsel_f64(svcmpge_f64(svptrue_b64(), svdup_f64(0), vec),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecSXA isPositiveStrict() const {
        return svsel_f64(svcmplt_f64(svptrue_b64(), svdup_f64(0), vec),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecSXA isNegativeStrict() const {
        return svsel_f64(svcmpgt_f64(svptrue_b64(), svdup_f64(0), vec),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecSXA isZero() const {
        return svsel_f64(svcmpeq_f64(svptrue_b64(), svdup_f64(0), vec),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecSXA isNotZero() const {
        return svsel_f64(svcmpne_f64(svptrue_b64(), svdup_f64(0), vec),
                         svdup_f64(1), svdup_f64(0));
    }

    inline InaVecMaskSXA<double> isPositiveMask() const {
        return svorr_z(svptrue_b64(), svcmpeq_f64(svptrue_b64(), svdup_f64(0), vec),
                       svcmple_f64(svptrue_b64(), svdup_f64(0), vec));
    }

    inline InaVecMaskSXA<double> isNegativeMask() const {
        return svorr_z(svptrue_b64(), svcmpeq_f64(svptrue_b64(), svdup_f64(0), vec),
                       svcmpge_f64(svptrue_b64(), svdup_f64(0), vec));
    }

    inline InaVecMaskSXA<double> isPositiveStrictMask() const {
        return svcmplt_f64(svptrue_b64(), svdup_f64(0), vec);
    }

    inline InaVecMaskSXA<double> isNegativeStrictMask() const {
        return svcmpgt_f64(svptrue_b64(), svdup_f64(0), vec);
    }

    inline InaVecMaskSXA<double> isZeroMask() const {
        return svcmpeq_f64(svptrue_b64(), svdup_f64(0), vec);
    }

    inline InaVecMaskSXA<double> isNotZeroMask() const {
        return svcmpne_f64(svptrue_b64(), svdup_f64(0), vec);
    }

    // Static basic methods
    inline static InaVecSXA GetZero() {
        return InaVecSXA(svdup_f64(0));
    }

    inline static InaVecSXA GetOne() {
        return InaVecSXA(svdup_f64(1));
    }

    inline static InaVecSXA Min(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svmin_f64_z(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSXA Max(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svmax_f64_z(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSXA IsLowerOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f64_z(svcmple_f64(svptrue_b64(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsLower(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f64_z(svcmplt_f64(svptrue_b64(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsGreaterOrEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f64_z(svcmpge_f64(svptrue_b64(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsGreater(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f64_z(svcmpgt_f64(svptrue_b64(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f64_z(svcmpeq_f64(svptrue_b64(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecSXA IsNotEqual(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svdup_f64_z(svcmpne_f64(svptrue_b64(), inVec1.vec, inVec2.vec), 1);
    }

    inline static InaVecMaskSXA<double> IsLowerOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmple_f64(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<double> IsLowerMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmplt_f64(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<double> IsGreaterOrEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpge_f64(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<double> IsGreaterMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpgt_f64(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<double> IsEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpeq_f64(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskSXA<double> IsNotEqualMask(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svcmpne_f64(svptrue_b64(), inVec1.vec, inVec2.vec);
    }

    inline static InaVecSXA BitsAnd(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f64_u64(svand_u64_z(svptrue_b64(), svreinterpret_u64_f64(inVec1.vec), svreinterpret_u64_f64(inVec2.vec)));
    }

    inline static InaVecSXA BitsNotAnd(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f64_u64(svbic_u64_z(svptrue_b64(), svreinterpret_u64_f64(inVec2.vec), svreinterpret_u64_f64(inVec1.vec)));
    }

    inline static InaVecSXA BitsOr(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f64_u64(svorr_u64_z(svptrue_b64(), svreinterpret_u64_f64(inVec1.vec), svreinterpret_u64_f64(inVec2.vec)));
    }

    inline static InaVecSXA BitsXor(const InaVecSXA& inVec1, const InaVecSXA& inVec2) {
        return svreinterpret_f64_u64(sveor_u64_z(svptrue_b64(), svreinterpret_u64_f64(inVec1.vec), svreinterpret_u64_f64(inVec2.vec)));
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
        return svsel_f64(__vr(inMask), inIfTrue.vec, svdup_f64(0));
    }

    inline static InaVecSXA IfFalse(const InaVecMaskSXA<double>& inMask, const InaVecSXA& inIfFalse) {
        return svsel_f64(__vr(inMask), svdup_f64(0), inIfFalse.vec);
    }

    // Inner operators
    inline InaVecSXA<double>& operator+=(const InaVecSXA<double>& inVec){
        vec = svadd_f64_z(svptrue_b64(),vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<double>& operator-=(const InaVecSXA<double>& inVec){
        vec = svsub_f64_z(svptrue_b64(),vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<double>& operator/=(const InaVecSXA<double>& inVec){
        vec = svdiv_f64_z(svptrue_b64(),vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<double>& operator*=(const InaVecSXA<double>& inVec){
        vec = svmul_f64_z(svptrue_b64(),vec,inVec.vec);
        return *this;
    }

    inline InaVecSXA<double> operator-() const {
        return svneg_f64_z(svptrue_b64(), vec);
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
    return svadd_f64_z(svptrue_b64(),inVec1.getVec(), inVec2.getVec());
}

inline InaVecSXA<double> operator-(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return svsub_f64_z(svptrue_b64(),inVec1.getVec(), inVec2.getVec());
}

inline InaVecSXA<double> operator/(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return svdiv_f64_z(svptrue_b64(),inVec1.getVec(), inVec2.getVec());
}

inline InaVecSXA<double> operator*(const InaVecSXA<double>& inVec1, const InaVecSXA<double>& inVec2){
    return svmul_f64_z(svptrue_b64(),inVec1.getVec(), inVec2.getVec());
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
