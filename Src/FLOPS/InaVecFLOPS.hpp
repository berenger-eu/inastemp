///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECFLOPS_HPP
#define INAVECFLOPS_HPP

#include <type_traits>

#include "Common/InaUtils.hpp"


class InaVecFLOPSStats {
    size_t mulop;
    size_t divop;
    size_t addop;
    size_t subop;
    size_t rsqrt;
    size_t sqrt;

public:
    InaVecFLOPSStats()
        : mulop(0), divop(0), addop(0),
          subop(0), rsqrt(0), sqrt(0){
    }

    void incMulOp(const size_t inNbEvents = 1){
        mulop += inNbEvents;
    }
    void incDivOp(const size_t inNbEvents = 1){
        divop += inNbEvents;
    }
    void incAddOp(const size_t inNbEvents = 1){
        addop += inNbEvents;
    }
    void incSubOp(const size_t inNbEvents = 1){
        subop += inNbEvents;
    }
    void incRsqrt(const size_t inNbEvents = 1){
        rsqrt += inNbEvents;
    }
    void incSqrt(const size_t inNbEvents = 1){
        sqrt += inNbEvents;
    }

    size_t getMulOp() const {
        return mulop;
    }
    size_t getDivOp() const {
        return divop;
    }
    size_t getAddOp() const {
        return addop;
    }
    size_t getSubOp() const {
        return subop;
    }
    size_t getRsqrt() const {
        return rsqrt;
    }
    size_t getSqrt() const {
        return sqrt;
    }
};


/**
 * This class defines the interface that each vector
 * type must provide.
 * It is not made to be use as an abstract data type
 * by having for example a pointer on InaVecFLOPS.
 * It simply help to implement each type.
 */
template < class VecType >
class InaVecFLOPS : public VecType{
    static InaVecFLOPSStats FlopsStats;

public:
    static const InaVecFLOPSStats& GetFlopsStats(){
        return FlopsStats;
    }

    static void ResetFlopsStats(){
        FlopsStats = InaVecFLOPSStats();
    }

    using Parent               = VecType;
    using VecRawType           = typename VecType::VecRawType;
    using MaskType             = typename VecType::MaskType;
    using RealType             = typename VecType::RealType;
    static const int Alignement= VecType::Alignement;
    static const bool IsOfFixedSize = VecType::IsOfFixedSize;

    using VecType::GetVecLength;

    using VecType::VecType;

    using VecType::IsRealFma;

    inline InaVecFLOPS(){}
    inline InaVecFLOPS(const InaVecFLOPS&) = default;
    inline InaVecFLOPS& operator = (const InaVecFLOPS&) = default;

    inline InaVecFLOPS(const VecType& inVec) : Parent(inVec){};
    inline InaVecFLOPS& operator = (const VecType& inVec){
        this->Parent::operator=(inVec);
        return *this;
    }

    inline InaVecFLOPS(const std::initializer_list<RealType> lst)
        : Parent(lst.begin()){
    }

    /** Load and store */

    inline /*not explicit*/ InaVecFLOPS(const VecRawType inVec)
        : Parent(inVec){
    }

    inline VecType& operator=(const VecRawType inVec){
        this->Parent::operator =(inVec);
        return *this;
    }

    //! Convert inValue into a VecType by duplicating it and
    //! setting all the values of the array equal to inValue
    //! @code VecType[0:last-val-idx] = inValue
    inline VecType setFromRawType(const VecRawType inValue) const {
        return Parent::setFromRawType(inValue);
    }

    inline explicit operator VecRawType() const{
        return this->Parent::operator VecRawType();
    }

    inline VecRawType getVec() const{
        return Parent::getVec();
    }

//    inline /*not explicit*/ InaVecFLOPS(const RealType inVal)
//        : Parent(inVal){
//    }

//    inline VecType& operator=(const RealType inVal){
//        this->Parent::operator =(inVal);
//        return *this;
//    }

    inline void setFromScalar(const RealType inVal){
        this->Parent::setFromScalar(inVal);
    }

    // Constructor from vec
    inline explicit InaVecFLOPS(const RealType ptr[])
        : Parent(ptr){
    }

    //! Convert values from inArray into a VecType
    //! inArray might not be aligned
    //! @code idx in [0:last-val-idx] => VecType[idx] = inArray[idx]
    inline VecType& setFromArray(const RealType ptr[]){
        Parent::setFromArray(ptr);
        return *this;
    }

    //! Convert values from inArray into a VecType
    //! inArray must be aligned
    //! @code idx in [0:last-val-idx] => VecType[idx] = inArray[idx]
    inline VecType& setFromAlignedArray(const RealType ptr[]){
        Parent::setFromAlignedArray(ptr);
        return *this;
    }

    //! Convert values at position inIndirection from inArray into a VecType
    //! inArray might not be aligned
    //! @code idx in [0:last-val-idx] => VecType[idx] = inArray[inIndirection[idx]]
    inline VecType& setFromIndirectArray(const RealType values[], const int inIndirection[]) {
        Parent::setFromIndirectArray(values, inIndirection);
        return *this;
    }

    //! Convert values at position inIndirection from inArray into a VecType
    //! inArray might not be aligned
    //! @code idx in [0:last-val-idx] => VecType[idx] = inArray[inIndirection1[idx]*inLeadingDimension
    //! @code                                                                         + inIndirection2[idx]]
    inline VecType& setFromIndirect2DArray(const RealType inArray[], const int inIndirection1[],
                                 const int inLeadingDimension, const int inIndirection2[]){
        Parent::setFromIndirect2DArray(inArray, inIndirection1, inLeadingDimension, inIndirection2);
        return *this;
    }

    //! inArray
    //! @code idx in [0:last-val-idx] => outArray[idx] = inVec[idx]
    inline void storeInArray(RealType ptr[]) const {
        Parent::storeInArray(ptr);
    }

    //! inArray might must be aligned
    //! @code idx in [0:last-val-idx] => outArray[idx] = inVec[idx]
    inline void storeInAlignedArray(RealType ptr[]) const {
        Parent::storeInAlignedArray(ptr);
    }

    // Acce to individual values
    //! Return the value at position inIdx in inVec
    //! @code return inVec[inIdx]
    inline RealType at(const int index) const {
        return Parent::at(index);
    }

    // Horizontal operation
    //! Sum all the values from inVec
    //! @code return inVec[0] + ... + inVec[last-val-idx]
    inline RealType horizontalSum() const {
        FlopsStats.incAddOp(size_t(GetVecLength())-1);
        return Parent::horizontalSum();
    }

    //! Multiply all the values from inVec
    //! @code return inVec[0] * ... * inVec[last-val-idx]
    inline RealType horizontalMul() const {
        FlopsStats.incMulOp(size_t(GetVecLength())-1);
        return Parent::horizontalMul();
    }

    //! Return the smallest value in the vector
    inline RealType minInVec() const {
        return Parent::minInVec();
    }

    //! Return the greatest value in the vector
    inline RealType maxInVec() const {
        return Parent::maxInVec();
    }

    //! Apply Sqrt to all values from inVec
    //! @code idx in [0:last-val-idx] => resVec[idx] = Sqrt(inVec[idx])
    inline VecType sqrt() const {
        FlopsStats.incSqrt((1)*size_t(GetVecLength()));
        return Parent::sqrt();
    }

    //! Apply exponential to all values from inVec
    //! @code idx in [0:last-val-idx] => resVec[idx] = Exp(inVec[idx])
    inline VecType exp() const {
        FlopsStats.incAddOp((1+5)*size_t(GetVecLength()));
        FlopsStats.incMulOp((1+1+5)*size_t(GetVecLength()));
        FlopsStats.incSubOp((1+1)*size_t(GetVecLength()));
        if(std::is_same<double,typename VecType::RealType>::value){
            FlopsStats.incAddOp((3)*size_t(GetVecLength()));
            FlopsStats.incMulOp((3)*size_t(GetVecLength()));
        }
        return Parent::exp();
    }

    //! Apply exponential to all values from inVec with low accuracy
    //! @code idx in [0:last-val-idx] => resVec[idx] = ExpLowExp(inVec[idx])
    inline VecType expLowAcc() const {
        FlopsStats.incAddOp((1+2)*size_t(GetVecLength()));
        FlopsStats.incMulOp((1+1+2)*size_t(GetVecLength()));
        FlopsStats.incSubOp((1+1)*size_t(GetVecLength()));
        if(std::is_same<double,typename VecType::RealType>::value){
            FlopsStats.incAddOp((1)*size_t(GetVecLength()));
            FlopsStats.incMulOp((1)*size_t(GetVecLength()));
        }
        return Parent::expLowAcc();
    }

    //! Apply exponential to all values from inVec
    //! @code idx in [0:last-val-idx] => resVec[idx] = Exp(inVec[idx])
    inline VecType exp10() const {
        FlopsStats.incAddOp((1+5)*size_t(GetVecLength()));
        FlopsStats.incMulOp((1+1+5)*size_t(GetVecLength()));
        FlopsStats.incSubOp((1+1)*size_t(GetVecLength()));
        if(std::is_same<double,typename VecType::RealType>::value){
            FlopsStats.incAddOp((3)*size_t(GetVecLength()));
            FlopsStats.incMulOp((3)*size_t(GetVecLength()));
        }
        return Parent::exp10();
    }

    //! Apply exponential to all values from inVec with low accuracy
    //! @code idx in [0:last-val-idx] => resVec[idx] = ExpLowExp(inVec[idx])
    inline VecType exp10LowAcc() const {
        FlopsStats.incAddOp((1+2)*size_t(GetVecLength()));
        FlopsStats.incMulOp((1+1+2)*size_t(GetVecLength()));
        FlopsStats.incSubOp((1+1)*size_t(GetVecLength()));
        if(std::is_same<double,typename VecType::RealType>::value){
            FlopsStats.incAddOp((1)*size_t(GetVecLength()));
            FlopsStats.incMulOp((1)*size_t(GetVecLength()));
        }
        return Parent::exp10LowAcc();
    }

    //! Apply exponential base 2 to all values from inVec
    //! @code idx in [0:last-val-idx] => resVec[idx] = Exp(inVec[idx])
    inline VecType exp2() const {
        FlopsStats.incAddOp((1+5)*size_t(GetVecLength()));
        FlopsStats.incMulOp((1+1+5)*size_t(GetVecLength()));
        FlopsStats.incSubOp((1+1)*size_t(GetVecLength()));
        if(std::is_same<double,typename VecType::RealType>::value){
            FlopsStats.incAddOp((3)*size_t(GetVecLength()));
            FlopsStats.incMulOp((3)*size_t(GetVecLength()));
        }
        return Parent::exp2();
    }

    //! Apply exponential base 2 to all values from inVec with low accuracy
    //! @code idx in [0:last-val-idx] => resVec[idx] = ExpLowExp(inVec[idx])
    inline VecType exp2LowAcc() const {
        FlopsStats.incAddOp((1+2)*size_t(GetVecLength()));
        FlopsStats.incMulOp((1+1+2)*size_t(GetVecLength()));
        FlopsStats.incSubOp((1+1)*size_t(GetVecLength()));
        if(std::is_same<double,typename VecType::RealType>::value){
            FlopsStats.incAddOp((1)*size_t(GetVecLength()));
            FlopsStats.incMulOp((1)*size_t(GetVecLength()));
        }
        return Parent::exp2LowAcc();
    }
    
    inline VecType log() const {
        FlopsStats.incAddOp((( std::numeric_limits<RealType>::max_digits10 * 2) + 1)*size_t(GetVecLength()));
        std::size_t pownumber = 0;
        for(int i = 0; i < ( std::numeric_limits<RealType>::max_digits10 * 2); i++)
            pownumber += InaUtils::FastPowNbMul(std::size_t(i));
        FlopsStats.incMulOp((pownumber + 1)*size_t(GetVecLength()));
        FlopsStats.incSubOp((1)*size_t(GetVecLength()));
        FlopsStats.incDivOp((( std::numeric_limits<RealType>::max_digits10 * 2) + 1)*size_t(GetVecLength()));
        return Parent::log();
    }

    inline VecType log2() const {
        FlopsStats.incAddOp((( std::numeric_limits<RealType>::max_digits10 * 2) + 1)*size_t(GetVecLength()));
        std::size_t pownumber = 0;
        for(int i = 0; i < ( std::numeric_limits<RealType>::max_digits10 * 2); i++)
            pownumber += InaUtils::FastPowNbMul(std::size_t(i));
        FlopsStats.incMulOp((pownumber + 1)*size_t(GetVecLength()));
        FlopsStats.incSubOp((1)*size_t(GetVecLength()));
        FlopsStats.incDivOp((( std::numeric_limits<RealType>::max_digits10 * 2) + 1 + 1)*size_t(GetVecLength()));
        return Parent::log2();
    }

    inline VecType log10() const {
        FlopsStats.incAddOp((( std::numeric_limits<RealType>::max_digits10 * 2) + 1)*size_t(GetVecLength()));
        std::size_t pownumber = 0;
        for(int i = 0; i < ( std::numeric_limits<RealType>::max_digits10 * 2); i++)
            pownumber += InaUtils::FastPowNbMul(std::size_t(i));
        FlopsStats.incMulOp((pownumber + 1)*size_t(GetVecLength()));
        FlopsStats.incSubOp((1)*size_t(GetVecLength()));
        FlopsStats.incDivOp((( std::numeric_limits<RealType>::max_digits10 * 2) + 1 + 1)*size_t(GetVecLength()));
        return Parent::log10();
    }

    inline VecType cos() const {
        FlopsStats.incAddOp(( std::numeric_limits<RealType>::max_digits10)*size_t(GetVecLength()));
        std::size_t pownumber = 0;
        for(int i = 0; i < ( std::numeric_limits<RealType>::max_digits10 * 2); i++)
            pownumber += InaUtils::FastPowNbMul(std::size_t(i));
        FlopsStats.incMulOp((pownumber + (std::numeric_limits<RealType>::max_digits10 * 2))*size_t(GetVecLength()));
        FlopsStats.incSubOp((std::numeric_limits<RealType>::max_digits10)*size_t(GetVecLength()));
        FlopsStats.incDivOp((std::numeric_limits<RealType>::max_digits10 * 2)*size_t(GetVecLength()));
        return Parent::cos();
    }

    inline VecType sin() const {
        FlopsStats.incAddOp(( std::numeric_limits<RealType>::max_digits10)*size_t(GetVecLength()));
        std::size_t pownumber = 0;
        for(int i = 0; i < ( std::numeric_limits<RealType>::max_digits10 * 2); i++)
            pownumber += InaUtils::FastPowNbMul(std::size_t(i));
        FlopsStats.incMulOp((pownumber + (std::numeric_limits<RealType>::max_digits10 * 2))*size_t(GetVecLength()));
        FlopsStats.incSubOp((std::numeric_limits<RealType>::max_digits10)*size_t(GetVecLength()));
        FlopsStats.incDivOp((std::numeric_limits<RealType>::max_digits10 * 2)*size_t(GetVecLength()));
        return Parent::sin();
    }

    inline VecType tan() const {
        FlopsStats.incAddOp(( std::numeric_limits<RealType>::max_digits10 * 2)*size_t(GetVecLength()));
        std::size_t pownumber = 0;
        for(int i = 0; i < ( std::numeric_limits<RealType>::max_digits10 * 2 * 2); i++)
            pownumber += InaUtils::FastPowNbMul(std::size_t(i));
        FlopsStats.incMulOp(((pownumber + (std::numeric_limits<RealType>::max_digits10 * 2)) * 2)*size_t(GetVecLength()));
        FlopsStats.incSubOp((std::numeric_limits<RealType>::max_digits10 * 2)*size_t(GetVecLength()));
        FlopsStats.incDivOp((std::numeric_limits<RealType>::max_digits10 * 2 * 2 + 1)*size_t(GetVecLength()));
        return Parent::tan();
    }
    //! Apply 1/Sqrt to all values from inVec
    //! @code idx in [0:last-val-idx] => resVec[idx] = 1/Sqrt(inVec[idx])
    inline VecType rsqrt() const {
        FlopsStats.incRsqrt((1)*size_t(GetVecLength()));
        return Parent::rsqrt();
    }

    //! Return the abs value of inVec
    //! @code idx in [0:last-val-idx] => resVec[idx] = abs(inVec[idx])
    inline VecType abs() const {
        return Parent::abs();
    }

    //! Convert inVal in integer and then to real
    //! @code idx in [0:last-val-idx] => floor(inVal[idx])
    inline VecType floor() const {
        return Parent::floor();
    }

    //! Return 1 or -1 for each value in inVec
    //! @code idx in [0:last-val-idx] => resVec[idx] = (inVec[idx]<0?-1:1)
    inline VecType signOf() const {
        return Parent::signOf();
    }

    //! Return 1 if positive or zero, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]>=0?1:0
    inline VecType isPositive() const {
        return Parent::isPositive();
    }

    //! Return 1 if negative or zero, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]<=0?1:0
    inline VecType isNegative() const {
        return Parent::isNegative();
    }

    //! Return 1 if positive, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]>0?1:0
    inline VecType isPositiveStrict() const {
        return Parent::isPositiveStrict();
    }

    //! Return 1 if negative, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]<0?1:0
    inline VecType isNegativeStrict() const {
        return Parent::isNegativeStrict();
    }

    //! Return 1 if equal to, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]==0?1:0
    inline VecType isZero() const {
        return Parent::isZero();
    }

    //! Return 1 if not equal to, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]!=0?1:0
    inline VecType isNotZero() const {
        return Parent::isNotZero();
    }

    //! Return ~0 if positive or zero, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]>=0?~0:0
    inline MaskType isPositiveMask() const {
        return Parent::isPositiveMask();
    }

    //! Return ~0 if negative or zero, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]<=0?~0:0
    inline MaskType isNegativeMask() const {
        return Parent::isNegativeMask();
    }

    //! Return ~0 if positive, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]>0?~0:0
    inline MaskType isPositiveStrictMask() const {
        return Parent::isPositiveStrictMask();
    }

    //! Return ~0 if negative, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]<0?~0:0
    inline MaskType isNegativeStrictMask() const {
        return Parent::isNegativeStrictMask();
    }

    //! Return ~0 if equal to, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]==0?~0:0
    inline MaskType isZeroMask() const {
        return Parent::isZeroMask();
    }

    //! Return ~0 if not equal to, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec[idx]!=0?~0:0
    inline MaskType isNotZeroMask() const {
        return Parent::isNotZeroMask();
    }

    // Static basic methods
    //! Return a vector where all values are 1
    //! @code VecType[0:last-val-idx] = 1
    //! Might be implemtented as ScalarToSimd(1)
    inline static VecType GetZero() {
        return Parent::GetZero();
    }

    //! Return a vector where all values are 0
    //! @code VecType[0:last-val-idx] = 0
    //! Might be implemtented as ScalarToSimd(0)
    inline static VecType GetOne() {
        return Parent::GetOne();
    }

    //! Return the minimum value between inVec1 and inVec2
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]<=inVec2[idx]?inVec1[idx]:inVec2[idx]
    inline static VecType Min(const VecType& inVec1, const VecType& inVec2) {
        return Parent::Min(inVec1, inVec2);
    }

    //! Return the maximum value between inVec1 and inVec2
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]>=inVec2[idx]?inVec1[idx]:inVec2[idx]
    inline static VecType Max(const VecType& inVec1, const VecType& inVec2) {
        return Parent::Max(inVec1, inVec2);
    }

    //! Return 1 if inVec1 <= inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]<=inVec2[idx]?1:0
    inline static VecType IsLowerOrEqual(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsLowerOrEqual(inVec1, inVec2);
    }

    //! Return 1 if inVec1 < inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]<inVec2[idx]?1:0
    inline static VecType IsLower(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsLower(inVec1, inVec2);
    }

    //! Return 1 if inVec1 >= inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]>=inVec2[idx]?1:0
    inline static VecType IsGreaterOrEqual(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsGreaterOrEqual(inVec1, inVec2);
    }

    //! Return 1 if inVec1 > inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]>inVec2[idx]?1:0
    inline static VecType IsGreater(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsGreater(inVec1, inVec2);
    }

    //! Return 1 if inVec1 == inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]==inVec2[idx]?1:0
    inline static VecType IsEqual(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsEqual(inVec1, inVec2);
    }

    //! Return 1 if inVec1 == inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]!=inVec2[idx]?1:0
    inline static VecType IsNotEqual(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsNotEqual(inVec1, inVec2);
    }

    //! Return ~0 if inVec1 <= inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]<=inVec2[idx]?~0:0
    inline static MaskType IsLowerOrEqualMask(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsLowerOrEqualMask(inVec1, inVec2);
    }

    //! Return ~0 if inVec1 < inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]<inVec2[idx]?~0:0
    inline static MaskType IsLowerMask(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsLowerMask(inVec1, inVec2);
    }

    //! Return ~0 if inVec1 >= inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]>=inVec2[idx]?~0:0
    inline static MaskType IsGreaterOrEqualMask(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsGreaterOrEqualMask(inVec1, inVec2);
    }

    //! Return ~0 if inVec1 > inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]>inVec2[idx]?~0:0
    inline static MaskType IsGreaterMask(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsGreaterMask(inVec1, inVec2);
    }

    //! Return ~0 if inVec1 == inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]==inVec2[idx]?~0:0
    inline static MaskType IsEqualMask(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsEqualMask(inVec1, inVec2);
    }

    //! Return ~0 if inVec1 == inVec2, else 0 for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]!=inVec2[idx]?~0:0
    inline static MaskType IsNotEqualMask(const VecType& inVec1, const VecType& inVec2) {
        return Parent::IsNotEqualMask(inVec1, inVec2);
    }

    /** Bits operations */
    //! Return inVec1 AND inVec2, for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]&inVec2[idx]
    inline static VecType BitsAnd(const VecType& inVec1, const VecType& inVec2) {
        return Parent::BitsAnd(inVec1, inVec2);
    }

    //! Return NOT inVec1 AND inVec2, for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = (~inVec1[idx])&inVec2[idx]
    inline static VecType BitsNotAnd(const VecType& inVec1, const VecType& inVec2) {
        return Parent::BitsNotAnd(inVec1, inVec2);
    }

    //! Return inVec1 OR inVec2, for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx]|inVec2[idx]
    inline static VecType BitsOr(const VecType& inVec1, const VecType& inVec2) {
        return Parent::BitsOr(inVec1, inVec2);
    }

    //! Return inVec1 XOR inVec2, for each value
    //! @code idx in [0:last-val-idx] => resVec[idx] = inVec1[idx] ^ inVec2[idx]
    inline static VecType BitsXor(const VecType& inVec1, const VecType& inVec2) {
        return Parent::BitsXor(inVec1, inVec2);
    }

    //! Return the name identifier of the vectorizer
    inline static  const char* GetName() {
        return Parent::GetName();
    }

    //! The method cannot be because the return type is unknown
    inline static  typename InaIfElse< VecType >::ThenClass If(const MaskType& inTest) {
        return Parent::If(inTest);
    }

    //! This is a ternary if/else test
    //! It can be implemented as : return (inMask & inIfTrue) | (~inMask & inIfFalse)
    //! @code idx in [0:last-val-idx] => resVec[idx] = inMask ? inIfTrue[idx] : inIfFalse[idx]
    inline static VecType IfElse(const MaskType& inMask, const VecType& inIfTrue, const VecType& inIfFalse) {
        return Parent::IfElse(inMask, inIfTrue, inIfFalse);
    }

    //! This is a ternary if/else test with 0 as else case
    //! It can be implemented as : return (inMask & inIfTrue)
    //! @code idx in [0:last-val-idx] => resVec[idx] = inMask ? inIfTrue[idx] : 0
    inline static VecType IfTrue(const MaskType& inMask, const VecType& inIfTrue) {
        return Parent::IfTrue(inMask, inIfTrue);
    }

    //! This is a ternary if/else test with 0 as if case
    //! It can be implemented as : return (!inMask & inIfFalse)
    //! @code idx in [0:last-val-idx] => resVec[idx] = inMask ? 0 : inIfFalse[idx]
    inline static VecType IfFalse(const MaskType& inMask, const VecType& inIfFalse) {
        return Parent::IfFalse(inMask, inIfFalse);
    }

    // Inner operators
    inline VecType& operator+=(const VecType& inVec){
        FlopsStats.incAddOp(size_t(GetVecLength()));
        this->Parent::operator +=(inVec);
        return *this;
    }

    inline VecType& operator-=(const VecType& inVec){
        FlopsStats.incSubOp(size_t(GetVecLength()));
        this->Parent::operator -=(inVec);
        return *this;
    }

    inline VecType& operator/=(const VecType& inVec){
        FlopsStats.incDivOp(size_t(GetVecLength()));
        this->Parent::operator /=(inVec);
        return *this;
    }

    inline VecType& operator*=(const VecType& inVec){
        FlopsStats.incMulOp(size_t(GetVecLength()));
        this->Parent::operator *=(inVec);
        return *this;
    }


    inline VecType pow(size_t power) const{
        FlopsStats.incMulOp(InaUtils::FastPowNbMul(power)*size_t(GetVecLength()));
        return this->Parent::pow(power);
    }


    template <class VecType2>
    friend InaVecFLOPS<VecType2> operator+(const InaVecFLOPS<VecType2>& inVec1, const InaVecFLOPS<VecType2>& inVec2);

    template <class VecType2>
    friend InaVecFLOPS<VecType2> operator-(const InaVecFLOPS<VecType2>& inVec1, const InaVecFLOPS<VecType2>& inVec2);

    template <class VecType2>
    friend InaVecFLOPS<VecType2> operator/(const InaVecFLOPS<VecType2>& inVec1, const InaVecFLOPS<VecType2>& inVec2);

    template <class VecType2>
    friend InaVecFLOPS<VecType2> operator*(const InaVecFLOPS<VecType2>& inVec1, const InaVecFLOPS<VecType2>& inVec2);


    // Multiple sum
    template <class ... Args>
    inline static void MultiHorizontalSum(RealType sumRes[], Args ...args){
        const size_t nbVecs = sizeof...(args);
        FlopsStats.incAddOp((size_t(GetVecLength())-1)*nbVecs);
        Parent::MultiHorizontalSum(sumRes, args...);
    }

    inline static VecType Fma(const VecType& inValAdd, const VecType& inValMul1, const VecType& inValMul2){
        FlopsStats.incAddOp(size_t(GetVecLength()));
        FlopsStats.incMulOp(size_t(GetVecLength()));
        return Parent::Fma(inValAdd, inValMul1, inValMul2);
    }
};

// Bits operators
template <class VecType>
inline InaVecFLOPS<VecType> operator&(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    return InaVecFLOPS<VecType>::BitsAnd(inVec1, inVec2);
}

template <class VecType>
inline InaVecFLOPS<VecType> operator|(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    return InaVecFLOPS<VecType>::BitsOr(inVec1, inVec2);
}

template <class VecType>
inline InaVecFLOPS<VecType> operator^(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    return InaVecFLOPS<VecType>::BitsXor(inVec1, inVec2);
}

// Dual operators
template <class VecType>
inline InaVecFLOPS<VecType> operator+(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    InaVecFLOPS<VecType>::FlopsStats.incAddOp(size_t(InaVecFLOPS<VecType>::GetVecLength()));
    return static_cast<const VecType&>(inVec1) + static_cast<const VecType&>(inVec2);
}

template <class VecType>
inline InaVecFLOPS<VecType> operator-(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    InaVecFLOPS<VecType>::FlopsStats.incSubOp(size_t(InaVecFLOPS<VecType>::GetVecLength()));
    return static_cast<const VecType&>(inVec1) - static_cast<const VecType&>(inVec2);
}

template <class VecType>
inline InaVecFLOPS<VecType> operator/(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    InaVecFLOPS<VecType>::FlopsStats.incDivOp(size_t(InaVecFLOPS<VecType>::GetVecLength()));
    return static_cast<const VecType&>(inVec1) / static_cast<const VecType&>(inVec2);
}

template <class VecType>
inline InaVecFLOPS<VecType> operator*(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    InaVecFLOPS<VecType>::FlopsStats.incMulOp(size_t(InaVecFLOPS<VecType>::GetVecLength()));
    return static_cast<const VecType&>(inVec1) * static_cast<const VecType&>(inVec2);
}

// Tests and comparions
template <class VecType>
inline typename InaVecFLOPS<VecType>::MaskType operator<(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    return InaVecFLOPS<VecType>::IsLowerMask(inVec1,inVec2);
}

template <class VecType>
inline typename InaVecFLOPS<VecType>::MaskType operator<=(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    return InaVecFLOPS<VecType>::IsLowerOrEqualMask(inVec1,inVec2);
}

template <class VecType>
inline typename InaVecFLOPS<VecType>::MaskType operator>(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    return InaVecFLOPS<VecType>::IsGreaterMask(inVec1,inVec2);
}

template <class VecType>
inline typename InaVecFLOPS<VecType>::MaskType operator>=(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    return InaVecFLOPS<VecType>::IsGreaterOrEqualMask(inVec1,inVec2);
}

template <class VecType>
inline typename InaVecFLOPS<VecType>::MaskType operator==(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    return InaVecFLOPS<VecType>::IsEqualMask(inVec1,inVec2);
}

template <class VecType>
inline typename InaVecFLOPS<VecType>::MaskType operator!=(const InaVecFLOPS<VecType>& inVec1, const InaVecFLOPS<VecType>& inVec2){
    return InaVecFLOPS<VecType>::IsNotEqualMask(inVec1,inVec2);
}


template < class VecType >
InaVecFLOPSStats InaVecFLOPS<VecType>::FlopsStats;

#endif // INAVECFLOPS_HPP
