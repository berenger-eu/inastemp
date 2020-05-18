///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INALOG_HPP
#define INALOG_HPP

#include <limits>


template < class VecType >
class InaMath : public VecType{

    inline constexpr static double CoeffLog2() {
        // return std::ln (2.);
        return 0.693147180559945309417232121458176568075500134360;
    }

    inline constexpr static double CoeffLog10() {
        // return std::ln(10.);
        return 2.302585092994045684017991454684364207601101488629;
    }

public:
    using Parent               = VecType;
    using VecRawType           = typename VecType::VecRawType;
    using MaskType             = typename VecType::MaskType;
    using RealType             = typename VecType::RealType;

    using VecType::VecType;

    inline VecType log(const VecType& inVec){
        VecType one = VecType(RealType(1));
        VecType res = VecType(RealType(2));
        VecType q = (inVec - one)/(inVec + one);
        VecType restmp = VecType(RealType(0));
        for(int i = 1; i < (std::numeric_limits<RealType>::max_digits10 * 4); i+=2){
            restmp += (one / VecType(RealType(i))) * q.pow(i);
        }
        res *= restmp;
        return res;
    }

    inline VecType log2(const VecType& inVec){
        return log(inVec) / VecType(RealType(CoeffLog2()));
    }

    inline VecType log10(const VecType& inVec){
        return log(inVec) / VecType(RealType(CoeffLog10()));
    }
};

#endif
