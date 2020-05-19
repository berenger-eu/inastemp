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
    constexpr static int precisionNumber = 4;
    using Parent               = VecType;
    using VecRawType           = typename VecType::VecRawType;
    using MaskType             = typename VecType::MaskType;
    using RealType             = typename VecType::RealType;

    using VecType::VecType;

    RealType factorials[(std::numeric_limits<RealType>::max_digits10 * precisionNumber * 2)];

    RealType factorial(int n)
    {
        if (n==1)
            return 1;
        else
            return RealType(n) * factorial(n - 1);
    }

    void precalcFactorials()
    {
        for (int i=1; i<(std::numeric_limits<RealType>::max_digits10 * precisionNumber * 2)+1; i++)
        {
            factorials[i-1] = factorial(i);
        }
    }
    VecType mod(VecType a1, RealType b)
    {
        RealType r[a1.GetVecLength()];
        RealType a[a1.GetVecLength()];
        a1.storeInArray(a);
        for(int i=0; i < a1.GetVecLength(); i++){
            r[i] = std::fmod(a[i], b);
        }
        return VecType(r);
    }
    inline VecType log(const VecType& inVec){
        VecType res = VecType(RealType(2));
        VecType q = (inVec - VecType(RealType(1)))/(inVec + VecType(RealType(1)));
        VecType restmp = VecType(RealType(0));
        for(int i = 1; i < (std::numeric_limits<RealType>::max_digits10 * precisionNumber); i+=2){
            restmp += (VecType(q).pow(std::size_t(i)) / VecType(RealType(i))) ;
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

    inline VecType cos(const VecType& inVec){
        precalcFactorials();
        VecType res = VecType(RealType(1));
        VecType curTermValue;
        for (int curTerm=1; curTerm < (std::numeric_limits<RealType>::max_digits10 * precisionNumber); curTerm++)
        {
            curTermValue = inVec.pow(std::size_t(curTerm*2));
            curTermValue /= VecType(factorials[ (curTerm*2) - 1 ]);
            if (curTerm & 0x01)
                res -= curTermValue;
            else
                res += curTermValue;
        }
        return res;
    }

    inline VecType sin(const VecType& inVec){
        precalcFactorials();
        VecType res = inVec;
        VecType curTermValue;
        for (int curTerm=1; curTerm < (std::numeric_limits<RealType>::max_digits10 * precisionNumber); curTerm++)
        {
            curTermValue = inVec.pow( std::size_t((curTerm*2) + 1) );
            curTermValue /= VecType(factorials[ (curTerm*2) ]);
            if (curTerm & 0x01)
                res -= curTermValue;
            else
                res += curTermValue;
        }
        return res;
    }

    inline VecType tan(const VecType& inVec){
        return sin(inVec)/cos(inVec);
    }

    inline VecType asin(const VecType& inVec){
        VecType res = inVec;
        VecType curTermValue;
        VecType curTermValue2 = VecType(RealType(0.5));
        for (int curTerm=1; curTerm < (std::numeric_limits<RealType>::max_digits10 * precisionNumber); curTerm++)
        {
            curTermValue = inVec.pow( std::size_t((curTerm*2) + 1) );
            curTermValue /= VecType(RealType((curTerm*2) + 1));
            curTermValue *= curTermValue2;
            res += curTermValue;
            curTermValue2 *= VecType(RealType((curTerm*2) + 1) / RealType(curTerm*2));
        }
        return res;
    }

    inline VecType acos(const VecType& inVec){
        return VecType(RealType(M_PI_2)) - asin(inVec);
    }

    inline VecType atan(const VecType& inVec){
        VecType res = inVec;
        VecType curTermValue;
        for (int curTerm=1; curTerm < (std::numeric_limits<RealType>::max_digits10 * precisionNumber); curTerm++)
        {
            curTermValue = inVec.pow( std::size_t((curTerm*2) + 1) );
            curTermValue /= VecType(RealType((curTerm*2) + 1));
            if (curTerm & 0x01)
                res -= curTermValue;
            else
                res += curTermValue;
        }
        return res;
    }

    inline VecType cosh(const VecType& inVec){
        precalcFactorials();
        VecType res = VecType(RealType(1));
        VecType curTermValue;
        for (int curTerm=1; curTerm < (std::numeric_limits<RealType>::max_digits10 * precisionNumber); curTerm++)
        {
            curTermValue = inVec.pow(std::size_t(curTerm*2));
            curTermValue /= VecType(factorials[ (curTerm*2) - 1 ]);
            res += curTermValue;
        }
        return res;
    }

    inline VecType sinh(const VecType& inVec){
        precalcFactorials();
        VecType res = inVec;
        VecType curTermValue;
        for (int curTerm=1; curTerm < (std::numeric_limits<RealType>::max_digits10 * precisionNumber); curTerm++)
        {
            curTermValue = inVec.pow( std::size_t((curTerm*2) + 1) );
            curTermValue /= VecType(factorials[ (curTerm*2) ]);
            res += curTermValue;
        }
        return res;
    }

    inline VecType tanh(const VecType& inVec){
        return sinh(inVec)/cosh(inVec);
    }

    inline VecType asinh(const VecType& inVec){
        return log(inVec + (VecType(RealType(1)) +  inVec.pow( std::size_t(2) )).sqrt() );
    }

    inline VecType acosh(const VecType& inVec){
        return log(inVec + ( (inVec + VecType(RealType(1))).sqrt() * (inVec - VecType(RealType(1))).sqrt() ) );
    }

    inline VecType atanh(const VecType& inVec){
        return VecType(RealType(0.5)) * log( (VecType(RealType(1)) + inVec) / (VecType(RealType(1)) - inVec) );
    }
};

#endif
