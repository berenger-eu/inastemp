///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INALOG_HPP
#define INALOG_HPP

#include <limits>


template < class VecType >
class InaMath {
public:
    using VecRawType           = typename VecType::VecRawType;
    using RealType             = typename VecType::RealType;

private:
    inline constexpr static double CoeffLog2() {
        // return std::ln (2.);
        return 0.693147180559945309417232121458176568075500134360;
    }

    inline constexpr static double CoeffLog10() {
        // return std::ln(10.);
        return 2.302585092994045684017991454684364207601101488629;
    }

    template <int MaxFactorial>
    inline const static std::array<RealType, MaxFactorial> ComputeFactorial(){
        std::array<RealType, MaxFactorial> factorial { 0 };
        if(MaxFactorial == 0){
            return factorial;
        }

        factorial[0] = 1;
        for(int idx = 1 ; idx < MaxFactorial ; ++idx){
            if (double(std::numeric_limits<RealType>::max()) > double(RealType(idx+1)) * double(factorial[idx-1])){
                factorial[idx] = RealType(idx+1) * factorial[idx-1];
            }
            else{
                factorial[idx] = std::numeric_limits<RealType>::max();
            }
        }
        return factorial;
    }

    constexpr static int MaxPrecisionNumber = (sizeof(RealType) == 4 ? 1 : 2);
    constexpr static int NbIterations = std::numeric_limits<RealType>::max_digits10 * MaxPrecisionNumber;
    constexpr static int MaxFactorial = (NbIterations * 2);
    const static inline std::array<RealType, MaxFactorial> Factorials { ComputeFactorial<MaxFactorial>() };

public:
    inline static VecType log(const VecType& inVec){
        VecType q = (inVec - VecType(RealType(1)))/(inVec + VecType(RealType(1)));
        VecType qpow2 = q*q;
        VecType qpowx = q;
        VecType res = VecType(RealType(0));
        for(int i = 1; i < NbIterations; i+=2){
            res += qpowx  / VecType(RealType(i)) ;
            qpowx *= qpow2;
        }
        return res * VecType(RealType(2));
    }

    inline static VecType log2(const VecType& inVec){
        return log(inVec) / VecType(RealType(CoeffLog2()));
    }

    inline static VecType log10(const VecType& inVec){
        return log(inVec) / VecType(RealType(CoeffLog10()));
    }

    inline static VecType cos(const VecType& inVec){
        VecType res = VecType(RealType(1));
        VecType vecPow2 = inVec*inVec;
        VecType vecPow2C = vecPow2;
        for (int idx=1; idx < NbIterations-1; idx += 2)
        {
            res -= vecPow2C / VecType(Factorials[ (idx*2) - 1 ]);
            vecPow2 *= vecPow2;
            res += vecPow2C / VecType(Factorials[ ((idx+1)*2) - 1 ]);
            vecPow2 *= vecPow2;
        }
        if (NbIterations%2 == 0){
            res -= vecPow2C / VecType(Factorials[ ((NbIterations-1)*2) - 1 ]);
        }
        return res;
    }

    inline static VecType sin(const VecType& inVec){
        VecType res = inVec;
        VecType vecPow2 = inVec*inVec;
        VecType vecPow2C = vecPow2 * inVec;
        for (int idx=1; idx < NbIterations-1 ; idx += 2)
        {
            res -= vecPow2C / VecType(Factorials[ (idx*2) ]);
            vecPow2 *= vecPow2;
            res += vecPow2C / VecType(Factorials[ ((idx+1)*2) ]);
            vecPow2 *= vecPow2;
        }
        if (NbIterations%2 == 0){
            res -= vecPow2C / VecType(Factorials[ ((NbIterations-1)*2) ]);
        }
        return res;
    }

    inline static VecType tan(const VecType& inVec){
        return sin(inVec)/cos(inVec);
    }

    inline static VecType asin(const VecType& inVec){
        VecType res = inVec;
        VecType coeff = VecType(RealType(0.5));
        VecType vecPow2 = inVec*inVec;
        VecType vecPow2C = vecPow2*inVec;
        for (int idx=1; idx < NbIterations; idx++)
        {
            res += (vecPow2C / VecType(RealType((idx*2) + 1))) * coeff;
            coeff *= VecType(RealType((idx*2) + 1) / RealType(idx*2));
            vecPow2C *= vecPow2;
        }
        return res;
    }

    inline static VecType acos(const VecType& inVec){
        return VecType(RealType(M_PI_2)) - asin(inVec);
    }

    inline static VecType atan(const VecType& inVec){
        VecType res = inVec;
        VecType vecPow2 = inVec*inVec;
        VecType vecPow2C = vecPow2 * inVec;
        for (int idx=1; idx < NbIterations-1 ; idx += 2)
        {
            res -= vecPow2C / VecType(RealType((idx*2) + 1));
            vecPow2C *= vecPow2;
            res += vecPow2C / VecType(RealType(((idx+1)*2) + 1));
            vecPow2C *= vecPow2;
        }
        if (NbIterations%2 == 0){
            res -= vecPow2C / VecType(RealType(((NbIterations-1)*2) + 1));
        }
        return res;
    }

    inline static VecType cosh(const VecType& inVec){
        VecType res = VecType(RealType(1));
        VecType vecPow2 = (inVec*inVec);
        VecType vecPow2C = vecPow2;
        for (int idx=1; idx < NbIterations; idx++)
        {
            res += vecPow2C / VecType(Factorials[ (idx*2) - 1 ]);
            vecPow2C *= vecPow2;
        }
        return res;
    }

    inline static VecType sinh(const VecType& inVec){
        VecType res = inVec;
        VecType vecPow2 = (inVec*inVec);
        VecType vecPow3C = vecPow2*inVec;
        for (int idx=1; idx < NbIterations; idx++)
        {
            res += vecPow3C / VecType(Factorials[ (idx*2) ]);
            vecPow3C *= vecPow2;
        }
        return res;
    }

    inline static VecType tanh(const VecType& inVec){
        return sinh(inVec)/cosh(inVec);
    }

    inline static VecType asinh(const VecType& inVec){
        return log(inVec + (VecType(RealType(1)) +  (inVec*inVec)).sqrt() );
    }

    inline static VecType acosh(const VecType& inVec){
        return log(inVec + ( (inVec + VecType(RealType(1))).sqrt() * (inVec - VecType(RealType(1))).sqrt() ) );
    }

    inline static VecType atanh(const VecType& inVec){
        return VecType(RealType(0.5)) * log( (VecType(RealType(1)) + inVec) / (VecType(RealType(1)) - inVec) );
    }
};

#endif
