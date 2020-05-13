///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef FLOPSTESTALL_HPP
#define FLOPSTESTALL_HPP

#include "InastempGlobal.h"
#include "UTester.hpp"

#include <cmath>
#include <cstring>

template < class VecType >
class FlopsTestAll : public UTester< FlopsTestAll< VecType > > {
    using Parent = UTester< FlopsTestAll< VecType > >;

    using RealType = typename VecType::RealType;
    using MaskType = typename VecType::MaskType;

    void TestBasic() {
        UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(0));
        UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(0));
        UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(0));
        UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(0));
        UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
        UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));


        VecType a = 1;
        {
            VecType res = a + a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(1));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));

            res += a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));
        }

        VecType::ResetFlopsStats();
        {
            VecType res = a * a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(1));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));

            res *= a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));
        }

        VecType::ResetFlopsStats();
        {
            VecType res = a / a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(1));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));

            res /= a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));
        }

        VecType::ResetFlopsStats();
        {
            VecType res = a - a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(1));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));

            res -= a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));
        }

        
        VecType::ResetFlopsStats();
        {
            VecType res = (a*a) + (a/a) - (a+a) * (a-a) / a;

            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));

            res = VecType(0.);      

            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(VecType::GetVecLength()) * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(VecType::GetVecLength()) * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(VecType::GetVecLength()) * size_t(0));
        }

        VecType::ResetFlopsStats();
        {
            VecType res = VecType::Fma(VecType(1),VecType(1),VecType(1));

            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(VecType::GetVecLength()));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(VecType::GetVecLength()));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(0));

            res = VecType::Fma(res, res, res);

            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp() , size_t(2 * VecType::GetVecLength()));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp() , size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp() , size_t(2 * VecType::GetVecLength()));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp() , size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt() , size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt() , size_t(0));
        }
    }

    void SetTests() {
        Parent::AddTest(&FlopsTestAll::TestBasic, "Basic test for vec type");
    }
};

#endif
