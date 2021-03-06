///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef CORETESTALLSTATIC_HPP
#define CORETESTALLSTATIC_HPP

#include "InastempStaticConfig.h"
#include "UTester.hpp"

#include <cmath>
#include <cstring>
#include <cassert>

#ifdef __NEC__
#define default_alignas
#else
#define default_alignas alignas(512)
#endif

inline constexpr size_t MaxTestValues(const size_t inSizeOfOneValue){
    return 2048/inSizeOfOneValue;
}

template <class VecType, int RemainingVec>
struct MultiHorizontalSumTester{
    template <class ... Params>
    static void addToTest(typename VecType::RealType res[], typename VecType::RealType good_res[],
                          const int my_idx, Params... params){
        std::array<typename VecType::RealType, MaxTestValues(sizeof(typename VecType::RealType))> values;
        for(size_t idx = 0 ; idx < size_t(VecType::GetVecLength()) ; ++idx){
            values[idx] = typename VecType::RealType(size_t(VecType::GetVecLength()*my_idx)+idx);
        }
        VecType v(&values[0]);

        MultiHorizontalSumTester<VecType, RemainingVec-1>::addToTest(res, good_res, my_idx+1, params..., v);

        good_res[my_idx] = v.horizontalSum();
    }
};

template <class VecType>
struct MultiHorizontalSumTester<VecType,0>{
    template <class ... Params>
    static void addToTest(typename VecType::RealType res[], typename VecType::RealType /*good_res*/[],
                          const int /*my_idx*/, Params... params){
        VecType::MultiHorizontalSum(res, params...);
    }
};


template < class VecType >
class TestAll : public UTester< TestAll< VecType > > {
    using Parent = UTester< TestAll< VecType > >;

    using RealType = typename VecType::RealType;
    using MaskType = typename VecType::MaskType;

    void equalToVecType(const VecType vec,
                       const VecType inReal) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()); ++idx) {
            UASSERTEEQUAL(vec.at(int(idx)), inReal.at(int(idx)));
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        RealType realsReals[VecType::GetVecLength()];
        inReal.storeInArray( realsReals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(reals[idx], realsReals[idx]);
        }

        default_alignas RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        default_alignas RealType realsRealsalign[VecType::GetVecLength()];
        inReal.storeInArray( realsRealsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(realsalign[idx], realsRealsalign[idx]);
        }
    }

    void equalToScalar(const VecType vec,
                       const RealType inReal) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(vec.at(int(idx)), inReal);
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(reals[idx], inReal);
        }

        default_alignas RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(realsalign[idx], inReal);
        }
    }

    void equalToScalarMask(const MaskType vec,
                           const MaskType inReal) {
        VecType vecval = VecType::IfTrue(vec, VecType(1));
        VecType realval = VecType::IfTrue(inReal, VecType(1));
        equalToVecType(vecval, realval);
    }

    void equalToArray(const VecType vec,
                      const RealType inReals[]) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(vec.at(int(idx)), inReals[idx]);
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(reals[idx], inReals[idx]);
        }

        default_alignas RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(realsalign[idx], inReals[idx]);
        }

        default_alignas char reals_forcena_buffer[VecType::GetVecLength()*sizeof(RealType)+1];
        RealType* reals_forcena = reinterpret_cast<RealType*>(&reals_forcena_buffer[1]);
        vec.storeInArray( reals_forcena);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(reals_forcena[idx], inReals[idx]);
        }
    }

    bool approxEqual(const float v1, const float v2) {
        return (std::abs(v1 - v2) / v2) <= 9.9999999999E-6f;
    }

    bool approxEqual(const double v1, const double v2) {
        return (std::abs(v1 - v2) / v2) <= 9.999999999999999E-12;
    }

    bool approxEqualLowAcc(const float v1, const float v2) {
        return (std::abs(v1 - v2) / v2) <= 9.9999999999E-2f;
    }

    bool approxEqualLowAcc(const double v1, const double v2) {
        return (std::abs(v1 - v2) / v2) <= 9.999999999999999E-5;
    }

    void approxEqualToScalar(const VecType vec,
                             const RealType inReal) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(vec.at(int(idx)), inReal));
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(reals[idx], inReal));
        }

        default_alignas RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(realsalign[idx], inReal));
        }
    }

    void approxEqualToArray(const VecType vec,
                            const RealType inReals[]) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(vec.at(int(idx)), inReals[idx]));
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(reals[idx], inReals[idx]));
        }

        default_alignas RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(realsalign[idx], inReals[idx]));
        }
    }

    void approxLowAccEqualToArray(const VecType vec,
                            const RealType inReals[]) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqualLowAcc(vec.at(int(idx)), inReals[idx]));
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqualLowAcc(reals[idx], inReals[idx]));
        }

        default_alignas RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqualLowAcc(realsalign[idx], inReals[idx]));
        }
    }

    void TestBasic() {
        equalToScalar(VecType(1), 1);
        equalToScalar(VecType(RealType(0)), 0);

        {
            RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = 1;
            }
            equalToScalar(VecType(reals), 1);
        }

        {
            assert(VecType::GetVecLength() <= 256);
            const RealType rv = 0;
            VecType vconstruct {{rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv}};
            VecType vcopyconstruct = {{rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv}};
            VecType vcopyop;
            vcopyop = VecType{{rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv,
                        rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv, rv}};
            vcopyop += vconstruct * vcopyconstruct; // unused
        }


        {
            default_alignas RealType reals[VecType::GetVecLength()];
#ifndef __NEC__
            default_alignas char buffer[VecType::GetVecLength()*sizeof(RealType)+8];
            RealType* realsna = reinterpret_cast<RealType*>(buffer+8);
#else
            default_alignas RealType buffer[VecType::GetVecLength()+1];
            RealType* realsna = reinterpret_cast<RealType*>(buffer+1);
#endif

            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(idx+1);
                realsna[idx] = reals[idx];
            }

            VecType vec_no_fal(reals);
            equalToArray(vec_no_fal, reals);
            equalToArray(vec_no_fal, realsna);

            VecType vec_no_fna(realsna);
            equalToArray(vec_no_fna, reals);
            equalToArray(vec_no_fna, realsna);

            VecType vec_no_fal2;
            vec_no_fal2.setFromArray(reals);
            equalToArray(vec_no_fal2, reals);
            equalToArray(vec_no_fal2, realsna);

            VecType vec_no_fna2;
            vec_no_fna2.setFromArray(realsna);
            equalToArray(vec_no_fna2, reals);
            equalToArray(vec_no_fna2, realsna);

            VecType vec_al_fal;
            vec_al_fal.setFromAlignedArray(reals);
            equalToArray(vec_al_fal, reals);
            equalToArray(vec_al_fal, realsna);

            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                   equalToScalar(vec_no_fal.at(int(idx)), RealType(idx+1));
                   equalToScalar(vec_no_fna.at(int(idx)), RealType(idx+1));
                   equalToScalar(vec_no_fal2.at(int(idx)), RealType(idx+1));
                   equalToScalar(vec_no_fna2.at(int(idx)), RealType(idx+1));
                   equalToScalar(vec_al_fal.at(int(idx)), RealType(idx+1));
            }
        }

        {
            RealType real = 1;
            int indirect[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                indirect[idx] = 0;
            }
            equalToScalar(VecType().setFromIndirectArray(&real, indirect), 1);
        }

        {
            RealType reals[VecType::GetVecLength()];
            int indirect[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                indirect[idx] = int(idx);
                reals[idx]    = RealType(idx);
            }
            equalToArray(VecType().setFromIndirectArray(reals, indirect), reals);
        }

        {
            const size_t limiteOffsetIn = std::min(32UL, sizeof(RealType)*size_t(VecType::GetVecLength()));
            for(size_t idxOffsetIn = 0 ; idxOffsetIn < limiteOffsetIn ; ++idxOffsetIn){
                unsigned char* bufferIn[sizeof(RealType)*VecType::GetVecLength()*2];
                RealType* realsIn = reinterpret_cast<RealType*>(&bufferIn[idxOffsetIn]);
                for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                    realsIn[idx]    = RealType(idx);
                }

                VecType vec(realsIn);
                equalToArray(vec, realsIn);

                vec.setFromArray(realsIn);
                equalToArray(vec, realsIn);

                for(size_t idxOffsetOut = 0 ; idxOffsetOut < sizeof(RealType)*size_t(VecType::GetVecLength()) ; ++idxOffsetOut){
                    unsigned char* bufferOut[sizeof(RealType)*VecType::GetVecLength()*2];
                    RealType* realsOut = reinterpret_cast<RealType*>(&bufferOut[idxOffsetOut]);

                    vec.storeInArray(realsOut);
                    for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                        UASSERTEEQUAL(realsOut[idx], realsIn[idx]);
                    }
                }
            }
        }

        {
            RealType reals[VecType::GetVecLength()];
            int indirect[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                indirect[idx] = int(idx);
                reals[idx]    = RealType(idx);
            }
            equalToArray(VecType().setFromIndirectArray(reals, indirect), reals);
        }

        {
            RealType real                     = 1;
            int indirect[VecType::GetVecLength()];
            for (int idx = 0; idx < (VecType::GetVecLength()) ; ++idx) {
                indirect[idx] = 0;
            }
            equalToScalar(VecType().setFromIndirect2DArray(&real, indirect, 0, indirect), 1);
        }

        {
            RealType reals[VecType::GetVecLength() * 2];
            int indirect1[VecType::GetVecLength()];
            int indirect2[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                indirect1[idx]                   = int(idx);
                indirect2[idx]                   = 1;
                reals[idx]                       = RealType(idx);
                reals[idx + size_t(VecType::GetVecLength())] = RealType(idx * 2);
            }
            equalToArray(VecType().setFromIndirect2DArray(reals, indirect1, 0, indirect1), reals);
            equalToArray(VecType().setFromIndirect2DArray(reals, indirect2, VecType::GetVecLength(), indirect1), &reals[VecType::GetVecLength()]);
        }

        {
            UASSERTEEQUAL(VecType(1).horizontalSum(), RealType(VecType::GetVecLength()));
            UASSERTEEQUAL(VecType(10).horizontalSum(), RealType(10 * VecType::GetVecLength()));

            UASSERTEEQUAL(VecType(1).horizontalMul(), RealType(1));
            UASSERTEEQUAL(VecType(10).horizontalMul(), RealType(pow(10, VecType::GetVecLength())));
        }

        {
            equalToScalar(VecType::Min(VecType(1),
                                         VecType(1)), RealType(1));
            equalToScalar(VecType::Min(VecType(RealType(0)),
                                         VecType(1)), RealType(0));
            equalToScalar(VecType::Min(VecType(1),
                                         VecType(RealType(0))), RealType(0));

            equalToScalar(VecType::Max(VecType(1),
                                         VecType(1)), RealType(1));
            equalToScalar(VecType::Max(VecType(RealType(0)),
                                         VecType(1)), RealType(1));
            equalToScalar(VecType::Max(VecType(1),
                                         VecType(RealType(0))), RealType(1));

            equalToScalar(VecType(1).abs(), RealType(1));
            equalToScalar(VecType(-1).abs(), RealType(1));
            equalToScalar(VecType(RealType(0)).abs(), RealType(0));
        }

        {
            RealType reals[VecType::GetVecLength()];
            RealType sum = 0;
            RealType mul = 1;
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(idx);
                sum += reals[idx];
                mul *= reals[idx];
            }
            UASSERTEEQUAL(VecType(reals).horizontalSum(), sum);
            UASSERTEEQUAL(VecType(reals).horizontalMul(), mul);
        }

        {
            equalToScalar(VecType::GetZero(), 0);
            equalToScalar(VecType::GetOne(), 1);
        }

        {
            RealType reals[VecType::GetVecLength()];
            RealType expres[VecType::GetVecLength()];
            RealType expreslowacc[VecType::GetVecLength()];
            RealType sqrtres[VecType::GetVecLength()];
            RealType rsqrtres[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx]    = RealType((idx%10) + 1);
                expres[idx]   = RealType(exp(reals[idx]));
                expreslowacc[idx]   = RealType(exp(reals[idx]));
                sqrtres[idx]  = RealType(sqrt(reals[idx]));
                rsqrtres[idx] = RealType(1 / sqrt(reals[idx]));
            }

            approxEqualToArray(VecType(reals).exp(), expres);
            approxLowAccEqualToArray(VecType(reals).expLowAcc(), expreslowacc);
            approxEqualToArray(VecType(reals).sqrt(), sqrtres);
            approxEqualToArray(VecType(reals).rsqrt(), rsqrtres);

            approxEqualToScalar(VecType(RealType(0)).exp(), std::exp(RealType(0)));
        }



        {
            default_alignas RealType reals[VecType::GetVecLength()];
            default_alignas RealType expres[VecType::GetVecLength()];
            default_alignas RealType expreslowacc[VecType::GetVecLength()];
            default_alignas RealType sqrtres[VecType::GetVecLength()];
            default_alignas RealType rsqrtres[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx]    = RealType((idx%10) + 1);
                expres[idx]   = RealType(exp(reals[idx]));
                expreslowacc[idx]   = RealType(exp(reals[idx]));
                sqrtres[idx]  = RealType(sqrt(reals[idx]));
                rsqrtres[idx] = RealType(1 / sqrt(reals[idx]));
            }

            approxEqualToArray(VecType().setFromAlignedArray(reals).exp(), expres);
            approxLowAccEqualToArray(VecType().setFromAlignedArray(reals).expLowAcc(), expreslowacc);
            approxEqualToArray(VecType().setFromAlignedArray(reals).sqrt(), sqrtres);
            approxEqualToArray(VecType().setFromAlignedArray(reals).rsqrt(), rsqrtres);
        }

        {
            equalToScalar(VecType(RealType(0)).signOf(), 0);
            equalToScalar(VecType(-1).signOf(), -1);
            equalToScalar(VecType(1).signOf(), 1);
            equalToScalar(VecType(-10).signOf(), -1);
            equalToScalar(VecType(10).signOf(), 1);

            equalToScalar(VecType(RealType(0)).isPositive(), 1);
            equalToScalar(VecType(-1).isPositive(), 0);
            equalToScalar(VecType(1).isPositive(), 1);
            equalToScalar(VecType(-10).isPositive(), 0);
            equalToScalar(VecType(10).isPositive(), 1);

            equalToScalar(VecType(RealType(0)).isNegative(), 1);
            equalToScalar(VecType(-1).isNegative(), 1);
            equalToScalar(VecType(1).isNegative(), 0);
            equalToScalar(VecType(-10).isNegative(), 1);
            equalToScalar(VecType(10).isNegative(), 0);

            equalToScalar(VecType(RealType(0)).isPositiveStrict(), 0);
            equalToScalar(VecType(-1).isPositiveStrict(), 0);
            equalToScalar(VecType(1).isPositiveStrict(), 1);
            equalToScalar(VecType(-10).isPositiveStrict(), 0);
            equalToScalar(VecType(10).isPositiveStrict(), 1);

            equalToScalar(VecType(RealType(0)).isNegativeStrict(), 0);
            equalToScalar(VecType(-1).isNegativeStrict(), 1);
            equalToScalar(VecType(1).isNegativeStrict(), 0);
            equalToScalar(VecType(-10).isNegativeStrict(), 1);
            equalToScalar(VecType(10).isNegativeStrict(), 0);


            equalToScalar(VecType(RealType(0)).isZero(), 1);
            equalToScalar(VecType(RealType(0)).isNotZero(), 0);

            equalToScalar(VecType(1).isZero(), 0);
            equalToScalar(VecType(1).isNotZero(), 1);
            equalToScalar(VecType(1).isZero(), 0);
            equalToScalar(VecType(1).isNotZero(), 1);
        }
        {
            equalToScalar(VecType::IsLowerOrEqual(VecType(RealType(0)),
                                                    VecType(RealType(0))),
                          1);
            equalToScalar(VecType::IsLowerOrEqual(VecType(-1),
                                                    VecType(-1)),
                          1);
            equalToScalar(VecType::IsLowerOrEqual(VecType(-1),
                                                    VecType(1)),
                          1);
            equalToScalar(VecType::IsLowerOrEqual(VecType(1),
                                                    VecType(-1)),
                          0);

            equalToScalar(VecType::IsLower(VecType(RealType(0)),
                                             VecType(RealType(0))),
                          0);
            equalToScalar(VecType::IsLower(VecType(-1),
                                             VecType(-1)),
                          0);
            equalToScalar(VecType::IsLower(VecType(-1),
                                             VecType(1)),
                          1);
            equalToScalar(VecType::IsLower(VecType(1),
                                             VecType(-1)),
                          0);

            equalToScalar(VecType::IsGreaterOrEqual(VecType(RealType(0)),
                                                      VecType(RealType(0))),
                          1);
            equalToScalar(VecType::IsGreaterOrEqual(VecType(-1),
                                                      VecType(-1)),
                          1);
            equalToScalar(VecType::IsGreaterOrEqual(VecType(-1),
                                                      VecType(1)),
                          0);
            equalToScalar(VecType::IsGreaterOrEqual(VecType(1),
                                                      VecType(-1)),
                          1);

            equalToScalar(VecType::IsGreater(VecType(RealType(0)),
                                               VecType(RealType(0))),
                          0);
            equalToScalar(VecType::IsGreater(VecType(-1),
                                               VecType(-1)),
                          0);
            equalToScalar(VecType::IsGreater(VecType(-1),
                                               VecType(1)),
                          0);
            equalToScalar(VecType::IsGreater(VecType(1),
                                               VecType(-1)),
                          1);

            equalToScalar(VecType::IsEqual(VecType(RealType(0)),
                                             VecType(RealType(0))),
                          1);
            equalToScalar(VecType::IsEqual(VecType(-1),
                                             VecType(-1)),
                          1);
            equalToScalar(VecType::IsEqual(VecType(-1),
                                             VecType(1)),
                          0);
            equalToScalar(VecType::IsEqual(VecType(1),
                                             VecType(-1)),
                          0);

            equalToScalar(VecType::IsNotEqual(VecType(RealType(0)),
                                                VecType(RealType(0))),
                          0);
            equalToScalar(VecType::IsNotEqual(VecType(-1),
                                                VecType(-1)),
                          0);
            equalToScalar(VecType::IsNotEqual(VecType(-1),
                                                VecType(1)),
                          1);
            equalToScalar(VecType::IsNotEqual(VecType(1),
                                                VecType(-1)),
                          1);
        }
        {
            const MaskType trueMask(true);
            const MaskType falseMask(false);

            equalToScalarMask(VecType(RealType(0)).isPositiveMask(), trueMask);
            equalToScalarMask(VecType(-1).isPositiveMask(), falseMask);
            equalToScalarMask(VecType(1).isPositiveMask(), trueMask);
            equalToScalarMask(VecType(-10).isPositiveMask(), falseMask);
            equalToScalarMask(VecType(10).isPositiveMask(), trueMask);

            equalToScalarMask(VecType(RealType(0)).isNegativeMask(), trueMask);
            equalToScalarMask(VecType(-1).isNegativeMask(), trueMask);
            equalToScalarMask(VecType(1).isNegativeMask(), falseMask);
            equalToScalarMask(VecType(-10).isNegativeMask(), trueMask);
            equalToScalarMask(VecType(10).isNegativeMask(), falseMask);

            equalToScalarMask(VecType(RealType(0)).isPositiveStrictMask(), falseMask);
            equalToScalarMask(VecType(-1).isPositiveStrictMask(), falseMask);
            equalToScalarMask(VecType(1).isPositiveStrictMask(), trueMask);
            equalToScalarMask(VecType(-10).isPositiveStrictMask(), falseMask);
            equalToScalarMask(VecType(10).isPositiveStrictMask(), trueMask);

            equalToScalarMask(VecType(RealType(0)).isNegativeStrictMask(), falseMask);
            equalToScalarMask(VecType(-1).isNegativeStrictMask(), trueMask);
            equalToScalarMask(VecType(1).isNegativeStrictMask(), falseMask);
            equalToScalarMask(VecType(-10).isNegativeStrictMask(), trueMask);
            equalToScalarMask(VecType(10).isNegativeStrictMask(), falseMask);


            equalToScalarMask(VecType(RealType(0)).isZeroMask(), trueMask);
            equalToScalarMask(VecType(RealType(0)).isNotZeroMask(), falseMask);

            equalToScalarMask(VecType(1).isZeroMask(), falseMask);
            equalToScalarMask(VecType(1).isNotZeroMask(), trueMask);
            equalToScalarMask(VecType(1).isZeroMask(), falseMask);
            equalToScalarMask(VecType(1).isNotZeroMask(), trueMask);
        }

        {
            const MaskType trueMask(true);
            const MaskType falseMask(false);

            equalToScalarMask(VecType::IsLowerOrEqualMask(VecType(RealType(0)),
                                                            VecType(RealType(0))),
                              trueMask);
            equalToScalarMask(VecType::IsLowerOrEqualMask(VecType(-1),
                                                            VecType(-1)),
                              trueMask);
            equalToScalarMask(VecType::IsLowerOrEqualMask(VecType(-1),
                                                            VecType(1)),
                              trueMask);
            equalToScalarMask(VecType::IsLowerOrEqualMask(VecType(1),
                                                            VecType(-1)),
                              falseMask);

            equalToScalarMask(VecType::IsLowerMask(VecType(RealType(0)),
                                                     VecType(RealType(0))),
                              falseMask);
            equalToScalarMask(VecType::IsLowerMask(VecType(-1),
                                                     VecType(-1)),
                              falseMask);
            equalToScalarMask(VecType::IsLowerMask(VecType(-1),
                                                     VecType(1)),
                              trueMask);
            equalToScalarMask(VecType::IsLowerMask(VecType(1),
                                                     VecType(-1)),
                              falseMask);

            equalToScalarMask(VecType::IsGreaterOrEqualMask(VecType(RealType(0)),
                                                              VecType(RealType(0))),
                              trueMask);
            equalToScalarMask(VecType::IsGreaterOrEqualMask(VecType(-1),
                                                              VecType(-1)),
                              trueMask);
            equalToScalarMask(VecType::IsGreaterOrEqualMask(VecType(-1),
                                                              VecType(1)),
                              falseMask);
            equalToScalarMask(VecType::IsGreaterOrEqualMask(VecType(1),
                                                              VecType(-1)),
                              trueMask);

            equalToScalarMask(VecType::IsGreaterMask(VecType(RealType(0)),
                                                       VecType(RealType(0))),
                              falseMask);
            equalToScalarMask(VecType::IsGreaterMask(VecType(-1),
                                                       VecType(-1)),
                              falseMask);
            equalToScalarMask(VecType::IsGreaterMask(VecType(-1),
                                                       VecType(1)),
                              falseMask);
            equalToScalarMask(VecType::IsGreaterMask(VecType(1),
                                                       VecType(-1)),
                              trueMask);

            equalToScalarMask(VecType::IsEqualMask(VecType(RealType(0)),
                                                     VecType(RealType(0))),
                              trueMask);
            equalToScalarMask(VecType::IsEqualMask(VecType(-1),
                                                     VecType(-1)),
                              trueMask);
            equalToScalarMask(VecType::IsEqualMask(VecType(-1),
                                                     VecType(1)),
                              falseMask);
            equalToScalarMask(VecType::IsEqualMask(VecType(1),
                                                     VecType(-1)),
                              falseMask);

            equalToScalarMask(VecType::IsNotEqualMask(VecType(RealType(0)),
                                                        VecType(RealType(0))),
                              falseMask);
            equalToScalarMask(VecType::IsNotEqualMask(VecType(-1),
                                                        VecType(-1)),
                              falseMask);
            equalToScalarMask(VecType::IsNotEqualMask(VecType(-1),
                                                        VecType(1)),
                              trueMask);
            equalToScalarMask(VecType::IsNotEqualMask(VecType(1),
                                                        VecType(-1)),
                              trueMask);
        }

        {
            equalToScalar(VecType(RealType(1)).floor(), std::floor(RealType(1)));
            equalToScalar(VecType(RealType(1.5)).floor(), std::floor(RealType(1.5)));
            equalToScalar(VecType(RealType(1.9)).floor(), std::floor(RealType(1.9)));
            equalToScalar(VecType(RealType(100000.9999)).floor(), std::floor(RealType(100000.9999)));
            equalToScalar(VecType(RealType(-1000.99)).floor(), std::floor(RealType(-1000.99)));
        }
        {
            const VecType trueMask(1);
            const VecType falseMask(RealType(0));

            equalToVecType(VecType::BitsOr(falseMask,
                                                falseMask),
                              falseMask);
            equalToVecType(VecType::BitsOr(trueMask,
                                                falseMask),
                              trueMask);
            equalToVecType(VecType::BitsOr(trueMask,
                                                trueMask),
                              trueMask);

            equalToVecType(VecType::BitsAnd(falseMask,
                                                 falseMask),
                              falseMask);
            equalToVecType(VecType::BitsAnd(trueMask,
                                                 falseMask),
                              falseMask);
            equalToVecType(VecType::BitsAnd(trueMask,
                                                 trueMask),
                              trueMask);


            equalToVecType(VecType::BitsXor(falseMask,
                                                 falseMask),
                              falseMask);
            equalToVecType(VecType::BitsXor(trueMask,
                                                 falseMask),
                              trueMask);
            equalToVecType(VecType::BitsXor(trueMask,
                                                 trueMask),
                              falseMask);


            equalToVecType(VecType::BitsNotAnd(falseMask,
                                                    falseMask),
                              falseMask);
            equalToVecType(VecType::BitsNotAnd(trueMask,
                                                    falseMask),
                              falseMask);
            equalToVecType(VecType::BitsNotAnd(trueMask,
                                                    trueMask),
                              falseMask);
            equalToVecType(VecType::BitsNotAnd(falseMask,trueMask),
                              trueMask);
        }
        {
            equalToScalar(VecType(0.) + VecType(0.),0);
            equalToScalar(VecType(0.) + VecType(10.),10);
            equalToScalar(VecType(10.) + VecType(10.),20);
            equalToScalar(VecType(10.) + VecType(0.),10);

            equalToScalar(VecType(0.) - VecType(0.),0);
            equalToScalar(VecType(0.) - VecType(10.),-10);
            equalToScalar(VecType(10.) - VecType(10.),0);
            equalToScalar(VecType(10.) - VecType(0.),10);

            equalToScalar(VecType(0.) * VecType(0.),0);
            equalToScalar(VecType(0.) * VecType(10.),0);
            equalToScalar(VecType(10.) * VecType(10.),100);
            equalToScalar(VecType(10.) * VecType(0.),0);

            equalToScalar(VecType(0.) / VecType(10.),0);
            equalToScalar(VecType(10.) / VecType(10.),1);
        }
        {
            equalToScalar(VecType(0.) += VecType(0.),0);
            equalToScalar(VecType(0.) += VecType(10.),10);
            equalToScalar(VecType(10.) += VecType(10.),20);
            equalToScalar(VecType(10.) += VecType(0.),10);

            equalToScalar(VecType(0.) -= VecType(0.),0);
            equalToScalar(VecType(0.) -= VecType(10.),-10);
            equalToScalar(VecType(10.) -= VecType(10.),0);
            equalToScalar(VecType(10.) -= VecType(0.),10);

            equalToScalar(VecType(0.) *= VecType(0.),0);
            equalToScalar(VecType(0.) *= VecType(10.),0);
            equalToScalar(VecType(10.) *= VecType(10.),100);
            equalToScalar(VecType(10.) *= VecType(0.),0);

            equalToScalar(VecType(0.) /= VecType(10.),0);
            equalToScalar(VecType(10.) /= VecType(10.),1);
        }
        {
            equalToScalar(-VecType(-4.),4.);
            equalToScalar(-VecType(4.),-4.);
            VecType a = 1.;
            equalToScalar(-a,-1.);
            equalToScalar((-1.) * a,-1.);
            a = VecType(-1.);
            equalToScalar(-a,1.);
            equalToScalar((-1.) * a,1.);
        }

        {
            equalToScalar(VecType(1.).pow(0),1.);
            equalToScalar(VecType(1.).pow(1),1.);
            equalToScalar(VecType(1.).pow(2),1.);

            equalToScalar(VecType(2.).pow(2),4.);

            equalToScalar(VecType(5.).pow(10),RealType(std::pow(5., 10)));
            equalToScalar(VecType(2.).pow(12),RealType(std::pow(2., 12)));
        }

        {
            VecType::MultiHorizontalSum(nullptr); // Should compile

            RealType res = 0;
            VecType::MultiHorizontalSum(&res, VecType(1));
            UASSERTEEQUAL(VecType(1).horizontalSum(), res);

            res = 0;
            VecType::MultiHorizontalSum(&res, VecType(10));
            UASSERTEEQUAL(VecType(10).horizontalSum(), res);
        }
        {
            RealType res[2] = {0};
            VecType v1(1);
            VecType v2(10);

            VecType::MultiHorizontalSum(res, v1, v2);

            UASSERTEEQUAL(v1.horizontalSum(), res[0]);
            UASSERTEEQUAL(v2.horizontalSum(), res[1]);
        }
        {
            RealType res[3] = {0};
            VecType v1(1);
            VecType v2(10);
            VecType v3(100);

            VecType::MultiHorizontalSum(res, v1, v2, v3);

            UASSERTEEQUAL(v1.horizontalSum(), res[0]);
            UASSERTEEQUAL(v2.horizontalSum(), res[1]);
            UASSERTEEQUAL(v3.horizontalSum(), res[2]);
        }
        {
            RealType res[4] = {0};
            VecType v1(1);
            VecType v2(10);
            VecType v3(100);
            VecType v4(1000);

            VecType::MultiHorizontalSum(res, v1, v2, v3, v4);

            UASSERTEEQUAL(v1.horizontalSum(), res[0]);
            UASSERTEEQUAL(v2.horizontalSum(), res[1]);
            UASSERTEEQUAL(v3.horizontalSum(), res[2]);
            UASSERTEEQUAL(v4.horizontalSum(), res[3]);
        }
        {
            RealType res[5] = {0};
            VecType v1(1);
            VecType v2(10);
            VecType v3(100);
            VecType v4(1000);
            VecType v5(10000);

            VecType::MultiHorizontalSum(res, v1, v2, v3, v4, v5);

            UASSERTEEQUAL(v1.horizontalSum(), res[0]);
            UASSERTEEQUAL(v2.horizontalSum(), res[1]);
            UASSERTEEQUAL(v3.horizontalSum(), res[2]);
            UASSERTEEQUAL(v4.horizontalSum(), res[3]);
            UASSERTEEQUAL(v5.horizontalSum(), res[4]);
        }
        {
            const int nb_vec_test = 3;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 4;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 5;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 6;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 7;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 8;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 9;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 10;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 11;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 12;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 13;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 14;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 15;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 16;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 17;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }
        {
            const int nb_vec_test = 18;
            RealType res[nb_vec_test] = {0};
            RealType good_res[nb_vec_test] = {0};
            MultiHorizontalSumTester<VecType,nb_vec_test>::addToTest(res, good_res, 0);
            for(int idx = 0 ; idx < nb_vec_test ; ++idx){
                UASSERTEEQUAL(good_res[idx], res[idx]);
            }
        }


        {
            RealType reals[VecType::GetVecLength()];
            for (size_t idx1 = 0; idx1 < size_t(VecType::GetVecLength()) ; ++idx1) {
                {
                    for (size_t idx2 = 0; idx2 < size_t(VecType::GetVecLength()) ; ++idx2) {
                        reals[idx2] = (idx1 == idx2 ? 10 : 0);
                    }
                    VecType vec(reals);
                    if(VecType::GetVecLength() == 1){
                        UASSERTEEQUAL(vec.minInVec(), RealType(10));
                        UASSERTEEQUAL(vec.maxInVec(), RealType(10));
                    }
                    else{
                        UASSERTEEQUAL(vec.minInVec(), RealType(0));
                        UASSERTEEQUAL(vec.maxInVec(), RealType(10));
                    }
                }
                {
                    for (size_t idx2 = 0; idx2 < size_t(VecType::GetVecLength()) ; ++idx2) {
                        reals[idx2] = (idx1 == idx2 ? -10 : 0);
                    }
                    VecType vec(reals);
                    if(VecType::GetVecLength() == 1){
                        UASSERTEEQUAL(vec.minInVec(), RealType(-10));
                        UASSERTEEQUAL(vec.maxInVec(), RealType(-10));
                    }
                    else{
                        UASSERTEEQUAL(vec.minInVec(), RealType(-10));
                        UASSERTEEQUAL(vec.maxInVec(), RealType(0));
                    }
                }
            }
        }

        {
            const VecType zero(RealType(0));
            const VecType one(RealType(1));
            equalToScalar(VecType::Fma(zero, zero, zero), RealType(0));
            equalToScalar(VecType::Fma(one, zero, zero), RealType(1));
            equalToScalar(VecType::Fma(zero, one, one), RealType(1));
            equalToScalar(VecType::Fma(one, one, one), RealType(2));

            RealType a[VecType::GetVecLength()];
            RealType b[VecType::GetVecLength()];
            RealType c[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                a[idx] = RealType(idx+1);
                b[idx] = RealType(idx*10);
                c[idx] = RealType(idx)+RealType(1.3);
            }

            RealType res[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                res[idx] = a[idx] + (b[idx] * c[idx]);
            }
            equalToArray(VecType::Fma(VecType(a), VecType(b), VecType(c)), res);
        }
    }

    void SetTests() {
        Parent::AddTest(&TestAll::TestBasic, "Basic test for vec type");
    }
};

#endif
