///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef CORETESTBLAS_HPP
#define CORETESTBLAS_HPP

#include "InastempGlobal.h"
#include "UTester.hpp"
#include "Common/InaBlas.hpp"
#include "Common/InaUtils.hpp"

#include <cmath>
#include <cstring>
#include <cassert>

template < class VecType >
class TestBlas : public UTester< TestBlas< VecType > > {
    using Parent = UTester< TestBlas< VecType > >;

    using RealType = typename VecType::RealType;
    using MaskType = typename VecType::MaskType;

    void equalToVecType(const VecType vec,
                       const VecType inReal) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(vec.at(int(idx)), inReal.at(int(idx)));
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        RealType realsReals[VecType::GetVecLength()];
        inReal.storeInArray( realsReals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(reals[idx], realsReals[idx]);
        }

        alignas(512) RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        alignas(512) RealType realsRealsalign[VecType::GetVecLength()];
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

        alignas(512) RealType realsalign[VecType::GetVecLength()];
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

        alignas(512) RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(realsalign[idx], inReals[idx]);
        }

        alignas(512) char reals_forcena_buffer[VecType::GetVecLength()*sizeof(RealType)+1];
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

        alignas(512) RealType realsalign[VecType::GetVecLength()];
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

        alignas(512) RealType realsalign[VecType::GetVecLength()];
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

        alignas(512) RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqualLowAcc(realsalign[idx], inReals[idx]));
        }
    }

    RealType* newArray(const long int value,
                       const long int size){
        alignas(512) char buffer[VecType::GetVecLength()*sizeof(RealType)+1];
        RealType* new_array = reinterpret_cast<RealType*>(&buffer[1]);
        for (long int idx = 0; idx < size ; ++idx) {
            new_array[idx] = RealType(value);
        }

        return new_array;
    }

    void equalArrayToArray(const RealType reals1[],
                           const RealType reals2[],
                           const long int size) {
       for (long int idx = 0; idx < size ; ++idx) {
           UASSERTEEQUAL(reals1[idx], reals2[idx]);
       }
    }

    void TestBasic() {

        // ZERO
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(0);
            }

            VecType vec(reals);

            RealType* reals2 = newArray(0, 37);

            alignas(512) char buff[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test = reinterpret_cast<RealType*>(&buff[1]);


            InaBlas<VecType> blas{};

            // NB VALUE = 37
            blas.setZero(ptr_test, 37);

            equalArrayToArray(reals2, ptr_test, 37);
            equalToArray(vec, ptr_test);
            equalToScalar(VecType::GetZero(), ptr_test[36]);
        }


        // SCALAR
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(5);
            }

            VecType vec(reals);

            RealType* reals2 = newArray(5, 45);

            InaBlas<VecType> blas{};

            alignas(512) char buff[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test = reinterpret_cast<RealType*>(&buff[1]);

            // NB VALUE = 45
            blas.setScalar(ptr_test, 5, 45);

            equalArrayToArray(reals2, ptr_test, 45);
            equalToArray(vec, ptr_test);
            equalToScalar(VecType(5), ptr_test[44]);
        }

        // MULT BY SCALAR
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(15);
            }

            VecType vec(reals);

            RealType* reals2 = newArray(15, 16);

            InaBlas<VecType> blas{};

            alignas(512) char buff[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test = reinterpret_cast<RealType*>(&buff[1]);

            // NB VALUE = 16
            blas.setScalar(ptr_test, 5, 16);

            blas.VecMultScalar(ptr_test, 3, 16);

            equalArrayToArray(reals2, ptr_test, 16);
            equalToArray(vec, ptr_test);
            equalToScalar(VecType(15), ptr_test[15]);
        }


        // ADD BY SCALAR
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(12);
            }

            VecType vec(reals);

            RealType* reals2 = newArray(12, 45);

            InaBlas<VecType> blas{};

            alignas(512) char buff[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test = reinterpret_cast<RealType*>(&buff[1]);

            // NB VALUE = 45
            blas.setScalar(ptr_test, 5, 45);

            blas.VecAddScalar(ptr_test, 7, 45);

            equalArrayToArray(reals2, ptr_test, 45);
            equalToArray(vec, ptr_test);
            equalToScalar(VecType(12), ptr_test[44]);
        }

        // MULT TERM TO TERM
        /*{
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(50);
            }

            VecType vec(reals);

            InaBlas<VecType> blas{};

            alignas(512) char buff[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test1 = reinterpret_cast<RealType*>(&buff[1]);

            alignas(512) char buff2[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test2 = reinterpret_cast<RealType*>(&buff2[1]);

            // NB VALUE = 33
            blas.setScalar(ptr_test1, 5, 33);
            blas.setScalar(ptr_test2, 10, 33);

            RealType* res = blas.MultTermToTerm(ptr_test1, ptr_test2, 33);

            RealType* reals2 = newArray(50, 33);

            equalArrayToArray(reals2, res, 33);
            equalToArray(vec, res);
            equalToScalar(VecType(50), res[32]);
        }

        // ADD TERM TO TERM
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(15);
            }

            VecType vec(reals);

            InaBlas<VecType> blas{};

            alignas(512) char buff[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test1 = reinterpret_cast<RealType*>(&buff[1]);

            alignas(512) char buff2[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test2 = reinterpret_cast<RealType*>(&buff2[1]);

            // NB VALUE = 47
            blas.setScalar(ptr_test1, 5, 47);
            blas.setScalar(ptr_test2, 10, 47);

            RealType* res = blas.AddTermToTerm(ptr_test1, ptr_test2, 47);

            RealType* reals2 = newArray(15, 47);

            equalArrayToArray(reals2, res, 47);
            equalToArray(vec, res);
            equalToScalar(VecType(15), res[46]);

            equalToScalar(VecType(15), reals2[46]);

        }*/

    }

    void SetTests() {
        Parent::AddTest(&TestBlas::TestBasic, "Basic test for blas functions");
    }

};

#endif
