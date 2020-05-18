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
#include "Common/InaTimer.hpp"

#include <unistd.h>
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
                       const unsigned long size){
        RealType* new_array = new RealType[size];
        for (unsigned long idx = 0; idx < size ; ++idx) {
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

            InaBlas<VecType> blas{};

            RealType* ptr_test = new RealType[32];

            // NB VALUE = 32
            blas.setZero(ptr_test, 32);

            RealType* reals2 = newArray(0, 32);

            equalArrayToArray(reals2, ptr_test, 32);
            equalToArray(vec, ptr_test);
            equalToScalar(VecType::GetZero(), ptr_test[31]);
        }


        // SCALAR
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(5);
            }

            VecType vec(reals);

            InaBlas<VecType> blas{};

            RealType* ptr_test = new RealType[16];

            // NB VALUE = 16
            blas.setScalar(ptr_test, 5, 16);

            RealType* reals2 = newArray(5, 16);

            equalArrayToArray(reals2, ptr_test, 16);
            equalToArray(vec, ptr_test);
            equalToScalar(VecType(5), ptr_test[15]);
        }

        // MULT BY SCALAR
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(20);
            }

            VecType vec(reals);

            InaBlas<VecType> blas{};

            RealType* ptr_test = new RealType[45];

            // NB VALUE = 45
            blas.setScalar(ptr_test, 5, 45);

            blas.VecMultScalar(ptr_test, 4, 45);

            RealType* reals2 = newArray(20, 45);

            equalArrayToArray(reals2, ptr_test, 45);
            equalToArray(vec, ptr_test);
            equalToScalar(VecType(20), ptr_test[44]);
        }


        // ADD BY SCALAR
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(12);
            }

            VecType vec(reals);

            InaBlas<VecType> blas{};

            RealType* ptr_test = new RealType[512];

            // NB VALUE = 512
            blas.setScalar(ptr_test, 5, 512);

            blas.VecAddScalar(ptr_test, 7, 512);

            RealType* reals2 = newArray(12, 512);

            equalArrayToArray(reals2, ptr_test, 512);
            equalToArray(vec, ptr_test);
            equalToScalar(VecType(12), ptr_test[411]);
        }

        // MULT TERM TO TERM
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(50);
            }

            VecType vec(reals);

            InaBlas<VecType> blas{};

            RealType* ptr_test1 = new RealType[32];
            RealType* ptr_test2 = new RealType[32];

            // NB VALUE = 33
            blas.setScalar(ptr_test1, 5, 32);
            blas.setScalar(ptr_test2, 10, 32);

            RealType* res = blas.MultTermToTerm(ptr_test1, ptr_test2, 32);

            RealType* reals2 = newArray(50, 32);

            equalArrayToArray(reals2, res, 32);
            equalToArray(vec, res);
            equalToScalar(VecType(50), res[31]);

            equalToScalar(VecType(50), reals2[31]);

        }

        // ADD TERM TO TERM
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(17);
            }

            VecType vec(reals);

            InaBlas<VecType> blas{};

            RealType* ptr_test1 = new RealType[256];
            RealType* ptr_test2 = new RealType[256];

            // NB VALUE = 256
            blas.setScalar(ptr_test1, 5, 256);
            blas.setScalar(ptr_test2, 12, 256);

            RealType* res = blas.AddTermToTerm(ptr_test1, ptr_test2, 256);

            RealType* reals2 = newArray(17, 256);

            equalArrayToArray(reals2, res, 256);
            equalToArray(vec, res);
            equalToScalar(VecType(17), res[255]);

            equalToScalar(VecType(17), reals2[31]);
        }

        // SCALAR PRODUCT
        {
            InaBlas<VecType> blas{};

            const unsigned long sizeVect=67;

            RealType* ptr_test1 = new RealType[sizeVect];
            RealType* ptr_test2 = new RealType[sizeVect];

            blas.setScalar(ptr_test1, 5, sizeVect);
            blas.setScalar(ptr_test2, 12, sizeVect);

            RealType scalarProduct = blas.ScalarProduct(ptr_test1, ptr_test2, 67);

            equalToScalar(VecType(4020), scalarProduct);
        }

        // MULT VECT MATRIX
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(17);
            }

            VecType vec(reals);

            InaBlas<VecType> blas{};

            // FIRST
            RealType* vectCalc = new RealType[3];
            vectCalc[0] = RealType(3);
            vectCalc[1] = RealType(4);
            vectCalc[2] = RealType(3);

            RealType* matCalc = new RealType[2*3];
            matCalc[0] = RealType(1);
            matCalc[1] = RealType(3);
            matCalc[2] = RealType(5);
            matCalc[3] = RealType(5);
            matCalc[4] = RealType(5);
            matCalc[5] = RealType(2);

            RealType* resCalc = blas.ProductVecMat(vectCalc, matCalc, 2, 3, 3);

            equalToScalar(VecType(30), resCalc[0]);
            equalToScalar(VecType(41), resCalc[1]);

            // SECOND
            RealType* vectCalc2 = new RealType[3];
            vectCalc2[0] = RealType(1);
            vectCalc2[1] = RealType(-2);
            vectCalc2[2] = RealType(0.5);

            RealType* matCalc2 = new RealType[3*3];
            matCalc2[0] = RealType(1);
            matCalc2[1] = RealType(2);
            matCalc2[2] = RealType(3);
            matCalc2[3] = RealType(4);
            matCalc2[4] = RealType(5);
            matCalc2[5] = RealType(6);
            matCalc2[6] = RealType(7);
            matCalc2[7] = RealType(8);
            matCalc2[8] = RealType(9);

            RealType* resCalc2 = blas.ProductVecMat(vectCalc2, matCalc2, 3, 3, 3);

            equalToScalar(VecType(-1.5), resCalc2[0]);
            equalToScalar(VecType(-3), resCalc2[1]);
            equalToScalar(VecType(-4.5), resCalc2[2]);


            // NB COLS = NB LINES VECTOR (NB Value)
            const unsigned long sizeVect=256;
            const unsigned long nbRows=128;
            const unsigned long nbCols=256;


            RealType* vect = new RealType[sizeVect];
            RealType* mat = new RealType[nbRows*nbCols];

            blas.setScalar(vect, 5, sizeVect);
            blas.setScalar(mat, 12, nbCols*nbRows);

            RealType* res = blas.ProductVecMat(vect, mat, nbRows, nbCols, sizeVect);

        }

    }

    void SetTests() {
        Parent::AddTest(&TestBlas::TestBasic, "Basic test for blas functions");
    }

};

#endif
