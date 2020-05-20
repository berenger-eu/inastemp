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


    void equalApprowLowAccArrayToArray(const RealType reals1[],
                           const RealType reals2[],
                           const long int size) {
       for (long int idx = 0; idx < size ; ++idx) {
           UASSERTETRUE(approxEqualLowAcc(reals1[idx], reals2[idx]));
       }
    }


    RealType** newMatrix(const unsigned long rows,
                       const unsigned long columns){

       RealType** newMat = new RealType* [rows];

       for (unsigned long i=0; i < rows; i++)
         newMat[i] = new RealType[columns];

       return newMat;
    }

    RealType* MultMatMat2(RealType* mat1, RealType* mat2,
                                    const unsigned long nbRows1, const unsigned long nbCols,
                                    const unsigned long nbCols2){

        RealType* matRes = new RealType[nbRows1*nbCols2];

        for(unsigned long idx_row = 0 ; idx_row < nbRows1; idx_row++){
            for(unsigned long idx_col = 0 ; idx_col < nbCols2; idx_col++){
                for(unsigned long idx_row2 = 0 ; idx_row2 < nbCols; idx_row2++){
                    matRes[idx_row*nbCols2+idx_col] += mat1[idx_row*nbCols+idx_row2] * mat2[idx_row2*nbCols2+idx_col];
                }
            }
        }

        return matRes;

    }

    RealType* MatTransposee(RealType* mat,
                       const unsigned long nbRows,
                       const unsigned long nbCols){

       RealType* matRes = new RealType[nbRows*nbCols];

        for(unsigned long idxCol = 0 ; idxCol < nbCols ; idxCol++) {
            for(unsigned long idxRow = 0; idxRow < nbRows ; idxRow++){
                matRes[idxCol*nbRows+idxRow] = mat[idxRow*nbCols+idxCol];
            }
        }

        return matRes;
    }


    void equalMatrixToMatrix(RealType** mat1, RealType** mat2,
                                    const unsigned long nbRows,
                                    const unsigned long nbCols) {

        for(unsigned long idx_row = 0 ; idx_row < nbRows; idx_row++){
            for(unsigned long idx_col = 0 ; idx_col < nbCols; idx_col++){
                UASSERTEEQUAL(mat1[idx_row][idx_col], mat2[idx_row][idx_col]);
            }
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


            RealType **matCalc = new RealType* [2];
            for (int i=0; i < 2; i++)
              matCalc[i] = new RealType[3];

            matCalc[0][0] = RealType(1);
            matCalc[0][1] = RealType(3);
            matCalc[0][2] = RealType(5);
            matCalc[1][0] = RealType(5);
            matCalc[1][1] = RealType(5);
            matCalc[1][2] = RealType(2);

            RealType* resCalc = blas.ProductVecMat(vectCalc, matCalc, 2, 3);

            equalToScalar(VecType(30), resCalc[0]);
            equalToScalar(VecType(41), resCalc[1]);


        }

        // MULT MATRIX MATRIX
        {
            InaBlas<VecType> blas{};

            RealType **mat1 = new RealType* [2];
            for (int i=0; i < 2; i++)
              mat1[i] = new RealType[3];

            mat1[0][0] = RealType(1);
            mat1[0][1] = RealType(2);
            mat1[0][2] = RealType(3);
            mat1[1][0] = RealType(4);
            mat1[1][1] = RealType(5);
            mat1[1][2] = RealType(6);


            RealType **mat2 = new RealType* [3];
            for (int i=0; i < 3; i++)
              mat2[i] = new RealType[2];

            mat2[0][0] = RealType(1);
            mat2[0][1] = RealType(2);
            mat2[1][0] = RealType(3);
            mat2[1][1] = RealType(4);
            mat2[2][0] = RealType(5);
            mat2[2][1] = RealType(6);


            RealType** resCalc = blas.ProductMatMat(mat1, mat2, 2, 3, 2);


            equalToScalar(VecType(22), resCalc[0][0]);
            equalToScalar(VecType(28), resCalc[0][1]);
            equalToScalar(VecType(49), resCalc[1][0]);
            equalToScalar(VecType(64), resCalc[1][1]);
        }

        // MULT MATRIX MATRIX 2
        {
            InaBlas<VecType> blas{};

            const unsigned long nbRows = 17;
            const unsigned long nbCols = 15;
            const unsigned long size = nbRows * nbCols;
            const unsigned long nbRows2 = nbCols;
            const unsigned long nbCols2 = 21;

            RealType* mat = new RealType[size];

            for(unsigned long idxRow = 0 ; idxRow < nbRows ; idxRow++) {
                for(unsigned long idxCol = 0 ; idxCol < nbCols ; idxCol++){
                    mat[(idxRow*nbCols)+idxCol]=(idxRow*idxCol)+idxCol+idxRow+RealType(10);
                }
            }

            RealType* mat2 = new RealType[nbRows2*nbCols2];

            for(unsigned long idxRow = 0 ; idxRow < nbRows2 ; idxRow++) {
                for(unsigned long idxCol = 0 ; idxCol < nbCols2 ; idxCol++){
                    mat2[(idxRow*nbCols2)+idxCol]=(idxRow*idxCol)+idxCol+idxRow+RealType(10);
                }
            }

            RealType* multMat = blas.ProductMatMat2(mat, mat2, nbRows, nbCols, nbCols2);
            RealType* multMatNoVect = MultMatMat2(mat, mat2, nbRows, nbCols, nbCols2);

            equalApprowLowAccArrayToArray(multMat, multMatNoVect, nbRows*nbCols2);
        }

        // TRANSPOSEE
        {
            InaBlas<VecType> blas{};

            const unsigned long nbRows = 37;
            const unsigned long nbCols = 23;
            const unsigned long size = nbRows * nbCols;

            RealType* mat = new RealType[nbRows*nbCols];

            for(unsigned long idxRow = 0 ; idxRow < nbRows ; idxRow++) {
                for(unsigned long idxCol = 0 ; idxCol < nbCols ; idxCol++){
                    mat[(idxRow*nbCols)+idxCol]=(idxRow*idxCol)+idxCol+idxRow;
                }
            }

            RealType* matTrans = blas.Transposee(mat, nbRows, nbCols);
            equalToScalar(mat[0], matTrans[0]);
            equalToScalar(mat[nbRows*nbCols-1], matTrans[nbRows*nbCols-1]);

            RealType* matTransNoVect = MatTransposee(mat, nbRows, nbCols);
            equalArrayToArray(matTransNoVect, matTrans, size);

            RealType* matTrans2 = blas.TransposeeOpti(mat, nbRows, nbCols);
            equalArrayToArray(matTransNoVect, matTrans2, size);

        }

        // MULT  TRANSPOSEE MATRIX
        {
            InaBlas<VecType> blas{};

            // cond: NB ROWS1 = NBROWS2
            const unsigned long nbRows = 23;
            const unsigned long nbCols = 12;
            const unsigned long size = nbRows * nbCols;
            const unsigned long nbRows2 = nbRows;
            const unsigned long nbCols2 = 37;

            RealType* mat = new RealType[size];

            for(unsigned long idxRow = 0 ; idxRow < nbRows ; idxRow++) {
                for(unsigned long idxCol = 0 ; idxCol < nbCols ; idxCol++){
                    mat[(idxRow*nbCols)+idxCol]=(idxRow*idxCol)+idxCol+idxRow+RealType(10);
                }
            }

            RealType* mat2 = new RealType[nbRows2*nbCols2];

            for(unsigned long idxRow = 0 ; idxRow < nbRows2 ; idxRow++) {
                for(unsigned long idxCol = 0 ; idxCol < nbCols2 ; idxCol++){
                    mat2[(idxRow*nbCols2)+idxCol]=(idxRow*idxCol)+idxCol+idxRow+RealType(10);
                }
            }

            RealType* multMat = blas.ProductTransposeeMat(mat, mat2, nbCols, nbRows, nbCols2);


            RealType* multMatNoVect = MultMatMat2(MatTransposee(mat, nbRows, nbCols), mat2, nbCols, nbRows, nbCols2);

            equalApprowLowAccArrayToArray(multMat, multMatNoVect, nbCols*nbCols2);
        }

        // MULT MATRIX TRANSPOSEE
        {

            InaBlas<VecType> blas{};

            // cond: NB COLS1 = NBCOLS2
            const unsigned long nbRows = 22;
            const unsigned long nbCols = 9;
            const unsigned long size = nbRows * nbCols;
            const unsigned long nbRows2 = 17;
            const unsigned long nbCols2 = nbCols;

            RealType* mat = new RealType[size];

            for(unsigned long idxRow = 0 ; idxRow < nbRows ; idxRow++) {
                for(unsigned long idxCol = 0 ; idxCol < nbCols ; idxCol++){
                    mat[(idxRow*nbCols)+idxCol]=(idxRow*idxCol)+idxCol+idxRow+RealType(10);
                }
            }

            RealType* mat2 = new RealType[nbRows2*nbCols2];

            for(unsigned long idxRow = 0 ; idxRow < nbRows2 ; idxRow++) {
                for(unsigned long idxCol = 0 ; idxCol < nbCols2 ; idxCol++){
                    mat2[(idxRow*nbCols2)+idxCol]=(idxRow*idxCol)+idxCol+idxRow+RealType(10);
                }
            }

            RealType* multMatMatTr = blas.ProductMatTransposee(mat, mat2, nbRows, nbCols, nbRows2);
            RealType* multMatNoVect = MultMatMat2(mat, MatTransposee(mat2, nbRows2, nbCols2), nbRows, nbCols, nbRows2);

            equalApprowLowAccArrayToArray(multMatMatTr, multMatNoVect, nbRows*nbRows2);
        }


        // MULT Trans Trans
        {
            InaBlas<VecType> blas{};

            // NBROWS = NBCOLS2
            const unsigned long nbRows = 55;
            const unsigned long nbCols = 45;
            const unsigned long size = nbRows * nbCols;
            const unsigned long nbRows2 = 49;
            const unsigned long nbCols2 = nbRows;

            RealType* mat = new RealType[size];

            for(unsigned long idxRow = 0 ; idxRow < nbRows ; idxRow++) {
                for(unsigned long idxCol = 0 ; idxCol < nbCols ; idxCol++){
                    mat[(idxRow*nbCols)+idxCol]=(idxRow*idxCol)+idxCol+idxRow+RealType(10);
                }
            }

            RealType* mat2 = new RealType[nbRows2*nbCols2];

            for(unsigned long idxRow = 0 ; idxRow < nbRows2 ; idxRow++) {
                for(unsigned long idxCol = 0 ; idxCol < nbCols2 ; idxCol++){
                    mat2[(idxRow*nbCols2)+idxCol]=(idxRow*idxCol)+idxCol+idxRow+RealType(10);
                }
            }

            RealType* multMat = blas.ProductTranspoTranspo(mat, mat2, nbCols, nbRows, nbRows2);
            RealType* multMatNoVect = MultMatMat2(MatTransposee(mat, nbRows, nbCols),
                                                  MatTransposee(mat2, nbRows2, nbCols2), nbCols, nbRows, nbRows2);

            equalApprowLowAccArrayToArray(multMat, multMatNoVect, nbCols*nbRows2);
        }

    }

    void SetTests() {
        Parent::AddTest(&TestBlas::TestBasic, "Basic test for blas functions");
    }

};

#endif
