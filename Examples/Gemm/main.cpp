///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

// In this example we time the duration to compute a gemm
// to simplify the example, the matrices are square
// and the size is a factor of the blocking coefficient

#include "InastempConfig.h"
#include "Common/InaTimer.hpp"

#include <memory>
#include <iostream>
#include <cstring>
#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
///
/// Utils
///
///////////////////////////////////////////////////////////////////////////////////

template <class RealType, size_t BlockSize>
void CopyMat(RealType* __restrict__ dest, const size_t leadingDest,
             const RealType* __restrict__ src, const size_t leadingSrc){

    for(size_t idxl = 0 ; idxl < leadingDest ; ++idxl){
        for(size_t idxb = 0 ; idxb < BlockSize ; ++idxb){
            dest[idxb*leadingDest + idxl] = src[idxb*leadingSrc+idxl];
        }
    }
}

template <class RealType, size_t BlockSize>
void CopyMatT(RealType* __restrict__ dest, const size_t leadingDest, const size_t lengthCopy,
             const RealType* __restrict__ src, const size_t leadingSrc){

    for(size_t idxl = 0 ; idxl < lengthCopy ; ++idxl){
        for(size_t idxb = 0 ; idxb < BlockSize ; ++idxb){
            dest[idxl*leadingDest + idxb] = src[idxb*leadingSrc+idxl];
        }
    }
}

template <class RealType>
void CheckEquality(const RealType goodMat[], const RealType testMat[], const size_t matDim){
    for(size_t idxRow = 0 ; idxRow < matDim ; ++idxRow){
        for(size_t idxCol = 0 ; idxCol < matDim ; ++idxCol){
            if( std::abs((testMat[idxCol * matDim + idxRow]-goodMat[idxCol * matDim + idxRow])
                         /goodMat[idxCol * matDim + idxRow])
                    > (sizeof(RealType) == 4 ? 1E-6 : 1E-10) ){
                std::cout << "Error at pos " << idxRow << " " << idxCol << "\n";
                std::cout << "    Should be " << goodMat[idxCol * matDim + idxRow] << "\n";
                std::cout << "    but is " << testMat[idxCol * matDim + idxRow] << "\n";
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////
///
/// Scalar basic implementation (very slow)
///
///////////////////////////////////////////////////////////////////////////////////


template <class RealType>
void ScalarGemmNoBlock(const RealType* __restrict__ A, const RealType* __restrict__ B,
                RealType* __restrict__ C, const size_t size){

    for(size_t j = 0 ; j < size ; ++j){
        for(size_t i = 0 ; i < size ; ++i){
            for(size_t k = 0 ; k < size ; ++k){
                C[k*size + j] += A[k*size + i] * B[j*size + i];
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////
///
/// Scalar with a blocking scheme
/// (similar as http://www.netlib.org/utk/papers/autoblock/node2.html)
///
///////////////////////////////////////////////////////////////////////////////////

template <class RealType, size_t BlockSize>
void ScalarGemm(const RealType* __restrict__ A, const RealType* __restrict__ B,
                RealType* __restrict__ C, const size_t size){

    for(size_t j_outer = 0 ; j_outer < size ; j_outer += BlockSize){
        for(size_t i_outer = 0 ; i_outer < size ; i_outer += BlockSize){
            for(size_t k = 0 ; k < size ; ++k){

                const RealType* __restrict__ ptrA = &A[k*size + i_outer];
                const RealType* __restrict__ ptrB = &B[j_outer*size + i_outer];
                RealType* __restrict__ ptrC = &C[k*size + j_outer];

                for(size_t idxBlockCol = 0 ; idxBlockCol < BlockSize ; ++idxBlockCol){
                    RealType v = 0;
                    for(size_t idxVal = 0 ; idxVal < BlockSize ; ++idxVal){
                        v += ptrA[idxVal] * ptrB[idxBlockCol*size + idxVal];
                    }
                    ptrC[idxBlockCol] += v;
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////
///
/// Scalar with a blocking scheme more close Goto-blas
///
///////////////////////////////////////////////////////////////////////////////////

template <class RealType, size_t PanelSizeK, size_t PanelSizeiA,
          size_t PanelSizejB, size_t BlockSize>
void ScalarGemmV2(const RealType* __restrict__ A, const RealType* __restrict__ B,
                RealType* __restrict__ C, const size_t size){

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK/BlockSize)*BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA/BlockSize)*BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB/BlockSize)*BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size/PanelSizeK)*PanelSizeK == size);
    assert((size/PanelSizeiA)*PanelSizeiA == size);
    assert((size/PanelSizejB)*PanelSizejB == size);

    for(size_t ip = 0 ; ip < size ; ip += PanelSizeiA){
        for(size_t jp = 0 ; jp < size ; jp += PanelSizejB){

            for(size_t kp = 0 ; kp < size ; kp += PanelSizeK){

                alignas(64) RealType panelA[PanelSizeiA*PanelSizeK];
                alignas(64) RealType panelB[PanelSizeK*BlockSize];

                for(size_t jb = 0 ; jb < PanelSizejB ; jb += BlockSize){

                    CopyMat<RealType, BlockSize>(panelB, PanelSizeK, &B[jp*size + kp], size);

                    for(size_t ib = 0 ; ib < PanelSizeiA ; ib += BlockSize){

                        if(jb == 0){
                            CopyMat<RealType, BlockSize>(&panelA[PanelSizeK*ib], PanelSizeK, &A[(ib+ip)*size + kp], size);
                        }

                        for(size_t idxRow = 0 ; idxRow < BlockSize ; ++idxRow){
                            for(size_t idxCol = 0 ; idxCol < BlockSize ; ++idxCol){
                                RealType sum = 0;
                                for(size_t idxK = 0 ; idxK < PanelSizeK ; ++idxK){
                                    sum += panelA[(idxRow+ib)*PanelSizeK+ idxK]
                                            * panelB[idxCol*PanelSizeK+ idxK];
                                }
                                C[(jp+jb+idxCol)*size + ip + ib + idxRow] += sum;
                            }
                        }
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////
///
/// InaVec with a blocking scheme more close Goto-blas
///
///////////////////////////////////////////////////////////////////////////////////

template <class RealType, size_t PanelSizeK, size_t PanelSizeiA,
          size_t PanelSizejB, class VecType>
void ScalarGemmIna(const RealType* __restrict__ A, const RealType* __restrict__ B,
                RealType* __restrict__ C, const size_t size){

    const int BlockSize = VecType::VecLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK/BlockSize)*BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA/BlockSize)*BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB/BlockSize)*BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size/PanelSizeK)*PanelSizeK == size);
    assert((size/PanelSizeiA)*PanelSizeiA == size);
    assert((size/PanelSizejB)*PanelSizejB == size);

    for(size_t ip = 0 ; ip < size ; ip += PanelSizeiA){
        for(size_t jp = 0 ; jp < size ; jp += PanelSizejB){

            for(size_t kp = 0 ; kp < size ; kp += PanelSizeK){

                alignas(64) RealType panelA[PanelSizeiA*PanelSizeK];
                alignas(64) RealType panelB[PanelSizeK*BlockSize];

                for(size_t jb = 0 ; jb < PanelSizejB ; jb += BlockSize){

                    CopyMat<RealType, BlockSize>(panelB, PanelSizeK, &B[jp*size + kp], size);

                    for(size_t ib = 0 ; ib < PanelSizeiA ; ib += BlockSize){

                        if(jb == 0){
                            CopyMat<RealType, BlockSize>(&panelA[PanelSizeK*ib], PanelSizeK, &A[(ib+ip)*size + kp], size);
                        }

                        for(size_t idxRow = 0 ; idxRow < BlockSize ; ++idxRow){
                            for(size_t idxCol = 0 ; idxCol < BlockSize ; ++idxCol){
                                VecType sum = 0;
                                for(size_t idxK = 0 ; idxK < PanelSizeK ; idxK += BlockSize){
                                    sum += VecType(&panelA[(idxRow+ib)*PanelSizeK+ idxK])
                                            * VecType(&panelB[idxCol*PanelSizeK+ idxK]);
                                }
                                C[(jp+jb+idxCol)*size + ip + ib + idxRow] += sum.horizontalSum();
                            }
                        }
                    }
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////
///
/// InaVec with a blocking scheme more close Goto-blas
/// but with more instructions
///
///////////////////////////////////////////////////////////////////////////////////

template <class RealType, size_t PanelSizeK, size_t PanelSizeiA,
          size_t PanelSizejB, class VecType>
void ScalarGemmInaV2(const RealType* __restrict__ A, const RealType* __restrict__ B,
                RealType* __restrict__ C, const size_t size){

    const int BlockSize = VecType::VecLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK/BlockSize)*BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA/BlockSize)*BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB/BlockSize)*BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size/PanelSizeK)*PanelSizeK == size);
    assert((size/PanelSizeiA)*PanelSizeiA == size);
    assert((size/PanelSizejB)*PanelSizejB == size);

    for(size_t ip = 0 ; ip < size ; ip += PanelSizeiA){
        for(size_t jp = 0 ; jp < size ; jp += PanelSizejB){

            for(size_t kp = 0 ; kp < size ; kp += PanelSizeK){

                alignas(64) RealType panelA[PanelSizeiA*PanelSizeK];
                alignas(64) RealType panelB[PanelSizeK*BlockSize];

                for(size_t jb = 0 ; jb < PanelSizejB ; jb += BlockSize){

                    CopyMat<RealType, BlockSize>(panelB, PanelSizeK, &B[jp*size + kp], size);

                    for(size_t ib = 0 ; ib < PanelSizeiA ; ib += BlockSize){

                        if(jb == 0){
                            CopyMatT<RealType, BlockSize>(&panelA[ib], PanelSizeiA, PanelSizeK,
                                                          &A[(ib+ip)*size + kp], size);
                        }

                        VecType sum[BlockSize];

                        for(size_t idxCol = 0 ; idxCol < BlockSize ; ++idxCol){
                            sum[idxCol] = 0.;
                        }

                        for(size_t idxK = 0 ; idxK < PanelSizeK ; ++idxK){
                            const VecType valA(&panelA[idxK*PanelSizeiA + ib]);
                            for(size_t idxCol = 0 ; idxCol < BlockSize ; ++idxCol){
                                sum[idxCol] += valA * VecType(panelB[idxCol*PanelSizeK + idxK]);
                            }
                        }

                        RealType* __restrict__ ptrC = &C[(jp+jb)*size + ip + ib];
                        for(size_t idxCol = 0 ; idxCol < BlockSize ; ++idxCol){
                            VecType res = sum[idxCol] + VecType(&ptrC[idxCol*size]);
                            res.storeInArray(&ptrC[idxCol*size]);
                        }
                    }
                }
            }
        }
    }
}

template <class RealType>
void compareGemmTime(){
    const size_t matDim = InaVecBestType<RealType>::VecLength * 128 * 2;
    const size_t nbFlops = matDim * matDim * matDim * 2;

    std::unique_ptr< RealType[] > A(new RealType[matDim*matDim]);
    std::unique_ptr< RealType[] > B(new RealType[matDim*matDim]);

    for(size_t idxVal = 0 ; idxVal < matDim*matDim ; ++idxVal){
        A[idxVal] = B[idxVal] = RealType((idxVal+1)%matDim);
    }

    /////////////////////////////////////////////////////////////

    std::unique_ptr< RealType[] > CScalarNoBlock(new RealType[matDim*matDim]);
    memset(CScalarNoBlock.get(), 0, sizeof(RealType)*matDim*matDim);
    if(matDim < 1000){

        InaTimer timer;

        ScalarGemmNoBlock<RealType>(A.get(), B.get(), CScalarNoBlock.get(), matDim);

        timer.stop();
        std::cout << "Scalar no block for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(nbFlops)/timer.getElapsed())/1E9 << "GFlop/s)\n";
    }
    else{
        std::cout << "Scalar no block for size skipped - matrix too large which makes this algorithm too slow!\n";
    }

    /////////////////////////////////////////////////////////////

    {
        std::unique_ptr< RealType[] > CScalar(new RealType[matDim*matDim]);
        memset(CScalar.get(), 0, sizeof(RealType)*matDim*matDim);

        InaTimer timer;

        ScalarGemm<RealType, InaVecBestType<RealType>::VecLength>(A.get(), B.get(), CScalar.get(), matDim);

        timer.stop();
        std::cout << "Scalar for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(nbFlops)/timer.getElapsed())/1E9 << "GFlop/s)\n";

    }

    /////////////////////////////////////////////////////////////

    const size_t PanelSizeA = 32;
    const size_t PanelSizeB = 32;
    const size_t PanelSizeK = 64;
    {
        std::unique_ptr< RealType[] > CScalar(new RealType[matDim*matDim]);
        memset(CScalar.get(), 0, sizeof(RealType)*matDim*matDim);

        InaTimer timer;

        ScalarGemmV2<RealType, PanelSizeK, PanelSizeA, PanelSizeB,
                InaVecBestType<RealType>::VecLength>(A.get(), B.get(), CScalar.get(), matDim);

        timer.stop();
        std::cout << "Scalar for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(nbFlops)/timer.getElapsed())/1E9 << "GFlop/s)\n";
    }

    /////////////////////////////////////////////////////////////

    {
        std::unique_ptr< RealType[] > CIna(new RealType[matDim*matDim]);
        memset(CIna.get(), 0, sizeof(RealType)*matDim*matDim);

        InaTimer timer;

        ScalarGemmIna<RealType, PanelSizeK, PanelSizeA, PanelSizeB,
                InaVecBestType<RealType>>(A.get(), B.get(), CIna.get(), matDim);

        timer.stop();
        std::cout << "Vector " << InaVecBestType<RealType>::GetName() << " for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(nbFlops)/timer.getElapsed())/1E9 << "GFlop/s)\n";
    }

    /////////////////////////////////////////////////////////////

    {
        std::unique_ptr< RealType[] > CIna(new RealType[matDim*matDim]);
        memset(CIna.get(), 0, sizeof(RealType)*matDim*matDim);

        InaTimer timer;

        ScalarGemmInaV2<RealType, PanelSizeK, PanelSizeA, PanelSizeB,
                InaVecBestType<RealType>>(A.get(), B.get(), CIna.get(), matDim);

        timer.stop();
        std::cout << "Vector V2 " << InaVecBestType<RealType>::GetName() << " for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(nbFlops)/timer.getElapsed())/1E9 << "GFlop/s)\n";
    }
}

int main(int /*argc*/, char* /*argv*/ []) {
    std::cout << "In Float:" << std::endl;
    compareGemmTime<float>();

    std::cout << "In Double:" << std::endl;
    compareGemmTime<double>();

    return 0;
}
