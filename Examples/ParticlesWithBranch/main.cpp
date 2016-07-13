///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

// In this example we time the duration to compute a given number of Exponential
// We compare the scalar version and the BestType

#include "InastempConfig.h"
#include "SCALAR/InaVecSCALARDouble.hpp"
#include "SCALAR/InaVecSCALARFloat.hpp"
#include "Common/InaTimer.hpp"

#ifdef INASTEMP_USE_SSE3
#include "SSE3/InaVecSSE3Double.hpp"
#include "SSE3/InaVecSSE3Float.hpp"
#endif

#ifdef INASTEMP_USE_AVX2
#include "AVX2/InaVecAVX2Double.hpp"
#include "AVX2/InaVecAVX2Float.hpp"
#include <immintrin.h>
#endif

#include <memory>
#include <iostream>
#include <random>

//! Here it is a dummy example where we compute the interactions between particles.
//! The Kernel is (m/r) but to show how it goes with a branch
//! we ask that when r < cutDistance the kernel becomes (m-constantIfCut)/r
//!
//! We compare :
//! - a scalar kernel
//! - an AVX2 by hand using intrinsics - if AVX available
//! - a kernel based on the AVX2 class from inastemp  - if AVX available
//! - an SSE3 by hand using intrinsics - if SSE3 available
//! - a kernel based on the SSE3 class from inastemp  - if SSE3 available
//! - a kernel based on the best available vectorizer from inastemp (might be AVX2!)



//////////////////////////////////////////////////////////////////////////////////////////////
/// SCALAR usual functions
//////////////////////////////////////////////////////////////////////////////////////////////

void ScalarFunction(const size_t nbParticles, const double* __restrict__ positionsX,
                        const double* __restrict__ positionsY,const double* __restrict__ positionsZ,
                        const double* __restrict__ physicalValues, double* __restrict__ potentials,
                        const double cutDistance, const double constantIfCut){
    for(size_t idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){
        for(size_t idxSource = idxTarget+1 ; idxSource < nbParticles ; ++idxSource){
            const double dx = positionsX[idxTarget] - positionsX[idxSource];
            const double dy = positionsY[idxTarget] - positionsY[idxSource];
            const double dz = positionsZ[idxTarget] - positionsZ[idxSource];

            const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            const double inv_distance = 1/distance;

            if(distance < cutDistance){
                potentials[idxTarget]  += ( inv_distance * physicalValues[idxSource] );
                potentials[idxSource] += ( inv_distance * physicalValues[idxTarget] );
            }
            else{
                potentials[idxTarget]  += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
                potentials[idxSource] += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
            }
        }
    }
}

void ScalarFunction(const size_t nbParticles, const float* __restrict__ positionsX,
                        const float* __restrict__ positionsY,const float* __restrict__ positionsZ,
                        const float* __restrict__ physicalValues, float* __restrict__ potentials,
                        const float cutDistance, const float constantIfCut){
    for(size_t idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){
        for(size_t idxSource = idxTarget+1 ; idxSource < nbParticles ; ++idxSource){
            const float dx = positionsX[idxTarget] - positionsX[idxSource];
            const float dy = positionsY[idxTarget] - positionsY[idxSource];
            const float dz = positionsZ[idxTarget] - positionsZ[idxSource];

            const float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            const float inv_distance = 1/distance;

            if(distance < cutDistance){
                potentials[idxTarget]  += ( inv_distance * physicalValues[idxSource] );
                potentials[idxSource] += ( inv_distance * physicalValues[idxTarget] );
            }
            else{
                potentials[idxTarget]  += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
                potentials[idxSource] += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////
/// AVX2 functions
//////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INASTEMP_USE_AVX2

void HandVectorizedFunctionAVX(const size_t nbParticles, const double* __restrict__ positionsX,
                        const double* __restrict__ positionsY,const double* __restrict__ positionsZ,
                        const double* __restrict__ physicalValues, double* __restrict__ potentials,
                        const double cutDistance, const double constantIfCut){

    const size_t VecLength = 4; // sizeof(__m256d)/sizeof(double)

    const __m256d VecOne = _mm256_set1_pd(1);
    const __m256d VecConstantIfCut = _mm256_set1_pd(constantIfCut);
    const __m256d VecCutDistance = _mm256_set1_pd(cutDistance);

    for(size_t idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){

        const __m256d targetX = _mm256_set1_pd(positionsX[idxTarget]);
        const __m256d targetY = _mm256_set1_pd(positionsY[idxTarget]);
        const __m256d targetZ = _mm256_set1_pd(positionsZ[idxTarget]);

        const __m256d targetPhysicalValue = _mm256_set1_pd(physicalValues[idxTarget]);
        __m256d targetPotential = _mm256_setzero_pd();

        const size_t lastToCompute = ((nbParticles-(idxTarget+1))/VecLength)*VecLength+(idxTarget+1);

        for(size_t idxSource = idxTarget+1 ; idxSource < lastToCompute ; idxSource += VecLength){
            const __m256d dx = targetX - _mm256_loadu_pd(&positionsX[idxSource]);
            const __m256d dy = targetY - _mm256_loadu_pd(&positionsY[idxSource]);
            const __m256d dz = targetZ - _mm256_loadu_pd(&positionsZ[idxSource]);

            const __m256d distance = _mm256_sqrt_pd(dx*dx + dy*dy + dz*dz);
            const __m256d inv_distance = VecOne/distance;

            const __m256d testRes = _mm256_cmp_pd(distance, VecCutDistance, _CMP_LT_OQ);

            const __m256d sourcesPhysicalValue = _mm256_loadu_pd(&physicalValues[idxSource]);

            targetPotential += inv_distance * _mm256_or_pd(_mm256_and_pd(testRes, sourcesPhysicalValue),
                                                           _mm256_andnot_pd(testRes, sourcesPhysicalValue-VecConstantIfCut));
            const __m256d resSource = inv_distance * _mm256_or_pd(_mm256_and_pd(testRes, targetPhysicalValue),
                                                                  _mm256_andnot_pd(testRes, targetPhysicalValue-VecConstantIfCut));
            const __m256d currentSource = _mm256_loadu_pd(&potentials[idxSource]);
            _mm256_storeu_pd(&potentials[idxSource], resSource+currentSource);
        }

        potentials[idxTarget] += InaVecAVX2<double>(targetPotential).horizontalSum();

        /////////////////////////////////////////////////////////////////////////////////

        for(size_t idxSource = lastToCompute ; idxSource < nbParticles ; ++idxSource){
            const double dx = positionsX[idxTarget] - positionsX[idxSource];
            const double dy = positionsY[idxTarget] - positionsY[idxSource];
            const double dz = positionsZ[idxTarget] - positionsZ[idxSource];

            const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            const double inv_distance = 1/distance;

            if(distance < cutDistance){
                potentials[idxTarget] += ( inv_distance * physicalValues[idxSource] );
                potentials[idxSource] += ( inv_distance * physicalValues[idxTarget] );
            }
            else{
                potentials[idxTarget] += ( inv_distance * (physicalValues[idxSource]-constantIfCut) );
                potentials[idxSource] += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
            }
        }
    }
}


void HandVectorizedFunctionAVX(const size_t nbParticles, const float* __restrict__ positionsX,
                        const float* __restrict__ positionsY,const float* __restrict__ positionsZ,
                        const float* __restrict__ physicalValues, float* __restrict__ potentials,
                        const float cutDistance, const float constantIfCut){

    const size_t VecLength = 8; // sizeof(__m256)/sizeof(float)

    const __m256 VecOne = _mm256_set1_ps(1);
    const __m256 VecConstantIfCut = _mm256_set1_ps(constantIfCut);
    const __m256 VecCutDistance = _mm256_set1_ps(cutDistance);

    for(size_t idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){

        const __m256 targetX = _mm256_set1_ps(positionsX[idxTarget]);
        const __m256 targetY = _mm256_set1_ps(positionsY[idxTarget]);
        const __m256 targetZ = _mm256_set1_ps(positionsZ[idxTarget]);

        const __m256 targetPhysicalValue = _mm256_set1_ps(physicalValues[idxTarget]);
        __m256 targetPotential = _mm256_setzero_ps();

        const size_t lastToCompute = ((nbParticles-(idxTarget+1))/VecLength)*VecLength+(idxTarget+1);

        for(size_t idxSource = idxTarget+1 ; idxSource < lastToCompute ; idxSource += VecLength){
            const __m256 dx = targetX - _mm256_loadu_ps(&positionsX[idxSource]);
            const __m256 dy = targetY - _mm256_loadu_ps(&positionsY[idxSource]);
            const __m256 dz = targetZ - _mm256_loadu_ps(&positionsZ[idxSource]);

            const __m256 distance = _mm256_sqrt_ps(dx*dx + dy*dy + dz*dz);
            const __m256 inv_distance = VecOne/distance;

            const __m256 testRes = _mm256_cmp_ps(distance, VecCutDistance, _CMP_LT_OQ);

            const __m256 sourcesPhysicalValue = _mm256_loadu_ps(&physicalValues[idxSource]);

            targetPotential += inv_distance * _mm256_or_ps(_mm256_and_ps(testRes, sourcesPhysicalValue),
                                                           _mm256_andnot_ps(testRes, sourcesPhysicalValue-VecConstantIfCut));
            const __m256 resSource = inv_distance * _mm256_or_ps(_mm256_and_ps(testRes, targetPhysicalValue),
                                                                  _mm256_andnot_ps(testRes, targetPhysicalValue-VecConstantIfCut));
            const __m256 currentSource = _mm256_loadu_ps(&potentials[idxSource]);
            _mm256_storeu_ps(&potentials[idxSource], resSource+currentSource);
        }

        potentials[idxTarget] += InaVecAVX2<float>(targetPotential).horizontalSum();

        /////////////////////////////////////////////////////////////////////////////////

        for(size_t idxSource = lastToCompute ; idxSource < nbParticles ; ++idxSource){
            const float dx = positionsX[idxTarget] - positionsX[idxSource];
            const float dy = positionsY[idxTarget] - positionsY[idxSource];
            const float dz = positionsZ[idxTarget] - positionsZ[idxSource];

            const float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            const float inv_distance = 1/distance;

            if(distance < cutDistance){
                potentials[idxTarget] += ( inv_distance * physicalValues[idxSource] );
                potentials[idxSource] += ( inv_distance * physicalValues[idxTarget] );
            }
            else{
                potentials[idxTarget] += ( inv_distance * (physicalValues[idxSource]-constantIfCut) );
                potentials[idxSource] += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
            }
        }
    }
}

#endif



//////////////////////////////////////////////////////////////////////////////////////////////
/// SSE3 functions
//////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INASTEMP_USE_SSE3

void HandVectorizedFunctionSSE3(const size_t nbParticles, const double* __restrict__ positionsX,
                        const double* __restrict__ positionsY,const double* __restrict__ positionsZ,
                        const double* __restrict__ physicalValues, double* __restrict__ potentials,
                        const double cutDistance, const double constantIfCut){

    const size_t VecLength = 2; // sizeof(__m128)/sizeof(double)

    const __m128d VecOne = _mm_set1_pd(1);
    const __m128d VecConstantIfCut = _mm_set1_pd(constantIfCut);
    const __m128d VecCutDistance = _mm_set1_pd(cutDistance);

    for(size_t idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){

        const __m128d targetX = _mm_set1_pd(positionsX[idxTarget]);
        const __m128d targetY = _mm_set1_pd(positionsY[idxTarget]);
        const __m128d targetZ = _mm_set1_pd(positionsZ[idxTarget]);

        const __m128d targetPhysicalValue = _mm_set1_pd(physicalValues[idxTarget]);
        __m128d targetPotential = _mm_setzero_pd();

        const size_t lastToCompute = ((nbParticles-(idxTarget+1))/VecLength)*VecLength+(idxTarget+1);

        for(size_t idxSource = idxTarget+1 ; idxSource < lastToCompute ; idxSource += VecLength){
            const __m128d dx = targetX - _mm_loadu_pd(&positionsX[idxSource]);
            const __m128d dy = targetY - _mm_loadu_pd(&positionsY[idxSource]);
            const __m128d dz = targetZ - _mm_loadu_pd(&positionsZ[idxSource]);

            const __m128d distance = _mm_sqrt_pd(dx*dx + dy*dy + dz*dz);
            const __m128d inv_distance = VecOne/distance;

            const __m128d testRes = _mm_cmp_pd(distance, VecCutDistance, _CMP_LT_OQ);

            const __m128d sourcesPhysicalValue = _mm_loadu_pd(&physicalValues[idxSource]);

            targetPotential += inv_distance * _mm_or_pd(_mm_and_pd(testRes, sourcesPhysicalValue),
                                                           _mm_andnot_pd(testRes, sourcesPhysicalValue-VecConstantIfCut));
            const __m128d resSource = inv_distance * _mm_or_pd(_mm_and_pd(testRes, targetPhysicalValue),
                                                                  _mm_andnot_pd(testRes, targetPhysicalValue-VecConstantIfCut));
            const __m128d currentSource = _mm_loadu_pd(&potentials[idxSource]);
            _mm_storeu_pd(&potentials[idxSource], resSource+currentSource);
        }

        potentials[idxTarget] += InaVecSSE3<double>(targetPotential).horizontalSum();

        /////////////////////////////////////////////////////////////////////////////////

        for(size_t idxSource = lastToCompute ; idxSource < nbParticles ; ++idxSource){
            const double dx = positionsX[idxTarget] - positionsX[idxSource];
            const double dy = positionsY[idxTarget] - positionsY[idxSource];
            const double dz = positionsZ[idxTarget] - positionsZ[idxSource];

            const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            const double inv_distance = 1/distance;

            if(distance < cutDistance){
                potentials[idxTarget] += ( inv_distance * physicalValues[idxSource] );
                potentials[idxSource] += ( inv_distance * physicalValues[idxTarget] );
            }
            else{
                potentials[idxTarget] += ( inv_distance * (physicalValues[idxSource]-constantIfCut) );
                potentials[idxSource] += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
            }
        }
    }
}

void HandVectorizedFunctionSSE3(const size_t nbParticles, const float* __restrict__ positionsX,
                        const float* __restrict__ positionsY,const float* __restrict__ positionsZ,
                        const float* __restrict__ physicalValues, float* __restrict__ potentials,
                        const float cutDistance, const float constantIfCut){

    const size_t VecLength = 4; // sizeof(__m128)/sizeof(float)

    const __m128 VecOne = _mm_set1_ps(1);
    const __m128 VecConstantIfCut = _mm_set1_ps(constantIfCut);
    const __m128 VecCutDistance = _mm_set1_ps(cutDistance);

    for(size_t idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){

        const __m128 targetX = _mm_set1_ps(positionsX[idxTarget]);
        const __m128 targetY = _mm_set1_ps(positionsY[idxTarget]);
        const __m128 targetZ = _mm_set1_ps(positionsZ[idxTarget]);

        const __m128 targetPhysicalValue = _mm_set1_ps(physicalValues[idxTarget]);
        __m128 targetPotential = _mm_setzero_ps();

        const size_t lastToCompute = ((nbParticles-(idxTarget+1))/VecLength)*VecLength+(idxTarget+1);

        for(size_t idxSource = idxTarget+1 ; idxSource < lastToCompute ; idxSource += VecLength){
            const __m128 dx = targetX - _mm_loadu_ps(&positionsX[idxSource]);
            const __m128 dy = targetY - _mm_loadu_ps(&positionsY[idxSource]);
            const __m128 dz = targetZ - _mm_loadu_ps(&positionsZ[idxSource]);

            const __m128 distance = _mm_sqrt_ps(dx*dx + dy*dy + dz*dz);
            const __m128 inv_distance = VecOne/distance;

            const __m128 testRes = _mm_cmp_ps(distance, VecCutDistance, _CMP_LT_OQ);

            const __m128 sourcesPhysicalValue = _mm_loadu_ps(&physicalValues[idxSource]);

            targetPotential += inv_distance * _mm_or_ps(_mm_and_ps(testRes, sourcesPhysicalValue),
                                                           _mm_andnot_ps(testRes, sourcesPhysicalValue-VecConstantIfCut));
            const __m128 resSource = inv_distance * _mm_or_ps(_mm_and_ps(testRes, targetPhysicalValue),
                                                                  _mm_andnot_ps(testRes, targetPhysicalValue-VecConstantIfCut));
            const __m128 currentSource = _mm_loadu_ps(&potentials[idxSource]);
            _mm_storeu_ps(&potentials[idxSource], resSource+currentSource);
        }

        potentials[idxTarget] += InaVecSSE3<float>(targetPotential).horizontalSum();

        /////////////////////////////////////////////////////////////////////////////////

        for(size_t idxSource = lastToCompute ; idxSource < nbParticles ; ++idxSource){
            const float dx = positionsX[idxTarget] - positionsX[idxSource];
            const float dy = positionsY[idxTarget] - positionsY[idxSource];
            const float dz = positionsZ[idxTarget] - positionsZ[idxSource];

            const float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            const float inv_distance = 1/distance;

            if(distance < cutDistance){
                potentials[idxTarget] += ( inv_distance * physicalValues[idxSource] );
                potentials[idxSource] += ( inv_distance * physicalValues[idxTarget] );
            }
            else{
                potentials[idxTarget] += ( inv_distance * (physicalValues[idxSource]-constantIfCut) );
                potentials[idxSource] += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
            }
        }
    }
}

#endif


//////////////////////////////////////////////////////////////////////////////////////////////
/// inastemp unique function
//////////////////////////////////////////////////////////////////////////////////////////////

template <class VecType, class RealType>
void VectorizedFunction(const size_t nbParticles, const RealType* __restrict__ positionsX,
                        const RealType* __restrict__ positionsY,const RealType* __restrict__ positionsZ,
                        const RealType* __restrict__ physicalValues, RealType* __restrict__ potentials,
                        const RealType cutDistance, const RealType constantIfCut){

    const VecType VecOne = 1;
    const VecType VecConstantIfCut = constantIfCut;
    const VecType VecCutDistance = cutDistance;

    for(size_t idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){

        const VecType targetX = positionsX[idxTarget];
        const VecType targetY = positionsY[idxTarget];
        const VecType targetZ = positionsZ[idxTarget];

        const VecType targetPhysicalValue = physicalValues[idxTarget];
        VecType targetPotential = VecType::GetZero();

        const size_t lastToCompute = ((nbParticles-(idxTarget+1))/VecType::VecLength)*VecType::VecLength+(idxTarget+1);

        for(size_t idxSource = idxTarget+1 ; idxSource < lastToCompute ; idxSource += VecType::VecLength){
            const VecType dx = targetX - VecType(&positionsX[idxSource]);
            const VecType dy = targetY - VecType(&positionsY[idxSource]);
            const VecType dz = targetZ - VecType(&positionsZ[idxSource]);

            const VecType distance = VecType(dx*dx + dy*dy + dz*dz).sqrt();
            const VecType inv_distance = VecOne/distance;

            const typename VecType::MaskType testRes = (distance < VecCutDistance);

            const VecType sourcesPhysicalValue = VecType(&physicalValues[idxSource]);

            targetPotential += inv_distance * VecType::IfElse(testRes, sourcesPhysicalValue,
                                                                      sourcesPhysicalValue-VecConstantIfCut);
            const VecType resSource = inv_distance * VecType::IfElse(testRes, targetPhysicalValue,
                                                                       targetPhysicalValue-VecConstantIfCut);

            /// The previous line are equivalent (in results and performance) to:
            // targetPotential += inv_distance * VecType::If(testRes).Then([&](){
            //     return sourcesPhysicalValue;
            // }).Else([&](){
            //     return sourcesPhysicalValue-VecConstantIfCut;
            // });
            // const VecType resSource = inv_distance * VecType::If(testRes).Then([&](){
            //     return targetPhysicalValue;
            // }).Else([&](){
            //     return targetPhysicalValue-VecConstantIfCut;
            // });

            const VecType currentSource = VecType(&potentials[idxSource]);
            (resSource+currentSource).storeInArray(&potentials[idxSource]);
        }

        potentials[idxTarget] += targetPotential.horizontalSum();

        /////////////////////////////////////////////////////////////////////////////////

        for(size_t idxSource = lastToCompute ; idxSource < nbParticles ; ++idxSource){
            const RealType dx = positionsX[idxTarget] - positionsX[idxSource];
            const RealType dy = positionsY[idxTarget] - positionsY[idxSource];
            const RealType dz = positionsZ[idxTarget] - positionsZ[idxSource];

            const RealType distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            const RealType inv_distance = 1/distance;

            if(distance < cutDistance){
                potentials[idxTarget] += ( inv_distance * physicalValues[idxSource] );
                potentials[idxSource] += ( inv_distance * physicalValues[idxTarget] );
            }
            else{
                potentials[idxTarget] += ( inv_distance * (physicalValues[idxSource]-constantIfCut) );
                potentials[idxSource] += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////
/// Main
//////////////////////////////////////////////////////////////////////////////////////////////

template <class RealType>
void Test(){
    const size_t NbParticles = 20000;
    const RealType cutDistance = RealType(0.5);
    const RealType constantIfCut = RealType(0.05);

    std::cout << "Program starts for N = " << NbParticles << " particles." << std::endl;

    /////////////////////////////////////////////////////////////

    std::mt19937_64 rdEngine;
    std::uniform_real_distribution<RealType> rdDist(0, 1);

    std::unique_ptr< RealType[] > positionsX(new RealType[NbParticles]);
    std::unique_ptr< RealType[] > positionsY(new RealType[NbParticles]);
    std::unique_ptr< RealType[] > positionsZ(new RealType[NbParticles]);
    std::unique_ptr< RealType[] > physicalValues(new RealType[NbParticles]);
    std::unique_ptr< RealType[] > potentials(new RealType[NbParticles]);

    for(size_t idxPart = 0 ; idxPart < NbParticles ; ++idxPart){
        positionsX[idxPart] = rdDist(rdEngine);
        positionsY[idxPart] = rdDist(rdEngine);
        positionsZ[idxPart] = rdDist(rdEngine);
        physicalValues[idxPart] = RealType(0.1);
        potentials[idxPart] = 0;
    }

    /////////////////////////////////////////////////////////////

    {
        InaTimer timer;

        ScalarFunction(NbParticles, positionsX.get(), positionsY.get(), positionsZ.get(),
                           physicalValues.get(), potentials.get(), cutDistance, constantIfCut);
        timer.stop();
        std::cout << "Scalar for " << (NbParticles*NbParticles/2) << " interactions took " << timer.getElapsed() << "s\n";
    }

    /////////////////////////////////////////////////////////////

    {
        InaTimer timer;

        VectorizedFunction<InaVecSCALAR<RealType>>(NbParticles, positionsX.get(), positionsY.get(), positionsZ.get(),
                           physicalValues.get(), potentials.get(), cutDistance, constantIfCut);

        timer.stop();
        std::cout << "InaVecSCALAR<RealType> for " << (NbParticles*NbParticles/2) << " interactions took " << timer.getElapsed() << "s\n";
    }

    /////////////////////////////////////////////////////////////


#ifdef INASTEMP_USE_AVX2
    {
        InaTimer timer;

        HandVectorizedFunctionAVX(NbParticles, positionsX.get(), positionsY.get(), positionsZ.get(),
                           physicalValues.get(), potentials.get(), cutDistance, constantIfCut);

        timer.stop();
        std::cout << "HandVectorizedFunctionAVX for " << (NbParticles*NbParticles/2) << " interactions took " << timer.getElapsed() << "s\n";
    }

    {
        InaTimer timer;

        VectorizedFunction<InaVecAVX2<RealType>>(NbParticles, positionsX.get(), positionsY.get(), positionsZ.get(),
                           physicalValues.get(), potentials.get(), cutDistance, constantIfCut);

        timer.stop();
        std::cout << InaVecAVX2<RealType>().GetName() << " for " << (NbParticles*NbParticles/2) << " interactions took " << timer.getElapsed() << "s\n";
    }
#endif

    /////////////////////////////////////////////////////////////


#ifdef INASTEMP_USE_SSE3
    {
        InaTimer timer;

        HandVectorizedFunctionSSE3(NbParticles, positionsX.get(), positionsY.get(), positionsZ.get(),
                           physicalValues.get(), potentials.get(), cutDistance, constantIfCut);

        timer.stop();
        std::cout << "HandVectorizedFunctionSSE3 for " << (NbParticles*NbParticles/2) << " interactions took " << timer.getElapsed() << "s\n";
    }
    {
        InaTimer timer;

        VectorizedFunction<InaVecSSE3<RealType>>(NbParticles, positionsX.get(), positionsY.get(), positionsZ.get(),
                           physicalValues.get(), potentials.get(), cutDistance, constantIfCut);

        timer.stop();
        std::cout << InaVecSSE3<RealType>().GetName() << " for " << (NbParticles*NbParticles/2) << " interactions took " << timer.getElapsed() << "s\n";
    }
#endif

    {
        InaTimer timer;

        VectorizedFunction<InaVecBestType<RealType>, RealType>(NbParticles, positionsX.get(), positionsY.get(), positionsZ.get(),
                           physicalValues.get(), potentials.get(), cutDistance, constantIfCut);

        timer.stop();
        std::cout << InaVecBestType<RealType>().GetName() << " for " << (NbParticles*NbParticles/2) << " interactions took " << timer.getElapsed() << "s\n";
    }

    /////////////////////////////////////////////////////////////
}

int main(int /*argc*/, char* /*argv*/ []) {

    std::cout << "Test FLOAT:" << std::endl;
    Test<float>();

    std::cout << "Test DOUBLE:" << std::endl;
    Test<double>();

    return 0;
}
