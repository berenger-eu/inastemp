
#ifndef __ARM_FEATURE_SVE
#error  __ARM_FEATURE_SVE is undefined
#endif


#include <arm_sve.h>


void daxpy_1_1(int64_t n, double da, double *dx, double *dy)
{
        int64_t i = 0;
        svbool_t pg = svwhilelt_b64(i, n);					     // [1]
        do
                {
                        svfloat64_t dx_vec = svld1(pg, &dx[i]);                     // [2]
                        svfloat64_t dy_vec = svld1(pg, &dy[i]);                     // [2]
                        svst1(pg, &dy[i], svmla_x(pg, dy_vec, dx_vec, da));         // [3]
                        i += svcntd();                                              // [4]
                        pg = svwhilelt_b64(i, n);                                   // [1]
                }
        while (svptest_any(svptrue_b64(), pg));                                   // [5]
}

int main(){
	return 0;
}
