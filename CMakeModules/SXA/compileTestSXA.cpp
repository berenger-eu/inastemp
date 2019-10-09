#include <velintrin.h>

/*
To test -S will ensure not to link :
/home/berenger/Projects/llvm-dev-install/bin/clang -v -target ve -O3 -fno-vectorize -fno-slp-vectorize -fno-crash-diagnostics main.cpp -S -o /tmp/test.s

To compile:
/home/berenger/Projects/llvm-dev-install/bin/clang -v -target ve -O3 -fno-vectorize -fno-slp-vectorize -fno-crash-diagnostics main.cpp -o /tmp/test.exe
*/


int main(){
    __vr vec = _ve_vbrd_vs_f64(4);
    __vr res = _ve_vfmuld_vvv(vec, vec);
    __vr resand = _vel_vand_vvvl(res, res, 10);
    
    return 0;
}

