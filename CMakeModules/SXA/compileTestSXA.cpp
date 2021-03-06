#include <velintrin.h>

/*
To test -S will ensure not to link :
/home/berenger/Projects/llvm-dev-install/bin/clang -v -target ve -O3 -fno-vectorize -fno-slp-vectorize -fno-crash-diagnostics main.cpp -S -o /tmp/test.s

To compile:
/home/berenger/Projects/llvm-dev-install/bin/clang -v -target ve -O3 -fno-vectorize -fno-slp-vectorize -fno-crash-diagnostics main.cpp -o /tmp/test.exe
*/


int main(){
    __vr vec = _vel_pvbrd_vsl(4,256);
    __vr res = _vel_vfmuld_vvvl(vec, vec, 256);
    __vr resand = _vel_vand_vvvl(res, res, 10);
    
    return 0;
}

