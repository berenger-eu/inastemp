#include <riscv_vector.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>


int main(){
    vint32m1_t acc;
    vint32m1_t op1;
    vint32m1_t op2;
    vint32m1_t res = vmacc_vv_i32m1(acc, op1, op2);

    return 0;
}
