///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
//
//
// This file ask the cpuid to get access to CPU properties.
// The file contains 3 mains parts:
// × First part is a wrapper in case we are on Windows or Linux
// × Second part is several call to the function and fill of an list
// × Third is out of scope, it prints the state and the properties
//   in a strict format in order to post process
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <inttypes.h>

int main(void) {
// https://github.com/cirosantilli/linux-kernel-module-cheat/blob/29fd625f3fda79f5e0ee6cac43517ba74340d513/baremetal/arch/aarch64/dump_regs.c#L17
    uint64_t id_aa64pfr0_el1;
    __asm__ ("mrs %0, id_aa64pfr0_el1" : "=r" (id_aa64pfr0_el1) : :);
    std::cout << "SVE=" << ((id_aa64pfr0_el1 >> 32) & 0xF?"TRUE":"FALSE") << ";";

    return 0;    
}

