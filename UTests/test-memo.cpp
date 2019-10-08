///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

#include "InastempGlobal.h"

#include "Common/InaMemory.hpp"

#include "UTester.hpp"

#include <cmath>
#include <cstring>
#include <vector>

class TestMemo : public UTester< TestMemo > {
    using Parent = UTester< TestMemo >;

    template <long int alignement>
    void TestBasic() {
        {
        
        }
        for(size_t size = 1 ; size < 10 ; ++size){
            {
                char* vecChar = InaMemory::_new<alignement, char>(size);
                UASSERTETRUE(vecChar != nullptr);
                UASSERTETRUE((reinterpret_cast<unsigned long>(vecChar) & (alignement-1)) == 0);
                InaMemory::_delete(vecChar);   
            }
            {
                double* vecDouble = InaMemory::_new<alignement, double>(size);
                UASSERTETRUE(vecDouble != nullptr);
                UASSERTETRUE((reinterpret_cast<unsigned long>(vecDouble) & (alignement-1)) == 0);
                InaMemory::_delete(vecDouble);   
            }
            {
                std::vector<double>* vecDouble = InaMemory::_new<alignement, std::vector<double> >(size);
                UASSERTETRUE(vecDouble != nullptr);
                UASSERTETRUE((reinterpret_cast<unsigned long>(vecDouble) & (alignement-1)) == 0);
                InaMemory::_delete(vecDouble);   
            }
        }
    }

    void SetTests() {
        Parent::AddTest(&TestMemo::TestBasic<1>, "Basic memo test 1");
        Parent::AddTest(&TestMemo::TestBasic<2>, "Basic memo test 2");
        Parent::AddTest(&TestMemo::TestBasic<4>, "Basic memo test 4");
        Parent::AddTest(&TestMemo::TestBasic<8>, "Basic memo test 8");
        Parent::AddTest(&TestMemo::TestBasic<16>, "Basic memo test 16");
        Parent::AddTest(&TestMemo::TestBasic<32>, "Basic memo test 32");
    }
};


int main() {
    return TestMemo().Run();
}
