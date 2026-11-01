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
#include <limits>
#include <stdexcept>

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

    struct alignas(128) AlignedObject {
        double value = 42.;
    };

    struct TrackedObject {
        inline static int attempts = 0;
        inline static int failAt = -1;
        inline static std::vector<int> destroyed;
        int id;

        TrackedObject() : id(attempts++) {
            if (id == failAt) {
                throw std::runtime_error("constructor failure");
            }
        }
        ~TrackedObject() { destroyed.push_back(id); }
    };

    void TestEmptyAndAlignment() {
        UASSERTETRUE((InaMemory::_new<1, double>(0) == nullptr));
        InaMemory::_delete<double>(nullptr);
        auto* objects = InaMemory::_new<1, AlignedObject>(3);
        UASSERTETRUE(reinterpret_cast<std::uintptr_t>(objects) % alignof(AlignedObject) == 0);
        for (int idx = 0; idx < 3; ++idx) {
            UASSERTEEQUAL(objects[idx].value, 42.);
        }
        InaMemory::_delete(objects);
    }

    void TestOverflow() {
        bool rejected = false;
        try {
            auto* objects = InaMemory::_new<32, double>(
                std::numeric_limits<std::size_t>::max() / sizeof(double) + 1);
            InaMemory::_delete(objects);
        } catch (const std::bad_array_new_length&) {
            rejected = true;
        }
        UASSERTETRUE(rejected);
        rejected = false;
        try {
            auto* objects = InaMemory::_new<32, char>(std::numeric_limits<std::size_t>::max());
            InaMemory::_delete(objects);
        } catch (const std::bad_array_new_length&) {
            rejected = true;
        }
        UASSERTETRUE(rejected);
    }

    void TestObjectLifetime() {
        TrackedObject::attempts = 0;
        TrackedObject::failAt = -1;
        TrackedObject::destroyed.clear();
        auto* objects = InaMemory::_new<1, TrackedObject>(4);
        InaMemory::_delete(objects);
        UASSERTETRUE((TrackedObject::destroyed == std::vector<int>{3, 2, 1, 0}));

        for (int failAt : {0, 2}) {
            TrackedObject::attempts = 0;
            TrackedObject::failAt = failAt;
            TrackedObject::destroyed.clear();
            bool threw = false;
            try {
                objects = InaMemory::_new<1, TrackedObject>(4);
                InaMemory::_delete(objects);
            } catch (const std::runtime_error&) {
                threw = true;
            }
            UASSERTETRUE(threw);
            UASSERTEEQUAL(TrackedObject::attempts, failAt + 1);
            const std::vector<int> expected = failAt == 0 ? std::vector<int>{} : std::vector<int>{1, 0};
            UASSERTETRUE(TrackedObject::destroyed == expected);
        }
    }

    void SetTests() {
        Parent::AddTest(&TestMemo::TestEmptyAndAlignment, "Empty and over-aligned allocations");
        Parent::AddTest(&TestMemo::TestOverflow, "Allocation size overflow");
        Parent::AddTest(&TestMemo::TestObjectLifetime, "Object destruction and constructor exceptions");
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
