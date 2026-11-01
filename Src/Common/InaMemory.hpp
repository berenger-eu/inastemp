///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAMEMORY_HPP
#define INAMEMORY_HPP

#include <cstddef>
#include <cstdint>
#include <cassert>
#include <limits>
#include <memory>
#include <new>

/**
 * Proposes methods to allocate aligned memory.
 * Any memory from this class must be dellocated with
 * the Delete method.
 * Note: the alignement must be a power of 2
 */
class InaMemory {
    struct Header{
        unsigned char* originPtr;
        std::size_t nbElements;
    };
    const static std::size_t HeaderSize = sizeof(Header);

    template < std::size_t Alignement >
    static void* Allocate(const std::size_t inSize, const std::size_t inNbelements) {
        // Return null for empty allocation
        if (inSize == 0) {
            return nullptr;
        }

        // Ensure it is a power of 2 > 0
        static_assert(Alignement != 0 && ((Alignement - 1) & Alignement) == 0, "Alignement must be a power of 2");

        // Keep the header aligned even when the requested alignment is small.
        constexpr std::size_t EffectiveAlignment =
            Alignement < alignof(Header) ? alignof(Header) : Alignement;
        constexpr std::size_t Offset = HeaderSize + EffectiveAlignment - 1;
        if (inSize > std::numeric_limits<std::size_t>::max() - Offset) {
            throw std::bad_array_new_length();
        }

        unsigned char* allocatedMemoryPtr = new unsigned char[inSize + Offset];
        void* alignedMemoryPtr = allocatedMemoryPtr + HeaderSize;
        std::size_t space = inSize + EffectiveAlignment - 1;
        std::align(EffectiveAlignment, inSize, alignedMemoryPtr, space);
        assert((reinterpret_cast<std::uintptr_t>(alignedMemoryPtr) & (Alignement - 1)) == 0);

        // Start the header's lifetime in the reserved storage.
        new (static_cast<unsigned char*>(alignedMemoryPtr) - HeaderSize)
            Header{allocatedMemoryPtr, inNbelements};

        // Return aligned address
        return reinterpret_cast< void* >(alignedMemoryPtr);
    }

public:
    template < std::size_t Alignement, class ObjectType >
    static ObjectType* _new(const std::size_t inNbElementsInArray) {
        if (inNbElementsInArray == 0) {
            return nullptr;
        }

        if (inNbElementsInArray > std::numeric_limits<std::size_t>::max() / sizeof(ObjectType)) {
            throw std::bad_array_new_length();
        }
        const std::size_t sizeInBytes = inNbElementsInArray * sizeof(ObjectType);
        constexpr std::size_t EffectiveAlignment =
            Alignement < alignof(ObjectType) ? alignof(ObjectType) : Alignement;

        ObjectType* alignedArray = reinterpret_cast< ObjectType* >(Allocate<EffectiveAlignment>(sizeInBytes, inNbElementsInArray));

        static_assert(Alignement != 0 && ((Alignement - 1) & Alignement) == 0,
                      "Alignement must be a power of 2");
        std::size_t constructed = 0;
        try {
            for (; constructed < inNbElementsInArray; ++constructed) {
                new (static_cast<void*>(&alignedArray[constructed])) ObjectType();
            }
        } catch (...) {
            while (constructed != 0) {
                alignedArray[--constructed].~ObjectType();
            }
            const Header* headerPtr = reinterpret_cast<const Header*>(
                reinterpret_cast<const unsigned char*>(alignedArray) - HeaderSize);
            delete[] headerPtr->originPtr;
            throw;
        }
        return alignedArray;
    }

    template < class ObjectType >
    static void _delete(const ObjectType* ptrToFree) {
        if (ptrToFree) {
            const Header* headerPtr = reinterpret_cast<const Header*>(reinterpret_cast< const unsigned char* >(ptrToFree) - HeaderSize);
            const std::size_t numberOfElements = headerPtr->nbElements;

            for (std::size_t idx = numberOfElements; idx != 0; --idx) {
                ptrToFree[idx - 1].~ObjectType();
            }

            delete[] headerPtr->originPtr;
        }
    }
};


#endif // INAMEMORY_HPP
