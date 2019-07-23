///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAMEMORY_HPP
#define INAMEMORY_HPP

#include <cstdint>
#include <cassert>

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

        // We will need to store the adress of the real blocks
        const std::size_t Offset = (Alignement < HeaderSize ? (HeaderSize / Alignement) + 1 : 1) * Alignement;

        unsigned char* allocatedMemoryPtr = new unsigned char[inSize + Offset];
        unsigned char* alignedMemoryPtr   = reinterpret_cast< unsigned char* >((reinterpret_cast< std::size_t >(allocatedMemoryPtr) + Offset) & ~(Alignement - 1));
        assert((reinterpret_cast< std::size_t >(alignedMemoryPtr) & (Alignement-1)) == 0);
        
        Header* headerPtr          = reinterpret_cast<Header*>(alignedMemoryPtr - HeaderSize);

        // Save real address then header value then alignement
        headerPtr->originPtr = allocatedMemoryPtr;
        headerPtr->nbElements      = inNbelements;

        // Return aligned address
        return reinterpret_cast< void* >(alignedMemoryPtr);
    }

public:
    template < std::size_t Alignement, class ObjectType >
    static ObjectType* _new(const std::size_t inNbElementsInArray) {
        if (inNbElementsInArray == 0) {
            return nullptr;
        }

        const std::size_t sizeInBytes = (inNbElementsInArray * sizeof(ObjectType));

        ObjectType* alignedArray = reinterpret_cast< ObjectType* >(Allocate<Alignement>(sizeInBytes, inNbElementsInArray));

        for (std::size_t idx = 0; idx < inNbElementsInArray; ++idx) {
            new (&alignedArray[idx]) ObjectType();
        }
        return alignedArray;
    }

    template < class ObjectType >
    static void _delete(const ObjectType* ptrToFree) {
        if (ptrToFree) {
            const Header* headerPtr = reinterpret_cast<const Header*>(reinterpret_cast< const unsigned char* >(ptrToFree) - HeaderSize);
            const std::size_t numberOfElements = headerPtr->nbElements;

            for (std::size_t idx = 0; idx < numberOfElements; ++idx) {
                ptrToFree[idx].~ObjectType();
            }

            delete[] headerPtr->originPtr;
        }
    }
};


#endif // INAMEMORY_HPP
