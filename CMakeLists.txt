###########################################################################
# Inastemp - Berenger Bramas MPCDF - 2016
# Under MIT Licence, please you must read the LICENCE file.
###########################################################################
cmake_minimum_required (VERSION 2.8.3 FATAL_ERROR)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

#===========================================================================
# Project Declaration
#===========================================================================
project(INASTEMP C CXX)

# check if compiling into source directories
string(COMPARE EQUAL "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" insource)
if(insource)
  message(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/)
set(INASTEMP_CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

# Adds the CMAKE_DEPENDENT_OPTION command
INCLUDE(CMakeDependentOption)

# Release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo" FORCE)
endif()

if($ENV{VERBOSE})
    MESSAGE(STATUS "Main -- CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}")
endif()

# Test if inastemp is included as a subproject
if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(INASTEMP_AS_SUBPROJECT FALSE)
else()
    set(INASTEMP_AS_SUBPROJECT TRUE)
endif()

if($ENV{VERBOSE})
    MESSAGE(STATUS "Main -- INASTEMP_AS_SUBPROJECT ${INASTEMP_AS_SUBPROJECT}")
endif()

# Build only the lib if subproject
if(NOT DEFINED INASTEMP_JUST_LIB)
    set(INASTEMP_JUST_LIB ${INASTEMP_AS_SUBPROJECT} CACHE BOOL
        "Set to on to build the lib for external use (without test and extra flags)" FORCE)
endif()

if($ENV{VERBOSE})
    MESSAGE(STATUS "Main -- INASTEMP_JUST_LIB ${INASTEMP_JUST_LIB}")
endif()

#===========================================================================
# Version Number
#===========================================================================

set(INASTEMP_MAJOR_VERSION 0)
set(INASTEMP_MINOR_VERSION 1)
set(INASTEMP_PATCH_VERSION 0)
set(INASTEMP_VERSION  "${INASTEMP_MAJOR_VERSION}.${INASTEMP_MINOR_VERSION}.${INASTEMP_PATCH_VERSION}" )

#===========================================================================
# Options
#===========================================================================
#  Ask CPU capacities
include(GetCpuInfos)
GetCpuInfos()

#  Ask compiler capacities
include(GetCompilerInfos)
GetCompilerInfos()

# All types from worse to best (ADD-NEW-HERE)
set(ALL_TYPES "SSE3;SSSE3;SSE41;SSE42;AVX;AVX2;AVX512KNL")
set(ALL_USED_TYPES "SCALAR")
set(INASTEMP_USE_SCALAR ON)
set(INASTEMP_CXX_FLAGS "-std=c++11")

# Set custom cpu <=> vec rules (maybe ADD-NEW-HERE if needed)
set(AVX512KNL_CPU_RULES "AVX512ER;AVX512PF;AVX512ER")

foreach(TYPE ${ALL_TYPES})
    if( COMPILER_INFO_${TYPE} )
        if(NOT DEFINED ${TYPE}_CPU_RULES)
            set( ${TYPE}_CPU_RULES "${TYPE}")
        endif()

        set(${TYPE}_SUPPORTED_BY_CPU TRUE)
        foreach(RULE ${${TYPE}_CPU_RULES})
            if(NOT DEFINED CPU_INFO_${RULE})
                MESSAGE(WARNING "Main -- The variable CPU_INFO_${RULE} cannot be read from the current CPU which might be too old. It turns off ${TYPE}.")
                set(${TYPE}_SUPPORTED_BY_CPU FALSE)                
                break()
            elseif(NOT ${CPU_INFO_${RULE}})
                set(${TYPE}_SUPPORTED_BY_CPU FALSE)                
                break()
            endif()
        endforeach()

        OPTION( INASTEMP_USE_${TYPE}  "Set to ON to compile with ${TYPE} support (default based on CPU capacities)" ${${TYPE}_SUPPORTED_BY_CPU} )

        if(${INASTEMP_USE_${TYPE}})
            set(INASTEMP_CXX_FLAGS "${INASTEMP_CXX_FLAGS} ${${TYPE}_FLAGS}")
            set(ALL_USED_TYPES "${ALL_USED_TYPES};${TYPE}")

            if(COMPILER_INFO_${TYPE}_NOOP)
                set(INASTEMP_USE_${TYPE}_OPERATORS ON)
            else()
                set(INASTEMP_USE_${TYPE}_OPERATORS OFF)
            endif()
        endif()
    else()
        if($ENV{VERBOSE})
            MESSAGE(STATUS "Main -- The compiler cannot compile ${TYPE} intrinsics")
        endif()
    endif()
endforeach()

list(LENGTH ALL_USED_TYPES INASTEMP_NB_USED_TYPES)
if($ENV{VERBOSE})
    MESSAGE(STATUS "Main -- INASTEMP_NB_USED_TYPES ${INASTEMP_NB_USED_TYPES}")
endif()
list(GET ALL_USED_TYPES -1 INASTEMP_BESTTYPE )

# Extra-flags
SET( INASTEMP_EXTRA_CXX_FLAGS "-funroll-loops" CACHE STRING "These flags are used for all configurations"  )
set(INASTEMP_CXX_FLAGS  "${INASTEMP_CXX_FLAGS} ${INASTEMP_EXTRA_CXX_FLAGS}")

# Options
if(NOT INASTEMP_AS_SUBPROJECT)
    SET( INASTEMP_EXTRA_CXX_FLAGS_RELEASE "-O3" CACHE STRING "These flags are added to the usual ones for RELEASE"  )
    set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} ${INASTEMP_EXTRA_CXX_FLAGS_RELEASE}")

    SET( INASTEMP_EXTRA_CXX_FLAGS_DEBUG "-O0" CACHE STRING "These flags are added to the usual ones for DEBUG"  )
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} ${INASTEMP_EXTRA_CXX_FLAGS_DEBUG}")

    set(INASTEMP_CXX_FLAGS "${INASTEMP_CXX_FLAGS} -fpic -Wextra -Wnon-virtual-dtor -Wshadow -Wpointer-arith -Wcast-qual -Wconversion  -Wall -Wno-sign-conversion -pedantic -Woverloaded-virtual -Wpointer-arith -Wcast-qual -Wconversion -Wno-error")

    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(INASTEMP_CXX_FLAGS  "${INASTEMP_CXX_FLAGS} -m64")
    endif()

    OPTION( INASTEMP_ATTACHE_SOURCE "Set to on to add -g flag" OFF )
    if( INASTEMP_ATTACHE_SOURCE )
        set(INASTEMP_CXX_FLAGS  "${INASTEMP_CXX_FLAGS} -g")
    endif(INASTEMP_ATTACHE_SOURCE)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "icpc.*$")
  set(INASTEMP_CXX_FLAGS "${INASTEMP_CXX_FLAGS} -diag-disable 2304 -diag-disable 10121 -diag-disable 10120")
else()
  set(INASTEMP_CXX_FLAGS "${INASTEMP_CXX_FLAGS} -Wold-style-cast")
endif()

#===========================================================================
# Generate config
#===========================================================================
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Src/InastempConfig.h.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/Src/InastempConfig.h  )

#===========================================================================
# Build lib
#===========================================================================
add_subdirectory(Src)


if(NOT INASTEMP_JUST_LIB)
    #===========================================================================
    # Build patterns
    #===========================================================================
    add_subdirectory(Patterns)

    #===========================================================================
    # Build unit tests
    #===========================================================================
    enable_testing()
    add_subdirectory(UTests)

    #===========================================================================
    # Build examples
    #===========================================================================
    add_subdirectory(Examples)
endif()

#===========================================================================
# Export cmake info
#===========================================================================
configure_file(${INASTEMP_SOURCE_DIR}/CMakeModules/InastempInfo.cmake.in
    ${INASTEMP_BINARY_DIR}/InastempInfo.cmake
    @ONLY IMMEDIATE
)
install(FILES  ${INASTEMP_BINARY_DIR}/InastempInfo.cmake  DESTINATION lib/  )

if(INASTEMP_AS_SUBPROJECT)
    set(INASTEMP_CXX_FLAGS  "${INASTEMP_CXX_FLAGS}" PARENT_SCOPE)
    set(INASTEMP_BINARY_DIR  "${INASTEMP_BINARY_DIR}" PARENT_SCOPE)
    set(INASTEMP_SOURCE_DIR  "${INASTEMP_SOURCE_DIR}" PARENT_SCOPE)
    set(INASTEMP_INCLUDE_DIR  "${INASTEMP_INCLUDE_DIR}" PARENT_SCOPE)
else()
    # nothing - we are not used by another lib here
endif()



