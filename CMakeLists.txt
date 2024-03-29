###########################################################################
# Inastemp - Berenger Bramas MPCDF - 2016
# Under MIT Licence, please you must read the LICENCE file.
###########################################################################
cmake_minimum_required (VERSION 3.0 FATAL_ERROR)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_POLICY_DEFAULT_CMP0077 OLD)

#===========================================================================
# Project Declaration
#===========================================================================
project(INASTEMP CXX)

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

if(NOT DEFINED INASTEMP_DO_INSTALL)
    set(INASTEMP_DO_INSTALL OFF CACHE BOOL
        "Set to on to install the headers file of inastemp" FORCE)
endif()

if($ENV{VERBOSE})
    MESSAGE(STATUS "Main -- INASTEMP_JUST_LIB ${INASTEMP_JUST_LIB}")
    MESSAGE(STATUS "Main -- INASTEMP_DO_INSTALL ${INASTEMP_DO_INSTALL}")
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

# Available from version 3.1
if(NOT ${CUSTOM_SYSTEM_PROCESSOR} STREQUAL "AURORA")
  if($ENV{VERBOSE})
    MESSAGE(STATUS "Main -- CMAKE_CXX_STANDARD is set to 17")
  endif()
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

if((${CMAKE_MAJOR_VERSION} STREQUAL 3 AND CMAKE_MINOR_VERSION LESS 6) OR (${CMAKE_MAJOR_VERSION} LESS 3) AND (CMAKE_CXX_COMPILER_ID STREQUAL "Intel"))
    MESSAGE(WARNING " Main -- CMAKE_CXX_STANDARD does not work with Intel compiler for cmake version < 3.6, we add -std=c++17 flag by hand which might not be compatible with your project")
    set(INASTEMP_CXX_FLAGS "-std=c++17")
elseif((${CMAKE_MAJOR_VERSION} STREQUAL 3 AND CMAKE_MINOR_VERSION LESS 1) OR (${CMAKE_MAJOR_VERSION} LESS 3))
    MESSAGE(WARNING " Main -- CMAKE_CXX_STANDARD does not work for cmake version < 3.1, we insert -std=c++17 flag by hand which might not be compatible with your project")
    set(INASTEMP_CXX_FLAGS "-std=c++17")
endif()

if(NOT CUSTOM_SYSTEM_PROCESSOR)
    set(CUSTOM_SYSTEM_PROCESSOR ${CMAKE_SYSTEM_PROCESSOR})
endif()

if($ENV{VERBOSE})
    MESSAGE(STATUS "Main -- CMAKE_SYSTEM_PROCESSOR = ${CMAKE_SYSTEM_PROCESSOR}")
    MESSAGE(STATUS "Main -- CUSTOM_SYSTEM_PROCESSOR = ${CUSTOM_SYSTEM_PROCESSOR}")
endif()

if(${CUSTOM_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
    if($ENV{VERBOSE})
        MESSAGE(STATUS "Main -- compile for ppc64le architecture")
    endif()

    #  Ask compiler capacities
    include(GetCompilerInfos)
    GetCompilerInfos()

    # All types from worse to best (ADD-NEW-HERE)
    set(ALL_TYPES "ALTIVEC")
    set(ALTIVEC_CPU_RULES "")

    set(INASTEMP_USE_SCALAR ON)
elseif(${CUSTOM_SYSTEM_PROCESSOR} STREQUAL "ARM" OR ${CUSTOM_SYSTEM_PROCESSOR} STREQUAL "aarch64")
    if($ENV{VERBOSE})
        MESSAGE(STATUS "Main -- compile for arm architecture")
    endif()
    set(CUSTOM_SYSTEM_PROCESSOR "ARM")
    
    #  Ask CPU capacities
    include(GetARMCpuInfos)
    GetARMCpuInfos()

    #  Ask compiler capacities
    include(GetCompilerInfos)
    GetCompilerInfos()

    # All types from worse to best (ADD-NEW-HERE)
    set(ALL_TYPES "SVE")
    set(SVE_CPU_RULES "")

    set(INASTEMP_USE_SCALAR ON)
elseif(${CUSTOM_SYSTEM_PROCESSOR} STREQUAL "AURORA")
    if($ENV{VERBOSE})
        MESSAGE(STATUS "Main -- compile for Nec AURORA architecture")
    endif()
    
    #  Ask compiler capacities
    include(GetCompilerInfos)
    GetCompilerInfos()

    # All types from worse to best (ADD-NEW-HERE)
    set(ALL_TYPES "SXA")
    set(SXA_CPU_RULES "")

    set(INASTEMP_USE_SCALAR ON)
else()
    if($ENV{VERBOSE})
        MESSAGE(STATUS "Main -- compile for x86 architecture")
    endif()
    #  Ask CPU capacities
    include(GetX86CpuInfos)
    GetX86CpuInfos()

    #  Ask compiler capacities
    include(GetCompilerInfos)
    GetCompilerInfos()

    # All types from worse to best (ADD-NEW-HERE)
    set(ALL_TYPES "SSE3;SSSE3;SSE41;SSE42;AVX;AVX2;AVX512COMMON;AVX512KNL;AVX512SKL")
    set(INASTEMP_USE_SCALAR ON)

    # Set custom cpu <=> vec rules (maybe ADD-NEW-HERE if needed)
    set(AVX512COMMON_CPU_RULES "AVX512F")
    set(AVX512KNL_CPU_RULES "AVX512F;AVX512ER;AVX512PF")
    set(AVX512SKL_CPU_RULES "AVX512F;AVX512CD;AVX512VL;AVX512BW;AVX512DQ")


    # Dependencies between types (maybe ADD-NEW-HERE if needed)
    set(SSSE3_DEP "SSE3")
    set(SSE41_DEP "SSSE3")
    set(SSE42_DEP "SSE41")
    set(AVX2_DEP "AVX")
    set(AVX512KNL_DEP "AVX512COMMON")
    set(AVX512SKL_DEP "AVX512COMMON")
    
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(INASTEMP_CXX_FLAGS  "${INASTEMP_CXX_FLAGS} -m64")
    endif()
endif()

# Enforce rules
set(ALL_TYPES_REVERSE ${ALL_TYPES})
list(REVERSE ALL_TYPES_REVERSE)
foreach(TYPE ${ALL_TYPES_REVERSE})   
    if( COMPILER_INFO_${TYPE} )
        if(NOT DEFINED ${TYPE}_CPU_RULES)
            set( ${TYPE}_CPU_RULES "${TYPE}")
        endif()

        # Value by default is based on CPU capacities
        set(${TYPE}_SUPPORTED_BY_CPU TRUE)
        foreach(RULE ${${TYPE}_CPU_RULES})
            if(NOT DEFINED CPU_INFO_${RULE})
                MESSAGE(WARNING "Main -- The variable CPU_INFO_${RULE} cannot be read from the current CPU which might be too old. It turns off ${TYPE}.")
                set(${TYPE}_SUPPORTED_BY_CPU FALSE)                
                break()
            elseif(NOT CPU_INFO_${RULE})
                set(${TYPE}_SUPPORTED_BY_CPU FALSE)                
                break()
            endif()
        endforeach()

        if(NOT DEFINED INASTEMP_USE_${TYPE})
            OPTION( INASTEMP_USE_${TYPE} "Set to ON to compile with ${TYPE} support (default based on CPU capacities)" ${${TYPE}_SUPPORTED_BY_CPU} )
        endif()

        if(INASTEMP_USE_${TYPE})
            set(INASTEMP_CXX_FLAGS "${INASTEMP_CXX_FLAGS} ${${TYPE}_FLAGS}")
            list(APPEND ALL_USED_TYPES "${TYPE}")

            if(COMPILER_INFO_${TYPE}_NOOP)
                set(INASTEMP_USE_${TYPE}_OPERATORS ON)
            else()
                set(INASTEMP_USE_${TYPE}_OPERATORS OFF)
            endif()

            # We turn it on if it is needed by another type
            foreach(DEPTYPE ${${TYPE}_DEP})
                if( NOT COMPILER_INFO_${TYPE} )
                    MESSAGE(ERROR "Main -- Type ${TYPE} depends on type ${DEPTYPE} which is not supported by the compiler.")
                endif()
                if( NOT INASTEMP_USE_${DEPTYPE} )
                    if($ENV{VERBOSE})
                        MESSAGE(STATUS "Main -- ${DEPTYPE} is turn on because it is needed by ${TYPE}")
                    endif()
                    set(INASTEMP_USE_${DEPTYPE} ON)
                endif()
            endforeach() 
        endif()
    else()
        if($ENV{VERBOSE})
            MESSAGE(STATUS "Main -- The compiler cannot compile ${TYPE} intrinsics")
        endif()
        set(INASTEMP_USE_${TYPE} FALSE)
    endif()
endforeach()
# The original order should be used as slow-to-fast types
list(APPEND ALL_USED_TYPES "SCALAR")
list(REVERSE ALL_USED_TYPES)
list(LENGTH ALL_USED_TYPES INASTEMP_NB_USED_TYPES)
if($ENV{VERBOSE})
    MESSAGE(STATUS "Main -- INASTEMP_NB_USED_TYPES ${INASTEMP_NB_USED_TYPES}")
    MESSAGE(STATUS "Main -- ALL_USED_TYPES ${ALL_USED_TYPES}")
endif()
list(GET ALL_USED_TYPES -1 INASTEMP_BESTTYPE )

# For test coverage
option(USE_COVERAGE "Enable profiling and coverage report analysis" OFF)

if( USE_COVERAGE )
    message(STATUS "Use code coverage")
    set(CMAKE_CXX_FLAGS "${INASTEMP_CXX_FLAGS} -g -O0 --coverage")
    list(APPEND INASTEMP_LIBRARIES gcov)
    # Then one can use these commande:
    # lcov --directory CMakeFiles/ -c -o coverage-rapport.info
    # genhtml -o ./coverage -t "Coverage example" coverage-rapport.info
    message(STATUS "USE_COVERAGE has been turned to ${USE_COVERAGE}")
endif()

# Options
if(NOT INASTEMP_AS_SUBPROJECT)
    SET( INASTEMP_EXTRA_CXX_FLAGS_RELEASE "-O3" CACHE STRING "These flags are added to the usual ones for RELEASE"  )
    set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} ${INASTEMP_EXTRA_CXX_FLAGS_RELEASE}")

    SET( INASTEMP_EXTRA_CXX_FLAGS_DEBUG "-O0" CACHE STRING "These flags are added to the usual ones for DEBUG"  )
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} ${INASTEMP_EXTRA_CXX_FLAGS_DEBUG}")

    set(INASTEMP_CXX_FLAGS "${INASTEMP_CXX_FLAGS} -fpic -Wextra -Wnon-virtual-dtor -Wshadow -Wpointer-arith -Wcast-qual -Wconversion  -Wall -Wno-sign-conversion -pedantic -Woverloaded-virtual -Wpointer-arith -Wcast-qual -Wconversion -Wno-error")

    if(${CUSTOM_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
        # Builting functions are not correctly making variables/parameters used
        set(INASTEMP_CXX_FLAGS "${INASTEMP_CXX_FLAGS} -Wno-unused-but-set-variable -Wno-unused-but-set-parameter")
    endif()
    if(${CUSTOM_SYSTEM_PROCESSOR} STREQUAL "ARM")
        # Builting functions are not correctly making variables/parameters used
        set(INASTEMP_CXX_FLAGS "${INASTEMP_CXX_FLAGS} -Wno-vla-extension")
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
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Src/InastempCompileConfig.h.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/Src/InastempCompileConfig.h  )

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
if(INASTEMP_DO_INSTALL)
    install(FILES  ${INASTEMP_BINARY_DIR}/InastempInfo.cmake  DESTINATION lib/  )
endif()

if(INASTEMP_AS_SUBPROJECT)
    set(INASTEMP_CXX_FLAGS  "${INASTEMP_CXX_FLAGS}" PARENT_SCOPE)
    set(INASTEMP_BINARY_DIR  "${INASTEMP_BINARY_DIR}" PARENT_SCOPE)
    set(INASTEMP_SOURCE_DIR  "${INASTEMP_SOURCE_DIR}" PARENT_SCOPE)
    set(INASTEMP_INCLUDE_DIR  "${INASTEMP_INCLUDE_DIR}" PARENT_SCOPE)
else()
    # nothing - we are not used by another lib here
endif()



