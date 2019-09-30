###########################################################################
# Inastemp - Berenger Bramas MPCDF - 2016
# Under MIT Licence, please you must read the LICENCE file.
###########################################################################
# This goes with the getARMCpuInfos.cpp
# This will create one CMAKE value per output option from the cpp file.
# For example the output of the CPP file can be:
# SSE3=TRUE;AVX=FALSE
# Then it will create:
# CPU_INFO_SSE3 = TRUE
# CPU_INFO_AVX = FALSE
#
# The binary should return 0 on success.
###########################################################################################
macro(GetARMCpuInfos)
# The original CPP file
set(GetARMCpuInfosFile "${PROJECT_SOURCE_DIR}/CMakeModules/getARMCpuInfos.cpp")

# Fatal error if the file does not exist
if(NOT EXISTS ${GetARMCpuInfosFile})
	message(FATAL_ERROR "The GetARMCpuInfosFile does not exist (${GetARMCpuInfosFile})")
endif()

OPTION( INASTEMP_ARMIE_CPU  "Set to ON to run the CPU detection over armie" OFF )

# Compile and execute the file
if(INASTEMP_ARMIE_CPU)
    SET( INASTEMP_ARMIE_CPU_ARGS "" CACHE STRING "Arguments for armie"  )

    if($ENV{VERBOSE})
    	message(STATUS "GetARMCpuInfosFile -- use armie")
    	message(STATUS "GetARMCpuInfosFile -- INASTEMP_ARMIE_CPU_ARGS = ${INASTEMP_ARMIE_CPU_ARGS}")
    endif()
    
    get_filename_component(
		GetARMCpuInfosFileExec ${GetARMCpuInfosFile}
		NAME_WE
	)

    set(GetARMCpuInfosFileExec "${CMAKE_CURRENT_BINARY_DIR}/${GetARMCpuInfosFileExec}")

    try_compile(COMPILE_RESULT_VAR
            ${CMAKE_CURRENT_BINARY_DIR}/GetARMCpuInfos
            ${GetARMCpuInfosFile}
            OUTPUT_VARIABLE comp
            COPY_FILE ${GetARMCpuInfosFileExec})

    if(COMPILE_RESULT_VAR)
        if(NOT EXISTS ${GetARMCpuInfosFileExec})
	        message(FATAL_ERROR "The GetARMCpuInfosFile compiled file does not exist (${GetARMCpuInfosFileExec})")
        endif()

        exec_program("armie ${INASTEMP_ARMIE_CPU_ARGS} -- ${GetARMCpuInfosFileExec}" ${CMAKE_CURRENT_BINARY_DIR}
                OUTPUT_VARIABLE run
                RETURN_VALUE RUN_RESULT_VAR)
    endif() 
else()
    # Simply try and run
    try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
          ${CMAKE_CURRENT_BINARY_DIR} ${GetARMCpuInfosFile}
          COMPILE_OUTPUT_VARIABLE comp
          RUN_OUTPUT_VARIABLE run)
endif()

# If it has successfuly compiled an run
if(COMPILE_RESULT_VAR AND (RUN_RESULT_VAR EQUAL 0) )
	set( CPU_OPTIONS ${run} )
	# For each value
	foreach(optionNode ${run})
		# Get name and value
		string(REPLACE "=" ";" optionNameAndValue ${optionNode})
		list(LENGTH optionNameAndValue optionLength)
		# If we get both
		if(optionLength EQUAL 2)
			list(GET optionNameAndValue 0 optionName)
			list(GET optionNameAndValue 1 optionValue)
			# create cmake variable
			set(CPU_INFO_${optionName} ${optionValue})
		else()
			message(WARNING "GetARMCpuInfosFile wrong format for ${optionNode}.")
		endif()
	endforeach()
	# output the sentence from the binrary
    if($ENV{VERBOSE})
    	message(STATUS "GetARMCpuInfosFile -- results : ${CPU_OPTIONS}")
    endif()
else()
	message(WARNING "GetARMCpuInfosFile -- did not return correctly.")
	message(WARNING "GetARMCpuInfosFile -- compilation output : ${comp}.")
	message(WARNING "GetARMCpuInfosFile -- execution output : ${run}.")
endif()

endmacro(GetARMCpuInfos)
