# CMakeLists.txt for OrbWeaver
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of OrbWeaver.

cmake_minimum_required(VERSION 3.0)

# We have the version here
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION ORB_VERSION)
string(REPLACE "." ";" VERSION_LIST ${ORB_VERSION})
list(GET VERSION_LIST 0 ORB_VERSION_MAJOR)
list(GET VERSION_LIST 1 ORB_VERSION_MINOR)

# We have only C++ sources, but until cmake 3.4 the check_symbol_exists
# macro didn't work properly if only defining C++ as the project language
if( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 4 )
	set(ORB_LANGS C CXX)
else()
	set(ORB_LANGS CXX)
endif()

project(OrbWeaver VERSION ${ORB_VERSION} LANGUAGES ${ORB_LANGS})
set(CMAKE_CXX_STANDARD 11)

set(ORB_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(ORB_DEFINES "")
set(ORB_LIBS "")
set(ORB_CXX_FLAGS "")
set(ORB_LINK_FLAGS "")

option(ORB_HDF5 "Desc" "ON")


macro(find_gsl)
	find_package(GSL REQUIRED)
	list(APPEND VR_INCLUDE_DIRS ${GSL_INCLUDE_DIRS})
	list(APPEND VR_LIBS ${GSL_LIBRARIES})
endmacro()


find_gsl()

macro(find_hdf5)
	if (HDF5_ROOT)
			set(ENV{HDF5_ROOT} ${HDF5_ROOT})
	endif()
	find_package(HDF5 COMPONENTS CXX)
	if (HDF5_FOUND)
		list(APPEND ORB_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
		list(APPEND ORB_LIBS ${HDF5_LIBRARIES})
		list(APPEND ORB_DEFINES USEHDF)
		set(ORB_HAS_HDF5 Yes)
	endif()
endmacro()



set(ORB_HAS_HDF5 No)
if (ORB_HDF5)
	find_hdf5()
endif()


# Include nbodylib, which is definitely needed
add_subdirectory(NBodylib NBodylib)
if (NBODYLIB_VERSION VERSION_LESS "1.10")
	message(FATAL_ERROR "NBodyLib version ${NBODYLIB_VERSION} unsupported, OrbWeaver requires >= 1.10")
endif()

list(APPEND ORB_INCLUDE_DIRS ${NBODYLIB_INCLUDE_DIRS})
list(APPEND ORB_DEFINES ${NBODYLIB_DEFINES})
list(APPEND ORB_LIBS ${NBODYLIB_LIBS})
list(APPEND ORB_CXX_FLAGS ${NBODYLIB_CXX_FLAGS})
list(APPEND ORB_LINK_FLAGS ${NBODYLIB_LINK_FLAGS})


set(ORB_HAS_OPENMP No)
list(FIND ORB_DEFINES "USEOMP" pos)
if(${pos} GREATER -1)
	set(ORB_HAS_OPENMP Yes)
endif()

# cleanup duplicates and convert to single string
list(REMOVE_DUPLICATES ORB_CXX_FLAGS)
string(REPLACE ";" " " ORB_CXX_FLAGS "${ORB_CXX_FLAGS}")
list(REMOVE_DUPLICATES ORB_LINK_FLAGS)
string(REPLACE ";" " " ORB_LINK_FLAGS "${ORB_LINK_FLAGS}")

include_directories(${ORB_INCLUDE_DIRS})


#
# Tell the world what what we are doing
#
macro(orbweaver_report feature)

	# Output feature name and underscore it in the next line
	message("\n${feature}")
	string(REGEX REPLACE "." "-" _underscores ${feature})
	message("${_underscores}\n")

	set(_args "${ARGN}")
	list(LENGTH _args _nargs)
	math(EXPR _nargs "${_nargs} - 1")
	foreach(_idx RANGE 0 ${_nargs} 2)

		# Items in the list come with a message first, then the variable name
		list(GET _args ${_idx} _msg)
		math(EXPR _idx2 "${_idx} + 1")
		list(GET _args ${_idx2} _varname)

		# We try to keep things up to 80 cols
		string(LENGTH ${_msg} _len)
		math(EXPR _nspaces "75 - ${_len}")
		string(RANDOM LENGTH ${_nspaces} _spaces)
		string(REGEX REPLACE "." " " _spaces "${_spaces}")
		string(CONCAT _msg "${_msg}" ${_spaces})
		message(" ${_msg} ${ORB_HAS_${_varname}}")
	endforeach()
endmacro()

message("\nOrbWeaver successfully configured with the following settings:")
orbweaver_report(Dependencies "OpenMP" OPENMP "HDF5" HDF5)


message("")
message("Compilation")
message("-----------")
message("")
message(" Include dirs: ${ORB_INCLUDE_DIRS}")
message(" Macros defined: ${ORB_DEFINES}")
message(" Libs: ${ORB_LIBS}")
message(" C++ flags: ${ORB_CXX_FLAGS}")
message(" Link flags: ${ORB_LINK_FLAGS}")
message("")

add_subdirectory(src)

# Export the include directories, if necessary
# If building on our own, add the "doc" target
if (_export)
        set(ORB_VERSION "${ORB_VERSION_MAJOR}.${ORB_VERSION_MINOR}" PARENT_SCOPE)
        set(ORB_VERSION_MAJOR "${ORB_VERSION_MAJOR}" PARENT_SCOPE)
        set(ORB_VERSION_MINOR "${ORB_VERSION_MINOR}" PARENT_SCOPE)
        set(ORB_INCLUDE_DIRS "${ORB_INCLUDE_DIRS}" PARENT_SCOPE)
        set(ORB_DEFINES "${ORB_DEFINES}" PARENT_SCOPE)
        set(ORB_CXX_FLAGS "${ORB_CXX_FLAGS}" PARENT_SCOPE)
        set(ORB_LINK_FLAGS "${ORB_LINK_FLAGS}" PARENT_SCOPE)
        set(ORB_LIBS "${ORB_LIBS}" PARENT_SCOPE)
        set(ORB_DOC_DIRS "${ORB_DOC_DIRS}" PARENT_SCOPE)
else()
        try_add_doc_target("${ORB_DOC_DIRS}")
endif()
