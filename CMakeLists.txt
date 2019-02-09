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
set(ORB_LIBS "")
set(ORB_CXX_FLAGS "")
set(ORB_LINK_FLAGS "")


macro(find_gsl)
	find_package(GSL REQUIRED)
	list(APPEND ORB_INCLUDE_DIRS ${GSL_INCLUDE_DIRS})
	list(APPEND ORB_LIBS ${GSL_LIBRARIES})
endmacro()


find_gsl()

macro(find_hdf5)
	if (HDF5_ROOT)
			set(ENV{HDF5_ROOT} ${HDF5_ROOT})
	endif()
	find_package(HDF5 COMPONENTS CXX REQUIRED)
	list(APPEND ORB_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
	list(APPEND ORB_LIBS ${HDF5_LIBRARIES})
endmacro()


find_hdf5()

# cleanup duplicates and convert to single string
list(REMOVE_DUPLICATES ORB_CXX_FLAGS)
string(REPLACE ";" " " ORB_CXX_FLAGS "${ORB_CXX_FLAGS}")
list(REMOVE_DUPLICATES ORB_LINK_FLAGS)
string(REPLACE ";" " " ORB_LINK_FLAGS "${ORB_LINK_FLAGS}")

include_directories(${ORB_INCLUDE_DIRS})


#
# Tell the world what what we are doing
#

message("\nOrbWeaver successfully configured with the following settings:")

message("")
message("Compilation")
message("-----------")
message("")
message(" Include dirs: ${ORB_INCLUDE_DIRS}")
message(" Libs: ${ORB_LIBS}")
message(" C++ flags: ${ORB_CXX_FLAGS}")
message(" Link flags: ${ORB_LINK_FLAGS}")
message("")

add_subdirectory(src)
