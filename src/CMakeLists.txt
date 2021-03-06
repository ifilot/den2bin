 #*************************************************************************
 #   CMakeLists.txt  --  This file is part of den2bin.                    *
 #                                                                        *
 #   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 #                                                                        *
 #   den2bin is free software: you can redistribute it and/or modify      *
 #   it under the terms of the GNU General Public License as published    *
 #   by the Free Software Foundation, either version 3 of the License,    *
 #   or (at your option) any later version.                               *
 #                                                                        *
 #   den2bin is distributed in the hope that it will be useful,           *
 #   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 #   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 #   See the GNU General Public License for more details.                 *
 #                                                                        *
 #   You should have received a copy of the GNU General Public License    *
 #   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 #                                                                        *
 #*************************************************************************/

# set minimum cmake requirements
cmake_minimum_required(VERSION 2.8)
project (den2bin)

# add custom directory to look for .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules )

# prepare configuration file
SET(VERSION_MAJOR "1")
SET(VERSION_MINOR "3")
SET(VERSION_MICRO "4")
configure_file(config.h.in ../src/config.h @ONLY)

# Enable release build
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# add OS specific
if(APPLE)
    add_definitions(-D_APPLE)
    SET(BOOST_INCLUDEDIR "/usr/local/include")
    SET(BOOST_LIBRARYDIR "/usr/local/lib")
else()
    SET(BOOST_INCLUDEDIR "/usr/include")
    SET(BOOST_LIBRARYDIR "/usr/lib/x86_64-linux-gnu")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    option(HAS_OPENMP "OpenMP enabled" ON)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# set Boost
set (Boost_NO_SYSTEM_PATHS ON)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_STATIC_RUNTIME OFF)
set (BOOST_ALL_DYN_LINK OFF)

# Include libraries
find_package(PkgConfig REQUIRED)
find_package(BZIP2 REQUIRED)
find_package(GLM REQUIRED)
find_package(Boost COMPONENTS regex iostreams filesystem REQUIRED)
pkg_check_modules(TCLAP tclap REQUIRED)

# Set include folders
include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_BINARY_DIR}
                    ${GLM_INCLUDE_DIRS}
                    ${Boost_INCLUDE_DIRS})

# Add sources
file(GLOB_RECURSE SOURCES "*.cpp")
add_executable(den2bin ${SOURCES})

# Set C++14
add_definitions(-std=c++14)

# Link libraries
if(APPLE)
    SET(CMAKE_MACOSX_RPATH TRUE)
    SET_TARGET_PROPERTIES(den2bin PROPERTIES INSTALL_RPATH "@executable_path/lib")
    target_link_libraries(den2bin ${Boost_LIBRARIES} ${BZIP2_LIBRARY})
else()
    target_link_libraries(den2bin ${Boost_LIBRARIES} ${BZIP2_LIBRARY})
endif()

# add Wno-literal-suffix to suppress warning messages
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}")

###
# Installing
##
install (TARGETS den2bin DESTINATION bin)
