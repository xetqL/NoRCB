# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( norcb )

set(CMAKE_CXX_STANDARD 17)
# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
set(CMAKE_BUILD_TYPE Debug)
if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory
include_directories(include)
include_directories(useful-algorithms/include)

file(GLOB SRC src/*)
# include for local package

message(${SRC})
# Creating entries for target: norcb
# ############################

#if (${CMAKE_BUILD_TYPE} EQUAL "Debug")
set(CMAKE_CXX_FLAGS "-g -O0")
#endif ()

add_executable( norcb ${SRC})

add_to_cached_list( CGAL_EXECUTABLE_TARGETS norcb )

# Link the executable to CGAL and third-party libraries
target_link_libraries(norcb PRIVATE CGAL::CGAL )

