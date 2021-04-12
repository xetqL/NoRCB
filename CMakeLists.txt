# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( norcb )
set(CMAKE_CXX_STANDARD 17)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

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

find_package(MPI REQUIRED)

# include for local directory
include_directories(include useful-algorithms/include ${MPI_INCLUDE_PATH})

file(GLOB_RECURSE SRC src/*)
# include for local package

# Creating entries for target: norcb
# ############################

#if (${CMAKE_BUILD_TYPE} EQUAL "Debug")
#set(CMAKE_CXX_FLAGS "-g -O0")
#endif ()

add_executable( norcb ${SRC} include/algorithm.hpp)
add_library(libnorcb ${SRC})

if(MPI_COMPILE_FLAGS)
  set_target_properties(norcb PROPERTIES
          COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(libnorcb PROPERTIES
          COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(norcb PROPERTIES
          LINK_FLAGS "${MPI_LINK_FLAGS}")
  set_target_properties(libnorcb PROPERTIES
          LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

add_to_cached_list( CGAL_EXECUTABLE_TARGETS norcb )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS libnorcb )

# Link the executable to CGAL and third-party libraries
target_link_libraries(norcb PRIVATE CGAL::CGAL ${MPI_LIBRARIES})
target_link_libraries(libnorcb PRIVATE CGAL::CGAL ${MPI_LIBRARIES})

