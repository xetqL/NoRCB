# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.15)

project( norcb )
set(CMAKE_CXX_STANDARD 17)

# CGAL and its components

if("${CMAKE_BUILD_TYPE}" STREQUAL Sanitize)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
  message("Flags for build type \"${CMAKE_BUILD_TYPE}\" are: ${CMAKE_CXX_FLAGS}")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL Profile)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -p ${CMAKE_CXX_FLAGS_DEBUG} -frounding-math")
  message("Flags for build type \"${CMAKE_BUILD_TYPE}\" are: ${CMAKE_CXX_FLAGS}")
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
file(GLOB_RECURSE HEADERS include/*)
list(REMOVE_ITEM SRC src/main.cpp)

add_library(norcb STATIC ${HEADERS})
set_target_properties(norcb PROPERTIES LINKER_LANGUAGE CXX)

if(MPI_COMPILE_FLAGS)
  set_target_properties(norcb PROPERTIES
          COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(norcb PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

  find_package( CGAL QUIET COMPONENTS  )
  if ( NOT CGAL_FOUND )
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()
  endif()
  target_link_libraries(norcb PUBLIC CGAL::CGAL)
  target_compile_definitions(norcb PUBLIC USE_CGAL)

target_link_libraries(norcb PUBLIC ${MPI_LIBRARIES})
