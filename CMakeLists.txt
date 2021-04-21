# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( norcb )
set(CMAKE_CXX_STANDARD 17)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if("${CMAKE_BUILD_TYPE}" STREQUAL Sanitize)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
  message("Flags for build type \"${CMAKE_BUILD_TYPE}\" are: ${CMAKE_CXX_FLAGS}")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL Profile)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -p ${CMAKE_CXX_FLAGS_DEBUG} -frounding-math")
  message("Flags for build type \"${CMAKE_BUILD_TYPE}\" are: ${CMAKE_CXX_FLAGS}")
endif()

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
list(REMOVE_ITEM SRC src/main.cpp)

# include for local package

# Creating entries for target: norcb
# ############################

#if (${CMAKE_BUILD_TYPE} EQUAL "Debug")
#set(CMAKE_CXX_FLAGS "-g -O0")
#endif ()

add_executable( norcb_test src/main.cpp ${SRC})
add_library(norcb src/norcb.cpp src/geom/utils.cpp)

if(MPI_COMPILE_FLAGS)
  set_target_properties(norcb_test PROPERTIES
          COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(norcb PROPERTIES
          COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(norcb_test PROPERTIES
          LINK_FLAGS "${MPI_LINK_FLAGS}")
  set_target_properties(norcb PROPERTIES
          LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

add_to_cached_list( CGAL_EXECUTABLE_TARGETS norcb_test )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS norcb )

# Link the executable to CGAL and third-party libraries
target_link_libraries(norcb_test PRIVATE CGAL::CGAL ${MPI_LIBRARIES})
target_link_libraries(norcb PRIVATE CGAL::CGAL ${MPI_LIBRARIES})

