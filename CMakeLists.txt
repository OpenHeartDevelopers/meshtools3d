# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( meshtools3d )


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS ImageIO  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

find_package( TBB QUIET ) 

if( TBB_FOUND ) 
  include(${TBB_USE_FILE}) 
  list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES}) 
endif() 


# include for local directory

# include for local package


# Creating entries for target: meshtools3d
# ############################

add_executable( meshtools3d  GetPot.cpp main.cpp Mesh.cpp Chrono.cpp ThicknessEvaluation.cpp  LaplaceSolver.cpp  mgmres.cpp VtkWriter.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS meshtools3d )

# Link the executable to CGAL and third-party libraries
target_link_libraries(meshtools3d   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

