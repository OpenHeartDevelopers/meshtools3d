cmake_minimum_required(VERSION 3.14)
project(meshtools3d LANGUAGES CXX)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Relocatable-install rpath: binary looks for bundled libs in ../lib next to itself.
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)
if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
  set(CMAKE_INSTALL_RPATH "@loader_path/../lib")
else()
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()

find_package(CGAL REQUIRED COMPONENTS ImageIO)
find_package(Boost REQUIRED)
find_package(TBB QUIET)

set(MESHTOOLS3D_INCLUDES "./m3d/include")

add_subdirectory(m3d/src)
include_directories(${MESHTOOLS3D_INCLUDES})

add_executable(meshtools3d applications/main.cpp)
target_link_libraries(meshtools3d PRIVATE m3d CGAL::CGAL CGAL::CGAL_ImageIO)

add_executable(laplace_solver applications/laplace_solver.cpp)
target_link_libraries(laplace_solver PRIVATE m3d)

if(TBB_FOUND)
  target_link_libraries(meshtools3d PRIVATE TBB::tbb TBB::tbbmalloc)
  target_compile_definitions(meshtools3d PRIVATE CGAL_LINKED_WITH_TBB)
  target_link_libraries(laplace_solver PRIVATE TBB::tbb TBB::tbbmalloc)
endif()

install(TARGETS meshtools3d laplace_solver
        RUNTIME DESTINATION bin)
