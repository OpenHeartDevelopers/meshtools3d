cmake_minimum_required(VERSION 3.21)
project(meshtools3d VERSION 2.0.0 LANGUAGES CXX)

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

# parfile_builder is intentionally stdlib-only — no m3d, no CGAL.
add_executable(parfile_builder applications/parfile_builder.cpp)

if(TBB_FOUND)
  target_link_libraries(meshtools3d PRIVATE TBB::tbb TBB::tbbmalloc)
  target_compile_definitions(meshtools3d PRIVATE CGAL_LINKED_WITH_TBB)
  target_link_libraries(laplace_solver PRIVATE TBB::tbb TBB::tbbmalloc)
endif()

# ── Install + dependency bundling ─────────────────────────────────────────────
# macOS: bundle Homebrew libs (not in /usr/lib or /System/Library) into lib/.
# Windows: bundle vcpkg DLLs alongside the executables in bin/.
# Linux: GMP, MPFR, TBB are system packages — not bundled, documented as deps.
if(APPLE)
  install(TARGETS meshtools3d laplace_solver parfile_builder
          RUNTIME_DEPENDENCY_SET _m3d_deps
          RUNTIME DESTINATION bin)
  install(RUNTIME_DEPENDENCY_SET _m3d_deps
          PRE_EXCLUDE_REGEXES  "^/usr/lib" "^/System/Library"
          POST_EXCLUDE_REGEXES "^/usr/lib" "^/System/Library"
          DESTINATION lib)
  # RUNTIME_DEPENDENCY_SET copies dylibs with their original LC_LOAD_DYLIB entries,
  # which encode absolute Homebrew paths. Rewrite every non-system dep to a
  # @loader_path-relative path so the package is self-contained.
  install(CODE [[
    set(_base "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}")
    file(GLOB _bins "${_base}/bin/meshtools3d" "${_base}/bin/laplace_solver" "${_base}/bin/parfile_builder")
    file(GLOB _libs "${_base}/lib/*.dylib")

    foreach(_file IN LISTS _bins _libs)
      cmake_path(GET _file PARENT_PATH _dir)
      cmake_path(GET _dir FILENAME _dirname)

      execute_process(COMMAND chmod u+w "${_file}")
      execute_process(COMMAND otool -L "${_file}" OUTPUT_VARIABLE _otool_out)
      string(REGEX MATCHALL "\t/[^\t\n (]+" _deps "${_otool_out}")

      foreach(_dep IN LISTS _deps)
        string(STRIP "${_dep}" _dep)
        if(_dep MATCHES "^/usr/lib|^/System/Library")
          continue()
        endif()
        cmake_path(GET _dep FILENAME _dep_name)
        if(NOT EXISTS "${_base}/lib/${_dep_name}")
          continue()
        endif()
        if(_dirname STREQUAL "bin")
          execute_process(COMMAND install_name_tool -change "${_dep}" "@loader_path/../lib/${_dep_name}" "${_file}")
        else()
          execute_process(COMMAND install_name_tool -change "${_dep}" "@loader_path/${_dep_name}" "${_file}")
        endif()
      endforeach()

      if(_dirname STREQUAL "lib" AND NOT IS_SYMLINK "${_file}")
        cmake_path(GET _file FILENAME _file_name)
        execute_process(COMMAND install_name_tool -id "@loader_path/${_file_name}" "${_file}")
      endif()
    endforeach()

    # The install_name_tool edits above rewrite Mach-O load commands, which
    # invalidates each file's code signature. On Apple Silicon, AMFI kills any
    # process that loads an invalidly-signed binary with SIGKILL ("killed: 9")
    # before it prints anything. Re-sign ad-hoc — libraries first, then the
    # executables — so the package launches on a fresh download.
    # (Runs last, after CPACK_STRIP_FILES has already stripped these files.)
    foreach(_file IN LISTS _libs _bins)
      if(IS_SYMLINK "${_file}")
        continue()
      endif()
      execute_process(COMMAND codesign --force --sign - "${_file}")
    endforeach()
  ]])
elseif(WIN32)
  install(TARGETS meshtools3d laplace_solver parfile_builder
          RUNTIME_DEPENDENCY_SET _m3d_deps
          RUNTIME DESTINATION bin)
  install(RUNTIME_DEPENDENCY_SET _m3d_deps
          PRE_EXCLUDE_REGEXES  "api-ms-win-.*" "ext-ms-.*"
          POST_EXCLUDE_REGEXES ".*[Ss]ystem32[/\\\\].*\\.dll"
          DESTINATION bin)
else()
  install(TARGETS meshtools3d laplace_solver parfile_builder
          RUNTIME DESTINATION bin)
endif()

# ── CPack ─────────────────────────────────────────────────────────────────────
string(TOLOWER "${CMAKE_SYSTEM_NAME}"      _cpack_sys)
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _cpack_arch)
string(REPLACE "darwin"  "macos"  _cpack_sys  "${_cpack_sys}")
string(REPLACE "aarch64" "arm64"  _cpack_arch "${_cpack_arch}")
string(REPLACE "amd64"   "x86_64" _cpack_arch "${_cpack_arch}")

set(CPACK_PACKAGE_NAME             "meshtools3d")
# CI overrides this with -DCPACK_PACKAGE_VERSION=<tag> so the tarball filename
# tracks the release tag; local builds fall back to PROJECT_VERSION.
if(NOT DEFINED CPACK_PACKAGE_VERSION)
  set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
endif()
set(CPACK_PACKAGE_VENDOR           "CEMRG")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "3D mesh generation tools for cardiac volumetric models")
set(CPACK_PACKAGE_FILE_NAME
    "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${_cpack_sys}-${_cpack_arch}")
set(CPACK_STRIP_FILES TRUE)
if(WIN32)
  set(CPACK_GENERATOR "ZIP")
else()
  set(CPACK_GENERATOR "TGZ")
endif()
include(CPack)
