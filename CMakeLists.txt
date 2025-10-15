cmake_minimum_required(VERSION 3.16)
project(Volume_Ascent_Edge_Algorithm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Put all artifacts in the build dir for predictable paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Third-party headers (vendored/submodules)
set(EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/cpp/third_party/eigen)
set(NANOFLANN_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/cpp/third_party/nanoflann/include)

include_directories(${EIGEN3_INCLUDE_DIR} ${NANOFLANN_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/cpp/include)

# Common compile flags
add_compile_definitions(NDEBUG)
add_compile_options(-O3 -ffast-math)

# ===== Libraries =====
add_library(ppp_ball SHARED cpp/ppp_ball.cpp)
target_include_directories(ppp_ball PUBLIC ${EIGEN3_INCLUDE_DIR})
target_compile_definitions(ppp_ball PRIVATE BUILD_LIB)
set_target_properties(ppp_ball PROPERTIES OUTPUT_NAME "ppp_ball")

add_library(initial_proj SHARED cpp/initial_proj.cpp)
target_include_directories(initial_proj PUBLIC ${EIGEN3_INCLUDE_DIR})
target_compile_definitions(initial_proj PRIVATE BUILD_LIB)
set_target_properties(initial_proj PROPERTIES OUTPUT_NAME "initial_proj")

# Keep output name "VA_edge" to match legacy loaders; env aliases map to wrappers
add_library(VA_edge SHARED cpp/VA_edge.cpp)
target_include_directories(VA_edge PUBLIC ${EIGEN3_INCLUDE_DIR} ${NANOFLANN_INCLUDE_DIR})
target_compile_definitions(VA_edge PRIVATE BUILD_LIB)
set_target_properties(VA_edge PROPERTIES OUTPUT_NAME "VA_edge")

# ===== Executables (optional) =====
add_executable(ppp_ball_cli cpp/ppp_ball.cpp)
target_compile_definitions(ppp_ball_cli PRIVATE BUILD_EXECUTABLE)

add_executable(initial_proj_cli cpp/initial_proj.cpp)
target_compile_definitions(initial_proj_cli PRIVATE BUILD_EXECUTABLE)

add_executable(VA_edge_cli cpp/VA_edge.cpp)
target_compile_definitions(VA_edge_cli PRIVATE BUILD_EXECUTABLE)
set_target_properties(VA_edge_cli PROPERTIES OUTPUT_NAME "VA_edge_cli")

# ===== Compute absolute lib paths for wrappers (cross‑platform) =====
# Use CMake's shared lib prefix/suffix to build the expected filenames.
function(_compute_lib_path OUT_VAR TARGET_NAME)
  get_target_property(_out ${TARGET_NAME} OUTPUT_NAME)
  if(NOT _out)
    set(_out ${TARGET_NAME})
  endif()
  set(_fname "${CMAKE_SHARED_LIBRARY_PREFIX}${_out}${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(_path "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${_fname}")
  set(${OUT_VAR} "${_path}" PARENT_SCOPE)
endfunction()

_compute_lib_path(_PPP_BALL_PATH ppp_ball)
_compute_lib_path(_PPP_PROJ_PATH initial_proj)
_compute_lib_path(_PPP_VA_EDGE_PATH VA_edge)

# Cache variables you can read from CMake or IDEs
set(PPP_BALL_LIB     "${_PPP_BALL_PATH}"     CACHE STRING "Path to ppp_ball shared library")
set(PPP_PROJ_LIB     "${_PPP_PROJ_PATH}"     CACHE STRING "Path to initial_proj shared library")
set(PPP_VA_EDGE_LIB  "${_PPP_VA_EDGE_PATH}"  CACHE STRING "Path to VA_edge shared library")
# Back‑compat alias for the Python wrapper that uses PPP_VAEDGE_LIB (no extra underscore)
set(PPP_VAEDGE_LIB   "${_PPP_VA_EDGE_PATH}"  CACHE STRING "Alias env var path to VA_edge shared library")

# Write helper scripts to export env vars after build
set(EXPORT_SH   ${CMAKE_BINARY_DIR}/export_lib_paths.sh)
set(EXPORT_PS1  ${CMAKE_BINARY_DIR}/export_lib_paths.ps1)
set(EXPORT_BAT  ${CMAKE_BINARY_DIR}/export_lib_paths.bat)

file(WRITE ${EXPORT_SH}  "export PPP_BALL_LIB=\"${PPP_BALL_LIB}\"
export PPP_PROJ_LIB=\"${PPP_PROJ_LIB}\"
export PPP_VA_EDGE_LIB=\"${PPP_VA_EDGE_LIB}\"
export PPP_VAEDGE_LIB=\"${PPP_VA_EDGE_LIB}\"
")
file(WRITE ${EXPORT_PS1} "$Env:PPP_BALL_LIB='${PPP_BALL_LIB}'
$Env:PPP_PROJ_LIB='${PPP_PROJ_LIB}'
$Env:PPP_VA_EDGE_LIB='${PPP_VA_EDGE_LIB}'
$Env:PPP_VAEDGE_LIB='${PPP_VA_EDGE_LIB}'
")
file(WRITE ${EXPORT_BAT} "@echo off
set PPP_BALL_LIB=${PPP_BALL_LIB}
set PPP_PROJ_LIB=${PPP_PROJ_LIB}
set PPP_VA_EDGE_LIB=${PPP_VA_EDGE_LIB}
set PPP_VAEDGE_LIB=${PPP_VA_EDGE_LIB}
")

add_custom_target(export_lib_paths ALL
  DEPENDS ppp_ball initial_proj VA_edge
  COMMENT "Writing env export scripts for libraries"
)

# Friendly configure messages
message(STATUS "PPP_BALL_LIB = ${PPP_BALL_LIB}")
message(STATUS "PPP_PROJ_LIB = ${PPP_PROJ_LIB}")
message(STATUS "PPP_VA_EDGE_LIB = ${PPP_VA_EDGE_LIB}")
message(STATUS "PPP_VAEDGE_LIB = ${PPP_VAEDGE_LIB}")
message(STATUS "Wrote: ${EXPORT_SH}; ${EXPORT_PS1}; ${EXPORT_BAT}")

# ===== Install (optional) =====
install(TARGETS ppp_ball initial_proj VA_edge
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib)
install(TARGETS ppp_ball_cli initial_proj_cli VA_edge_cli
        RUNTIME DESTINATION bin)
install(DIRECTORY python/VA DESTINATION .
        PATTERN "__pycache__" EXCLUDE)
