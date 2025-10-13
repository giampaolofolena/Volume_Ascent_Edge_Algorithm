cmake_minimum_required(VERSION 3.16)
project(ppp_vaedge LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Third-party headers (vendored/submodules)
set(EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/cpp/third_party/eigen)
set(NANOFLANN_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/cpp/third_party/nanoflann)

include_directories(${EIGEN3_INCLUDE_DIR} ${NANOFLANN_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/cpp/include)

# Common compile flags
add_compile_definitions(NDEBUG)
add_compile_options(-O3 -ffast-math)

# ===== Libraries =====
add_library(ppp_ball SHARED cpp/ppp_ball.cpp)
target_include_directories(ppp_ball PUBLIC ${EIGEN3_INCLUDE_DIR})

add_library(initial_proj SHARED cpp/initial_proj.cpp)
target_include_directories(initial_proj PUBLIC ${EIGEN3_INCLUDE_DIR})

add_library(VAedge SHARED cpp/VAedge.cpp)
target_include_directories(VAedge PUBLIC ${EIGEN3_INCLUDE_DIR} ${NANOFLANN_INCLUDE_DIR})
target_compile_definitions(VAedge PRIVATE BUILD_LIB)

# ===== Executables (optional) =====
add_executable(ppp_ball_cli cpp/ppp_ball.cpp)
target_compile_definitions(ppp_ball_cli PRIVATE BUILD_EXECUTABLE)

add_executable(initial_proj_cli cpp/initial_proj.cpp)
target_compile_definitions(initial_proj_cli PRIVATE BUILD_EXECUTABLE)

add_executable(VAedge_cli cpp/VAedge.cpp)
target_compile_definitions(VAedge_cli PRIVATE BUILD_EXECUTABLE)

# Output names/platform niceties
set_target_properties(ppp_ball PROPERTIES OUTPUT_NAME "ppp_ball")
set_target_properties(initial_proj PROPERTIES OUTPUT_NAME "initial_proj")
set_target_properties(VAedge PROPERTIES OUTPUT_NAME "VAedge")

# Install (optional)
install(TARGETS ppp_ball initial_proj VAedge
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib)
install(TARGETS ppp_ball_cli initial_proj_cli VAedge_cli
        RUNTIME DESTINATION bin)
install(DIRECTORY python/ppp DESTINATION .
        PATTERN "__pycache__" EXCLUDE)

