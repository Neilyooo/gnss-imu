cmake_minimum_required(VERSION 3.16)

project(gnss_imu_fusion LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Eigen3 REQUIRED)

add_library(gnss_imu_fusion STATIC
  src/benchmark.cpp
  src/data_manager.cpp
  src/dataset_io.cpp
  src/fusion.cpp
  src/internal/estimator_common.cpp
  src/estimators/eskf_estimator.cpp
  src/estimators/ukf_estimator.cpp
  src/estimators/fgo_estimator.cpp
)

target_include_directories(gnss_imu_fusion
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(gnss_imu_fusion
  PUBLIC
    Eigen3::Eigen
)

target_compile_options(gnss_imu_fusion
  PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)

add_executable(gnss_imu_bench
  src/main.cpp
)

target_link_libraries(gnss_imu_bench
  PRIVATE
    gnss_imu_fusion
)
