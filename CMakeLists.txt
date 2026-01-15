cmake_minimum_required(VERSION 3.20)

project(SoarAI
	VERSION 1.0.0
	LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Do not set for deterministic training
option(SOAR_BUILD_FAST "Build SoarAI with more aggressive compiler optimizations for performance boost" OFF)
option(SOAR_BUILD_API "Build C++ binaries for the REST API" OFF)

set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

find_package(CNum 0.2.2 EXACT REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
add_subdirectory(${CMAKE_SOURCE_DIR}/src)
