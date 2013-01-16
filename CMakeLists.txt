cmake_minimum_required(VERSION 2.8) 

project(math.jsx)
set(CMAKE_BUILD_TYPE debug)
file(GLOB SRC "src/*.cpp")

set(JSX ${PROJECT_SOURCE_DIR}/../jsx)
set(EXTERN ${JSX}/extern)

set(EXECUTABLE_OUTPUT_PATH ${JSX}/bin)
set(LIBRARY_OUTPUT_PATH ${JSX}/bin)
set(CMAKE_CURRENT_BINARY_DIR ${JSX}/bin)

set(CMAKE_SHARED_LINKER_FLAG "${CMAKE_SHARED_LINKER_FLAGS} lt_sys_lib_dlsearch_path_spec=${EXECUTABLE_OUTPUT_PATH}")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(CMAKE_CXX_FLAGS "-Wno-parentheses-equality -DNDEBUG")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set(CMAKE_CXX_FLAGS "-std=c++0x -fPIC -DNDEBUG")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

include_directories(${EXTERN})
include_directories(${JSX}/sdk/core/src)

link_directories(${JSX}/bin)
add_library(math.jsx SHARED ${SRC})
target_link_libraries(math.jsx core.jsx)
