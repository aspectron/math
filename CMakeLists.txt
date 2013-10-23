cmake_minimum_required(VERSION 2.8) 

project(math.jsx)

set(JSX_ROOT $ENV{JSX_ROOT})
if(NOT DEFINED JSX_ROOT)
	set(JSX_ROOT ~/dev/jsx)
	message("Environment variable JSX_ROOT not defined, using ${JSX_ROOT}")
endif(NOT DEFINED JSX_ROOT)

include(${JSX_ROOT}/CMakeLists.txt.global)

file(GLOB SRC "src/*.cpp")

add_library(math.jsx SHARED ${SRC})
target_link_libraries(math.jsx core.jsx)
