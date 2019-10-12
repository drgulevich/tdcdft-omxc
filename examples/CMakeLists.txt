# cmake version to be used
cmake_minimum_required(VERSION 2.6)

# project name
project(infqw)

# set relative path policy
cmake_policy(SET CMP0015 NEW)

# flags
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-DARMA_NO_DEBUG -march=native -Ofast -std=c++14") # -pg for gprof profiler
set(CMAKE_CXX_FLAGS_DEBUG "-DARMA_EXTRA_DEBUG -g -Wall -std=c++14")

# display info
message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
if(${CMAKE_BUILD_TYPE} STREQUAL Release)
    message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if(${CMAKE_BUILD_TYPE} STREQUAL Debug)
    message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

include_directories( ../include/ )

link_directories( ../ )

add_executable(infqw infqw.cpp)

target_link_libraries(infqw PRIVATE libtdcdft-omxc.a)
target_link_libraries(infqw PRIVATE armadillo)

add_custom_target(distclean
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt cmake_install.cmake Makefile
#    COMMAND ${CMAKE_COMMAND} -E remove install_manifest.txt *.so.*
    COMMENT "distclean invoked"
    )
