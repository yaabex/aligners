cmake_minimum_required(VERSION 2.8.4)
project(dsa C CXX)

add_subdirectory(gtest)
enable_testing(true)

include_directories(. ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Default build type
set(CMAKE_BUILD_TYPE Debug)

# Common CXX flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -std=c++14")

# Alignment sources
set(INCLUDE_FILES
    alignment/basic_score.hpp
    alignment/hirschberg.hpp
    alignment/needleman_wunsch.hpp)

# Testing sources
set(TEST_FILES
    tests/needlemanwunsch.cpp)

add_executable(tests ${INCLUDE_FILES} ${TEST_FILES})
target_link_libraries(tests gtest gtest_main)
add_test(tests tests)
