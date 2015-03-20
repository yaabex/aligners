cmake_minimum_required(VERSION 2.8.4)
project(dsa CXX)

enable_testing(true)

include_directories(. alignment)

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
    alignment)

set(TEST_SRC_DIR tests)

foreach (TEST_FILE ${TEST_FILES})
  add_executable(${TEST_FILE} ${TEST_SRC_DIR}/${TEST_FILE}.cpp ${INCLUDE_FILES})
  target_link_libraries(${TEST_FILE} gtest gtest_main)
  add_test(${TEST_FILE} ${TEST_FILE})
endforeach()

