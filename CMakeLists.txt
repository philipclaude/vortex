# setup the project
cmake_minimum_required(VERSION 3.5)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
project(vortex)

# setup libraries and build type
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "-std=c++17")
set_property(GLOBAL PROPERTY TARGET_MESSAGES OFF)

include(cmake/dependencies.cmake)
include(cmake/setup.cmake)
include(cmake/flags.cmake)
include(cmake/debug.cmake)
include(cmake/coverage.cmake)

# setup include directories
include_directories(src extern)
include_directories(SYSTEM ${VORTEX_INCLUDE_DIRS})

# build the libraries, tests and executables
add_subdirectory(src)
