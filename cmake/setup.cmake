# get the name of the binary directory
STRING(TOUPPER ${CMAKE_BINARY_DIR} BIN_DIR_NAME)
STRING(FIND ${BIN_DIR_NAME} "/" LAST_DIR_IDX REVERSE)
STRING(SUBSTRING ${BIN_DIR_NAME} LAST_DIR_IDX -1 BIN_DIR_NAME)
set(BUILD_TYPE_STRING "Choose the type of build, options are: Debug Release Coverage Memcheck.")

# determine the build type: default is debug
if(NOT CMAKE_BUILD_TYPE)
  if(BIN_DIR_NAME MATCHES "DEBUG")
    set(CMAKE_BUILD_TYPE "debug" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
  elseif(BIN_DIR_NAME MATCHES "COVERAGE" )
    set(CMAKE_BUILD_TYPE "coverage" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
  elseif(BIN_DIR_NAME MATCHES "RELEASE" )
    set(CMAKE_BUILD_TYPE "release" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
  elseif(BIN_DIR_NAME MATCHES "MEMCHECK" )
    set(CMAKE_BUILD_TYPE "memcheck" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
  else()
    set(CMAKE_BUILD_TYPE "release" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
  endif()
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
message(STATUS "build type = " ${CMAKE_BUILD_TYPE})

# setup binary and lib directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# setup include and link directories
include_directories(${CMAKE_SOURCE_DIR}/src)
link_directories(${CMAKE_BINARY_DIR}/lib)

# static versus dynamic builds
if(VORTEX_BUILD_DYNAMIC)
    set(BUILD_SHARED_LIBS TRUE)
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
    add_definitions(-DVORTEX_DYNAMIC_LIBS)
else()
    set(BUILD_SHARED_LIBS FALSE)
endif()
