# add all the source files
set(SOURCES
  array2d.cpp
  elements.cpp
  field.cpp
  graphics.cpp
  halfedges.cpp
  io.cpp
  library.cpp
  log.cpp
  math/linalg.cpp
  math/mat.cpp
  math/spmat.cpp
  math/sym.cpp
  math/vec.cpp
  mesh.cpp
  mesher.cpp
  numerics.cpp
  #predicates.c
  quadrature.cpp
  texture.cpp
  triangulate.cpp
  util.cpp
  voronoi.cpp
  optimaltransport.cpp
  )

# create an object library for all the sources
add_library(vortex_src OBJECT ${SOURCES})
add_definitions(-DVORTEX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_definitions(-DROBUST_PREDICATES=0)

target_compile_features(vortex_src PRIVATE cxx_std_17)
#target_compile_definitions(vortex_src PUBLIC -DVORTEX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_target_properties(vortex_src PROPERTIES
            CXX_STANDARD 17
            CXX_EXTENSIONS OFF)

# create both shared and static libraries with same name
set(VORTEX_LIBS vortex_static vortex_shared )
add_library(vortex_shared SHARED $<TARGET_OBJECTS:vortex_src> $<TARGET_OBJECTS:vortex_external>)
set_target_properties(vortex_shared PROPERTIES OUTPUT_NAME vortex CLEAN_DIRECT_OUTPUT 1)
target_link_libraries(vortex_shared dl)
if (UNIX AND NOT APPLE)
  target_link_libraries(vortex_shared rt dl)
endif()
add_library(vortex_static STATIC $<TARGET_OBJECTS:vortex_src> $<TARGET_OBJECTS:vortex_external>)
set_target_properties(vortex_static PROPERTIES OUTPUT_NAME vortex CLEAN_DIRECT_OUTPUT 1)

add_custom_target(vortex_lib DEPENDS ${VORTEX_LIBS})
target_link_libraries(vortex_shared ${VORTEX_EXTERNAL_LIBRARIES})
target_link_libraries(vortex_static ${VORTEX_EXTERNAL_LIBRARIES})

# create executables
add_executable(vortex vortex.cpp)
add_dependencies(vortex vortex_lib)
target_link_libraries(vortex vortex_static dl)

# function to add unit test files with a suffix in directories
function(add_test_files var type)
  set(sources)
  foreach(dir ${ARGN})
      file(GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${dir}/*_${type}.cpp")
      list(APPEND sources ${_sources})
  endforeach()
  set(${var} ${${var}} ${sources} PARENT_SCOPE)
endfunction()

add_test_files(TEST_FILES "test" ${CMAKE_CURRENT_SOURCE_DIR})
add_test_files(TEST_FILES "test" ${CMAKE_CURRENT_SOURCE_DIR}/math)

# set the directory for the test executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)

# create a target for each test source
foreach(test_src ${TEST_FILES})

  # replace the '/' with '_' and "../" with ""
  string(REPLACE "../" "" test_src_underscored ${test_src})
  string(REPLACE "/" "_" test_src_underscored ${test_src_underscored})

  # remove the extension
  get_filename_component(test_bin ${test_src_underscored} NAME_WE)

  # assign the name of the executable
  set(test_exe ${test_bin}_exe)

  # the actual executable: target for simply building
  add_executable(${test_exe} ${test_src})
  target_link_libraries(${test_exe} vortex_shared ${VORTEX_EXTERNAL_LIBRARIES})
  target_compile_definitions(${test_exe} PUBLIC -DSTANDALONE -DVORTEX_FULL_UNIT_TEST=0)

  # create the target
  add_custom_target(${test_bin} command $<TARGET_FILE:${test_exe}> 1 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data)

  # targets for memory & optionally stack checking & gdb
  ADD_DEBUG_TARGETS(${test_bin} ${CMAKE_SOURCE_DIR}/data)

endforeach()

# create the full unit executable
add_executable(unit_exe tester.cpp ${TEST_FILES})
target_link_libraries(unit_exe vortex_shared)
target_compile_definitions(unit_exe PUBLIC -DSTDOUT_REDIRECT="${CMAKE_BINARY_DIR}/unit_tests_output.txt" -DVORTEX_FULL_UNIT_TEST=1)

# target for running all unit tests from a single executable
add_custom_target(unit command $<TARGET_FILE:unit_exe> WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data)

ADD_DEBUG_TARGETS(unit ${CMAKE_SOURCE_DIR})
ADD_COVERAGE_UT(unit_coverage unit)

# reset the directory for the executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
