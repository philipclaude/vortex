
if(BUILD_TYPE AND BUILD_TYPE MATCHES "COVERAGE")

  # find lcov
  find_program(LCOV lcov)
  if(NOT LCOV)
    message(FATAL_ERROR "could not find lcov")
  endif()

  # find genhtml
  find_program(GENHTML genhtml)
  if(NOT GENHTML )
    message(FATAL_ERROR "could not find genhtml")
  endif()

  # find gcov
  find_program(GCOV gcov)
  if(NOT GCOV)
    message(FATAL_ERROR "could not find gcov")
  endif()

  # set the 'open' command for displaying coverage in a browser
 if(APPLE)
    set(OPEN open)
  else()
    set(OPEN xdg-open)
  endif()

  # set coverage flags
  set(COVERAGE_INFO coverage.info)
  set(HTMLDIR coverage_html)
  set(LCOV_FLAGS --capture -q --gcov-tool ${GCOV} --no-external --base-directory ${CMAKE_SOURCE_DIR} --directory . --output-file ${COVERAGE_INFO} --ignore mismatch --ignore-errors unused,mismatch,inconsistent)
  set(GENHTML_FLAGS ${COVERAGE_INFO} -q --legend --frames --show-details --demangle-cpp --output-directory ${HTMLDIR} -t "vortex code coverage")
  set(LCOV_REMOVES \"_test.*\" \"${CMAKE_BINARY_DIR}*\" \"${CMAKE_SOURCE_DIR}/extern/*\")

  # branch coverage slows down lcov, especially when looking at predicates
  set(BRANCH_COVERAGE "")

  # target for writing coverage information to output
  add_custom_target(coverage_info COMMAND ${CMAKE_COMMAND} -E echo "Generating tracefile ${COVERAGE_INFO}..."
    COMMAND ${LCOV} ${LCOV_FLAGS}
    COMMAND ${LCOV} --remove ${COVERAGE_INFO} ${LCOV_REMOVES} -q --output-file ${COVERAGE_INFO} ${BRANCH_COVERAGE}
    COMMAND ${LCOV} --summary ${COVERAGE_INFO} ${BRANCH_COVERAGE}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    # macro for adding coverage information for a specific unit test
    macro(ADD_COVERAGE_UT COVERAGE_NAME TARGET_NAME)
      add_custom_target(${COVERAGE_NAME}
                        COMMAND ${CMAKE_COMMAND} -E echo "Generating tracefile ${COVERAGE_INFO}..."
                        COMMAND ${LCOV} ${LCOV_FLAGS}
                        COMMAND ${LCOV} --remove ${COVERAGE_INFO} ${LCOV_REMOVES} -q -o ${COVERAGE_INFO} ${BRANCH_COVERAGE}
                        COMMAND ${CMAKE_COMMAND} -E echo "Generating html documents in ${HTMLDIR}..."
                        COMMAND ${GENHTML}  ${GENHTML_FLAGS}
                        COMMAND ${LCOV} --summary ${COVERAGE_INFO} ${BRANCH_COVERAGE}
                        DEPENDS ${TARGET_NAME}
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    endmacro()

    # coverage information with unit test files
    ADD_COVERAGE_UT(coverage_test "")

    # target for displaying the coverage information in a browser
    add_custom_target(coverage_show
                      COMMAND ${OPEN} ${HTMLDIR}/index.html
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    # target for removing coverage files
    add_custom_target(coverage_clean
                      COMMAND ${CMAKE_COMMAND} -E remove_directory ${HTMLDIR}
                      COMMAND ${CMAKE_COMMAND} -E remove ${COVERAGE_INFO}
                      COMMAND find ${CMAKE_SOURCE_DIR} -name "*.gcda" | xargs rm -f
                      COMMAND find ${CMAKE_SOURCE_DIR} -name "*.gcov" | xargs rm -f
                      COMMAND find ${CMAKE_SOURCE_DIR} -name "*.gcno" | xargs rm -f
                      COMMAND ${CMAKE_COMMAND} -E echo "-- Removed all coverage files"
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
else()
    #
    # === build type is not meant for coverage, define the targets which do nothing ===
    #
    set(NO_COVERAGE_MESSAGE "Please set the CMAKE_BUILD_TYPE to \\'coverage\\' or \\'coverage_release\\' to generate coverage information.")

    macro(ADD_COVERAGE_UT COVERAGE_NAME TARGET_NAME)
      add_custom_target(${COVERAGE_NAME} COMMAND ${CMAKE_COMMAND} -E echo ${NO_COVERAGE_MESSAGE})
    endmacro()

    # override all the targets
    ADD_COVERAGE_UT(coverage_test "" )
    add_custom_target(coverage_info    COMMAND ${CMAKE_COMMAND} -E echo ${NO_COVERAGE_MESSAGE})
    add_custom_target(coverage_show    COMMAND ${CMAKE_COMMAND} -E echo ${NO_COVERAGE_MESSAGE})
    add_custom_target(coverage_clean   COMMAND ${CMAKE_COMMAND} -E echo ${NO_COVERAGE_MESSAGE})
    add_custom_target(coverage_cleaner COMMAND ${CMAKE_COMMAND} -E echo ${NO_COVERAGE_MESSAGE})
endif()
