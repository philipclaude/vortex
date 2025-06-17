
# set main flags for memcheck
set(MEMCHECK_FLAGS --tool=memcheck --leak-check=yes --show-leak-kinds=all --errors-for-leak-kinds=all --num-callers=50 --error-exitcode=1 --gen-suppressions=all)

# add the suppressions list to the memcheck flags
#list(APPEND MEMCHECK_FLAGS --suppressions=${CMAKE_SOURCE_DIR}/cmake/suppressions/osx)
#list(APPEND MEMCHECK_FLAGS --suppressions=${CMAKE_SOURCE_DIR}/cmake/suppressions/ubuntu)

# find gcov and valgrind
find_program(COVERAGE_COMMAND NAMES gcov)
find_program(MEMORYCHECK_COMMAND NAMES valgrind)

# cmake options for memcheck flags
set(VALGRIND_MEMCHECK_FLAGS CACHE STRING "additonal valgrind memcheck flags: Usefull ones are --track-origins=yes")

# set the main memcheck command
set(MEMCHECK_COMMAND ${MEMORYCHECK_COMMAND} ${MEMCHECK_FLAGS} ${VALGRIND_MEMCHECK_FLAGS})

if(BUILD_TYPE AND BUILD_TYPE MATCHES "DEBUG")
  set(STACKCHECK_COMMAND ${MEMORYCHECK_COMMAND} --tool=exp-sgcheck --gen-suppressions=all)
  set(STACKCHECK_COMMAND ${STACKCHECK_COMMAND} --suppressions=${CMAKE_SOURCE_DIR}/cmake/suppressions/stackcheck)
else()
  set(STACKCHECK_COMMAND ${CMAKE_COMMAND} -E echo "Please set CMAKE_BUILD_TYPE to \\'debug\\' for stack checking")
endif()

# set callgrind and massif commands
set(CALLGRIND_COMMAND ${MEMORYCHECK_COMMAND} --tool=callgrind)
set(MASSIF_COMMAND ${MEMORYCHECK_COMMAND} --tool=massif)

#===============================================================================
# GDB
#===============================================================================
if(APPLE)
  find_program(GDB_COMMAND NAMES lldb gdb)
else()
  find_program(GDB_COMMAND NAMES gdb)
endif()

if(NOT GDB_COMMAND)
  message(STATUS "Please install gdb (lldb for Apple) for debugging")
  set(USE_GDB OFF)
else()
  option(USE_GDB "debugger" ON)
endif()

#===============================================================================
# GUI for CallGrind
#===============================================================================
if(APPLE)
  set(KCACHEGRINDEXEC qcachegrind)
else()
  set(KCACHEGRINDEXEC kcachegrind)
endif()

find_program(KCACHEGRIND ${KCACHEGRINDEXEC})
if(NOT KCACHEGRIND)
  message(STATUS "Please install ${KCACHEGRINDEXEC} for graphical profiling.")
  set(USE_KCACHEGRIND OFF)
else()
  option(USE_KCACHEGRIND "Use ${KCACHEGRINDEXEC} with valgrinds callgrind." ON)
endif()

#===============================================================================
# GUI for massif
#===============================================================================
set(MASSIFGUIEXEC massif-visualizer)

find_program(MASSIFGUI ${MASSIFGUIEXEC})
if(NOT MASSIFGUI)
  message(STATUS "Please install ${MASSIFGUIEXEC} for graphical heap memory usage.")
  set(USE_MASSIFGUI OFF)
else()
  option(USE_MASSIFGUI "Use ${MASSIFGUIEXEC} with valgrinds massif." ON)
endif()

if(BIN_DIR_NAME MATCHES "DEPLOY")
  set(VALGRIND_INIT OFF)
else()
  set(VALGRIND_INIT ON)
endif()

option(VALGRIND_CALLGRIND "enable make targets for valgrinds callgrind tool." OFF)
option(VALGRIND_MEMCHECK "enable make targets for valgrinds memcheck tool." ${VALGRIND_INIT})
option(VALGRIND_STACKCHECK "enable make targets for valgrinds stackcheck tool." OFF)
option(VALGRIND_MASSIF "enable make targets for valgrinds massif heap memory useage tool." OFF)

#===============================================================================
# function for adding debugging targets
#===============================================================================
function(ADD_DEBUG_TARGETS UNIT_TEST WORKDIR)

  if(VALGRIND_MEMCHECK)
    # add the memory checking target: assumes ${UNIT_TEST}_exe exists
    add_custom_target(${UNIT_TEST}_memcheck COMMAND ${MEMCHECK_COMMAND} $<TARGET_FILE:${UNIT_TEST}_exe> $(UNITARGS)
                      WORKING_DIRECTORY ${WORKDIR})
  endif()

  if(VALGRIND_STACKCHECK)
    # add the stack checking target
    add_custom_target(${UNIT_TEST}_stackcheck COMMAND ${STACKCHECK_COMMAND} $<TARGET_FILE:${UNIT_TEST}_exe> $(UNITARGS)
                       WORKING_DIRECTORY ${WORKDIR})
  endif()

  set(CALLGRIND_FILE $<TARGET_FILE_DIR:${UNIT_TEST}_exe>/callgrind.out)
  set(CALLGRIND_MESSAGE "Callgrind file: ${CALLGRIND_FILE}")

  if(VALGRIND_CALLGRIND)
    if(KCACHEGRIND AND USE_KCACHEGRIND)
      # add the callgrind target with kcachegrind gui
      add_custom_target(${UNIT_TEST}_callgrind
                         COMMAND ${CALLGRIND_COMMAND} --callgrind-out-file=${CALLGRIND_FILE} $<TARGET_FILE:${UNIT_TEST}_exe>
                         COMMAND ${KCACHEGRIND} ${CALLGRIND_FILE}
                         COMMAND ${CMAKE_COMMAND} -E echo
                         COMMAND ${CMAKE_COMMAND} -E echo ${CALLGRIND_MESSAGE}
                         COMMAND ${CMAKE_COMMAND} -E echo
                         WORKING_DIRECTORY ${WORKDIR})
    else()
      # add the callgrind target
      add_custom_target(${UNIT_TEST}_callgrind
                         COMMAND ${CALLGRIND_COMMAND} --callgrind-out-file=${CALLGRIND_FILE} $<TARGET_FILE:${UNIT_TEST}_exe>
                         COMMAND ${CMAKE_COMMAND} -E echo
                         COMMAND ${CMAKE_COMMAND} -E echo ${CALLGRIND_MESSAGE}
                         COMMAND ${CMAKE_COMMAND} -E echo
                         WORKING_DIRECTORY ${WORKDIR})
    endif()
  endif()

  set(MASSIF_FILE $<TARGET_FILE_DIR:${UNIT_TEST}_exe>/massif.out)
  set(MASSIF_MESSAGE "Massif file: ${MASSIF_FILE}")

  if(VALGRIND_MASSIF)
    if(MASSIFGUI AND USE_MASSIFGUI)
      # add the massif target with MASSIF_GUI
      add_custom_target(${UNIT_TEST}_massif
                         COMMAND ${MASSIF_COMMAND} --massif-out-file=${MASSIF_FILE} $<TARGET_FILE:${UNIT_TEST}_exe>
                         COMMAND ${MASSIFGUI} ${MASSIF_FILE}
                         COMMAND ${CMAKE_COMMAND} -E echo
                         COMMAND ${CMAKE_COMMAND} -E echo ${MASSIF_MESSAGE}
                         COMMAND ${CMAKE_COMMAND} -E echo
                         WORKING_DIRECTORY ${WORKDIR})
    else()
      # add the massif target
      add_custom_target(${UNIT_TEST}_massif
                         COMMAND ${MASSIF_COMMAND} --massif-out-file=${MASSIF_FILE} $<TARGET_FILE:${UNIT_TEST}_exe>
                         COMMAND ${CMAKE_COMMAND} -E echo
                         COMMAND ${CMAKE_COMMAND} -E echo ${MASSIF_MESSAGE}
                         COMMAND ${CMAKE_COMMAND} -E echo
                         WORKING_DIRECTORY ${WORKDIR})
    endif()
  endif()

  if(USE_GDB)
    add_custom_target(${UNIT_TEST}_gdb
      COMMAND ${GDB_COMMAND} $<TARGET_FILE:${UNIT_TEST}_exe>
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data)
  endif()
endfunction()
