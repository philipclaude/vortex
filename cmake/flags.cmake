# normal warnings
set(NORMAL_WARNINGS -Wall -Wextra)

#
# switch over compilers
# supported: gnu, clang, intel
#
if(CMAKE_COMPILER_IS_GNUCXX)
  #
  # === gnu compiler flags
  #
  message(STATUS "setting gnu compiler flags")

  # warnings
  set(GNU_WARNING_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wunused-result -Winit-self -Wno-variadic-macros -Wno-vla -Wno-strict-overflow")
  set(FULL_WARNINGS ${GNU_WARNING_FLAGS})

  # basic flags
  set(CMAKE_CXX_FLAGS "${GNU_WARNING_FLAGS} -std=c++11 -fopenmp -fstrict-aliasing -Wstrict-aliasing -pedantic -Wnon-virtual-dtor" CACHE STRING "C++ Flags" FORCE)
	if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	  # stack size for predicates (per thread) is too small, short term solution is to use discontiguous stack but this is only supported for g++ on Linux
		#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsplit-stack")
	endif()
  set(CMAKE_C_FLAGS "${GNU_WARNING_FLAGS} -fopenmp -fstrict-aliasing -Wstrict-aliasing" CACHE STRING "C flags" FORCE)
  if(NOT CYGWIN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" CACHE STRING "C++ flags" FORCE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC" CACHE STRING "C flags" FORCE)
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -U__STRICT_ANSI__ -Wa,-mbig-obj -Og" CACHE STRING "C++ flags" FORCE)
  endif()

  # debug flags
  set(CMAKE_CXX_FLAGS_DEBUG "-g -ftrapv -fbounds-check" CACHE STRING "C++ debug flags" FORCE)
  if(NOT CYGWIN)
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0" CACHE STRING "C++ debug flags" FORCE)
  else()
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og" CACHE STRING "C++ debug flags" FORCE)
  endif()
  set(CMAKE_C_FLAGS_DEBUG "-g -O0 -ftrapv -fbounds-check" CACHE STRING "C debug flags" FORCE)

  # release flags
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops" CACHE STRING "C++ release flags" FORCE)
  set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG  -funroll-loops" CACHE STRING "C release flags" FORCE)

  # Look for the gold linker and use it if it's available. Needed for sanitizer on ubuntu
  find_program(GOLD gold)
  if(GOLD)
    set(GOLD_FLAG "-fuse-ld=gold")
  endif()
  unset(GOLD CACHE)

  # memcheck flags
  set(CMAKE_CXX_FLAGS_MEMCHECK "-g -Os -fsanitize=address -fno-omit-frame-pointer" CACHE STRING "C++ memcheck flags" FORCE)
  set(CMAKE_C_FLAGS_MEMCHECK "-g -Os -fsanitize=address -fno-omit-frame-pointer" CACHE STRING "C memcheck flags" FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS_MEMCHECK "-fsanitize=address ${GOLD_FLAG}" CACHE STRING "Flags used by the linker during the creation of .so's." FORCE)

  # coverage flags
  set(GNU_NO_INLINE_FLAGS "-DALWAYS_INLINE=inline -fno-inline -fno-inline-functions -fno-inline-small-functions -fno-inline-functions-called-once -fno-default-inline -fno-implicit-inline-templates")
  set(CMAKE_CXX_FLAGS_COVERAGE "-g -O0 --coverage -ftest-coverage ${GNU_NO_INLINE_FLAGS}" CACHE STRING "C++ coverage flags" FORCE)

  # linker flags
  set(CMAKE_SHARED_LINKER_FLAGS "" CACHE STRING "Flags used by the linker during the creation of .so's." FORCE)
  set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "--coverage" CACHE STRING "Executable linker flags for coverage testing" FORCE)

  # profiler flags
  if(USE_PROFILER)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pg")
  endif()

elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
  #
  # === intel compiler flags
  #
  message(STATUS "setting intel compiler flags")

  # warnings
  set(INTEL_WARNING_FLAGS "-wd3415")
  set(FULL_WARNINGS ${INTEL_WARNING_FLAGS})

  # basic flags
  set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -fPIC -fstrict-aliasing -ansi-alias-check ${INTEL_WARINNGS}" CACHE STRING "C++ flags" FORCE)
  #set(CMAKE_CXX_FLAGS "-Wall -fPIC -fstrict-aliasing -ansi-alias-check ${INTEL_WARINNGS}" CACHE STRING "C++ flags" FORCE)
  set(CMAKE_C_FLAGS "-Wall -fPIC -fstrict-aliasing -ansi-alias-check" CACHE STRING "C flags" FORCE)

  # debug flags
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "C++ debug flags" FORCE)
  set(CMAKE_C_FLAGS_DEBUG "-g -O0" CACHE STRING "C debug flags" FORCE)

  # release flags
  set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "C++ release flags" FORCE)
  set(CMAKE_C_FLAGS_RELEASE "-O3" CACHE STRING "C release flags" FORCE)

  # --- no coverage with intel ---

  # linker flags
  #set(CMAKE_SHARED_LINKER_FLAGS "-stdlib=libc++" CACHE STRING "Flags used by the linker during the creation of .so's." FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS "" CACHE STRING "Flags used by the linker during the creation of .so's." FORCE)
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
  #
  #  === clang compiler flags ===
  #
  message(STATUS "setting clang compiler flags")

  # warnings
  set(CLANG_WARNING_FLAGS "-Wall -Wstrict-aliasing -Wnon-virtual-dtor -pedantic -Wno-variadic-macros -Wno-vla-extension")
  if(CLANG_VERSION VERSION_GREATER 3.4 AND NOT APPLE)
    set(CLANG_WARNING_FLAGS "${CLANG_WARNING_FLAGS} -Wno-absolute-value")
  endif()
  set(FULL_WARNINGS ${CLANG_WARNING_FLAGS})

  # basic flags
  set(CMAKE_CXX_FLAGS "${CLANG_WARNING_FLAGS} -fstrict-aliasing" CACHE STRING "C++ flags" FORCE)
  set(CMAKE_C_FLAGS "-Wall -fstrict-aliasing -Wstrict-aliasing" CACHE STRING "C flags" FORCE)
  if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" CACHE STRING "C++ flags" FORCE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC" CACHE STRING "C flags" FORCE)
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" CACHE STRING "C++ flags" FORCE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC" CACHE STRING "C flags" FORCE)
  endif()

  # debug flags
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fsanitize=address -ftrapv -fno-omit-frame-pointer" CACHE STRING "C++ debug flags" FORCE)
  set(CMAKE_C_FLAGS_DEBUG "-g -O0 -ftrapv -fno-omit-frame-pointer" CACHE STRING "C debug flags" FORCE)

  # release flags
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops" CACHE STRING "C++ release flags" FORCE)
  set(CMAKE_C_FLAGS_RELEASE "-O3 -funroll-loops" CACHE STRING "C release flags" FORCE)

  # memcheck flags
  set(CMAKE_CXX_FLAGS_MEMCHECK "-g -O0 -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls" CACHE STRING "C++ debug flags" FORCE)

  # coverage flags
  set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_DEBUG} --coverage" CACHE STRING "C++ coverage flags" FORCE)
  set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_DEBUG} --coverage" CACHE STRING "C++ coverage flags" FORCE)

  # linker flags
  if(NOT CYGWIN)
    set(CMAKE_EXE_LINKER_FLAGS "-lc++" CACHE STRING "Executable linker flags." FORCE)
  endif()
  set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "--coverage" CACHE STRING "Executable linker flags for coverage testing" FORCE)
  if(APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-undefined,error" CACHE STRING "Flags used by the linker during the creation of .dylib's" FORCE)
  else()
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined" CACHE STRING "Flags used by the linker during the creation of .so's" FORCE)
  endif()
else()
  #message(FATAL_ERROR "unknown compiler: " ${CMAKE_CXX_COMPILER_ID})
endif()
