set(extern_repositories)
function(add_extern_repository name)
	set(options "")
	set(one_value_args GIT_REPOSITORY FULL_HISTORY SKIP_CONFIG)
	set(multi_value_args "")
	cmake_parse_arguments(ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})
	if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/${name})	
		if(ARG_GIT_REPOSITORY)
			set(clone_opts --recursive --depth=1 --single-branch)
			if(ARG_FULL_HISTORY)
				set(clone_opts --recursive)
			endif()
			set(fetch_get git)
			set(fetch_url ${ARG_GIT_REPOSITORY})
			set(fetch_arg clone ${clone_opts} ${ARG_GIT_REPOSITORY} ${CMAKE_CURRENT_SOURCE_DIR}/extern/${name})
		else()
			message(FATAL_ERROR "unknown repository type")
		endif()
		message(STATUS "fetching ${name} from ${fetch_url}")
		execute_process(COMMAND ${fetch_get} ${fetch_arg} RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
	endif()
	if(NOT ARG_SKIP_CONFIG)
		add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/${name} EXCLUDE_FROM_ALL)
	endif()
	list(APPEND extern_repositories "${CMAKE_CURRENT_SOURCE_DIR}/extern/${name}")
	set(extern_repositories ${extern_repositories} PARENT_SCOPE)
endfunction()

set(WITH_GMF_AIO TRUE)
set(WITH_GMF_FORTRAN FALSE)
set(WINGS_BUILD_APPS FALSE)
set(ABSL_PROPAGATE_CXX_STD ON)
set(ABSL_USE_SYSTEM_INCLUDES ON)

add_extern_repository(fmt GIT_REPOSITORY "https://github.com/fmtlib/fmt")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/fmt/include)
add_extern_repository(libmeshb GIT_REPOSITORY "https://github.com/LoicMarechal/libMeshb" SKIP_CONFIG TRUE)
add_extern_repository(stb GIT_REPOSITORY "https://github.com/nothings/stb" SKIP_CONFIG TRUE)
add_extern_repository(tinyobjloader GIT_REPOSITORY "https://github.com/tinyobjloader/tinyobjloader")
add_extern_repository(argparse GIT_REPOSITORY "https://github.com/p-ranav/argparse")
add_extern_repository(morton GIT_REPOSITORY "https://github.com/morton-nd/morton-nd")
add_extern_repository(wings GIT_REPOSITORY "https://github.com/middleburygcl/wings" SKIP_CONFIG TRUE)
add_extern_repository(OpenNL GIT_REPOSITORY "https://github.com/middleburygcl/geogram.psm.OpenNL" SKIP_CONFIG TRUE)
add_extern_repository(PCK GIT_REPOSITORY "https://github.com/middleburygcl/geogram.psm.Predicates" SKIP_CONFIG TRUE)
add_extern_repository(trees GIT_REPOSITORY "https://github.com/middleburygcl/trees.git" SKIP_CONFIG TRUE)
add_extern_repository(stlext GIT_REPOSITORY "https://github.com/middleburygcl/stlext.git" SKIP_CONFIG TRUE)
add_extern_repository(abseil GIT_REPOSITORY "https://github.com/abseil/abseil-cpp")
add_extern_repository(nlopt GIT_REPOSITORY "https://github.com/stevengj/nlopt")
add_extern_repository(json GIT_REPOSITORY "https://github.com/nlohmann/json")

# utilities to clean up and update repositories
add_custom_target(vortex_clean_extern COMMAND rm -rf ${extern_repositories})

set(WINGS_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/extern/wings/wings.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/extern/wings/util/log.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/extern/wings/util/shader.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/extern/wings/util/util.cpp
)
add_library(vortex_wings ${WINGS_SOURCES})
target_compile_definitions(vortex_wings PRIVATE WINGS_COMPILE_STB)

# external repositories
set(external_libraries fmt argparse vortex_wings nlopt)
set(external_libraries ${external_libraries} absl::hash absl::container_memory absl::flat_hash_set absl::memory)

# OpenGL
set(GL_LIBRARIES)
if (APPLE)
        find_library(OpenGL_LIBRARY OpenGL)
        set(GL_LIBRARIES ${LIBRARIES} ${OpenGL_LIBRARY})
else()
        find_package(OpenGL COMPONENTS REQUIRED OpenGL EGL)
        set(GL_LIBRARIES ${LIBRARIES} OpenGL::EGL pthread)
endif()
target_link_libraries(vortex_wings ${GL_LIBRARIES} fmt)

# TBB
find_package(TBB)
if(TBB_FOUND)
	message(STATUS "found TBB")
	set(external_libraries ${external_libraries} TBB::tbb)
	add_definitions(-DVORTEX_WITH_TBB=1)
else()
	message(WARNING "did not find TBB")
	add_definitions(-DVORTEX_WITH_TBB=0)
endif()

# OpenMP
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
	message(STATUS "found OpenMP")
	set(external_libraries ${external_libraries} OpenMP::OpenMP_CXX)
	add_definitions(-DVORTEX_WITH_OMP=1)
else()
	message(STATUS "did not find OpenMP")
	add_definitions(-DVORTEX_WITH_OMP=0)
endif()

set(VORTEX_EXTERNAL_LIBRARIES ${external_libraries} ${GL_LIBRARIES})

# set all include directories
set(VORTEX_INCLUDE_DIRS
	${CMAKE_CURRENT_SOURCE_DIR}/extern/abseil
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/libmeshb/sources
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/OpenNL/OpenNL_psm
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/PCK
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/stb
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/wings
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/morton/include
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/shapelib
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/trees
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/stlext
	${CMAKE_CURRENT_SOURCE_DIR}/extern/json/include
)


set(EXTERN_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/extern/OpenNL/OpenNL_psm/OpenNL_psm.c
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/PCK/Predicates_psm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/extern/tinyobjloader/tiny_obj_loader.cc
	${CMAKE_CURRENT_SOURCE_DIR}/extern/libmeshb/sources/libmeshb7.c
)
add_library(vortex_external OBJECT ${EXTERN_SOURCES})

set_target_properties(vortex_external PROPERTIES COMPILE_FLAGS "-w")
