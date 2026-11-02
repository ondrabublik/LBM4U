cmake_minimum_required(VERSION 3.16)
project(cLBMnet CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release)
endif()

option(LBM_FP32 "Store the distribution functions in single precision" OFF)
option(LBM_BUILD_TESTS "Build the unit tests" ON)

set(CORE_SOURCES
	Collision.cpp
	Deflate.cpp
	Block.cpp
	BlockBasic.cpp
	BlockFreeSurface.cpp
	BlockFreeSurfaceLagrangeParticles.cpp
	BlockFreeSurfaceLagrangeFibers.cpp
	BlockFreeSurfaceEulerParticlesFVM.cpp
	BlockFreeSurfaceEulerParticlesLBM.cpp
	LagrangeParticles.cpp
	LagrangeFibers.cpp
	MarchingCubes.cpp
	Mesh.cpp
	SimulationIO.cpp)

add_library(lbmcore STATIC ${CORE_SOURCES})
target_include_directories(lbmcore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(LBM_FP32)
	target_compile_definitions(lbmcore PUBLIC LBM_FP32)
endif()

if(MSVC)
	# OpenMP 3.1 features (atomic write) need the LLVM runtime of MSVC
	target_compile_options(lbmcore PUBLIC /openmp:llvm /W3 /permissive- /utf-8)
	target_compile_definitions(lbmcore PUBLIC _CRT_SECURE_NO_WARNINGS)
else()
	find_package(OpenMP REQUIRED)
	target_link_libraries(lbmcore PUBLIC OpenMP::OpenMP_CXX)
	target_compile_options(lbmcore PRIVATE -Wall -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable)
	if(MINGW)
		# self contained executables (no libgomp/libstdc++ DLLs needed)
		set(LBM_STATIC_LINK -static)
	endif()
endif()

add_executable(cLBMnet main.cpp)
target_link_libraries(cLBMnet PRIVATE lbmcore ${LBM_STATIC_LINK})

if(LBM_BUILD_TESTS)
	enable_testing()
	add_executable(lbm_tests
		tests/main.cpp
		tests/TestCase.cpp
		tests/test_core.cpp
		tests/test_single_phase.cpp
		tests/test_free_surface.cpp
		tests/test_refinement.cpp)
	target_link_libraries(lbm_tests PRIVATE lbmcore ${LBM_STATIC_LINK})
	add_test(NAME lbm_tests COMMAND lbm_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
