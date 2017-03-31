#CMakeLists.txt to generate SlimCodeML (simplified)
cmake_minimum_required(VERSION 2.6)
project(SlimCodeML)

#SET(CMAKE_VERBOSE_MAKEFILE 1) #For verbose output
SET(CMAKE_BUILD_TYPE RELEASE) #Use extra optimization flags; set to DEBUG if undesired
#Set g++ specific compiler settings
if(CMAKE_COMPILER_IS_GNUCXX) 
	set(CMAKE_CXX_FLAGS_RELEASE "-Wno-unused-result -O3 -ffast-math -mtune=native -minline-stringops-dynamically -funroll-loops -funsafe-loop-optimizations" CACHE "Release mode options" STRING FORCE)
	set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-Wl,-O1" CACHE "Release mode linker options" STRING FORCE)
endif(CMAKE_COMPILER_IS_GNUCXX)

# Set BLAS/LAPACK directories and libraries
if(WIN32) #Untested
	set(BLAS_LINK_DIR ${CMAKE_SOURCE_DIR}/lib CACHE PATH "BLAS lib dir" )
	set(LAPACK_LINK_DIR ${CMAKE_SOURCE_DIR}/lib CACHE PATH "LAPACK lib dir" )
	set(MATH_LIB_NAMES blas_win32_MD.lib lapack_win32_MD.lib CACHE STRING "Math libraries" )
else() #You need libblas.a and liblapack.a in your source directory; otherweise change paths below
	set(BLAS_LINK_DIR ${CMAKE_SOURCE_DIR}  CACHE PATH "BLAS lib dir")
	set(LAPACK_LINK_DIR ${CMAKE_SOURCE_DIR} CACHE PATH "LAPACK lib dir")
	set(MATH_LIB_NAMES blas lapack blas lapack gfortranbegin gfortran CACHE STRING "Math libraries" )
endif()

link_directories( ${BLAS_LINK_DIR} )
link_directories( ${LAPACK_LINK_DIR} )

#Sources (ordered as suggested by gprof -R)
set(SRCS
	tools.cpp
	treesub.cpp
	codeml.cpp
	eigenqrev.cpp
	treespace.cpp
	global.cpp
)

if(WIN32)
	add_definitions(/D_CRT_SECURE_NO_WARNINGS)
	add_definitions(/D_SECURE_SCL=0)
	add_definitions(/Oi)
endif(WIN32)

#Executables
add_executable(slimcodeml ${SRCS})

#Link libraries
target_link_libraries(slimcodeml ${MATH_LIB_NAMES})


