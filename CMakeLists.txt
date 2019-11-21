# CMake powered build system for CMVS
# Create by Pierre Moulon
# First draft the 4 june 2010
# Second draft 24 sept 2010
# Third draft (linux) 12 oct 2010
# 28 july 2014 Make VisualStudio Win32/Win64 detection automatic
# 15.01.2015 Added tinycthread, nlopt, CImg, removed BLAS, LAPACK, GSL, lmfit, winpthread (R. Hiestand)

cmake_minimum_required(VERSION 3.15.0)

project(CMVS-PMVS2 C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# By default build in Release mode
if(NOT CMAKE_BUILD_TYPE AND NOT MSVC)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE AND NOT MSVC)

add_subdirectory(lib)

include_directories(
	${CMAKE_SOURCE_DIR}/program/base/image
	${CMAKE_SOURCE_DIR}/program/base/numeric
	${CMAKE_SOURCE_DIR}/program/base/pmvs
	${CMAKE_SOURCE_DIR}/thirdParty/
	${CMAKE_SOURCE_DIR}/thirdParty/tinycthread
)

# Threads
find_package(Threads REQUIRED)
set(PMVS_LIBRARIES ${PMVS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Eigen
include_directories(${CMAKE_SOURCE_DIR}/lib/eigen)

# nlopt
include_directories(${CMAKE_SOURCE_DIR}/lib/nlopt/src/api)
include_directories(${CMAKE_BINARY_DIR}/lib/nlopt/src/api)

# cimg
include_directories(${CMAKE_SOURCE_DIR}/lib/cimg)

# graclus
include_directories(${CMAKE_SOURCE_DIR}/lib/graclus/metisLib)

# jpeg-turbo
include_directories(${CMAKE_SOURCE_DIR}/lib/libjpeg-turbo)
include_directories(${CMAKE_BINARY_DIR}/lib/libjpeg-turbo)

# ZLIB
#find_package(ZLIB REQUIRED)
#include_directories(${ZLIB_INCLUDE_DIR})

# PNG
set(PMVS_USE_PNG OFF)
if(PMVS_USE_PNG)
  find_package(PNG)
  if(PNG_FOUND)
	  message("PNG library found, enabling support: ${PNG_LIBRARIES}")
	  ADD_DEFINITIONS(-DPMVS_HAVE_PNG)
	  include_directories(${PNG_INCLUDE_DIRS})
	  set(PMVS_LIBRARIES ${PMVS_LIBRARIES} ${PNG_LIBRARIES})
  endif(PNG_FOUND)
endif(PMVS_USE_PNG)

# TIFF
set(PMVS_USE_TIFF OFF)
if(PMVS_USE_TIFF)
  find_package(TIFF)
  if(TIFF_FOUND)
	  message("TIFF library found, enabling support: ${TIFF_LIBRARIES}")
	  ADD_DEFINITIONS(-DPMVS_HAVE_TIFF)
	  include_directories(${TIFF_INCLUDE_DIR})
	  set(PMVS_LIBRARIES ${PMVS_LIBRARIES} ${TIFF_LIBRARIES})
  endif(TIFF_FOUND)
endif(PMVS_USE_TIFF)

add_subdirectory(base)
add_subdirectory(main)
add_subdirectory(thirdParty)
