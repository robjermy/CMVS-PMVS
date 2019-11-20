# CMake powered build system for CMVS
# Create by Pierre Moulon
# First draft the 4 june 2010
# Second draft 24 sept 2010
# Third draft (linux) 12 oct 2010
# 28 july 2014 Make VisualStudio Win32/Win64 detection automatic
# 15.01.2015 Added tinycthread, nlopt, CImg, removed BLAS, LAPACK, GSL, lmfit, winpthread (R. Hiestand)

cmake_minimum_required(VERSION 2.8)

project(CMVS-PMVS2 C CXX)

set(CMAKE_CXX_STANDARD 17)

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
	${CMAKE_SOURCE_DIR}/thirdParty/graclus1.2/metisLib
	${CMAKE_SOURCE_DIR}/thirdParty/nlopt-2.4.2/api
	${CMAKE_SOURCE_DIR}/thirdParty/cimg
	${CMAKE_SOURCE_DIR}/thirdParty/tinycthread
)

# Threads
find_package(Threads REQUIRED)
set(PMVS_LIBRARIES ${PMVS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Eigen
set( EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/lib/eigen )
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${EIGEN3_INCLUDE_DIR}/cmake" )
find_package(Eigen3 3.3.7 REQUIRED)
message("Found Eigen ${EIGEN3_VERSION}")
include_directories( ${EIGEN3_INCLUDE_DIR} )

# Boost
set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     OFF)
find_package(Boost)
if(Boost_FOUND)
	message("Using system boost library: ${Boost_INCLUDE_DIRS}")
	include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
	message("No system boost library found, using internal version")
	set(PMVS_USE_INTERNAL_BOOST TRUE)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/assert/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/concept_check/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/config/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/core/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/detail/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/disjoint_sets/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/graph/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/iterator/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/mpl/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/predef/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/preprocessor/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/property_map/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/smart_ptr/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/static_assert/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/throw_exception/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/tuple/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/type_traits/include)
	include_directories(${CMAKE_SOURCE_DIR}/lib/boost/libs/utility/include)
endif(Boost_FOUND)

# JPEG
find_package(JPEG)
if(JPEG_FOUND)
	message("Using system JPEG library: ${JPEG_LIBRARIES}")
	include_directories(${JPEG_INCLUDE_DIR})
	set(PMVS_LIBRARIES ${PMVS_LIBRARIES} ${JPEG_LIBRARIES})
else(JPEG_FOUND)
	message("No system JPEG library found, using internal version")
	set(PMVS_USE_INTERNAL_JPEG TRUE)
	include_directories(${CMAKE_SOURCE_DIR}/thirdParty/jpeg)
	include_directories(${CMAKE_BINARY_DIR}/thirdParty/jpeg)	# For jconfig.h
	set(PMVS_LIBRARIES ${PMVS_LIBRARIES} jpeg)
endif(JPEG_FOUND)

# ZLIB
#find_package(ZLIB REQUIRED)
#include_directories(${ZLIB_INCLUDE_DIR})

# PNG
set(PMVS_USE_PNG        OFF)
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
set(PMVS_USE_TIFF        OFF)
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

