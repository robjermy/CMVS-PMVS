cmake_minimum_required(VERSION 3.6)

project(cmvs3 VERSION 3.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Use gold linker to speed up linking time, see cmake/useGoldLinker.cmake
include(useGoldLinker)

add_subdirectory(lib)

include_directories(
	${CMAKE_SOURCE_DIR}/program/base/image
	${CMAKE_SOURCE_DIR}/program/base/numeric
	${CMAKE_SOURCE_DIR}/program/base/pmvs
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
