cmake_minimum_required(VERSION 3.6)

project(cmvs3 VERSION 3.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Use gold linker to speed up linking time, see cmake/useGoldLinker.cmake
include(useGoldLinker)

# threads
find_package(Threads REQUIRED)

# eigen
set(EIGEN_INCLUDE_DIR 
	${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen
)

# nlopt
set(NLOPT_INCLUDE_DIR
	${CMAKE_CURRENT_SOURCE_DIR}/lib/nlopt/src/api
	${CMAKE_CURRENT_BINARY_DIR}/lib/nlopt/src/api
)

# cimg
set(CIMG_INCLUDE_DIR
	${CMAKE_CURRENT_SOURCE_DIR}/lib/cimg
) 

# graclus
set(GRACLUS_INCLUDE_DIR
	${CMAKE_CURRENT_SOURCE_DIR}/lib/graclus/metisLib
)

# libjpeg-turbo
set(JPEG_INCLUDE_DIR
	${CMAKE_CURRENT_SOURCE_DIR}/lib/libjpeg-turbo
	${CMAKE_CURRENT_BINARY_DIR}/lib/libjpeg-turbo
)

# PNG
set(PMVS_USE_PNG OFF)
if(PMVS_USE_PNG)
  find_package(PNG)
  if(PNG_FOUND)
	  message("PNG library found, enabling support: ${PNG_LIBRARIES}")
	  ADD_DEFINITIONS(-DPMVS_HAVE_PNG)
	  include_directories(${PNG_INCLUDE_DIRS})
	  set(PMVS_LIBRARIES ${PMVS_LIBRARIES} ${PNG_LIBRARIES})
  endif()
endif()

# TIFF
set(PMVS_USE_TIFF OFF)
if(PMVS_USE_TIFF)
  find_package(TIFF)
  if(TIFF_FOUND)
	  message("TIFF library found, enabling support: ${TIFF_LIBRARIES}")
	  ADD_DEFINITIONS(-DPMVS_HAVE_TIFF)
	  include_directories(${TIFF_INCLUDE_DIR})
	  set(PMVS_LIBRARIES ${PMVS_LIBRARIES} ${TIFF_LIBRARIES})
  endif()
endif()

set(CMVS_BASE_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/base
)

# external libraries
add_subdirectory(lib)

# base functionality
add_subdirectory(base)

# main executables
add_subdirectory(main)
