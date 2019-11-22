if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	add_definitions(-DNUMBITS=64)
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

add_library(cmvs_lib
  bundle.cpp
  graclus.cpp
)

target_include_directories(cmvs_lib PUBLIC
  ${CMVS_BASE_INCLUDE_DIR}
  ${GRACLUS_INCLUDE_DIR}
)