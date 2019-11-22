if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	add_definitions(-DNUMBITS=64)
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

file(GLOB _SRCS "*.[hc]pp")
add_library(cmvs_lib ${_SRCS})

target_include_directories(cmvs_lib PUBLIC
  ${CMVS_BASE_INCLUDE_DIR}
)

target_include_directories(cmvs_lib PRIVATE
  ${GRACLUS_INCLUDE_DIR}
)