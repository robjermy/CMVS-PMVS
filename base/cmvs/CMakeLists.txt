if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	add_definitions(-DNUMBITS=64)
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

add_library(cmvs_lib
  bundle.cpp
  graclus.cpp
)
