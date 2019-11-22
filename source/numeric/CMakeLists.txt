file(GLOB _SRCS "*.cpp")
add_library(numeric_lib ${_SRCS})

target_include_directories(numeric_lib PUBLIC 
  ${CMVS_INCLUDE_DIR}
)

target_include_directories(numeric_lib PRIVATE 
  ${EIGEN_INCLUDE_DIR}
)