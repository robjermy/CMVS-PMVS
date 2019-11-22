file(GLOB _SRCS "*.[hc]pp")
add_library(numeric_lib ${_SRCS})

target_include_directories(numeric_lib PRIVATE 
  ${EIGEN_INCLUDE_DIR}
)