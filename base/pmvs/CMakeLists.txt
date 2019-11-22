file(GLOB _SRCS "*.[hc]pp")
add_library(pmvs_lib ${_SRCS})

target_include_directories(pmvs_lib PUBLIC
  ${CMVS_BASE_INCLUDE_DIR}
)

target_include_directories(pmvs_lib PRIVATE
  ${CIMG_INCLUDE_DIR}
  ${NLOPT_INCLUDE_DIR}
)