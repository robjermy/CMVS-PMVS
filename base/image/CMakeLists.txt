file(GLOB _SRCS "*.[hc]pp")
add_library(image_lib ${_SRCS})

target_link_libraries(image_lib jpeg-static)

target_include_directories(image_lib PUBLIC
  ${CMVS_BASE_INCLUDE_DIR}
)

target_include_directories(image_lib PRIVATE
  ${CIMG_INCLUDE_DIR}
  ${JPEG_INCLUDE_DIR}
)