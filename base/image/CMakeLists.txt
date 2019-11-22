add_library(image_lib
  camera.cpp
  image.cpp
  photo.cpp
  photoSetS.cpp
)

target_link_libraries(image_lib jpeg-static)
target_include_directories(image_lib PRIVATE
  ${CIMG_INCLUDE_DIR}
  ${JPEG_INCLUDE_DIR}
)

target_include_directories(image_lib PUBLIC
  ${CMVS_BASE_INCLUDE_DIR}
)