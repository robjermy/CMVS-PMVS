add_library(image_lib
  camera.cpp
  image.cpp
  photo.cpp
  photoSetS.cpp
)

set_property(TARGET image_lib PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(image_lib jpeg-static)
target_include_directories(image_lib PRIVATE ${CMAKE_SOURCE_DIR}/lib/libjpeg-turbo)