add_library(numeric_lib
  mat2.hpp
  mat3.hpp
  mat4.hpp
  mylapack.hpp
  mylapack.cpp
  vec2.hpp
  vec3.hpp
  vec4.hpp
)

target_include_directories(numeric_lib PRIVATE ${EIGEN_INCLUDE_DIR})