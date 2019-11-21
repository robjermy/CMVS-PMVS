# Add libjpeg-turbo directly to our build.
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/lib/libjpeg-turbo
  EXCLUDE_FROM_ALL
)

# Add nlopt directly to our build.
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/lib/nlopt
  EXCLUDE_FROM_ALL
)

# Add graclus directly to our build.
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/lib/graclus
  EXCLUDE_FROM_ALL
)

# Add eigen directly to our build.
file(COPY ${CMAKE_SOURCE_DIR}/lib/eigen-custom/CMakeLists.txt DESTINATION ${CMAKE_SOURCE_DIR}/lib/eigen)
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/lib/eigen
  EXCLUDE_FROM_ALL
)