# Add eigen directly to our build.
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/lib/eigen
  ${CMAKE_SOURCE_DIR}/lib/eigen-build
  EXCLUDE_FROM_ALL
)

# Add nlopt directly to our build.
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/lib/nlopt
  ${CMAKE_SOURCE_DIR}/lib/nlopt-build
  EXCLUDE_FROM_ALL
)
