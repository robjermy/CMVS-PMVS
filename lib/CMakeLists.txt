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

# Add graclus directly to our build.
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/lib/graclus
  ${CMAKE_SOURCE_DIR}/lib/graclus-build
  EXCLUDE_FROM_ALL
)
