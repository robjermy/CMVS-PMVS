# Download and unpack eigen at configure time
# Copy file with external library definition to download folder
configure_file(cmake/eigen.cmake eigen-download/CMakeLists.txt)

# Generate CMake files that will be used to perform the download
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
RESULT_VARIABLE result
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/eigen-download)

if(result)
message(FATAL_ERROR "CMake step for eigen failed: ${result}")
endif()

# Download eigen from git here, checkout branch
execute_process(COMMAND ${CMAKE_COMMAND} --build .
RESULT_VARIABLE result
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/eigen-download )

if(result)
message(FATAL_ERROR "Build step for eigen failed: ${result}")
endif()

# Add eigen directly to our build.
add_subdirectory(
  ${CMAKE_SOURCE_DIR}/lib/eigen
  ${CMAKE_SOURCE_DIR}/lib/eigen/build
  EXCLUDE_FROM_ALL
)

# Download and unpack boost at configure time
# Copy file with external library definition to download folder
configure_file(cmake/boost.cmake boost-download/CMakeLists.txt)

# Generate CMake files that will be used to perform the download
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
RESULT_VARIABLE result
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/boost-download)

if(result)
message(FATAL_ERROR "CMake step for boost failed: ${result}")
endif()

# Download boost from git here, checkout branch
execute_process(COMMAND ${CMAKE_COMMAND} --build .
RESULT_VARIABLE result
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/boost-download )

if(result)
message(FATAL_ERROR "Build step for boost failed: ${result}")
endif()
