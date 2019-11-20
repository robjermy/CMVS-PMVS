cmake_minimum_required(VERSION 3.8.0)

project(eigen-download NONE)

include(ExternalProject)
ExternalProject_Add(eigen
  GIT_REPOSITORY    https://github.com/eigenteam/eigen-git-mirror.git
  GIT_TAG           3.3.7
  GIT_SHALLOW       TRUE
  SOURCE_DIR        "${CMAKE_SOURCE_DIR}/lib/eigen"
  BINARY_DIR        ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)
