cmake_minimum_required(VERSION 3.8.0)

project(boost-download NONE)

include(ExternalProject)
ExternalProject_Add(boost
  GIT_REPOSITORY    https://github.com/boostorg/boost.git
  GIT_TAG           afb333b7c5101041f0280b2edf155c55114c9c95
  GIT_SHALLOW       TRUE
  SOURCE_DIR        "${CMAKE_SOURCE_DIR}/lib/boost"
  BINARY_DIR        ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)
