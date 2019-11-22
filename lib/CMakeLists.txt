include(ExternalProject)
find_package(Git REQUIRED)
include(cmvsExtensions)

ExternalGitProject(libjpeg-turbo https://github.com/libjpeg-turbo/libjpeg-turbo.git 2.0.3)
add_subdirectory(
  ${libjpeg-turbo_SOURCE_DIR}
  EXCLUDE_FROM_ALL
)

# Header only
ExternalGitProject(cimg https://github.com/dtschump/CImg.git v.2.7.5)

# Eigen lib with modified CMakeLists to not build uninstall target
ExternalGitProject(eigen https://github.com/eigenteam/eigen-git-mirror.git 3.3.7)
file(COPY ${CMAKE_SOURCE_DIR}/lib/eigen-custom/CMakeLists.txt DESTINATION ${eigen_SOURCE_DIR})
add_subdirectory(
  ${eigen_SOURCE_DIR}
  EXCLUDE_FROM_ALL
)

ExternalGitProject(nlopt https://github.com/stevengj/nlopt.git v2.6.1)
add_subdirectory(
  ${nlopt_SOURCE_DIR}
  EXCLUDE_FROM_ALL
)

ExternalGitProject(graclus https://github.com/robjermy/Graclus.git 1.3)
add_subdirectory(
  ${graclus_SOURCE_DIR}
  EXCLUDE_FROM_ALL
)
