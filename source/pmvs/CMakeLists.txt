file(GLOB _SRCS "*.cpp")
add_library(pmvs_lib ${_SRCS})

target_include_directories(pmvs_lib PUBLIC
  ${CMVS_INCLUDE_DIR}
  ${JSON_INCLUDE_DIR}
)

target_include_directories(pmvs_lib PRIVATE
  ${CIMG_INCLUDE_DIR}
  ${NLOPT_INCLUDE_DIR}
)

target_link_libraries(pmvs_lib PRIVATE
  nlohmann_json::nlohmann_json
)