add_subdirectory(numeric)
add_subdirectory(image)
add_subdirectory(pmvs)
add_subdirectory(cmvs)

add_executable(pmvs3 pmvs.cpp)
target_link_libraries(pmvs3 PRIVATE
  pmvs_lib
  image_lib
  numeric_lib
  nlopt
  Threads::Threads
  nlohmann_json::nlohmann_json
)

target_include_directories(pmvs3 PUBLIC
  ${JSON_INCLUDE_DIR}
)

add_executable(genOption genOption.cpp)

add_executable(cmvs3 cmvs.cpp)

target_link_libraries(cmvs3 PRIVATE
  cmvs_lib
  image_lib
  graclus
  Threads::Threads
)

target_include_directories(cmvs3 PUBLIC
  ${CMVS_INCLUDE_DIR}
)

install(TARGETS genOption DESTINATION bin/)
install(TARGETS cmvs3 DESTINATION bin/)
install(TARGETS pmvs3 DESTINATION bin/)
