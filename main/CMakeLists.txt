
add_executable(pmvs2 pmvs2.cpp)
target_link_libraries(pmvs2 PRIVATE
  pmvs_lib
  image_lib
  numeric_lib
  nlopt
  Threads::Threads
)

target_include_directories(pmvs2 PUBLIC
  ${CMVS_BASE_INCLUDE_DIR}
)

add_executable(genOption genOption.cpp)

add_executable(cmvs cmvs.cpp)

target_link_libraries(cmvs PRIVATE
  cmvs_lib
  image_lib
  graclus
  Threads::Threads
)

install(TARGETS genOption DESTINATION bin/)
install(TARGETS cmvs DESTINATION bin/)
install(TARGETS pmvs2 DESTINATION bin/)
