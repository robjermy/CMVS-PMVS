
add_executable(pmvs2 pmvs2.cpp)
target_link_libraries(pmvs2 pmvs_lib image_lib numeric_lib  nlopt ${PMVS_LIBRARIES})

add_executable(genOption genOption.cpp)

add_executable(cmvs cmvs.cpp)
target_link_libraries(cmvs cmvs_lib image_lib graclus ${PMVS_LIBRARIES})

install(TARGETS genOption DESTINATION bin/)
install(TARGETS cmvs DESTINATION bin/)
install(TARGETS pmvs2 DESTINATION bin/)
