
add_executable(pmvs2 pmvs2.cpp)

target_link_libraries(pmvs2 pmvs_lib image_lib numeric_lib ${PMVS_LIBRARIES} nlopt tinycthread)

add_executable(genOption genOption.cpp)

add_executable(cmvs cmvs.cpp)
target_link_libraries(cmvs cmvs_lib image_lib graclus_lib ${PMVS_LIBRARIES} tinycthread)

install(TARGETS genOption DESTINATION bin/)
install(TARGETS cmvs DESTINATION bin/)
install(TARGETS pmvs2 DESTINATION bin/)
