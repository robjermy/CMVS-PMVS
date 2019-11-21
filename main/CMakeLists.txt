
add_executable(pmvs2 pmvs2.cpp)
target_link_libraries(pmvs2 pmvs_lib image_lib numeric_lib ${PMVS_LIBRARIES} nlopt tinycthread)

add_executable(genOption genOption.cpp)

add_executable(cmvs cmvs.cpp)
target_link_libraries(cmvs cmvs_lib image_lib graclus ${PMVS_LIBRARIES} tinycthread)

# cimg
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/cimg/CImg.h)

install(TARGETS genOption DESTINATION bin/)
install(TARGETS cmvs DESTINATION bin/)
install(TARGETS pmvs2 DESTINATION bin/)
