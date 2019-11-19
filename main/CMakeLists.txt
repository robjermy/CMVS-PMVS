
ADD_EXECUTABLE(pmvs2 pmvs2.cpp)

TARGET_LINK_LIBRARIES(pmvs2 pmvs_lib image_lib numeric_lib ${PMVS_LIBRARIES} nlopt tinycthread)

ADD_EXECUTABLE(genOption genOption.cpp)

ADD_EXECUTABLE(cmvs cmvs.cpp)
TARGET_LINK_LIBRARIES(cmvs cmvs_lib image_lib graclus_lib ${PMVS_LIBRARIES} tinycthread)

INSTALL(TARGETS genOption DESTINATION bin/)
INSTALL(TARGETS cmvs DESTINATION bin/)
INSTALL(TARGETS pmvs2 DESTINATION bin/)
