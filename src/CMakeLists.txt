set(SOURCES lib.cpp)
set(HEADERS MotionVec.h ForceVec.h PTransform.h RBInertia.h
						ABInertia.h EigenTypedef.h EigenUtility.h Operators.h
						MathFunc.h
						SpaceVecAlg)

include_directories(${EIGEN_INCLUDE_DIRS})

# !!! SpaceVecAlg is an empty library !!!
# !!! Its aim is to test header file and have an IDE integration !!!
# !!! Don't link against it !!!
add_library(SpaceVecAlg ${SOURCES} ${HEADERS})

set(INSTALL_PATH include/SpaceVecAlg)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})

