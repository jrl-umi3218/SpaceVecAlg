set(SOURCES Operators.cpp)
set(HEADERS MotionVec.h ForceVec.h PTransform.h RBInertia.h
						ABInertia.h EigenTypedef.h EigenUtility.h Operators.h)

find_path(EIGEN_INCLUDE_DIRS Eigen/Array
                             /usr/include/eigen3
                             /usr/local/include/eigen3
                             /opt/local/include/eigen3)

include_directories(${EIGEN_INCLUDE_DIRS})
add_library(SpaceVecAlg ${SOURCES} ${HEADERS})
