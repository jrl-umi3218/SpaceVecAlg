# Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
#
# This file is part of SpaceVecAlg.
#
# SpaceVecAlg is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SpaceVecAlg is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SpaceVecAlg.  If not, see <http://www.gnu.org/licenses/>.

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

