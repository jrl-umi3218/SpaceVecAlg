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
set(HEADERS SpaceVecAlg/MotionVec.h
            SpaceVecAlg/ForceVec.h
            SpaceVecAlg/ImpedanceVec.h
            SpaceVecAlg/AdmittanceVec.h
            SpaceVecAlg/PTransform.h
            SpaceVecAlg/RBInertia.h
            SpaceVecAlg/ABInertia.h
            SpaceVecAlg/EigenTypedef.h
            SpaceVecAlg/EigenUtility.h
            SpaceVecAlg/Operators.h
            SpaceVecAlg/MathFunc.h
            SpaceVecAlg/Conversions.h
            SpaceVecAlg/SpaceVecAlg)



# !!! SpaceVecAlg is an empty library !!!
# !!! Its aim is to test header file and have an IDE integration !!!
# !!! Don't link against it !!!
add_library(SpaceVecAlg ${SOURCES} ${HEADERS})

if(DEFINED Eigen_INCLUDE_DIR)
  target_include_directories(SpaceVecAlg SYSTEM INTERFACE "${Eigen_INCLUDE_DIR}")
else()
  target_include_directories(SpaceVecAlg SYSTEM INTERFACE "${_Eigen_INCLUDE_DIRS}")
endif()

# Targets:
#   * <prefix>/lib/libbar.a
#   * <prefix>/lib/libbaz.a
#   * header location after install: <prefix>/include/foo/Bar.hpp
#   * headers can be included by C++ code `#include <foo/Bar.hpp>`
install(
    TARGETS SpaceVecAlg
    EXPORT "${TARGETS_EXPORT_NAME}"
    LIBRARY DESTINATION "lib"
    ARCHIVE DESTINATION "lib"
    RUNTIME DESTINATION "bin"
    INCLUDES DESTINATION "${INCLUDE_INSTALL_DIR}"
)
install(FILES ${HEADERS} DESTINATION ${INCLUDE_INSTALL_DESTINATION})

