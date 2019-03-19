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

set(HEADERS_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include")

set(HEADERS_INCLUDE_DIR
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:${HEADERS_INSTALL_DIR}>)

set(HEADERS
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/MotionVec.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/ForceVec.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/ImpedanceVec.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/AdmittanceVec.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/PTransform.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/RBInertia.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/ABInertia.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/EigenTypedef.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/EigenUtility.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/Operators.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/MathFunc.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/Conversions.h
  ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/SpaceVecAlg
)

add_library(SpaceVecAlg INTERFACE)
if(COMMAND target_sources)
  target_sources(SpaceVecAlg INTERFACE ${HEADERS})
endif()
target_compile_features(SpaceVecAlg INTERFACE cxx_auto_type cxx_constexpr)
target_include_directories(SpaceVecAlg SYSTEM INTERFACE "${Eigen_INCLUDE_DIRS}")
target_include_directories(SpaceVecAlg INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)

install(
    TARGETS SpaceVecAlg
    EXPORT "${TARGETS_EXPORT_NAME}"
)
install(DIRECTORY SpaceVecAlg DESTINATION ${INCLUDE_INSTALL_DIR})

