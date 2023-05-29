#
# Copyright 2012-2020 CNRS-UM LIRMM, CNRS-AIST JRL
#

set(HEADERS_INCLUDE_DIR
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>)

set(HEADERS
    ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/fwd.h
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
    ${HEADERS_INCLUDE_DIR}/SpaceVecAlg/SpaceVecAlg)

add_library(SpaceVecAlg INTERFACE)
if(COMMAND target_sources)
  target_sources(SpaceVecAlg INTERFACE ${HEADERS})
endif()
target_compile_features(SpaceVecAlg INTERFACE cxx_auto_type cxx_constexpr)
if(MSVC)
  target_compile_options(SpaceVecAlg INTERFACE /source-charset:utf-8)
endif()
if(TARGET Eigen3::Eigen)
  target_link_libraries(SpaceVecAlg INTERFACE Eigen3::Eigen)
else()
  target_include_directories(SpaceVecAlg SYSTEM
                             INTERFACE "${EIGEN3_INCLUDE_DIR}")
endif()
target_include_directories(
  SpaceVecAlg INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:include>)

install(TARGETS SpaceVecAlg EXPORT "${TARGETS_EXPORT_NAME}")
install(DIRECTORY SpaceVecAlg DESTINATION ${INCLUDE_INSTALL_DIR})
