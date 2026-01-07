#
# Copyright 2012-2025 CNRS-UM LIRMM, CNRS-AIST JRL
#

set(HEADERS
    src/SpaceVecAlg/fwd.h
    src/SpaceVecAlg/MotionVec.h
    src/SpaceVecAlg/ForceVec.h
    src/SpaceVecAlg/ImpedanceVec.h
    src/SpaceVecAlg/AdmittanceVec.h
    src/SpaceVecAlg/PTransform.h
    src/SpaceVecAlg/RBInertia.h
    src/SpaceVecAlg/ABInertia.h
    src/SpaceVecAlg/EigenTypedef.h
    src/SpaceVecAlg/EigenUtility.h
    src/SpaceVecAlg/Operators.h
    src/SpaceVecAlg/MathFunc.h
    src/SpaceVecAlg/Conversions.h
    src/SpaceVecAlg/SpaceVecAlg)

add_library(SpaceVecAlg INTERFACE)
jrl_target_enforce_msvc_conformance(SpaceVecAlg PRIVATE)
target_link_libraries(SpaceVecAlg INTERFACE Eigen3::Eigen)
target_compile_features(SpaceVecAlg INTERFACE cxx_auto_type cxx_constexpr)
if(MSVC)
  target_compile_options(SpaceVecAlg INTERFACE /source-charset:utf-8)
endif()
target_include_directories(
  SpaceVecAlg INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:include>)

jrl_target_headers(SpaceVecAlg PUBLIC HEADERS ${HEADERS})
set_target_properties(SpaceVecAlg PROPERTIES VERSION ${PROJECT_VERSION})

# Declare export components
jrl_add_export_component(NAME SpaceVecAlg TARGETS SpaceVecAlg)
