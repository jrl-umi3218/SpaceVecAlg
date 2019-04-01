#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

CONFIGURE_FILE(setup.in.py "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")

SET(PYTHON "python")
SET(PIP "pip")
IF(${PYTHON_BINDING_FORCE_PYTHON3})
  SET(PYTHON "python3")
  SET(PIP "pip3")
ENDIF()

# Build the bindings locally at build time for test purposes
ADD_CUSTOM_TARGET(sva-python-bindings ALL
  COMMAND ${PYTHON} setup.py build_ext --inplace
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Generating local SpaceVecAlg Python bindings"
  DEPENDS include/sva_wrapper.hpp sva/c_sva.pxd sva/sva.pxd sva/sva.pyx
  SOURCES include/sva_wrapper.hpp sva/c_sva.pxd sva/sva.pxd sva/sva.pyx
)

IF(NOT ${DISABLE_TESTS})
  ADD_CUSTOM_TARGET(test-sva-python-bindings
    COMMAND nosetests
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Test SpaceVecAlg Python bindings"
    DEPENDS tests/test_sva_containers.py tests/test_sva_inertia.py tests/test_sva_pickle.py tests/test_sva_ptransform.py tests/test_sva_vector.py
    SOURCES tests/test_sva_containers.py tests/test_sva_inertia.py tests/test_sva_pickle.py tests/test_sva_ptransform.py tests/test_sva_vector.py
  )
  ADD_TEST(NAME PythonBindingsTest
    COMMAND nosetests
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endif()

SET(PIP_EXTRA_OPTIONS "")
IF(${PYTHON_BINDING_USER_INSTALL})
  SET(PIP_EXTRA_OPTIONS "--user")
ENDIF()
# Install the bindings
IF(DEFINED PYTHON_DEB_ROOT)
  ADD_CUSTOM_TARGET(install-sva-python-bindings
    COMMAND python setup.py install --root=${PYTHON_DEB_ROOT} --install-layout=deb && python3 setup.py install --root=${PYTHON_DEB_ROOT} --install-layout=deb
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Install SpaceVecAlg Python bindings (Debian package)"
  )
ELSE()
  ADD_CUSTOM_TARGET(install-sva-python-bindings
    COMMAND ${PIP} install ${PIP_EXTRA_OPTIONS} .
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Install SpaceVecAlg Python bindings"
  )
ENDIF()
INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install-sva-python-bindings)")
