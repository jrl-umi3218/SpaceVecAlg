#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

configure_file(setup.in.py "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")

set(PYTHON "python")
set(PIP "pip")
if(${PYTHON_BINDING_FORCE_PYTHON3})
  set(PYTHON "python3")
  set(PIP "pip3")
endif()

# Build the bindings locally at build time for test purposes
add_custom_target(sva-python-bindings ALL
  COMMAND ${PYTHON} setup.py build_ext --inplace
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Generating local SpaceVecAlg Python bindings"
  DEPENDS include/sva_wrapper.hpp sva/c_sva.pxd sva/sva.pxd sva/sva.pyx
  SOURCES include/sva_wrapper.hpp sva/c_sva.pxd sva/sva.pxd sva/sva.pyx
)

if(NOT ${DISABLE_TESTS})
  add_custom_target(test-sva-python-bindings
    COMMAND nosetests
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Test SpaceVecAlg Python bindings"
    DEPENDS tests/test_sva_containers.py tests/test_sva_inertia.py tests/test_sva_pickle.py tests/test_sva_ptransform.py tests/test_sva_vector.py
    SOURCES tests/test_sva_containers.py tests/test_sva_inertia.py tests/test_sva_pickle.py tests/test_sva_ptransform.py tests/test_sva_vector.py
  )
  add_test(NAME PythonBindingsTest
    COMMAND nosetests
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endif()

set(PIP_EXTRA_OPTIONS "")
if(${PYTHON_BINDING_USER_INSTALL})
  set(PIP_EXTRA_OPTIONS "--user")
endif()
# Install the bindings
if(DEFINED PYTHON_DEB_ROOT)
  add_custom_target(install-sva-python-bindings
    COMMAND python setup.py install --root=${PYTHON_DEB_ROOT} --install-layout=deb && python3 setup.py install --root=${PYTHON_DEB_ROOT} --install-layout=deb
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Install SpaceVecAlg Python bindings (Debian package)"
  )
else()
  add_custom_target(install-sva-python-bindings
    COMMAND ${PIP} install ${PIP_EXTRA_OPTIONS} .
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Install SpaceVecAlg Python bindings"
  )
endif()
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install-sva-python-bindings)")
