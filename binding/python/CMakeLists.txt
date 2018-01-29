# Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
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

MACRO(GET_SVA_PROPERTY PROPERTY)
  GET_TARGET_PROPERTY(SVA_${PROPERTY} SpaceVecAlg ${PROPERTY})
  IF(NOT SVA_${PROPERTY})
    SET(SVA_${PROPERTY} "")
  ENDIF()
ENDMACRO()
GET_SVA_PROPERTY(COMPILE_FLAGS)
GET_SVA_PROPERTY(INCLUDE_DIRECTORIES)
GET_SVA_PROPERTY(LINK_FLAGS)
GET_SVA_PROPERTY(LINK_LIBRARIES)
GET_SVA_PROPERTY(LOCATION)
CONFIGURE_FILE(setup.in.py ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)

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
)

IF(NOT ${DISABLE_TESTS})
  ADD_CUSTOM_TARGET(test-sva-python-bindings
    COMMAND nosetests
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Test SpaceVecAlg Python bindings"
    DEPENDS tests/test_sva_containers.py tests/test_sva_inertia.py tests/test_sva_pickle.py tests/test_sva_ptransform.py tests/test_sva_vector.py
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
INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_MAKE_PROGRAM} install-sva-python-bindings)")
