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

INCLUDE(../../cmake/python.cmake)

# Look for Python 2.7
SET(Python_ADDITIONAL_VERSIONS 2.7)
FINDPYTHON()

# Create the package in build dir for testing purpose
SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/spacevecalg)
CONFIGURE_FILE(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/spacevecalg/__init__.py COPYONLY)
SET(OUTPUT_BINDING ${CMAKE_CURRENT_BINARY_DIR}/spacevecalg.cpp)

# Generate Python binding code
ADD_CUSTOM_COMMAND (
  OUTPUT ${OUTPUT_BINDING}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate.py ${OUTPUT_BINDING}
  DEPENDS generate.py
)

# Build the library
SET(SOURCES ${OUTPUT_BINDING})

INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(SYSTEM ${EIGEN3_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../../src)

ADD_LIBRARY(_spacevecalg SHARED ${SOURCES})

target_link_libraries(_spacevecalg ${PYTHON_LIBRARY})

SET_TARGET_PROPERTIES(_spacevecalg PROPERTIES PREFIX "")

# Install rules
INSTALL(TARGETS _spacevecalg DESTINATION "${PYTHON_SITELIB}/spacevecalg")
PYTHON_INSTALL_BUILD(spacevecalg __init__.py "${PYTHON_SITELIB}")
