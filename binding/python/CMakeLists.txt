INCLUDE(../../cmake/python.cmake)

# Look for Python 2.7
SET(Python_ADDITIONAL_VERSIONS 2.7)
FINDPYTHON()

## Define PYTHON_DISTLIB
EXECUTE_PROCESS(
  COMMAND "${PYTHON_EXECUTABLE}" "-c"
  "import sys, os; print os.sep.join(['lib', 'python' + sys.version[:3], 'dist-packages'])"
  OUTPUT_VARIABLE PYTHON_DISTLIB
  ERROR_QUIET)
# Remove final \n of the variable PYTHON_DISTLIB
STRING(REPLACE "\n" "" PYTHON_DISTLIB "${PYTHON_DISTLIB}")


# create the package in build dir for testing purpose
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/spacevecalg)
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/spacevecalg/__init__.py COPYONLY)

set(OUTPUT_BINDING ${CMAKE_CURRENT_BINARY_DIR}/spacevecalg.cpp)

# generate python binding code
add_custom_command (
  OUTPUT ${OUTPUT_BINDING}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate.py ${OUTPUT_BINDING}
  DEPENDS generate.py
)

# build the library
set(SOURCES ${OUTPUT_BINDING})
include_directories(.)
include_directories(../../src)
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})
add_library(_spacevecalg SHARED ${SOURCES})

set_target_properties(_spacevecalg PROPERTIES PREFIX "")

# install rules
set(INSTALL_PATH "${PYTHON_DISTLIB}/spacevecalg/")
install(TARGETS _spacevecalg DESTINATION ${INSTALL_PATH})
install(FILES __init__.py DESTINATION ${INSTALL_PATH})

