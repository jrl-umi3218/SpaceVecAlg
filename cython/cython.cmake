# Copyright (C) 2019 LAAS-CNRS, JRL AIST-CNRS, INRIA.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

option(PYTHON_BINDING "Generate Python binding" ON)
set(CYTHON_SETUP_IN_PY_LOCATION "${CMAKE_CURRENT_LIST_DIR}/setup.in.py")
set(CYTHON_DUMMY_CPP_LOCATION "${CMAKE_CURRENT_LIST_DIR}/dummy.cpp")
set(PYTHON_EXTRA_CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/python")

# Find the Python packages required depending on binding options
macro(_setup_python_for_cython)
  # FindPython3.cmake only exists from CMake 3.12
  if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
    list(APPEND CMAKE_MODULE_PATH ${PYTHON_EXTRA_CMAKE_MODULE_PATH})
  endif()
  set(PYTHON_VERSION Python3)
  if(PYTHON_BINDING)
    list(APPEND PYTHON_BINDING_VERSIONS Python3)
    if(NOT DEFINED Python3_EXECUTABLE)
      find_program(DEFAULT_PYTHON3_EXECUTABLE python3)
      if(DEFAULT_PYTHON3_EXECUTABLE)
        set(Python3_EXECUTABLE ${DEFAULT_PYTHON3_EXECUTABLE})
      endif()
    endif()

    # --- Standard CMake Virtual Env Management ---
    # set(Python3_FIND_VIRTUALENV "FIRST")
    # ---------------------------------------------

    find_package(${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development
                                                       NumPy)

    # --- Print out the detected paths ---
    message(STATUS "---------------------------------------------------")
    message(STATUS "Python3 Interpreter: ${Python3_EXECUTABLE}")
    message(STATUS "Python3 Include Dirs: ${Python3_INCLUDE_DIRS}")
    message(STATUS "Python3 NumPy Include Dirs: ${Python3_NumPy_INCLUDE_DIRS}")
    message(STATUS "---------------------------------------------------")
  endif()
endmacro()

# This macro adds a dummy shared library target to extract compilation flags
# from an interface library
macro(_CYTHON_DUMMY_TARGET TARGET)
  if(NOT TARGET _cython_dummy_${TARGET})
    add_library(_cython_dummy_${TARGET} SHARED EXCLUDE_FROM_ALL
                "${CYTHON_DUMMY_CPP_LOCATION}")
    target_link_libraries(_cython_dummy_${TARGET} PUBLIC ${TARGET})
    set_target_properties(_cython_dummy_${TARGET} PROPERTIES FOLDER
                                                             "bindings/details")
  endif()
endmacro()

# Check wheter a target is an interface library or not
macro(_is_interface_library TARGET OUT)
  get_target_property(target_type ${TARGET} TYPE)
  if(${target_type} STREQUAL "INTERFACE_LIBRARY")
    set(${OUT} True)
  else()
    set(${OUT} False)
  endif()
endmacro()

# Check whether a target is a static library or not
macro(_is_static_library TARGET OUT)
  get_target_property(target_type ${TARGET} TYPE)
  if(${target_type} STREQUAL "STATIC_LIBRARY")
    set(${OUT} True)
  else()
    set(${OUT} False)
  endif()
endmacro()

# Copy bindings source to build directories and create appropriate target for
# building, installing and testing
macro(
  _ADD_CYTHON_BINDINGS_TARGETS
  PYTHON
  PYTHON_EXECUTABLE
  PACKAGE
  SOURCES
  GENERATE_SOURCES
  TARGETS
  WITH_TESTS)
  set(SETUP_LOCATION
      "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}/${PYTHON}/$<CONFIGURATION>")
  set(${PACKAGE}_${PYTHON}_SETUP_LOCATION
      "${SETUP_LOCATION}"
      CACHE INTERNAL "")
  if(TARGET cython_${PYTHON}_${PACKAGE})
    target_include_directories(cython_${PYTHON}_${PACKAGE}
                               INTERFACE "${SETUP_LOCATION}")
  endif()
  if(DEFINED CMAKE_BUILD_TYPE)
    file(MAKE_DIRECTORY
         "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}/${PYTHON}/${CMAKE_BUILD_TYPE}")
  else()
    foreach(CFG ${CMAKE_CONFIGURATION_TYPES})
      file(MAKE_DIRECTORY
           "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}/${PYTHON}/${CFG}")
    endforeach()
  endif()
  file(
    GENERATE
    OUTPUT "${SETUP_LOCATION}/setup.py"
    INPUT "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}/setup.in.py")
  # Target to build the bindings
  set(TARGET_NAME ${PACKAGE}-${PYTHON}-bindings)
  add_custom_target(
    ${TARGET_NAME} ALL
    COMMAND ${CMAKE_COMMAND} -E chdir "${SETUP_LOCATION}" ${PYTHON_EXECUTABLE}
            setup.py build_ext --inplace
    COMMENT "Generating local ${PACKAGE} ${PYTHON} bindings"
    DEPENDS ${SOURCES} ${GENERATE_SOURCES}
    SOURCES ${SOURCES} ${GENERATE_SOURCES})
  set_target_properties(${TARGET_NAME} PROPERTIES FOLDER "bindings")
  add_dependencies(${TARGET_NAME} ${TARGETS})
  # Copy sources
  foreach(F ${GENERATE_SOURCES})
    file(
      GENERATE
      OUTPUT "${SETUP_LOCATION}/${F}"
      INPUT "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}/configured/${F}")
  endforeach()
  set(I 0)
  foreach(SRC ${SOURCES})
    if(IS_ABSOLUTE ${SRC})
      if(NOT ${SRC} MATCHES "^${CMAKE_CURRENT_BINARY_DIR}")
        message(
          FATAL_ERROR
            "Source provided to ADD_CYTHON_BINDINGS must have a relative path or an absolute path in CMAKE_CURRENT_BINARY_DIR (${CMAKE_CURRENT_BINARY_DIR})"
        )
      endif()
      file(RELATIVE_PATH REL_SRC "${CMAKE_CURRENT_BINARY_DIR}" "${SRC}")
      set(FILE_IN "${SRC}")
      set(FILE_OUT "${SETUP_LOCATION}/${REL_SRC}")
    else()
      set(FILE_IN "${CMAKE_CURRENT_SOURCE_DIR}/${SRC}")
      set(FILE_OUT "${SETUP_LOCATION}/${SRC}")
    endif()
    add_custom_target(
      copy-sources-${I}-${TARGET_NAME}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${FILE_IN} ${FILE_OUT}
      DEPENDS ${FILE_IN})
    set_target_properties(copy-sources-${I}-${TARGET_NAME}
                          PROPERTIES FOLDER "bindings/details")
    add_dependencies(${TARGET_NAME} copy-sources-${I}-${TARGET_NAME})
    math(EXPR I "${I} + 1")
  endforeach()
  # Manual target to force regeneration
  add_custom_target(
    force-${TARGET_NAME}
    COMMAND ${CMAKE_COMMAND} -E chdir "${SETUP_LOCATION}" ${PYTHON_EXECUTABLE}
            setup.py build_ext --inplace --force
    COMMENT "Generating local ${PACKAGE} ${PYTHON} bindings (forced)")
  set_target_properties(force-${TARGET_NAME} PROPERTIES FOLDER "bindings")
  # Tests
  if(${WITH_TESTS} AND ${BUILD_TESTING})
    if(WIN32)
      set(ENV_VAR "PATH")
      set(PATH_SEP ";")
    else()
      set(ENV_VAR "LD_LIBRARY_PATH")
      set(PATH_SEP ":")
    endif()
    set(EXTRA_LD_PATH "")
    foreach(TGT ${TARGETS})
      _is_interface_library(${TGT} IS_INTERFACE)
      if(NOT ${IS_INTERFACE})
        set(EXTRA_LD_PATH
            "$<TARGET_FILE_DIR:${TGT}>${PATH_SEP}${EXTRA_LD_PATH}")
      endif()
    endforeach()
    if(${WITH_TESTS})
      add_test(
        NAME test-${TARGET_NAME}
        COMMAND
          ${CMAKE_COMMAND} -E env "${ENV_VAR}=${EXTRA_LD_PATH}$ENV{${ENV_VAR}}"
          ${CMAKE_COMMAND} -E env "PYTHONPATH=.${PATH_SEP}$ENV{PYTHONPATH}"
          ${CMAKE_COMMAND} -E chdir "${SETUP_LOCATION}" ${PYTHON_EXECUTABLE} -m
          pytest)
    endif()
  endif()
  # Install targets
  if(DEFINED PYTHON_DEB_ROOT)
    add_custom_target(
      install-${TARGET_NAME}
      COMMAND ${CMAKE_COMMAND} -E chdir "${SETUP_LOCATION}" ${PYTHON_EXECUTABLE}
              setup.py install --root=${PYTHON_DEB_ROOT} --install-layout=deb
      COMMENT "Install ${PACKAGE} ${PYTHON} bindings (Debian layout)")
  else()
    set(PIP_EXTRA_OPTIONS "")
    # set(PIP_EXTRA_OPTIONS "--no-system-install")
    add_custom_target(
      install-${TARGET_NAME}
      COMMAND ${CMAKE_COMMAND} -E chdir "${SETUP_LOCATION}" ${PYTHON_EXECUTABLE}
              -m pip install . ${PIP_EXTRA_OPTIONS} --upgrade
      COMMENT "Install ${PACKAGE} ${PYTHON} bindings")
    set_target_properties(install-${TARGET_NAME} PROPERTIES FOLDER "bindings")
    add_custom_target(
      uninstall-${TARGET_NAME}
      COMMAND ${CMAKE_COMMAND} -E chdir "${SETUP_LOCATION}" ${PYTHON_EXECUTABLE}
              -m pip uninstall -y ${PACKAGE}
      COMMENT "Removing ${PACKAGE} ${PYTHON} bindings")
    set_target_properties(uninstall-${TARGET_NAME} PROPERTIES FOLDER "bindings")
    add_dependencies(uninstall uninstall-${TARGET_NAME})
  endif()
  install(
    CODE "EXECUTE_PROCESS(COMMAND \"${CMAKE_COMMAND}\" --build \"${CMAKE_BINARY_DIR}\" --config \${CMAKE_INSTALL_CONFIG_NAME} --target install-${TARGET_NAME})"
  )
endmacro()

# .rst: .. command:: ADD_CYTHON_BINDINGS(PACKAGE TARGETS targets... [VERSION
# version] [MODULES modules...] [EXPORT_SOURCES sources...] [PRIVATE_SOURCES
# ...] [GENERATE_SOURCES ...])
#
# This macro add cython bindings using one or more libraries built by the
# project.
#
# :PACKAGE:          Name of the Python package
#
# :TARGETS:          Name of the targets that the bindings should link to
#
# :VERSION:          Version of the bindings, defaults to ``PROJECT_VERSION``
#
# :MODULES:          Python modules built by this macro call. Defaults to
# ``PACKAGE.PACKAGE``
#
# :EXPORT_SOURCES:   Sources that will be installed along with the package
# (typically, public pxd files and __init__.py)
#
# :PRIVATE_SOURCES:  Sources that are needed to built the package but will not
# be installed
#
# :GENERATE_SOURCES: Sources that will be configured and then generated in the
# correct location, the generated files are then considered as PRIVATE_SOURCES
#
# The macro will generate a setup.py script in
# ``$CMAKE_CURRENT_BINARY_DIR/$PACKAGE/$PYTHON/$<CONFIGURATION>`` and copy the
# provided sources in this location. Relative paths are preferred to provide
# sources but one can use absolute paths if and only if the absolute path starts
# with ``$CMAKE_CURRENT_BINARY_DIR``
#
macro(ADD_CYTHON_BINDINGS PACKAGE)
  set(options)
  set(oneValueArgs VERSION)
  set(multiValueArgs MODULES TARGETS EXPORT_SOURCES PRIVATE_SOURCES
                     GENERATE_SOURCES)
  cmake_parse_arguments(CYTHON_BINDINGS "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" ${ARGN})
  if(NOT DEFINED CYTHON_BINDINGS_VERSION)
    set(CYTHON_BINDINGS_VERSION ${PROJECT_VERSION})
  endif()
  if(NOT DEFINED CYTHON_BINDINGS_EXPORT_SOURCES)
    set(CYTHON_BINDINGS_EXPORT_SOURCES)
  endif()
  if(NOT DEFINED CYTHON_BINDINGS_PRIVATE_SOURCES)
    set(CYTHON_BINDINGS_PRIVATE_SOURCES)
  endif()
  if(NOT DEFINED CYTHON_BINDINGS_GENERATE_SOURCES)
    set(CYTHON_BINDINGS_GENERATE_SOURCES)
  endif()
  if(NOT DEFINED CYTHON_BINDINGS_MODULES)
    set(CYTHON_BINDINGS_MODULES "${PACKAGE}.${PACKAGE}")
  endif()
  if(NOT DEFINED CYTHON_BINDINGS_TARGETS)
    message(
      FATAL_ERROR
        "Error in ADD_CYTHON_BINDINGS, bindings should depend on at least one target"
    )
  endif()
  # Setup the basic setup script
  set(CYTHON_BINDINGS_SOURCES)
  list(APPEND CYTHON_BINDINGS_SOURCES ${CYTHON_BINDINGS_EXPORT_SOURCES})
  list(APPEND CYTHON_BINDINGS_SOURCES ${CYTHON_BINDINGS_PRIVATE_SOURCES})
  set(WITH_TESTS False)
  foreach(SRC ${CYTHON_BINDINGS_SOURCES})
    if(${SRC} MATCHES "^tests/")
      set(WITH_TESTS True)
    endif()
  endforeach()
  set(CYTHON_BINDINGS_PACKAGE_NAME ${PACKAGE})
  set(CYTHON_BINDINGS_COMPILE_DEFINITIONS)
  set(CYTHON_BINDINGS_CXX_STANDARD)
  set(CYTHON_BINDINGS_INCLUDE_DIRECTORIES)
  set(CYTHON_BINDINGS_LINK_FLAGS)
  set(CYTHON_BINDINGS_LIBRARIES)
  set(CYTHON_BINDINGS_STATIC_LIBRARIES)
  set(CYTHON_BINDINGS_TARGET_FILES)
  foreach(TGT ${CYTHON_BINDINGS_TARGETS})
    _is_interface_library(${TGT} IS_INTERFACE)
    if(${IS_INTERFACE})
      _cython_dummy_target(${TGT})
      list(APPEND CYTHON_BINDINGS_COMPILE_DEFINITIONS
           "$<TARGET_PROPERTY:_cython_dummy_${TGT},COMPILE_DEFINITIONS>")
      list(
        APPEND CYTHON_BINDINGS_COMPILE_DEFINITIONS
        "$<TARGET_PROPERTY:_cython_dummy_${TGT},INTERFACE_COMPILE_DEFINITIONS>")
      list(APPEND CYTHON_BINDINGS_CXX_STANDARD
           "$<TARGET_PROPERTY:_cython_dummy_${TGT},CXX_STANDARD>")
      list(APPEND CYTHON_BINDINGS_INCLUDE_DIRECTORIES
           "$<TARGET_PROPERTY:_cython_dummy_${TGT},INCLUDE_DIRECTORIES>")
      list(
        APPEND CYTHON_BINDINGS_INCLUDE_DIRECTORIES
        "$<TARGET_PROPERTY:_cython_dummy_${TGT},INTERFACE_INCLUDE_DIRECTORIES>")
      list(APPEND CYTHON_BINDINGS_LINK_FLAGS
           "$<TARGET_PROPERTY:_cython_dummy_${TGT},LINK_FLAGS>")
    else()
      _is_static_library(${TGT} IS_STATIC)
      list(APPEND CYTHON_BINDINGS_COMPILE_DEFINITIONS
           "$<TARGET_PROPERTY:${TGT},COMPILE_DEFINITIONS>")
      list(APPEND CYTHON_BINDINGS_COMPILE_DEFINITIONS
           "$<TARGET_PROPERTY:${TGT},INTERFACE_COMPILE_DEFINITIONS>")
      list(APPEND CYTHON_BINDINGS_CXX_STANDARD
           "$<TARGET_PROPERTY:${TGT},CXX_STANDARD>")
      list(APPEND CYTHON_BINDINGS_INCLUDE_DIRECTORIES
           "$<TARGET_PROPERTY:${TGT},INCLUDE_DIRECTORIES>")
      list(APPEND CYTHON_BINDINGS_INCLUDE_DIRECTORIES
           "$<TARGET_PROPERTY:${TGT},INTERFACE_INCLUDE_DIRECTORIES>")
      list(APPEND CYTHON_BINDINGS_LINK_FLAGS
           "$<TARGET_PROPERTY:${TGT},LINK_FLAGS>")
      list(APPEND CYTHON_BINDINGS_LIBRARIES "$<TARGET_LINKER_FILE:${TGT}>")
      list(APPEND CYTHON_BINDINGS_TARGET_FILES "$<TARGET_LINKER_FILE:${TGT}>")
      if(${IS_STATIC})
        list(APPEND CYTHON_BINDINGS_STATIC_LIBRARIES
             "$<TARGET_LINKER_FILE:${TGT}>")
      endif()
    endif()
  endforeach()
  configure_file("${CYTHON_SETUP_IN_PY_LOCATION}"
                 "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}/setup.in.py")
  foreach(F ${CYTHON_BINDINGS_GENERATE_SOURCES})
    configure_file(${F}
                   "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}/configured/${F}")
  endforeach()
  _add_cython_bindings_targets(
    "python3"
    ${Python3_EXECUTABLE}
    ${PACKAGE}
    "${CYTHON_BINDINGS_SOURCES}"
    "${CYTHON_BINDINGS_GENERATE_SOURCES}"
    "${CYTHON_BINDINGS_TARGETS}"
    ${WITH_TESTS})
endmacro()

# In this macro PYTHON is the module we should search and PYTHON_B is the name
# for the bindings
macro(_MAKE_CYTHON_LIBRARY PACKAGE PYTHON PYTHON_B OUT)
  set(SETUP_LOCATION_VAR ${PACKAGE}_${PYTHON_B}_SETUP_LOCATION)
  set(TGT_NAME cython_${PYTHON_B}_${PACKAGE})
  set(${OUT} ${TGT_NAME})
  if(NOT TARGET ${TGT_NAME})
    set(${PYTHON}_FIND_VERSION_COUNT 3)
    set(${PYTHON}_FIND_VERSION_EXACT TRUE)
    execute_process(
      COMMAND ${PYTHON_B} -c "import sys; print(sys.version_info.major);"
      OUTPUT_VARIABLE ${PYTHON}_FIND_VERSION_MAJOR
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
      COMMAND ${PYTHON_B} -c "import sys; print(sys.version_info.minor);"
      OUTPUT_VARIABLE ${PYTHON}_FIND_VERSION_MINOR
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    find_package(${PYTHON} REQUIRED COMPONENTS Interpreter Development)
    add_library(${TGT_NAME} INTERFACE)
    target_link_libraries(${TGT_NAME} INTERFACE ${${PYTHON}_LIBRARIES})
    target_include_directories(${TGT_NAME}
                               INTERFACE "${${PYTHON}_INCLUDE_DIRS}")
    if(DEFINED ${SETUP_LOCATION_VAR})
      set(SETUP_LOCATION "${${SETUP_LOCATION_VAR}}")
      target_include_directories(${TGT_NAME} INTERFACE "${SETUP_LOCATION}")
    endif()
    add_dependencies(${TGT_NAME} ${PACKAGE}-${PYTHON_B}-bindings)
  endif()
endmacro()

macro(_APPEND_CYTHON_LIBRARY PACKAGE PYTHON PYTHON_B OUT)
  _make_cython_library(${PACKAGE} ${PYTHON} ${PYTHON_B} LIB)
  list(APPEND ${OUT} ${LIB})
endmacro()

# .rst: .. command:: GET_CYTHON_LIBRARIES(PACKAGE VAR)
#
# This macro search Python versions according to the specified bindings settings
# then returns appropriate targets in the provided VAR variable
#
# It creates interface targets that include the generated bindings directory and
# link to the correct Python version
#
macro(GET_CYTHON_LIBRARIES PACKAGE VAR)
  # FindPython(2|3).cmake only exists from CMake 3.12
  if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
    list(APPEND CMAKE_MODULE_PATH ${PYTHON_EXTRA_CMAKE_MODULE_PATH})
  endif()
  set(${VAR})
  _append_cython_library(${PACKAGE} Python3 python3 ${VAR})
endmacro()

# .rst: .. command:: GET_PYTHON_NAMES(VAR)
#
# This macro returns the names of Python versions according to the specified
# bindings
#
macro(GET_PYTHON_NAMES VAR)
  set(${VAR})
  list(APPEND ${VAR} Python3)
endmacro()

# .rst: .. command:: MAKE_CYTHON_BINDINGS(PACKAGE TARGETS targets... [VERSION
# version] [MODULES modules...] [EXPORT_SOURCES sources...] [PRIVATE_SOURCES
# ...] [GENERATE_SOURCES ...])
#
# This function adds cython bindings using one or more libraries built by the
# project. It is similar to ADD_CYTHON_BINDINGS but the process is entirely
# handled by CMake which gives us better incremental builds.
#
# For each module ``a.b.c`` it creates a target name ``c_${PYTHON_VERSION}`` for
# each ``PYTHON_VERSION`` in ``PYTHON_BINDING_VERSIONS``
#
# :PACKAGE:          Name of the Python package
#
# :TARGETS:          Name of the targets that the bindings should link to
#
# :VERSION:          Version of the bindings, defaults to ``PROJECT_VERSION``
#
# :MODULES:          Python modules built by this macro call. Defaults to
# ``PACKAGE.PACKAGE``
#
# :EXPORT_SOURCES:   Sources that will be installed along with the package
# (typically, public pxd files and __init__.py)
#
# :PRIVATE_SOURCES:  Sources that are needed to built the package but will not
# be installed
#
# :GENERATE_SOURCES: Sources that will be configured and then generated in the
# correct location, the generated files are installed unless they are test files
#
# The macro will generate a setup.py script in
# ``$CMAKE_CURRENT_BINARY_DIR/$PACKAGE/$PYTHON/$<CONFIGURATION>`` and copy the
# provided sources in this location. Relative paths are preferred to provide
# sources but one can use absolute paths if and only if the absolute path starts
# with ``$CMAKE_CURRENT_BINARY_DIR``
#
function(MAKE_CYTHON_BINDINGS PACKAGE)
  set(options)
  set(oneValueArgs VERSION)
  set(multiValueArgs MODULES TARGETS EXPORT_SOURCES PRIVATE_SOURCES
                     GENERATE_SOURCES)
  cmake_parse_arguments(CYTHON_BINDINGS "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" ${ARGN})
  if(NOT DEFINED CYTHON_BINDINGS_VERSION)
    set(CYTHON_BINDINGS_VERSION ${PROJECT_VERSION})
  endif()
  if(NOT DEFINED CYTHON_BINDINGS_EXPORT_SOURCES)
    set(CYTHON_BINDINGS_EXPORT_SOURCES)
  endif()
  if(NOT DEFINED CYTHON_BINDINGS_PRIVATE_SOURCES)
    set(CYTHON_BINDINGS_PRIVATE_SOURCES)
  endif()
  if(NOT DEFINED CYTHON_BINDINGS_GENERATE_SOURCES)
    set(CYTHON_BINDINGS_GENERATE_SOURCES)
  endif()
  if(NOT DEFINED CYTHON_BINDINGS_MODULES)
    set(CYTHON_BINDINGS_MODULES "${PACKAGE}.${PACKAGE}")
  endif()
  if(NOT DEFINED CYTHON_BINDINGS_TARGETS)
    message(
      FATAL_ERROR
        "Error in ADD_CYTHON_BINDINGS, bindings should depend on at least one target"
    )
  endif()
  set(CYTHON_BINDINGS_SOURCES)
  list(APPEND CYTHON_BINDINGS_SOURCES ${CYTHON_BINDINGS_EXPORT_SOURCES})
  list(APPEND CYTHON_BINDINGS_SOURCES ${CYTHON_BINDINGS_PRIVATE_SOURCES})
  list(APPEND CYTHON_BINDINGS_SOURCES ${CYTHON_BINDINGS_GENERATE_SOURCES})
  set(CYTHON_BINDINGS_COMPILE_SOURCES)
  set(WITH_TESTS False)
  foreach(SRC ${CYTHON_BINDINGS_SOURCES})
    if(${SRC} MATCHES "^tests/")
      set(WITH_TESTS True)
    endif()
    if(${SRC} MATCHES ".pyx$")
      list(APPEND CYTHON_BINDINGS_COMPILE_SOURCES ${SRC})
    endif()
  endforeach()
  add_library(_cython_dummy_${PACKAGE} SHARED EXCLUDE_FROM_ALL
              "${CYTHON_DUMMY_CPP_LOCATION}")
  target_link_libraries(_cython_dummy_${PACKAGE}
                        INTERFACE ${CYTHON_BINDINGS_TARGETS})
  set_target_properties(_cython_dummy_${PACKAGE} PROPERTIES FOLDER
                                                            "bindings/details")
  foreach(PYTHON ${PYTHON_BINDING_VERSIONS})
    set(PACKAGE_OUTPUT_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON}/$<CONFIG>/${PACKAGE})
    if(DEFINED PYTHON_DEB_ROOT)
      execute_process(
        COMMAND
          ${${PYTHON}_EXECUTABLE} -c
          "import sysconfig; print(sysconfig.get_path('platlib', vars={'base': ''}))"
        RESULT_VARIABLE PYTHON_INSTALL_DESTINATION_FOUND
        OUTPUT_VARIABLE PYTHON_INSTALL_DESTINATION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    else()
      execute_process(
        COMMAND ${${PYTHON}_EXECUTABLE} -c
                "import sysconfig; print(sysconfig.get_path('purelib'))"
        RESULT_VARIABLE PYTHON_INSTALL_DESTINATION_FOUND
        OUTPUT_VARIABLE PYTHON_INSTALL_DESTINATION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      # Debian/Ubuntu has a specific problem here See
      # https://github.com/mesonbuild/meson/issues/8739 for an overview of the
      # problem
      if(EXISTS /etc/debian_version)
        execute_process(
          COMMAND
            ${Python3_EXECUTABLE} -c
            "import sys; print(\"python{}.{}\".format(sys.version_info.major, sys.version_info.minor));"
          OUTPUT_VARIABLE PYTHON_VERSION
          OUTPUT_STRIP_TRAILING_WHITESPACE)
        string(REPLACE "python3/" "${PYTHON_VERSION}/"
                       PYTHON_INSTALL_DESTINATION
                       "${PYTHON_INSTALL_DESTINATION}")
      endif()
    endif()
    foreach(F ${CYTHON_BINDINGS_GENERATE_SOURCES})
      configure_file(${F} ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON}/cmake/${F})
      file(
        GENERATE
        OUTPUT ${PACKAGE_OUTPUT_DIRECTORY}/${F}
        INPUT ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON}/cmake/${F})
    endforeach()
    foreach(F ${CYTHON_BINDINGS_EXPORT_SOURCES})
      file(
        GENERATE
        OUTPUT ${PACKAGE_OUTPUT_DIRECTORY}/${F}
        INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${F})
    endforeach()
    foreach(F ${CYTHON_BINDINGS_PRIVATE_SOURCES})
      if(${F} MATCHES "^tests/")
        file(
          GENERATE
          OUTPUT ${PACKAGE_OUTPUT_DIRECTORY}/${F}
          INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${F})
      endif()
    endforeach()
    install(
      DIRECTORY ${PACKAGE_OUTPUT_DIRECTORY}/
      DESTINATION ${PYTHON_INSTALL_DESTINATION}
      # We can't use PACKAGE_OUTPUT_DIRECTORY because it contains a
      # generator-expression
      REGEX "^${CMAKE_CURRENT_BINARY_DIR}/${PYTHON}/[A-z]*/${PACKAGE}/tests.*"
            EXCLUDE
      PATTERN ".pytest_cache/*" EXCLUDE
      PATTERN "__pycache__/*" EXCLUDE)
    # Make an uninstall rule that:
    #
    # * Remove the installed module fully (including the empty directory)
    # * Remove trace of bindings that could have been installed by
    #   add_cython_bindings in the past
    set(UNINSTALL_TARGET_NAME uninstall-${PACKAGE}-${PYTHON}-bindings)
    add_custom_target(
      ${UNINSTALL_TARGET_NAME}
      COMMAND ${CMAKE_COMMAND} -E rm -rf
              ${PYTHON_INSTALL_DESTINATION}/${PACKAGE}*.dist-info)
    add_dependencies(uninstall ${UNINSTALL_TARGET_NAME})
    if(WITH_TESTS AND BUILD_TESTING)
      if(WIN32)
        set(ENV_VAR "PATH")
        set(PATH_SEP ";")
      else()
        set(ENV_VAR "LD_LIBRARY_PATH")
        set(PATH_SEP ":")
      endif()
      set(EXTRA_LD_PATH "")
      foreach(TGT ${CYTHON_BINDINGS_TARGETS})
        _is_interface_library(${TGT} IS_INTERFACE)
        if(NOT ${IS_INTERFACE})
          set(EXTRA_LD_PATH
              "$<TARGET_FILE_DIR:${TGT}>${PATH_SEP}${EXTRA_LD_PATH}")
        endif()
      endforeach()
      add_test(
        NAME test-${PACKAGE}-${PYTHON}-bindings
        COMMAND
          ${CMAKE_COMMAND} -E env "${ENV_VAR}=${EXTRA_LD_PATH}$ENV{${ENV_VAR}}"
          ${CMAKE_COMMAND} -E env "PYTHONPATH=.${PATH_SEP}$ENV{PYTHONPATH}"
          ${CMAKE_COMMAND} -E chdir "${PACKAGE_OUTPUT_DIRECTORY}"
          ${${PYTHON}_EXECUTABLE} -m pytest)
    endif()
    foreach(MOD ${CYTHON_BINDINGS_MODULES})
      string(REPLACE "." "/" SRC ${MOD})
      set(SRC "${SRC}.pyx")
      if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${SRC}")
        message(
          FATAL_ERROR "Expected to find ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}")
      endif()
      string(REGEX REPLACE ".pyx$" ".cpp" SRC_CPP ${SRC})
      string(REGEX REPLACE "/[^/]*$" "" SRC_DIR ${SRC})
      string(REGEX REPLACE "^(.*)\\..*$" "\\1" LIB_NAME ${MOD})
      string(REGEX REPLACE "\\." "_" LIB_NAME ${LIB_NAME})
      string(REGEX REPLACE "^.*\\.(.*)$" "\\1" LIB_OUTPUT_NAME ${MOD})
      string(REGEX REPLACE "^([^/]*)/.*$" "\\1" MOD_FOLDER ${SRC})
      set(MOD_OUTPUT_DIRECTORY ${PACKAGE_OUTPUT_DIRECTORY}/${SRC_DIR})
      set(CPP_OUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON}/${SRC_DIR})
      set(CPP_OUT ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON}/${SRC_CPP})
      file(MAKE_DIRECTORY ${CPP_OUT_DIR})
      add_custom_command(
        OUTPUT ${CPP_OUT}
        COMMAND
          ${${PYTHON}_EXECUTABLE} -m cython --cplus -o ${CPP_OUT}
          "-I$<JOIN:$<REMOVE_DUPLICATES:$<TARGET_PROPERTY:_cython_dummy_${PACKAGE},INCLUDE_DIRECTORIES>>,;-I>"
          -I${CMAKE_CURRENT_SOURCE_DIR}/include
          ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}
        DEPENDS ${CYTHON_BINDINGS_SOURCES} ${CYTHON_BINDINGS_TARGETS}
        COMMAND_EXPAND_LISTS)
      set(TARGET_NAME ${LIB_NAME}_${PYTHON})
      if(${PYTHON} STREQUAL "Python3")
        python3_add_library(${TARGET_NAME} MODULE ${CPP_OUT})
      else()
        message(FATAL_ERROR "Unknown Python value: ${PYTHON}")
      endif()
      # Cython is likely to generate code that won't compile without warnings
      if(UNIX AND NOT DEFINED CXX_DISABLE_WERROR)
        target_compile_options(${TARGET_NAME} PRIVATE -Wno-error)
      endif()
      if(UNIX)
        # Cython does a lot of casts that remove the const qualifier
        target_compile_options(${TARGET_NAME} PRIVATE -Wno-cast-qual)
        # Cython usually includes the deprecated NumPy API
        target_compile_options(${TARGET_NAME} PRIVATE -Wno-cpp)
        # Cython does some fishy conversions
        target_compile_options(${TARGET_NAME} PRIVATE -Wno-conversion
                                                      -Wno-overflow)
        # Generating API might look like unusued variables
        target_compile_options(${TARGET_NAME} PRIVATE -Wno-unused-variable
                                                      -Wno-unused-function)
      endif()
      target_include_directories(${TARGET_NAME}
                                 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
      target_include_directories(
        ${TARGET_NAME} INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON})
      target_link_libraries(
        ${TARGET_NAME} PUBLIC ${CYTHON_BINDINGS_TARGETS} ${PYTHON}::Python
                              ${PYTHON}::NumPy)
      set_target_properties(
        ${TARGET_NAME}
        PROPERTIES CXX_VISIBILITY_PRESET default
                   PREFIX ""
                   DEBUG_POSTFIX ""
                   OUTPUT_NAME ${LIB_OUTPUT_NAME}
                   LIBRARY_OUTPUT_DIRECTORY ${MOD_OUTPUT_DIRECTORY}
                   RUNTIME_OUTPUT_DIRECTORY ${MOD_OUTPUT_DIRECTORY})
      if(NOT TARGET ${UNINSTALL_TARGET_NAME}-${MOD_FOLDER})
        add_custom_target(
          ${UNINSTALL_TARGET_NAME}-${MOD_FOLDER}
          COMMAND ${CMAKE_COMMAND} -E rm -rf
                  ${PYTHON_INSTALL_DESTINATION}/${MOD_FOLDER})
        add_dependencies(${UNINSTALL_TARGET_NAME}
                         ${UNINSTALL_TARGET_NAME}-${MOD_FOLDER})
      endif()
    endforeach()
  endforeach()
endfunction()
