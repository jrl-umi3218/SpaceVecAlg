# We are now ready to compile the actual extension module
nanobind_add_module(
  # Name of the extension
  sva_pywrap
  # Target the stable ABI for Python 3.12+, which reduces the number of binary
  # wheels that must be built. This does nothing on older Python versions
  STABLE_ABI
  # Build libnanobind statically and merge it into the extension (which itself
  # remains a shared library)
  #
  # If your project builds multiple extensions, you can replace this flag by
  # NB_SHARED to conserve space by reusing a shared libnanobind across libraries
  NB_STATIC
  LTO
  NB_SUPPRESS_WARNINGS
  # Source code goes here
  src/sva/sva_module.cpp
  src/sva/bind_PTransformd.cpp
  src/sva/bind_ForceVecd.cpp
  src/sva/bind_MotionVecd.cpp
  src/sva/bind_ABInertiad.cpp
  src/sva/bind_RBInertiad.cpp)
jrl_check_python_module_name(sva_pywrap)
target_link_libraries(sva_pywrap PRIVATE SpaceVecAlg)
target_link_libraries(sva_pywrap PRIVATE nanoeigenpy::nanoeigenpy_headers)
jrl_target_set_default_compile_options(sva_pywrap PRIVATE)
jrl_target_set_output_directory(sva_pywrap OUTPUT_DIRECTORY
                                ${CMAKE_BINARY_DIR}/lib/site-packages/sva)

# XXX: as-is the default generated __init__.py includes the package as
# sva.sva_pywrap.<classname>. I prefer all symbols to be in the public namespace
# sva.<classname> thus we install our own jrl_python_generate_init_py(
# sva_pywrap OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib/site-packages/sva/__init__.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.conf.py
               ${CMAKE_BINARY_DIR}/lib/site-packages/sva/__init__.py COPYONLY)

# This is needed for LSP to recognise sva.<symbol> instead of
# sva.sva_pywrap.<symbol>
file(WRITE ${CMAKE_BINARY_DIR}/lib/site-packages/sva/__init__.pyi
     "from .sva_pywrap import *\n")

# Display output of stub generation
set(CMAKE_EXECUTE_PROCESS_COMMAND_ECHO STDOUT)

nanobind_add_stub(
  sva_pywrap_stub
  MODULE
  sva.sva_pywrap
  OUTPUT
  ${CMAKE_BINARY_DIR}/lib/site-packages/sva/sva_pywrap.pyi
  PYTHON_PATH
  ${CMAKE_BINARY_DIR}/lib/site-packages
  MARKER_FILE
  ${CMAKE_BINARY_DIR}/lib/site-packages/sva/py.typed
  DEPENDS
  sva_pywrap
  VERBOSE)

jrl_python_compile_all(DIRECTORY ${CMAKE_BINARY_DIR}/lib/site-packages/sva
                       VERBOSE)

jrl_python_compute_install_dir(python_install_dir)
install(
  DIRECTORY ${CMAKE_BINARY_DIR}/lib/site-packages/sva
  DESTINATION ${python_install_dir}
  FILES_MATCHING
  PATTERN "*.py"
  PATTERN "*.pyc"
  PATTERN "*.pyi"
  PATTERN "*.typed"
  PATTERN "*.so"
  PATTERN "*.pyd")
