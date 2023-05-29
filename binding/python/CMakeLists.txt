#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

make_cython_bindings(
  sva
  TARGETS
  SpaceVecAlg
  EXPORT_SOURCES
  sva/c_sva.pxd
  sva/sva.pxd
  sva/__init__.py
  PRIVATE_SOURCES
  include/sva_wrapper.hpp
  sva/sva.pyx
  sva/c_sva_private.pxd
  tests/test_sva_containers.py
  tests/test_sva_inertia.py
  tests/test_sva_pickle.py
  tests/test_sva_ptransform.py
  tests/test_sva_vector.py)
