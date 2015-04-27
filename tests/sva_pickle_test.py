#! /usr/bin/env python

import sys
import os
sys.path.insert(0, os.path.join("@CMAKE_CURRENT_BINARY_DIR@", '../binding/python'))

import pickle

import eigen3 as e3
import spacevecalg as sva

if __name__ == '__main__':
  # regiter custom pickle function
  sva.copy_reg_pickle()

  mv = sva.MotionVecd(e3.Vector6d.Random())
  fv = sva.ForceVecd(e3.Vector6d.Random())
  pt = sva.PTransformd(e3.Quaterniond(e3.Vector4d.Random().normalized()),
                                      e3.Vector3d.Random())
  rb = sva.RBInertiad(3., e3.Vector3d.Random(), e3.Matrix3d.Random())
  ab = sva.ABInertiad(e3.Matrix3d.Random(), e3.Matrix3d.Random(),
                      e3.Matrix3d.Random())

  def test(v):
    pickled = pickle.dumps(v)
    v2 = pickle.loads(pickled)
    assert(v == v2)

  test(mv)
  test(fv)
  test(pt)
  test(rb)
  test(ab)
