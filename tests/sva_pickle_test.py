#! /usr/bin/env python2

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
