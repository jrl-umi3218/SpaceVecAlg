#! /usr/bin/env python

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

import unittest

import pickle

import eigen as e3
import sva

class TestSVAPickle(unittest.TestCase):
  def test(self):
    # register custom pickle function
    sva.copy_reg_pickle()

    mv = sva.MotionVecd(e3.Vector6d.Random())
    fv = sva.ForceVecd(e3.Vector6d.Random())
    pt = sva.PTransformd(e3.Quaterniond(e3.Vector4d.Random().normalized()),
                                        e3.Vector3d.Random())
    rb = sva.RBInertiad(3., e3.Vector3d.Random(), e3.Matrix3d.Random())
    ab = sva.ABInertiad(e3.Matrix3d.Random(), e3.Matrix3d.Random(),
                        e3.Matrix3d.Random())

    def test_pickle(v):
      pickled = pickle.dumps(v)
      v2 = pickle.loads(pickled)
      self.assertEqual(v, v2)

    test_pickle(mv)
    test_pickle(fv)
    test_pickle(pt)
    test_pickle(rb)
    test_pickle(ab)
