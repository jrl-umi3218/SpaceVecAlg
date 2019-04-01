#! /usr/bin/env python

#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

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
