#!/usr/bin/env python

#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

import unittest

import eigen
import sva

import numpy as np

import sys

TOL = 0.00001

class TestSVAPTransformdVector(unittest.TestCase):
  def test(self):
    create_random_pt = lambda: sva.PTransformd(eigen.Quaterniond(eigen.Vector4d.Random().normalized()), eigen.Vector3d().Random()*100)

    # Check empty creation and filling
    v1 = sva.PTransformdVector()
    for i in range(100):
      v1.append(create_random_pt())
    assert(len(v1) == 100)

    # Check creation from a single PTransformd object
    v2 = sva.PTransformdVector(create_random_pt())
    assert(len(v2) == 1)

    # Check creation from a list
    v3 = sva.PTransformdVector([create_random_pt() for i in range(100)])
    assert(len(v3) == 100)

    # Check creation from a lot of PTransformd
    v3 = sva.PTransformdVector(*[create_random_pt() for i in range(100)])
    assert(len(v3) == 100)

    # Check access
    pt = create_random_pt()
    v4 = sva.PTransformdVector([pt]*100)
    assert(all([pt == vi for vi in v4]))

class TestSVAMotionVecdVector(unittest.TestCase):
  def test(self):
    create_random_mv = lambda: sva.MotionVecd(eigen.Vector3d().Random()*100, eigen.Vector3d().Random()*100)

    # Check empty creation and filling
    v1 = sva.MotionVecdVector()
    for i in range(100):
      v1.append(create_random_mv())
    assert(len(v1) == 100)

    # Check creation from a single MotionVecd object
    v2 = sva.MotionVecdVector(create_random_mv())
    assert(len(v2) == 1)

    # Check creation from a list
    v3 = sva.MotionVecdVector([create_random_mv() for i in range(100)])
    assert(len(v3) == 100)

    # Check creation from a lot of MotionVecd
    v3 = sva.MotionVecdVector(*[create_random_mv() for i in range(100)])
    assert(len(v3) == 100)

    # Check access
    mv = create_random_mv()
    v4 = sva.MotionVecdVector([mv]*100)
    assert(all([mv == vi for vi in v4]))

if __name__ == "__main__":
  suite = unittest.TestSuite()
  suite.addTest(TestSVAPTransformdVector('test'))
  suite.addTest(TestSVAMotionVecdVector('test'))
  unittest.TextTestRunner(verbosity=2).run(suite)
