#!/usr/bin/env python

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

import eigen
import sva

import numpy as np

import sys

TOL = 0.00001

def isUpperNull(m):
  assert(isinstance(m, eigen.Matrix3d))
  return m[0,1] == m[0,2] == m[1,2] == 0

class TestRBInertiad(unittest.TestCase):
  def test(self):
    mass = 1.
    I = eigen.Matrix3d([[1, 2, 3],[2, 1, 4], [3, 4, 1]])
    h = eigen.Vector3d.Random()

    rb2 = sva.RBInertiad(mass, h, I)
    self.assertEqual(rb2.mass(), mass)
    self.assertEqual(rb2.momentum(), h)
    self.assertEqual(rb2.inertia(), I)
    self.assertTrue(isUpperNull(rb2.lowerTriangularInertia()))

    rb3 = sva.RBInertiad(mass, h, rb2.lowerTriangularInertia())
    self.assertEqual(rb3.mass(), mass)
    self.assertEqual(rb3.momentum(), h)
    self.assertEqual(rb3.inertia(), I)
    self.assertTrue(isUpperNull(rb3.lowerTriangularInertia()))

    rb4 = rb2 + rb3
    self.assertEqual(rb4.mass(), mass + mass)
    self.assertEqual(rb4.momentum(), h + h)
    self.assertEqual(rb4.inertia(), I + I)
    self.assertTrue(isUpperNull(rb4.lowerTriangularInertia()))

    rb5 = 2.*rb2
    self.assertEqual(rb5.mass(), 2.*mass)
    self.assertEqual(rb5.momentum(), 2.*h)
    self.assertEqual(rb5.inertia(), 2.*I)
    self.assertTrue(isUpperNull(rb5.lowerTriangularInertia()))

    rb6 = rb2*2.
    self.assertEqual(rb6.mass(), 2.*mass)
    self.assertEqual(rb6.momentum(), 2.*h)
    self.assertEqual(rb6.inertia(), 2.*I)
    self.assertTrue(isUpperNull(rb6.lowerTriangularInertia()))

    rb7 = rb2 - rb3
    self.assertEqual(rb7.mass(), mass - mass)
    self.assertEqual(rb7.momentum(), h - h)
    self.assertEqual(rb7.inertia(), I - I)
    self.assertTrue(isUpperNull(rb7.lowerTriangularInertia()))

    rb8 = -rb2
    self.assertEqual(rb8, rb2*-1)
    self.assertTrue(isUpperNull(rb8.lowerTriangularInertia()))

    rb9 = sva.RBInertiad(rb2)
    rb9 += rb3
    self.assertEqual(rb9, rb2 + rb3)
    self.assertTrue(isUpperNull(rb9.lowerTriangularInertia()))

    rb10 = sva.RBInertiad(rb2)
    rb10 -= rb3
    self.assertEqual(rb10, rb2 - rb3)
    self.assertTrue(isUpperNull(rb10.lowerTriangularInertia()))

    self.assertEqual(rb2, rb2)
    self.assertNotEqual(rb2, rb6)

    self.assertTrue(rb2 != rb6)
    self.assertTrue(not(rb2 != rb2))

class TestABInertiad(unittest.TestCase):
  def test(self):
    M = eigen.Matrix3d([[1, 2, 3],[2, 1, 4], [3, 4, 1]])
    H = eigen.Matrix3d.Random()
    I = eigen.Matrix3d([[1, 2, 3],[2, 1, 4], [3, 4, 1]])

    ab1 = sva.ABInertiad(M, H, I)
    self.assertEqual(ab1.massMatrix(), M)
    self.assertTrue(isUpperNull(ab1.lowerTriangularMassMatrix()))
    self.assertEqual(ab1.gInertia(), H)
    self.assertEqual(ab1.inertia(), I)
    self.assertTrue(isUpperNull(ab1.lowerTriangularInertia()))

    ab2 = sva.ABInertiad(ab1.lowerTriangularMassMatrix(), H, ab1.lowerTriangularInertia())
    self.assertEqual(ab2.massMatrix(), M)
    self.assertTrue(isUpperNull(ab2.lowerTriangularMassMatrix()))
    self.assertEqual(ab2.gInertia(), H)
    self.assertEqual(ab2.inertia(), I)
    self.assertTrue(isUpperNull(ab2.lowerTriangularInertia()))

    ab3 = ab1 + ab2
    self.assertEqual(ab3.massMatrix(), M + M)
    self.assertTrue(isUpperNull(ab3.lowerTriangularMassMatrix()))
    self.assertEqual(ab3.gInertia(), H + H)
    self.assertEqual(ab3.inertia(), I + I)
    self.assertTrue(isUpperNull(ab3.lowerTriangularInertia()))

    ab4 = 2.*ab2
    self.assertEqual(ab4.massMatrix(), 2.*M)
    self.assertTrue(isUpperNull(ab4.lowerTriangularMassMatrix()))
    self.assertEqual(ab4.gInertia(), 2.*H)
    self.assertEqual(ab4.inertia(), 2.*I)
    self.assertTrue(isUpperNull(ab4.lowerTriangularInertia()))

    ab5 = ab2*2.
    self.assertEqual(ab5.massMatrix(), 2.*M)
    self.assertTrue(isUpperNull(ab5.lowerTriangularMassMatrix()))
    self.assertEqual(ab5.gInertia(), 2.*H)
    self.assertEqual(ab5.inertia(), 2.*I)
    self.assertTrue(isUpperNull(ab5.lowerTriangularInertia()))

    ab6 = ab1 - ab2
    self.assertEqual(ab6.massMatrix(), M - M)
    self.assertTrue(isUpperNull(ab6.lowerTriangularMassMatrix()))
    self.assertEqual(ab6.gInertia(), H - H)
    self.assertEqual(ab6.inertia(), I - I)
    self.assertTrue(isUpperNull(ab6.lowerTriangularInertia()))

    ab7 = -ab1
    self.assertEqual(ab7, -1.*ab1)
    self.assertTrue(isUpperNull(ab7.lowerTriangularMassMatrix()))
    self.assertTrue(isUpperNull(ab7.lowerTriangularInertia()))

    ab8 = sva.ABInertiad(ab1)
    ab8 += ab2
    self.assertEqual(ab8, ab1 + ab2)
    self.assertTrue(isUpperNull(ab8.lowerTriangularMassMatrix()))
    self.assertTrue(isUpperNull(ab8.lowerTriangularInertia()))

    ab9 = sva.ABInertiad(ab1)
    ab9 -= ab2
    self.assertEqual(ab9, ab1 - ab2)
    self.assertTrue(isUpperNull(ab9.lowerTriangularMassMatrix()))
    self.assertTrue(isUpperNull(ab9.lowerTriangularInertia()))

    self.assertEqual(ab2, ab2)
    self.assertNotEqual(ab2, ab5)

    self.assertTrue(ab2 != ab5)
    self.assertTrue(not(ab2 != ab2))

class TestRBInertiadLeftOperatorsTest(unittest.TestCase):
  def test(self):
    mass = 1.
    I = eigen.Matrix3d([[1, 2, 3],[2, 1, 4], [3, 4, 1]])
    h = eigen.Vector3d.Random()
    rb = sva.RBInertiad(mass, h, I)
    rb6d = rb.matrix()

    w = eigen.Vector3d.Random()
    v = eigen.Vector3d.Random()
    mVec = sva.MotionVecd(w, v)
    mVec6d = mVec.vector()

    fVec = rb*mVec
    fVec6d = rb6d*mVec6d
    self.assertTrue(isinstance(fVec, sva.ForceVecd))
    self.assertAlmostEqual((fVec6d - fVec.vector()).norm(), 0, delta = TOL)

    # Skip vectorized version test

class TestABInertiadLeftOperatorsTest(unittest.TestCase):
  def test(self):
    M = eigen.Matrix3d([[1, 2, 3],[2, 1, 4], [3, 4, 1]])
    H = eigen.Matrix3d.Random()
    I = eigen.Matrix3d([[1, 2, 3],[2, 1, 4], [3, 4, 1]])

    ab = sva.ABInertiad(M, H, I)
    ab6d = ab.matrix()

    mass = 1.
    h = eigen.Vector3d.Random() * 100
    rb = sva.RBInertiad(mass, h, I)
    rb6d = rb.matrix()

    w = eigen.Vector3d.Random()
    v = eigen.Vector3d.Random()
    mVec = sva.MotionVecd(w, v)
    mVec6d = mVec.vector()

    abRes = ab + rb
    abRes6d = ab6d + rb6d
    self.assertTrue(isinstance(abRes, sva.ABInertiad))
    self.assertAlmostEqual((abRes6d - abRes.matrix()).norm(), 0, delta = TOL)
    self.assertTrue(isUpperNull(abRes.lowerTriangularMassMatrix()))
    self.assertTrue(isUpperNull(abRes.lowerTriangularInertia()))

    fVec = ab*mVec
    fVec6d = ab6d*mVec6d
    self.assertTrue(isinstance(fVec, sva.ForceVecd))
    self.assertAlmostEqual((fVec6d - fVec.vector()).norm(), 0, delta = TOL)

    # Skip vectorized version test

if __name__ == "__main__":
  suite = unittest.TestSuite()
  suite.addTest(TestRBInertiad('test'))
  suite.addTest(TestABInertiad('test'))
  suite.addTest(TestRBInertiadLeftOperatorsTest('test'))
  suite.addTest(TestABInertiadLeftOperatorsTest('test'))
  unittest.TextTestRunner(verbosity=2).run(suite)
