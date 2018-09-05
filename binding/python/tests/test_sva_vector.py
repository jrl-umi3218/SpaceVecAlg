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

class TestMotionVecd(unittest.TestCase):
  def test(self):
    w = eigen.Vector3d.Random()
    v = eigen.Vector3d.Random()

    vec = sva.MotionVecd(w, v)
    m = vec.vector()

    self.assertEqual(w, vec.angular())
    self.assertEqual(v, vec.linear())

    self.assertEqual(m, eigen.Vector6d(list(w) + list(v)))

    self.assertEqual((5.*vec).vector(), 5.*m)
    self.assertEqual((vec*5.).vector(), 5.*m)
    self.assertEqual((vec/5.).vector(), m/5)
    self.assertEqual((-vec).vector(), -m)

    w2 = eigen.Vector3d.Random()
    v2 = eigen.Vector3d.Random()
    vec2 = sva.MotionVecd(w2, v2)
    m2 = vec2.vector()

    self.assertEqual((vec + vec2).vector(), (m + m2))
    self.assertEqual((vec - vec2).vector(), (m - m2))

    vec_pluseq = sva.MotionVecd(vec)
    vec_pluseq += vec2
    self.assertEqual(vec_pluseq, vec + vec2)

    vec_minuseq = sva.MotionVecd(vec)
    vec_minuseq -= vec2
    self.assertEqual(vec_minuseq, vec - vec2)

    self.assertEqual(vec, vec)
    self.assertNotEqual(vec, -vec)

    self.assertTrue(vec != (-vec))
    self.assertTrue(not(vec != vec))

    z = sva.MotionVecd([0, 0, 0], [0, 0, 0])
    self.assertEqual(sva.MotionVecd.Zero(), z)

class TestForceVecd(unittest.TestCase):
  def test(self):
    n = eigen.Vector3d.Random()
    f = eigen.Vector3d.Random()

    vec = sva.ForceVecd(n, f)
    m = vec.vector()

    self.assertEqual(n, vec.couple())
    self.assertEqual(f, vec.force())

    self.assertEqual(m, eigen.Vector6d(list(n) + list(f)))

    self.assertEqual((5.*vec).vector(), 5.*m)
    self.assertEqual((vec*5.).vector(), 5.*m)
    self.assertEqual((vec/5.).vector(), m/5)
    self.assertEqual((-vec).vector(), -m)

    n2 = eigen.Vector3d.Random()
    f2 = eigen.Vector3d.Random()
    vec2 = sva.ForceVecd(n2, f2)
    m2 = vec2.vector()

    self.assertEqual((vec + vec2).vector(), (m + m2))
    self.assertEqual((vec - vec2).vector(), (m - m2))

    vec_pluseq = sva.ForceVecd(vec)
    vec_pluseq += vec2
    self.assertEqual(vec_pluseq, vec + vec2)

    vec_minuseq = sva.ForceVecd(vec)
    vec_minuseq -= vec2
    self.assertEqual(vec_minuseq, vec - vec2)

    self.assertEqual(vec, vec)
    self.assertNotEqual(vec, -vec)

    self.assertTrue(vec != (-vec))
    self.assertTrue(not(vec != vec))

    z = sva.ForceVecd([0, 0, 0], [0, 0, 0])
    self.assertEqual(sva.ForceVecd.Zero(), z)

class TestMotionVecdLeftOperatorsTest(unittest.TestCase):
  def test(self):
    w = eigen.Vector3d.Random()*100
    v = eigen.Vector3d.Random()*100
    n = eigen.Vector3d.Random()*100
    f = eigen.Vector3d.Random()*100

    mVec = sva.MotionVecd(w, v)
    fVec = sva.ForceVecd(n, f)

    mm = mVec.vector()
    mf = fVec.vector()

    self.assertAlmostEqual(mVec.dot(fVec) - (mm.transpose()*mf)[0], 0, delta = TOL)

    w2 = eigen.Vector3d.Random()*100
    v2 = eigen.Vector3d.Random()*100

    mVec2 = sva.MotionVecd(w2, v2)
    mm2 = mVec2.vector()

    crossM = mVec.cross(mVec2)
    self.assertTrue(isinstance(crossM, sva.MotionVecd))
    self.assertAlmostEqual((crossM.vector() - sva.vector6ToCrossMatrix(mm)*mm2).norm(), 0, delta = TOL)

    crossF = mVec.crossDual(fVec)
    self.assertTrue(isinstance(crossF, sva.ForceVecd))
    self.assertAlmostEqual((crossF.vector() - sva.vector6ToCrossDualMatrix(mm)*mf).norm(), 0, delta = TOL)

class TestImpedanceVecd(unittest.TestCase):
  def test(self):
    w = eigen.Vector3d.Random()
    v = eigen.Vector3d.Random()
    vec = sva.ImpedanceVecd(w, v)
    z = vec.vector()

    self.assertEqual(w, vec.angular())
    self.assertEqual(v, vec.linear())
    self.assertEqual(z, eigen.Vector6d(w.x(), w.y(), w.z(), v.x(), v.y(), v.z()))

    self.assertEqual((5.*vec).vector(), 5.*z)
    self.assertEqual((vec*5.).vector(), 5.*z)
    self.assertEqual((vec/5.).vector(), z/5.)

    vec *= 5.
    self.assertEqual(vec.vector(), 5.*z)
    vec /= 5.
    [ self.assertAlmostEqual(x, 0, delta = TOL) for x in vec.vector() - z ]

    w2 = eigen.Vector3d.Random()
    v2 = eigen.Vector3d.Random()
    vec2 = sva.ImpedanceVecd(w2, v2)
    z2 = vec2.vector()

    self.assertEqual((vec + vec2).vector(), z + z2)

    vec_pluseq = sva.ImpedanceVecd(vec)
    self.assertEqual(vec_pluseq, vec)
    vec_pluseq += vec2
    self.assertEqual(vec_pluseq, vec + vec2)

    self.assertEqual(vec, vec)
    self.assertTrue(not (vec != vec))

    w = eigen.Vector3d.Random()
    v = eigen.Vector3d.Random()
    mv = sva.MotionVecd(eigen.Vector3d.Random(), eigen.Vector3d.Random())

    fv = vec * mv
    self.assertTrue(isinstance(fv, sva.ForceVecd))
    res = eigen.Vector6d([ x*y for x,y in zip(vec.vector(), mv.vector()) ])
    self.assertEqual(res, fv.vector())

    fv2 = mv * vec
    self.assertEqual(fv, fv2)

    hv = sva.ImpedanceVecd(11., 42.)
    self.assertEqual(hv.angular(), eigen.Vector3d(11., 11., 11.))
    self.assertEqual(hv.linear(), eigen.Vector3d(42., 42., 42.))

    z = sva.ImpedanceVecd(0., 0.)
    self.assertEqual(sva.ImpedanceVecd.Zero(), z)

class TestAdmittanceVecd(unittest.TestCase):
  def test(self):
    w = eigen.Vector3d.Random()
    v = eigen.Vector3d.Random()
    vec = sva.AdmittanceVecd(w, v)
    z = vec.vector()

    self.assertEqual(w, vec.angular())
    self.assertEqual(v, vec.linear())
    self.assertEqual(z, eigen.Vector6d(w.x(), w.y(), w.z(), v.x(), v.y(), v.z()))

    self.assertEqual((5.*vec).vector(), 5.*z)
    self.assertEqual((vec*5.).vector(), 5.*z)
    self.assertEqual((vec/5.).vector(), z/5.)

    vec *= 5.
    self.assertEqual(vec.vector(), 5.*z)
    vec /= 5.
    [ self.assertAlmostEqual(x, 0, delta = TOL) for x in vec.vector() - z ]

    w2 = eigen.Vector3d.Random()
    v2 = eigen.Vector3d.Random()
    vec2 = sva.AdmittanceVecd(w2, v2)
    z2 = vec2.vector()

    self.assertEqual((vec + vec2).vector(), z + z2)

    vec_pluseq = sva.AdmittanceVecd(vec)
    self.assertEqual(vec_pluseq, vec)
    vec_pluseq += vec2
    self.assertEqual(vec_pluseq, vec + vec2)

    self.assertEqual(vec, vec)
    self.assertTrue(not (vec != vec))

    w = eigen.Vector3d.Random()
    v = eigen.Vector3d.Random()
    fv = sva.ForceVecd(eigen.Vector3d.Random(), eigen.Vector3d.Random())

    mv = vec * fv
    self.assertTrue(isinstance(mv, sva.MotionVecd))
    res = eigen.Vector6d([ x*y for x,y in zip(vec.vector(), fv.vector()) ])
    self.assertEqual(res, mv.vector())

    mv2 = fv * vec
    self.assertEqual(mv, mv2)

    hv = sva.AdmittanceVecd(11., 42.)
    self.assertEqual(hv.angular(), eigen.Vector3d(11., 11., 11.))
    self.assertEqual(hv.linear(), eigen.Vector3d(42., 42., 42.))

    z = sva.AdmittanceVecd(0., 0.)
    self.assertEqual(sva.AdmittanceVecd.Zero(), z)


if __name__ == "__main__":
  suite = unittest.TestSuite()
  suite.addTest(TestMotionVecd('test'))
  suite.addTest(TestForceVecd('test'))
  suite.addTest(TestMotionVecdLeftOperatorsTest('test'))
  suite.addTest(TestImpedanceVecd('test'))
  suite.addTest(TestAdmittanceVecd('test'))
  unittest.TextTestRunner(verbosity=2).run(suite)
