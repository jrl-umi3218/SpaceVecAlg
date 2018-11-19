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

class TestRotationMatrix(unittest.TestCase):
  def test(self):
    theta2d = eigen.Vector2d.Random()*10
    theta = theta2d[0]

    self.assertAlmostEqual((sva.RotX(theta) - eigen.AngleAxisd(-theta, eigen.Vector3d.UnitX()).matrix()).norm(), 0, delta = TOL)
    self.assertAlmostEqual((sva.RotY(theta) - eigen.AngleAxisd(-theta, eigen.Vector3d.UnitY()).matrix()).norm(), 0, delta = TOL)
    self.assertAlmostEqual((sva.RotZ(theta) - eigen.AngleAxisd(-theta, eigen.Vector3d.UnitZ()).matrix()).norm(), 0, delta = TOL)

class TestPTransformd(unittest.TestCase):
  def test(self):
    Em = eigen.AngleAxisd(np.pi/2, eigen.Vector3d.UnitX()).inverse().matrix()
    Eq = eigen.Quaterniond(eigen.AngleAxisd(np.pi/2, eigen.Vector3d.UnitX()).inverse())
    r = eigen.Vector3d.Random()*100

    pt1 = sva.PTransformd.Identity()
    self.assertEqual(pt1.rotation(), eigen.Matrix3d.Identity())
    self.assertEqual(pt1.translation(), eigen.Vector3d.Zero())

    pt2 = sva.PTransformd(Em, r)
    self.assertEqual(pt2.rotation(), Em)
    self.assertEqual(pt2.translation(), r)

    pt3 = sva.PTransformd(Eq, r)
    self.assertEqual(pt3.rotation(), Eq.toRotationMatrix())
    self.assertEqual(pt3.translation(), r)

    pt4 = sva.PTransformd(Eq)
    self.assertEqual(pt4.rotation(), Eq.toRotationMatrix())
    self.assertEqual(pt4.translation(), eigen.Vector3d.Zero())

    pt5 = sva.PTransformd(Em)
    self.assertEqual(pt5.rotation(), Em)
    self.assertEqual(pt5.translation(), eigen.Vector3d.Zero())

    pt6 = sva.PTransformd(r)
    self.assertEqual(pt6.rotation(), eigen.Matrix3d.Identity())
    self.assertEqual(pt6.translation(), r)

    pttmp = sva.PTransformd(eigen.AngleAxisd(np.pi/4, eigen.Vector3d.UnitY()).matrix(), eigen.Vector3d.Random()*100.)

    pt7 = pt2*pttmp
    ptm = pt2.matrix() * pttmp.matrix()
    self.assertAlmostEqual((pt7.matrix() - ptm).norm(), 0, delta = TOL)

    pt8 = pt2.inv()
    self.assertAlmostEqual((pt8.matrix() - pt2.matrix().inverse()).norm(), 0, delta = TOL)

    self.assertEqual(pt2, pt2)
    self.assertNotEqual(pt2, pt8)

    self.assertTrue(pt2 != pt8)
    self.assertTrue(not(pt2 != pt2))

class TestPTransformdLeftOperator(unittest.TestCase):
  def test(self):
    Eq = eigen.AngleAxisd(np.pi/2, eigen.Vector3d.UnitX())
    r = eigen.Vector3d.Random() * 100
    pt = sva.PTransformd(Eq, r)
    ptInv = pt.inv()
    pt6d = pt.matrix()
    ptInv6d = ptInv.matrix()
    ptDual6d = pt.dualMatrix()

    M = eigen.Matrix3d([[1,2,3],[2,1,4],[3,4,1]])
    H = eigen.Matrix3d.Random() * 100
    I = eigen.Matrix3d([[1,2,3],[2,1,4],[3,4,1]])
    ab = sva.ABInertiad(M, H, I)
    ab6d = ab.matrix()

    mass = 1.
    h = eigen.Vector3d.Random()*100
    rb = sva.RBInertiad(mass, h, I)
    rb6d = rb.matrix()

    w = eigen.Vector3d.Random()*100
    v = eigen.Vector3d.Random()*100
    n = eigen.Vector3d.Random()*100
    f = eigen.Vector3d.Random()*100

    mVec = sva.MotionVecd(w, v)
    mVec6d = mVec.vector()

    fVec = sva.ForceVecd(n, f)
    fVec6d = fVec.vector()

    mvRes1 = pt*mVec
    mvRes16d = pt6d*mVec6d

    self.assertAlmostEqual((mvRes1.vector() - mvRes16d).norm(), 0, delta = TOL)
    self.assertAlmostEqual((mvRes1.angular() - pt.angularMul(mVec)).norm(), 0, delta = TOL)
    self.assertAlmostEqual((mvRes1.linear() - pt.linearMul(mVec)).norm(), 0, delta = TOL)

    # Note: the C++ version would test the vectorized version here but this is not supported by the Python bindings

    mvRes2 = pt.invMul(mVec)
    mvRes26d = ptInv6d*mVec6d
    self.assertAlmostEqual((mvRes2.vector() - mvRes26d).norm(), 0, delta = TOL)
    self.assertAlmostEqual((mvRes2.angular() - pt.angularInvMul(mVec)).norm(), 0, delta = TOL)
    self.assertAlmostEqual((mvRes2.linear() - pt.linearInvMul(mVec)).norm(), 0, delta = TOL)

    # Same note about the vectorized version

    fvRes1 = pt.dualMul(fVec)
    fvRes16d = ptDual6d*fVec6d

    self.assertAlmostEqual((fvRes1.vector() - fvRes16d).norm(), 0, delta = TOL)
    self.assertAlmostEqual((fvRes1.couple() - pt.coupleDualMul(fVec)).norm(), 0, delta = TOL)
    self.assertAlmostEqual((fvRes1.force() - pt.forceDualMul(fVec)).norm(), 0, delta = TOL)

    # Same note about the vectorized version

    fvRes2 = pt.transMul(fVec)
    fvRes26d = pt6d.transpose() * fVec6d

    self.assertAlmostEqual((fvRes2.vector() - fvRes26d).norm(), 0, delta= TOL)
    self.assertAlmostEqual((fvRes2.couple() - pt.coupleTransMul(fVec)).norm(), 0, delta = TOL)
    self.assertAlmostEqual((fvRes2.force() - pt.forceTransMul(fVec)).norm(), 0, delta = TOL)

    # Same note about the vectorized version

    rbRes1 = pt.dualMul(rb)
    rbRes16d = ptDual6d * rb6d * ptInv6d
    self.assertAlmostEqual((rbRes1.matrix() - rbRes16d).norm(), 0, delta = TOL)

    rbRes2 = pt.transMul(rb);
    rbRes26d = pt6d.transpose()*rb6d*pt6d
    self.assertAlmostEqual((rbRes2.matrix() - rbRes26d).norm(), 0, delta = TOL);

    abRes1 = pt.dualMul(ab);
    abRes16d = ptDual6d*ab6d*ptInv6d
    self.assertAlmostEqual((abRes1.matrix() - abRes16d).norm(), 0, delta = TOL)

    abRes2 = pt.transMul(ab);
    abRes26d = pt6d.transpose()*ab6d*pt6d
    self.assertAlmostEqual((abRes2.matrix() - abRes26d).norm(), 0, delta = TOL)

class TestEulerAngle(unittest.TestCase):
  def test(self):
    res = eigen.Vector3d()

    res = sva.rotationError(eigen.Matrix3d.Identity(), sva.RotX(np.pi/2))
    self.assertAlmostEqual((res - eigen.Vector3d(np.pi/2, 0, 0)).norm(), 0, delta = TOL)

    res = sva.rotationError(eigen.Matrix3d.Identity(), sva.RotY(np.pi/2))
    self.assertAlmostEqual((res - eigen.Vector3d(0, np.pi/2, 0)).norm(), 0, delta = TOL)

    res = sva.rotationError(eigen.Matrix3d.Identity(), sva.RotZ(np.pi/2))
    self.assertAlmostEqual((res - eigen.Vector3d(0, 0, np.pi/2)).norm(), 0, delta = TOL)

    res = sva.rotationError(sva.RotZ(np.pi/4), sva.RotZ(np.pi/2))
    self.assertAlmostEqual((res - eigen.Vector3d(0, 0, np.pi/4)).norm(), 0, delta = TOL)

class TestInterpolate(unittest.TestCase):
  def test(self):
    from_ = sva.PTransformd(eigen.Matrix3d.Identity(), eigen.Vector3d.Zero())
    to = sva.PTransformd(eigen.AngleAxisd(np.pi, eigen.Vector3d.UnitZ()), eigen.Vector3d(1, 2, -3))

    res = sva.interpolate(from_, to, 0.5)
    self.assertAlmostEqual((res.rotation() - eigen.AngleAxisd(np.pi/2, eigen.Vector3d.UnitZ()).matrix()).norm(), 0, delta = TOL)
    self.assertAlmostEqual((res.translation() - eigen.Vector3d(0.5, 1., -1.5)).norm(), 0,  delta = TOL)

class TestTransformError(unittest.TestCase):
  def test(self):
    X_a_b = sva.PTransformd(eigen.Quaterniond(eigen.Vector4d.Random()).normalized(), eigen.Vector3d.Random())
    X_a_c = sva.PTransformd(eigen.Quaterniond(eigen.Vector4d.Random()).normalized(), eigen.Vector3d.Random())

    V_a_b = sva.transformVelocity(X_a_b)
    V_a_c = sva.transformVelocity(X_a_c)

    self.assertAlmostEqual((V_a_b.angular() - sva.rotationVelocity(X_a_b.rotation())).norm(), 0, delta = TOL)
    self.assertAlmostEqual((V_a_b.linear() - X_a_b.translation()).norm(), 0, delta = TOL)
    self.assertAlmostEqual((V_a_c.angular() - sva.rotationVelocity(X_a_c.rotation())).norm(), 0, delta = TOL)
    self.assertAlmostEqual((V_a_c.linear() - X_a_c.translation()).norm(), 0, delta = TOL)

    V_b_c_a = sva.transformError(X_a_b, X_a_c)
    w_b_c_a = sva.rotationError(X_a_b.rotation(), X_a_c.rotation())
    v_b_c_a = X_a_c.translation() - X_a_b.translation()

    self.assertAlmostEqual((V_b_c_a.angular() - w_b_c_a).norm(), 0, delta = TOL)
    self.assertAlmostEqual((V_b_c_a.linear() - v_b_c_a).norm(), 0, delta = TOL)

class Testsinc_inv(unittest.TestCase):
  def test(self):
    dummy_sinc_inv = lambda x : x/np.sin(x)
    eps = sys.float_info.epsilon

    t = -1.
    nrIter = 333
    for i in range(nrIter):
      self.assertEqual(dummy_sinc_inv(t), sva.sinc_inv(t))
      t += 2./nrIter

    self.assertTrue(np.isnan(dummy_sinc_inv(0)))
    self.assertEqual(sva.sinc_inv(0), 1)

    self.assertEqual(dummy_sinc_inv(eps), sva.sinc_inv(eps))
    self.assertEqual(dummy_sinc_inv(np.sqrt(eps)),
                     sva.sinc_inv(np.sqrt(eps)))
    self.assertEqual(dummy_sinc_inv(np.sqrt(np.sqrt(eps))),
                     sva.sinc_inv(np.sqrt(np.sqrt(eps))))

if __name__ == "__main__":
  suite = unittest.TestSuite()
  suite.addTest(TestRotationMatrix('test'))
  suite.addTest(TestPTransformd('test'))
  suite.addTest(TestPTransformdLeftOperator('test'))
  suite.addTest(TestEulerAngle('test'))
  suite.addTest(TestInterpolate('test'))
  suite.addTest(TestTransformError('test'))
  suite.addTest(Testsinc_inv('test'))
  unittest.TextTestRunner(verbosity=2).run(suite)
