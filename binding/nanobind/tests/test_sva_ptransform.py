import pytest
import numpy as np
import nanoeigenpy as eigen
import sva
import sys

TOL = 1e-5

def test_rotation_matrix():
    theta = np.random.randn() * 10

    # X axis
    aa_x = eigen.AngleAxis(-theta, np.array([1.0, 0.0, 0.0]))
    rot_x = aa_x.matrix()
    assert np.allclose(sva.RotX(theta), rot_x, atol=TOL)

    # Y axis
    aa_y = eigen.AngleAxis(-theta, np.array([0.0, 1.0, 0.0]))
    rot_y = aa_y.matrix()
    assert np.allclose(sva.RotY(theta), rot_y, atol=TOL)

    # Z axis
    aa_z = eigen.AngleAxis(-theta, np.array([0.0, 0.0, 1.0]))
    rot_z = aa_z.matrix()
    assert np.allclose(sva.RotZ(theta), rot_z, atol=TOL)

def test_ptransformd():
    # Rotation matrix and quaternion from AngleAxis
    Em = eigen.AngleAxis(np.pi/2, np.array([1.0, 0.0, 0.0])).inverse().matrix()
    Eq = eigen.Quaternion(eigen.AngleAxis(np.pi/2, np.array([1.0, 0.0, 0.0])).inverse())
    r = np.random.randn(3) * 100

    pt1 = sva.PTransformd.Identity()
    assert np.allclose(pt1.rotation(), np.eye(3), atol=TOL)
    assert np.allclose(pt1.translation(), np.zeros(3), atol=TOL)

    pt2 = sva.PTransformd(Em, r)
    assert np.allclose(pt2.rotation(), Em, atol=TOL)
    assert np.allclose(pt2.translation(), r, atol=TOL)

    pt3 = sva.PTransformd(Eq, r)
    assert np.allclose(pt3.rotation(), Eq.matrix(), atol=TOL)
    assert np.allclose(pt3.translation(), r, atol=TOL)

    pt4 = sva.PTransformd(Eq)
    assert np.allclose(pt4.rotation(), Eq.matrix(), atol=TOL)
    assert np.allclose(pt4.translation(), np.zeros(3), atol=TOL)

    pt5 = sva.PTransformd(Em)
    assert np.allclose(pt5.rotation(), Em, atol=TOL)
    assert np.allclose(pt5.translation(), np.zeros(3), atol=TOL)

    pt6 = sva.PTransformd(r)
    assert np.allclose(pt6.rotation(), np.eye(3), atol=TOL)
    assert np.allclose(pt6.translation(), r, atol=TOL)

    pttmp = sva.PTransformd(
        eigen.AngleAxis(np.pi/4, np.array([0.0, 1.0, 0.0])).matrix(),
        np.random.randn(3) * 100
    )

    pt7 = pt2 * pttmp
    ptm = pt2.matrix() @ pttmp.matrix()
    assert np.allclose(pt7.matrix(), ptm, atol=TOL)

    pt8 = pt2.inv()
    assert np.allclose(pt8.matrix(), np.linalg.inv(pt2.matrix()), atol=TOL)

    assert pt2 == pt2
    assert pt2 != pt8

    assert pt2 != pt8
    assert not (pt2 != pt2)

def test_ptransformd_left_operator():
    Eq = eigen.AngleAxis(np.pi/2, np.array([1.0, 0.0, 0.0]))
    r = np.random.randn(3) * 100
    pt = sva.PTransformd(Eq.matrix(), r)
    ptInv = pt.inv()
    pt6d = pt.matrix()
    ptInv6d = ptInv.matrix()
    ptDual6d = pt.dualMatrix()

    M = np.array([[1,2,3],[2,1,4],[3,4,1]], dtype=float)
    H = np.random.randn(3,3) * 100
    I = np.array([[1,2,3],[2,1,4],[3,4,1]], dtype=float)
    ab = sva.ABInertiad(M, H, I)
    ab6d = ab.matrix()

    mass = 1.
    h = np.random.randn(3) * 100
    rb = sva.RBInertiad(mass, h, I)
    rb6d = rb.matrix()

    w = np.random.randn(3) * 100
    v = np.random.randn(3) * 100
    n = np.random.randn(3) * 100
    f = np.random.randn(3) * 100

    mVec = sva.MotionVecd(w, v)
    mVec6d = mVec.vector()

    fVec = sva.ForceVecd(n, f)
    fVec6d = fVec.vector()

    mvRes1 = pt * mVec
    mvRes16d = pt6d @ mVec6d
    assert np.linalg.norm(mvRes1.vector() - mvRes16d) < TOL
    assert np.linalg.norm(mvRes1.angular() - pt.angularMul(mVec)) < TOL
    assert np.linalg.norm(mvRes1.linear() - pt.linearMul(mVec)) < TOL

    mvRes2 = pt.invMul(mVec)
    mvRes26d = ptInv6d @ mVec6d
    assert np.linalg.norm(mvRes2.vector() - mvRes26d) < TOL
    assert np.linalg.norm(mvRes2.angular() - pt.angularInvMul(mVec)) < TOL
    assert np.linalg.norm(mvRes2.linear() - pt.linearInvMul(mVec)) < TOL

    fvRes1 = pt.dualMul(fVec)
    fvRes16d = ptDual6d @ fVec6d
    assert np.linalg.norm(fvRes1.vector() - fvRes16d) < TOL
    assert np.linalg.norm(fvRes1.couple() - pt.coupleDualMul(fVec)) < TOL
    assert np.linalg.norm(fvRes1.force() - pt.forceDualMul(fVec)) < TOL

    fvRes2 = pt.transMul(fVec)
    fvRes26d = pt6d.T @ fVec6d
    assert np.linalg.norm(fvRes2.vector() - fvRes26d) < TOL
    assert np.linalg.norm(fvRes2.couple() - pt.coupleTransMul(fVec)) < TOL
    assert np.linalg.norm(fvRes2.force() - pt.forceTransMul(fVec)) < TOL

    rbRes1 = pt.dualMul(rb)
    rbRes16d = ptDual6d @ rb6d @ ptInv6d
    assert np.linalg.norm(rbRes1.matrix() - rbRes16d) < TOL

    rbRes2 = pt.transMul(rb)
    rbRes26d = pt6d.T @ rb6d @ pt6d
    assert np.linalg.norm(rbRes2.matrix() - rbRes26d) < TOL

    abRes1 = pt.dualMul(ab)
    abRes16d = ptDual6d @ ab6d @ ptInv6d
    assert np.linalg.norm(abRes1.matrix() - abRes16d) < TOL

    abRes2 = pt.transMul(ab)
    abRes26d = pt6d.T @ ab6d @ pt6d
    assert np.linalg.norm(abRes2.matrix() - abRes26d) < TOL



def test_euler_angle():
    # X axis
    res = sva.rotationError(np.eye(3), sva.RotX(np.pi/2))
    expected = np.array([np.pi/2, 0, 0])
    assert np.isclose(np.linalg.norm(res - expected), 0, atol=TOL)

    # Y axis
    res = sva.rotationError(np.eye(3), sva.RotY(np.pi/2))
    expected = np.array([0, np.pi/2, 0])
    assert np.isclose(np.linalg.norm(res - expected), 0, atol=TOL)

    # Z axis
    res = sva.rotationError(np.eye(3), sva.RotZ(np.pi/2))
    expected = np.array([0, 0, np.pi/2])
    assert np.isclose(np.linalg.norm(res - expected), 0, atol=TOL)

    # Z axis, relative rotation
    res = sva.rotationError(sva.RotZ(np.pi/4), sva.RotZ(np.pi/2))
    expected = np.array([0, 0, np.pi/4])
    assert np.isclose(np.linalg.norm(res - expected), 0, atol=TOL)

def test_interpolate():
    from_ = sva.PTransformd(np.eye(3), np.zeros(3))
    to = sva.PTransformd(sva.RotZ(np.pi), np.array([1, 2, -3]))
    res = sva.interpolate(from_, to, 0.5)
    expected_rot = sva.RotZ(np.pi/2.)
    expected_trans = np.array([0.5, 1., -1.5])
    assert np.allclose(res.rotation(), expected_rot, atol=TOL)
    assert np.allclose(res.translation(), expected_trans, atol=TOL)

def test_transform_error():
    X_a_b = sva.PTransformd(eigen.Quaternion(np.random.randn(4)).normalized(), np.random.randn(3))
    X_a_c = sva.PTransformd(eigen.Quaternion(np.random.randn(4)).normalized(), np.random.randn(3))

    V_a_b = sva.transformVelocity(X_a_b)
    V_a_c = sva.transformVelocity(X_a_c)

    assert np.isclose(np.linalg.norm(V_a_b.angular() - sva.rotationVelocity(X_a_b.rotation())), 0, atol=TOL)
    assert np.isclose(np.linalg.norm(V_a_b.linear() - X_a_b.translation()), 0, atol=TOL)
    assert np.isclose(np.linalg.norm(V_a_c.angular() - sva.rotationVelocity(X_a_c.rotation())), 0, atol=TOL)
    assert np.isclose(np.linalg.norm(V_a_c.linear() - X_a_c.translation()), 0, atol=TOL)

    V_b_c_a = sva.transformError(X_a_b, X_a_c)
    w_b_c_a = sva.rotationError(X_a_b.rotation(), X_a_c.rotation())
    v_b_c_a = X_a_c.translation() - X_a_b.translation()

    assert np.isclose(np.linalg.norm(V_b_c_a.angular() - w_b_c_a), 0, atol=TOL)
    assert np.isclose(np.linalg.norm(V_b_c_a.linear() - v_b_c_a), 0, atol=TOL)

def test_sinc_inv():
    dummy_sinc_inv = lambda x : x/np.sin(x)
    eps = sys.float_info.epsilon

    t = -1.
    nrIter = 333
    for i in range(nrIter):
        assert np.isclose(dummy_sinc_inv(t), sva.sinc_inv(t), atol=TOL)
        t += 2./nrIter

    # Test for nan and warning at zero
    with pytest.warns(RuntimeWarning):
        assert np.isnan(dummy_sinc_inv(0))
    assert sva.sinc_inv(0) == 1

    assert np.isclose(dummy_sinc_inv(eps), sva.sinc_inv(eps), atol=TOL)
    assert np.isclose(dummy_sinc_inv(np.sqrt(eps)), sva.sinc_inv(np.sqrt(eps)), atol=TOL)
    assert np.isclose(dummy_sinc_inv(np.sqrt(np.sqrt(eps))), sva.sinc_inv(np.sqrt(np.sqrt(eps))), atol=TOL)
