import pytest
import numpy as np
import sva

TOL = 1e-5

def is_upper_null(m):
    # m is a 3x3 numpy array
    return m[0,1] == m[0,2] == m[1,2] == 0

def test_rbinertiad():
    mass = 1.
    I = np.array([[1, 2, 3],[2, 1, 4], [3, 4, 1]], dtype=float)
    h = np.random.randn(3)

    rb2 = sva.RBInertiad(mass, h, I)
    assert rb2.mass() == mass
    assert np.allclose(rb2.momentum(), h, atol=TOL)
    assert np.allclose(rb2.inertia(), I, atol=TOL)
    assert is_upper_null(rb2.lowerTriangularInertia())

    rb3 = sva.RBInertiad(mass, h, rb2.lowerTriangularInertia())
    assert rb3.mass() == mass
    assert np.allclose(rb3.momentum(), h, atol=TOL)
    assert np.allclose(rb3.inertia(), I, atol=TOL)
    assert is_upper_null(rb3.lowerTriangularInertia())

    rb4 = rb2 + rb3
    assert rb4.mass() == mass + mass
    assert np.allclose(rb4.momentum(), h + h, atol=TOL)
    assert np.allclose(rb4.inertia(), I + I, atol=TOL)
    assert is_upper_null(rb4.lowerTriangularInertia())

    rb5 = 2. * rb2
    assert rb5.mass() == 2. * mass
    assert np.allclose(rb5.momentum(), 2. * h, atol=TOL)
    assert np.allclose(rb5.inertia(), 2. * I, atol=TOL)
    assert is_upper_null(rb5.lowerTriangularInertia())

    rb6 = rb2 * 2.
    assert rb6.mass() == 2. * mass
    assert np.allclose(rb6.momentum(), 2. * h, atol=TOL)
    assert np.allclose(rb6.inertia(), 2. * I, atol=TOL)
    assert is_upper_null(rb6.lowerTriangularInertia())

    rb7 = rb2 - rb3
    assert rb7.mass() == mass - mass
    assert np.allclose(rb7.momentum(), h - h, atol=TOL)
    assert np.allclose(rb7.inertia(), I - I, atol=TOL)
    assert is_upper_null(rb7.lowerTriangularInertia())

    rb8 = -rb2
    assert np.allclose(rb8.mass(), -mass, atol=TOL)
    assert np.allclose(rb8.momentum(), -h, atol=TOL)
    assert np.allclose(rb8.inertia(), -I, atol=TOL)
    assert is_upper_null(rb8.lowerTriangularInertia())

    rb9 = sva.RBInertiad(rb2)
    rb9 += rb3
    assert np.allclose(rb9.mass(), rb2.mass() + rb3.mass(), atol=TOL)
    assert np.allclose(rb9.momentum(), rb2.momentum() + rb3.momentum(), atol=TOL)
    assert np.allclose(rb9.inertia(), rb2.inertia() + rb3.inertia(), atol=TOL)
    assert is_upper_null(rb9.lowerTriangularInertia())

    rb10 = sva.RBInertiad(rb2)
    rb10 -= rb3
    assert np.allclose(rb10.mass(), rb2.mass() - rb3.mass(), atol=TOL)
    assert np.allclose(rb10.momentum(), rb2.momentum() - rb3.momentum(), atol=TOL)
    assert np.allclose(rb10.inertia(), rb2.inertia() - rb3.inertia(), atol=TOL)
    assert is_upper_null(rb10.lowerTriangularInertia())

    assert rb2 == rb2
    assert rb2 != rb6
    assert rb2 != rb6
    assert not (rb2 != rb2)

def test_abinertiad():
    M = np.array([[1, 2, 3],[2, 1, 4], [3, 4, 1]], dtype=float)
    H = np.random.randn(3,3)
    I = np.array([[1, 2, 3],[2, 1, 4], [3, 4, 1]], dtype=float)

    ab1 = sva.ABInertiad(M, H, I)
    assert np.allclose(ab1.massMatrix(), M, atol=TOL)
    assert is_upper_null(ab1.lowerTriangularMassMatrix())
    assert np.allclose(ab1.gInertia(), H, atol=TOL)
    assert np.allclose(ab1.inertia(), I, atol=TOL)
    assert is_upper_null(ab1.lowerTriangularInertia())

    ab2 = sva.ABInertiad(ab1.lowerTriangularMassMatrix(), H, ab1.lowerTriangularInertia())
    assert np.allclose(ab2.massMatrix(), M, atol=TOL)
    assert is_upper_null(ab2.lowerTriangularMassMatrix())
    assert np.allclose(ab2.gInertia(), H, atol=TOL)
    assert np.allclose(ab2.inertia(), I, atol=TOL)
    assert is_upper_null(ab2.lowerTriangularInertia())

    ab3 = ab1 + ab2
    assert np.allclose(ab3.massMatrix(), M + M, atol=TOL)
    assert is_upper_null(ab3.lowerTriangularMassMatrix())
    assert np.allclose(ab3.gInertia(), H + H, atol=TOL)
    assert np.allclose(ab3.inertia(), I + I, atol=TOL)
    assert is_upper_null(ab3.lowerTriangularInertia())

    ab4 = 2. * ab2
    assert np.allclose(ab4.massMatrix(), 2. * M, atol=TOL)
    assert is_upper_null(ab4.lowerTriangularMassMatrix())
    assert np.allclose(ab4.gInertia(), 2. * H, atol=TOL)
    assert np.allclose(ab4.inertia(), 2. * I, atol=TOL)
    assert is_upper_null(ab4.lowerTriangularInertia())

    ab5 = ab2 * 2.
    assert np.allclose(ab5.massMatrix(), 2. * M, atol=TOL)
    assert is_upper_null(ab5.lowerTriangularMassMatrix())
    assert np.allclose(ab5.gInertia(), 2. * H, atol=TOL)
    assert np.allclose(ab5.inertia(), 2. * I, atol=TOL)
    assert is_upper_null(ab5.lowerTriangularInertia())

    ab6 = ab1 - ab2
    assert np.allclose(ab6.massMatrix(), M - M, atol=TOL)
    assert is_upper_null(ab6.lowerTriangularMassMatrix())
    assert np.allclose(ab6.gInertia(), H - H, atol=TOL)
    assert np.allclose(ab6.inertia(), I - I, atol=TOL)
    assert is_upper_null(ab6.lowerTriangularInertia())

    ab7 = -ab1
    assert np.allclose(ab7.massMatrix(), -M, atol=TOL)
    assert np.allclose(ab7.gInertia(), -H, atol=TOL)
    assert np.allclose(ab7.inertia(), -I, atol=TOL)
    assert is_upper_null(ab7.lowerTriangularMassMatrix())
    assert is_upper_null(ab7.lowerTriangularInertia())

    ab8 = sva.ABInertiad(ab1)
    ab8 += ab2
    assert np.allclose(ab8.massMatrix(), ab1.massMatrix() + ab2.massMatrix(), atol=TOL)
    assert np.allclose(ab8.gInertia(), ab1.gInertia() + ab2.gInertia(), atol=TOL)
    assert np.allclose(ab8.inertia(), ab1.inertia() + ab2.inertia(), atol=TOL)
    assert is_upper_null(ab8.lowerTriangularMassMatrix())
    assert is_upper_null(ab8.lowerTriangularInertia())

    ab9 = sva.ABInertiad(ab1)
    ab9 -= ab2
    assert np.allclose(ab9.massMatrix(), ab1.massMatrix() - ab2.massMatrix(), atol=TOL)
    assert np.allclose(ab9.gInertia(), ab1.gInertia() - ab2.gInertia(), atol=TOL)
    assert np.allclose(ab9.inertia(), ab1.inertia() - ab2.inertia(), atol=TOL)
    assert is_upper_null(ab9.lowerTriangularMassMatrix())
    assert is_upper_null(ab9.lowerTriangularInertia())

    assert ab2 == ab2
    assert ab2 != ab5
    assert ab2 != ab5
    assert not (ab2 != ab2)

def test_rbinertiad_left_operators():
    mass = 1.
    I = np.array([[1, 2, 3],[2, 1, 4], [3, 4, 1]], dtype=float)
    h = np.random.randn(3)
    rb = sva.RBInertiad(mass, h, I)
    rb6d = rb.matrix()

    w = np.random.randn(3)
    v = np.random.randn(3)
    mVec = sva.MotionVecd(w, v)
    mVec6d = mVec.vector()

    fVec = rb * mVec
    fVec6d = rb6d @ mVec6d
    assert isinstance(fVec, sva.ForceVecd)
    assert np.linalg.norm(fVec6d - fVec.vector()) < TOL

def test_abinertiad_left_operators():
    M = np.array([[1, 2, 3],[2, 1, 4], [3, 4, 1]], dtype=float)
    H = np.random.randn(3,3)
    I = np.array([[1, 2, 3],[2, 1, 4], [3, 4, 1]], dtype=float)

    ab = sva.ABInertiad(M, H, I)
    ab6d = ab.matrix()

    mass = 1.
    h = np.random.randn(3) * 100
    rb = sva.RBInertiad(mass, h, I)
    rb6d = rb.matrix()

    w = np.random.randn(3)
    v = np.random.randn(3)
    mVec = sva.MotionVecd(w, v)
    mVec6d = mVec.vector()

    abRes = ab + rb
    abRes6d = ab6d + rb6d
    assert isinstance(abRes, sva.ABInertiad)
    assert np.linalg.norm(abRes6d - abRes.matrix()) < TOL
    assert is_upper_null(abRes.lowerTriangularMassMatrix())
    assert is_upper_null(abRes.lowerTriangularInertia())

    fVec = ab * mVec
    fVec6d = ab6d @ mVec6d
    assert isinstance(fVec, sva.ForceVecd)
    assert np.linalg.norm(fVec6d - fVec.vector()) < TOL
