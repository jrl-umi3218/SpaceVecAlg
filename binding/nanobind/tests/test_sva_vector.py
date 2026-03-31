import numpy as np
import sva

TOL = 1e-5

def test_motionvecd_basic():
    w = np.random.randn(3)
    v = np.random.randn(3)
    vec = sva.MotionVecd(w, v)
    m = vec.vector()

    assert np.allclose(w, vec.angular(), atol=TOL)
    assert np.allclose(v, vec.linear(), atol=TOL)
    assert np.allclose(m, np.concatenate([w, v]), atol=TOL)
    assert np.allclose((5.*vec).vector(), 5.*m, atol=TOL)
    assert np.allclose((vec*5.).vector(), 5.*m, atol=TOL)
    assert np.allclose((vec/5.).vector(), m/5, atol=TOL)
    assert np.allclose((-vec).vector(), -m, atol=TOL)

def test_motionvecd_operators():
    w = np.random.randn(3)
    v = np.random.randn(3)
    vec = sva.MotionVecd(w, v)
    w2 = np.random.randn(3)
    v2 = np.random.randn(3)
    vec2 = sva.MotionVecd(w2, v2)
    m = vec.vector()
    m2 = vec2.vector()

    assert np.allclose((vec + vec2).vector(), (m + m2), atol=TOL)
    assert np.allclose((vec - vec2).vector(), (m - m2), atol=TOL)

    vec_pluseq = sva.MotionVecd(vec)
    vec_pluseq += vec2
    assert np.allclose(vec_pluseq.vector(), (vec + vec2).vector(), atol=TOL)

    vec_minuseq = sva.MotionVecd(vec)
    vec_minuseq -= vec2
    assert np.allclose(vec_minuseq.vector(), (vec - vec2).vector(), atol=TOL)

    assert np.allclose(vec.vector(), vec.vector(), atol=TOL)
    assert not np.allclose(vec.vector(), (-vec).vector(), atol=TOL)

    z = sva.MotionVecd(np.zeros(3), np.zeros(3))
    assert np.allclose(sva.MotionVecd.Zero().vector(), z.vector(), atol=TOL)

def test_forcevecd_basic():
    n = np.random.randn(3)
    f = np.random.randn(3)
    vec = sva.ForceVecd(n, f)
    m = vec.vector()

    assert np.allclose(n, vec.couple(), atol=TOL)
    assert np.allclose(f, vec.force(), atol=TOL)
    assert np.allclose(m, np.concatenate([n, f]), atol=TOL)
    assert np.allclose((5.*vec).vector(), 5.*m, atol=TOL)
    assert np.allclose((vec*5.).vector(), 5.*m, atol=TOL)
    assert np.allclose((vec/5.).vector(), m/5, atol=TOL)
    assert np.allclose((-vec).vector(), -m, atol=TOL)

def test_forcevecd_operators():
    n = np.random.randn(3)
    f = np.random.randn(3)
    vec = sva.ForceVecd(n, f)
    n2 = np.random.randn(3)
    f2 = np.random.randn(3)
    vec2 = sva.ForceVecd(n2, f2)
    m = vec.vector()
    m2 = vec2.vector()

    assert np.allclose((vec + vec2).vector(), (m + m2), atol=TOL)
    assert np.allclose((vec - vec2).vector(), (m - m2), atol=TOL)

    vec_pluseq = sva.ForceVecd(vec)
    vec_pluseq += vec2
    assert np.allclose(vec_pluseq.vector(), (vec + vec2).vector(), atol=TOL)

    vec_minuseq = sva.ForceVecd(vec)
    vec_minuseq -= vec2
    assert np.allclose(vec_minuseq.vector(), (vec - vec2).vector(), atol=TOL)

    assert np.allclose(vec.vector(), vec.vector(), atol=TOL)
    assert not np.allclose(vec.vector(), (-vec).vector(), atol=TOL)

    z = sva.ForceVecd(np.zeros(3), np.zeros(3))
    assert np.allclose(sva.ForceVecd.Zero().vector(), z.vector(), atol=TOL)



def test_motionvecd_left_operators():
    w = np.random.randn(3) * 100
    v = np.random.randn(3) * 100
    n = np.random.randn(3) * 100
    f = np.random.randn(3) * 100

    mVec = sva.MotionVecd(w, v)
    fVec = sva.ForceVecd(n, f)

    mm = mVec.vector()
    mf = fVec.vector()

    # Dot product
    dot_val = mVec.dot(fVec)
    assert np.isclose(dot_val, np.dot(mm, mf), atol=TOL)

    # Cross product
    w2 = np.random.randn(3) * 100
    v2 = np.random.randn(3) * 100
    mVec2 = sva.MotionVecd(w2, v2)
    mm2 = mVec2.vector()

    crossM = mVec.cross(mVec2)
    assert isinstance(crossM, sva.MotionVecd)
    # Assuming sva.vector6ToCrossMatrix exists and returns a 6x6 numpy array
    cross_matrix = sva.vector6ToCrossMatrix(mm)
    assert np.allclose(crossM.vector(), cross_matrix @ mm2, atol=TOL)

    crossF = mVec.crossDual(fVec)
    assert isinstance(crossF, sva.ForceVecd)
    # Assuming sva.vector6ToCrossDualMatrix exists and returns a 6x6 numpy array
    cross_dual_matrix = sva.vector6ToCrossDualMatrix(mm)
    assert np.allclose(crossF.vector(), cross_dual_matrix @ mf, atol=TOL)


def test_impedancevecd_basic_and_operators():
    w = np.random.randn(3)
    v = np.random.randn(3)
    vec = sva.ImpedanceVecd(w, v)
    z = vec.vector()

    assert np.allclose(w, vec.angular(), atol=TOL)
    assert np.allclose(v, vec.linear(), atol=TOL)
    assert np.allclose(z, np.concatenate([w, v]), atol=TOL)

    assert np.allclose((5.*vec).vector(), 5.*z, atol=TOL)
    assert np.allclose((vec*5.).vector(), 5.*z, atol=TOL)
    assert np.allclose((vec/5.).vector(), z/5., atol=TOL)

    vec2 = sva.ImpedanceVecd(np.random.randn(3), np.random.randn(3))
    z2 = vec2.vector()

    assert np.allclose((vec + vec2).vector(), (z + z2), atol=TOL)

    vec_pluseq = sva.ImpedanceVecd(vec)
    assert np.allclose(vec_pluseq.vector(), vec.vector(), atol=TOL)
    vec_pluseq += vec2
    assert np.allclose(vec_pluseq.vector(), (vec + vec2).vector(), atol=TOL)

    assert np.allclose(vec.vector(), vec.vector(), atol=TOL)
    assert not np.any(vec.vector() != vec.vector())

    # *= and /=
    vec_cpy = sva.ImpedanceVecd(w, v)
    vec_cpy *= 5.
    assert np.allclose(vec_cpy.vector(), 5.*z, atol=TOL)
    vec_cpy /= 5.
    assert np.allclose(vec_cpy.vector(), z, atol=TOL)

    # Multiplication with MotionVecd
    mv = sva.MotionVecd(np.random.randn(3), np.random.randn(3))
    fv = vec * mv
    assert isinstance(fv, sva.ForceVecd)
    res = np.array([x*y for x, y in zip(vec.vector(), mv.vector())])
    assert np.allclose(res, fv.vector(), atol=TOL)
    fv2 = mv * vec
    assert np.allclose(fv.vector(), fv2.vector(), atol=TOL)

    # Scalar constructor
    hv = sva.ImpedanceVecd(11., 42.)
    assert np.allclose(hv.angular(), np.full(3, 11.), atol=TOL)
    assert np.allclose(hv.linear(), np.full(3, 42.), atol=TOL)

    z0 = sva.ImpedanceVecd(0., 0.)
    assert np.allclose(sva.ImpedanceVecd.Zero().vector(), z0.vector(), atol=TOL)


def test_admittancevecd_basic_and_operators():
    w = np.random.randn(3)
    v = np.random.randn(3)
    vec = sva.AdmittanceVecd(w, v)
    z = vec.vector()

    assert np.allclose(w, vec.angular(), atol=TOL)
    assert np.allclose(v, vec.linear(), atol=TOL)
    assert np.allclose(z, np.concatenate([w, v]), atol=TOL)

    assert np.allclose((5.*vec).vector(), 5.*z, atol=TOL)
    assert np.allclose((vec*5.).vector(), 5.*z, atol=TOL)
    assert np.allclose((vec/5.).vector(), z/5., atol=TOL)

    vec2 = sva.AdmittanceVecd(np.random.randn(3), np.random.randn(3))
    z2 = vec2.vector()

    assert np.allclose((vec + vec2).vector(), (z + z2), atol=TOL)

    vec_pluseq = sva.AdmittanceVecd(vec)
    assert np.allclose(vec_pluseq.vector(), vec.vector(), atol=TOL)
    vec_pluseq += vec2
    assert np.allclose(vec_pluseq.vector(), (vec + vec2).vector(), atol=TOL)

    assert np.allclose(vec.vector(), vec.vector(), atol=TOL)
    assert not np.any(vec.vector() != vec.vector())

    # *= and /=
    vec_cpy = sva.AdmittanceVecd(w, v)
    vec_cpy *= 5.
    assert np.allclose(vec_cpy.vector(), 5.*z, atol=TOL)
    vec_cpy /= 5.
    assert np.allclose(vec_cpy.vector(), z, atol=TOL)

    # Multiplication with ForceVecd
    fv = sva.ForceVecd(np.random.randn(3), np.random.randn(3))
    mv = vec * fv
    assert isinstance(mv, sva.MotionVecd)
    res = np.array([x*y for x, y in zip(vec.vector(), fv.vector())])
    assert np.allclose(res, mv.vector(), atol=TOL)
    mv2 = fv * vec
    assert np.allclose(mv.vector(), mv2.vector(), atol=TOL)

    # Scalar constructor
    hv = sva.AdmittanceVecd(11., 42.)
    assert np.allclose(hv.angular(), np.full(3, 11.), atol=TOL)
    assert np.allclose(hv.linear(), np.full(3, 42.), atol=TOL)

    z0 = sva.AdmittanceVecd(0., 0.)
    assert np.allclose(sva.AdmittanceVecd.Zero().vector(), z0.vector(), atol=TOL)
