# Copyright 2026 CNRS-UM LIRMM, CNRS-AIST JRL

import pickle
import numpy as np
import sva

def check_pickle(obj, atol=1e-8, rtol=1e-5):
    pickled = pickle.dumps(obj)
    obj2 = pickle.loads(pickled)
    # Try numpy allclose for objects with .vector() or .matrix()
    if hasattr(obj, "vector"):
        assert np.allclose(obj.vector(), obj2.vector(), atol=atol, rtol=rtol)
    elif hasattr(obj, "matrix"):
        assert np.allclose(obj.matrix(), obj2.matrix(), atol=atol, rtol=rtol)
    else:
        # Fallback to equality
        assert obj == obj2

def test_motionvecd_pickle():
    mv = sva.MotionVecd(np.random.randn(6))
    check_pickle(mv)

def test_forcevecd_pickle():
    fv = sva.ForceVecd(np.random.randn(6))
    check_pickle(fv)

def test_ptransformd_pickle():
    quat_np = np.random.randn(4)
    quat_np /= np.linalg.norm(quat_np)
    quat = sva.Quaternion(quat_np)
    trans = np.random.randn(3)
    pt = sva.PTransformd(quat, trans)
    check_pickle(pt)

def test_rbinertiad_pickle():
    mass = np.random.rand()
    h = np.random.randn(3)
    I = np.random.randn(3, 3)
    rb = sva.RBInertiad(mass, h, I)
    check_pickle(rb)

def test_abinertiad_pickle():
    M = np.random.randn(3, 3)
    H = np.random.randn(3, 3)
    I = np.random.randn(3, 3)
    ab = sva.ABInertiad(M, H, I)
    check_pickle(ab)
