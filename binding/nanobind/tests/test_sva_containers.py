import numpy as np
import sva

TOL = 1e-5

def random_quaternion():
    q = np.random.randn(4)
    q /= np.linalg.norm(q)
    return sva.Quaternion(q)

def random_vec3():
    return np.random.randn(3) * 100

def test_ptransformd_list():
    create_random_pt = lambda: sva.PTransformd(random_quaternion(), random_vec3())

    v1 = []
    for _ in range(100):
        v1.append(create_random_pt())
    assert len(v1) == 100

    pt = create_random_pt()
    v2 = [pt]
    assert len(v2) == 1

    v3 = [create_random_pt() for _ in range(100)]
    assert len(v3) == 100

    v4 = [*([create_random_pt() for _ in range(100)])]
    assert len(v4) == 100

    pt = create_random_pt()
    v5 = [pt] * 100
    assert all(pt == vi for vi in v5)

def test_motionvecd_list():
    create_random_mv = lambda: sva.MotionVecd(random_vec3(), random_vec3())

    v1 = []
    for _ in range(100):
        v1.append(create_random_mv())
    assert len(v1) == 100

    mv = create_random_mv()
    v2 = [mv]
    assert len(v2) == 1

    v3 = [create_random_mv() for _ in range(100)]
    assert len(v3) == 100

    v4 = [*([create_random_mv() for _ in range(100)])]
    assert len(v4) == 100

    mv = create_random_mv()
    v5 = [mv] * 100
    assert all(mv == vi for vi in v5)

def test_forcevecd_list():
    create_random_fv = lambda: sva.ForceVecd(random_vec3(), random_vec3())

    v1 = []
    for _ in range(100):
        v1.append(create_random_fv())
    assert len(v1) == 100

    fv = create_random_fv()
    v2 = [fv]
    assert len(v2) == 1

    v3 = [create_random_fv() for _ in range(100)]
    assert len(v3) == 100

    v4 = [*([create_random_fv() for _ in range(100)])]
    assert len(v4) == 100

    fv = create_random_fv()
    v5 = [fv] * 100
    assert all(fv == vi for vi in v5)

def test_impedancevecd_list():
    create_random_iv = lambda: sva.ImpedanceVecd(random_vec3(), random_vec3())

    v1 = []
    for _ in range(100):
        v1.append(create_random_iv())
    assert len(v1) == 100

    iv = create_random_iv()
    v2 = [iv]
    assert len(v2) == 1

    v3 = [create_random_iv() for _ in range(100)]
    assert len(v3) == 100

    v4 = [*([create_random_iv() for _ in range(100)])]
    assert len(v4) == 100

    iv = create_random_iv()
    v5 = [iv] * 100
    assert all(iv == vi for vi in v5)

def test_admittancevecd_list():
    create_random_av = lambda: sva.AdmittanceVecd(random_vec3(), random_vec3())

    v1 = []
    for _ in range(100):
        v1.append(create_random_av())
    assert len(v1) == 100

    av = create_random_av()
    v2 = [av]
    assert len(v2) == 1

    v3 = [create_random_av() for _ in range(100)]
    assert len(v3) == 100

    v4 = [*([create_random_av() for _ in range(100)])]
    assert len(v4) == 100

    av = create_random_av()
    v5 = [av] * 100
    assert all(av == vi for vi in v5)
