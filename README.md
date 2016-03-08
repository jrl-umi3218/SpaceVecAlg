[![License LGPL 3](https://img.shields.io/badge/license-LGPLv3-green.svg)](http://www.gnu.org/licenses/lgpl-3.0.txt)
[![Build Status](https://travis-ci.org/jrl-umi3218/SpaceVecAlg.svg?branch=master)](https://travis-ci.org/jrl-umi3218/SpaceVecAlg)

SpaceVecAlg
========

SpaceVecAlg aim to implement Spatial Vector Algebra with the Eigen3 linear algebra library.

All this implementation is based on appendix A of [Roy Featherstone Rigid Body Dynamics Algorithms book](http://www.springer.com/fr/book/9780387743141).

Documentation
-----

Features:
 * Featherstone Spatial Vector Algebra C++11 implementation
 * Header only
 * Use Eigen3 as linear algebra library
 * Python binding

A short tutorial is available [here](https://github.com/jorisv/sva_rbdyn_presentation/blob/master/presentation_release.pdf).

To learn more about Spatial Vector Algebra you can find some presentations on the following [page](http://royfeatherstone.org/spatial/).

The [SpaceVecAlg and RBDyn tutorial](https://github.com/jorisv/sva_rbdyn_tutorials) is also a big ressource to understand [how to use SpaceVecAlg](http://nbviewer.ipython.org/github/jorisv/sva_rbdyn_tutorials/blob/master/SpaceVecAlg.ipynb).
Also you will find a lot of IPython Notebook that will present real use case.

Finally you can build a Doxygen documentation by typing `make doc` in the build directory. After a `make install` the documentation will be in `CMAKE_INSTALL_PREFIX/share/doc/SpaceVecAlg` (see the Installing section).

### Appendix A table transcription to C++

In this section `a` stand for a double, `v` for a motion vector, `f` for a force vector,
`I` for a rigid body inertia, `I^a` for a articulated body inertia and
`X` for a plücker transfrom.

`.` stand for the dot product, `x`for the cross product and `x^{\*}` for the cross product dual.

`r` stand for a 3d translation vector, `E` for a 3d rotation matrix,
`m` for a mass, `c` for the center of mass 3d vector from the body origin,
`I_c` for the 3d rotational inertia matrix at CoM frame.

#### Table A.2 transcription

Operation                    | C++
-----------------------------|----
rx(theta)                    | `sva::RotX(theta)`
ry(theta)                    | `sva::RotY(theta)`
rz(theta)                    | `sva::RotZ(theta)`
X = rotx(theta)              | `sva::PTransformd(sva::RotX(theta))`
X = roty(theta)              | `sva::PTransformd(sva::RotY(theta))`
X = rotz(theta)              | `sva::PTransformd(sva::RotZ(theta))`
X = xlt(r)                   | `sva::PTransformd(r)`
x = crm(v)                   | `sva::vector6ToCrossMatrix(v)`
v x^{\*} = crf(v)            | `sva::vector6ToCrossDualMatrix(v)`
I = E\*mcI(m, c, I_c)\*E{^T} | `inertiaToOrigin(I_c, m, c, E)`
v = XtoV(X)                  | `sva::transformVelocity(X)`


#### Table A.4 transcription

Operation           | C++
--------------------|----
a v                 | `a*sva::MotionVecd()`
a f                 | `a*sva::ForceVecd()`
a I                 | `a*sva::RBInertiad()`
a I^a               | `a*sva::ABInertiad()`
v_1 + v_2           | `sva::MotionVecd() + sva::MotionVecd()`
f_1 + f_2           | `sva::ForceVecd() + sva::ForceVecd()`
I_1 + I_1           | `sva::RBInertiad() + sva::RBInertiad()`
I_1^a + I_2^a       | `sva::ABInertiad() + sva::ABInertiad()`
I_1^a + I_2^a       | `sva::ABInertiad() + sva::RBInertiad()`
v . f               | `sva::MotionVecd().dot(sva::ForceVecd())`
v_1 x v_2           | `sva::MotionVecd().cross(sva::MotionVecd())`
v x^\* f            | `sva::MotionVecd().crossDual(sva::ForceVecd())`
I v                 | `sva::RBInertiad()\*sva::MotionVecd()`
I^a v               | `sva::ABInertiad()\*sva::MotionVecd()`
X_1 X_2             | `sva::PTransformd()\*sva::PTransformd()`
X^{-1}              | `sva::PTransformd().inv()`
X v                 | `sva::PTransformd()\*sva::MotionVecd()`
X^{-1} v            | `sva::PTransformd().invMul(sva::MotionVecd())`
X^{\*} f            | `sva::PTransformd().dualMul(sva::ForceVecd())`
X^{T} f             | `sva::PTransformd().transMul(sva::ForceVecd())`
X^{\*} I X^{-1}     | `sva::PTransformd().dualMul(sva::RBInertiad())`
X^{T} I X           | `sva::PTransformd().transMul(sva::RBInertiad())`
X^{\*} I^a X^{-1}   | `sva::PTransformd().dualMul(sva::ABInertiad())`
X^{T} I^a X         | `sva::PTransformd().transMul(sva::ABInertiad())`

#### Table A.3 transcription

Here `w` stand for the 3d angular velocity, `v` for the 3d linear velocity,
`n` for the 3d torque, `f` for the 3d force, `E` for the 3d rotation matrix,
`r` for the 3d translation vector, `q` for a unit quaternion, `m` for a mass,
`h` for the first moment of mass (h = m c) at body frame,
`I` for the 3d rotational inertia at body frame, `M` for the 3d mass matrix,
and `H` for the 3d generalized inertia matrix.

Operation                    | C++
-----------------------------|----
mv(w, v)                     | `sva::MotionVecd(w, v)`
fv(n, f)                     | `sva::ForceVecd(n, f)`
plx(E, r)                    | `sva::PTransform(E, r)`
plx(q, r)                    | `sva::PTransform(q, r)`
rbi(m, h, I)                 | `sva::RBInertia(m, h, I)`
abi(M, H, I)                 | `sva::ABInertia(M, H, I)`


Installing
------

### Manual

#### Dependencies

To compile you need the following tools:

 * [Git]()
 * [CMake]() >= 2.8
 * [pkg-config]()
 * [doxygen]()
 * [g++]() >= 4.7 (for C++11 support)
 * [Boost](http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html) >= 1.49
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2

For Python bindings:

 * [PyBindGen](https://launchpad.net/pybindgen) = 0.16
 * [Eigen3ToPython](https://github.com/jorisv/Eigen3ToPython) (to use the Python binding)

#### Building

```sh
git clone --recursive https://github.com/jrl-umi3218/SpaceVecAlg
cd SpaceVecAlg
mkdir _build
cd _build
cmake [options] ..
make && make intall
```

Where the main options are:

 * `-DCMAKE_BUIlD_TYPE=Release` Build in Release mode
 * `-DCMAKE_INSTALL_PREFIX=some/path/to/install` default is `/usr/local`
 * `-DPYTHON_BINDING=ON` Build the python binding
 * `-DUNIT_TESTS=ON` Build unit tests.
 * `-DPYTHON_DEB_LAYOUT=OFF` install python library in `site-packages` (ON will install in `dist-packages`)

### Arch Linux

You can use the following [AUR package](https://aur.archlinux.org/packages/spacevecalg-git).


Pulling git subtree
------

To update sync cmake or .travis directory with their upstream git repository:

	git fetch git://github.com/jrl-umi3218/jrl-cmakemodules.git master
	git subtree pull --prefix cmake git://github.com/jrl-umi3218/jrl-cmakemodules.git master --squash

	git fetch git://github.com/jrl-umi3218/jrl-travis.git master
	git subtree pull --prefix .travis git://github.com/jrl-umi3218/jrl-travis.git master --squash
