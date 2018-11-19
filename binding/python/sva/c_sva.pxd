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

cimport eigen.c_eigen as c_eigen
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "<SpaceVecAlg/SpaceVecAlg>" namespace "sva":
  cdef cppclass AdmittanceVec[T]:
    AdmittanceVec()
    AdmittanceVec(const AdmittanceVec[T] &)
    AdmittanceVec(const c_eigen.Matrix[T,c_eigen.six,c_eigen.one]&)
    AdmittanceVec(const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&,
                  const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&)

    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] angular() const
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] linear() const
    c_eigen.Matrix[T,c_eigen.six,c_eigen.one] vector() const

    AdmittanceVec[T] operator+(const AdmittanceVec[T] &)
    AdmittanceVec[T] operator*(const T &)
    AdmittanceVec[T] operator/(const T &)

    bool operator==(const AdmittanceVec[T] &)
    bool operator!=(const AdmittanceVec[T] &)

  cdef cppclass ImpedanceVec[T]:
    ImpedanceVec()
    ImpedanceVec(const ImpedanceVec[T] &)
    ImpedanceVec(const c_eigen.Matrix[T,c_eigen.six,c_eigen.one]&)
    ImpedanceVec(const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&,
                 const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&)

    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] angular() const
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] linear() const
    c_eigen.Matrix[T,c_eigen.six,c_eigen.one] vector() const

    ImpedanceVec[T] operator+(const ImpedanceVec[T] &)
    ImpedanceVec[T] operator*(const T &)
    ImpedanceVec[T] operator/(const T &)

    bool operator==(const ImpedanceVec[T] &)
    bool operator!=(const ImpedanceVec[T] &)

  cdef cppclass ForceVec[T]:
    ForceVec()
    ForceVec(const ForceVec[T] &)
    ForceVec(const c_eigen.Matrix[T,c_eigen.six,c_eigen.one]&)
    ForceVec(const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&,
             const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&)

    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] couple() const
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] moment() const
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] force() const
    c_eigen.Matrix[T,c_eigen.six,c_eigen.one] vector() const

    ForceVec[T] operator+(const ForceVec[T]&)
    ForceVec[T] operator-()
    ForceVec[T] operator-(const ForceVec[T]&)
    ForceVec[T] operator*(const T&)
    ForceVec[T] operator/(const T&)

    bool operator==(const ForceVec[T]&)
    bool operator!=(const ForceVec[T]&)

  cdef cppclass MotionVec[T]:
    MotionVec()
    MotionVec(const MotionVec[T]&)
    MotionVec(const c_eigen.Matrix[T,c_eigen.six,c_eigen.one] &)
    MotionVec(const c_eigen.Matrix[T,c_eigen.three,c_eigen.one] &,
              const c_eigen.Matrix[T,c_eigen.three,c_eigen.one] &)

    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] angular() const
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] linear() const
    c_eigen.Matrix[T,c_eigen.six,c_eigen.one] vector() const

    MotionVec[T] cross(const MotionVec[T] &)
    ForceVec[T] crossDual(const ForceVec[T] &)
    T dot(const ForceVec[T] &)

    MotionVec[T] operator+(const MotionVec[T] &)
    MotionVec[T] operator-()
    MotionVec[T] operator-(const MotionVec[T] &)
    MotionVec[T] operator*(const T &)
    MotionVec[T] operator/(const T &)

    bool operator==(const MotionVec[T] &)
    bool operator!=(const MotionVec[T] &)

  cdef cppclass RBInertia[T]:
    RBInertia()
    RBInertia(const RBInertia[T]&)
    RBInertia(double, const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&,
                const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&)

    T mass() const
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] momentum() const
    c_eigen.Matrix[T,c_eigen.three,c_eigen.three] inertia() const
    c_eigen.Matrix[T,c_eigen.six,c_eigen.six] matrix() const
    c_eigen.Matrix[T,c_eigen.three,c_eigen.three] lowerTriangularInertia() const

    RBInertia[T] operator+(const RBInertia[T]&)
    RBInertia[T] operator-()
    RBInertia[T] operator-(const RBInertia[T]&)
    RBInertia[T] scalar_mul "operator*"(const T&)
    ForceVec[T] mvec_mul "operator*"(const MotionVec[T]&)

    bool operator==(const RBInertia[T]&)
    bool operator!=(const RBInertia[T]&)

  cdef cppclass ABInertia[T]:
    ABInertia()
    ABInertia(const ABInertia[T]&)
    ABInertia(const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&,
              const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&,
              const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&)

    c_eigen.Matrix[T,c_eigen.three,c_eigen.three] massMatrix()
    c_eigen.Matrix[T,c_eigen.three,c_eigen.three] gInertia()
    c_eigen.Matrix[T,c_eigen.three,c_eigen.three] inertia()
    c_eigen.Matrix[T,c_eigen.six,c_eigen.six] matrix()
    c_eigen.Matrix[T,c_eigen.three,c_eigen.three] lowerTriangularMassMatrix() const
    c_eigen.Matrix[T,c_eigen.three,c_eigen.three] lowerTriangularInertia() const

    ABInertia[T] operator+(const ABInertia[T]&)
    ABInertia[T] operator+(const RBInertia[T]&)
    ABInertia[T] operator-()
    ABInertia[T] operator-(const ABInertia[T]&)
    ABInertia[T] scalar_mul "operator*"(const T&)
    ForceVec[T] mvec_mul "operator*"(const MotionVec[T]&)

    bool operator==(const ABInertia[T]&)
    bool operator!=(const ABInertia[T]&)

  cdef cppclass PTransform[T]:
    PTransform()
    PTransform(const PTransform[T]&)
    PTransform(const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&)
    PTransform(const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&,
               const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&)
    PTransform(const c_eigen.Quaternion[T]&,
               const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&)
    PTransform(const c_eigen.Quaternion[T]&)

    c_eigen.Matrix[T,c_eigen.three,c_eigen.three] rotation() const
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] translation() const
    c_eigen.Matrix[T,c_eigen.six,c_eigen.six] matrix() const
    c_eigen.Matrix[T,c_eigen.six,c_eigen.six] dualMatrix() const

    PTransform[T] inv()

    PTransform[T] operator*(const PTransform[T]&)
    MotionVec[T] operator*(const MotionVec[T]&)

    MotionVec[T] invMul(const MotionVec[T]&)

    ForceVec[T] dualMul(const ForceVec[T]&)
    ForceVec[T] transMul(const ForceVec[T]&)

    RBInertia[T] dualMul(const RBInertia[T]&)
    RBInertia[T] transMul(const RBInertia[T]&)

    ABInertia[T] dualMul(const ABInertia[T]&)
    ABInertia[T] transMul(const ABInertia[T]&)

    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] angularMul(const MotionVec[T]&)
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] linearMul(const MotionVec[T]&)

    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] angularInvMul(const MotionVec[T]&)
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] linearInvMul(const MotionVec[T]&)

    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] coupleDualMul(const ForceVec[T]&)
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] forceDualMul(const ForceVec[T]&)

    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] coupleTransMul(const ForceVec[T]&)
    c_eigen.Matrix[T,c_eigen.three,c_eigen.one] forceTransMul(const ForceVec[T]&)

    bool operator==(const PTransform[T]&)
    bool operator!=(const PTransform[T]&)

  ctypedef AdmittanceVec[double] AdmittanceVecd
  ctypedef ImpedanceVec[double] ImpedanceVecd
  ctypedef ForceVec[double] ForceVecd
  ctypedef MotionVec[double] MotionVecd
  ctypedef RBInertia[double] RBInertiad
  ctypedef ABInertia[double] ABInertiad
  ctypedef PTransform[double] PTransformd

  c_eigen.Matrix[T,c_eigen.three,c_eigen.three] RotX[T](T)
  c_eigen.Matrix[T,c_eigen.three,c_eigen.three] RotY[T](T)
  c_eigen.Matrix[T,c_eigen.three,c_eigen.three] RotZ[T](T)

  c_eigen.Matrix[T,c_eigen.three,c_eigen.three] rotationError[T](const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&, const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&)
  c_eigen.Matrix[T,c_eigen.three,c_eigen.three] rotationVelocity[T](const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&)

  MotionVec[T] transformError[T](const PTransform[T]&, const PTransform[T]&)
  MotionVec[T] transformVelocity[T](const PTransform[T]&)

  c_eigen.Matrix[T,c_eigen.three,c_eigen.three] vector3ToCrossMatrix[T](const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&)
  c_eigen.Matrix[T,c_eigen.six,c_eigen.six] vector6ToCrossMatrix[T](const c_eigen.Matrix[T,c_eigen.six,c_eigen.one]&)
  c_eigen.Matrix[T,c_eigen.six,c_eigen.six] vector6ToCrossDualMatrix[T](const c_eigen.Matrix[T,c_eigen.six,c_eigen.one]&)

  c_eigen.Matrix[T,c_eigen.three,c_eigen.three] inertiaToOrigin[T](const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&, T, const c_eigen.Matrix[T,c_eigen.three,c_eigen.one]&, const c_eigen.Matrix[T,c_eigen.three,c_eigen.three]&)

  PTransform[T] interpolate[T](const PTransform[T]&, const PTransform[T]&, T)

  T sinc_inv[T](T)
