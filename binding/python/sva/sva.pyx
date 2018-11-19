# distutils: language = c++

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
cimport sva.c_sva as c_sva
cimport sva.c_sva_private as c_sva_private
cimport eigen.eigen as eigen
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

cdef class AdmittanceVecd(object):
  def __copyctor__(self, AdmittanceVecd other):
    self.impl = c_sva.AdmittanceVecd(other.impl)
  def __v6ctor__(self, eigen.Vector6d other):
    self.impl = c_sva.AdmittanceVecd(other.impl)
  def __v3v3ctor__(self, eigen.Vector3d angular, eigen.Vector3d linear):
    self.impl = c_sva.AdmittanceVecd(angular.impl, linear.impl)
  def __homoctor__(self, double angular, double linear):
    self.impl = c_sva_private.AdmittanceVecdHomo(angular, linear)
  def __cinit__(self, *args, skip_alloc = False):
    if len(args) == 0:
      self.impl = c_sva.AdmittanceVecd()
    elif len(args) == 1 and isinstance(args[0], AdmittanceVecd):
      self.__copyctor__(args[0])
    elif len(args) == 1 and isinstance(args[0], eigen.Vector6d):
      self.__v6ctor__(args[0])
    elif len(args) == 2 and isinstance(args[0], eigen.Vector3d) and isinstance(args[1], eigen.Vector3d):
      self.__v3v3ctor__(args[0], args[1])
    elif len(args) == 2 and isinstance(args[0], list) and isinstance(args[1], list):
      self.__v3v3ctor__(eigen.Vector3d(args[0]), eigen.Vector3d(args[1]))
    elif len(args) == 2:
      self.__homoctor__(args[0], args[1])
    else:
      raise TypeError("Invalid arguments passed to AdmittanceVecd ctor")

  def angular(self):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>self.impl.angular())
  def linear(self):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>self.impl.linear())
  def vector(self):
    return eigen.Vector6dFromC(<c_eigen.Vector6d>self.impl.vector())

  def __repr__(self):
    return "sva.AdmittanceVecd"
  def __str__(self):
    return c_sva_private.AdmittanceVecdToString(self.impl)

  def __add__(AdmittanceVecd self, AdmittanceVecd other):
    return AdmittanceVecdFromC(self.impl + other.impl)
  def __iadd__(self, AdmittanceVecd other):
    c_sva_private.av_iadd(self.impl, other.impl)
    return self

  def __mul(self, double s):
    return AdmittanceVecdFromC(self.impl*s)
  def __mul_fv(self, ForceVecd fv):
    return MotionVecdFromC(c_sva_private.mul_av_fv(self.impl, deref(fv.impl)))
  def __mul__(self, other):
    if isinstance(self, AdmittanceVecd):
      if isinstance(other, ForceVecd):
        return self.__mul_fv(other)
      return self.__mul(other)
    else:
      return other.__mul__(self)
  def __imul__(self, double other):
    c_sva_private.av_imul(self.impl, other)
    return self

  def __div__(AdmittanceVecd self, double other):
    return self.__truediv__(other)
  def __truediv__(AdmittanceVecd self, double other):
    return AdmittanceVecdFromC(self.impl/other)
  def __idiv__(self, double other):
    return self.__itruediv__(other)
  def __itruediv__(self, double other):
    c_sva_private.av_idiv(self.impl, other)
    return self

  def __richcmp__(AdmittanceVecd self, AdmittanceVecd other, int op):
    if op == 2:
      return self.impl == other.impl
    elif op == 3:
      return self.impl != other.impl
    else:
      raise NotImplementedError("This comparison is not supported")

  @staticmethod
  def Zero():
    return AdmittanceVecdFromC(c_sva_private.AdmittanceVecdZero())
  @staticmethod
  def pickle(mv):
    return AdmittanceVecd, (list(mv.angular()), list(mv.linear()))

cdef AdmittanceVecd AdmittanceVecdFromC(const c_sva.AdmittanceVecd& av):
  cdef AdmittanceVecd ret = AdmittanceVecd()
  ret.impl = av
  return ret

cdef class ImpedanceVecd(object):
  def __copyctor__(self, ImpedanceVecd other):
    self.impl = c_sva.ImpedanceVecd(other.impl)
  def __v6ctor__(self, eigen.Vector6d other):
    self.impl = c_sva.ImpedanceVecd(other.impl)
  def __v3v3ctor__(self, eigen.Vector3d angular, eigen.Vector3d linear):
    self.impl = c_sva.ImpedanceVecd(angular.impl, linear.impl)
  def __homoctor__(self, double angular, double linear):
    self.impl = c_sva_private.ImpedanceVecdHomo(angular, linear)
  def __cinit__(self, *args, skip_alloc = False):
    if len(args) == 0:
      self.impl = c_sva.ImpedanceVecd()
    elif len(args) == 1 and isinstance(args[0], ImpedanceVecd):
      self.__copyctor__(args[0])
    elif len(args) == 1 and isinstance(args[0], eigen.Vector6d):
      self.__v6ctor__(args[0])
    elif len(args) == 2 and isinstance(args[0], eigen.Vector3d) and isinstance(args[1], eigen.Vector3d):
      self.__v3v3ctor__(args[0], args[1])
    elif len(args) == 2 and isinstance(args[0], list) and isinstance(args[1], list):
      self.__v3v3ctor__(eigen.Vector3d(args[0]), eigen.Vector3d(args[1]))
    elif len(args) == 2:
      self.__homoctor__(args[0], args[1])
    else:
      raise TypeError("Invalid arguments passed to ImpedanceVecd ctor")

  def angular(self):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>self.impl.angular())
  def linear(self):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>self.impl.linear())
  def vector(self):
    return eigen.Vector6dFromC(<c_eigen.Vector6d>self.impl.vector())

  def __repr__(self):
    return "sva.ImpedanceVecd"
  def __str__(self):
    return c_sva_private.ImpedanceVecdToString(self.impl)

  def __add__(ImpedanceVecd self, ImpedanceVecd other):
    return ImpedanceVecdFromC(self.impl + other.impl)
  def __iadd__(self, ImpedanceVecd other):
    c_sva_private.iv_iadd(self.impl, other.impl)
    return self

  def __mul(self, double s):
    return ImpedanceVecdFromC(self.impl*s)
  def __mul_mv(self, MotionVecd mv):
    return ForceVecdFromC(c_sva_private.mul_iv_mv(self.impl, deref(mv.impl)))
  def __mul__(self, other):
    if isinstance(self, ImpedanceVecd):
      if isinstance(other, MotionVecd):
        return self.__mul_mv(other)
      return self.__mul(other)
    else:
      return other.__mul__(self)
  def __imul__(self, double other):
    c_sva_private.iv_imul(self.impl, other)
    return self

  def __div__(ImpedanceVecd self, double other):
    return self.__truediv__(other)
  def __truediv__(ImpedanceVecd self, double other):
    return ImpedanceVecdFromC(self.impl/other)
  def __idiv__(self, double other):
    return self.__itruediv__(other)
  def __itruediv__(self, double other):
    c_sva_private.iv_idiv(self.impl, other)
    return self

  def __richcmp__(ImpedanceVecd self, ImpedanceVecd other, int op):
    if op == 2:
      return self.impl == other.impl
    elif op == 3:
      return self.impl != other.impl
    else:
      raise NotImplementedError("This comparison is not supported")

  @staticmethod
  def Zero():
    return ImpedanceVecdFromC(c_sva_private.ImpedanceVecdZero())
  @staticmethod
  def pickle(mv):
    return ImpedanceVecd, (list(mv.angular()), list(mv.linear()))

cdef ImpedanceVecd ImpedanceVecdFromC(const c_sva.ImpedanceVecd& iv):
  cdef ImpedanceVecd ret = ImpedanceVecd()
  ret.impl = iv
  return ret

cdef class ForceVecd(object):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __copyctor__(self, ForceVecd other):
    self.impl = new c_sva.ForceVecd(deref(other.impl))
  def __v6ctor__(self, eigen.Vector6d other):
    self.impl = new c_sva.ForceVecd(other.impl)
  def __v3v3ctor__(self, eigen.Vector3d couple, eigen.Vector3d force):
    self.impl = new c_sva.ForceVecd(couple.impl,force.impl)
  def __cinit__(self, *args, skip_alloc = False):
    self.__own_impl = True
    if len(args) == 0:
      if not skip_alloc:
        self.impl = new c_sva.ForceVecd()
    elif len(args) == 1 and isinstance(args[0], ForceVecd):
      self.__copyctor__(args[0])
    elif len(args) == 1 and isinstance(args[0], eigen.Vector6d):
      self.__v6ctor__(args[0])
    elif len(args) == 2 and isinstance(args[0], eigen.Vector3d) and isinstance(args[1], eigen.Vector3d):
      self.__v3v3ctor__(args[0], args[1])
    elif len(args) == 2 and isinstance(args[0], list) and isinstance(args[1], list):
      self.__v3v3ctor__(eigen.Vector3d(args[0]), eigen.Vector3d(args[1]))
    else:
      raise TypeError("Invalid arguments passed to ForceVecd ctor")

  def couple(self):
    cdef eigen.Vector3d ret = eigen.Vector3d()
    ret.impl = <c_eigen.Vector3d>(self.impl.couple())
    return ret
  def moment(self):
    cdef eigen.Vector3d ret = eigen.Vector3d()
    ret.impl = <c_eigen.Vector3d>(self.impl.moment())
    return ret
  def force(self):
    cdef eigen.Vector3d ret = eigen.Vector3d()
    ret.impl = <c_eigen.Vector3d>(self.impl.force())
    return ret
  def vector(self):
    cdef eigen.Vector6d ret = eigen.Vector6d()
    ret.impl = <c_eigen.Vector6d>(self.impl.vector())
    return ret

  def __repr__(self):
    return "sva.ForceVecd"
  def __str__(self):
    return c_sva_private.ForceVecToString[double](deref(self.impl))

  def __add__(ForceVecd self, ForceVecd other):
    return ForceVecdFromC(deref(self.impl) + deref(other.impl))
  def __iadd__(self, ForceVecd other):
    c_sva_private.fv_iadd[double](self.impl, other.impl)
    return self

  def __neg__(self):
    return ForceVecdFromC(-deref(self.impl))

  def __sub__(ForceVecd self, ForceVecd other):
    return ForceVecdFromC(deref(self.impl) - deref(other.impl))
  def __isub__(self, ForceVecd other):
    c_sva_private.fv_isub[double](self.impl, other.impl)
    return self

  def __mul(self, double s):
    return ForceVecdFromC(deref(self.impl)*s)
  def __mul__(self, other):
    if isinstance(self, ForceVecd):
      if isinstance(other, AdmittanceVecd):
        return other.__mul__(self)
      return self.__mul(other)
    else:
      return other.__mul__(self)
  def __imul__(self, double other):
    c_sva_private.fv_imul(self.impl, other)
    return self

  def __div__(ForceVecd self, double other):
    return self.__truediv__(other)
  def __truediv__(ForceVecd self, double other):
    return ForceVecdFromC(deref(self.impl)/other)
  def __idiv__(ForceVecd self, double other):
    return self.__itruediv__(other)
  def __itruediv__(self, double other):
    c_sva_private.fv_idiv(self.impl, other)
    return self

  def __richcmp__(ForceVecd self, ForceVecd other, int op):
    if op == 2:
      return deref(self.impl) == deref(other.impl)
    elif op == 3:
      return deref(self.impl) != deref(other.impl)
    else:
      raise NotImplementedError("This comparison is not supported")

  @staticmethod
  def Zero():
    return ForceVecdFromC(c_sva_private.ForceVecdZero())
  @staticmethod
  def pickle(fv):
    return ForceVecd, (list(fv.couple()), list(fv.force()))

cdef ForceVecd ForceVecdFromC(const c_sva.ForceVecd& fv, cppbool copy = True):
  cdef ForceVecd ret = ForceVecd(skip_alloc = True)
  if copy:
    ret.impl = new c_sva.ForceVecd(fv)
  else:
    ret.__own_impl = False
    ret.impl = &(c_sva_private.const_cast_fvd(fv))
  return ret

cdef class ForceVecdVector(object):
  def __addForceVecd(self, ForceVecd pt):
    self.v.push_back(deref(pt.impl))
  def __cinit__(self, *args):
    if len(args) == 1 and isinstance(args[0], list):
      for pt in args[0]:
        self.__addForceVecd(pt)
    elif len(args) == 1 and isinstance(args[0], ForceVecd):
      self.__addForceVecd(args[0])
    else:
      for pt in args:
        self.__addForceVecd(pt)

cdef class MotionVecd(object):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __copyctor__(self, MotionVecd other):
    self.impl = new c_sva.MotionVecd(deref(other.impl))
  def __v6ctor__(self, eigen.Vector6d other):
    self.impl = new c_sva.MotionVecd(other.impl)
  def __v3v3ctor__(self, eigen.Vector3d angular, eigen.Vector3d linear):
    self.impl = new c_sva.MotionVecd(angular.impl, linear.impl)
  def __cinit__(self, *args, skip_alloc = False):
    self.__own_impl = True
    if len(args) == 0:
      if not skip_alloc:
        self.impl = new c_sva.MotionVecd()
    elif len(args) == 1 and isinstance(args[0], MotionVecd):
      self.__copyctor__(args[0])
    elif len(args) == 1 and isinstance(args[0], eigen.Vector6d):
      self.__v6ctor__(args[0])
    elif len(args) == 2 and isinstance(args[0], eigen.Vector3d) and isinstance(args[1], eigen.Vector3d):
      self.__v3v3ctor__(args[0], args[1])
    elif len(args) == 2 and isinstance(args[0], list) and isinstance(args[1], list):
      self.__v3v3ctor__(eigen.Vector3d(args[0]), eigen.Vector3d(args[1]))
    else:
      raise TypeError("Invalid arguments passed to MotionVecd ctor")

  def angular(self):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>self.impl.angular())
  def linear(self):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>self.impl.linear())
  def vector(self):
    return eigen.Vector6dFromC(<c_eigen.Vector6d>self.impl.vector())

  def cross(self, MotionVecd other):
    return MotionVecdFromC(self.impl.cross(deref(other.impl)))
  def crossDual(self, ForceVecd other):
    return ForceVecdFromC(self.impl.crossDual(deref(other.impl)))
  def dot(self, ForceVecd other):
    return self.impl.dot(deref(other.impl))

  def __repr__(self):
    return "sva.MotionVecd"
  def __str__(self):
    return c_sva_private.MotionVecToString[double](deref(self.impl))

  def __add__(MotionVecd self, MotionVecd other):
    return MotionVecdFromC(deref(self.impl) + deref(other.impl))
  def __iadd__(self, MotionVecd other):
    c_sva_private.mv_iadd[double](self.impl, other.impl)
    return self

  def __neg__(self):
    return MotionVecdFromC(-deref(self.impl))

  def __sub__(MotionVecd self, MotionVecd other):
    return MotionVecdFromC(deref(self.impl) - deref(other.impl))
  def __isub__(self, MotionVecd other):
    c_sva_private.mv_isub(self.impl, other.impl)
    return self

  def __mul(self, double s):
    return MotionVecdFromC(deref(self.impl)*s)
  def __mul__(self, other):
    if isinstance(self, MotionVecd):
      if isinstance(other, ImpedanceVecd):
        return other.__mul__(self)
      return self.__mul(other)
    else:
      return other.__mul__(self)
  def __imul__(self, double other):
    c_sva_private.mv_imul(self.impl, other)
    return self

  def __div__(MotionVecd self, double other):
    return self.__truediv__(other)
  def __truediv__(MotionVecd self, double other):
    return MotionVecdFromC(deref(self.impl)/other)
  def __idiv__(self, double other):
    return self.__itruediv__(other)
  def __itruediv__(self, double other):
    c_sva_private.mv_idiv(self.impl, other)
    return self

  def __richcmp__(MotionVecd self, MotionVecd other, int op):
    if op == 2:
      return deref(self.impl) == deref(other.impl)
    elif op == 3:
      return deref(self.impl) != deref(other.impl)
    else:
      raise NotImplementedError("This comparison is not supported")

  @staticmethod
  def Zero():
    return MotionVecdFromC(c_sva_private.MotionVecdZero())
  @staticmethod
  def pickle(mv):
    return MotionVecd, (list(mv.angular()), list(mv.linear()))

cdef MotionVecd MotionVecdFromC(const c_sva.MotionVecd& mv, cppbool copy = True):
  cdef MotionVecd ret = MotionVecd(skip_alloc = True)
  if copy:
      ret.impl = new c_sva.MotionVecd(mv)
  else:
      ret.__own_impl = False
      ret.impl = &(c_sva_private.const_cast_mvd(mv))
  return ret

cdef class MotionVecdVector(object):
  def __addMotionVecd(self, MotionVecd pt):
    self.v.push_back(deref(pt.impl))
  def __cinit__(self, *args):
    if len(args) == 1 and isinstance(args[0], list):
      for pt in args[0]:
        self.__addMotionVecd(pt)
    elif len(args) == 1 and isinstance(args[0], MotionVecd):
      self.__addMotionVecd(args[0])
    else:
      for pt in args:
        self.__addMotionVecd(pt)
  def append(self, MotionVecd pt):
    self.__addMotionVecd(pt)
  def __iter__(self):
    for mv in self.v:
      yield MotionVecdFromC(mv, copy = False)
  def __len__(self):
    return self.v.size()

cdef class MotionVecdVectorVector(object):
  def __addMotionVecd(self, int i, MotionVecd pt):
    self.v[i].push_back(deref(pt.impl))
  def __cinit__(self, *args):
    self.v.resize(len(args))
    for i in xrange(len(args)):
      for pt in args[i]:
        self.__addMotionVecd(i, pt)

cdef class RBInertiad(object):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __copyctor__(self, RBInertiad other):
    self.impl = new c_sva.RBInertiad(deref(other.impl))
  def __ctor__(self, double mass, eigen.Vector3d momentum, eigen.Matrix3d inertia_m):
    self.impl = new c_sva.RBInertiad(mass, momentum.impl, inertia_m.impl)
  def __cinit__(self, *args, skip_alloc = False):
    self.__own_impl = True
    if len(args) == 0:
      if not skip_alloc:
        self.impl = new c_sva.RBInertiad()
    elif len(args) == 1 and isinstance(args[0], RBInertiad):
      self.__copyctor__(args[0])
    elif len(args) == 3 and isinstance(args[1], eigen.Vector3d) and isinstance(args[2], eigen.Matrix3d):
      self.__ctor__(args[0], args[1], args[2])
    elif len(args) == 3 and isinstance(args[1], list) and isinstance(args[2], list):
      self.__ctor__(args[0], eigen.Vector3d(args[1]), eigen.Matrix3d(args[2]))
    else:
      raise TypeError("Invalid arguments passed to RBInertiad ctor")

  def mass(self):
    return self.impl.mass()
  def momentum(self):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>(self.impl.momentum()))
  def inertia(self):
    return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(self.impl.inertia()))
  def matrix(self):
    return eigen.Matrix6dFromC(<c_eigen.Matrix6d>(self.impl.matrix()))
  def lowerTriangularInertia(self):
    return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(self.impl.lowerTriangularInertia()))

  def __repr__(self):
    return "sva.RBInertiad"
  def __str__(self):
    return c_sva_private.RBInertiaToString[double](deref(self.impl))

  def __add__(RBInertiad self, RBInertiad other):
    return RBInertiadFromC(deref(self.impl) + deref(other.impl))
  def __iadd__(self, RBInertiad other):
    c_sva_private.rbi_iadd(self.impl, other.impl)
    return self

  def __neg__(self):
    return RBInertiadFromC(-deref(self.impl))

  def __sub__(RBInertiad self, RBInertiad other):
    return RBInertiadFromC(deref(self.impl) - deref(other.impl))
  def __isub__(self, RBInertiad other):
    c_sva_private.rbi_isub(self.impl, other.impl)
    return self

  def __scalar_mul(self, double s):
    return RBInertiadFromC(self.impl.scalar_mul(s))
  def __mvec_mul(self, MotionVecd mv):
    return ForceVecdFromC(self.impl.mvec_mul(deref(mv.impl)))
  def __mul__(self, other):
    if isinstance(self, RBInertiad):
      if isinstance(other, MotionVecd):
        return self.__mvec_mul(other)
      else:
        return self.__scalar_mul(other)
    else:
      return other.__mul__(self)
  def __imul__(self, double other):
    c_sva_private.rbi_imul(self.impl, other)
    return self

  def __richcmp__(RBInertiad self, RBInertiad other, int op):
    if op == 2:
      return deref(self.impl) == deref(other.impl)
    elif op == 3:
      return deref(self.impl) != deref(other.impl)
    else:
      raise NotImplementedError("This comparison is not supported")

  @staticmethod
  def pickle(rb):
    return RBInertiad, (rb.mass(), list(rb.momentum()), list(rb.inertia()))

cdef RBInertiad RBInertiadFromC(const c_sva.RBInertiad& fv, cppbool copy=True):
  cdef RBInertiad ret = RBInertiad(skip_alloc = True)
  if copy:
      ret.impl = new c_sva.RBInertiad(fv)
  else:
      ret.__own_impl = False
      ret.impl = &(c_sva_private.const_cast_rbid(fv))
  return ret

cdef class ABInertiad(object):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __copyctor__(self, ABInertiad other):
    self.impl = new c_sva.ABInertiad(deref(other.impl))
  def __ctor__(self, eigen.Matrix3d mass_m, eigen.Matrix3d generalized_inertia_m, eigen.Matrix3d inertia_m):
    self.impl = new c_sva.ABInertiad(mass_m.impl, generalized_inertia_m.impl, inertia_m.impl)
  def __cinit__(self, *args, skip_alloc = False):
    self.__own_impl = True
    if len(args) == 0:
      if not skip_alloc:
        self.impl = new c_sva.ABInertiad()
    elif len(args) == 1 and isinstance(args[0], ABInertiad):
      self.__copyctor__(args[0])
    elif len(args) == 3 and isinstance(args[0], eigen.Matrix3d) and isinstance(args[1], eigen.Matrix3d) and isinstance(args[2], eigen.Matrix3d):
      self.__ctor__(args[0], args[1], args[2])
    elif len(args) == 3 and isinstance(args[0], list) and isinstance(args[1], list) and isinstance(args[2], list):
      self.__ctor__(eigen.Matrix3d(args[0]), eigen.Matrix3d(args[1]), eigen.Matrix3d(args[2]))
    else:
      raise TypeError("Invalid arguments passed to ABInertiad ctor")

  def massMatrix(self):
    return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(self.impl.massMatrix()))
  def gInertia(self):
    return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(self.impl.gInertia()))
  def inertia(self):
    return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(self.impl.inertia()))
  def matrix(self):
    return eigen.Matrix6dFromC(<c_eigen.Matrix6d>(self.impl.matrix()))
  def lowerTriangularInertia(self):
    return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(self.impl.lowerTriangularInertia()))
  def lowerTriangularMassMatrix(self):
    return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(self.impl.lowerTriangularMassMatrix()))

  def __repr__(self):
    return "sva.ABInertiad"
  def __str__(self):
    return c_sva_private.ABInertiaToString[double](deref(self.impl))

  def __abinertia_add(self, ABInertiad other):
    return ABInertiadFromC(deref(self.impl) + deref(other.impl))
  def __rbinertia_add(self, RBInertiad other):
    return ABInertiadFromC(deref(self.impl) + deref(other.impl))
  def __add__(self, other):
    if isinstance(self, ABInertiad):
      if isinstance(other, ABInertiad):
        return self.__abinertia_add(other)
      elif isinstance(other, RBInertiad):
        return self.__rbinertia_add(other)
      else:
        raise TypeError("ABInertiad or RBInertiad required")
    else:
      return other.__add__(self)
  def __abinertia_iadd(self, ABInertiad other):
    c_sva_private.abi_iadd(self.impl, other.impl)
  def __rbinertia_iadd(self, RBInertiad other):
    c_sva_private.abi_rbi_iadd(self.impl, other.impl)
  def __iadd__(self, other):
    if isinstance(other, ABInertiad):
      self.__abinertia_iadd(other)
    elif isinstance(other, RBInertiad):
      self.__rbinertia_iadd(other)
    else:
      raise TypeError("ABInertiad or RBInertiad required")
    return self

  def __neg__(self):
    return ABInertiadFromC(-deref(self.impl))

  def __sub__(ABInertiad self, ABInertiad other):
    return ABInertiadFromC(deref(self.impl) - deref(other.impl))
  def __isub__(self, ABInertiad other):
    c_sva_private.abi_isub(self.impl, other.impl)
    return self

  def __scalar_mul(self, double s):
    return ABInertiadFromC(self.impl.scalar_mul(s))
  def __mvec_mul(self, MotionVecd mv):
    return ForceVecdFromC(self.impl.mvec_mul(deref(mv.impl)))
  def __mul__(self, other):
    if isinstance(self, ABInertiad):
      if isinstance(other, MotionVecd):
        return self.__mvec_mul(other)
      else:
        return self.__scalar_mul(other)
    else:
      return other.__mul__(self)
  def __imul__(self, double other):
    c_sva_private.abi_imul(self.impl, other)
    return self

  def __richcmp__(ABInertiad self, ABInertiad other, int op):
    if op == 2:
      return deref(self.impl) == deref(other.impl)
    elif op == 3:
      return deref(self.impl) != deref(other.impl)
    else:
      raise NotImplementedError("This comparison is not supported")

  @staticmethod
  def pickle(ab):
    return ABInertiad, (list(ab.massMatrix()), list(ab.gInertia()), list(ab.inertia()))

cdef ABInertiad ABInertiadFromC(const c_sva.ABInertiad& fv, cppbool copy=True):
  cdef ABInertiad ret = ABInertiad(skip_alloc = True)
  if copy:
      ret.impl = new c_sva.ABInertiad(fv)
  else:
      ret.__own_impl = False
      ret.impl = &(c_sva_private.const_cast_abid(fv))
  return ret

cdef class PTransformd(object):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __copyctor__(self, PTransformd other):
    self.impl = new c_sva.PTransformd(deref(other.impl))
  def __rtctor__(self, eigen.Matrix3d rot, eigen.Vector3d trans):
    self.impl = new c_sva.PTransformd(rot.impl, trans.impl)
  def __qtctor__(self, eigen.Quaterniond rot, eigen.Vector3d trans):
    self.impl = new c_sva.PTransformd(rot.impl, trans.impl)
  def __aatctor__(self, eigen.AngleAxisd rot, eigen.Vector3d trans):
    self.impl = new c_sva.PTransformd(rot.impl.matrix(), trans.impl)
  def __mctor__(self, eigen.Matrix3d rot):
    self.impl = new c_sva.PTransformd(rot.impl)
  def __qctor__(self, eigen.Quaterniond rot):
    self.impl = new c_sva.PTransformd(rot.impl)
  def __aactor__(self, eigen.AngleAxisd rot):
    self.impl = new c_sva.PTransformd(rot.impl.matrix())
  def __tctor__(self, eigen.Vector3d trans):
    self.impl = c_sva_private.NewPTransformdFromV3(trans.impl)
  def __cinit__(self, *args, skip_alloc = False):
    self.__own_impl = True
    if len(args) == 0:
      if not skip_alloc:
        self.impl = new c_sva.PTransformd()
    elif len(args) == 1 and isinstance(args[0], PTransformd):
      self.__copyctor__(args[0])
    elif len(args) == 1 and isinstance(args[0], eigen.Quaterniond):
      self.__qctor__(args[0])
    elif len(args) == 1 and isinstance(args[0], eigen.AngleAxisd):
      self.__aactor__(args[0])
    elif len(args) == 1 and isinstance(args[0], eigen.Matrix3d):
      self.__mctor__(args[0])
    elif len(args) == 1 and isinstance(args[0], eigen.Vector3d):
      self.__tctor__(args[0])
    elif len(args) == 2 and isinstance(args[1], eigen.Vector3d):
      if isinstance(args[0], eigen.Quaterniond):
        self.__qtctor__(args[0], args[1])
      elif isinstance(args[0], eigen.Matrix3d):
        self.__rtctor__(args[0], args[1])
      elif isinstance(args[0], eigen.AngleAxisd):
        self.__aatctor__(*args)
      else:
        raise TypeError("Invalid first argument passed to PTransformd ctor")
    elif len(args) == 2 and isinstance(args[0], list) and isinstance(args[1], list):
      self.__rtctor__(eigen.Matrix3d(args[0]), eigen.Vector3d(args[1]))
    else:
      raise TypeError("Invalid arguments passed to PTransformd ctor")

  def rotation(self):
    return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(self.impl.rotation()))
  def translation(self):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>(self.impl.translation()))
  def matrix(self):
    return eigen.Matrix6dFromC(<c_eigen.Matrix6d>(self.impl.matrix()))
  def dualMatrix(self):
    return eigen.Matrix6dFromC(<c_eigen.Matrix6d>(self.impl.dualMatrix()))

  def inv(self):
    return PTransformdFromC(self.impl.inv())

  def __repr__(self):
    return "sva.PTransformd"
  def __str__(self):
    return c_sva_private.PTransformToString[double](deref(self.impl))

  def __pt_mul(self, PTransformd pt):
    return PTransformdFromC(deref(self.impl)*deref(pt.impl))
  def __mvec_mul(self, MotionVecd mv):
    return MotionVecdFromC(deref(self.impl)*deref(mv.impl))
  def __mul__(self, other):
    if isinstance(self, PTransformd):
      if isinstance(other, MotionVecd):
        return self.__mvec_mul(other)
      elif isinstance(other, PTransformd):
        return self.__pt_mul(other)
      else:
        raise TypeError("Unsupported operands PTransformd and {0}".format(type(other)))
    else:
      return other.__mul__(self)

  def invMul(self, MotionVecd other):
    return MotionVecdFromC(self.impl.invMul(deref(other.impl)))
  def __fv_dualMul(self, ForceVecd other):
    return ForceVecdFromC(self.impl.dualMul(deref(other.impl)))
  def __rbi_dualMul(self, RBInertiad other):
    return RBInertiadFromC(self.impl.dualMul(deref(other.impl)))
  def __abi_dualMul(self, ABInertiad other):
    return ABInertiadFromC(self.impl.dualMul(deref(other.impl)))
  def dualMul(self, other):
    if isinstance(other, ForceVecd):
      return self.__fv_dualMul(other)
    elif isinstance(other, RBInertiad):
      return self.__rbi_dualMul(other)
    elif isinstance(other, ABInertiad):
      return self.__abi_dualMul(other)
    else:
      raise TypeError("Cannot call dualMul with type {0}".format(type(other)))

  def __fv_transMul(self, ForceVecd other):
    return ForceVecdFromC(self.impl.transMul(deref(other.impl)))
  def __rbi_transMul(self, RBInertiad other):
    return RBInertiadFromC(self.impl.transMul(deref(other.impl)))
  def __abi_transMul(self, ABInertiad other):
    return ABInertiadFromC(self.impl.transMul(deref(other.impl)))
  def transMul(self, other):
    if isinstance(other, ForceVecd):
      return self.__fv_transMul(other)
    elif isinstance(other, RBInertiad):
      return self.__rbi_transMul(other)
    elif isinstance(other, ABInertiad):
      return self.__abi_transMul(other)
    else:
      raise TypeError("Cannot call transMul with type {0}".format(type(other)))

  def angularMul(self, MotionVecd mv):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>(self.impl.angularMul(deref(mv.impl))))

  def linearMul(self, MotionVecd mv):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>(self.impl.linearMul(deref(mv.impl))))

  def angularInvMul(self, MotionVecd mv):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>(self.impl.angularInvMul(deref(mv.impl))))

  def linearInvMul(self, MotionVecd mv):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>(self.impl.linearInvMul(deref(mv.impl))))

  def coupleDualMul(self, ForceVecd fv):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>(self.impl.coupleDualMul(deref(fv.impl))))

  def forceDualMul(self, ForceVecd fv):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>(self.impl.forceDualMul(deref(fv.impl))))

  def coupleTransMul(self, ForceVecd fv):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>(self.impl.coupleTransMul(deref(fv.impl))))

  def forceTransMul(self, ForceVecd fv):
    return eigen.Vector3dFromC(<c_eigen.Vector3d>(self.impl.forceTransMul(deref(fv.impl))))

  def __richcmp__(PTransformd self, PTransformd other, int op):
    if op == 2:
      return deref(self.impl) == deref(other.impl)
    elif op == 3:
      return deref(self.impl) != deref(other.impl)
    else:
      raise NotImplementedError("This comparison is not supported")

  @staticmethod
  def Identity():
    return PTransformdFromC(c_sva_private.PTransformIdentity[double]())
  @staticmethod
  def pickle(pt):
    return PTransformd, (list(pt.rotation()), list(pt.translation()))

cdef PTransformd PTransformdFromC(const c_sva.PTransformd & pt, cppbool
        copy=True):
  cdef PTransformd ret = PTransformd(skip_alloc = True)
  if copy:
      ret.impl = new c_sva.PTransformd(pt)
  else:
      ret.__own_impl = False
      ret.impl = &(c_sva_private.const_cast_ptd(pt))
  return ret

cdef class PTransformdVector(object):
  def __dealloc__(self):
    if self.__own_impl:
      del self.v
  def __addPTransformd(self, PTransformd pt):
    self.v.push_back(deref(pt.impl))
  def __cinit__(self, *args, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.v = new vector[c_sva.PTransformd]()
    if len(args) == 1 and isinstance(args[0], (list, PTransformdVector)):
      for pt in args[0]:
        self.__addPTransformd(pt)
    elif len(args) == 1 and isinstance(args[0], PTransformd):
      self.__addPTransformd(args[0])
    else:
      for pt in args:
        self.__addPTransformd(pt)
  def append(self, PTransformd pt):
    self.__addPTransformd(pt)
  def __iter__(self):
    for pt in deref(self.v):
      yield PTransformdFromC(pt, copy = False)
  def __len__(self):
    return self.v.size()

cdef PTransformdVector PTransformdVectorFromC(const vector[c_sva.PTransformd]&v,
    cppbool copy=True):
  cdef PTransformdVector ret = PTransformdVector(skip_alloc = False)
  if copy:
    ret.v = new vector[c_sva.PTransformd](v)
  else:
    ret.__own_impl = False
    ret.v = &(c_sva_private.const_cast_pt_vec(v))
  return ret

def RotX(double theta):
  return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(c_sva.RotX[double](theta)))

def RotY(double theta):
  return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(c_sva.RotY[double](theta)))

def RotZ(double theta):
  return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(c_sva.RotZ[double](theta)))

def rotationError(eigen.Matrix3d E_a_b, eigen.Matrix3d E_a_c):
  return eigen.Vector3dFromC(<c_eigen.Vector3d>(c_sva.rotationError[double](E_a_b.impl, E_a_c.impl)))

def rotationVelocity(eigen.Matrix3d E_a_b):
  return eigen.Vector3dFromC(<c_eigen.Vector3d>(c_sva.rotationVelocity[double](E_a_b.impl)))

def transformError(PTransformd X_a_b, PTransformd X_a_c):
  return MotionVecdFromC(c_sva.transformError[double](deref(X_a_b.impl),
      deref(X_a_c.impl)))

def transformVelocity(PTransformd X_a_b):
  return MotionVecdFromC(c_sva.transformVelocity[double](deref(X_a_b.impl)))

def vector3ToCrossMatrix(eigen.Vector3d v):
  return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(c_sva.vector3ToCrossMatrix[double](v.impl)))

def vector6ToCrossMatrix(eigen.Vector6d v):
  return eigen.Matrix6dFromC(<c_eigen.Matrix6d>(c_sva.vector6ToCrossMatrix[double](v.impl)))

def vector6ToCrossDualMatrix(eigen.Vector6d v):
  return eigen.Matrix6dFromC(<c_eigen.Matrix6d>(c_sva.vector6ToCrossDualMatrix[double](v.impl)))

def inertiaToOrigin(eigen.Matrix3d inertia, double mass, eigen.Vector3d com, eigen.Matrix3d rotation):
  return eigen.Matrix3dFromC(<c_eigen.Matrix3d>(c_sva.inertiaToOrigin[double](inertia.impl,
      mass, com.impl, rotation.impl)))

def interpolate(PTransformd frm, PTransformd to, double t = 0.5):
  return PTransformdFromC(c_sva.interpolate[double](deref(frm.impl),
      deref(to.impl), t))

def sinc_inv(double x):
  return c_sva.sinc_inv[double](x)

def copy_reg_pickle():
  # Python 2/3 support
  # Try to import copyreg first (Python 3)
  # Otherwise import copy_reg (Python 2)
  try:
    import copyreg
  except ImportError:
    import copy_reg as copyreg
  for c in [MotionVecd, ForceVecd, PTransformd, RBInertiad, ABInertiad]:
    copyreg.pickle(c, c.pickle)
