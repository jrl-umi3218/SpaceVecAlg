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

cimport eigen.eigen as eigen
cimport c_sva
from libcpp.vector cimport vector
from libcpp cimport bool as cppbool

cdef class AdmittanceVecd(object):
  cdef c_sva.AdmittanceVecd impl

cdef AdmittanceVecd AdmittanceVecdFromC(const c_sva.AdmittanceVecd&)

cdef class ImpedanceVecd(object):
  cdef c_sva.ImpedanceVecd impl

cdef ImpedanceVecd ImpedanceVecdFromC(const c_sva.ImpedanceVecd&)

cdef class ForceVecd(object):
  cdef c_sva.ForceVecd * impl
  cdef cppbool __own_impl

cdef ForceVecd ForceVecdFromC(const c_sva.ForceVecd&, cppbool copy=?)

cdef class ForceVecdVector(object):
  cdef vector[c_sva.ForceVecd] v

cdef class MotionVecd(object):
  cdef c_sva.MotionVecd * impl
  cdef cppbool __own_impl

cdef MotionVecd MotionVecdFromC(const c_sva.MotionVecd&, cppbool copy=?)

cdef class MotionVecdVector(object):
  cdef vector[c_sva.MotionVecd] v

cdef class MotionVecdVectorVector(object):
  cdef vector[vector[c_sva.MotionVecd]] v

cdef class RBInertiad(object):
  cdef c_sva.RBInertiad * impl
  cdef cppbool __own_impl

cdef RBInertiad RBInertiadFromC(const c_sva.RBInertiad&, cppbool copy=?)

cdef class ABInertiad(object):
  cdef c_sva.ABInertiad * impl
  cdef cppbool __own_impl

cdef ABInertiad ABInertiadFromC(const c_sva.ABInertiad&, cppbool copy=?)

cdef class PTransformd(object):
  cdef c_sva.PTransformd * impl
  cdef cppbool __own_impl

cdef PTransformd PTransformdFromC(const c_sva.PTransformd&, cppbool copy=?)

cdef class PTransformdVector(object):
  cdef vector[c_sva.PTransformd] * v
  cdef cppbool __own_impl

cdef PTransformdVector PTransformdVectorFromC(const vector[c_sva.PTransformd]&,
    cppbool copy=?)
