#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from sva.c_sva cimport *
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "sva_wrapper.hpp" namespace "sva":
  string AdmittanceVecdToString(const AdmittanceVecd &)
  string ImpedanceVecdToString(const ImpedanceVecd &)
  string ForceVecToString[T](const ForceVec[T] &)
  string MotionVecToString[T](const MotionVec[T] &)
  string RBInertiaToString[T](const RBInertia[T] &)
  string ABInertiaToString[T](const ABInertia[T] &)
  string PTransformToString[T](const PTransform[T] &)

  PTransform[T] PTransformIdentity[T]()

  PTransformd * NewPTransformdFromV3(const c_eigen.Matrix[double, c_eigen.three, c_eigen.one]&)

  MotionVecd mul_av_fv(AdmittanceVecd & lhs, ForceVecd & rhs)
  void av_iadd(AdmittanceVecd & lhs, AdmittanceVecd & rhs)
  void av_imul(AdmittanceVecd & lhs, double s)
  void av_idiv(AdmittanceVecd & lhs, double s)
  AdmittanceVecd AdmittanceVecdZero()
  AdmittanceVecd AdmittanceVecdHomo(double, double)

  ForceVecd mul_iv_mv(ImpedanceVecd & lhs, MotionVecd & rhs)
  void iv_iadd(ImpedanceVecd & lhs, ImpedanceVecd & rhs)
  void iv_imul(ImpedanceVecd & lhs, double s)
  void iv_idiv(ImpedanceVecd & lhs, double s)
  ImpedanceVecd ImpedanceVecdZero()
  ImpedanceVecd ImpedanceVecdHomo(double, double)

  ForceVecd & const_cast_fvd(const ForceVecd &)
  void fv_iadd[T](ForceVec[T] * lhs, ForceVec[T] * rhs)
  void fv_isub[T](ForceVec[T] * lhs, ForceVec[T] * rhs)
  void fv_imul[T](ForceVec[T] * lhs, double s)
  void fv_idiv[T](ForceVec[T] * lhs, double s)
  ForceVecd ForceVecdZero()

  MotionVecd & const_cast_mvd(const MotionVecd &)
  void mv_iadd[T](MotionVec[T] * lhs, MotionVec[T] * rhs)
  void mv_isub[T](MotionVec[T] * lhs, MotionVec[T] * rhs)
  void mv_imul[T](MotionVec[T] * lhs, double s)
  void mv_idiv[T](MotionVec[T] * lhs, double s)
  MotionVecd MotionVecdZero()

  RBInertiad & const_cast_rbid(const RBInertiad &)
  void rbi_iadd[T](RBInertia[T] * lhs, RBInertia[T] * rhs)
  void rbi_isub[T](RBInertia[T] * lhs, RBInertia[T] * rhs)
  void rbi_imul[T](RBInertia[T] * lhs, double s)

  ABInertiad & const_cast_abid(const ABInertiad &)
  void abi_iadd[T](ABInertia[T] * lhs, ABInertia[T] * rhs)
  void abi_rbi_iadd[T](ABInertia[T] * lhs, RBInertia[T] * rhs)
  void abi_isub[T](ABInertia[T] * lhs, ABInertia[T] * rhs)
  void abi_imul[T](ABInertia[T] * lhs, double s)

  PTransformd & const_cast_ptd(const PTransformd &)
  vector[PTransformd]& const_cast_pt_vec(const vector[PTransformd]&)
