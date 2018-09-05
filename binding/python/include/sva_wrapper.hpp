/* Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
 *
 * This file is part of SpaceVecAlg.
 *
 * SpaceVecAlg is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SpaceVecAlg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with SpaceVecAlg.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <SpaceVecAlg/SpaceVecAlg>
#include <sstream>

namespace sva
{

std::string AdmittanceVecdToString(const sva::AdmittanceVecd & av)
{
  std::stringstream ss;
  ss << av;
  return ss.str();
}

std::string ImpedanceVecdToString(const sva::ImpedanceVecd & iv)
{
  std::stringstream ss;
  ss << iv;
  return ss.str();
}

template<typename T>
std::string ForceVecToString(const sva::ForceVec<T> & m)
{
  std::stringstream ss;
  ss << m;
  return ss.str();
}

template<typename T>
std::string MotionVecToString(const sva::MotionVec<T> & m)
{
  std::stringstream ss;
  ss << m;
  return ss.str();
}

template<typename T>
std::string RBInertiaToString(const sva::RBInertia<T> & m)
{
  std::stringstream ss;
  ss << m;
  return ss.str();
}

template<typename T>
std::string ABInertiaToString(const sva::ABInertia<T> & m)
{
  std::stringstream ss;
  ss << m;
  return ss.str();
}

template<typename T>
std::string PTransformToString(const sva::PTransform<T> & m)
{
  std::stringstream ss;
  ss << m;
  return ss.str();
}

template<typename T>
PTransform<T> PTransformIdentity()
{
  return PTransform<T>::Identity();
}

sva::PTransformd * NewPTransformdFromV3(const Eigen::Matrix<double, 3, 1> & v)
{
  return new sva::PTransformd(v);
}

sva::ForceVec<double>& const_cast_fvd(const sva::ForceVec<double> & rhs)
{
  return const_cast<sva::ForceVec<double>&>(rhs);
}

template<typename T>
void fv_iadd(sva::ForceVec<T> * lhs, sva::ForceVec<T> * rhs)
{
  *lhs += *rhs;
}

template<typename T>
void fv_isub(sva::ForceVec<T> * lhs, sva::ForceVec<T> * rhs)
{
  *lhs -= *rhs;
}

template<typename T>
void fv_imul(sva::ForceVec<T> * lhs, double s)
{
  *lhs = (*lhs)*s;
}

template<typename T>
void fv_idiv(sva::ForceVec<T> * lhs, double s)
{
  *lhs = (*lhs)/s;
}

sva::MotionVec<double>& const_cast_mvd(const sva::MotionVec<double> & rhs)
{
  return const_cast<sva::MotionVec<double>&>(rhs);
}

template<typename T>
void mv_iadd(sva::MotionVec<T> * lhs, sva::MotionVec<T> * rhs)
{
  *lhs += *rhs;
}

template<typename T>
void mv_isub(sva::MotionVec<T> * lhs, sva::MotionVec<T> * rhs)
{
  *lhs -= *rhs;
}

template<typename T>
void mv_imul(sva::MotionVec<T> * lhs, double s)
{
  *lhs = (*lhs)*s;
}

template<typename T>
void mv_idiv(sva::MotionVec<T> * lhs, double s)
{
  *lhs = (*lhs)/s;
}

sva::RBInertia<double>& const_cast_rbid(const sva::RBInertia<double> & rhs)
{
  return const_cast<sva::RBInertia<double>&>(rhs);
}

template<typename T>
void rbi_iadd(sva::RBInertia<T> * lhs, sva::RBInertia<T> * rhs)
{
  *lhs += *rhs;
}

template<typename T>
void rbi_isub(sva::RBInertia<T> * lhs, sva::RBInertia<T> * rhs)
{
  *lhs -= *rhs;
}

template<typename T>
void rbi_imul(sva::RBInertia<T> * lhs, double s)
{
  *lhs = (*lhs)*s;
}

sva::ABInertia<double>& const_cast_abid(const sva::ABInertia<double> & rhs)
{
  return const_cast<sva::ABInertia<double>&>(rhs);
}

template<typename T>
void abi_iadd(sva::ABInertia<T> * lhs, sva::ABInertia<T> * rhs)
{
  *lhs += *rhs;
}
template<typename T>
void abi_rbi_iadd(sva::ABInertia<T> * lhs, sva::RBInertia<T> * rhs)
{
  *lhs = *lhs + *rhs;
}

template<typename T>
void abi_isub(sva::ABInertia<T> * lhs, sva::ABInertia<T> * rhs)
{
  *lhs -= *rhs;
}

template<typename T>
void abi_imul(sva::ABInertia<T> * lhs, double s)
{
  *lhs = (*lhs)*s;
}

sva::PTransform<double>& const_cast_ptd(const sva::PTransform<double> & rhs)
{
  return const_cast<sva::PTransform<double>&>(rhs);
}

std::vector<sva::PTransformd>& const_cast_pt_vec(const std::vector<sva::PTransformd> & rhs)
{
  return const_cast<std::vector<sva::PTransformd>&>(rhs);
}

ForceVecd ForceVecdZero()
{
  return ForceVecd::Zero();
}

MotionVecd MotionVecdZero()
{
  return MotionVecd::Zero();
}

AdmittanceVecd AdmittanceVecdZero()
{
  return AdmittanceVecd::Zero();
}

ImpedanceVecd ImpedanceVecdZero()
{
  return ImpedanceVecd::Zero();
}

AdmittanceVecd AdmittanceVecdHomo(double a, double l)
{
  return AdmittanceVecd(a, l);
}

ImpedanceVecd ImpedanceVecdHomo(double a, double l)
{
  return ImpedanceVecd(a, l);
}

MotionVecd mul_av_fv(AdmittanceVecd & lhs, ForceVecd & rhs)
{
  return lhs * rhs;
}
void av_iadd(AdmittanceVecd & lhs, AdmittanceVecd & rhs)
{
  lhs += rhs;
}
void av_imul(AdmittanceVecd & lhs, double s)
{
  lhs *= s;
}
void av_idiv(AdmittanceVecd & lhs, double s)
{
  lhs /= s;
}

ForceVecd mul_iv_mv(ImpedanceVecd & lhs, MotionVecd & rhs)
{
  return lhs * rhs;
}
void iv_iadd(ImpedanceVecd & lhs, ImpedanceVecd & rhs)
{
  lhs += rhs;
}
void iv_imul(ImpedanceVecd & lhs, double s)
{
  lhs *= s;
}
void iv_idiv(ImpedanceVecd & lhs, double s)
{
  lhs /= s;
}

}
