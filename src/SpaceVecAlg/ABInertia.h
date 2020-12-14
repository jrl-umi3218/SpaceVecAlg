/*
 * Copyright 2012-2020 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "EigenTypedef.h"
#include "fwd.h"

namespace sva
{

/**
 * Spatial Articulated Body Inertia compact representation.
 * See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
 */
template<typename T>
class ABInertia
{
public:
  typedef Eigen::Vector3<T> vector3_t;
  typedef Eigen::Matrix3<T> matrix3_t;
  typedef Eigen::Matrix6<T> matrix6_t;

public:
  ABInertia() : M_(), H_(), I_() {}

  /**
   * @param M Mass matrix.
   * @param H Generalized inertia matrix.
   * @param I Inertia matrix.
   */
  ABInertia(const matrix3_t & M, const matrix3_t & H, const matrix3_t & I)
  : M_(matrix3_t::Zero()), H_(H), I_(matrix3_t::Zero())
  {
    M_.template triangularView<Eigen::Lower>() = M;
    I_.template triangularView<Eigen::Lower>() = I;
  }

  /**
   * @param M Lower triangular view of Mass matrix.
   * @param H Generalized inertia matrix.
   * @param I Lower triangular view Inertia matrix.
   */
  ABInertia(const Eigen::TriangularView<matrix3_t, Eigen::Lower> & ltM,
            const matrix3_t & H,
            const Eigen::TriangularView<matrix3_t, Eigen::Lower> & ltI)
  : M_(ltM), H_(H), I_(ltI)
  {
  }

  // Accessor
  /// @return Mass matrix with a zero upper part.
  const matrix3_t & lowerTriangularMassMatrix() const
  {
    return M_;
  }

  /// @return Mass matrix.
  matrix3_t massMatrix() const
  {
    matrix3_t M;
    M.template triangularView<Eigen::Upper>() = M_.transpose();
    M.template triangularView<Eigen::StrictlyLower>() = M_;
    return M;
  }

  /// @return Generalized inertia matrix.
  const matrix3_t & gInertia() const
  {
    return H_;
  }

  /// @return Inertia matrix with a zero upper part.
  const matrix3_t & lowerTriangularInertia() const
  {
    return I_;
  }

  /// @return Inertia matrix.
  matrix3_t inertia() const
  {
    matrix3_t I;
    I.template triangularView<Eigen::Upper>() = I_.transpose();
    I.template triangularView<Eigen::StrictlyLower>() = I_;
    return I;
  }

  /// @retrun Non compact spatial articulated body inertia matrix.
  matrix6_t matrix() const
  {
    matrix6_t m;
    m << inertia(), H_, H_.transpose(), massMatrix();
    return m;
  }

  template<typename T2>
  ABInertia<T2> cast() const
  {
    return ABInertia<T2>(M_.template cast<T2>(), H_.template cast<T2>(), I_.template cast<T2>());
  }

  // Operators
  ABInertia<T> operator+(const ABInertia<T> & rbI) const
  {
    matrix3_t M, I;
    M.template triangularView<Eigen::Lower>() = M_ + rbI.M_;
    I.template triangularView<Eigen::Lower>() = I_ + rbI.I_;
    return ABInertia<T>(M, H_ + rbI.H_, I);
  }

  ABInertia<T> operator-(const ABInertia<T> & rbI) const
  {
    matrix3_t M, I;
    M.template triangularView<Eigen::Lower>() = M_ - rbI.M_;
    I.template triangularView<Eigen::Lower>() = I_ - rbI.I_;
    return ABInertia<T>(M, H_ - rbI.H_, I);
  }

  ABInertia<T> operator-() const
  {
    return ABInertia<T>(-M_, -H_, -I_);
  }

  ABInertia<T> & operator+=(const ABInertia<T> & rbI)
  {
    M_.template triangularView<Eigen::Lower>() += rbI.M_;
    H_ += rbI.H_;
    I_.template triangularView<Eigen::Lower>() += rbI.I_;
    return *this;
  }

  ABInertia<T> & operator-=(const ABInertia<T> & rbI)
  {
    M_.template triangularView<Eigen::Lower>() -= rbI.M_;
    H_ -= rbI.H_;
    I_.template triangularView<Eigen::Lower>() -= rbI.I_;
    return *this;
  }

  template<typename T2>
  ABInertia<T> operator*(T2 scalar) const
  {
    matrix3_t M, I;
    M.template triangularView<Eigen::Lower>() = scalar * M_;
    I.template triangularView<Eigen::Lower>() = scalar * I_;
    return ABInertia<T>(M, scalar * H_, I);
  }

  /// @return Ia + I
  ABInertia<T> operator+(const RBInertia<T> & rbI) const;

  /// @return Ia * v
  ForceVec<T> operator*(const MotionVec<T> & mv) const;

  /// @see operator*(const MotionVec<T>& mv) const
  template<typename Derived>
  void mul(const Eigen::MatrixBase<Derived> & mv, Eigen::MatrixBase<Derived> & result) const;

  bool operator==(const ABInertia<T> & abI) const
  {
    return M_ == abI.M_ && H_ == abI.H_ && I_ == abI.I_;
  }

  bool operator!=(const ABInertia<T> & abI) const
  {
    return M_ != abI.M_ || H_ != abI.H_ || I_ != abI.I_;
  }

private:
  matrix3_t M_;
  matrix3_t H_;
  matrix3_t I_;
};

template<typename T, typename T2>
inline ABInertia<T> operator*(T2 scalar, const ABInertia<T> & abI)
{
  return abI * scalar;
}

template<typename T>
inline std::ostream & operator<<(std::ostream & out, const ABInertia<T> & abI)
{
  out << abI.matrix();
  return out;
}

} // namespace sva
