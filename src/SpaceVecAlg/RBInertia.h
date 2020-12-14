/*
 * Copyright 2012-2020 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "EigenUtility.h"
#include "fwd.h"

namespace sva
{

template<typename T>
Eigen::Matrix3<T> inertiaToOrigin(const Eigen::Matrix3<T> & inertia,
                                  T mass,
                                  const Eigen::Vector3<T> & com,
                                  const Eigen::Matrix3<T> & rotation);

/**
 * Spatial Rigid Body Inertia compact representation.
 * See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
 */
template<typename T>
class RBInertia
{
public:
  typedef Eigen::Vector3<T> vector3_t;
  typedef Eigen::Matrix3<T> matrix3_t;
  typedef Eigen::Matrix6<T> matrix6_t;

public:
  RBInertia() : m_(), h_(), I_() {}

  /**
   * @param m Mass.
   * @param h Spatial momentum.
   * @param I Inertia matrix at body origin.
   */
  RBInertia(T m, const vector3_t & h, const matrix3_t & I) : m_(m), h_(h), I_(matrix3_t::Zero())
  {
    I_.template triangularView<Eigen::Lower>() = I;
  }

  /**
   * @param m Mass.
   * @param h Spatial momentum.
   * @param I Eigen::Lower triangular view of Inertia matrix at body origin.
   */
  RBInertia(T m, const vector3_t & h, const Eigen::TriangularView<matrix3_t, Eigen::Lower> & ltI)
  : m_(m), h_(h), I_(ltI)
  {
  }

  // Accessor
  /// @return Mass.
  T mass() const
  {
    return m_;
  }

  /// @return Spatial momentum.
  const vector3_t & momentum() const
  {
    return h_;
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

  /// @retrun Non compact spatial rigid body inertia matrix.
  matrix6_t matrix() const
  {
    matrix6_t m;
    matrix3_t hCross = vector3ToCrossMatrix(h_);
    m << inertia(), hCross, hCross.transpose(), matrix3_t::Identity() * m_;
    return m;
  }

  template<typename T2>
  RBInertia<T2> cast() const
  {
    return RBInertia<T2>(T2(m_), h_.template cast<T2>(), I_.template cast<T2>());
  }

  // Operators
  RBInertia<T> operator+(const RBInertia<T> & rbI) const
  {
    matrix3_t I;
    I.template triangularView<Eigen::Lower>() = I_ + rbI.I_;
    return RBInertia<T>(m_ + rbI.m_, h_ + rbI.h_, I);
  }

  RBInertia<T> operator-(const RBInertia<T> & rbI) const
  {
    matrix3_t I;
    I.template triangularView<Eigen::Lower>() = I_ - rbI.I_;
    return RBInertia<T>(m_ - rbI.m_, h_ - rbI.h_, I);
  }

  RBInertia<T> operator-() const
  {
    return RBInertia<T>(-m_, -h_, -I_);
  }

  RBInertia<T> & operator+=(const RBInertia<T> & rbI)
  {
    I_.template triangularView<Eigen::Lower>() += rbI.I_;
    m_ += rbI.m_;
    h_ += rbI.h_;
    return *this;
  }

  RBInertia<T> & operator-=(const RBInertia<T> & rbI)
  {
    I_.template triangularView<Eigen::Lower>() -= rbI.I_;
    m_ -= rbI.m_;
    h_ -= rbI.h_;
    return *this;
  }

  template<typename T2>
  RBInertia<T> operator*(T2 scalar) const
  {
    matrix3_t I;
    I.template triangularView<Eigen::Lower>() = scalar * I_;
    return RBInertia<T>(scalar * m_, scalar * h_, I);
  }

  /// @return I*v
  ForceVec<T> operator*(const MotionVec<T> & mv) const;

  /// @see operator*(const MotionVec<T>& mv) const
  template<typename Derived>
  void mul(const Eigen::MatrixBase<Derived> & mv, Eigen::MatrixBase<Derived> & result) const;

  bool operator==(const RBInertia<T> & rbI) const
  {
    return m_ == rbI.m_ && h_ == rbI.h_ && I_ == rbI.I_;
  }

  bool operator!=(const RBInertia<T> & rbI) const
  {
    return m_ != rbI.m_ || h_ != rbI.h_ || I_ != rbI.I_;
  }

private:
  T m_;
  vector3_t h_;
  matrix3_t I_;
};

template<typename T, typename T2>
inline RBInertia<T> operator*(T2 scalar, const RBInertia<T> & rbI)
{
  return rbI * scalar;
}

template<typename T>
inline std::ostream & operator<<(std::ostream & out, const RBInertia<T> & rbI)
{
  out << rbI.matrix();
  return out;
}

template<typename T>
Eigen::Matrix3<T> inertiaToOrigin(const Eigen::Matrix3<T> & inertia,
                                  T mass,
                                  const Eigen::Vector3<T> & com,
                                  const Eigen::Matrix3<T> & rotation)
{
  Eigen::Matrix3<T> trans =
      Eigen::vector3ToCrossMatrix<T>(mass * com) * Eigen::vector3ToCrossMatrix<T>(com).transpose();
  return rotation * (inertia + trans) * rotation.transpose();
}

} // namespace sva
