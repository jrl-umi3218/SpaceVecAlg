/*
 * Copyright 2012-2020 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "EigenTypedef.h"
#include "fwd.h"
#include <type_traits>

namespace sva
{

/**
 * Spatial Motion Vector compact representations.
 * See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
 */
template<typename T>
class MotionVec
{
public:
  typedef Eigen::Vector3<T> vector3_t;
  typedef Eigen::Vector6<T> vector6_t;

  friend class PTransform<T>;

public:
  /// Zero motion vector
  static MotionVec<T> Zero()
  {
    return MotionVec<T>(vector3_t::Zero(), vector3_t::Zero());
  }

public:
  MotionVec() : angular_(), linear_() {}

  /**
   * @param vec Spatial motion vector with angular motion in head
   * and linear motion in tail.
   */
  MotionVec(const vector6_t & vec) : angular_(vec.template head<3>()), linear_(vec.template tail<3>()) {}

  /**
   * @param angular Angular motion.
   * @param linear Linear motion.
   */
  MotionVec(const vector3_t & angular, const vector3_t & linear) : angular_(angular), linear_(linear) {}

  // Accessor
  /// @return Angular motion
  const vector3_t & angular() const
  {
    return angular_;
  }

  /// @return Angular motion
  vector3_t & angular()
  {
    return angular_;
  }

  /// @return Linear motion
  const vector3_t & linear() const
  {
    return linear_;
  }

  /// @return Linear motion
  vector3_t & linear()
  {
    return linear_;
  }

  /// @return Non compact spatial motion vector.
  vector6_t vector() const
  {
    return ((vector6_t() << angular_, linear_).finished());
  }

  template<typename T2>
  MotionVec<T2> cast() const
  {
    return MotionVec<T2>(angular_.template cast<T2>(), linear_.template cast<T2>());
  }

  // Operators
  MotionVec<T> operator+(const MotionVec<T> & mv) const
  {
    return MotionVec<T>(angular_ + mv.angular_, linear_ + mv.linear_);
  }

  MotionVec<T> operator-(const MotionVec<T> & mv) const
  {
    return MotionVec<T>(angular_ - mv.angular_, linear_ - mv.linear_);
  }

  MotionVec<T> operator-() const
  {
    return MotionVec<T>(-angular_, -linear_);
  }

  MotionVec<T> & operator+=(const MotionVec<T> & mv)
  {
    angular_ += mv.angular_;
    linear_ += mv.linear_;
    return *this;
  }

  MotionVec<T> & operator-=(const MotionVec<T> & mv)
  {
    angular_ -= mv.angular_;
    linear_ -= mv.linear_;
    return *this;
  }

  template<typename T2, typename std::enable_if<std::is_arithmetic<T2>::value, int>::type = 0>
  MotionVec<T> operator*(T2 scalar) const
  {
    return MotionVec<T>(scalar * angular_, scalar * linear_);
  }

  template<typename T2, typename std::enable_if<std::is_arithmetic<T2>::value, int>::type = 0>
  MotionVec<T> & operator*=(T2 scalar)
  {
    angular_ *= scalar;
    linear_ *= scalar;
    return *this;
  }

  template<typename T2>
  MotionVec<T> & operator/=(T2 scalar)
  {
    angular_ /= scalar;
    linear_ /= scalar;
    return *this;
  }

  template<typename T2>
  MotionVec<T> operator/(T2 scalar) const
  {
    return MotionVec<T>(angular_ / scalar, linear_ / scalar);
  }

  /// @return v x v
  MotionVec<T> cross(const MotionVec<T> & mv2) const;

  /// @see cross
  template<typename Derived>
  void cross(const Eigen::MatrixBase<Derived> & mv2, Eigen::MatrixBase<Derived> & result) const;

  /// @return v x* f
  ForceVec<T> crossDual(const ForceVec<T> & fv2) const;

  /// @see crossDual
  template<typename Derived>
  void crossDual(const Eigen::MatrixBase<Derived> & fv2, Eigen::MatrixBase<Derived> & result) const;

  /// @return v.v
  T dot(const ForceVec<T> & fv2) const;

  bool operator==(const MotionVec<T> & mv) const
  {
    return (angular_ == mv.angular_) && (linear_ == mv.linear_);
  }

  bool operator!=(const MotionVec<T> & mv) const
  {
    return (angular_ != mv.angular_) || (linear_ != mv.linear_);
  }

private:
  vector3_t angular_;
  vector3_t linear_;
};

template<typename T, typename T2>
inline MotionVec<T> operator*(T2 scalar, const MotionVec<T> & mv)
{
  return mv * scalar;
}

template<typename T>
inline std::ostream & operator<<(std::ostream & out, const MotionVec<T> & mv)
{
  out << mv.vector().transpose();
  return out;
}

} // namespace sva
