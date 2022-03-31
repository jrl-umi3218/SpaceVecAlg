/*
 * Copyright 2012-2020 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "EigenTypedef.h"
#include "fwd.h"

namespace sva
{

/**
 * Create a rotation matrix about the X axis.
 * The rotation is exprimed in successor frame.
 * @param theta rotation in radian.
 */
template<typename T>
Eigen::Matrix3<T> RotX(T theta);

/**
 * Create a rotation matrix about the Y axis.
 * The rotation is exprimed in successor frame.
 * @param theta rotation in radian.
 */
template<typename T>
Eigen::Matrix3<T> RotY(T theta);

/**
 * Create a rotation matrix about the Z axis.
 * The rotation is exprimed in successor frame.
 * @param theta rotation in radian.
 */
template<typename T>
Eigen::Matrix3<T> RotZ(T theta);

/**
 * Compute the 3D rotation error between two matrix E_a_b and E_a_c in the 'a' frame.
 * This method convert the 3D rotation matrix E_b_c into a rotation vector.
 * The matrix E_b_c is computed as follow E_a_c = E_b_c*E_a_b.
 * Then the error is computed with E_b_a*rotationVelocity(E_b_c).
 * @return XYZ rotation in radian.
 */
template<typename T>
Eigen::Vector3<T> rotationError(const Eigen::Matrix3<T> & E_a_b, const Eigen::Matrix3<T> & E_a_c);

/**
 * Compute the 3D rotation vector of the rotation matrix E_a_b in the 'a' frame.
 * If we integrate this rotation vector for 1 second we must
 * have the rotation matrix E_a_b.
 * (see exponential matrix and logarithmic matrix).
 */
template<typename T>
Eigen::Vector3<T> rotationVelocity(const Eigen::Matrix3<T> & E_a_b);

/**
 * Compute the 6D error between two PTransform in the 'a' frame.
 * This method convert the 6D transformation matrix X_b_c into a motion vector.
 * The matrix X_b_c is computed as follow X_a_c = X_b_c*X_a_b.
 * Then the error is computed with PTransform(E_b_a)*transformVelocity(X_b_c).
 * @return XYZ rotation in radian.
 */
template<typename T>
MotionVec<T> transformError(const PTransform<T> & X_a_b, const PTransform<T> & X_a_c);

/**
 * Compute the motion vector of the matrix X_a_b in the 'a' frame.
 * If we integrate this motion vector for 1 second we must
 * have the transformation matrix X_a_b.
 * This function can be see as an implementation of the function XtoV
 * (see Featherstone appendix) but with the use of logarithmic
 * matrix to compute the rotational error.
 */
template<typename T>
MotionVec<T> transformVelocity(const PTransform<T> & X_a_b);

/**
 * Plücker transform compact representation.
 * Use 3D matrix as rotation internal representation.
 * Quaternion are inversed as they must be expressed in successor frame.
 * See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
 */
template<typename T>
class PTransform
{
  typedef Eigen::Vector3<T> vector3_t;
  typedef Eigen::Matrix3<T> matrix3_t;
  typedef Eigen::Matrix6<T> matrix6_t;
  typedef Eigen::Quaternion<T> quaternion_t;

public:
  /// Identity transformation.
  static PTransform<T> Identity()
  {
    return PTransform<T>(matrix3_t::Identity(), vector3_t::Zero());
  }

public:
  // Constructors
  /// Default constructor. Rotation and translation are uninitialized.
  PTransform() : E_(), r_() {}

  /// Copy constructor.
  template<typename T2>
  PTransform(const PTransform<T2> & pt) : E_(pt.rotation().template cast<T>()), r_(pt.translation().template cast<T>())
  {
  }

  /**
   * @param rot Rotation matrix.
   * @param trans Translation vector.
   */
  PTransform(const matrix3_t & rot, const vector3_t & trans) : E_(rot), r_(trans) {}

  /**
   * @param rot Rotation quaternion.
   * @param trans Translation vector.
   */
  PTransform(const quaternion_t & rot, const vector3_t & trans) : E_(rot.matrix()), r_(trans) {}

  /**
   * Rotation only transform.
   * @param rot Rotation quaternion.
   */
  PTransform(const quaternion_t & rot) : E_(rot.matrix()), r_(vector3_t::Zero()) {}

  /**
   * Rotation only transform.
   * @param rot Rotation matrix.
   */
  PTransform(const matrix3_t & rot) : E_(rot), r_(vector3_t::Zero()) {}

  /**
   * Translation only transform.
   * @param trans Translation vector.
   */
  PTransform(const vector3_t & trans) : E_(matrix3_t::Identity()), r_(trans) {}

  // Accessor
  /// @return Rotation matrix.
  const matrix3_t & rotation() const
  {
    return E_;
  }

  /// @return Rotation matrix.
  matrix3_t & rotation()
  {
    return E_;
  }

  /// @return Translation vector.
  const vector3_t & translation() const
  {
    return r_;
  }

  /// @return Translation vector.
  vector3_t & translation()
  {
    return r_;
  }

  /// @return Non compact Plücker transformation matrix.
  matrix6_t matrix() const
  {
    matrix6_t m;
    m << E_, matrix3_t::Zero(), -E_ * vector3ToCrossMatrix(r_), E_;
    return m;
  }

  /// @return Non compact dual Plücker transformation matrix.
  matrix6_t dualMatrix() const
  {
    matrix6_t m;
    m << E_, -E_ * vector3ToCrossMatrix(r_), matrix3_t::Zero(), E_;
    return m;
  }

  template<typename T2>
  PTransform<T2> cast() const
  {
    return PTransform<T2>(E_.template cast<T2>(), r_.template cast<T2>());
  }

  // Operators
  /// @return X*X
  PTransform<T> operator*(const PTransform<T> & pt) const
  {
    return PTransform<T>(E_ * pt.E_, pt.r_ + pt.E_.transpose() * r_);
  }

  /// @return Xv
  MotionVec<T> operator*(const MotionVec<T> & mv) const;

  /// @return compute angular part of @see operator*(const MotionVec<T>& mv) const
  Eigen::Vector3<T> angularMul(const MotionVec<T> & mv) const;
  /// @return compute linear part of @see operator*(const MotionVec<T>& mv) const
  Eigen::Vector3<T> linearMul(const MotionVec<T> & mv) const;
  /// @see operator*(const MotionVec<T>& mv) const;
  template<typename Derived>
  void mul(const Eigen::MatrixBase<Derived> & mv, Eigen::MatrixBase<Derived> & result) const;

  /// @return X^-1 v
  MotionVec<T> invMul(const MotionVec<T> & mv) const;
  /// @return compute angular part of @see invMul
  Eigen::Vector3<T> angularInvMul(const MotionVec<T> & mv) const;
  /// @return compute linear part of @see invMul
  Eigen::Vector3<T> linearInvMul(const MotionVec<T> & mv) const;
  /// @see invMul
  template<typename Derived>
  void invMul(const Eigen::MatrixBase<Derived> & mv, Eigen::MatrixBase<Derived> & result) const;

  /// @return X*v
  ForceVec<T> dualMul(const ForceVec<T> & fv) const;
  /// @return compute couple part of @see dualMul
  Eigen::Vector3<T> coupleDualMul(const ForceVec<T> & fv) const;
  /// @return compute force part of @see dualMul
  Eigen::Vector3<T> forceDualMul(const ForceVec<T> & fv) const;
  /// @see dualMul
  template<typename Derived>
  void dualMul(const Eigen::MatrixBase<Derived> & fv, Eigen::MatrixBase<Derived> & result) const;

  /// @return Xtv
  ForceVec<T> transMul(const ForceVec<T> & fv) const;
  /// @return compute couple part of @see transMul
  Eigen::Vector3<T> coupleTransMul(const ForceVec<T> & fv) const;
  /// @return compute force part of @see transMul
  Eigen::Vector3<T> forceTransMul(const ForceVec<T> & fv) const;
  /// @see transMul
  template<typename Derived>
  void transMul(const Eigen::MatrixBase<Derived> & fv, Eigen::MatrixBase<Derived> & result) const;

  /// @return X*IX^-1
  RBInertia<T> dualMul(const RBInertia<T> & rbI) const;
  /// @return XtIX
  RBInertia<T> transMul(const RBInertia<T> & rbI) const;

  /// @return X*IX^-1
  ABInertia<T> dualMul(const ABInertia<T> & rbI) const;
  /// @return XtIX
  ABInertia<T> transMul(const ABInertia<T> & rbI) const;

  /// @return Inverse Plücker transformation.
  PTransform<T> inv() const
  {
    return PTransform<T>(E_.transpose(), -E_ * r_);
  }

  bool operator==(const PTransform<T> & pt) const
  {
    return E_ == pt.E_ && r_ == pt.r_;
  }

  bool operator!=(const PTransform<T> & pt) const
  {
    return E_ != pt.E_ || r_ != pt.r_;
  }

private:
  matrix3_t E_;
  vector3_t r_;
};

template<typename T>
inline Eigen::Matrix3<T> RotX(T theta)
{
  T s = std::sin(theta), c = std::cos(theta);
  return (Eigen::Matrix3<T>() << 1., 0., 0., 0., c, s, 0., -s, c).finished();
}

template<typename T>
inline Eigen::Matrix3<T> RotY(T theta)
{
  T s = std::sin(theta), c = std::cos(theta);
  return (Eigen::Matrix3<T>() << c, 0., -s, 0., 1., 0., s, 0., c).finished();
}

template<typename T>
inline Eigen::Matrix3<T> RotZ(T theta)
{
  T s = std::sin(theta), c = std::cos(theta);
  return (Eigen::Matrix3<T>() << c, s, 0., -s, c, 0., 0., 0., 1.).finished();
}

template<typename T>
inline Eigen::Vector3<T> rotationError(const Eigen::Matrix3<T> & E_a_b, const Eigen::Matrix3<T> & E_a_c)
{
  Eigen::Matrix3<T> E_b_c = E_a_c * E_a_b.transpose();
  return Eigen::Vector3<T>(E_a_b.transpose() * rotationVelocity(E_b_c));
}

template<typename T>
inline Eigen::Vector3<T> rotationVelocity(const Eigen::Matrix3<T> & E_a_b)
{
  constexpr T eps = std::numeric_limits<T>::epsilon();
  constexpr T sqeps = details::sqrt(eps);
  constexpr T sqsqeps = details::sqrt(sqeps);
  constexpr T pi = 3.1415926535897932;

  T trace = E_a_b(0, 0) + E_a_b(1, 1) + E_a_b(2, 2);
  T acosV = (trace - 1.) * 0.5;
  T theta = std::acos(std::min(std::max(acosV, -1.), 1.));

  Eigen::Vector3<T> w(-E_a_b(2, 1) + E_a_b(1, 2), -E_a_b(0, 2) + E_a_b(2, 0), -E_a_b(1, 0) + E_a_b(0, 1));

  if(1 + trace < sqsqeps) // in case of angle close to pi
  {
    // adapted from https://vision.in.tum.de/_media/members/demmeln/nurlanov2021so3log.pdf, sec2.2
    Eigen::Vector3<T> s = (2 * E_a_b.diagonal() + Eigen::Vector3<T>::Constant(1 - trace)) / (3 - trace);
    Eigen::Vector3<T> tn2 = theta * s.cwiseSqrt();
    // If theta is really close to pi the sign of n is derived from (13)
    // We choose the first non-zero coefficient to be positive arbitrarly
    if(theta > pi - 1e-7)
    {
      if(tn2(0) > 0)
      {
        if(E_a_b(0, 1) + E_a_b(1, 0) < 0)
        {
          tn2(1) = -tn2(1);
        }
        if(E_a_b(0, 2) + E_a_b(2, 0) < 0)
        {
          tn2(2) = -tn2(2);
        }
      }
      else if(tn2(1) > 0) // only if tn2(0) == 0
      {
        if(E_a_b(1, 2) + E_a_b(2, 1) < 0)
        {
          tn2(2) = -tn2(2);
        }
      }
      return tn2;
    }
    // The sign is derived from (9)
    return (w.array() >= 0).select(tn2, -tn2);
  }
  else
  {
    w *= sinc_inv(theta) * 0.5;
  }

  return w;
}

template<typename T>
inline MotionVec<T> transformError(const PTransform<T> & X_a_b, const PTransform<T> & X_a_c)
{
  PTransform<T> X_b_c = X_a_c * X_a_b.inv();
  return PTransform<T>(Eigen::Matrix3<T>(X_a_b.rotation().transpose())) * transformVelocity(X_b_c);
}

template<typename T>
inline MotionVec<T> transformVelocity(const PTransform<T> & X_a_b)
{
  return MotionVec<T>(rotationVelocity(X_a_b.rotation()), X_a_b.translation());
}

// interpolate between transformations, t must be between 0 and 1
template<typename T>
PTransform<T> interpolate(const PTransform<T> & from, const PTransform<T> & to, double t)
{
  Eigen::Quaternion<T> qfrom(from.rotation());
  Eigen::Quaternion<T> qto(to.rotation());
  PTransform<T> result(qfrom.slerp(t, qto), (from.translation() * (1. - t) + to.translation() * t));
  return result;
}

template<typename T>
inline std::ostream & operator<<(std::ostream & out, const PTransform<T> & pt)
{
  out << pt.matrix();
  return out;
}

} // namespace sva
