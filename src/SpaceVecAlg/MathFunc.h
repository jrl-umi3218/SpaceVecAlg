/*
 * Copyright 2012-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <cmath>
#include <limits>

namespace sva
{

namespace details
{
template<typename T>
T constexpr sqrtNewtonRaphson(T x, T curr, T prev)
{
  return curr == prev ? curr : sqrtNewtonRaphson(x, static_cast<T>(0.5) * (curr + x / curr), curr);
}

/**
 * Constexpr version of the square root
 * Return value:
 *   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
 *   - Otherwise, returns NaN
 * Copied from https://stackoverflow.com/a/34134071
 *
 * \note Should not be used for value x such that x+1==x or x+1==1 in floating point arithmetic.
 * \internal This limitation could be alleviated by giving a better initial guess. Finding the closest 100^n to x and
 * starting with 10^n could be a way to go.
 */
template<typename T>
T constexpr sqrt(T x)
{
  return x >= static_cast<T>(0) && x < std::numeric_limits<T>::infinity() ? sqrtNewtonRaphson(x, x, static_cast<T>(0))
                                                                          : std::numeric_limits<T>::quiet_NaN();
}

template<typename T>
bool constexpr eq(T x, T y)
{
  return x < y ? (y - x) <= y * std::numeric_limits<T>::epsilon() : (x - y) <= x * std::numeric_limits<T>::epsilon();
}

template<typename T>
T constexpr cbrtNewtonRaphson(T x, T curr, T prev)
{
  return eq(curr, prev)
             ? curr
             : cbrtNewtonRaphson(x, (static_cast<T>(2) * curr + x / (curr * curr)) / static_cast<T>(3), curr);
}

template<typename T>
T constexpr cbrtSub(T x)
{
  return x < std::numeric_limits<T>::infinity() ? cbrtNewtonRaphson(x, sqrt(x), static_cast<T>(0))
                                                : std::numeric_limits<T>::quiet_NaN();
}

/**
 * Constexpr version of the cubic root
 * Return value:
 *   - For a finite of "x", returns an approximation for the cubic root of "x"
 *   - Otherwise, returns NaN
 * \note Should not be used for value x such that x+1==x or x+1==1 in floating point arithmetic.
 * \internal The limitation derives from the one on sqrt.
 */
template<typename T>
T constexpr cbrt(T x)
{
  return x >= static_cast<T>(0) ? cbrtSub(x) : -cbrtSub(-x);
}

/** Compute the value \f$ \frac{1}{x^2} - \frac{1+\cos(x)}{2 x \sin(x)} \f$.
 */
template<typename T>
inline T SO3JacF2(const T & x)
{
  // Taylor expansion at 0 is 1/360 * (1 + x^2/60 + x^4/2520 + x^6/100800 + x^8/3991680)
  using details::cbrt;
  using details::sqrt;
  constexpr T ulp = std::numeric_limits<T>::epsilon();
  constexpr T taylor_2_bound = sqrt(static_cast<T>(60) * ulp);
  constexpr T taylor_4_bound = sqrt(sqrt(static_cast<T>(2520) * ulp));
  constexpr T taylor_6_bound = sqrt(cbrt(static_cast<T>(100800) * ulp));
  constexpr T taylor_8_bound = sqrt(sqrt(sqrt(static_cast<T>(3991680) * ulp)));

  double absx = std::abs(x);
  if(absx >= taylor_8_bound)
  {
    return static_cast<T>(1) / (x * x) - (static_cast<T>(1) + std::cos(x)) / (static_cast<T>(2) * x * std::sin(x));
  }
  else
  {
    // approximation by taylor series in x at 0 up to order 0
    T result = static_cast<T>(1);

    if(absx >= taylor_2_bound)
    {
      T x2 = x * x;
      // approximation by taylor series in x at 0 up to order 2
      result += x2 / static_cast<T>(60);

      if(absx >= taylor_4_bound)
      {
        T x4 = x2 * x2;
        // approximation by taylor series in x at 0 up to order 4
        result += x4 / static_cast<T>(2520);
        if(absx >= taylor_6_bound)
        {
          // approximation by taylor series in x at 0 up to order 6
          result += (x2 * x4) / static_cast<T>(100800);
        }
      }
    }

    return result / static_cast<T>(12);
  }
}

/** Compute the value \f$ \frac{x + \sin(x)}{2x^2(1-\cos(x)} - \frac{2}{x^3}\f$
 * which is the derivative of \f$ \frac{1}{x^2} - \frac{1+\cos(x)}{2 x \sin(x)} \f$.
 */
template<typename T>
inline T dSO3JacF2(const T & x)
{
  // Taylor expansion at 0 is x/360 * (1 + x^2/21 + x^4/560 + x^6/16632 + (691*x^8)/363242880) + x^10/17297280 +
  // (3617*x^12)/2117187072000 We need to go this far in the Taylor expansion, because for numbers around
  // taylor_8_bound, the analytical formula gives a relative error of about 1e-6 with doubles, which is too high. Going
  // up to 13th order gives us a max error of about 3e-10.
  using details::cbrt;
  using details::sqrt;
  constexpr T ulp = std::numeric_limits<T>::epsilon();
  constexpr T taylor_2_bound = sqrt(static_cast<T>(21) * ulp);
  constexpr T taylor_4_bound = sqrt(sqrt(static_cast<T>(560) * ulp));
  constexpr T taylor_6_bound = sqrt(cbrt(static_cast<T>(16632) * ulp));
  constexpr T taylor_8_bound = sqrt(sqrt(sqrt(static_cast<T>(363242880) * ulp / static_cast<T>(691))));
  constexpr T taylor_12_bound = sqrt(sqrt(cbrt(static_cast<T>(2117187072000) * ulp / static_cast<T>(3617))));

  double absx = std::abs(x);
  if(absx >= taylor_12_bound)
  {
    T x2 = x * x;
    return (x + std::sin(x)) / (static_cast<T>(2) * x2 * (1 - std::cos(x))) - static_cast<T>(2) / (x2 * x);
  }
  else
  {
    // approximation by taylor series in x at 0 up to order 1
    T result = static_cast<T>(1);

    if(absx >= taylor_2_bound)
    {
      T x2 = x * x;
      // approximation by taylor series in x at 0 up to order 3
      result += x2 / static_cast<T>(21);

      if(absx >= taylor_4_bound)
      {
        T x4 = x2 * x2;
        // approximation by taylor series in x at 0 up to order 5
        result += x4 / static_cast<T>(560);
        if(absx >= taylor_6_bound)
        {
          // approximation by taylor series in x at 0 up to order 7
          result += (x2 * x4) / static_cast<T>(16632);

          if(absx >= taylor_8_bound)
          {
            T x8 = x4 * x4;
            // approximation by taylor series in x at 0 up to order 11
            result += (static_cast<T>(691) * x8) / static_cast<T>(363242880) + (x2 * x8) / static_cast<T>(17297280);
          }
        }
      }
    }

    return x * result / static_cast<T>(360);
  }
}

} // namespace details

/** sinus cardinal: sin(x)/x
 * Code adapted from boost::math::detail::sinc
 */
template<typename T>
T sinc(const T x)
{
  constexpr T taylor_0_bound = std::numeric_limits<double>::epsilon();
  constexpr T taylor_2_bound = details::sqrt(taylor_0_bound);
  constexpr T taylor_n_bound = details::sqrt(taylor_2_bound);

  if(std::abs(x) >= taylor_n_bound)
  {
    return (std::sin(x) / x);
  }
  else
  {
    // approximation by taylor series in x at 0 up to order 0
    T result = static_cast<T>(1);

    if(std::abs(x) >= taylor_0_bound)
    {
      T x2 = x * x;

      // approximation by taylor series in x at 0 up to order 2
      result -= x2 / static_cast<T>(6);

      if(std::abs(x) >= taylor_2_bound)
      {
        // approximation by taylor series in x at 0 up to order 4
        result += (x2 * x2) / static_cast<T>(120);
      }
    }

    return (result);
  }
}

/**
 * Compute 1/sinc(x).
 * This code is inspired by boost/math/special_functions/sinc.hpp.
 */
template<typename T>
T sinc_inv(const T x)
{
  constexpr T taylor_0_bound = std::numeric_limits<T>::epsilon();
  constexpr T taylor_2_bound = details::sqrt(taylor_0_bound);
  constexpr T taylor_n_bound = details::sqrt(taylor_2_bound);

  // We use the 4th order taylor series around 0 of x/sin(x) to compute
  // this function:
  //
  //     x^2  7x^4
  // 1 + ── + ──── + O(x^6)
  //     6    360
  // this approximation is valid around 0.
  // if x is far from 0, our approximation is not valid
  // since x^6 becomes non negligible we use the normal computation of the function
  // (i.e. taylor_2_bound^6 + taylor_0_bound == taylor_0_bound but
  //       taylor_n_bound^6 + taylor_0_bound != taylor_0).

  if(std::abs(x) >= taylor_n_bound)
  {
    return (x / std::sin(x));
  }
  else
  {
    // x is below taylor_n_bound so we don't care about the 6th order term of
    // the taylor series.
    // We set the 0 order term.
    T result = static_cast<T>(1);

    if(std::abs(x) >= taylor_0_bound)
    {
      // x is above the machine epsilon so x^2 is meaningful.
      T x2 = x * x;
      result += x2 / static_cast<T>(6);

      if(std::abs(x) >= taylor_2_bound)
      {
        // x is above the machine sqrt(epsilon) so x^4 is meaningful.
        result += static_cast<T>(7) * (x2 * x2) / static_cast<T>(360);
      }
    }

    return (result);
  }
}

/** Inverse of the right Jacobian of SO(3) as defined in
 * "A micro Lie theory for state estimation in robotics" by Solà et al. (see in particular eq. 144)
 *
 * Inverse of left Jacobian is simply obtained by transposing (see eq. 147)
 *
 * sva::rotationError(X,Y) is such that is derivative wrt X is \f$ J_r^{-1} \f$ as given by this function and
 * its derivative wrt Y is \f$ -J_r^{-T} \f$.
 *
 * \param u Point of so(3) at which to compute the matrix.
 */
template<typename T>
Eigen::Matrix3<T> SO3RightJacInv(const Eigen::Vector3<T> & u)
{
  auto nu = u.norm();

  Eigen::Matrix3<T> C2;
  auto f2 = details::SO3JacF2(nu);
  auto xx = f2 * u.x() * u.x();
  auto xy = f2 * u.x() * u.y();
  auto xz = f2 * u.x() * u.z();
  auto yy = f2 * u.y() * u.y();
  auto yz = f2 * u.y() * u.z();
  auto zz = f2 * u.z() * u.z();
  // clang-format off
  C2 << -yy - zz,    xy   ,    xz,
           xy   , -xx - zz,    yz,
           xz   ,    yz   , -xx - yy;
  // clang-format on

  Eigen::Matrix3<T> C = vector3ToCrossMatrix(Eigen::Vector3<T>(u / 2));

  return Eigen::Matrix3<T>::Identity() + C + C2;
}

/** Derivative w.r.t. time of Inverse of the right Jacobian of SO(3)
 *
 * \param u Point of so(3) at which to compute the matrix.
 * \param du Derivative w.r.t. time of \p u.
 */
template<typename T>
Eigen::Matrix3<T> SO3RightJacInvDot(const Eigen::Vector3<T> & u, const Eigen::Vector3<T> & du)
{
  auto nu = u.norm();
  Eigen::Matrix3<T> C2;
  Eigen::Matrix3<T> C3;
  auto f2 = details::SO3JacF2(nu);
  auto df2 = (u.dot(du) / nu) * details::dSO3JacF2(nu);
  auto xx = df2 * u.x() * u.x();
  auto xy = df2 * u.x() * u.y();
  auto xz = df2 * u.x() * u.z();
  auto yy = df2 * u.y() * u.y();
  auto yz = df2 * u.y() * u.z();
  auto zz = df2 * u.z() * u.z();
  auto dxx = 2 * f2 * du.x() * u.x();
  auto dyy = 2 * f2 * du.y() * u.y();
  auto dzz = 2 * f2 * du.z() * u.z();
  auto c12 = f2 * (du.y() * u.x() + du.x() * u.y());
  auto c13 = f2 * (du.z() * u.x() + du.x() * u.z());
  auto c23 = f2 * (du.z() * u.y() + du.y() * u.z());
  // clang-format off
  C2 << -yy - zz,    xy   ,    xz,
           xy   , -xx - zz,    yz,
           xz   ,    yz   , -xx - yy;
  C3 << -dyy - dzz,    c12    ,     c13,
           c12    , -dxx - dzz,     c23,
           c13    ,    c23    , -dxx - dyy;
  // clang-format on

  Eigen::Matrix3<T> C = vector3ToCrossMatrix(Eigen::Vector3<T>(du / 2));

  return C + C2 + C3;
}

} // namespace sva
