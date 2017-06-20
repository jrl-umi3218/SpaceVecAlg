// Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
//
// This file is part of SpaceVecAlg.
//
// SpaceVecAlg is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SpaceVecAlg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with SpaceVecAlg.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include "PTransform.h"

namespace sva
{

/**
* \addtogroup Conversions Convert to and from sva types
*
* @{
*/

/**
* @brief Generic functions to convert to/from sva types
*
* Each function in the conversions namespace allows to convert to and from
* another type. An important point to note is that sva internally uses a
* left-handed convention. Hence, all functions will assume that the non-sva
* objects use a right-handed convention and convert accordingly. This is
* overridable by the "rightHandedness" argument.
*/
namespace conversions
{
  //! Alias for right handedness (default)
  constexpr bool RightHanded = true;
  //! Alias for left handedness
  constexpr bool LeftHanded = false;

  /**
   * @brief Convert an homogeneous matrix into a Plucker Transform
   *
   * @return Plucker Transform equivalent to the input homogeneous matrix
   * @param m A 4x4 Eigen Matrix
   * @param RightHandedness Handedness of the input homogeneous matrix.
   * Defaults to right handedness.
   */
  template<typename Derived>
  PTransform<typename Derived::Scalar> fromHomogeneous(const Eigen::MatrixBase<Derived>& m,
                                            bool rightHandedness = RightHanded)
  {
    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived, 4, 4);
    if(rightHandedness)
    {
      return PTransform<typename Derived::Scalar>(m.template block<3,3>(0,0).transpose(), m.template block<3,1>(0,3));
    }
    else
    {
      return PTransform<typename Derived::Scalar>(m.template block<3,3>(0,0), m.template block<3,1>(0,3));
    }
  }

  /**
   * @brief Convert a Plucker Transform into an homogeneous matrix
   *
   * @return An Eigen 4x4 homogeneous matrix equivalent to the input Plucker Transform
   * @param pt Plucker Transform
   * @param RightHandedness Handedness of the output homogeneous matrix.
   * Defaults to right handedness.
   */
  template<typename T>
  Eigen::Matrix<T, 4, 4> toHomogeneous(const PTransform<T>& pt,
                                       bool rightHandedness = RightHanded)
  {
    Eigen::Matrix<T, 4, 4> res = Eigen::Matrix<T, 4, 4>::Zero();
    res(3, 3) = 1.0;

    res.template block<3, 1>(0, 3) = pt.translation();

    if(rightHandedness)
    {
      res.template block<3, 3>(0, 0) = pt.rotation().transpose();
    }
    else
    {
      res.template block<3, 3>(0, 0) = pt.rotation();
    }
    return res;
  }

  //! Define an Eigen::Affine3<T>
  template<typename T>
  using affine3_t = Eigen::Transform<T, 3, Eigen::TransformTraits::Affine>;

  /**
   * @brief Convert an Eigen::Affine3<T> into a Plucker Transform
   *
   * @return A Plucker Transform of the same type as the input Transform
   * @param a Input Eigen Transform. Only rotation and translation components
   * are considered.
   * @param RightHandedness Handedness of the input Eigen Transform.
   * Defaults to right handedness.
   */

  template<typename T>
  PTransform<T> fromAffine(const affine3_t<T>& a,
                           bool rightHandedness = RightHanded)
  {
    if(rightHandedness)
    {
      return PTransform<T>(a.rotation().transpose(), a.translation());
    }
    else
    {
      return PTransform<T>(a.rotation(), a.translation());
    }
  }

  /**
   * @brief Convert a Plucker Transform into an Eigen::Affine3<T>
   *
   * @return An Eigen::Transform of dimension 3, and same scalar type as the input PTransform.
   * @param pt Input Plucker Transform.
   * @param RightHandedness Handedness of the output Eigen Transform.
   * Defaults to right handedness.
   */

  template<typename T>
  affine3_t<T> toAffine(const PTransform<T>& pt,
                        bool rightHandedness = RightHanded)
  {
    affine3_t<T> ret;
    ret.setIdentity();
    ret.translation() = pt.translation();
    if(rightHandedness)
    {
      ret.matrix().template block<3, 3>(0, 0) = pt.rotation().transpose();
    }
    else
    {
      ret.matrix().template block<3, 3>(0, 0) = pt.rotation();
    }
    return ret;
  }
}

/** @} */

}
