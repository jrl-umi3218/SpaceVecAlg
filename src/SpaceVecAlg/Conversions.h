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

}

/** @} */

}
