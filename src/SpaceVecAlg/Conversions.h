// Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
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

namespace conversions
{
  constexpr bool RightHanded = true;
  constexpr bool LeftHanded = false;

  template<typename T>
  PTransform<T> fromHomogeneous(const Eigen::Ref<Eigen::Matrix<T, 4, 4>>& m,
                                bool rightHandedness = RightHanded)
  {
    if(rightHandedness)
    {
      return PTransform<T>(m.block<3,3>(0,0).transpose(), m.block<3,1>(0,3));
    }
    else
    {
      return PTransform<T>(m.block<3,3>(0,0), m.block<3,1>(0,3));
    }
  }

  template<typename T>
  Eigen::Matrix<T, 4, 4> toHomogeneous(const PTransform<T>& pt,
                                       bool rightHandedness = RightHanded)
  {
    Eigen::Matrix<T, 4, 4> res = Eigen::Matrix<T, 4, 4>::Identity();
    res.block<3, 1>(0, 3) = pt.translation();

    if(rightHandedness)
    {
      res.block<3, 3>(0, 0) = pt.rotation().transpose();
    }
    else
    {
      res.block<3, 3>(0, 0) = pt.rotation();
    }
    return res;
  }

}

}
