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

namespace Eigen
{

template<typename T>
using Vector6 = Matrix<T, 6, 1>;
template<typename T>
using Matrix6 = Matrix<T, 6, 6>;

template<typename T>
using Vector3 = Matrix<T, 3, 1>;
template<typename T>
using Matrix3 = Matrix<T, 3, 3>;

typedef Vector6<double> Vector6d;
typedef Matrix6<double> Matrix6d;

} // namespace Eigen
