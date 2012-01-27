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

// includes
// Eigen
#include <Eigen/Core>

// SpaceVecAlg
#include "EigenTypedef.h"

namespace Eigen
{

inline Matrix3d vector3ToCrossMatrix(const Vector3d& vec)
{
	Matrix3d mat;
	mat << 0.,-vec(2),vec(1),
				 vec(2),0.,-vec(0),
				 -vec(1),vec(0),0.;
	return mat;
}

inline Matrix6d vector6ToCrossMatrix(const Vector6d& vec)
{
	Matrix6d mat;
	Matrix3d c13 = vector3ToCrossMatrix(vec.head<3>());
	mat << c13, Matrix3d::Zero(),
				 vector3ToCrossMatrix(vec.tail<3>()), c13;
	return mat;
}

inline Matrix6d vector6ToCrossDualMatrix(const Vector6d& vec)
{
	return -vector6ToCrossMatrix(vec).transpose();
}

}
