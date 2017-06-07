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

/**
	* Convert a 3D vector to a cross product matrix.
	*/
template<typename T>
inline Matrix3<T> vector3ToCrossMatrix(const Vector3<T>& vec)
{
	Matrix3<T> mat;
	mat << 0.,-vec(2),vec(1),
				 vec(2),0.,-vec(0),
				 -vec(1),vec(0),0.;
	return mat;
}

/**
	* Convert a 6D vector to a spatial cross product matrix.
	*/
template<typename T>
inline Matrix6<T> vector6ToCrossMatrix(const Vector6<T>& vec)
{
	Matrix6<T> mat;
	Matrix3<T> c13 = vector3ToCrossMatrix<T>(vec.template head<3>());
	mat << c13, Matrix3d::Zero(),
				 vector3ToCrossMatrix<T>(vec.template tail<3>()), c13;
	return mat;
}

/**
	* Convert a 6D vector to a spatial dual cross product matrix.
	*/
template<typename T>
inline Matrix6<T> vector6ToCrossDualMatrix(const Vector6<T>& vec)
{
	return -vector6ToCrossMatrix<T>(vec).transpose();
}

} // namespace Eigen
