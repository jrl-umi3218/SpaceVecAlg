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
#include <Eigen/Geometry>

// SpaceVecAlg
#include "EigenTypedef.h"
#include "EigenUtility.h"

namespace sva
{

using namespace Eigen;

/**
	* Plücker transform compact representation.
	* Use 3D matrix as rotation internal representation.
	* See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
	*/
class PTransform
{
public:
	// Constructor
	/// Identity transformation.
	PTransform():
		E_(Matrix3d::Identity()),
		r_(Vector3d::Zero())
	{}

	/**
		* @param rot Rotation matrix.
		* @param trans Translation vector.
		*/
	PTransform(const Matrix3d& rot, const Vector3d& trans):
		E_(rot),
		r_(trans)
	{}

	/**
		* @param rot Rotation quaternion.
		* @param trans Translation vector.
		*/
	PTransform(const Quaterniond& rot, const Vector3d& trans):
		E_(),
		r_(trans)
	{
		E_ = rot.toRotationMatrix();
	}

	/**
		* Rotation only transform.
		* @param rot Rotation quaternion.
		*/
	PTransform(const Quaterniond& rot):
		E_(),
		r_(Vector3d::Zero())
	{
		E_ = rot.toRotationMatrix();
	}

	/**
		* Rotation only transform.
		* @param rot Rotation matrix.
		*/
	PTransform(const Matrix3d& rot):
		E_(rot),
		r_(Vector3d::Zero())
	{}

	/**
		* Translation only transform.
		* @param trans Translation vector.
		*/
	PTransform(const Vector3d& trans):
		E_(Matrix3d::Identity()),
		r_(trans)
	{}

	// Accessor
	/// @return Rotation matrix.
	const Matrix3d& rotation() const
	{
		return E_;
	}

	/// @return Rotation matrix.
	Matrix3d& rotation()
	{
		return E_;
	}

	/// @return Translation vector.
	const Vector3d& translation() const
	{
		return r_;
	}

	/// @return Translation vector.
	Vector3d& translation()
	{
		return r_;
	}

	/// @return Non compact Plücker transformation matrix.
	Matrix6d matrix() const
	{
		Matrix6d m;
		m << E_, Matrix3d::Zero(),
				 -E_*vector3ToCrossMatrix(r_), E_;
		return m;
	}

	/// @return Non compact dual Plücker transformation matrix.
	Matrix6d dualMatrix() const
	{
		Matrix6d m;
		m << E_, -E_*vector3ToCrossMatrix(r_),
				 Matrix3d::Zero(), E_;
		return m;
	}

	// Operators
	/// @return X*X
	PTransform operator*(const PTransform& pt) const
	{
		return PTransform(E_*pt.E_, pt.r_ + pt.E_.transpose()*r_);
	}

	/// @return Xv
	MotionVec operator*(const MotionVec& mv);
	/// @return X⁻¹v
	MotionVec invMul(const MotionVec& mv);

	/// @return X*v
	ForceVec dualMul(const ForceVec& fv);
	/// @return Xtv
	ForceVec transMul(const ForceVec& fv);

	/// @return X*IX⁻¹
	RBInertia dualMul(const RBInertia& rbI);
	/// @return XtIX
	RBInertia transMul(const RBInertia& rbI);

	/// @return X*IX⁻¹
	ABInertia dualMul(const ABInertia& rbI);
	/// @return XtIX
	ABInertia transMul(const ABInertia& rbI);

	/// @return Inverse Plücker transformation.
	PTransform inv() const
	{
		return PTransform(E_.transpose(), -E_*r_);
	}

private:
	Matrix3d E_;
	Vector3d r_;
};

} // namespace sva
