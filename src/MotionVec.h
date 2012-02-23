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

namespace sva
{

using namespace Eigen;

/**
	* Spatial Motion Vector compact representations.
	* See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
	*/
class MotionVec
{
public:
	MotionVec():
		mv_()
	{}

	/**
		* @param vec Spatial motion vector with angular motion in head
		* and linear motion in tail.
		*/
	MotionVec(const Vector6d& vec):
		mv_(vec)
	{}

	/**
		* @param angular Angular motion.
		* @param linear Linear motion.
		*/
	MotionVec(const Vector3d& angular, const Vector3d& linear):
		mv_((Vector6d() << angular, linear).finished())
	{}

	// Accessor
	/// @return Angular motion part (3 first parameters).
	Vector3d angular() const
	{
		return mv_.head<3>();
	}

	/// @return Linear motion part (3 last parameters).
	Vector3d linear() const
	{
		return mv_.tail<3>();
	}

	/// @return Non compact spatial motion vector.
	const Vector6d& vector() const
	{
		return mv_;
	}

	/// @return Non compact spatial motion vector.
	Vector6d& vector()
	{
		return mv_;
	}

	// Operators
	MotionVec operator+(const MotionVec& mv) const
	{
		return MotionVec(mv_ + mv.mv_);
	}

	MotionVec operator-(const MotionVec& mv) const
	{
		return MotionVec(mv_ - mv.mv_);
	}

	friend MotionVec operator*(double scalar, const MotionVec& fv);

private:
	Vector6d mv_;
};

inline MotionVec operator*(double scalar, const MotionVec& mv)
{
	return MotionVec(scalar * mv.mv_);
}

} // namespace sva
