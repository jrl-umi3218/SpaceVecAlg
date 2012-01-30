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

class PTransform
{
public:
	// Constructor
	PTransform():
		E_(Matrix3d::Identity()),
		r_(Vector3d::Zero())
	{}

	PTransform(const Matrix3d& rot, const Vector3d& trans):
		E_(rot),
		r_(trans)
	{}

	PTransform(const Quaterniond& rot, const Vector3d& trans):
		E_(),
		r_(trans)
	{
		E_ = rot.toRotationMatrix();
	}

	PTransform(const Quaterniond& rot):
		E_(),
		r_(Vector3d::Zero())
	{
		E_ = rot.toRotationMatrix();
	}

	PTransform(const Matrix3d& rot):
		E_(rot),
		r_(Vector3d::Zero())
	{}

	PTransform(const Vector3d& trans):
		E_(Matrix3d::Identity()),
		r_(trans)
	{}

	// Accessor
	const Matrix3d& rotation() const
	{
		return E_;
	}

	Matrix3d& rotation()
	{
		return E_;
	}

	const Vector3d& translation() const
	{
		return r_;
	}

	Vector3d& translation()
	{
		return r_;
	}

	Matrix6d matrix() const
	{
		Matrix6d m;
		m << E_, Matrix3d::Zero(),
				 -E_*vector3ToCrossMatrix(r_), E_;
		return m;
	}

	Matrix6d dualMatrix() const
	{
		Matrix6d m;
		m << E_, -E_*vector3ToCrossMatrix(r_),
				 Matrix3d::Zero(), E_;
		return m;
	}

	// Operators
	PTransform operator*(const PTransform& pt) const
	{
		return PTransform(E_*pt.E_, pt.r_ + pt.E_.transpose()*r_);
	}

	PTransform inv() const
	{
		return PTransform(E_.transpose(), -E_*r_);
	}

private:
	Matrix3d E_;
	Vector3d r_;
};

}
