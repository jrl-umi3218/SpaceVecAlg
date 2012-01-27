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
#include "EigenUtility.h"

namespace sva
{

using namespace Eigen;

class RBInertia
{
public:
	RBInertia():
		m_(0.),
		h_(),
		I_()
	{}

	RBInertia(double m, const Vector3d& h, const Matrix3d& I):
		m_(m),
		h_(h),
		I_(Matrix3d::Zero())
	{
		I_.triangularView<Lower>() = I;
	}

	RBInertia(double m, const Vector3d& h,
						const TriangularView<Matrix3d, Lower>& ltI):
		m_(m),
		h_(h),
		I_(ltI)
	{}

	// Accessor
	double mass() const
	{
		return m_;
	}

	const Vector3d& momentum() const
	{
		return h_;
	}

	const Matrix3d& lowerTriangularInertia() const
	{
		return I_;
	}

	Matrix3d inertia() const
	{
		Matrix3d I;
		I.triangularView<Upper>() = I_.transpose();
		I.triangularView<StrictlyLower>() = I_;
		return I;
	}

	Matrix6d matrix() const
	{
		Matrix6d m;
		Matrix3d hCross = vector3ToCrossMatrix(h_);
		m << inertia(), hCross,
				 hCross.transpose(), Matrix3d::Identity()*m_;
		return m;
	}

	// Operators
	RBInertia operator+(const RBInertia& rbI) const
	{
		Matrix3d I;
		I.triangularView<Lower>() = I_ + rbI.I_;
		return RBInertia(m_ + rbI.m_,
										 h_ + rbI.h_,
										 I);
	}

	friend RBInertia operator*(double scalar, const RBInertia& rbI);

private:
	double m_;
	Vector3d h_;
	Matrix3d I_;
};

inline RBInertia operator*(double scalar, const RBInertia& rbI)
{
	Matrix3d I;
	I.triangularView<Lower>() = scalar * rbI.I_;
	return RBInertia(scalar * rbI.m_, scalar * rbI.h_, I);
}

}
