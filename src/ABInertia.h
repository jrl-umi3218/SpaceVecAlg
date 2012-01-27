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

class ABInertia
{
public:
	ABInertia():
		M_(),
		H_(),
		I_()
	{}

	ABInertia(const Matrix3d& M, const Matrix3d& H, const Matrix3d& I):
		M_(),
		H_(H),
		I_()
	{
		M_.triangularView<Lower>() = M;
		I_.triangularView<Lower>() = I;
	}

	ABInertia(const TriangularView<Matrix3d, Lower>& ltM,
						const Matrix3d& H,
						const TriangularView<Matrix3d, Lower>& ltI):
		M_(ltM),
		H_(H),
		I_(ltI)
	{}

	// Accessor
	const Matrix3d& lowerTriangularMassMatrix() const
	{
		return M_;
	}

	Matrix3d massMatrix() const
	{
		Matrix3d M;
		M.triangularView<Upper>() = M_;
		M.triangularView<StrictlyLower>() = M_;
		return M;
	}

	const Matrix3d& gInertia() const
	{
		return H_;
	}

	const Matrix3d& lowerTriangularInertia() const
	{
		return I_;
	}

	Matrix3d inertia() const
	{
		Matrix3d I;
		I.triangularView<Upper>() = I_;
		I.triangularView<StrictlyLower>() = I_;
		return I;
	}

	Matrix6d matrix() const
	{
		Matrix6d m;
		m << inertia(), H_,
				 H_.transpose(), massMatrix();
		return m;
	}

	// Operators
	ABInertia operator+(const ABInertia& rbI) const
	{
		Matrix3d M, I;
		M.triangularView<Lower>() = M_ + rbI.M_;
		I.triangularView<Lower>() = I_ + rbI.I_;
		return ABInertia(M, H_ + rbI.H_, I);
	}

	friend ABInertia operator*(double scalar, const ABInertia& rbI);

private:
	Matrix3d M_;
	Matrix3d H_;
	Matrix3d I_;
};

inline ABInertia operator*(double scalar, const ABInertia& rbI)
{
	Matrix3d M, I;
	M.triangularView<Lower>() = scalar*rbI.M_;
	I.triangularView<Lower>() = scalar*rbI.I_;
	return ABInertia(M, scalar*rbI.H_, I);
}

}
