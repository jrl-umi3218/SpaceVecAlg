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

namespace sva
{

using namespace Eigen;

/**
	* Spatial Articulated Body Inertia compact representation.
	* See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
	*/
class ABInertia
{
public:
	ABInertia():
		M_(),
		H_(),
		I_()
	{}

	/**
		* @param M Mass matrix.
		* @param H Generalized inertia matrix.
		* @param I Inertia matrix.
		*/
	ABInertia(const Matrix3d& M, const Matrix3d& H, const Matrix3d& I):
		M_(Matrix3d::Zero()),
		H_(H),
		I_(Matrix3d::Zero())
	{
		M_.triangularView<Lower>() = M;
		I_.triangularView<Lower>() = I;
	}

	/**
		* @param M Lower triangular view of Mass matrix.
		* @param H Generalized inertia matrix.
		* @param I Lower triangular view Inertia matrix.
		*/
	ABInertia(const TriangularView<Matrix3d, Lower>& ltM,
						const Matrix3d& H,
						const TriangularView<Matrix3d, Lower>& ltI):
		M_(ltM),
		H_(H),
		I_(ltI)
	{}

	// Accessor
	/// @return Mass matrix with a zero upper part.
	const Matrix3d& lowerTriangularMassMatrix() const
	{
		return M_;
	}

	/// @return Mass matrix.
	Matrix3d massMatrix() const
	{
		Matrix3d M;
		M.triangularView<Upper>() = M_.transpose();
		M.triangularView<StrictlyLower>() = M_;
		return M;
	}

	/// @return Generalized inertia matrix.
	const Matrix3d& gInertia() const
	{
		return H_;
	}

	/// @return Inertia matrix with a zero upper part.
	const Matrix3d& lowerTriangularInertia() const
	{
		return I_;
	}

	/// @return Inertia matrix.
	Matrix3d inertia() const
	{
		Matrix3d I;
		I.triangularView<Upper>() = I_.transpose();
		I.triangularView<StrictlyLower>() = I_;
		return I;
	}

	/// @retrun Non compact spatial articulated body inertia matrix.
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

	ABInertia operator*(double scalar) const
	{
		Matrix3d M, I;
		M.triangularView<Lower>() = scalar*M_;
		I.triangularView<Lower>() = scalar*I_;
		return ABInertia(M, scalar*H_, I);
	}

	/// @return Ia + I
	ABInertia operator+(const RBInertia& rbI) const;

	/// @return Ia * v
	ForceVec operator*(const MotionVec& mv) const;

	bool operator==(const ABInertia& abI) const
	{
		return M_ == abI.M_ && H_ == abI.H_ && I_ == abI.I_;
	}

	bool operator!=(const ABInertia& abI) const
	{
		return M_ != abI.M_ || H_ != abI.H_ || I_ != abI.I_;
	}

private:
	Matrix3d M_;
	Matrix3d H_;
	Matrix3d I_;
};

inline ABInertia operator*(double scalar, const ABInertia& abI)
{
	return abI*scalar;
}

inline std::ostream& operator<<(std::ostream& out, const ABInertia& abI)
{
	out << abI.matrix();
	return out;
}

} // namespace sva
