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
	* Spatial Rigid Body Inertia compact representation.
	* See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
	*/
class RBInertia
{
public:
	RBInertia():
		m_(),
		h_(),
		I_()
	{}

	/**
		* @param m Mass.
		* @param h Spatial momentum.
		* @param I Inertia matrix.
		*/
	RBInertia(double m, const Vector3d& h, const Matrix3d& I):
		m_(m),
		h_(h),
		I_(Matrix3d::Zero())
	{
		I_.triangularView<Lower>() = I;
	}

	/**
		* @param m Mass.
		* @param h Spatial momentum.
		* @param I Lower triangular view of Inertia matrix.
		*/
	RBInertia(double m, const Vector3d& h,
						const TriangularView<Matrix3d, Lower>& ltI):
		m_(m),
		h_(h),
		I_(ltI)
	{}

	// Accessor
	/// @return Mass.
	double mass() const
	{
		return m_;
	}

	/// @return Spatial momentum.
	const Vector3d& momentum() const
	{
		return h_;
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

	/// @retrun Non compact spatial rigid body inertia matrix.
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

	RBInertia operator*(double scalar) const
	{
		Matrix3d I;
		I.triangularView<Lower>() = scalar *I_;
		return RBInertia(scalar * m_, scalar * h_, I);
	}

	/// @return I*v
	ForceVec operator*(const MotionVec& mv);

private:
	double m_;
	Vector3d h_;
	Matrix3d I_;
};

inline RBInertia operator*(double scalar, const RBInertia& rbI)
{
	return rbI*scalar;
}

} // namespace sva
