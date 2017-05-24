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

namespace sva
{

using namespace Eigen;

/**
	* Spatial Articulated Body Inertia compact representation.
	* See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
	*/
template<typename T>
class ABInertia
{
public:
	typedef Vector3<T> vector3_t;
	typedef Matrix3<T> matrix3_t;
	typedef Matrix6<T> matrix6_t;

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
	ABInertia(const matrix3_t& M, const matrix3_t& H, const matrix3_t& I):
		M_(matrix3_t::Zero()),
		H_(H),
		I_(matrix3_t::Zero())
	{
		M_.template triangularView<Lower>() = M;
		I_.template triangularView<Lower>() = I;
	}

	/**
		* @param M Lower triangular view of Mass matrix.
		* @param H Generalized inertia matrix.
		* @param I Lower triangular view Inertia matrix.
		*/
	ABInertia(const TriangularView<matrix3_t, Lower>& ltM,
						const matrix3_t& H,
						const TriangularView<matrix3_t, Lower>& ltI):
		M_(ltM),
		H_(H),
		I_(ltI)
	{}

	// Accessor
	/// @return Mass matrix with a zero upper part.
	const matrix3_t& lowerTriangularMassMatrix() const
	{
		return M_;
	}

	/// @return Mass matrix.
	matrix3_t massMatrix() const
	{
		matrix3_t M;
		M.template triangularView<Upper>() = M_.transpose();
		M.template triangularView<StrictlyLower>() = M_;
		return M;
	}

	/// @return Generalized inertia matrix.
	const matrix3_t& gInertia() const
	{
		return H_;
	}

	/// @return Inertia matrix with a zero upper part.
	const matrix3_t& lowerTriangularInertia() const
	{
		return I_;
	}

	/// @return Inertia matrix.
	matrix3_t inertia() const
	{
		matrix3_t I;
		I.template triangularView<Upper>() = I_.transpose();
		I.template triangularView<StrictlyLower>() = I_;
		return I;
	}

	/// @retrun Non compact spatial articulated body inertia matrix.
	matrix6_t matrix() const
	{
		matrix6_t m;
		m << inertia(), H_,
				 H_.transpose(), massMatrix();
		return m;
	}

	template<typename T2>
	ABInertia<T2> cast() const
	{
		return ABInertia<T2>(M_.template cast<T2>(), H_.template cast<T2>(),
												I_.template cast<T2>());
	}

	// Operators
	ABInertia<T> operator+(const ABInertia<T>& rbI) const
	{
		matrix3_t M, I;
		M.template triangularView<Lower>() = M_ + rbI.M_;
		I.template triangularView<Lower>() = I_ + rbI.I_;
		return ABInertia<T>(M, H_ + rbI.H_, I);
	}

	ABInertia<T> operator-(const ABInertia<T>& rbI) const
	{
		matrix3_t M, I;
		M.template triangularView<Lower>() = M_ - rbI.M_;
		I.template triangularView<Lower>() = I_ - rbI.I_;
		return ABInertia<T>(M, H_ - rbI.H_, I);
	}

	ABInertia<T> operator-() const
	{
		return ABInertia<T>(-M_, -H_, -I_);
	}

	ABInertia<T>& operator+=(const ABInertia<T>& rbI)
	{
		M_.template triangularView<Lower>() += rbI.M_;
		H_ += rbI.H_;
		I_.template triangularView<Lower>() +=  rbI.I_;
		return *this;
	}

	ABInertia<T>& operator-=(const ABInertia<T>& rbI)
	{
		M_.template triangularView<Lower>() -= rbI.M_;
		H_ -= rbI.H_;
		I_.template triangularView<Lower>() -=  rbI.I_;
		return *this;
	}

	template<typename T2>
	ABInertia<T> operator*(T2 scalar) const
	{
		matrix3_t M, I;
		M.template triangularView<Lower>() = scalar*M_;
		I.template triangularView<Lower>() = scalar*I_;
		return ABInertia<T>(M, scalar*H_, I);
	}

	/// @return Ia + I
	ABInertia<T> operator+(const RBInertia<T>& rbI) const;

	/// @return Ia * v
	ForceVec<T> operator*(const MotionVec<T>& mv) const;

	/// @see operator*(const MotionVec<T>& mv) const
	template<typename Derived>
	void mul(const Eigen::MatrixBase<Derived>& mv,
		Eigen::MatrixBase<Derived>& result) const;

	bool operator==(const ABInertia<T>& abI) const
	{
		return M_ == abI.M_ && H_ == abI.H_ && I_ == abI.I_;
	}

	bool operator!=(const ABInertia<T>& abI) const
	{
		return M_ != abI.M_ || H_ != abI.H_ || I_ != abI.I_;
	}

private:
	matrix3_t M_;
	matrix3_t H_;
	matrix3_t I_;
};

template<typename T, typename T2>
inline ABInertia<T> operator*(T2 scalar, const ABInertia<T>& abI)
{
	return abI*scalar;
}

template<typename T>
inline std::ostream& operator<<(std::ostream& out, const ABInertia<T>& abI)
{
	out << abI.matrix();
	return out;
}

} // namespace sva
