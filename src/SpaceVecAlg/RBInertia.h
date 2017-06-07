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

template<typename T>
Matrix3<T> inertiaToOrigin(const Matrix3<T>& inertia,
												 T mass, const Vector3<T>& com,
												 const Matrix3<T>& rotation);

/**
	* Spatial Rigid Body Inertia compact representation.
	* See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
	*/
template<typename T>
class RBInertia
{
public:
	typedef Vector3<T> vector3_t;
	typedef Matrix3<T> matrix3_t;
	typedef Matrix6<T> matrix6_t;

public:
	RBInertia():
		m_(),
		h_(),
		I_()
	{}

	/**
		* @param m Mass.
		* @param h Spatial momentum.
		* @param I Inertia matrix at body origin.
		*/
	RBInertia(T m, const vector3_t& h, const matrix3_t& I):
		m_(m),
		h_(h),
		I_(matrix3_t::Zero())
	{
		I_.template triangularView<Lower>() = I;
	}

	/**
		* @param m Mass.
		* @param h Spatial momentum.
		* @param I Lower triangular view of Inertia matrix at body origin.
		*/
	RBInertia(T m, const vector3_t& h,
						const TriangularView<matrix3_t, Lower>& ltI):
		m_(m),
		h_(h),
		I_(ltI)
	{}

	// Accessor
	/// @return Mass.
	T mass() const
	{
		return m_;
	}

	/// @return Spatial momentum.
	const vector3_t& momentum() const
	{
		return h_;
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

	/// @retrun Non compact spatial rigid body inertia matrix.
	matrix6_t matrix() const
	{
		matrix6_t m;
		matrix3_t hCross = vector3ToCrossMatrix(h_);
		m << inertia(), hCross,
				 hCross.transpose(), matrix3_t::Identity()*m_;
		return m;
	}

	template<typename T2>
	RBInertia<T2> cast() const
	{
		return RBInertia<T2>(T2(m_), h_.template cast<T2>(), I_.template cast<T2>());
	}

	// Operators
	RBInertia<T> operator+(const RBInertia<T>& rbI) const
	{
		matrix3_t I;
		I.template triangularView<Lower>() = I_ + rbI.I_;
		return RBInertia<T>(m_ + rbI.m_,
											 h_ + rbI.h_,
											 I);
	}

	RBInertia<T> operator-(const RBInertia<T>& rbI) const
	{
		matrix3_t I;
		I.template triangularView<Lower>() = I_ - rbI.I_;
		return RBInertia<T>(m_ - rbI.m_,
											 h_ - rbI.h_,
											 I);
	}

	RBInertia<T> operator-() const
	{
		return RBInertia<T>(-m_, -h_, -I_);
	}

	RBInertia<T>& operator+=(const RBInertia<T>& rbI)
	{
		I_.template triangularView<Lower>() += rbI.I_;
		m_ += rbI.m_;
		h_ += rbI.h_;
		return *this;
	}

	RBInertia<T>& operator-=(const RBInertia<T>& rbI)
	{
		I_.template triangularView<Lower>() -= rbI.I_;
		m_ -= rbI.m_;
		h_ -= rbI.h_;
		return *this;
	}

	template<typename T2>
	RBInertia<T> operator*(T2 scalar) const
	{
		matrix3_t I;
		I.template triangularView<Lower>() = scalar *I_;
		return RBInertia<T>(scalar * m_, scalar * h_, I);
	}

	/// @return I*v
	ForceVec<T> operator*(const MotionVec<T>& mv) const;

	/// @see operator*(const MotionVec<T>& mv) const
	template<typename Derived>
	void mul(const Eigen::MatrixBase<Derived>& mv,
		Eigen::MatrixBase<Derived>& result) const;

	bool operator==(const RBInertia<T>& rbI) const
	{
		return m_ == rbI.m_ && h_ == rbI.h_ && I_ == rbI.I_;
	}

	bool operator!=(const RBInertia<T>& rbI) const
	{
		return m_ != rbI.m_ || h_ != rbI.h_ || I_ != rbI.I_;
	}

private:
	T m_;
	vector3_t h_;
	matrix3_t I_;
};

template <typename T, typename T2>
inline RBInertia<T> operator*(T2 scalar, const RBInertia<T>& rbI)
{
	return rbI*scalar;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& out, const RBInertia<T>& rbI)
{
	out << rbI.matrix();
	return out;
}

template<typename T>
Matrix3<T> inertiaToOrigin(const Matrix3<T>& inertia,
												 T mass, const Vector3<T>& com,
												 const Matrix3<T>& rotation)
{
	Matrix3<T> trans = vector3ToCrossMatrix<T>(mass*com)*vector3ToCrossMatrix<T>(com).transpose();
	return rotation*(inertia + trans)*rotation.transpose();
}

} // namespace sva
