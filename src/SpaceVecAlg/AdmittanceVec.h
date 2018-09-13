// Copyright 2012-2018 CNRS-UM LIRMM, CNRS-AIST JRL
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
// You should have receaved a copy of the GNU Lesser General Public License
// along with SpaceVecAlg.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <type_traits>

namespace sva
{

using namespace Eigen;

/**
	* Admittance Vector.
	* Mechanical admittance is a measure of how much a structure moves
	* in response to a given force:
	*
	*     v = A * F
	*
	* where F is the exerted force, A the admittance and v the resulting output
	* velocity. Admittance is the reciprocal of mechanical impedance; see
	* <https://en.wikipedia.org/wiki/Mechanical_impedance>.
	*/
template<typename T>
class AdmittanceVec
{
public:
	typedef Vector3<T> vector3_t;
	typedef Vector6<T> vector6_t;

public:
	/// Zero admittance vector
	static AdmittanceVec<T> Zero()
	{
		return AdmittanceVec<T>(vector3_t::Zero(), vector3_t::Zero());
	}

public:
	/** Empty constructor */
	AdmittanceVec():
		angular_(),
		linear_()
	{}

	/** Define admittance from 6D vector.
		* @param vec Admittance vector with angular motion in head
		* and linear motion in tail.
		*/
	AdmittanceVec(const vector6_t& vec):
		angular_(vec.template head<3>()),
		linear_(vec.template tail<3>())
	{}

	/** Define admittance from angular and linear components.
		* @param angular Angular admittance.
		* @param linear Linear admittance.
		*/
	AdmittanceVec(const vector3_t& angular, const vector3_t& linear):
		angular_(angular),
		linear_(linear)
	{}

	/** Homogeneous admittance constructor.
		* @param angular Angular admittance.
		* @param linear Linear admittance.
		*/
	AdmittanceVec(T angular, T linear):
		angular_(angular, angular, angular),
		linear_(linear, linear, linear)
	{}

	// Accessor
	/// @return Angular admittance
	const vector3_t& angular() const
	{
		return angular_;
	}

	/// @return Angular admittance
	vector3_t& angular()
	{
		return angular_;
	}

	/// @return Linear admittance
	const vector3_t& linear() const
	{
		return linear_;
	}

	/// @return Linear admittance
	vector3_t& linear()
	{
		return linear_;
	}

	/// @return Non-compact admittance vector.
	vector6_t vector() const
	{
		return ((vector6_t() << angular_, linear_).finished());
	}

	template<typename T2>
	AdmittanceVec<T2> cast() const
	{
		return AdmittanceVec<T2>(angular_.template cast<T2>(), linear_.template cast<T2>());
	}

	// Operators
	AdmittanceVec<T> operator+(const AdmittanceVec<T>& av) const
	{
		return AdmittanceVec<T>(angular_ + av.angular_, linear_ + av.linear_);
	}

	AdmittanceVec<T>& operator+=(const AdmittanceVec<T>& av)
	{
		angular_ += av.angular_;
		linear_ += av.linear_;
		return *this;
	}

	template<typename T2, typename std::enable_if<std::is_arithmetic<T2>::value, int>::type = 0>
	AdmittanceVec<T> operator*(T2 scalar) const
	{
		return AdmittanceVec<T>(scalar*angular_, scalar*linear_);
	}

	template<typename T2>
	AdmittanceVec<T> & operator*=(T2 scalar)
	{
		angular_ *= scalar;
		linear_ *= scalar;
		return *this;
	}

	template<typename T2>
	AdmittanceVec<T> operator/(T2 scalar) const
	{
		return AdmittanceVec<T>(angular_/scalar, linear_/scalar);
	}

	template<typename T2>
	AdmittanceVec<T> & operator/=(T2 scalar)
	{
		angular_ /= scalar;
		linear_ /= scalar;
		return *this;
	}

	bool operator==(const AdmittanceVec<T>& av) const
	{
		return (angular_ == av.angular_) && (linear_ == av.linear_);
	}

	bool operator!=(const AdmittanceVec<T>& av) const
	{
		return !(*this == av);
	}

private:
	vector3_t angular_;
	vector3_t linear_;
};

template<typename T, typename T2>
inline AdmittanceVec<T> operator*(T2 scalar, const AdmittanceVec<T>& av)
{
	return av*scalar;
}

template<typename T>
inline MotionVec<T> operator*(const AdmittanceVec<T>& av, const ForceVec<T>& fv)
{
	return MotionVec<T>(av.angular().cwiseProduct(fv.couple()), av.linear().cwiseProduct(fv.force()));
}

template<typename T>
inline MotionVec<T> operator*(const ForceVec<T>& fv, const AdmittanceVec<T>& av)
{
	return av * fv;
}

template<typename T>
inline std::ostream& operator<<(std::ostream& out, const AdmittanceVec<T>& av)
{
	out << av.vector().transpose();
	return out;
}

} // namespace sva
