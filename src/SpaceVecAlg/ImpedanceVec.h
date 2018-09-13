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
// You should have received a copy of the GNU Lesser General Public License
// along with SpaceVecAlg.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <type_traits>

namespace sva
{

using namespace Eigen;

/**
	* Impedance Vector.
	* Mechanical impedance is a measure of how much a structure resists
	* motion when subjected to a given force:
	*
	*     F = Z * v
	*
	* where F is the exerted force, Z the impedance and v the resulting output
	* velocity. See <https://en.wikipedia.org/wiki/Mechanical_impedance>.
	*/
template<typename T>
class ImpedanceVec
{
public:
	typedef Vector3<T> vector3_t;
	typedef Vector6<T> vector6_t;

public:
	/// Zero impedance vector
	static ImpedanceVec<T> Zero()
	{
		return ImpedanceVec<T>(vector3_t::Zero(), vector3_t::Zero());
	}

public:
	/** Empty constructor */
	ImpedanceVec():
		angular_(),
		linear_()
	{}

	/** Define impedance from 6D vector.
		* @param vec Impedance vector with angular motion in head
		* and linear motion in tail.
		*/
	ImpedanceVec(const vector6_t& vec):
		angular_(vec.template head<3>()),
		linear_(vec.template tail<3>())
	{}

	/** Define impedance from angular and linear components.
		* @param angular Angular impedance.
		* @param linear Linear impedance.
		*/
	ImpedanceVec(const vector3_t& angular, const vector3_t& linear):
		angular_(angular),
		linear_(linear)
	{}

	/** Homogeneous impedance constructor.
		* @param angular Angular impedance.
		* @param linear Linear impedance.
		*/
	ImpedanceVec(T angular, T linear):
		angular_(angular, angular, angular),
		linear_(linear, linear, linear)
	{}

	// Accessor
	/// @return Angular impedance
	const vector3_t& angular() const
	{
		return angular_;
	}

	/// @return Angular impedance
	vector3_t& angular()
	{
		return angular_;
	}

	/// @return Linear impedance
	const vector3_t& linear() const
	{
		return linear_;
	}

	/// @return Linear impedance
	vector3_t& linear()
	{
		return linear_;
	}

	/// @return Non compact impedance vector.
	vector6_t vector() const
	{
		return ((vector6_t() << angular_, linear_).finished());
	}

	template<typename T2>
	ImpedanceVec<T2> cast() const
	{
		return ImpedanceVec<T2>(angular_.template cast<T2>(), linear_.template cast<T2>());
	}

	// Operators
	ImpedanceVec<T> operator+(const ImpedanceVec<T>& iv) const
	{
		return ImpedanceVec<T>(angular_ + iv.angular_, linear_ + iv.linear_);
	}

	ImpedanceVec<T>& operator+=(const ImpedanceVec<T>& iv)
	{
		angular_ += iv.angular_;
		linear_ += iv.linear_;
		return *this;
	}

	template<typename T2, typename std::enable_if<std::is_arithmetic<T2>::value, int>::type = 0>
	ImpedanceVec<T> operator*(T2 scalar) const
	{
		return ImpedanceVec<T>(scalar*angular_, scalar*linear_);
	}

	template<typename T2>
	ImpedanceVec<T> & operator*=(T2 scalar)
	{
		angular_ *= scalar;
		linear_ *= scalar;
		return *this;
	}

	template<typename T2>
	ImpedanceVec<T> operator/(T2 scalar) const
	{
		return ImpedanceVec<T>(angular_/scalar, linear_/scalar);
	}

	template<typename T2>
	ImpedanceVec<T> & operator/=(T2 scalar)
	{
		angular_ /= scalar;
		linear_ /= scalar;
		return *this;
	}

	bool operator==(const ImpedanceVec<T>& iv) const
	{
		return (angular_ == iv.angular_) && (linear_ == iv.linear_);
	}

	bool operator!=(const ImpedanceVec<T>& iv) const
	{
		return !(*this == iv);
	}

private:
	vector3_t angular_;
	vector3_t linear_;
};

template<typename T, typename T2>
inline ImpedanceVec<T> operator*(T2 scalar, const ImpedanceVec<T>& iv)
{
	return iv*scalar;
}

template<typename T>
inline ForceVec<T> operator*(const ImpedanceVec<T>& iv, const MotionVec<T>& mv)
{
	return ForceVec<T>(iv.angular().cwiseProduct(mv.angular()), iv.linear().cwiseProduct(mv.linear()));
}

template<typename T>
inline ForceVec<T> operator*(const MotionVec<T> & mv, const ImpedanceVec<T> & iv)
{
	return iv * mv;
}

template<typename T>
inline std::ostream& operator<<(std::ostream& out, const ImpedanceVec<T>& iv)
{
	out << iv.vector().transpose();
	return out;
}

} // namespace sva
