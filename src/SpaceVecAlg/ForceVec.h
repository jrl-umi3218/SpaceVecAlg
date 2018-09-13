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

#include <type_traits>

namespace sva
{

using namespace Eigen;

/**
	* Spatial Force Vector compact representations.
	* See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
	*/
template<typename T>
class ForceVec
{
public:
	typedef Vector3<T> vector3_t;
	typedef Vector6<T> vector6_t;

	friend class PTransform<T>;

public:
	/// Zero force vector
	static ForceVec<T> Zero()
	{
		return ForceVec<T>(vector3_t::Zero(), vector3_t::Zero());
	}

public:
	ForceVec():
		couple_(),
		force_()
	{}

	/// @param vec Spatial force vector with couple in head and force in tail.
	ForceVec(const vector6_t& vec):
		couple_(vec.template head<3>()),
		force_(vec.template tail<3>())
	{}

	/**
		* @param couple Couple.
		* @param force Force.
		*/
	ForceVec(const vector3_t& couple, const vector3_t& force):
		couple_(couple),
		force_(force)
	{}

	// Accessor
	/// @return Couple
	///
	/// @note The term "couple" is used in SpaceVecAlg with the general meaning
	/// of "moment" of a force vector. It should not to be confused with the
	/// more precise meaning of a "pure moment" that is also found in mechanics
	/// <https://en.wikipedia.org/wiki/Couple_(mechanics)>.
	///
	/// @sa moment()
	vector3_t& couple()
	{
		return couple_;
	}

	const vector3_t& couple() const
	{
		return couple_;
	}

	/// @return Moment
	vector3_t& moment()
	{
		return couple_;
	}

	const vector3_t& moment() const
	{
		return couple_;
	}

	/// @return Force
	vector3_t& force()
	{
		return force_;
	}

	const vector3_t& force() const
	{
		return force_;
	}

	/// @return Non compact spatial force vector.
	vector6_t vector() const
	{
		return ((vector6_t() << couple_, force_).finished());
	}

	template<typename T2>
	ForceVec<T2> cast() const
	{
		return ForceVec<T2>(couple_.template cast<T2>(), force_.template cast<T2>());
	}

	// Operators
	ForceVec<T> operator+(const ForceVec<T>& fv) const
	{
		return ForceVec<T>(couple_ + fv.couple_, force_ + fv.force_);
	}

	ForceVec<T> operator-(const ForceVec<T>& fv) const
	{
		return ForceVec<T>(couple_ - fv.couple_, force_ - fv.force_);
	}

	ForceVec<T> operator-() const
	{
		return ForceVec<T>(-couple_, -force_);
	}

	ForceVec<T>& operator+=(const ForceVec<T>& fv)
	{
		couple_ += fv.couple_;
		force_ += fv.force_;
		return *this;
	}

	ForceVec<T>& operator-=(const ForceVec<T>& fv)
	{
		couple_ -= fv.couple_;
		force_ -= fv.force_;
		return *this;
	}

	template<typename T2, typename std::enable_if<std::is_arithmetic<T2>::value, int>::type = 0>
	ForceVec<T> operator*(T2 scalar) const
	{
		return ForceVec<T>(scalar*couple_, scalar*force_);
	}

	template<typename T2>
	ForceVec<T> operator/(T2 scalar) const
	{
		return ForceVec<T>(couple_/scalar, force_/scalar);
	}

	bool operator==(const ForceVec<T>& fv) const
	{
		return (couple_ == fv.couple_) && (force_ == fv.force_);
	}

	bool operator!=(const ForceVec<T>& fv) const
	{
		return (couple_ != fv.couple_) || (force_ != fv.force_);
	}

private:
	vector3_t couple_;
	vector3_t force_;
};

template<typename T, typename T2>
inline ForceVec<T> operator*(T2 scalar, const ForceVec<T>& fv)
{
	return fv*scalar;
}

template<typename T>
inline std::ostream& operator<<(std::ostream& out, const ForceVec<T>& fv)
{
	out << fv.vector().transpose();
	return out;
}

} // namespace sva
