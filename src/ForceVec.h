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
	* Spatial Force Vector compact representations.
	* See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
	*/
template<typename T>
class ForceVec
{
public:
	typedef Vector3<T> vector3_t;
	typedef Vector6<T> vector6_t;

public:
	ForceVec():
		fv_()
	{}

	/// @param vec Spatial force vector with couple in head and force in tail.
	ForceVec(const vector6_t& vec):
		fv_(vec)
	{}

	/**
		* @param couple Couple.
		* @param force Force.
		*/
	ForceVec(const vector3_t& couple, const vector3_t& force):
		fv_((vector6_t() << couple, force).finished())
	{}

	// Accessor
	/// @return Couple part (3 first parameters).
	vector3_t couple() const
	{
		return fv_.template head<3>();
	}

	/// @return Force part (3 last parameters).
	vector3_t force() const
	{
		return fv_.template tail<3>();
	}

	/// @return Non compact spatial force vector.
	const vector6_t& vector() const
	{
		return fv_;
	}

	/// @return Non compact spatial force vector.
	vector6_t& vector()
	{
		return fv_;
	}

	template<typename T2>
	ForceVec<T2> cast() const
	{
		return ForceVec<T2>(fv_.cast<T2>());
	}

	// Operators
	ForceVec<T> operator+(const ForceVec<T>& fv) const
	{
		return ForceVec<T>(fv_ + fv.fv_);
	}

	ForceVec<T> operator-(const ForceVec<T>& fv) const
	{
		return ForceVec<T>(fv_ - fv.fv_);
	}

	ForceVec<T> operator-() const
	{
		return ForceVec<T>(-fv_);
	}

	template<typename T2>
	ForceVec<T> operator*(T2 scalar) const
	{
		return ForceVec<T>(scalar * fv_);
	}

	bool operator==(const ForceVec<T>& fv) const
	{
		return fv_ == fv.fv_;
	}

	bool operator!=(const ForceVec<T>& fv) const
	{
		return fv_ != fv.fv_;
	}

private:
	vector6_t fv_;
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
