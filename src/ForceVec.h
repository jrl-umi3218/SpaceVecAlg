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
class ForceVec
{
public:
	ForceVec():
		fv_()
	{}

	/// @param vec Spatial force vector with couple in head and force in tail.
	ForceVec(const Vector6d& vec):
		fv_(vec)
	{}

	/**
		* @param couple Couple.
		* @param force Force.
		*/
	ForceVec(const Vector3d& couple, const Vector3d& force):
		fv_((Vector6d() << couple, force).finished())
	{}

	// Accessor
	/// @return Couple part (3 first parameters).
	Vector3d couple() const
	{
		return fv_.head<3>();
	}

	/// @return Force part (3 last parameters).
	Vector3d force() const
	{
		return fv_.tail<3>();
	}

	/// @return Non compact spatial force vector.
	const Vector6d& vector() const
	{
		return fv_;
	}

	/// @return Non compact spatial force vector.
	Vector6d& vector()
	{
		return fv_;
	}

	// Operators
	ForceVec operator+(const ForceVec& fv) const
	{
		return ForceVec(fv_ + fv.fv_);
	}

	ForceVec operator-(const ForceVec& fv) const
	{
		return ForceVec(fv_ - fv.fv_);
	}

	ForceVec operator-() const
	{
		return ForceVec(-fv_);
	}

	ForceVec operator*(double scalar) const
	{
		return ForceVec(scalar * fv_);
	}

private:
	Vector6d fv_;
};

inline ForceVec operator*(double scalar, const ForceVec& fv)
{
	return fv*scalar;
}

} // namespace sva
