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
	* Create a rotation matrix about the X axis.
	* The rotation is exprimed in successor frame.
	* @param theta rotation in radian.
	*/
template<typename T>
Matrix3<T> RotX(T theta);

/**
	* Create a rotation matrix about the Y axis.
	* The rotation is exprimed in successor frame.
	* @param theta rotation in radian.
	*/
template<typename T>
Matrix3<T> RotY(T theta);

/**
	* Create a rotation matrix about the Z axis.
	* The rotation is exprimed in successor frame.
	* @param theta rotation in radian.
	*/
template<typename T>
Matrix3<T> RotZ(T theta);

/**
	* Compute the 3D rotation error between two matrix E_a_b and E_a_c in the 'a' frame.
	* This method convert the 3D rotation matrix E_b_c into a rotation vector.
	* The matrix E_b_c is computed as follow E_a_c = E_b_c*E_a_b.
	* Then the error is computed with E_b_a*rotationVelocity(E_b_c).
	* @return XYZ rotation in radian.
	*/
template<typename T>
Vector3<T> rotationError(const Matrix3<T>& E_a_b, const Matrix3<T>& E_a_c);

/**
	* Compute the 3D rotation vector of the rotation matrix E_a_b in the 'a' frame.
	* If we integrate this rotation vector for 1 second we must
	* have the rotation matrix E_a_b.
	* (see exponential matrix and logarithmic matrix).
	*/
template<typename T>
Vector3<T> rotationVelocity(const Matrix3<T>& E_a_b);

/**
	* Compute the 6D error between two PTransform in the 'a' frame.
	* This method convert the 6D transformation matrix X_b_c into a motion vector.
	* The matrix X_b_c is computed as follow X_a_c = X_b_c*X_a_b.
	* Then the error is computed with PTransform(E_b_a)*transformVelocity(X_b_c).
	* @return XYZ rotation in radian.
	*/
template<typename T>
MotionVec<T> transformError(const PTransform<T>& X_a_b, const PTransform<T>& X_a_c);

/**
	* Compute the motion vector of the matrix X_a_b in the 'a' frame.
	* If we integrate this motion vector for 1 second we must
	* have the transformation matrix X_a_b.
	* This function can be see as an implementation of the function XtoV
	* (see Featherstone appendix) but with the use of logarithmic
	* matrix to compute the rotational error.
	*/
template<typename T>
MotionVec<T> transformVelocity(const PTransform<T>& X_a_b);


/**
	* Plücker transform compact representation.
	* Use 3D matrix as rotation internal representation.
	* Quaternion are inversed as they must be expressed in successor frame.
	* See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
	*/
template<typename T>
class PTransform
{
	typedef Vector3<T> vector3_t;
	typedef Matrix3<T> matrix3_t;
	typedef Matrix6<T> matrix6_t;
	typedef Quaternion<T> quaternion_t;

public:
	/// Identity transformation.
	static PTransform<T> Identity()
	{
		return PTransform<T>(matrix3_t::Identity(), vector3_t::Zero());
	}

public:
	// Constructors
	/// Default constructor. Rotation and translation are uninitialized.
	PTransform():
		E_(),
		r_()
	{}

	/// Copy constructor.
	template<typename T2>
	PTransform(const PTransform<T2>& pt):
		E_(pt.rotation().template cast<T>()),
		r_(pt.translation().template cast<T>())
	{}

	/**
		* @param rot Rotation matrix.
		* @param trans Translation vector.
		*/
	PTransform(const matrix3_t& rot, const vector3_t& trans):
		E_(rot),
		r_(trans)
	{}

	/**
		* @param rot Rotation quaternion.
		* @param trans Translation vector.
		*/
	PTransform(const quaternion_t& rot, const vector3_t& trans):
		E_(rot.matrix()),
		r_(trans)
	{
	}

	/**
		* Rotation only transform.
		* @param rot Rotation quaternion.
		*/
	PTransform(const quaternion_t& rot):
		E_(rot.matrix()),
		r_(vector3_t::Zero())
	{
	}

	/**
		* Rotation only transform.
		* @param rot Rotation matrix.
		*/
	PTransform(const matrix3_t& rot):
		E_(rot),
		r_(vector3_t::Zero())
	{}

	/**
		* Translation only transform.
		* @param trans Translation vector.
		*/
	PTransform(const vector3_t& trans):
		E_(matrix3_t::Identity()),
		r_(trans)
	{}

	// Accessor
	/// @return Rotation matrix.
	const matrix3_t& rotation() const
	{
		return E_;
	}

	/// @return Rotation matrix.
	matrix3_t& rotation()
	{
		return E_;
	}

	/// @return Translation vector.
	const vector3_t& translation() const
	{
		return r_;
	}

	/// @return Translation vector.
	vector3_t& translation()
	{
		return r_;
	}

	/// @return Non compact Plücker transformation matrix.
	matrix6_t matrix() const
	{
		matrix6_t m;
		m << E_, matrix3_t::Zero(),
				 -E_*vector3ToCrossMatrix(r_), E_;
		return m;
	}

	/// @return Non compact dual Plücker transformation matrix.
	matrix6_t dualMatrix() const
	{
		matrix6_t m;
		m << E_, -E_*vector3ToCrossMatrix(r_),
				 matrix3_t::Zero(), E_;
		return m;
	}

	template<typename T2>
	PTransform<T2> cast() const
	{
		return PTransform<T2>(E_.template cast<T2>(), r_.template cast<T2>());
	}

	// Operators
	/// @return X*X
	PTransform<T> operator*(const PTransform<T>& pt) const
	{
		return PTransform<T>(E_*pt.E_, pt.r_ + pt.E_.transpose()*r_);
	}

	/// @return Xv
	MotionVec<T> operator*(const MotionVec<T>& mv) const;

	/// @return compute angular part of @see operator*(const MotionVec<T>& mv) const
	Eigen::Vector3<T> angularMul(const MotionVec<T>& mv) const;
	/// @return compute linear part of @see operator*(const MotionVec<T>& mv) const
	Eigen::Vector3<T> linearMul(const MotionVec<T>& mv) const;
	/// @see operator*(const MotionVec<T>& mv) const;
	template<typename Derived>
	void mul(const Eigen::MatrixBase<Derived>& mv,
		Eigen::MatrixBase<Derived>& result) const;


	/// @return X⁻¹v
	MotionVec<T> invMul(const MotionVec<T>& mv) const;
	/// @return compute angular part of @see invMul
	Eigen::Vector3<T> angularInvMul(const MotionVec<T>& mv) const;
	/// @return compute linear part of @see invMul
	Eigen::Vector3<T> linearInvMul(const MotionVec<T>& mv) const;
	/// @see invMul
	template<typename Derived>
	void invMul(const Eigen::MatrixBase<Derived>& mv,
		Eigen::MatrixBase<Derived>& result) const;


	/// @return X*v
	ForceVec<T> dualMul(const ForceVec<T>& fv) const;
	/// @return compute couple part of @see dualMul
	Eigen::Vector3<T> coupleDualMul(const ForceVec<T>& fv) const;
	/// @return compute force part of @see dualMul
	Eigen::Vector3<T> forceDualMul(const ForceVec<T>& fv) const;
	/// @see dualMul
	template<typename Derived>
	void dualMul(const Eigen::MatrixBase<Derived>& fv,
		Eigen::MatrixBase<Derived>& result) const;


	/// @return Xtv
	ForceVec<T> transMul(const ForceVec<T>& fv) const;
	/// @return compute couple part of @see transMul
	Eigen::Vector3<T> coupleTransMul(const ForceVec<T>& fv) const;
	/// @return compute force part of @see transMul
	Eigen::Vector3<T> forceTransMul(const ForceVec<T>& fv) const;
	/// @see transMul
	template<typename Derived>
	void transMul(const Eigen::MatrixBase<Derived>& fv,
		Eigen::MatrixBase<Derived>& result) const;


	/// @return X*IX⁻¹
	RBInertia<T> dualMul(const RBInertia<T>& rbI) const;
	/// @return XtIX
	RBInertia<T> transMul(const RBInertia<T>& rbI) const;

	/// @return X*IX⁻¹
	ABInertia<T> dualMul(const ABInertia<T>& rbI) const;
	/// @return XtIX
	ABInertia<T> transMul(const ABInertia<T>& rbI) const;

	/// @return Inverse Plücker transformation.
	PTransform<T> inv() const
	{
		return PTransform<T>(E_.transpose(), -E_*r_);
	}

	bool operator==(const PTransform<T>& pt) const
	{
		return E_ == pt.E_ && r_ == pt.r_;
	}

	bool operator!=(const PTransform<T>& pt) const
	{
		return E_ != pt.E_ || r_ != pt.r_;
	}

private:
	matrix3_t E_;
	vector3_t r_;
};


template<typename T>
inline Matrix3<T> RotX(T theta)
{
	T s = std::sin(theta), c = std::cos(theta);
	return (Matrix3<T>() << 1., 0., 0.,
												0., c, s,
												0., -s, c).finished();
}


template<typename T>
inline Matrix3<T> RotY(T theta)
{
	T s = std::sin(theta), c = std::cos(theta);
	return (Matrix3<T>() << c, 0., -s,
												0., 1., 0.,
												s, 0., c).finished();
}


template<typename T>
inline Matrix3<T> RotZ(T theta)
{
	T s = std::sin(theta), c = std::cos(theta);
	return (Matrix3<T>() << c, s, 0.,
												-s, c, 0.,
												0., 0., 1.).finished();
}


template<typename T>
inline Vector3<T> rotationError(const Matrix3<T>& E_a_b, const Matrix3<T>& E_a_c)
{
	Matrix3<T> E_b_c = E_a_c*E_a_b.transpose();
	return Vector3<T>(E_a_b.transpose()*rotationVelocity(E_b_c));
}


template<typename T>
inline Vector3<T> rotationVelocity(const Matrix3<T>& E_a_b)
{
	Vector3<T> w;
	T acosV = (E_a_b(0,0) + E_a_b(1,1) + E_a_b(2,2) - 1.)*0.5;
	T theta = std::acos(std::min(std::max(acosV,-1.),1.));

	w = Vector3<T>(-E_a_b(2,1) + E_a_b(1,2),
								-E_a_b(0,2) + E_a_b(2,0),
								-E_a_b(1,0) + E_a_b(0,1));
	w *= sinc_inv(theta)*0.5;

	return w;
}


template<typename T>
inline MotionVec<T> transformError(const PTransform<T>& X_a_b, const PTransform<T>& X_a_c)
{
	PTransform<T> X_b_c = X_a_c*X_a_b.inv();
	return PTransform<T>(Matrix3<T>(X_a_b.rotation().transpose()))*transformVelocity(X_b_c);
}


template<typename T>
inline MotionVec<T> transformVelocity(const PTransform<T>& X_a_b)
{
	return MotionVec<T>(rotationVelocity(X_a_b.rotation()),
		X_a_b.translation());
}


// interpolate between transformations, t must be between 0 and 1
template<typename T>
PTransform<T> interpolate(const PTransform<T>& from, const PTransform<T>& to, double t)
{
	Eigen::Quaternion<T> qfrom(from.rotation());
	Eigen::Quaternion<T> qto(to.rotation());
	PTransform<T> result(qfrom.slerp(t, qto), (from.translation()*(1.-t) + to.translation()*t));
	return result;
}


template<typename T>
inline std::ostream& operator<<(std::ostream& out, const PTransform<T>& pt)
{
	out << pt.matrix();
	return out;
}

} // namespace sva
