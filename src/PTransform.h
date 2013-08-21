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
	* Compute the rotation error between two matrix in world frame.
	* Compute the XYZ rotation of rotTo rotation matrix
	* calculated has follow rotDes = rotTo*rotCur
	* @param prec Precision to know if the rotTo matrix is identity.
	* @return XYZ rotation in radian.
	*/
template<typename T>
Vector3<T> rotationError(const Matrix3<T>& rotCur, const Matrix3<T>& rotDes,
	double prec=1e-8);

/**
	* Compute the rotation vector of the rotation matrix.
	* If we integrate this rotation vector for 1 second we must
	* have the rotation matrix rot.
	* (see exponential matrix and logarithmic matrix).
	* @param prec Precision to know if the rot matrix is identity.
	*/
template<typename T>
Vector3<T> rotationVelocity(const Matrix3<T>& rot, double prec=1e-8);


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
		E_(),
		r_(trans)
	{
		E_ = rot.inverse().toRotationMatrix();
	}

	/**
		* Rotation only transform.
		* @param rot Rotation quaternion.
		*/
	PTransform(const quaternion_t& rot):
		E_(),
		r_(vector3_t::Zero())
	{
		E_ = rot.inverse().toRotationMatrix();
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
		return PTransform<T2>(E_.cast<T2>(), r_.cast<T2>());
	}

	// Operators
	/// @return X*X
	PTransform<T> operator*(const PTransform<T>& pt) const
	{
		return PTransform<T>(E_*pt.E_, pt.r_ + pt.E_.transpose()*r_);
	}

	/// @return Xv
	MotionVec<T> operator*(const MotionVec<T>& mv) const;
	/// @return X⁻¹v
	MotionVec<T> invMul(const MotionVec<T>& mv) const;

	/// @return X*v
	ForceVec<T> dualMul(const ForceVec<T>& fv) const;
	/// @return Xtv
	ForceVec<T> transMul(const ForceVec<T>& fv) const;

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
inline Vector3<T> rotationError(const Matrix3<T>& rotCur, const Matrix3<T>& rotDes,
	double prec)
{
	Matrix3<T> rotTo = rotDes*rotCur.transpose();
	return Vector3<T>(rotCur.transpose()*rotationVelocity(rotTo, prec));
}


template<typename T>
inline Vector3<T> rotationVelocity(const Matrix3<T>& rot, double prec)
{
	Vector3<T> w;
	T theta = std::acos((rot(0,0) + rot(1,1) + rot(2,2) - 1.)*0.5);

	if(rot.isIdentity(prec))
	{
		w.setZero();
	}
	else
	{
		w = Vector3<T>(-rot(2,1) + rot(1,2),
									-rot(0,2) + rot(2,0),
									-rot(1,0) + rot(0,1));
		w *= theta/(2.*std::sin(theta));
	}

	return w;
}


template<typename T>
inline std::ostream& operator<<(std::ostream& out, const PTransform<T>& pt)
{
	out << pt.matrix();
	return out;
}

} // namespace sva
