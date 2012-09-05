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
Matrix3d RotX(double theta);

/**
	* Create a rotation matrix about the Y axis.
	* The rotation is exprimed in successor frame.
	* @param theta rotation in radian.
	*/
Matrix3d RotY(double theta);

/**
	* Create a rotation matrix about the Z axis.
	* The rotation is exprimed in successor frame.
	* @param theta rotation in radian.
	*/
Matrix3d RotZ(double theta);

/**
	* Compute the rotation error between two matrix.
	* Compute the XYZ rotation of rotTo rotation matrix
	* calculated has follow rotDes = rotTo*rotCur
	* @param prec Precision to know if the rotTo matrix is identity.
	* @return XYZ rotation in radian.
	*/
Vector3d rotationError(const Matrix3d& rotCur, const Matrix3d& rotDes,
	double prec=1e-8);

/**
	* Compute the rotation vector of the rotation matrix.
	* If we integrate this rotation vector for 1 second we must
	* have the rotation matrix rot.
	* (see exponential matrix and logarithmic matrix).
	* @param prec Precision to know if the rot matrix is identity.
	*/
Vector3d rotationVelocity(const Matrix3d& rot, double prec=1e-8);


/**
	* Plücker transform compact representation.
	* Use 3D matrix as rotation internal representation.
	* Quaternion are inversed as they must be expressed in successor frame.
	* See Roy Featherstone «Rigid Body Dynamics Algorithms» page 247.
	*/
class PTransform
{

public:
	/// Identity transformation.
	static PTransform Identity()
	{
		return PTransform(Matrix3d::Identity(), Vector3d::Zero());
	}

public:
	// Constructors
	/// Default constructor. Rotation and translation are uninitialized.
	PTransform():
		E_(),
		r_()
	{}

	/**
		* @param rot Rotation matrix.
		* @param trans Translation vector.
		*/
	PTransform(const Matrix3d& rot, const Vector3d& trans):
		E_(rot),
		r_(trans)
	{}

	/**
		* @param rot Rotation quaternion.
		* @param trans Translation vector.
		*/
	PTransform(const Quaterniond& rot, const Vector3d& trans):
		E_(),
		r_(trans)
	{
		E_ = rot.inverse().toRotationMatrix();
	}

	/**
		* Rotation only transform.
		* @param rot Rotation quaternion.
		*/
	PTransform(const Quaterniond& rot):
		E_(),
		r_(Vector3d::Zero())
	{
		E_ = rot.inverse().toRotationMatrix();
	}

	/**
		* Rotation only transform.
		* @param rot Rotation matrix.
		*/
	PTransform(const Matrix3d& rot):
		E_(rot),
		r_(Vector3d::Zero())
	{}

	/**
		* Translation only transform.
		* @param trans Translation vector.
		*/
	PTransform(const Vector3d& trans):
		E_(Matrix3d::Identity()),
		r_(trans)
	{}

	// Accessor
	/// @return Rotation matrix.
	const Matrix3d& rotation() const
	{
		return E_;
	}

	/// @return Rotation matrix.
	Matrix3d& rotation()
	{
		return E_;
	}

	/// @return Translation vector.
	const Vector3d& translation() const
	{
		return r_;
	}

	/// @return Translation vector.
	Vector3d& translation()
	{
		return r_;
	}

	/// @return Non compact Plücker transformation matrix.
	Matrix6d matrix() const
	{
		Matrix6d m;
		m << E_, Matrix3d::Zero(),
				 -E_*vector3ToCrossMatrix(r_), E_;
		return m;
	}

	/// @return Non compact dual Plücker transformation matrix.
	Matrix6d dualMatrix() const
	{
		Matrix6d m;
		m << E_, -E_*vector3ToCrossMatrix(r_),
				 Matrix3d::Zero(), E_;
		return m;
	}

	// Operators
	/// @return X*X
	PTransform operator*(const PTransform& pt) const
	{
		return PTransform(E_*pt.E_, pt.r_ + pt.E_.transpose()*r_);
	}

	/// @return Xv
	MotionVec operator*(const MotionVec& mv) const;
	/// @return X⁻¹v
	MotionVec invMul(const MotionVec& mv) const;

	/// @return X*v
	ForceVec dualMul(const ForceVec& fv) const;
	/// @return Xtv
	ForceVec transMul(const ForceVec& fv) const;

	/// @return X*IX⁻¹
	RBInertia dualMul(const RBInertia& rbI) const;
	/// @return XtIX
	RBInertia transMul(const RBInertia& rbI) const;

	/// @return X*IX⁻¹
	ABInertia dualMul(const ABInertia& rbI) const;
	/// @return XtIX
	ABInertia transMul(const ABInertia& rbI) const;

	/// @return Inverse Plücker transformation.
	PTransform inv() const
	{
		return PTransform(E_.transpose(), -E_*r_);
	}

	bool operator==(const PTransform& pt) const
	{
		return E_ == pt.E_ && r_ == pt.r_;
	}

	bool operator!=(const PTransform& pt) const
	{
		return E_ != pt.E_ || r_ != pt.r_;
	}

private:
	Matrix3d E_;
	Vector3d r_;
};


inline Matrix3d RotX(double theta)
{
	double s = std::sin(theta), c = std::cos(theta);
	return (Matrix3d() << 1., 0., 0.,
												0., c, s,
												0., -s, c).finished();
}


inline Matrix3d RotY(double theta)
{
	double s = std::sin(theta), c = std::cos(theta);
	return (Matrix3d() << c, 0., -s,
												0., 1., 0.,
												s, 0., c).finished();
}


inline Matrix3d RotZ(double theta)
{
	double s = std::sin(theta), c = std::cos(theta);
	return (Matrix3d() << c, s, 0.,
												-s, c, 0.,
												0., 0., 1.).finished();
}


inline Vector3d rotationError(const Matrix3d& rotCur, const Matrix3d& rotDes,
	double prec)
{
	Matrix3d rotTo = rotDes*rotCur.transpose();
	return rotCur*rotationVelocity(rotTo, prec);
}


inline Vector3d rotationVelocity(const Matrix3d& rot, double prec)
{
	Vector3d w;
	double theta = std::acos((rot(0,0) + rot(1,1) + rot(2,2) - 1.)*0.5);

	if(rot.isIdentity(prec))
	{
		w.setZero();
	}
	else
	{
		w = Vector3d(-rot(2,1) + rot(1,2),
								 -rot(0,2) + rot(2,0),
								 -rot(1,0) + rot(0,1));
		w *= theta/(2.*std::sin(theta));
	}

	return w;
}


inline std::ostream& operator<<(std::ostream& out, const PTransform& pt)
{
	out << pt.matrix();
	return out;
}

} // namespace sva
