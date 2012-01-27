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

// associated header
#include "Operators.h"

// includes
// Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

// SpaceVecAlg
#include "EigenUtility.h"

using namespace sva;

MotionVec cross(const MotionVec& mv1, const MotionVec& mv2)
{
	return MotionVec(mv1.angular().cross(mv2.angular()),
									 mv1.angular().cross(mv2.linear()) +
									 mv1.linear().cross(mv2.angular()));
}

ForceVec crossDual(const MotionVec& mv1, const ForceVec& fv2)
{
	return ForceVec(mv1.angular().cross(fv2.couple()) +
									mv1.linear().cross(fv2.force()),
									mv1.angular().cross(fv2.force()));
}

double dot(const MotionVec& mv1, const ForceVec& fv2)
{
	return mv1.angular().dot(fv2.couple()) + mv1.linear().dot(fv2.force());
}

ForceVec operator*(const RBInertia& rbI, const MotionVec& mv)
{
	return ForceVec(rbI.inertia()*mv.angular() + rbI.momentum().cross(mv.linear()),
									rbI.mass()*mv.linear() - rbI.momentum().cross(mv.angular()));
}

ABInertia operator+(const ABInertia& abI, const RBInertia& rbI)
{
	Matrix3d M, I;
	M.triangularView<Lower>() = abI.massMatrix() + Matrix3d::Identity()*rbI.mass();
	I.triangularView<Lower>() = abI.inertia() + rbI.inertia();
	return ABInertia(M, abI.gInertia() + vector3ToCrossMatrix(rbI.momentum()), I);
}

ForceVec operator*(const ABInertia& rbI, const MotionVec& mv)
{
	return ForceVec(rbI.inertia()*mv.angular() + rbI.gInertia()*mv.linear(),
									rbI.massMatrix()*mv.linear() +
									rbI.gInertia().transpose()*mv.angular());
}

MotionVec operator*(const PTransform& pt, const MotionVec& mv)
{
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	return MotionVec(E*mv.angular(),
									 E*(mv.linear() - r.cross(mv.angular())));
}

MotionVec invMul(const PTransform& pt, const MotionVec& mv)
{
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	return MotionVec(E.transpose()*mv.angular(),
									 E.transpose()*mv.linear() +
									 r.cross(E.transpose()*mv.angular()));
}

ForceVec dualMul(const PTransform& pt, const ForceVec& fv)
{
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	return ForceVec(E*(fv.couple() - r.cross(fv.force())),
									E*fv.force());
}

ForceVec transMul(const PTransform& pt, const ForceVec& fv)
{
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	Vector3d n = E.transpose()*fv.couple();
	n.noalias() += r.cross(E.transpose()*fv.force());
	return ForceVec(n,
									E.transpose()*fv.force());
}

RBInertia dualMul(const PTransform& pt, const RBInertia& rbI)
{
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	Matrix3d I;
	I.triangularView<Lower>() = E*(rbI.inertia() + vector3ToCrossMatrix(r)*
			vector3ToCrossMatrix(rbI.momentum()) +
			vector3ToCrossMatrix(rbI.momentum() - rbI.mass()*r)*
			vector3ToCrossMatrix(r))*E.transpose();
	return RBInertia(rbI.mass(),
									 E*(rbI.momentum() - rbI.mass()*r),
									 I);
}

RBInertia transMul(const PTransform& pt, const RBInertia& rbI)
{
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	Matrix3d I;
	I.triangularView<Lower>() = E.transpose()*rbI.inertia()*E -
			vector3ToCrossMatrix(r)*vector3ToCrossMatrix(E.transpose()*rbI.momentum()) -
			vector3ToCrossMatrix(E.transpose()*rbI.momentum() + rbI.mass()*r)*
			vector3ToCrossMatrix(r);
	return RBInertia(rbI.mass(), E.transpose()*rbI.momentum() + rbI.mass()*r, I);
}

ABInertia dualMul(const PTransform& pt, const ABInertia& rbI)
{
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();

	Matrix3d massM =rbI.massMatrix();
	Matrix3d rCross = vector3ToCrossMatrix(r);
	Matrix3d tmpI = rbI.gInertia() - rCross*massM;

	Matrix3d M, I;
	M.triangularView<Lower>() = E*massM*E.transpose();
	I.triangularView<Lower>() = E*(rbI.inertia() - rCross*rbI.gInertia().transpose() +
			(tmpI*rCross))*E.transpose();

	return ABInertia(M,
									 E*(tmpI)*E.transpose(),
									 I);
}

ABInertia transMul(const PTransform& pt, const ABInertia& rbI)
{
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();

	Matrix3d Mp = E.transpose()*rbI.massMatrix()*E;
	Matrix3d Hp = E.transpose()*rbI.gInertia()*E;
	Matrix3d rCross = vector3ToCrossMatrix(r);

	Matrix3d M, I;
	M.triangularView<Lower>() = Mp;
	I.triangularView<Lower>() = (E.transpose()*rbI.inertia()*E +
															 rCross*Hp.transpose() -
															 (Hp + rCross*Mp)*rCross);
	return ABInertia(M,
									 Hp + rCross*Mp,
									 I);

}

