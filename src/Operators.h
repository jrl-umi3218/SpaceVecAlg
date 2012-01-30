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

// Don't include it directly, include SpaceVecAlg instead

namespace sva
{

// sva::MotionVec Left
sva::MotionVec cross(const sva::MotionVec& mv1, const sva::MotionVec& mv2);
sva::ForceVec crossDual(const sva::MotionVec& mv1, const sva::ForceVec& fv2);
double dot(const sva::MotionVec& mv1, const sva::ForceVec& fv2);

// sva::RBInertia Left
sva::ForceVec operator*(const sva::RBInertia& rbI, const sva::MotionVec& mv);

// sva::ABInertia Left
sva::ABInertia operator+(const sva::ABInertia& abI, const sva::RBInertia& rbI);
sva::ForceVec operator*(const sva::ABInertia& rbI, const sva::MotionVec& mv);

// sva::PTransform Left
sva::MotionVec operator*(const sva::PTransform& pt, const sva::MotionVec& mv);
sva::MotionVec invMul(const sva::PTransform& pt, const sva::MotionVec& mv);

sva::ForceVec dualMul(const sva::PTransform& pt, const sva::ForceVec& fv);
sva::ForceVec transMul(const sva::PTransform& pt, const sva::ForceVec& fv);

sva::RBInertia dualMul(const sva::PTransform& pt, const sva::RBInertia& rbI);
sva::RBInertia transMul(const sva::PTransform& pt, const sva::RBInertia& rbI);

sva::ABInertia dualMul(const sva::PTransform& pt, const sva::ABInertia& rbI);
sva::ABInertia transMul(const sva::PTransform& pt, const sva::ABInertia& rbI);

// Implementation

inline MotionVec cross(const MotionVec& mv1, const MotionVec& mv2)
{
	return MotionVec(mv1.angular().cross(mv2.angular()),
									 mv1.angular().cross(mv2.linear()) +
									 mv1.linear().cross(mv2.angular()));
}

inline ForceVec crossDual(const MotionVec& mv1, const ForceVec& fv2)
{
	return ForceVec(mv1.angular().cross(fv2.couple()) +
									mv1.linear().cross(fv2.force()),
									mv1.angular().cross(fv2.force()));
}

inline double dot(const sva::MotionVec& mv1, const sva::ForceVec& fv2)
{
	return mv1.angular().dot(fv2.couple()) + mv1.linear().dot(fv2.force());
}

inline ForceVec operator*(const RBInertia& rbI, const MotionVec& mv)
{
	return ForceVec(rbI.inertia()*mv.angular() + rbI.momentum().cross(mv.linear()),
									rbI.mass()*mv.linear() - rbI.momentum().cross(mv.angular()));
}

inline ABInertia operator+(const ABInertia& abI, const RBInertia& rbI)
{
	using namespace Eigen;
	Matrix3d M, I;
	M.triangularView<Lower>() = abI.massMatrix() + Matrix3d::Identity()*rbI.mass();
	I.triangularView<Lower>() = abI.inertia() + rbI.inertia();
	return ABInertia(M, abI.gInertia() + vector3ToCrossMatrix(rbI.momentum()), I);
}

inline ForceVec operator*(const ABInertia& rbI, const MotionVec& mv)
{
	return ForceVec(rbI.inertia()*mv.angular() + rbI.gInertia()*mv.linear(),
									rbI.massMatrix()*mv.linear() +
									rbI.gInertia().transpose()*mv.angular());
}

inline MotionVec operator*(const PTransform& pt, const MotionVec& mv)
{
	using namespace Eigen;
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	return MotionVec(E*mv.angular(),
									 E*(mv.linear() - r.cross(mv.angular())));
}

inline MotionVec invMul(const PTransform& pt, const MotionVec& mv)
{
	using namespace Eigen;
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	return MotionVec(E.transpose()*mv.angular(),
									 E.transpose()*mv.linear() +
									 r.cross(E.transpose()*mv.angular()));
}

inline ForceVec dualMul(const PTransform& pt, const ForceVec& fv)
{
	using namespace Eigen;
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	return ForceVec(E*(fv.couple() - r.cross(fv.force())),
									E*fv.force());
}

inline ForceVec transMul(const PTransform& pt, const ForceVec& fv)
{
	using namespace Eigen;
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	Vector3d n = E.transpose()*fv.couple();
	n.noalias() += r.cross(E.transpose()*fv.force());
	return ForceVec(n,
									E.transpose()*fv.force());
}

inline RBInertia dualMul(const PTransform& pt, const RBInertia& rbI)
{
	using namespace Eigen;
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

inline RBInertia transMul(const PTransform& pt, const RBInertia& rbI)
{
	using namespace Eigen;
	const Matrix3d& E = pt.rotation();
	const Vector3d& r = pt.translation();
	Matrix3d I;
	I.triangularView<Lower>() = E.transpose()*rbI.inertia()*E -
			vector3ToCrossMatrix(r)*vector3ToCrossMatrix(E.transpose()*rbI.momentum()) -
			vector3ToCrossMatrix(E.transpose()*rbI.momentum() + rbI.mass()*r)*
			vector3ToCrossMatrix(r);
	return RBInertia(rbI.mass(), E.transpose()*rbI.momentum() + rbI.mass()*r, I);
}

inline ABInertia dualMul(const PTransform& pt, const ABInertia& rbI)
{
	using namespace Eigen;
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

inline ABInertia transMul(const PTransform& pt, const ABInertia& rbI)
{
	using namespace Eigen;
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

} // namespace sva
