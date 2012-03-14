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

// Operators implementation

inline MotionVec MotionVec::cross(const MotionVec& mv2)
{
	return MotionVec(angular().cross(mv2.angular()),
									 angular().cross(mv2.linear()) +
									 linear().cross(mv2.angular()));
}

inline ForceVec MotionVec::crossDual(const ForceVec& fv2)
{
	return ForceVec(angular().cross(fv2.couple()) +
									linear().cross(fv2.force()),
									angular().cross(fv2.force()));
}

inline double MotionVec::dot(const sva::ForceVec& fv2)
{
	return angular().dot(fv2.couple()) + linear().dot(fv2.force());
}

inline ForceVec RBInertia::operator*(const MotionVec& mv)
{
	return ForceVec(inertia()*mv.angular() + momentum().cross(mv.linear()),
									mass()*mv.linear() - momentum().cross(mv.angular()));
}

inline ABInertia ABInertia::operator+(const RBInertia& rbI)
{
	using namespace Eigen;
	Matrix3d M, I;
	M.triangularView<Lower>() = massMatrix() + Matrix3d::Identity()*rbI.mass();
	I.triangularView<Lower>() = inertia() + rbI.inertia();
	return ABInertia(M, gInertia() + vector3ToCrossMatrix(rbI.momentum()), I);
}

inline ForceVec ABInertia::operator*(const MotionVec& mv)
{
	return ForceVec(inertia()*mv.angular() + gInertia()*mv.linear(),
									massMatrix()*mv.linear() +
									gInertia().transpose()*mv.angular());
}

inline MotionVec PTransform::operator*(const MotionVec& mv)
{
	using namespace Eigen;
	const Matrix3d& E = rotation();
	const Vector3d& r = translation();
	return MotionVec(E*mv.angular(),
									 E*(mv.linear() - r.cross(mv.angular())));
}

inline MotionVec PTransform::invMul(const MotionVec& mv)
{
	using namespace Eigen;
	const Matrix3d& E = rotation();
	const Vector3d& r = translation();
	return MotionVec(E.transpose()*mv.angular(),
									 E.transpose()*mv.linear() +
									 r.cross(E.transpose()*mv.angular()));
}

inline ForceVec PTransform::dualMul(const ForceVec& fv)
{
	using namespace Eigen;
	const Matrix3d& E = rotation();
	const Vector3d& r = translation();
	return ForceVec(E*(fv.couple() - r.cross(fv.force())),
									E*fv.force());
}

inline ForceVec PTransform::transMul(const ForceVec& fv)
{
	using namespace Eigen;
	const Matrix3d& E = rotation();
	const Vector3d& r = translation();
	Vector3d n = E.transpose()*fv.couple();
	n.noalias() += r.cross(E.transpose()*fv.force());
	return ForceVec(n,
									E.transpose()*fv.force());
}

inline RBInertia PTransform::dualMul(const RBInertia& rbI)
{
	using namespace Eigen;
	const Matrix3d& E = rotation();
	const Vector3d& r = translation();
	Matrix3d I;
	I.triangularView<Lower>() = E*(rbI.inertia() + vector3ToCrossMatrix(r)*
			vector3ToCrossMatrix(rbI.momentum()) +
			vector3ToCrossMatrix(rbI.momentum() - rbI.mass()*r)*
			vector3ToCrossMatrix(r))*E.transpose();
	return RBInertia(rbI.mass(),
									 E*(rbI.momentum() - rbI.mass()*r),
									 I);
}

inline RBInertia PTransform::transMul(const RBInertia& rbI)
{
	using namespace Eigen;
	const Matrix3d& E = rotation();
	const Vector3d& r = translation();
	Matrix3d I;
	I.triangularView<Lower>() = E.transpose()*rbI.inertia()*E -
			vector3ToCrossMatrix(r)*vector3ToCrossMatrix(E.transpose()*rbI.momentum()) -
			vector3ToCrossMatrix(E.transpose()*rbI.momentum() + rbI.mass()*r)*
			vector3ToCrossMatrix(r);
	return RBInertia(rbI.mass(), E.transpose()*rbI.momentum() + rbI.mass()*r, I);
}

inline ABInertia PTransform::dualMul(const ABInertia& rbI)
{
	using namespace Eigen;
	const Matrix3d& E = rotation();
	const Vector3d& r = translation();

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

inline ABInertia PTransform::transMul(const ABInertia& rbI)
{
	using namespace Eigen;
	const Matrix3d& E = rotation();
	const Vector3d& r = translation();

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
