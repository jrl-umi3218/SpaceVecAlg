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

// Don't include it directly, include SpaceVecAlg instead

namespace sva
{

// Operators implementation

namespace sva_internal
{


template<typename Derived1, typename Derived2, typename Derived3>
inline void colwiseCrossEq(const Eigen::MatrixBase<Derived1>& m1,
													const Eigen::MatrixBase<Derived2>& m2,
													Eigen::MatrixBase<Derived3> const & result)
{
	Eigen::MatrixBase<Derived3>& result_nc = const_cast<Eigen::MatrixBase<Derived3> &>(result);
	result_nc.row(0) = (m1.row(1)*m2.coeff(2) - m1.row(2)*m2.coeff(1));
	result_nc.row(1) = (m1.row(2)*m2.coeff(0) - m1.row(0)*m2.coeff(2));
	result_nc.row(2) = (m1.row(0)*m2.coeff(1) - m1.row(1)*m2.coeff(0));
}

template<typename Derived1, typename Derived2, typename Derived3>
inline void colwiseCrossPlusEq(const Eigen::MatrixBase<Derived1>& m1,
														 const Eigen::MatrixBase<Derived2>& m2,
														 Eigen::MatrixBase<Derived3> const & result)
{
	Eigen::MatrixBase<Derived3>& result_nc = const_cast<Eigen::MatrixBase<Derived3> &>(result);
	result_nc.row(0) += (m1.row(1)*m2.coeff(2) - m1.row(2)*m2.coeff(1));
	result_nc.row(1) += (m1.row(2)*m2.coeff(0) - m1.row(0)*m2.coeff(2));
	result_nc.row(2) += (m1.row(0)*m2.coeff(1) - m1.row(1)*m2.coeff(0));
}

template<typename Derived1, typename Derived2, typename Derived3>
inline void colwiseCrossMinusEq(const Eigen::MatrixBase<Derived1>& m1,
															const Eigen::MatrixBase<Derived2>& m2,
															Eigen::MatrixBase<Derived3> const & result)
{
	Eigen::MatrixBase<Derived3>& result_nc = const_cast<Eigen::MatrixBase<Derived3> &>(result);
	result_nc.row(0) -= (m1.row(1)*m2.coeff(2) - m1.row(2)*m2.coeff(1));
	result_nc.row(1) -= (m1.row(2)*m2.coeff(0) - m1.row(0)*m2.coeff(2));
	result_nc.row(2) -= (m1.row(0)*m2.coeff(1) - m1.row(1)*m2.coeff(0));
}

template<typename Derived1, typename Derived2, typename Derived3>
inline void colwiseLeftMultEq(const Eigen::MatrixBase<Derived1>& m1,
														const Eigen::MatrixBase<Derived2>& m2,
														Eigen::MatrixBase<Derived3> const & result)
{
	Eigen::MatrixBase<Derived3>& result_nc = const_cast<Eigen::MatrixBase<Derived3> &>(result);
	for(typename Derived1::Index i = 0; i < m1.cols(); ++i)
	{
		result_nc.col(i) = m2*m1.col(i);
	}
}

} // internal


template<typename Derived>
inline Eigen::Block<Derived, 3, Dynamic>
motionAngular(Eigen::MatrixBase<Derived>& mv)
{
	return Eigen::Block<Derived, 3, Dynamic>(mv.derived(), 0, 0, 3, mv.cols());
}

template<typename Derived>
inline Eigen::Block<const Derived, 3, Dynamic>
motionAngular(const Eigen::MatrixBase<Derived>& mv)
{
	return Eigen::Block<const Derived, 3, Dynamic>(mv.derived(), 0, 0, 3, mv.cols());
}

template<typename Derived>
inline Eigen::Block<Derived, 3, Dynamic>
motionLinear(Eigen::MatrixBase<Derived>& mv)
{
	return Eigen::Block<Derived, 3, Dynamic>(mv.derived(), 3, 0, 3, mv.cols());
}

template<typename Derived>
inline Eigen::Block<const Derived, 3, Dynamic>
motionLinear(const Eigen::MatrixBase<Derived>& mv)
{
	return Eigen::Block<const Derived, 3, Dynamic>(mv.derived(), 3, 0, 3, mv.cols());
}

template<typename Derived>
inline Eigen::Block<Derived, 3, Dynamic>
forceCouple(Eigen::MatrixBase<Derived>& mv)
{
	return Eigen::Block<Derived, 3, Dynamic>(mv.derived(), 0, 0, 3, mv.cols());
}

template<typename Derived>
inline Eigen::Block<const Derived, 3, Dynamic>
forceCouple(const Eigen::MatrixBase<Derived>& mv)
{
	return Eigen::Block<const Derived, 3, Dynamic>(mv.derived(), 0, 0, 3, mv.cols());
}

template<typename Derived>
inline Eigen::Block<Derived, 3, Dynamic>
forceForce(Eigen::MatrixBase<Derived>& mv)
{
	return Eigen::Block<Derived, 3, Dynamic>(mv.derived(), 3, 0, 3, mv.cols());
}

template<typename Derived>
inline Eigen::Block<const Derived, 3, Dynamic>
forceForce(const Eigen::MatrixBase<Derived>& mv)
{
	return Eigen::Block<const Derived, 3, Dynamic>(mv.derived(), 3, 0, 3, mv.cols());
}

template<typename T>
inline MotionVec<T> MotionVec<T>::cross(const MotionVec<T>& mv2) const
{
	return MotionVec<T>(angular().cross(mv2.angular()),
										 angular().cross(mv2.linear()) +
										 linear().cross(mv2.angular()));
}

template<typename T>
template<typename Derived>
inline void MotionVec<T>::cross(const Eigen::MatrixBase<Derived>& mv2,
	Eigen::MatrixBase<Derived>& result) const
{
	using namespace sva_internal;
	static_assert(Derived::RowsAtCompileTime == 6,
							 "the matrix must have exactly 6 rows");
	static_assert(std::is_same<typename Derived::Scalar, T>::value,
							 "motion vec and matrix must be the same type");

	colwiseCrossEq(motionAngular(mv2), -angular_, motionAngular(result));

	colwiseCrossEq(motionLinear(mv2), -angular_, motionLinear(result));
	colwiseCrossMinusEq(motionAngular(mv2), linear_, motionLinear(result));
}

template<typename T>
inline ForceVec<T> MotionVec<T>::crossDual(const ForceVec<T>& fv2) const
{
	return ForceVec<T>(angular().cross(fv2.couple()) +
										linear().cross(fv2.force()),
										angular().cross(fv2.force()));
}

template<typename T>
template<typename Derived>
inline void MotionVec<T>::crossDual(const Eigen::MatrixBase<Derived>& fv2,
	Eigen::MatrixBase<Derived>& result) const
{
	using namespace sva_internal;
	static_assert(Derived::RowsAtCompileTime == 6,
							 "the matrix must have exactly 6 rows");
	static_assert(std::is_same<typename Derived::Scalar, T>::value,
							 "motion vec and matrix must be the same type");

	colwiseCrossEq(forceCouple(fv2), -angular_, forceCouple(result));
	colwiseCrossMinusEq(forceForce(fv2), linear_, forceCouple(result));

	colwiseCrossEq(forceForce(fv2), -angular_, forceForce(result));
}

template<typename T>
inline T MotionVec<T>::dot(const sva::ForceVec<T>& fv2) const
{
	return angular().dot(fv2.couple()) + linear().dot(fv2.force());
}

template<typename T>
inline ForceVec<T> RBInertia<T>::operator*(const MotionVec<T>& mv) const
{
	return ForceVec<T>(inertia()*mv.angular() + momentum().cross(mv.linear()),
										mass()*mv.linear() - momentum().cross(mv.angular()));
}

template<typename T>
template<typename Derived>
inline void RBInertia<T>::mul(const Eigen::MatrixBase<Derived>& mv,
	Eigen::MatrixBase<Derived>& result) const
{
	using namespace sva_internal;
	static_assert(Derived::RowsAtCompileTime == 6,
							 "the matrix must have exactly 6 rows");
	static_assert(std::is_same<typename Derived::Scalar, T>::value,
							 "motion vec and matrix must be the same type");

	forceCouple(result).noalias() = inertia()*motionAngular(mv);
	colwiseCrossMinusEq(motionLinear(mv), momentum(), forceCouple(result));

	forceForce(result).noalias() = motionLinear(mv)*mass();
	colwiseCrossPlusEq(motionAngular(mv), momentum(), forceForce(result));
}

template<typename T>
inline ABInertia<T> ABInertia<T>::operator+(const RBInertia<T>& rbI) const
{
	using namespace Eigen;
	Matrix3<T> M, I;
	M.template triangularView<Lower>() = massMatrix() + Matrix3<T>::Identity()*rbI.mass();
	I.template triangularView<Lower>() = inertia() + rbI.inertia();
	return ABInertia<T>(M, gInertia() + vector3ToCrossMatrix(rbI.momentum()), I);
}

template<typename T>
inline ForceVec<T> ABInertia<T>::operator*(const MotionVec<T>& mv) const
{
	return ForceVec<T>(inertia()*mv.angular() + gInertia()*mv.linear(),
									massMatrix()*mv.linear() +
									gInertia().transpose()*mv.angular());
}

template<typename T>
template<typename Derived>
inline void ABInertia<T>::mul(const Eigen::MatrixBase<Derived>& mv,
	Eigen::MatrixBase<Derived>& result) const
{
	static_assert(Derived::RowsAtCompileTime == 6,
							 "the matrix must have exactly 6 rows");
	static_assert(std::is_same<typename Derived::Scalar, T>::value,
							 "motion vec and matrix must be the same type");

	forceCouple(result).noalias() = inertia()*motionAngular(mv);
	forceCouple(result).noalias() += gInertia()*motionLinear(mv);

	forceForce(result).noalias() = inertia()*motionLinear(mv);
	forceForce(result).noalias() += gInertia().transpose()*motionAngular(mv);
}

template<typename T>
inline MotionVec<T> PTransform<T>::operator*(const MotionVec<T>& mv) const
{
	using namespace Eigen;
	return MotionVec<T>(angularMul(mv), linearMul(mv));
}

template<typename T>
inline Eigen::Vector3<T> PTransform<T>::angularMul(const MotionVec<T>& mv) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	return (E*mv.angular()).eval();
}

template<typename T>
inline Eigen::Vector3<T> PTransform<T>::linearMul(const MotionVec<T>& mv) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();
	return (E*(mv.linear() - r.cross(mv.angular()))).eval();
}

template<typename T>
template<typename Derived>
inline void PTransform<T>::mul(const Eigen::MatrixBase<Derived>& mv,
	Eigen::MatrixBase<Derived>& result) const
{
	using namespace Eigen;
	using namespace sva_internal;
	static_assert(Derived::RowsAtCompileTime == 6,
							 "the matrix must have exactly 6 rows");
	static_assert(std::is_same<typename Derived::Scalar, T>::value,
							 "motion vec and matrix must be the same type");

	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();

	motionAngular(result).noalias() = E*motionAngular(mv);

	motionLinear(result).noalias() = motionLinear(mv);
	colwiseCrossPlusEq(motionAngular(mv), r, motionLinear(result));
	colwiseLeftMultEq(motionLinear(result), E, motionLinear(result));
}

template<typename T>
inline MotionVec<T> PTransform<T>::invMul(const MotionVec<T>& mv) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();
	MotionVec<T> ret;
	ret.angular_.noalias() = E.transpose()*mv.angular_;
	ret.linear_.noalias() = E.transpose()*mv.linear_;
	ret.linear_.noalias() += r.cross(ret.angular_);
	return ret;
}

template<typename T>
inline Eigen::Vector3<T> PTransform<T>::angularInvMul(const MotionVec<T>& mv) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	return (E.transpose()*mv.angular()).eval();
}

template<typename T>
inline Eigen::Vector3<T> PTransform<T>::linearInvMul(const MotionVec<T>& mv) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();
	return (E.transpose()*mv.linear() + r.cross(E.transpose()*mv.angular())).eval();
}

template<typename T>
template<typename Derived>
inline void PTransform<T>::invMul(const Eigen::MatrixBase<Derived>& mv,
	Eigen::MatrixBase<Derived>& result) const
{
	using namespace Eigen;
	using namespace sva_internal;
	static_assert(Derived::RowsAtCompileTime == 6,
							 "the matrix must have exactly 6 rows");
	static_assert(std::is_same<typename Derived::Scalar, T>::value,
							 "motion vec and matrix must be the same type");

	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();

	motionAngular(result).noalias() = E.transpose()*motionAngular(mv);

	motionLinear(result).noalias() = E.transpose()*motionLinear(mv);
	colwiseCrossMinusEq(motionAngular(result), r, motionLinear(result));
}

template<typename T>
inline ForceVec<T> PTransform<T>::dualMul(const ForceVec<T>& fv) const
{
	using namespace Eigen;
	return ForceVec<T>(coupleDualMul(fv), forceDualMul(fv));
}

template<typename T>
inline Eigen::Vector3<T> PTransform<T>::coupleDualMul(const ForceVec<T>& fv) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();
	return (E*(fv.couple() - r.cross(fv.force()))).eval();
}

template<typename T>
inline Eigen::Vector3<T> PTransform<T>::forceDualMul(const ForceVec<T>& fv) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	return (E*fv.force()).eval();
}

template<typename T>
template<typename Derived>
inline void PTransform<T>::dualMul(const Eigen::MatrixBase<Derived>& fv,
	Eigen::MatrixBase<Derived>& result) const
{
	using namespace Eigen;
	using namespace sva_internal;
	static_assert(Derived::RowsAtCompileTime == 6,
							 "the matrix must have exactly 6 rows");
	static_assert(std::is_same<typename Derived::Scalar, T>::value,
							 "force vec and matrix must be the same type");

	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();

	forceCouple(result).noalias() = forceCouple(fv);
	colwiseCrossPlusEq(forceForce(fv), r, forceCouple(result));
	colwiseLeftMultEq(forceCouple(result), E, forceCouple(result));

	forceForce(result).noalias() = E*forceForce(fv);
}

template<typename T>
inline ForceVec<T> PTransform<T>::transMul(const ForceVec<T>& fv) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();
	ForceVec<T> ret;
	ret.force_.noalias() = E.transpose()*fv.force_;
	ret.couple_.noalias() = E.transpose()*fv.couple_;
	ret.couple_.noalias() += r.cross(ret.force_);
	return ret;
}

template<typename T>
inline Eigen::Vector3<T> PTransform<T>::coupleTransMul(const ForceVec<T>& fv) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();

	return (E.transpose()*fv.couple_ + r.cross(E.transpose()*fv.force_)).eval();
}

template<typename T>
inline Eigen::Vector3<T> PTransform<T>::forceTransMul(const ForceVec<T>& fv) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	return (E.transpose()*fv.force_).eval();
}

template<typename T>
template<typename Derived>
inline void PTransform<T>::transMul(const Eigen::MatrixBase<Derived>& fv,
	Eigen::MatrixBase<Derived>& result) const
{
	using namespace Eigen;
	using namespace sva_internal;
	static_assert(Derived::RowsAtCompileTime == 6,
							 "the matrix must have exactly 6 rows");
	static_assert(std::is_same<typename Derived::Scalar, T>::value,
							 "force vec and matrix must be the same type");

	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();

	forceForce(result).noalias() = E.transpose()*forceForce(fv);

	forceCouple(result).noalias() = E.transpose()*forceCouple(fv);
	colwiseCrossMinusEq(forceForce(result), r, forceCouple(result));
}

template<typename T>
inline RBInertia<T> PTransform<T>::dualMul(const RBInertia<T>& rbI) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();
	Matrix3<T> I;
	I.template triangularView<Lower>() = E*(rbI.inertia() + vector3ToCrossMatrix(r)*
			vector3ToCrossMatrix<T>(rbI.momentum()) +
			vector3ToCrossMatrix<T>(rbI.momentum() - rbI.mass()*r)*
			vector3ToCrossMatrix<T>(r))*E.transpose();
	return RBInertia<T>(rbI.mass(),
										 E*(rbI.momentum() - rbI.mass()*r),
										 I);
}

template<typename T>
inline RBInertia<T> PTransform<T>::transMul(const RBInertia<T>& rbI) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();
	Matrix3<T> I;
	I.template triangularView<Lower>() = E.transpose()*rbI.inertia()*E -
			vector3ToCrossMatrix<T>(r)*vector3ToCrossMatrix<T>(E.transpose()*rbI.momentum()) -
			vector3ToCrossMatrix<T>(E.transpose()*rbI.momentum() + rbI.mass()*r)*
			vector3ToCrossMatrix<T>(r);
	return RBInertia<T>(rbI.mass(), E.transpose()*rbI.momentum() + rbI.mass()*r, I);
}

template<typename T>
inline ABInertia<T> PTransform<T>::dualMul(const ABInertia<T>& rbI) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();

	Matrix3<T> massM =rbI.massMatrix();
	Matrix3<T> rCross = vector3ToCrossMatrix(r);
	Matrix3<T> tmpI = rbI.gInertia() - rCross*massM;

	Matrix3<T> M, I;
	M.template triangularView<Lower>() = E*massM*E.transpose();
	I.template triangularView<Lower>() = E*(rbI.inertia() - rCross*rbI.gInertia().transpose() +
			(tmpI*rCross))*E.transpose();

	return ABInertia<T>(M,
									 E*(tmpI)*E.transpose(),
									 I);
}

template<typename T>
inline ABInertia<T> PTransform<T>::transMul(const ABInertia<T>& rbI) const
{
	using namespace Eigen;
	const Matrix3<T>& E = rotation();
	const Vector3<T>& r = translation();

	Matrix3<T> Mp(E.transpose()*rbI.massMatrix()*E);
	Matrix3<T> Hp(E.transpose()*rbI.gInertia()*E);
	Matrix3<T> rCross(vector3ToCrossMatrix(r));

	Matrix3<T> M, I;
	M.template triangularView<Lower>() = Mp;
	I.template triangularView<Lower>() = (E.transpose()*rbI.inertia()*E +
																			rCross*Hp.transpose() -
																			(Hp + rCross*Mp)*rCross);
	return ABInertia<T>(M,
									 Hp + rCross*Mp,
									 I);

}

} // namespace sva
