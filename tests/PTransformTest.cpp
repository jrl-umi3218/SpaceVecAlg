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

// check memory allocation in some method
#define EIGEN_RUNTIME_NO_MALLOC

// includes
// std
#include <iostream>

// boost
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PTransformd test
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

typedef Eigen::Matrix<double, 6, Eigen::Dynamic> Matrix6Xd;

const double TOL = 0.00001;

bool isUpperNull(const Eigen::Matrix3d& m)
{
	using namespace Eigen;
	return (Matrix3d(m.triangularView<StrictlyUpper>()).array() == 0.).all();
}

BOOST_AUTO_TEST_CASE(RotationMatrixTest)
{
	using namespace Eigen;
	using namespace sva;

	Vector2d theta2d = Vector2d::Random()*10;
	double theta = theta2d(0);

	BOOST_CHECK_SMALL((RotX(theta) - AngleAxisd(-theta, Vector3d::UnitX()).matrix()).array().abs().sum(), TOL);
	BOOST_CHECK_SMALL((RotY(theta) - AngleAxisd(-theta, Vector3d::UnitY()).matrix()).array().abs().sum(), TOL);
	BOOST_CHECK_SMALL((RotZ(theta) - AngleAxisd(-theta, Vector3d::UnitZ()).matrix()).array().abs().sum(), TOL);
}

BOOST_AUTO_TEST_CASE(PTransformdTest)
{
	using namespace Eigen;
	using namespace sva;
	namespace constants = boost::math::constants;

	Matrix3d Em = AngleAxisd(constants::pi<double>()/2., Vector3d(1., 0., 0.)).inverse().toRotationMatrix();
	Quaterniond Eq;
	Eq = AngleAxisd(constants::pi<double>()/2., Vector3d(1., 0., 0.)).inverse();
	Vector3d r = Vector3d::Random()*100.;

	// Identity
	PTransformd pt1 = PTransformd::Identity();

	BOOST_CHECK_EQUAL(pt1.rotation(), Matrix3d::Identity());
	BOOST_CHECK_EQUAL(pt1.translation(), Vector3d::Zero());

	// Matrix3d Vector3d constructor
	PTransformd pt2(Em, r);

	BOOST_CHECK_EQUAL(pt2.rotation(), Em);
	BOOST_CHECK_EQUAL(pt2.translation(), r);

	// Quaternion Vector3d constructor
	PTransformd pt3(Eq, r);

	BOOST_CHECK_EQUAL(pt3.rotation(), Eq.toRotationMatrix());
	BOOST_CHECK_EQUAL(pt3.translation(), r);

	// Quaternion constructor
	PTransformd pt4(Eq);

	BOOST_CHECK_EQUAL(pt4.rotation(), Eq.toRotationMatrix());
	BOOST_CHECK_EQUAL(pt4.translation(), Vector3d::Zero());

	// Matrix3d constructor
	PTransformd pt5(Em);

	BOOST_CHECK_EQUAL(pt5.rotation(), Em);
	BOOST_CHECK_EQUAL(pt5.translation(), Vector3d::Zero());

	// Vector3d constructor
	PTransformd pt6(r);

	BOOST_CHECK_EQUAL(pt6.rotation(), Matrix3d::Identity());
	BOOST_CHECK_EQUAL(pt6.translation(), r);

	// operator*(PTransformd)
	PTransformd pttmp(AngleAxisd(constants::pi<double>()/4., Vector3d(0.,1.,0.)).
															toRotationMatrix(),
									 Vector3d::Random()*100.);

	PTransformd pt7 = pt2*pttmp;
	Matrix6d ptm(pt2.matrix()*pttmp.matrix());

	BOOST_CHECK_SMALL((pt7.matrix() - ptm).array().abs().sum(), TOL);

	// inv
	PTransformd pt8 = pt2.inv();
	BOOST_CHECK_SMALL((pt8.matrix() - pt2.matrix().inverse()).array().abs().sum(), TOL);

	// ==
	BOOST_CHECK_EQUAL(pt2, pt2);
	BOOST_CHECK_NE(pt2, pt8);

	// !=
	BOOST_CHECK(pt2 != pt8);
	BOOST_CHECK(!(pt2 != pt2));
}

BOOST_AUTO_TEST_CASE(PTransformdLeftOperatorsTest)
{
	using namespace Eigen;
	using namespace sva;
	namespace constants = boost::math::constants;

	Quaterniond Eq;
	Eq = AngleAxisd(constants::pi<double>()/2., Vector3d(1., 0., 0.));
	Vector3d r = Vector3d::Random()*100.;
	PTransformd pt(Eq, r);
	PTransformd ptInv= pt.inv();
	Matrix6d pt6d = pt.matrix();
	Matrix6d ptInv6d = ptInv.matrix();
	Matrix6d ptDual6d = pt.dualMatrix();

	Matrix3d M, H, I;
	M << 1., 2., 3.,
			 2., 1., 4.,
			 3., 4., 1.;
	H = Matrix3d::Random()*100.;
	I << 1., 2., 3.,
			 2., 1., 4.,
			 3., 4., 1.;

	ABInertiad ab(M, H, I);
	Matrix6d ab6d = ab.matrix();

	double mass = 1.;
	Vector3d h = Vector3d::Random()*100.;
	RBInertiad rb(mass, h, I);
	Matrix6d rb6d = rb.matrix();

	Vector3d w, v, n, f;
	w = Vector3d::Random()*100.;
	v = Vector3d::Random()*100.;
	n = Vector3d::Random()*100.;
	f = Vector3d::Random()*100.;

	sva::MotionVecd mVec(w, v);
	Vector6d mVec6d = mVec.vector();

	sva::ForceVecd fVec(n, f);
	Vector6d fVec6d = fVec.vector();

	// PTransformd * MotionVecd
	MotionVecd mvRes1 = pt*mVec;
	Vector6d mvRes16d(pt6d*mVec6d);

	BOOST_CHECK_SMALL((mvRes1.vector() - mvRes16d).array().abs().sum(), TOL);

	// test the angular and linear version
	BOOST_CHECK_SMALL((mvRes1.angular() - pt.angularMul(mVec)).array().abs().sum(), TOL);
	BOOST_CHECK_SMALL((mvRes1.linear() - pt.linearMul(mVec)).array().abs().sum(), TOL);

	// test the vectorized version
	Matrix6Xd mv1Vec6Xd(6, 2);
	Matrix6Xd mvRes1Vec6Xd(6, 2);
	mv1Vec6Xd << mVec.vector(), mVec.vector();

	internal::set_is_malloc_allowed(false);
	pt.mul(mv1Vec6Xd, mvRes1Vec6Xd);
	internal::set_is_malloc_allowed(true);

	BOOST_CHECK_SMALL((mvRes1.vector() - mvRes1Vec6Xd.col(0)).norm(), TOL);
	BOOST_CHECK_EQUAL(mvRes1Vec6Xd.col(0), mvRes1Vec6Xd.col(1));

	// PTransformd^-1 * MotionVecd
	MotionVecd mvRes2 = pt.invMul(mVec);
	Vector6d mvRes26d(ptInv6d*mVec6d);

	BOOST_CHECK_SMALL((mvRes2.vector() - mvRes26d).array().abs().sum(), TOL);

	// test the angular and linear version
	BOOST_CHECK_SMALL((mvRes2.angular() - pt.angularInvMul(mVec)).array().abs().sum(), TOL);
	BOOST_CHECK_SMALL((mvRes2.linear() - pt.linearInvMul(mVec)).array().abs().sum(), TOL);

	// test the vectorized version
	Matrix6Xd mv2Vec6Xd(6, 2);
	Matrix6Xd mvRes2Vec6Xd(6, 2);
	mv2Vec6Xd << mVec.vector(), mVec.vector();

	internal::set_is_malloc_allowed(false);
	pt.invMul(mv2Vec6Xd, mvRes2Vec6Xd);
	internal::set_is_malloc_allowed(true);

	BOOST_CHECK_SMALL((mvRes2.vector() - mvRes2Vec6Xd.col(0)).norm(), TOL);
	BOOST_CHECK_EQUAL(mvRes2Vec6Xd.col(0), mvRes2Vec6Xd.col(1));

	// PTransformd* * ForceVecd
	ForceVecd fvRes1 = pt.dualMul(fVec);
	Vector6d fvRes16d(ptDual6d*fVec6d);

	BOOST_CHECK_SMALL((fvRes1.vector() - fvRes16d).array().abs().sum(), TOL);

	// test the couple and force version
	BOOST_CHECK_SMALL((fvRes1.couple() - pt.coupleDualMul(fVec)).array().abs().sum(), TOL);
	BOOST_CHECK_SMALL((fvRes1.force() - pt.forceDualMul(fVec)).array().abs().sum(), TOL);

	// test the vectorized version
	Matrix6Xd fv1Vec6Xd(6, 2);
	Matrix6Xd fvRes1Vec6Xd(6, 2);
	fv1Vec6Xd << fVec.vector(), fVec.vector();

	internal::set_is_malloc_allowed(false);
	pt.dualMul(fv1Vec6Xd, fvRes1Vec6Xd);
	internal::set_is_malloc_allowed(true);

	BOOST_CHECK_SMALL((fvRes1.vector() - fvRes1Vec6Xd.col(0)).norm(), TOL);
	BOOST_CHECK_EQUAL(fvRes1Vec6Xd.col(0), fvRes1Vec6Xd.col(1));

	// PTransformd T * ForceVecd
	ForceVecd fvRes2 = pt.transMul(fVec);
	Vector6d fvRes26d(pt6d.transpose()*fVec6d);

	BOOST_CHECK_SMALL((fvRes2.vector() - fvRes26d).array().abs().sum(), TOL);

	// test the couple and force version
	BOOST_CHECK_SMALL((fvRes2.couple() - pt.coupleTransMul(fVec)).array().abs().sum(), TOL);
	BOOST_CHECK_SMALL((fvRes2.force() - pt.forceTransMul(fVec)).array().abs().sum(), TOL);

	// test the vectorized version
	Matrix6Xd fv2Vec6Xd(6, 2);
	Matrix6Xd fvRes2Vec6Xd(6, 2);
	fv2Vec6Xd << fVec.vector(), fVec.vector();

	internal::set_is_malloc_allowed(false);
	pt.transMul(fv2Vec6Xd, fvRes2Vec6Xd);
	internal::set_is_malloc_allowed(true);

	BOOST_CHECK_SMALL((fvRes2.vector() - fvRes2Vec6Xd.col(0)).norm(), TOL);
	BOOST_CHECK_EQUAL(fvRes2Vec6Xd.col(0), fvRes2Vec6Xd.col(1));

	// PTransformd* * RBInertiad * PTransformd^-1
	RBInertiad rbRes1 = pt.dualMul(rb);
	Matrix6d rbRes16d(ptDual6d*rb6d*ptInv6d);

	BOOST_CHECK_SMALL((rbRes1.matrix() - rbRes16d).array().abs().sum(), TOL);

	// PTransformd T * RBInertiad * PTransformd
	RBInertiad rbRes2 = pt.transMul(rb);
	Matrix6d rbRes26d(pt6d.transpose()*rb6d*pt6d);

	BOOST_CHECK_SMALL((rbRes2.matrix() - rbRes26d).array().abs().sum(), TOL);

	// PTransformd* * ABInertiad * PTransformd^-1
	ABInertiad abRes1 = pt.dualMul(ab);
	Matrix6d abRes16d(ptDual6d*ab6d*ptInv6d);

	BOOST_CHECK_SMALL((abRes1.matrix() - abRes16d).array().abs().sum(), TOL);

	// PTransformd T * ABInertiad * PTransformd
	ABInertiad abRes2 = pt.transMul(ab);
	Matrix6d abRes26d(pt6d.transpose()*ab6d*pt6d);

	BOOST_CHECK_SMALL((abRes2.matrix() - abRes26d).array().abs().sum(), TOL);
}

BOOST_AUTO_TEST_CASE(EulerAngleTest)
{
	using namespace Eigen;
	using namespace sva;
	namespace cst = boost::math::constants;

	Vector3d res;

	res = rotationError<double>(Matrix3d::Identity(), RotX(cst::pi<double>()/2.));
	BOOST_CHECK_SMALL((res - Vector3d(cst::pi<double>()/2., 0., 0.)).norm(), TOL);

	res = rotationError<double>(Matrix3d::Identity(), RotY(cst::pi<double>()/2.));
	BOOST_CHECK_SMALL((res - Vector3d(0., cst::pi<double>()/2., 0.)).norm(), TOL);

	res = rotationError<double>(Matrix3d::Identity(), RotZ(cst::pi<double>()/2.));
	BOOST_CHECK_SMALL((res - Vector3d(0., 0., cst::pi<double>()/2.)).norm(), TOL);

	res = rotationError<double>(RotZ(cst::pi<double>()/4.), RotZ(cst::pi<double>()/2.));
	BOOST_CHECK_SMALL((res - Vector3d(0., 0., cst::pi<double>()/4.)).norm(), TOL);
}

BOOST_AUTO_TEST_CASE(InterpolateTest)
{
	using namespace Eigen;
	using namespace sva;
	namespace cst = boost::math::constants;

	PTransformd from(Matrix3d::Identity(), Vector3d(0., 0., 0.));
	PTransformd to(AngleAxisd(cst::pi<double>(), Vector3d::UnitZ()).toRotationMatrix(),
								 Vector3d(1., 2., -3.));

	PTransformd res = interpolate<double>(from, to, 0.5);

	BOOST_CHECK_SMALL((res.rotation() -
		AngleAxisd(cst::pi<double>()/2., Vector3d::UnitZ()).toRotationMatrix()).norm(), TOL);
	BOOST_CHECK_SMALL((res.translation() - Vector3d(0.5, 1., -1.5)).norm(), TOL);

	res = interpolate<double>(from, to, 0);
	BOOST_CHECK_SMALL((res.rotation() - from.rotation()).norm(), TOL);
	BOOST_CHECK_SMALL((res.translation() - from.translation()).norm(), TOL);

	res = interpolate<double>(from, to, 1);
	BOOST_CHECK_SMALL((res.rotation() - to.rotation()).norm(), TOL);
	BOOST_CHECK_SMALL((res.translation() - to.translation()).norm(), TOL);
}

BOOST_AUTO_TEST_CASE(TransformError)
{
	using namespace Eigen;
	using namespace sva;
	namespace cst = boost::math::constants;

	PTransformd X_a_b(Quaterniond(Vector4d::Random()).normalized(), Vector3d::Random());
	PTransformd X_a_c(Quaterniond(Vector4d::Random()).normalized(), Vector3d::Random());

	MotionVecd V_a_b = transformVelocity(X_a_b);
	MotionVecd V_a_c = transformVelocity(X_a_c);

	BOOST_CHECK_SMALL((V_a_b.angular() - rotationVelocity(X_a_b.rotation())).norm(), TOL);
	BOOST_CHECK_SMALL((V_a_b.linear() - X_a_b.translation()).norm(), TOL);
	BOOST_CHECK_SMALL((V_a_c.angular() - rotationVelocity(X_a_c.rotation())).norm(), TOL);
	BOOST_CHECK_SMALL((V_a_c.linear() - X_a_c.translation()).norm(), TOL);

	MotionVecd V_b_c_a = transformError(X_a_b, X_a_c);
	Vector3d w_b_c_a = rotationError(X_a_b.rotation(), X_a_c.rotation());
	Vector3d v_b_c_a = X_a_c.translation() - X_a_b.translation();

	BOOST_CHECK_SMALL((V_b_c_a.angular() - w_b_c_a).norm(), TOL);
	BOOST_CHECK_SMALL((V_b_c_a.linear() - v_b_c_a).norm(), TOL);
}


BOOST_AUTO_TEST_CASE(sinc_invTest)
{
	auto dummy_sinc_inv = [](double x){return x/std::sin(x);};
	double eps = std::numeric_limits<double>::epsilon();

	// test equality between -1 and 1 (avoid 0)
	double t = -1.;
	const int nrIter = 333;
	for(int i = 0; i < nrIter; ++i)
	{
		BOOST_CHECK_EQUAL(dummy_sinc_inv(t), sva::sinc_inv(t));
		t += 2./nrIter;
	}
	BOOST_CHECK(std::isnan(dummy_sinc_inv(0.)));
	BOOST_CHECK_EQUAL(sva::sinc_inv(0.), 1.);

	// not sure thoses test will work on all architectures
	BOOST_CHECK_EQUAL(dummy_sinc_inv(eps), sva::sinc_inv(eps));
	BOOST_CHECK_EQUAL(dummy_sinc_inv(std::sqrt(eps)),
									 sva::sinc_inv(std::sqrt(eps)));
	BOOST_CHECK_EQUAL(dummy_sinc_inv(std::sqrt(std::sqrt(eps))),
									 sva::sinc_inv(std::sqrt(std::sqrt(eps))));
}


template<typename T>
inline Eigen::Vector3<T> oldRotationVelocity(const Eigen::Matrix3<T>& E_a_b, double prec)
{
	Eigen::Vector3<T> w;
	T acosV = (E_a_b(0,0) + E_a_b(1,1) + E_a_b(2,2) - 1.)*0.5;
	T theta = std::acos(acosV);

	if(E_a_b.isIdentity(prec))
	{
		w.setZero();
	}
	else
	{
		w = Eigen::Vector3<T>(-E_a_b(2,1) + E_a_b(1,2),
									-E_a_b(0,2) + E_a_b(2,0),
									-E_a_b(1,0) + E_a_b(0,1));
		w *= theta/(2.*std::sin(theta));
}

	return w;
}


BOOST_AUTO_TEST_CASE(oldVsNewRotationVelocity)
{
	using namespace Eigen;
	using namespace sva;

	Matrix3d r1(RotZ(1.)*RotX(1.5)*RotZ(2.));
	Matrix3d r2(RotZ(1.043)*RotY(0.3422)*RotX(-0.30943));
	Matrix3d r3(RotX(-0.8348)*RotY(-0.2344)*RotZ(0.2344));
	Matrix3d r4(Matrix3d::Identity());

	BOOST_CHECK_EQUAL(oldRotationVelocity(r1, 1e-7), rotationVelocity(r1));
	BOOST_CHECK_EQUAL(oldRotationVelocity(r2, 1e-7), rotationVelocity(r2));
	BOOST_CHECK_EQUAL(oldRotationVelocity(r3, 1e-7), rotationVelocity(r3));
	BOOST_CHECK_EQUAL(oldRotationVelocity(r4, 1e-7), rotationVelocity(r4));
}
