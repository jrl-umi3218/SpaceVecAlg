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
#include <SpaceVecAlg>

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

	BOOST_CHECK_EQUAL(RotX(theta), AngleAxisd(-theta, Vector3d::UnitX()).matrix());
	BOOST_CHECK_EQUAL(RotY(theta), AngleAxisd(-theta, Vector3d::UnitY()).matrix());
	BOOST_CHECK_EQUAL(RotZ(theta), AngleAxisd(-theta, Vector3d::UnitZ()).matrix());
}

BOOST_AUTO_TEST_CASE(PTransformdTest)
{
	using namespace Eigen;
	using namespace sva;
	namespace constants = boost::math::constants;

	Matrix3d Em = AngleAxisd(constants::pi<double>()/2., Vector3d(1., 0., 0.)).inverse().toRotationMatrix();
	Quaterniond Eq;
	Eq = AngleAxisd(constants::pi<double>()/2., Vector3d(1., 0., 0.));
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

	BOOST_CHECK_EQUAL(pt3.rotation(), Eq.inverse().toRotationMatrix());
	BOOST_CHECK_EQUAL(pt3.translation(), r);

	// Quaternion constructor
	PTransformd pt4(Eq);

	BOOST_CHECK_EQUAL(pt4.rotation(), Eq.inverse().toRotationMatrix());
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
