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

// includes
// std
#include <iostream>

// boost
#define BOOST_TEST_MODULE MotionVec ForceVec test
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// SpaceVecAlg
#include <SpaceVecAlg>

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

BOOST_AUTO_TEST_CASE(PTransformTest)
{
	using namespace Eigen;
	using namespace sva;
	namespace constants = boost::math::constants;

	Matrix3d Em = AngleAxisd(constants::pi<double>()/2., Vector3d(1., 0., 0.)).inverse().toRotationMatrix();
	Quaterniond Eq;
	Eq = AngleAxisd(constants::pi<double>()/2., Vector3d(1., 0., 0.));
	Vector3d r = Vector3d::Random()*100.;

	// Identity
	PTransform pt1 = PTransform::Identity();

	BOOST_CHECK_EQUAL(pt1.rotation(), Matrix3d::Identity());
	BOOST_CHECK_EQUAL(pt1.translation(), Vector3d::Zero());

	// Matrix3d Vector3d constructor
	PTransform pt2(Em, r);

	BOOST_CHECK_EQUAL(pt2.rotation(), Em);
	BOOST_CHECK_EQUAL(pt2.translation(), r);

	// Quaternion Vector3d constructor
	PTransform pt3(Eq, r);

	BOOST_CHECK_EQUAL(pt3.rotation(), Eq.inverse().toRotationMatrix());
	BOOST_CHECK_EQUAL(pt3.translation(), r);

	// Quaternion constructor
	PTransform pt4(Eq);

	BOOST_CHECK_EQUAL(pt4.rotation(), Eq.inverse().toRotationMatrix());
	BOOST_CHECK_EQUAL(pt4.translation(), Vector3d::Zero());

	// Matrix3d constructor
	PTransform pt5(Em);

	BOOST_CHECK_EQUAL(pt5.rotation(), Em);
	BOOST_CHECK_EQUAL(pt5.translation(), Vector3d::Zero());

	// Vector3d constructor
	PTransform pt6(r);

	BOOST_CHECK_EQUAL(pt6.rotation(), Matrix3d::Identity());
	BOOST_CHECK_EQUAL(pt6.translation(), r);

	// operator*(PTransform)
	PTransform pttmp(AngleAxisd(constants::pi<double>()/4., Vector3d(0.,1.,0.)).
															toRotationMatrix(),
									 Vector3d::Random()*100.);

	PTransform pt7 = pt2*pttmp;
	Matrix6d ptm = pt2.matrix()*pttmp.matrix();

	BOOST_CHECK_SMALL((pt7.matrix() - ptm).array().abs().sum(), TOL);

	// inv
	PTransform pt8 = pt2.inv();
	BOOST_CHECK_SMALL((pt8.matrix() - pt2.matrix().inverse()).array().abs().sum(), TOL);

	// ==
	BOOST_CHECK_EQUAL(pt2, pt2);
	BOOST_CHECK_NE(pt2, pt8);

	// !=
	BOOST_CHECK(pt2 != pt8);
	BOOST_CHECK(!(pt2 != pt2));
}

BOOST_AUTO_TEST_CASE(PTransformLeftOperatorsTest)
{
	using namespace Eigen;
	using namespace sva;
	namespace constants = boost::math::constants;

	Quaterniond Eq;
	Eq = AngleAxisd(constants::pi<double>()/2., Vector3d(1., 0., 0.));
	Vector3d r = Vector3d::Random()*100.;
	PTransform pt(Eq, r);
	PTransform ptInv= pt.inv();
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

	ABInertia ab(M, H, I);
	Matrix6d ab6d = ab.matrix();

	double mass = 1.;
	Vector3d h = Vector3d::Random()*100.;
	RBInertia rb(mass, h, I);
	Matrix6d rb6d = rb.matrix();

	Vector3d w, v, n, f;
	w = Vector3d::Random()*100.;
	v = Vector3d::Random()*100.;
	n = Vector3d::Random()*100.;
	f = Vector3d::Random()*100.;

	sva::MotionVec mVec(w, v);
	Vector6d mVec6d = mVec.vector();

	sva::ForceVec fVec(n, f);
	Vector6d fVec6d = fVec.vector();

	// PTransform * MotionVec
	MotionVec mvRes1 = pt*mVec;
	Vector6d mvRes16d = pt6d*mVec6d;

	BOOST_CHECK_SMALL((mvRes1.vector() - mvRes16d).array().abs().sum(), TOL);

	// PTransform^-1 * MotionVec
	MotionVec mvRes2 = pt.invMul(mVec);
	Vector6d mvRes26d = ptInv6d*mVec6d;

	BOOST_CHECK_SMALL((mvRes2.vector() - mvRes26d).array().abs().sum(), TOL);

	// PTransform* * ForceVec
	ForceVec fvRes1 = pt.dualMul(fVec);
	Vector6d fvRes16d = ptDual6d*fVec6d;

	BOOST_CHECK_SMALL((fvRes1.vector() - fvRes16d).array().abs().sum(), TOL);

	// PTransform T * ForceVec
	ForceVec fvRes2 = pt.transMul(fVec);
	Vector6d fvRes26d = pt6d.transpose()*fVec6d;

	BOOST_CHECK_SMALL((fvRes2.vector() - fvRes26d).array().abs().sum(), TOL);

	// PTransform* * RBInertia * PTransform^-1
	RBInertia rbRes1 = pt.dualMul(rb);
	Matrix6d rbRes16d = ptDual6d*rb6d*ptInv6d;

	BOOST_CHECK_SMALL((rbRes1.matrix() - rbRes16d).array().abs().sum(), TOL);

	// PTransform T * RBInertia * PTransform
	RBInertia rbRes2 = pt.transMul(rb);
	Matrix6d rbRes26d = pt6d.transpose()*rb6d*pt6d;

	BOOST_CHECK_SMALL((rbRes2.matrix() - rbRes26d).array().abs().sum(), TOL);

	// PTransform* * ABInertia * PTransform^-1
	ABInertia abRes1 = pt.dualMul(ab);
	Matrix6d abRes16d = ptDual6d*ab6d*ptInv6d;

	BOOST_CHECK_SMALL((abRes1.matrix() - abRes16d).array().abs().sum(), TOL);

	// PTransform T * ABInertia * PTransform
	ABInertia abRes2 = pt.transMul(ab);
	Matrix6d abRes26d = pt6d.transpose()*ab6d*pt6d;

	BOOST_CHECK_SMALL((abRes2.matrix() - abRes26d).array().abs().sum(), TOL);
}
