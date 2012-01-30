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

// SpaceVecAlg
#include <ABInertia.h>
#include <MotionVec.h>
#include <Operators.h>
#include <RBInertia.h>

const double TOL = 0.00001;

bool isUpperNull(const Eigen::Matrix3d& m)
{
	using namespace Eigen;
	return (Matrix3d(m.triangularView<StrictlyUpper>()).array() == 0.).all();
}

BOOST_AUTO_TEST_CASE(RBInertiaTest)
{
	using namespace Eigen;
	using namespace sva;

	double mass = 1.;
	Matrix3d I;
	I << 1., 2., 3.,
			 2., 1., 4.,
			 3., 4., 1.;
	Vector3d h = Vector3d::Random()*100.;

	// Default constructor
	RBInertia rb1;

	BOOST_CHECK_EQUAL(rb1.mass(), 0.);

	// Parametrized constructor double Vector3d Matrix3d
	RBInertia rb2(mass, h, I);

	BOOST_CHECK_EQUAL(rb2.mass(), mass);
	BOOST_CHECK_EQUAL(rb2.momentum(), h);
	BOOST_CHECK_EQUAL(rb2.inertia(), I);
	BOOST_CHECK(isUpperNull(rb2.lowerTriangularInertia()));

	// Parametrized constructor double Vector3d Matrix3d
	RBInertia rb3(mass, h, I.triangularView<Lower>());

	BOOST_CHECK_EQUAL(rb3.mass(), mass);
	BOOST_CHECK_EQUAL(rb3.momentum(), h);
	BOOST_CHECK_EQUAL(rb3.inertia(), I);
	BOOST_CHECK(isUpperNull(rb3.lowerTriangularInertia()));

	// rbI + rbI
	RBInertia rb4 = rb2 + rb3;

	BOOST_CHECK_EQUAL(rb4.mass(), mass + mass);
	BOOST_CHECK_EQUAL(rb4.momentum(), h + h);
	BOOST_CHECK_EQUAL(rb4.inertia(), I + I);
	BOOST_CHECK(isUpperNull(rb4.lowerTriangularInertia()));

	// alpha * rbI
	RBInertia rb5 = 2.*rb2;

	BOOST_CHECK_EQUAL(rb5.mass(), 2.*mass);
	BOOST_CHECK_EQUAL(rb5.momentum(), 2.*h);
	BOOST_CHECK_EQUAL(rb5.inertia(), 2.*I);
	BOOST_CHECK(isUpperNull(rb5.lowerTriangularInertia()));
}

BOOST_AUTO_TEST_CASE(ABInertiaTest)
{
	using namespace Eigen;
	using namespace sva;

	double mass = 1.;
	Matrix3d M, H, I;
	M << 1., 2., 3.,
			 2., 1., 4.,
			 3., 4., 1.;
	H = Matrix3d::Random()*100.;
	I << 1., 2., 3.,
			 2., 1., 4.,
			 3., 4., 1.;

	// Parametrized constructor double Vector3d Matrix3d
	ABInertia ab1(M, H, I);

	BOOST_CHECK_EQUAL(ab1.massMatrix(), M);
	BOOST_CHECK(isUpperNull(ab1.lowerTriangularMassMatrix()));
	BOOST_CHECK_EQUAL(ab1.gInertia(), H);
	BOOST_CHECK_EQUAL(ab1.inertia(), I);
	BOOST_CHECK(isUpperNull(ab1.lowerTriangularInertia()));

	// Parametrized constructor double Vector3d Matrix3d
	ABInertia ab2(M.triangularView<Lower>(), H, I.triangularView<Lower>());

	BOOST_CHECK_EQUAL(ab2.massMatrix(), M);
	BOOST_CHECK(isUpperNull(ab2.lowerTriangularMassMatrix()));
	BOOST_CHECK_EQUAL(ab2.gInertia(), H);
	BOOST_CHECK_EQUAL(ab2.inertia(), I);
	BOOST_CHECK(isUpperNull(ab2.lowerTriangularInertia()));

	// abI + abI
	ABInertia ab3 = ab1 + ab2;

	BOOST_CHECK_EQUAL(ab3.massMatrix(), M + M);
	BOOST_CHECK(isUpperNull(ab3.lowerTriangularMassMatrix()));
	BOOST_CHECK_EQUAL(ab3.gInertia(), H + H);
	BOOST_CHECK_EQUAL(ab3.inertia(), I + I);
	BOOST_CHECK(isUpperNull(ab3.lowerTriangularInertia()));

	// alpha * rbI
	ABInertia ab4 = 2.*ab2;

	BOOST_CHECK_EQUAL(ab4.massMatrix(), 2.*M);
	BOOST_CHECK(isUpperNull(ab4.lowerTriangularMassMatrix()));
	BOOST_CHECK_EQUAL(ab4.gInertia(), 2.*H);
	BOOST_CHECK_EQUAL(ab4.inertia(), 2.*I);
	BOOST_CHECK(isUpperNull(ab4.lowerTriangularInertia()));
}

BOOST_AUTO_TEST_CASE(RBInertiaLeftOperatorsTest)
{
	using namespace Eigen;
	using namespace sva;
	double mass = 1.;
	Matrix3d I;
	I << 1., 2., 3.,
			 2., 1., 4.,
			 3., 4., 1.;
	Vector3d h = Vector3d::Random()*100.;
	RBInertia rb(mass, h, I);
	Matrix6d rb6d = rb.matrix();

	Vector3d w, v;
	w = Vector3d::Random()*100.;
	v = Vector3d::Random()*100.;
	sva::MotionVec mVec(w, v);
	Vector6d mVec6d = mVec.vector();

	// RBInertia * MotionVec
	ForceVec fVec = rb*mVec;
	Vector6d fVec6d = rb6d*mVec6d;

	BOOST_CHECK_SMALL((fVec6d - fVec.vector()).array().abs().sum(), TOL);
}

BOOST_AUTO_TEST_CASE(ABInertiaLeftOperatorsTest)
{
	using namespace Eigen;
	using namespace sva;
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

	Vector3d w, v;
	w = Vector3d::Random()*100.;
	v = Vector3d::Random()*100.;
	sva::MotionVec mVec(w, v);
	Vector6d mVec6d = mVec.vector();

	// ABInertia + RBInertia
	ABInertia abRes = ab + rb;
	Matrix6d abRes6d = ab6d + rb6d;

	BOOST_CHECK_SMALL((abRes6d - abRes.matrix()).array().abs().sum(), TOL);
	BOOST_CHECK(isUpperNull(abRes.lowerTriangularMassMatrix()));
	BOOST_CHECK(isUpperNull(abRes.lowerTriangularInertia()));

	// ABInertia * MotionVec
	ForceVec fVec = ab*mVec;
	Vector6d fVec6d = ab6d*mVec6d;

	BOOST_CHECK_SMALL((fVec6d - fVec.vector()).array().abs().sum(), TOL);
}

