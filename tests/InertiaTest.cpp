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
#define BOOST_TEST_MODULE RBInertiad ABinertiad test
#include <boost/test/unit_test.hpp>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

typedef Eigen::Matrix<double, 6, Eigen::Dynamic> Matrix6Xd;

const double TOL = 0.00001;

bool isUpperNull(const Eigen::Matrix3d& m)
{
	using namespace Eigen;
	return (Matrix3d(m.triangularView<StrictlyUpper>()).array() == 0.).all();
}

BOOST_AUTO_TEST_CASE(RBInertiadTest)
{
	using namespace Eigen;
	using namespace sva;

	double mass = 1.;
	Matrix3d I;
	I << 1., 2., 3.,
			 2., 1., 4.,
			 3., 4., 1.;
	Vector3d h = Vector3d::Random()*100.;

	// Parametrized constructor double Vector3d Matrix3d
	RBInertiad rb2(mass, h, I);

	BOOST_CHECK_EQUAL(rb2.mass(), mass);
	BOOST_CHECK_EQUAL(rb2.momentum(), h);
	BOOST_CHECK_EQUAL(rb2.inertia(), I);
	BOOST_CHECK(isUpperNull(rb2.lowerTriangularInertia()));

	// Parametrized constructor double Vector3d Matrix3d
	RBInertiad rb3(mass, h, I.triangularView<Lower>());

	BOOST_CHECK_EQUAL(rb3.mass(), mass);
	BOOST_CHECK_EQUAL(rb3.momentum(), h);
	BOOST_CHECK_EQUAL(rb3.inertia(), I);
	BOOST_CHECK(isUpperNull(rb3.lowerTriangularInertia()));

	// rbI + rbI
	RBInertiad rb4 = rb2 + rb3;

	BOOST_CHECK_EQUAL(rb4.mass(), mass + mass);
	BOOST_CHECK_EQUAL(rb4.momentum(), h + h);
	BOOST_CHECK_EQUAL(rb4.inertia(), I + I);
	BOOST_CHECK(isUpperNull(rb4.lowerTriangularInertia()));

	// alpha * rbI
	RBInertiad rb5 = 2.*rb2;

	BOOST_CHECK_EQUAL(rb5.mass(), 2.*mass);
	BOOST_CHECK_EQUAL(rb5.momentum(), 2.*h);
	BOOST_CHECK_EQUAL(rb5.inertia(), 2.*I);
	BOOST_CHECK(isUpperNull(rb5.lowerTriangularInertia()));

	// rbI * alpha
	RBInertiad rb6 = rb2*2.;

	BOOST_CHECK_EQUAL(rb6.mass(), 2.*mass);
	BOOST_CHECK_EQUAL(rb6.momentum(), 2.*h);
	BOOST_CHECK_EQUAL(rb6.inertia(), 2.*I);
	BOOST_CHECK(isUpperNull(rb6.lowerTriangularInertia()));

	// rbI - rbI
	RBInertiad rb7 = rb2 - rb3;

	BOOST_CHECK_EQUAL(rb7.mass(), mass - mass);
	BOOST_CHECK_EQUAL(rb7.momentum(), h - h);
	BOOST_CHECK_EQUAL(rb7.inertia(), I - I);
	BOOST_CHECK(isUpperNull(rb7.lowerTriangularInertia()));

	// -rbI
	RBInertiad rb8 = -rb2;

	BOOST_CHECK_EQUAL(rb8, rb2*-1.);
	BOOST_CHECK(isUpperNull(rb8.lowerTriangularInertia()));

	// rbI += rbI
	RBInertiad rb9(rb2);
	rb9 += rb3;

	BOOST_CHECK_EQUAL(rb9, rb2 + rb3);
	BOOST_CHECK(isUpperNull(rb9.lowerTriangularInertia()));

	// rbI -= rbI
	RBInertiad rb10(rb2);
	rb10 -= rb3;

	BOOST_CHECK_EQUAL(rb10, rb2 - rb3);
	BOOST_CHECK(isUpperNull(rb10.lowerTriangularInertia()));

	// ==
	BOOST_CHECK_EQUAL(rb2, rb2);
	BOOST_CHECK_NE(rb2, rb6);

	// !=
	BOOST_CHECK(rb2 != rb6);
	BOOST_CHECK(!(rb2 != rb2));
}

BOOST_AUTO_TEST_CASE(ABInertiadTest)
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

	// Parametrized constructor double Vector3d Matrix3d
	ABInertiad ab1(M, H, I);

	BOOST_CHECK_EQUAL(ab1.massMatrix(), M);
	BOOST_CHECK(isUpperNull(ab1.lowerTriangularMassMatrix()));
	BOOST_CHECK_EQUAL(ab1.gInertia(), H);
	BOOST_CHECK_EQUAL(ab1.inertia(), I);
	BOOST_CHECK(isUpperNull(ab1.lowerTriangularInertia()));

	// Parametrized constructor double Vector3d Matrix3d
	ABInertiad ab2(M.triangularView<Lower>(), H, I.triangularView<Lower>());

	BOOST_CHECK_EQUAL(ab2.massMatrix(), M);
	BOOST_CHECK(isUpperNull(ab2.lowerTriangularMassMatrix()));
	BOOST_CHECK_EQUAL(ab2.gInertia(), H);
	BOOST_CHECK_EQUAL(ab2.inertia(), I);
	BOOST_CHECK(isUpperNull(ab2.lowerTriangularInertia()));

	// abI + abI
	ABInertiad ab3 = ab1 + ab2;

	BOOST_CHECK_EQUAL(ab3.massMatrix(), M + M);
	BOOST_CHECK(isUpperNull(ab3.lowerTriangularMassMatrix()));
	BOOST_CHECK_EQUAL(ab3.gInertia(), H + H);
	BOOST_CHECK_EQUAL(ab3.inertia(), I + I);
	BOOST_CHECK(isUpperNull(ab3.lowerTriangularInertia()));

	// alpha * rbI
	ABInertiad ab4 = 2.*ab2;

	BOOST_CHECK_EQUAL(ab4.massMatrix(), 2.*M);
	BOOST_CHECK(isUpperNull(ab4.lowerTriangularMassMatrix()));
	BOOST_CHECK_EQUAL(ab4.gInertia(), 2.*H);
	BOOST_CHECK_EQUAL(ab4.inertia(), 2.*I);
	BOOST_CHECK(isUpperNull(ab4.lowerTriangularInertia()));

	// abI * alpha
	ABInertiad ab5 = ab2*2.;

	BOOST_CHECK_EQUAL(ab5.massMatrix(), 2.*M);
	BOOST_CHECK(isUpperNull(ab5.lowerTriangularMassMatrix()));
	BOOST_CHECK_EQUAL(ab5.gInertia(), 2.*H);
	BOOST_CHECK_EQUAL(ab5.inertia(), 2.*I);
	BOOST_CHECK(isUpperNull(ab5.lowerTriangularInertia()));

	// abI - abI
	ABInertiad ab6 = ab1 - ab2;

	BOOST_CHECK_EQUAL(ab6.massMatrix(), M - M);
	BOOST_CHECK(isUpperNull(ab6.lowerTriangularMassMatrix()));
	BOOST_CHECK_EQUAL(ab6.gInertia(), H - H);
	BOOST_CHECK_EQUAL(ab6.inertia(), I - I);
	BOOST_CHECK(isUpperNull(ab6.lowerTriangularInertia()));

	// -abI
	ABInertiad ab7 = -ab1;
	BOOST_CHECK_EQUAL(ab7, ab1*-1.);
	BOOST_CHECK(isUpperNull(ab7.lowerTriangularMassMatrix()));
	BOOST_CHECK(isUpperNull(ab7.lowerTriangularInertia()));

	// abI += abI
	ABInertiad ab8(ab1);
	ab8 += ab2;

	BOOST_CHECK_EQUAL(ab8, ab1 + ab2);
	BOOST_CHECK(isUpperNull(ab8.lowerTriangularMassMatrix()));
	BOOST_CHECK(isUpperNull(ab8.lowerTriangularInertia()));

	// abI -= abI
	ABInertiad ab9(ab1);
	ab9 -= ab2;

	BOOST_CHECK_EQUAL(ab9, ab1 - ab2);
	BOOST_CHECK(isUpperNull(ab9.lowerTriangularMassMatrix()));
	BOOST_CHECK(isUpperNull(ab9.lowerTriangularInertia()));

	// ==
	BOOST_CHECK_EQUAL(ab2, ab2);
	BOOST_CHECK_NE(ab2, ab5);

	// !=
	BOOST_CHECK(ab2 != ab5);
	BOOST_CHECK(!(ab2 != ab2));
}

BOOST_AUTO_TEST_CASE(RBInertiadLeftOperatorsTest)
{
	using namespace Eigen;
	using namespace sva;
	double mass = 1.;
	Matrix3d I;
	I << 1., 2., 3.,
			 2., 1., 4.,
			 3., 4., 1.;
	Vector3d h = Vector3d::Random()*100.;
	RBInertiad rb(mass, h, I);
	Matrix6d rb6d = rb.matrix();

	Vector3d w, v;
	w = Vector3d::Random()*100.;
	v = Vector3d::Random()*100.;
	sva::MotionVecd mVec(w, v);
	Vector6d mVec6d = mVec.vector();

	// RBInertiad * MotionVecd
	ForceVecd fVec = rb*mVec;
	Vector6d fVec6d(rb6d*mVec6d);

	BOOST_CHECK_SMALL((fVec6d - fVec.vector()).array().abs().sum(), TOL);

	// vectorized version
	Matrix6Xd mVec6Xd(6, 2);
	Matrix6Xd fVecRes6Xd(6, 2);
	mVec6Xd << mVec.vector(), mVec.vector();

	internal::set_is_malloc_allowed(false);
	rb.mul(mVec6Xd, fVecRes6Xd);
	internal::set_is_malloc_allowed(true);

#ifdef __i386__
	BOOST_CHECK_SMALL((fVec.vector() - fVecRes6Xd.col(0)).array().abs().sum(), TOL);
#else
	BOOST_CHECK_EQUAL(fVec.vector(), fVecRes6Xd.col(0));
#endif
	BOOST_CHECK_EQUAL(fVecRes6Xd.col(0), fVecRes6Xd.col(1));
}

BOOST_AUTO_TEST_CASE(ABInertiadLeftOperatorsTest)
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

	ABInertiad ab(M, H, I);
	Matrix6d ab6d = ab.matrix();

	double mass = 1.;
	Vector3d h = Vector3d::Random()*100.;
	RBInertiad rb(mass, h, I);
	Matrix6d rb6d = rb.matrix();

	Vector3d w, v;
	w = Vector3d::Random()*100.;
	v = Vector3d::Random()*100.;
	sva::MotionVecd mVec(w, v);
	Vector6d mVec6d = mVec.vector();

	// ABInertiad + RBInertiad
	ABInertiad abRes = ab + rb;
	Matrix6d abRes6d = ab6d + rb6d;

	BOOST_CHECK_SMALL((abRes6d - abRes.matrix()).array().abs().sum(), TOL);
	BOOST_CHECK(isUpperNull(abRes.lowerTriangularMassMatrix()));
	BOOST_CHECK(isUpperNull(abRes.lowerTriangularInertia()));

	// ABInertiad * MotionVecd
	ForceVecd fVec = ab*mVec;
	Vector6d fVec6d(ab6d*mVec6d);

	BOOST_CHECK_SMALL((fVec6d - fVec.vector()).array().abs().sum(), TOL);

	// vectorized version
	Matrix6Xd mVec6Xd(6, 2);
	Matrix6Xd fVecRes6Xd(6, 2);
	mVec6Xd << mVec.vector(), mVec.vector();

	internal::set_is_malloc_allowed(false);
	ab.mul(mVec6Xd, fVecRes6Xd);
	internal::set_is_malloc_allowed(true);

#ifdef __i386__
	BOOST_CHECK_SMALL((fVec.vector() - fVecRes6Xd.col(0)).array().abs().sum(), TOL);
#else
	BOOST_CHECK_EQUAL(fVec.vector(), fVecRes6Xd.col(0));
#endif
	BOOST_CHECK_EQUAL(fVecRes6Xd.col(0), fVecRes6Xd.col(1));
}

