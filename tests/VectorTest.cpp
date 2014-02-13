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
#define BOOST_TEST_MODULE MotionVecd ForceVecd test
#include <boost/test/unit_test.hpp>

// SpaceVecAlg
#include <SpaceVecAlg>

typedef Eigen::Matrix<double, 6, Eigen::Dynamic> Matrix6Xd;

const double TOL = 0.00001;

BOOST_AUTO_TEST_CASE(MotionVecdTest)
{
	using namespace Eigen;
	Vector3d w, v;
	Vector6d m;
	w = Vector3d::Random();
	v = Vector3d::Random();

	sva::MotionVecd vec(w, v);
	m = vec.vector();

	// angular
	BOOST_CHECK_EQUAL(w, vec.angular());

	// linear
	BOOST_CHECK_EQUAL(v, vec.linear());

	// vector
	BOOST_CHECK_EQUAL(m, (Vector6d() << w, v).finished());

	// alpha*M
	BOOST_CHECK_EQUAL((5.*vec).vector(), 5.*m);

	// M*alpha
	BOOST_CHECK_EQUAL((vec*5.).vector(), 5.*m);

	// -M
	BOOST_CHECK_EQUAL((-vec).vector(), -m);

	Vector3d w2, v2;
	w2 = Vector3d::Random();
	v2 = Vector3d::Random();
	Vector6d m2;
	sva::MotionVecd vec2(w2, v2);
	m2 = vec2.vector();

	// M + M
	BOOST_CHECK_EQUAL((vec + vec2).vector(), m + m2);

	// M - M
	BOOST_CHECK_EQUAL((vec - vec2).vector(), m - m2);

	// M += M
	sva::MotionVecd vec_pluseq(vec);
	vec_pluseq += vec2;
	BOOST_CHECK_EQUAL(vec_pluseq, vec + vec2);

	// M -= M
	sva::MotionVecd vec_minuseq(vec);
	vec_minuseq -= vec2;
	BOOST_CHECK_EQUAL(vec_minuseq, vec - vec2);

	// ==
	BOOST_CHECK_EQUAL(vec, vec);
	BOOST_CHECK_NE(vec, -vec);

	// !=
	BOOST_CHECK(vec != (-vec));
	BOOST_CHECK(!(vec != vec));
}

BOOST_AUTO_TEST_CASE(ForceVecdTest)
{
	using namespace Eigen;
	Vector3d n, f;
	Vector6d m;
	n = Vector3d::Random();
	f = Vector3d::Random();

	sva::ForceVecd vec(n, f);
	m = vec.vector();

	// couple
	BOOST_CHECK_EQUAL(n, vec.couple());

	// force
	BOOST_CHECK_EQUAL(f, vec.force());

	// vector
	BOOST_CHECK_EQUAL(m, (Vector6d() << n, f).finished());

	// alpha*F
	BOOST_CHECK_EQUAL((5.*vec).vector(), 5.*m);

	// F*alpha
	BOOST_CHECK_EQUAL((vec*5.).vector(), 5.*m);

	// -F
	BOOST_CHECK_EQUAL((-vec).vector(), -m);

	Vector3d n2, f2;
	n2 = Vector3d::Random();
	f2 = Vector3d::Random();
	Vector6d m2;
	sva::ForceVecd vec2(n2, f2);
	m2 = vec2.vector();

	// F + F
	BOOST_CHECK_EQUAL((vec + vec2).vector(), m + m2);

	// F - F
	BOOST_CHECK_EQUAL((vec - vec2).vector(), m - m2);

	// M += M
	sva::ForceVecd vec_pluseq(vec);
	vec_pluseq += vec2;
	BOOST_CHECK_EQUAL(vec_pluseq, vec + vec2);

	// M -= M
	sva::ForceVecd vec_minuseq(vec);
	vec_minuseq -= vec2;
	BOOST_CHECK_EQUAL(vec_minuseq, vec - vec2);

	// ==
	BOOST_CHECK_EQUAL(vec, vec);
	BOOST_CHECK_NE(vec, -vec);

	// !=
	BOOST_CHECK(vec != (-vec));
	BOOST_CHECK(!(vec != vec));
}

BOOST_AUTO_TEST_CASE(MotionVecdLeftOperatorsTest)
{
	using namespace Eigen;
	Vector3d w, v, n, f;
	w = Vector3d::Random()*100.;
	v = Vector3d::Random()*100.;
	n = Vector3d::Random()*100.;
	f = Vector3d::Random()*100.;

	sva::MotionVecd mVec(w, v);
	sva::ForceVecd fVec(n, f);

	Vector6d mm, mf;
	mm = mVec.vector();
	mf = fVec.vector();

	// dot(MotionVecd, ForceVecd)
	BOOST_CHECK_SMALL(mVec.dot(fVec) - mm.transpose()*mf, TOL);

	// cross(MotionVecd, MotionVecd)
	Vector3d w2, v2;
	w2 = Vector3d::Random()*100.;
	v2 = Vector3d::Random()*100.;

	sva::MotionVecd mVec2(w2, v2);
	Vector6d mm2;
	mm2 = mVec2.vector();

	sva::MotionVecd crossM = mVec.cross(mVec2);
	BOOST_CHECK_SMALL((crossM.vector() - vector6ToCrossMatrix(mm)*mm2).
										 array().abs().sum(), TOL);

	// crossDual(MotionVecd, ForceVecd)
	sva::ForceVecd crossF = mVec.crossDual(fVec);
	BOOST_CHECK_SMALL((crossF.vector() - vector6ToCrossDualMatrix(mm)*mf).
										 array().abs().sum(), TOL);

	// test the vectorized version
	Matrix6Xd crossMVec(6, 2);
	Matrix6Xd crossMVecRes(6, 2);
	crossMVec << mVec2.vector(), mVec2.vector();

	internal::set_is_malloc_allowed(false);
	mVec.cross(crossMVec, crossMVecRes);
	internal::set_is_malloc_allowed(true);

	BOOST_CHECK_EQUAL(crossM.vector(), crossMVecRes.col(0));
	BOOST_CHECK_EQUAL(crossMVecRes.col(0), crossMVecRes.col(1));

	Matrix6Xd crossFVec(6, 2);
	Matrix6Xd crossFVecRes(6, 2);
	crossFVec << fVec.vector(), fVec.vector();

	internal::set_is_malloc_allowed(false);
	mVec.crossDual(crossFVec, crossFVecRes);
	internal::set_is_malloc_allowed(true);

	BOOST_CHECK_EQUAL(crossF.vector(), crossFVecRes.col(0));
	BOOST_CHECK_EQUAL(crossFVecRes.col(0), crossFVecRes.col(1));
}

