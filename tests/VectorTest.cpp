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
#include <SpaceVecAlg>

const double TOL = 0.00001;

BOOST_AUTO_TEST_CASE(MotionVecTest)
{
	using namespace Eigen;
	Vector3d w, v;
	Vector6d m;
	w = Vector3d::Random();
	v = Vector3d::Random();

	sva::MotionVec vec(w, v);
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
	Vector6d m2;
	w2 *= 2.;
	v2 *= 3.;
	sva::MotionVec vec2(w2, v2);
	m2 = vec2.vector();

	// M + M
	BOOST_CHECK_EQUAL((vec + vec2).vector(), m + m2);

	// M - M
	BOOST_CHECK_EQUAL((vec - vec2).vector(), m - m2);

	// ==
	BOOST_CHECK_EQUAL(vec, vec);
	BOOST_CHECK_NE(vec, -vec);

	// !=
	BOOST_CHECK(vec != (-vec));
	BOOST_CHECK(!(vec != vec));
}

BOOST_AUTO_TEST_CASE(ForceVecTest)
{
	using namespace Eigen;
	Vector3d n, f;
	Vector6d m;
	n = Vector3d::Random();
	f = Vector3d::Random();

	sva::ForceVec vec(n, f);
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
	Vector6d m2;
	f2 *= 2.;
	n2 *= 3.;
	sva::ForceVec vec2(n2, f2);
	m2 = vec2.vector();

	// F + F
	BOOST_CHECK_EQUAL((vec + vec2).vector(), m + m2);

	// F - F
	BOOST_CHECK_EQUAL((vec - vec2).vector(), m - m2);

	// ==
	BOOST_CHECK_EQUAL(vec, vec);
	BOOST_CHECK_NE(vec, -vec);

	// !=
	BOOST_CHECK(vec != (-vec));
	BOOST_CHECK(!(vec != vec));
}

BOOST_AUTO_TEST_CASE(MotionVecLeftOperatorsTest)
{
	using namespace Eigen;
	Vector3d w, v, n, f;
	w = Vector3d::Random()*100.;
	v = Vector3d::Random()*100.;
	n = Vector3d::Random()*100.;
	f = Vector3d::Random()*100.;

	sva::MotionVec mVec(w, v);
	sva::ForceVec fVec(n, f);

	Vector6d mm, mf;
	mm = mVec.vector();
	mf = fVec.vector();

	// dot(MotionVec, ForceVec)
	BOOST_CHECK_SMALL(mVec.dot(fVec) - mm.transpose()*mf, TOL);

	// cross(MotionVec, MotionVec)
	Vector3d w2, v2;
	w2 = Vector3d::Random()*100.;
	v2 = Vector3d::Random()*100.;

	sva::MotionVec mVec2(w2, v2);
	Vector6d mm2;
	mm2 = mVec2.vector();

	sva::MotionVec crossM = mVec.cross(mVec2);
	BOOST_CHECK_SMALL((crossM.vector() - vector6ToCrossMatrix(mm)*mm2).
										 array().abs().sum(), TOL);

	// crossDual(MotionVec, ForceVec)
	sva::ForceVec crossF = mVec.crossDual(fVec);
	BOOST_CHECK_SMALL((crossF.vector() - vector6ToCrossDualMatrix(mm)*mf).
										 array().abs().sum(), TOL);
}

