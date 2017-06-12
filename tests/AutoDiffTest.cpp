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

// includes
// std
#include <iostream>

// boost
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE AutoDiff test
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// Eigen
#include <Eigen/Core>
#include <unsupported/Eigen/AutoDiff>

// AutoDiff
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> derivative_t;
typedef Eigen::AutoDiffScalar<derivative_t> scalar_t;

/*
namespace std
{
scalar_t sin(const scalar_t& t)
{
  return Eigen::sin(t);
}
scalar_t cos(const scalar_t& t)
{
  return Eigen::cos(t);
}
}
*/

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>
namespace sva
{
template<>
inline Matrix3<scalar_t> RotZ<scalar_t>(scalar_t theta)
{
	using namespace Eigen;

	Matrix3<scalar_t> ret;
	scalar_t z(0., derivative_t::Zero(theta.derivatives().rows()));
	scalar_t o(1., derivative_t::Zero(theta.derivatives().rows()));
	scalar_t s = sin(theta), c = cos(theta);

	ret << c, s, z,
				-s, c, z,
				z, z, o;

	return ret;
}
}

using boost::math::constants::pi;

BOOST_AUTO_TEST_CASE(PTransformVsPTransform)
{
	using namespace sva;
	using namespace Eigen;

	{
		Vector3<scalar_t> pt1T(0., 1., 0.);
		pt1T[0].derivatives().setZero(2);
		pt1T[1].derivatives().setZero(2);
		pt1T[1].derivatives()(0) = 1.;
		pt1T[2].derivatives().setZero(2);

		PTransform<scalar_t> pt1(pt1T);
		PTransformd pt2(Vector3d(1., 0., 0.));

		PTransform<scalar_t> res1 = pt1*pt1;
		PTransform<scalar_t> res2 = pt1*PTransform<scalar_t>{pt2};
		BOOST_CHECK(PTransform<scalar_t>{pt2} == pt2.cast<scalar_t>());

		std::cout << res1.translation()[0].derivatives().transpose() << std::endl;
		std::cout << res1.translation()[1].derivatives().transpose() << std::endl;
		std::cout << res1.translation()[2].derivatives().transpose() << std::endl;

		std::cout << std::endl;

		std::cout << res2.translation()[0].derivatives().transpose() << std::endl;
		std::cout << res2.translation()[1].derivatives().transpose() << std::endl;
		std::cout << res2.translation()[2].derivatives().transpose() << std::endl;

		std::cout << std::endl;
	}

	{
		Vector3<scalar_t> prismY(0., 1., 0.);
		prismY[0].derivatives().setZero(2, 1);
		prismY[1].derivatives().setZero(2, 1);
		prismY[1].derivatives()(0) = 1.;
		prismY[2].derivatives().setZero(2, 1);

		scalar_t rotZ = pi<double>()/4.;
		rotZ.derivatives().setZero(2, 1);
		rotZ.derivatives()(1) = 1.;

		// PTransform<scalar_t> Rot(RotZ(rotZ));
		PTransform<scalar_t> Rot(AngleAxis<scalar_t>(-rotZ, Vector3<scalar_t>::UnitZ()).matrix());
		PTransform<scalar_t> Prism(prismY);

		PTransform<scalar_t> res = Prism*Rot;

		std::cout << res.translation()[0].derivatives().transpose() << std::endl;
		std::cout << res.translation()[1].derivatives().transpose() << std::endl;
		std::cout << res.translation()[2].derivatives().transpose() << std::endl;

		std::cout << std::endl;
	}
}
