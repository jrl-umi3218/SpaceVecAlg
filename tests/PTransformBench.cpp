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
#define BOOST_TEST_MODULE PTranformd bench
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/timer/timer.hpp>
// The inclusion of boost chrono was commented in timer.hpp for boost >= 1.60.
// Because of this, the auto-link feature does not incude the chrono library
// anymore, what causes a link error. 
// (see also https://svn.boost.org/trac/boost/ticket/11862)
// We add manually the line.
// Possible alternative: include only for specific version of boost and 
// auto-link capable compiler
#include <boost/chrono/chrono.hpp>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

typedef Eigen::Matrix<double, 6, Eigen::Dynamic> Matrix6Xd;

BOOST_AUTO_TEST_CASE(PTransfromd_PTransformd)
{
	using namespace sva;

	const std::size_t size = 10000000;
	std::vector<PTransformd> pt1(size, PTransformd::Identity());
	std::vector<PTransformd> pt2(size, PTransformd::Identity());
	std::vector<PTransformd> ptRes(size);


	std::cout << "PTransform vs PTransform" << std::endl;
	{
		boost::timer::auto_cpu_timer t;
		for(std::size_t i = 0; i < size; ++i)
		{
			ptRes[i] = pt1[i]*pt2[i];
		}
	}
	std::cout << std::endl;
}


BOOST_AUTO_TEST_CASE(PTransfromd_MotionVec)
{
	using namespace sva;

	const std::size_t size = 10000000;
	std::vector<PTransformd> pt1(size, PTransformd::Identity());
	std::vector<MotionVecd> mv(size, MotionVecd(Eigen::Vector6d::Random()));
	std::vector<MotionVecd> mvRes(size);


	std::cout << "PTransform vs MotionVec" << std::endl;
	{
		boost::timer::auto_cpu_timer t;
		for(std::size_t i = 0; i < size; ++i)
		{
			mvRes[i] = pt1[i]*mv[i];
		}
	}
	std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(PTransfromd_MotionEigen)
{
	using namespace sva;

	const std::size_t size = 10000000;
	const std::size_t cols = 3;
	std::vector<PTransformd> pt1(size, PTransformd::Identity());
	std::vector<Matrix6Xd> mv(size, Matrix6Xd::Random(6,cols));
	std::vector<Matrix6Xd> mvRes(size, Matrix6Xd(6,cols));


	std::cout << "PTransform vs MotionEigen" << std::endl;
	{
		boost::timer::auto_cpu_timer t;
		for(std::size_t i = 0; i < size; ++i)
		{
			pt1[i].mul(mv[i], mvRes[i]);
		}
	}
	std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(PTransfromd_as_matrix_MotionEigen)
{
	using namespace sva;

	const std::size_t size = 10000000;
	const std::size_t cols = 3;
	std::vector<PTransformd> pt1(size, PTransformd::Identity());
	std::vector<Matrix6Xd> mv(size, Matrix6Xd::Random(6,cols));
	std::vector<Matrix6Xd> mvRes(size, Matrix6Xd(6,cols));

	std::cout << "PTransform as matrix vs MotionEigen" << std::endl;
	{
		boost::timer::auto_cpu_timer t;
		for(std::size_t i = 0; i < size; ++i)
		{
			mvRes[i].noalias() = pt1[i].matrix()*mv[i];
		}
	}
	std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(PTransfromd_MotionEigen_as_motion)
{
	using namespace sva;

	const std::size_t size = 10000000;
	const std::size_t cols = 3;
	std::vector<PTransformd> pt1(size, PTransformd::Identity());
	std::vector<Matrix6Xd> mv(size, Matrix6Xd::Random(6,cols));
	std::vector<Matrix6Xd> mvRes(size, Matrix6Xd(6,cols));

	std::cout << "PTransform vs MotionEigen as MotionVec" << std::endl;
	{
		boost::timer::auto_cpu_timer t;
		for(std::size_t i = 0; i < size; ++i)
		{
			for(std::size_t j = 0; j < cols; ++j)
			{
				mvRes[i].col(j).noalias() = (pt1[i]*MotionVecd(mv[i].col(j))).vector();
			}
		}
	}
	std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(PTransfromd_inv_MotionVec)
{
	using namespace sva;

	const std::size_t size = 10000000;
	std::vector<PTransformd> pt1(size, PTransformd::Identity());
	std::vector<MotionVecd> mv(size, MotionVecd(Eigen::Vector6d::Random()));
	std::vector<MotionVecd> mvRes(size);


	std::cout << "PTransform_inv vs MotionVec" << std::endl;
	{
		boost::timer::auto_cpu_timer t;
		for(std::size_t i = 0; i < size; ++i)
		{
			mvRes[i] = pt1[i].inv()*mv[i];
		}
	}
	std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(PTransfromd_invMul_MotionVec)
{
	using namespace sva;

	const std::size_t size = 10000000;
	std::vector<PTransformd> pt1(size, PTransformd::Identity());
	std::vector<MotionVecd> mv(size, MotionVecd(Eigen::Vector6d::Random()));
	std::vector<MotionVecd> mvRes(size);


	std::cout << "PTransform_invMul vs MotionVec" << std::endl;
	{
		boost::timer::auto_cpu_timer t;
		for(std::size_t i = 0; i < size; ++i)
		{
			mvRes[i] = pt1[i].invMul(mv[i]);
		}
	}
	std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(PTransfromd_dual_ForceVec)
{
	using namespace sva;

	const std::size_t size = 10000000;
	std::vector<PTransformd> pt1(size, PTransformd::Identity());
	std::vector<ForceVecd> mv(size, ForceVecd(Eigen::Vector6d::Random()));
	std::vector<ForceVecd> mvRes(size);


	std::cout << "PTransform dual vs ForceVec" << std::endl;
	{
		boost::timer::auto_cpu_timer t;
		for(std::size_t i = 0; i < size; ++i)
		{
			mvRes[i] = pt1[i].dualMul(mv[i]);
		}
	}
	std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(PTransfromd_trans_ForceVec)
{
	using namespace sva;

	const std::size_t size = 10000000;
	std::vector<PTransformd> pt1(size, PTransformd::Identity());
	std::vector<ForceVecd> mv(size, ForceVecd(Eigen::Vector6d::Random()));
	std::vector<ForceVecd> mvRes(size);


	std::cout << "PTransform trans vs ForceVec" << std::endl;
	{
		boost::timer::auto_cpu_timer t;
		for(std::size_t i = 0; i < size; ++i)
		{
			mvRes[i] = pt1[i].transMul(mv[i]);
		}
	}
	std::cout << std::endl;
}
