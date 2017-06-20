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
#include <SpaceVecAlg/Conversions.h>

BOOST_AUTO_TEST_CASE(ConversionsHomogeneous)
{
    using namespace Eigen;
    using namespace sva;

    Matrix4d hom = Matrix4d::Zero();
    hom(0,1) = -1;
    hom(1,0) = +1;
    hom(2,2) = +1;
    hom(3,3) = +1;
    PTransform<double> pt = conversions::fromHomogeneous(hom,conversions::RightHanded);

    const Vector4d vec(4,3,2,1);
    Matrix4d homVec = Matrix4d::Identity();
    homVec.block<4,1>(0,3) = vec;
    PTransform<double> ptVec = conversions::fromHomogeneous(homVec,conversions::RightHanded);

    Matrix4d hom2 = conversions::toHomogeneous(pt,conversions::RightHanded);

    const PTransform<double> rotatedPT= ptVec*pt;
    const Vector4d rotated = hom*vec;
    const Vector4d rotated2 = hom2*vec;

    BOOST_CHECK_EQUAL(rotated, rotated2);
    BOOST_CHECK_EQUAL((rotated.block<3,1>(0,0)), rotatedPT.translation());
}
