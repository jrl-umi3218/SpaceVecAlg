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

#pragma once

// includes
// SpaceVecAlg
#include "ABInertia.h"
#include "EigenTypedef.h"
#include "ForceVec.h"
#include "MotionVec.h"
#include "PTransform.h"
#include "RBInertia.h"


// sva::MotionVec Left
sva::MotionVec cross(const sva::MotionVec& mv1, const sva::MotionVec& mv2);
sva::ForceVec crossDual(const sva::MotionVec& mv1, const sva::ForceVec& fv2);
double dot(const sva::MotionVec& mv1, const sva::ForceVec& fv2);

// sva::RBInerita Left
sva::ForceVec operator*(const sva::RBInertia& rbI, const sva::MotionVec& mv);

// sva::ABInertia Left
sva::ABInertia operator+(const sva::ABInertia& abI, const sva::RBInertia& rbI);
sva::ForceVec operator*(const sva::ABInertia& rbI, const sva::MotionVec& mv);

// sva::PTransform Left
sva::MotionVec operator*(const sva::PTransform& pt, const sva::MotionVec& mv);
sva::MotionVec invMul(const sva::PTransform& pt, const sva::MotionVec& mv);

sva::ForceVec dualMul(const sva::PTransform& pt, const sva::ForceVec& fv);
sva::ForceVec transMul(const sva::PTransform& pt, const sva::ForceVec& fv);

sva::RBInertia dualMul(const sva::PTransform& pt, const sva::RBInertia& rbI);
sva::RBInertia transMul(const sva::PTransform& pt, const sva::RBInertia& rbI);

sva::ABInertia dualMul(const sva::PTransform& pt, const sva::ABInertia& rbI);
sva::ABInertia transMul(const sva::PTransform& pt, const sva::ABInertia& rbI);
