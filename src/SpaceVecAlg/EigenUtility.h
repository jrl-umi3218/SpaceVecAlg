/*
 * Copyright 2012-2020 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "EigenTypedef.h"

namespace Eigen
{

/**
 * Convert a 3D vector to a cross product matrix.
 */
template<typename T>
inline Matrix3<T> vector3ToCrossMatrix(const Vector3<T> & vec)
{
  Matrix3<T> mat;
  mat << 0., -vec(2), vec(1), vec(2), 0., -vec(0), -vec(1), vec(0), 0.;
  return mat;
}

/**
 * Convert a 6D vector to a spatial cross product matrix.
 */
template<typename T>
inline Matrix6<T> vector6ToCrossMatrix(const Vector6<T> & vec)
{
  Matrix6<T> mat;
  Matrix3<T> c13 = vector3ToCrossMatrix<T>(vec.template head<3>());
  mat << c13, Matrix3d::Zero(), vector3ToCrossMatrix<T>(vec.template tail<3>()), c13;
  return mat;
}

/**
 * Convert a 6D vector to a spatial dual cross product matrix.
 */
template<typename T>
inline Matrix6<T> vector6ToCrossDualMatrix(const Vector6<T> & vec)
{
  return -vector6ToCrossMatrix<T>(vec).transpose();
}

} // namespace Eigen

namespace sva
{

using Eigen::vector3ToCrossMatrix;
using Eigen::vector6ToCrossDualMatrix;
using Eigen::vector6ToCrossMatrix;

} // namespace sva
