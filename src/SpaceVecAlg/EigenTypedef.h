/*
 * Copyright 2012-2020 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <Eigen/Core>

namespace Eigen
{

template<typename T>
using Vector6 = Matrix<T, 6, 1>;
template<typename T>
using Matrix6 = Matrix<T, 6, 6>;

template<typename T>
using Vector3 = Matrix<T, 3, 1>;
template<typename T>
using Matrix3 = Matrix<T, 3, 3>;

typedef Vector6<double> Vector6d;
typedef Matrix6<double> Matrix6d;

} // namespace Eigen
