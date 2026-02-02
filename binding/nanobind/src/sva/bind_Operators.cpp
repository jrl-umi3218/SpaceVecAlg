// filepath: binding/nanobind/src/sva/bind_operators_api.cpp
#include <SpaceVecAlg/Operators.h>
#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_Operators(nb::module_ & m)
{
  using Mat3X = Eigen::Matrix<double, 3, Eigen::Dynamic>;

  m.def(
      "colwiseCrossEq",
      [](const Mat3X & m1, const Mat3X & m2, Mat3X & result)
      {
        // result is modified in-place, shape (3, N)
        sva::sva_internal::colwiseCrossEq(m1, m2, result);
      },
      nb::arg("m1"), nb::arg("m2"), nb::arg("result"),
      "Column-wise cross product: result[:,i] = m1[:,i] x m2[:,i] for each column");
}
