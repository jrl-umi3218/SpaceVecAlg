#include <SpaceVecAlg/EigenUtility.h>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_EigenUtility(nb::module_ & m)
{
  using Vec3 = Eigen::Vector3d;
  using Vec6 = Eigen::Matrix<double, 6, 1>;
  using Mat3 = Eigen::Matrix3d;
  using Mat6 = Eigen::Matrix<double, 6, 6>;

  m.def(
      "vector3ToCrossMatrix", [](const Vec3 & v) { return sva::vector3ToCrossMatrix(v); }, nb::arg("vec"),
      "Convert a 3D vector to a cross product matrix (3x3)");

  m.def(
      "vector6ToCrossMatrix", [](const Vec6 & v) { return sva::vector6ToCrossMatrix(v); }, nb::arg("vec"),
      "Convert a 6D vector to a spatial cross product matrix (6x6)");

  m.def(
      "vector6ToCrossDualMatrix", [](const Vec6 & v) { return sva::vector6ToCrossDualMatrix(v); }, nb::arg("vec"),
      "Convert a 6D vector to a spatial dual cross product matrix (6x6)");
}
