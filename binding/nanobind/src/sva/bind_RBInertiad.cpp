#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;

void bind_RBInertiad(nb::module_ & sva)
{
  using RBI = sva::RBInertiad;
  using Vec3 = Eigen::Vector3d;
  using Mat3 = Eigen::Matrix3d;
  using Mat6 = Eigen::Matrix<double, 6, 6>;
  using MV = sva::MotionVecd;
  using FV = sva::ForceVecd;

  auto rbi = nb::class_<RBI>(sva, "RBInertiad");
  rbi.def(nb::init(), "Default constructor, mass, momentum, and inertia are uninitialized")
      .def(nb::init<double, const Vec3 &, const Mat3 &>(),
           "Constructor from mass, spatial momentum, and inertia matrix", nb::arg("mass"), nb::arg("momentum"),
           nb::arg("inertia"))
      .def("mass", &RBI::mass, "Get the mass")
      .def("momentum", &RBI::momentum, nb::rv_policy::reference_internal, "Get the spatial momentum (const)")
      .def("lowerTriangularInertia", &RBI::lowerTriangularInertia, nb::rv_policy::reference_internal,
           "Get the lower triangular inertia matrix (const)")
      .def("inertia", &RBI::inertia, "Get the full inertia matrix")
      .def("matrix", &RBI::matrix, "Get the 6x6 rigid body inertia matrix")
      .def(
          "cast", [](const RBI & self) { return self.cast<double>(); },
          "Cast to another scalar type (only double supported in binding)")
      .def(
          "__add__", [](const RBI & self, const RBI & other) { return self + other; }, nb::arg("other"),
          "Add two RBInertiad")
      .def(
          "__sub__", [](const RBI & self, const RBI & other) { return self - other; }, nb::arg("other"),
          "Subtract two RBInertiad")
      .def(
          "__neg__", [](const RBI & self) { return -self; }, "Unary minus")
      .def(
          "__iadd__", [](RBI & self, const RBI & other) -> RBI & { return self += other; }, nb::arg("other"),
          "In-place addition")
      .def(
          "__isub__", [](RBI & self, const RBI & other) -> RBI & { return self -= other; }, nb::arg("other"),
          "In-place subtraction")
      .def(
          "__mul__", [](const RBI & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__rmul__", [](const RBI & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__eq__", [](const RBI & self, const RBI & other) { return self == other; }, nb::arg("other"),
          "Check equality")
      .def(
          "__ne__", [](const RBI & self, const RBI & other) { return self != other; }, nb::arg("other"),
          "Check inequality")
      .def("__repr__",
           [](const RBI & self)
           {
             std::ostringstream ss;
             ss << self;
             return ss.str();
           });
  // Optionally, bind RBI * MotionVecd if needed:
  // .def("__mul__", [](const RBI & self, const MV & mv) { return self * mv; }, nb::arg("motion_vec"), "Multiply
  // RBInertiad by MotionVecd")
}
