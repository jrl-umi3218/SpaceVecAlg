#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;

void bind_ABInertiad(nb::module_ & sva)
{
  using ABI = sva::ABInertiad;
  using Vec3 = Eigen::Vector3d;
  using Mat3 = Eigen::Matrix3d;
  using Mat6 = Eigen::Matrix<double, 6, 6>;
  using RBInertia = sva::RBInertiad;
  using MV = sva::MotionVecd;
  using FV = sva::ForceVecd;

  auto abi = nb::class_<ABI>(sva, "ABInertiad");
  abi.def(nb::init(), "Default constructor, matrices are uninitialized")
      .def(nb::init<const Mat3 &, const Mat3 &, const Mat3 &>(),
           "Constructor from mass, generalized inertia, and inertia matrices", nb::arg("mass"), nb::arg("g_inertia"),
           nb::arg("inertia"))
      .def(nb::init<const ABI &>(), "Copy constructor")
      .def("lowerTriangularMassMatrix", &ABI::lowerTriangularMassMatrix, nb::rv_policy::reference_internal,
           "Get the lower triangular mass matrix (const)")
      .def("massMatrix", &ABI::massMatrix, "Get the full mass matrix")
      .def("gInertia", &ABI::gInertia, nb::rv_policy::reference_internal, "Get the generalized inertia matrix (const)")
      .def("lowerTriangularInertia", &ABI::lowerTriangularInertia, nb::rv_policy::reference_internal,
           "Get the lower triangular inertia matrix (const)")
      .def("inertia", &ABI::inertia, "Get the full inertia matrix")
      .def("matrix", &ABI::matrix, "Get the 6x6 articulated body inertia matrix")
      .def(
          "cast", [](const ABI & self) { return self.cast<double>(); },
          "Cast to another scalar type (only double supported in binding)")
      .def(
          "__add__", [](const ABI & self, const ABI & other) { return self + other; }, nb::arg("other"),
          "Add two ABInertiad")
      .def(
          "__sub__", [](const ABI & self, const ABI & other) { return self - other; }, nb::arg("other"),
          "Subtract two ABInertiad")
      .def(
          "__neg__", [](const ABI & self) { return -self; }, "Unary minus")
      .def(
          "__iadd__", [](ABI & self, const ABI & other) -> ABI & { return self += other; }, nb::arg("other"),
          "In-place addition")
      .def(
          "__isub__", [](ABI & self, const ABI & other) -> ABI & { return self -= other; }, nb::arg("other"),
          "In-place subtraction")
      .def(
          "__mul__", [](const ABI & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__rmul__", [](const ABI & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__eq__", [](const ABI & self, const ABI & other) { return self == other; }, nb::arg("other"),
          "Check equality")
      .def(
          "__ne__", [](const ABI & self, const ABI & other) { return self != other; }, nb::arg("other"),
          "Check inequality")
      .def("__repr__",
           [](const ABI & self)
           {
             std::ostringstream ss;
             ss << self;
             return ss.str();
           })
      // Optionally, bind ABI + RBInertiad and ABI * MotionVecd if needed:
      .def(
          "__add__", [](const ABI & self, const RBInertia & rbI) { return self + rbI; }, nb::arg("rb_inertia"),
          "Add ABInertiad and RBInertiad")
      .def(
          "__mul__", [](const ABI & self, const MV & mv) { return self * mv; }, nb::arg("motion_vec"),
          "Multiply ABInertiad by MotionVecd");
}
// ...existing code...
