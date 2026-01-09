#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;

void bind_MotionVecd(nb::module_ & sva)
{
  using MV = sva::MotionVecd;
  using FV = sva::ForceVecd;
  using Vec3 = Eigen::Vector3d;
  using Vec6 = Eigen::Matrix<double, 6, 1>;

  auto mv = nb::class_<MV>(sva, "MotionVecd");
  mv.def(nb::init(), "Default constructor, angular and linear are uninitialized")
      .def(nb::init<const Vec6 &>(), "Constructor from 6D vector", nb::arg("vec"))
      .def(nb::init<const Vec3 &, const Vec3 &>(), "Constructor from angular and linear vectors", nb::arg("angular"),
           nb::arg("linear"))
      .def_static("Zero", &MV::Zero, "Return a zero motion vector")
      .def("angular", nb::overload_cast<>(&MV::angular, nb::const_), nb::rv_policy::reference_internal,
           "Get the angular part")
      .def("angular", nb::overload_cast<>(&MV::angular), nb::rv_policy::reference_internal,
           "Get the angular part (mutable)")
      .def("linear", nb::overload_cast<>(&MV::linear, nb::const_), nb::rv_policy::reference_internal,
           "Get the linear part")
      .def("linear", nb::overload_cast<>(&MV::linear), nb::rv_policy::reference_internal,
           "Get the linear part (mutable)")
      .def("vector", &MV::vector, "Get the 6D vector representation")
      .def(
          "cast", [](const MV & self) { return self.cast<double>(); },
          "Cast to another scalar type (only double supported in binding)")
      .def("dot", &MV::dot, nb::arg("force_vec"), "Dot product with a ForceVecd")
      .def(
          "cross", [](const MV & self, const MV & other) { return self.cross(other); }, nb::arg("other"),
          "Cross product with another MotionVecd")
      .def(
          "crossDual", [](const MV & self, const FV & fv) { return self.crossDual(fv); }, nb::arg("force_vec"),
          "Cross dual with a ForceVecd")
      .def(
          "__add__", [](const MV & self, const MV & other) { return self + other; }, nb::arg("other"),
          "Add two MotionVecd")
      .def(
          "__sub__", [](const MV & self, const MV & other) { return self - other; }, nb::arg("other"),
          "Subtract two MotionVecd")
      .def(
          "__neg__", [](const MV & self) { return -self; }, "Unary minus")
      .def(
          "__mul__", [](const MV & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__rmul__", [](const MV & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__truediv__", [](const MV & self, double scalar) { return self / scalar; }, nb::arg("scalar"),
          "Divide by scalar")
      .def(
          "__iadd__", [](MV & self, const MV & other) -> MV & { return self += other; }, nb::arg("other"),
          "In-place addition")
      .def(
          "__isub__", [](MV & self, const MV & other) -> MV & { return self -= other; }, nb::arg("other"),
          "In-place subtraction")
      .def(
          "__imul__", [](MV & self, double scalar) -> MV & { return self *= scalar; }, nb::arg("scalar"),
          "In-place multiplication")
      .def(
          "__itruediv__", [](MV & self, double scalar) -> MV & { return self /= scalar; }, nb::arg("scalar"),
          "In-place division")
      .def(
          "__eq__", [](const MV & self, const MV & other) { return self == other; }, nb::arg("other"), "Check equality")
      .def(
          "__ne__", [](const MV & self, const MV & other) { return self != other; }, nb::arg("other"),
          "Check inequality")
      .def("__repr__",
           [](const MV & self)
           {
             std::ostringstream ss;
             ss << self;
             return ss.str();
           });
}
// ...existing code...
