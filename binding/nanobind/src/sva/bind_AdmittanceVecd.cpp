#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;

void bind_AdmittanceVecd(nb::module_ & sva)
{
  using AV = sva::AdmittanceVecd;
  using FV = sva::ForceVecd;
  using MV = sva::MotionVecd;
  using Vec3 = Eigen::Vector3d;
  using Vec6 = Eigen::Matrix<double, 6, 1>;

  auto av = nb::class_<AV>(sva, "AdmittanceVecd");
  av.def(nb::init(), "Default constructor, angular and linear are uninitialized")
      .def(nb::init<const Vec6 &>(), "Constructor from 6D vector", nb::arg("vec"))
      .def(nb::init<const Vec3 &, const Vec3 &>(), "Constructor from angular and linear vectors", nb::arg("angular"),
           nb::arg("linear"))
      .def(nb::init<const AV &>(), "Copy constructor")
      .def(nb::init<double, double>(), "Homogeneous constructor", nb::arg("angular"), nb::arg("linear"))
      .def_static("Zero", &AV::Zero, "Return a zero admittance vector")
      .def("angular", nb::overload_cast<>(&AV::angular, nb::const_), nb::rv_policy::reference_internal,
           "Get the angular part")
      .def("angular", nb::overload_cast<>(&AV::angular), nb::rv_policy::reference_internal,
           "Get the angular part (mutable)")
      .def("linear", nb::overload_cast<>(&AV::linear, nb::const_), nb::rv_policy::reference_internal,
           "Get the linear part")
      .def("linear", nb::overload_cast<>(&AV::linear), nb::rv_policy::reference_internal,
           "Get the linear part (mutable)")
      .def("vector", &AV::vector, "Get the 6D vector representation")
      .def(
          "cast", [](const AV & self) { return self.cast<double>(); },
          "Cast to another scalar type (only double supported in binding)")
      .def(
          "__add__", [](const AV & self, const AV & other) { return self + other; }, nb::arg("other"),
          "Add two AdmittanceVecd")
      .def(
          "__iadd__", [](AV & self, const AV & other) -> AV & { return self += other; }, nb::arg("other"),
          "In-place addition")
      .def(
          "__mul__", [](const AV & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__rmul__", [](const AV & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__imul__", [](AV & self, double scalar) -> AV & { return self *= scalar; }, nb::arg("scalar"),
          "In-place multiplication")
      .def(
          "__truediv__", [](const AV & self, double scalar) { return self / scalar; }, nb::arg("scalar"),
          "Divide by scalar")
      .def(
          "__itruediv__", [](AV & self, double scalar) -> AV & { return self /= scalar; }, nb::arg("scalar"),
          "In-place division")
      .def(
          "__eq__", [](const AV & self, const AV & other) { return self == other; }, nb::arg("other"), "Check equality")
      .def(
          "__ne__", [](const AV & self, const AV & other) { return self != other; }, nb::arg("other"),
          "Check inequality")
      .def(
          "__sub__", [](const AV & self, const AV & other)
          { return AV(self.angular() - other.angular(), self.linear() - other.linear()); }, nb::arg("other"),
          "Subtract two AdmittanceVecd")
      .def(
          "__isub__",
          [](AV & self, const AV & other) -> AV &
          {
            self.angular() -= other.angular();
            self.linear() -= other.linear();
            return self;
          },
          nb::arg("other"), "In-place subtraction")
      // Multiplication with ForceVecd (returns MotionVecd)
      .def(
          "__mul__", [](const AV & self, const FV & fv) { return self * fv; }, nb::arg("force_vec"),
          "Multiply with ForceVecd (returns MotionVecd)")
      .def(
          "__rmul__", [](const AV & self, const FV & fv) { return self * fv; }, nb::arg("force_vec"),
          "Multiply with ForceVecd (returns MotionVecd)")
      .def("__repr__",
           [](const AV & self)
           {
             std::ostringstream ss;
             ss << self;
             return ss.str();
           })
      // pickle (serialization)
      .def("__getstate__", [](const AV & self) { return self.vector(); })
      .def("__setstate__", [](AV & self, nb::handle h) { self = AV(nb::cast<Vec6>(h)); });
}
