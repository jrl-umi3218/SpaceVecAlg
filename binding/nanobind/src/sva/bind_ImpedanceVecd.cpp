#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;

void bind_ImpedanceVecd(nb::module_ & sva)
{
  using IV = sva::ImpedanceVecd;
  using MV = sva::MotionVecd;
  using FV = sva::ForceVecd;
  using Vec3 = Eigen::Vector3d;
  using Vec6 = Eigen::Matrix<double, 6, 1>;

  auto iv = nb::class_<IV>(sva, "ImpedanceVecd");
  iv.def(nb::init(), "Default constructor, angular and linear are uninitialized")
      .def(nb::init<const Vec6 &>(), "Constructor from 6D vector", nb::arg("vec"))
      .def(nb::init<const Vec3 &, const Vec3 &>(), "Constructor from angular and linear vectors", nb::arg("angular"),
           nb::arg("linear"))
      .def(nb::init<const IV &>(), "Copy constructor")
      .def(nb::init<double, double>(), "Homogeneous constructor", nb::arg("angular"), nb::arg("linear"))
      .def_static("Zero", &IV::Zero, "Return a zero impedance vector")
      .def("angular", nb::overload_cast<>(&IV::angular, nb::const_), nb::rv_policy::reference_internal,
           "Get the angular part")
      .def("angular", nb::overload_cast<>(&IV::angular), nb::rv_policy::reference_internal,
           "Get the angular part (mutable)")
      .def("linear", nb::overload_cast<>(&IV::linear, nb::const_), nb::rv_policy::reference_internal,
           "Get the linear part")
      .def("linear", nb::overload_cast<>(&IV::linear), nb::rv_policy::reference_internal,
           "Get the linear part (mutable)")
      .def("vector", &IV::vector, "Get the 6D vector representation")
      .def(
          "cast", [](const IV & self) { return self.cast<double>(); },
          "Cast to another scalar type (only double supported in binding)")
      .def(
          "__add__", [](const IV & self, const IV & other) { return self + other; }, nb::arg("other"),
          "Add two ImpedanceVecd")
      .def(
          "__iadd__", [](IV & self, const IV & other) -> IV & { return self += other; }, nb::arg("other"),
          "In-place addition")
      .def(
          "__mul__", [](const IV & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__rmul__", [](const IV & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__imul__", [](IV & self, double scalar) -> IV & { return self *= scalar; }, nb::arg("scalar"),
          "In-place multiplication")
      .def(
          "__truediv__", [](const IV & self, double scalar) { return self / scalar; }, nb::arg("scalar"),
          "Divide by scalar")
      .def(
          "__itruediv__", [](IV & self, double scalar) -> IV & { return self /= scalar; }, nb::arg("scalar"),
          "In-place division")
      .def(
          "__eq__", [](const IV & self, const IV & other) { return self == other; }, nb::arg("other"), "Check equality")
      .def(
          "__ne__", [](const IV & self, const IV & other) { return self != other; }, nb::arg("other"),
          "Check inequality")
      .def(
          "__sub__", [](const IV & self, const IV & other)
          { return IV(self.angular() - other.angular(), self.linear() - other.linear()); }, nb::arg("other"),
          "Subtract two ImpedanceVecd")
      .def(
          "__isub__",
          [](IV & self, const IV & other) -> IV &
          {
            self.angular() -= other.angular();
            self.linear() -= other.linear();
            return self;
          },
          nb::arg("other"), "In-place subtraction")
      // Multiplication with MotionVecd (returns ForceVecd)
      .def(
          "__mul__", [](const IV & self, const MV & mv) { return self * mv; }, nb::arg("motion_vec"),
          "Multiply with MotionVecd (returns ForceVecd)")
      .def(
          "__rmul__", [](const MV & mv, const IV & self) { return self * mv; }, nb::arg("motion_vec"),
          "Multiply with MotionVecd (returns ForceVecd)")
      .def("__repr__",
           [](const IV & self)
           {
             std::ostringstream ss;
             ss << self;
             return ss.str();
           })
      // pickle (serialization)
      .def("__getstate__", [](const IV & self) { return self.vector(); })
      .def("__setstate__", [](IV & self, nb::handle h) { self = IV(nb::cast<Vec6>(h)); });
}
