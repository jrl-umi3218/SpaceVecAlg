#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;

void bind_ForceVecd(nb::module_ & sva)
{
  using FV = sva::ForceVecd;
  using Vec3 = Eigen::Vector3d;
  using Vec6 = Eigen::Matrix<double, 6, 1>;

  auto fv = nb::class_<FV>(sva, "ForceVecd");
  fv.def(nb::init(), "Default constructor, couple and force are uninitialized")
      .def(nb::init<const Vec6 &>(), "Constructor from 6D vector", nb::arg("vec"))
      .def(nb::init<const Vec3 &, const Vec3 &>(), "Constructor from couple and force vectors", nb::arg("couple"),
           nb::arg("force"))
      .def_static("Zero", &FV::Zero, "Return a zero force vector")
      .def("couple", nb::overload_cast<>(&FV::couple, nb::const_), nb::rv_policy::reference_internal,
           "Get the couple (const)")
      .def("couple", nb::overload_cast<>(&FV::couple), nb::rv_policy::reference_internal, "Get the couple (mutable)")
      .def("moment", nb::overload_cast<>(&FV::moment, nb::const_), nb::rv_policy::reference_internal,
           "Get the moment (const, same as couple)")
      .def("moment", nb::overload_cast<>(&FV::moment), nb::rv_policy::reference_internal,
           "Get the moment (mutable, same as couple)")
      .def("force", nb::overload_cast<>(&FV::force, nb::const_), nb::rv_policy::reference_internal,
           "Get the force (const)")
      .def("force", nb::overload_cast<>(&FV::force), nb::rv_policy::reference_internal, "Get the force (mutable)")
      .def("vector", &FV::vector, "Get the 6D vector representation")
      .def(
          "cast", [](const FV & self) { return self.cast<double>(); },
          "Cast to another scalar type (only double supported in binding)")
      .def(
          "__add__", [](const FV & self, const FV & other) { return self + other; }, nb::arg("other"),
          "Add two ForceVecd")
      .def(
          "__sub__", [](const FV & self, const FV & other) { return self - other; }, nb::arg("other"),
          "Subtract two ForceVecd")
      .def(
          "__neg__", [](const FV & self) { return -self; }, "Unary minus")
      .def(
          "__mul__", [](const FV & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__rmul__", [](const FV & self, double scalar) { return self * scalar; }, nb::arg("scalar"),
          "Multiply by scalar")
      .def(
          "__truediv__", [](const FV & self, double scalar) { return self / scalar; }, nb::arg("scalar"),
          "Divide by scalar")
      .def(
          "__iadd__", [](FV & self, const FV & other) -> FV & { return self += other; }, nb::arg("other"),
          "In-place addition")
      .def(
          "__isub__", [](FV & self, const FV & other) -> FV & { return self -= other; }, nb::arg("other"),
          "In-place subtraction")
      .def(
          "__imul__", [](FV & self, double scalar) -> FV & { return self *= scalar; }, nb::arg("scalar"),
          "In-place multiplication")
      .def(
          "__itruediv__", [](FV & self, double scalar) -> FV & { return self /= scalar; }, nb::arg("scalar"),
          "In-place division")
      .def(
          "__eq__", [](const FV & self, const FV & other) { return self == other; }, nb::arg("other"), "Check equality")
      .def(
          "__ne__", [](const FV & self, const FV & other) { return self != other; }, nb::arg("other"),
          "Check inequality")
      .def("__repr__",
           [](const FV & self)
           {
             std::ostringstream ss;
             ss << self;
             return ss.str();
           })
      // pickle (serialization)
      .def("__getstate__",
           [](const FV & self)
           {
             // Serialize as the internal 6D vector directly
             return self.vector();
           })
      .def("__setstate__",
           [](FV & self, nb::handle h)
           {
             // Accept a single vector directly
             self = FV(nb::cast<Vec6>(h));
           });
}
