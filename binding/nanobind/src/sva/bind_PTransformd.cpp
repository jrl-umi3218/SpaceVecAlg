#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanoeigenpy/geometry/quaternion.hpp>
#include <sstream>

namespace nb = nanobind;

void bind_PTransformd(nb::module_ & sva)
{
  nanoeigenpy::exposeQuaternion<double>(sva, "Quaternion");

  using PT = sva::PTransformd;
  using Vec3 = Eigen::Vector3d;
  using Mat3 = Eigen::Matrix3d;

  auto pt = nb::class_<PT>(sva, "PTransformd");
  pt.def(nb::init(), "Default constructor, rotation and translation are unitialized")
      .def(nb::init<const Mat3 &, const Vec3 &>(), "Constructor from rotation matrix and translation vector",
           nb::arg("rotation"), nb::arg("translation"))
      .def(nb::init<const Eigen::Quaterniond &, const Vec3 &>(), "Constructor from quaternion and translation vector",
           nb::arg("quaternion"), nb::arg("translation"))
      .def(nb::init<const Mat3 &>(), "Constructor from rotation matrix", nb::arg("rotation"))
      .def(nb::init<const Eigen::Quaterniond &>(), "Constructor from quaternion", nb::arg("quaternion"))
      .def(nb::init<const Vec3 &>(), "Constructor from translation vector", nb::arg("translation"))
      .def_static("Identity", &PT::Identity, "Creates an identity transformation")
      .def("rotation", nb::overload_cast<>(&PT::rotation, nb::const_), nb::rv_policy::reference_internal,
           "Get the rotation matrix")
      .def("translation", nb::overload_cast<>(&PT::translation, nb::const_), nb::rv_policy::reference_internal,
           "Get the translation vector")
      .def("inv", &PT::inv, "Return the inverse transformation")
      .def("matrix", &PT::matrix, "Return the 6x6 Plücker transformation matrix")
      .def("__repr__",
           [](const PT & self)
           {
             std::ostringstream ss;
             ss << self;
             return ss.str();
           });

  pt.def("dualMatrix", &PT::dualMatrix, "Return the 6x6 dual Plücker transformation matrix");

  // Bind templated free functions in the sva namespace
  sva.def("RotX", &sva::RotX<double>, nb::arg("theta"), "Create a rotation matrix about the X axis");
  sva.def("RotY", &sva::RotY<double>, nb::arg("theta"), "Create a rotation matrix about the Y axis");
  sva.def("RotZ", &sva::RotZ<double>, nb::arg("theta"), "Create a rotation matrix about the Z axis");
  sva.def("rotationError", &sva::rotationError<double>, nb::arg("E_a_b"), nb::arg("E_a_c"),
          "Compute the 3D rotation error between two matrices");
  sva.def("rotationVelocity", &sva::rotationVelocity<double>, nb::arg("E_a_b"),
          "Compute the 3D rotation vector of the rotation matrix");
  sva.def("transformError", &sva::transformError<double>, nb::arg("X_a_b"), nb::arg("X_a_c"),
          "Compute the 6D error between two PTransformd objects");
  sva.def("transformVelocity", &sva::transformVelocity<double>, nb::arg("X_a_b"),
          "Compute the motion vector of the matrix X_a_b");
  sva.def("interpolate", &sva::interpolate<double>, nb::arg("from"), nb::arg("to"), nb::arg("t"),
          "Interpolate between two PTransformd objects");
  sva.def("sinc_inv", &sva::sinc_inv<double>, nb::arg("x"), "numerically stable inverse sinc function");

  // bind operators
  pt.def(
      "__mul__", [](const PT & self, const PT & other) { return self * other; }, nb::arg("other"),
      "Multiply two PTransformd objects");

  pt.def(
      "__eq__", [](const PT & self, const PT & other) { return self == other; }, nb::arg("other"),
      "Check equality of two PTransformd objects");

  pt.def(
      "__ne__", [](const PT & self, const PT & other) { return self != other; }, nb::arg("other"),
      "Check inequality of two PTransformd objects");
}
// ...existing code...
