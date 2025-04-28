#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;

void bind_MotionVecd(nb::module_ & sva)
{
  using PT = sva::MotionVecd;

  // sva.def("hello", []() {}, "hello function documentation");
  auto pt = nb::class_<sva::MotionVecd>(sva, "MotionVecd");
  pt.def(nb::init(), "Default constructor, rotation and translation are unitialized");
  pt.def("__repr__",
         [](PT & self)
         {
           std::ostringstream ss;
           ss << self;
           return ss.str();
         });
}
