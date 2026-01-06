#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;

void bind_ForceVecd(nb::module_ & sva)
{
  using PT = sva::ForceVecd;

  // sva.def("hello", []() {}, "hello function documentation");
  auto pt = nb::class_<sva::ForceVecd>(sva, "ForceVecd");
  pt.def(nb::init(), "Default constructor, rotation and translation are unitialized");
  pt.def("__repr__",
         [](PT & self)
         {
           std::ostringstream ss;
           ss << self;
           return ss.str();
         });
}
