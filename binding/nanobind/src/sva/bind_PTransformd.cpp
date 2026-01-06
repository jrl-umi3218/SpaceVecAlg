#include <SpaceVecAlg/SpaceVecAlg>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;

void bind_PTransformd(nb::module_ & sva)
{
  using PT = sva::PTransformd;

  // sva.def("hello", []() {}, "hello function documentation");
  auto pt = nb::class_<sva::PTransformd>(sva, "PTransformd");
  pt.def(nb::init(), "Default constructor, rotation and translation are unitialized");
  pt.def_static("Identity", &PT::Identity, "Creates an identity transformation");
  pt.def("__repr__",
         [](PT & self)
         {
           std::ostringstream ss;
           ss << self;
           return ss.str();
         });
}
