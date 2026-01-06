#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_PTransformd(nb::module_ &);
void bind_MotionVecd(nb::module_ &);
void bind_ForceVecd(nb::module_ &);

NB_MODULE(sva_pywrap, m)
{
  // register functions and classes here
  m.doc() = "Python bindings for the SpaceVecAlg library";
  bind_PTransformd(m);
  bind_ForceVecd(m);
  bind_MotionVecd(m);
}
