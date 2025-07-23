#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_PTransformd(nb::module_ &);
void bind_MotionVecd(nb::module_ &);
void bind_ForceVecd(nb::module_ &);

NB_MODULE(_sva, unused_m)
{
  /**
   * We use a trick to load the symbols in the _sva namespace within the parent sva module's namespace.
   * This is done so that we expose the sva types as sva.Type instead of sva._sva.Type.
   *
   * In particular this makes types exported in stub files cleaner and more user-friendly, and makes sphinx autoapi
   * cross-references work accross projects.
   *
   * See https://github.com/wjakob/nanobind/issues/420#issuecomment-1950233151 for further details
   */
  (void)unused_m;
  nb::module_ m = nb::module_::import_("sva");
  // register functions and classes here
  m.doc() = "Python bindings for the SpaceVecAlg library";
  bind_PTransformd(m);
  bind_ForceVecd(m);
  bind_MotionVecd(m);
}
