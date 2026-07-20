#include <bla.hpp>

namespace ngsbem
{
  using namespace ngbla;

  NGS_DLL_HEADER double LaplaceSL_Polygon (FlatArray<Vec<3>> vertices, Vec<3> x);
  NGS_DLL_HEADER Vec<3> LaplaceGradSL_Polygon (FlatArray<Vec<3>> vertices, Vec<3> x);
  NGS_DLL_HEADER double LaplaceDL_Polygon (FlatArray<Vec<3>> vertices, Vec<3> x,
                                           Vec<3> source_normal);
}
