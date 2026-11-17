#include "fmmkernels.hpp"

namespace ngsfmm
{
  // Construct expansions here so their virtual methods are not emitted by every caller.
  template <typename T_MP, typename T_VALUE, typename T_Kappa>
  shared_ptr<BaseSingularMLExpansion> FMMInterface<T_MP,T_VALUE,T_Kappa> ::
  CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
  {
    return make_shared<SingularMLExpansion<T_MP,T_Kappa>> (c, r, kappa, fmm_params);
  }

  template <typename T_MP, typename T_VALUE, typename T_Kappa>
  shared_ptr<BaseRegularMLExpansion> FMMInterface<T_MP,T_VALUE,T_Kappa> ::
  CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
  {
    return make_shared<RegularMLExpansion<T_MP,T_Kappa>> (c, r, kappa, fmm_params);
  }

  // Export the supported factory combinations for callers outside this library.
  template class FMMInterface<Vec<1,Complex>,double,double>;
  template class FMMInterface<Vec<1,Complex>,Complex,double>;
  template class FMMInterface<Vec<1,Complex>,Complex,Complex>;
  template class FMMInterface<Vec<1,Complex32>,double,double>;
  template class FMMInterface<Vec<1,Complex32>,Complex,double>;
  template class FMMInterface<Vec<1,Complex32>,Complex,Complex>;
  template class FMMInterface<Vec<3,Complex>,double,double>;
  template class FMMInterface<Vec<3,Complex>,Complex,double>;
  template class FMMInterface<Vec<3,Complex>,Complex,Complex>;
  template class FMMInterface<Vec<3,Complex32>,double,double>;
  template class FMMInterface<Vec<3,Complex32>,Complex,double>;
  template class FMMInterface<Vec<3,Complex32>,Complex,Complex>;
  template class FMMInterface<Vec<6,Complex>,double,double>;
  template class FMMInterface<Vec<6,Complex32>,double,double>;
}
