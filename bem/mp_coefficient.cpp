#include "mp_coefficient.hpp"

namespace ngsbem
{
  template <typename entry_type>
  shared_ptr<RegularMLExpansionCF<entry_type>>
  SingularMLExpansionCF<entry_type> :: CreateRegularExpansion(Vec<3> center, double r) const
  {
    mlmp->CalcMP();
    return make_shared<RegularMLExpansionCF<entry_type>> (this->MLExpansion(), center, r);
  }

  template class SingularMLExpansionCF<Complex>;
  template class SingularMLExpansionCF<Vec<3,Complex>>;
}
