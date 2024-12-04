#include <fespace.hpp>
#include <elementbyelement.hpp>

using namespace ngcomp;

namespace ngcomp
{
  shared_ptr<BaseMatrix> CreateMatrixFreeOperator (shared_ptr<FESpace> source_fes,
                                                   shared_ptr<FESpace> target_fes,
                                                   std::function<tuple<Matrix<>,Array<int>,Array<int>>(ElementId)> creator,
                                                   LocalHeap & lh)
  {
    // bool mixed = target_fes != nullptr;
    if (!target_fes) target_fes = source_fes;
    auto ma = source_fes->GetMeshAccess();


    const Table<size_t> & table = ma->GetElementsOfClass();

    shared_ptr<BaseMatrix> sum;
  
    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;
      
        ElementId ei(VOL,elclass_inds[0]);
        // auto & trafo = ma->GetTrafo(ei, lh);
        // auto & felx = source_fes->GetFE (ei, lh);
        // auto & fely = GetTestSpace()->GetFE (ei, lh);
        // MixedFiniteElement fel(felx, fely);

        auto [elmat,sind,tind] = creator(ei);
        
        Table<DofId> sdofs(elclass_inds.Size(), sind.Size());
        Table<DofId> tdofs(elclass_inds.Size(), tind.Size());

        Array<DofId> dofs;
        for (auto i : Range(elclass_inds.Size()))
          {
            source_fes->GetDofNrs(ElementId(VOL,elclass_inds[i]), dofs);
            for (int j = 0; j < sind.Size(); j++)
              sdofs[i][j] = dofs[sind[j]];
            for (int j = 0; j < tind.Size(); j++)
              tdofs[i][j] = dofs[tind[j]];
          }

        auto mat = make_shared<ConstantElementByElementMatrix<>>
          (target_fes->GetNDof(), source_fes->GetNDof(),
           elmat, std::move(tdofs), std::move(sdofs));
      
        if (sum)
          sum = make_shared<SumMatrix>(sum, mat);
        else
          sum = mat;
      }

    return sum;
  }
                                                 
}
