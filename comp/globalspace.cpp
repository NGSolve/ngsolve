#include <comp.hpp>
#include "globalspace.hpp"

namespace ngcomp
{

  void GlobalSpace::VolDiffOp::CalcMatrix
    (const FiniteElement & bfel,
     const BaseMappedIntegrationPoint & mip,
     SliceMatrix<double,ColMajor> mat,
     LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<double> basisvec(basis->Dimension(), lh);
    basis -> Evaluate (mip, basisvec);

    for (int i = 0, ii = 0; i < dim; i++)
      for (int j = 0; j < vecdim; j++)
        mat(j,i) = basisvec(ii++);
  }


  GlobalSpace::GlobalSpace(shared_ptr<MeshAccess> ama,
                           const Flags& flags)
    : FESpace(ama, flags)
  {
    order = 5;
    basis = std::any_cast<shared_ptr<CoefficientFunction>>(flags.GetAnyFlag("basis"));

    auto dims = basis->Dimensions();
    if (dims.Size() == 1)
      {
        dim = dims[0];
        vecdim = 1;
      }
    else
      {
        dim = dims[1];
        vecdim = dims[0];
      }

    SetNDof(dim);
    evaluator[VOL] = make_shared<VolDiffOp>(basis, dim, vecdim);
  }


  FiniteElement & GlobalSpace :: GetFE (ElementId ei, Allocator & lh) const
  {
    return *new (lh) FE(dim, ma->GetElType(ei));
  }
  
  void GlobalSpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize0();
    for (int i = 0; i < dim; i++)
      dnums.Append(i);
  }
  
}
