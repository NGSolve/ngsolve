#include <comp.hpp>
#include "compressedfespace.hpp"
#include "plateaufespace.hpp"

namespace ngcomp
{

  PlateauFESpace ::
  PlateauFESpace (shared_ptr<FESpace> afes, Array<Region> aplateaus)
    : CompressedFESpace(afes), plateaus(aplateaus)
  {
    /*
    cout << "PlateauFESpace, plateaus: " << endl;
    for (auto p : plateaus)
      cout << "p = " << p.Mask() << endl;
    */
  }

  
  void PlateauFESpace :: Update()
  {
    // space->Update();   // rely on autoupdate
    FESpace::Update();
    size_t ndofall = space->GetNDof();

    // cout << "ndof = " << ndofall << endl;
    
    all2comp.SetSize(ndofall);
    for (size_t i = 0; i < all2comp.Size(); i++)
      all2comp[i] = i;

    ctofdof.SetSize(ndofall);
    for (int i : Range(ndofall))
      ctofdof[i] = space->GetDofCouplingType(i);
    

    Array<DofId> plateau_refdof(plateaus.Size());

    size_t nv = ma->GetNV();
    // cout << "nv = " << nv << endl;
    for (auto ri : Range(plateaus))
      {
        BitArray regdofs = space->GetDofs(plateaus[ri]);

        DofId first = -1;
        for (size_t i = 0; i < nv; i++)
          if (regdofs[i])
            {
              first = i;
              break;
            }
        for (size_t i = 0; i < nv; i++)
          if (regdofs[i])
            {
              all2comp[i] = first;
              if (i != first)
                SetDofCouplingType(i, UNUSED_DOF);
            }
        for (size_t i = nv; i < ndofall; i++)
          if (regdofs[i])
            {
              all2comp[i]=-1;
              SetDofCouplingType(i, UNUSED_DOF);              
            }
      }
    SetNDof(ndofall);
  }
  
}
