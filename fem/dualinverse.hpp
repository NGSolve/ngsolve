#ifndef FILE_DUALINVERSE
#define FILE_DUALINVERSE

/*
  Inverse duality matrix D_ij = <dual_i, phi_j> of a vector-valued element,
  computed once on the reference element and shared through a slot owned by
  the space (one per vertex-order class)
*/

#include "scalarfe.hpp"
#include "elementtransformation.hpp"

namespace ngfem
{
  template <ELEMENT_TYPE ET, int DIMS, typename FEL>
  double * GetDualInverse (const FEL & fel, std::atomic<double*> & slot, LocalHeap & lh)
  {
    double * pinv = slot.load (std::memory_order_acquire);
    if (pinv) return pinv;

    constexpr int DIM = ET_trait<ET>::DIM;
    int ndof = fel.GetNDof();
    FE_ElementTransformation<DIM,DIM> trafo(ET);
    Matrix<> D(ndof, ndof), shape(ndof, DIMS), dshape(ndof, DIMS);
    D = 0.0;
    for (int vb = VOL; vb < DIM; vb++)
      {
        Facet2ElementTrafo f2el(ET, VorB(vb));
        for (int f : Range(f2el.GetNFacets()))
          {
            HeapReset hr(lh);
            IntegrationRule irf(f2el.FacetType(f), 2*fel.Order()+2);
            auto & mir = trafo(f2el(f, irf, lh), lh);
            for (int i : Range(mir.Size()))
              {
                fel.CalcMappedShape (mir[i], shape);
                fel.CalcDualShape (mir[i], dshape);
                D += mir[i].GetWeight() * dshape * Trans(shape);
              }
          }
      }
    CalcInverse (D);
    pinv = new double[ndof*ndof];
    FlatMatrix<>(ndof, ndof, pinv) = D;
    double * expected = nullptr;
    if (!slot.compare_exchange_strong (expected, pinv, std::memory_order_acq_rel))
      {
        delete [] pinv;
        pinv = expected;
      }
    return pinv;
  }
}

#endif
