/*

You may want to include this file to keep your application compatible with 
different versions of ngsolve.

*/




#ifdef FILE_POSTPROC

namespace ngcomp
{

  // // argument MeshAccess ma  was removed since it is available in GridFunction u
  // // change: ngsolve-5.1 - r1177, 20130404
  // NGS_DLL_HEADER 
  // inline void SetValues (const MeshAccess & ma,
  //                        const CoefficientFunction & coef,
  //                        GridFunction & u,
  //                        bool vb,
  //                        DifferentialOperator * diffop,   // NULL is FESpace evaluator
  //                        LocalHeap & clh)
  // {
  //   SetValues (coef, u, vb, diffop, clh);
  // }
}
  

#endif



#define L2HighOrderFiniteElement DGFiniteElement






