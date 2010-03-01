#ifndef FILE_POSTPROC
#define FILE_POSTPROC

/*********************************************************************/
/* File:   postproc.hh                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngcomp
{

  /* 
     Postprocessing functions
  */

  template <class SCAL>
  extern NGS_DLL_HEADER void CalcFlux (const MeshAccess & ma, 
			const S_GridFunction<SCAL> & u,
			S_GridFunction<SCAL> & flux,
			const BilinearFormIntegrator & bli,
			bool applyd, bool add,
			int domain);

  template <class SCAL>
  extern NGS_DLL_HEADER void CalcFluxProject (const MeshAccess & ma, 
			       const S_GridFunction<SCAL> & u,
			       S_GridFunction<SCAL> & flux,
			       const BilinearFormIntegrator & bli,
			       bool applyd, int domain,
			       LocalHeap & lh);

  template <class SCAL>
  extern NGS_DLL_HEADER void CalcFluxProject (const MeshAccess & ma, 
			       const S_GridFunction<SCAL> & u,
			       S_GridFunction<SCAL> & flux,
			       const BilinearFormIntegrator & bli,
			       bool applyd, const BitArray & domains, LocalHeap & lh);


  template <class SCAL>
  extern NGS_DLL_HEADER void SetValues (const MeshAccess & ma, 
			 const CoefficientFunction & coef,
			 GridFunction & u,
			 bool bound,
			 DifferentialOperator * diffop,   // NULL is FESpace evaluator
			 LocalHeap & clh);
    


  template <class SCAL>
  extern NGS_DLL_HEADER int CalcPointFlux (const MeshAccess & ma, 
			    const GridFunction & u,
			    const FlatVector<double> & point,
			    FlatVector<SCAL> & flux,
			    const BilinearFormIntegrator & bli,
			    bool applyd,
			    LocalHeap & lh,
			    int component = 0);

  template <class SCAL>
  extern NGS_DLL_HEADER int CalcPointFlux (const MeshAccess & ma, 
			    const GridFunction & u,
			    const FlatVector<double> & point,
			    const Array<int> & domains,
			    FlatVector<SCAL> & flux,
			    const BilinearFormIntegrator & bli,
			    bool applyd,
			    LocalHeap & lh,
			    int component = 0);

  template <class SCAL>
  extern NGS_DLL_HEADER void CalcError (const MeshAccess & ma, 
			 const S_GridFunction<SCAL> & bu,
			 const S_GridFunction<SCAL> & bflux,
			 const BilinearFormIntegrator & bli,
			 FlatVector<double> & err,
			 int domain,
			 LocalHeap & lh);
  template <class SCAL>
  extern NGS_DLL_HEADER void CalcError (const MeshAccess & ma, 
			 const S_GridFunction<SCAL> & u,
			 const S_GridFunction<SCAL> & flux,
			 const BilinearFormIntegrator & bli,
			 FlatVector<double> & err,
			 const BitArray & domains, LocalHeap & lh);


  template <class SCAL>
  NGS_DLL_HEADER void CalcDifference (const MeshAccess & ma, 
		       const S_GridFunction<SCAL> & u1,
		       const S_GridFunction<SCAL> & u2,
		       const BilinearFormIntegrator & bli1,
		       const BilinearFormIntegrator & bli2,
		       FlatVector<double> & diff,
		       int domain, LocalHeap & lh);
  
  template <class SCAL>
  NGS_DLL_HEADER void CalcDifference (const MeshAccess & ma, 
		       const S_GridFunction<SCAL> & u1,
		       const BilinearFormIntegrator & bli1,
		       const CoefficientFunction * coef_real, 
		       const CoefficientFunction * coef_imag,
		       FlatVector<double> & diff,
		       int domain, LocalHeap & lh);



  template <class SCAL>
  extern NGS_DLL_HEADER void CalcGradient (const MeshAccess & ma,
			    const FESpace & fesh1,
			    const S_BaseVector<SCAL> & vech1,
			    const FESpace & feshcurl,
			    S_BaseVector<SCAL> & vechcurl);

  template <class SCAL>
  extern NGS_DLL_HEADER void CalcGradientT (const MeshAccess & ma,
			     const FESpace & feshcurl,
			     const S_BaseVector<SCAL> & vechcurl,
			     const FESpace & fesh1,
			     S_BaseVector<SCAL> & vech1);


  template <class SCAL>
  extern NGS_DLL_HEADER void CalcErrorHierarchical (const MeshAccess & ma, 
				     const S_BilinearForm<SCAL> & bfa,
				     const S_BilinearForm<SCAL> & bfa2,
				     const S_LinearForm<SCAL> & lff,
				     S_GridFunction<SCAL> & gfu,
				     const FESpace & festest,
				     FlatVector<double> & err,
				     LocalHeap & lh);


  /*
 /// calculate elementwise flux
 extern void CalcFlux (const MeshAccess & ma, const FESpace & fes,
 const BDBIntegrator<> & bli,
 const BaseVector & u, BaseVector & flux, int applyd = 1);

 ///


 extern void CalcFluxNodal (const MeshAccess & ma, 
 const FESpace & fespaceu,
 const FESpace & fespaceflux,
 const BDBIntegrator<> & bli,
 const BaseVector & u, BaseVector & flux, int applyd = 1,
 int dom = 0);

 ///


 extern void CalcFluxElement (const MeshAccess & ma, 
 const FESpace & fespaceu,
 const ElementFESpace & fespaceflux,
 const BDBIntegrator<> & bli,
 const BaseVector & u, BaseVector & flux, int applyd = 1,
 int dom = 0);


 ///


 extern void CalcBoundaryFlux (const MeshAccess & ma, 
 const FESpace & fespaceu,
 const BDBBoundaryIntegrator<> & bli,
 const BaseVector & u, BaseVector & flux, int applyd = 1);

 ///


 extern void CalcBoundaryFluxNodal (const MeshAccess & ma, 
 const FESpace & fespaceu,
 const FESpace & fespaceflux,
 const BDBBoundaryIntegrator<> & bli,
 const BaseVector & u, BaseVector & flux, int applyd = 1,
 int dom = 0);


 ///


 extern void ZZErrorEstimator2 (const MeshAccess & ma, 
 const FESpace & fespaceu,
 const FESpace & fespaceflux,
 const BDBIntegrator<> & bli,
 const BaseVector & u, const BaseVector & flux, 
 Vector & elerr, int dom = 0);

 ///


 extern void ZZBoundaryErrorEstimator2 (const MeshAccess & ma, 
 const FESpace & fespaceu,
 const FESpace & fespaceflux,
 const BDBBoundaryIntegrator<> & bli,
 const BaseVector & u, const BaseVector & flux, 
 Vector & elerr, int dom = 0);




 ///

 extern void AverageElementData (const MeshAccess & ma,
 const BaseSystemVector & eldata, 
 BaseSystemVector & nodaldata,
 int dom = 0);



 ///

 extern void ZZErrorEstimate (const MeshAccess & ma, 
 const BaseSystemVector & eldata, 
 const BaseSystemVector & nodaldata,
 Vector & elerr,
 int dom = 0);

 ///

 extern void Interpolate2Nodal (const MeshAccess & ma,
 const FESpace & fespace,
 const FESpace & fespacenodal,
 const BaseSystemVector & data, 
 BaseSystemVector & nodaldata,
 int dom = 0);

 ///

 extern void ZZErrorEstimate (const MeshAccess & ma, 
 const FESpace & fespace,
 const FESpace & fespacenodal,
 const BaseSystemVector & eldata, 
 const BaseSystemVector & nodaldata,
 Vector & elerr,
 int dom = 0);

 ///

 extern int PointEvaluation (const MeshAccess & ma,
 const FESpace & fespace,
 const BaseVector & vec,
 BDBIntegrator<> & bli,
 double * point,
 Vector & result,
 int applyd = 0);
  */

}

#endif

