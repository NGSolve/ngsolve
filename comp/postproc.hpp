#ifndef FILE_POSTPROC
#define FILE_POSTPROC

/*********************************************************************/
/* File:   postproc.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngcomp
{

  /* 
     Postprocessing functions
  */

  template <class SCAL>
  extern NGS_DLL_HEADER void CalcFlux (const S_GridFunction<SCAL> & u,
                                       S_GridFunction<SCAL> & flux,
                                       shared_ptr<BilinearFormIntegrator> bli,
                                       bool applyd, bool add,
                                       int domain);
  
  extern NGS_DLL_HEADER 
  void CalcFluxProject (const GridFunction & u,
			GridFunction & flux,
			shared_ptr<BilinearFormIntegrator> bli,
			bool applyd, int domain,
			LocalHeap & lh);

  template <class SCAL>
  extern NGS_DLL_HEADER void CalcFluxProject (const S_GridFunction<SCAL> & u,
                                              S_GridFunction<SCAL> & flux,
                                              shared_ptr<BilinearFormIntegrator> bli,
                                              bool applyd, const BitArray & domains, LocalHeap & lh);


  extern NGS_DLL_HEADER 
  void SetValues (shared_ptr<CoefficientFunction> coef,
		  GridFunction & u,
		  VorB vb,
		  DifferentialOperator * diffop,   // NULL is FESpace evaluator
		  LocalHeap & clh,
                  bool dualdiffop = false, bool use_simd = true, int mdcomp=0,
                  optional<shared_ptr<BitArray>> definedonelements = nullopt,
                  int bonus_intorder=0);
  
  extern NGS_DLL_HEADER 
  void SetValues (shared_ptr<CoefficientFunction> coef,
		  GridFunction & u,
		  const Region & region, 
		  DifferentialOperator * diffop,   // NULL is FESpace evaluator
		  LocalHeap & clh,
                  bool dualdiffop = false, bool use_simd = true, int mdcomp=0,
                  optional<shared_ptr<BitArray>> definedonelements = nullopt,
                  int bonus_intorder=0);
  

  template <class SCAL>
  extern NGS_DLL_HEADER
  int CalcPointFlux (const GridFunction & u,
		     const FlatVector<double> & point,
		     FlatVector<SCAL> & flux,
		     shared_ptr<BilinearFormIntegrator> bli,
		     bool applyd,
		     LocalHeap & lh,
		     int component = 0);

  
  template <class SCAL>
  extern NGS_DLL_HEADER 
  int CalcPointFlux (const GridFunction & u,
		     const FlatVector<double> & point,
		     const Array<int> & domains,
		     FlatVector<SCAL> & flux,
		     shared_ptr<BilinearFormIntegrator> bli,
		     bool applyd,
		     LocalHeap & lh,
		     int component = 0);


  extern NGS_DLL_HEADER 
  void CalcError (const GridFunction & bu,
		  const GridFunction & bflux,
		  shared_ptr<BilinearFormIntegrator> bli,
		  FlatVector<double> & err,
		  int domain,
		  LocalHeap & lh);

  template <class SCAL>
  extern NGS_DLL_HEADER void CalcError (const S_GridFunction<SCAL> & u,
                                        const S_GridFunction<SCAL> & flux,
                                        shared_ptr<BilinearFormIntegrator> bli,
                                        FlatVector<double> & err,
                                        const BitArray & domains, LocalHeap & lh);


  template <class SCAL>
  NGS_DLL_HEADER void CalcDifference (const S_GridFunction<SCAL> & u1,
                                      const S_GridFunction<SCAL> & u2,
                                      shared_ptr<BilinearFormIntegrator> bli1,
                                      shared_ptr<BilinearFormIntegrator> bli2,
                                      FlatVector<double> & diff,
                                      int domain, LocalHeap & lh);
  
  NGS_DLL_HEADER void CalcDifference (const GridFunction & u1,
                                      shared_ptr<BilinearFormIntegrator> bli1,
                                      shared_ptr<CoefficientFunction> coef, 
                                      FlatVector<double> & diff,
                                      int domain, LocalHeap & lh);



  template <class SCAL>
  extern NGS_DLL_HEADER void CalcGradient (shared_ptr<MeshAccess> ma,
			    const FESpace & fesh1,
			    const S_BaseVector<SCAL> & vech1,
			    const FESpace & feshcurl,
			    S_BaseVector<SCAL> & vechcurl);

  template <class SCAL>
  extern NGS_DLL_HEADER void CalcGradientT (shared_ptr<MeshAccess> ma,
			     const FESpace & feshcurl,
			     const S_BaseVector<SCAL> & vechcurl,
			     const FESpace & fesh1,
			     S_BaseVector<SCAL> & vech1);


  template <class SCAL>
  extern NGS_DLL_HEADER 
  void CalcErrorHierarchical (const S_BilinearForm<SCAL> & bfa,
                              const S_BilinearForm<SCAL> & bfa2,
                              const S_LinearForm<SCAL> & lff,
                              S_GridFunction<SCAL> & gfu,
                              const FESpace & festest,
                              FlatVector<double> & err,
                              LocalHeap & lh);

}

#endif

