/*********************************************************************/
/* File:   postproc.cpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   Postprocessing functions
*/

#include <comp.hpp>

#include <solve.hpp>
#include <fem.hpp>
// extern PDE * pde;
extern ngstd::AutoPtr<ngsolve::PDE> pde;


namespace ngcomp
{
  using namespace ngcomp;



  template <class SCAL>
  void CalcFluxProject (const MeshAccess & ma, 
			const S_GridFunction<SCAL> & u,
			S_GridFunction<SCAL> & flux,
			const BilinearFormIntegrator & bli,
			bool applyd, const BitArray & domains, LocalHeap & clh)
  {
    static int timer = NgProfiler::CreateTimer ("CalcFluxProject");
    NgProfiler::RegionTimer reg (timer);
    
    ma.PushStatus ("Post-processing");

    const FESpace & fes = u.GetFESpace();
    const FESpace & fesflux = flux.GetFESpace();

    bool bound = bli.BoundaryForm();

    int ne      = bound ? ma.GetNSE() : ma.GetNE();
    int dim     = fes.GetDimension();
    int dimflux = fesflux.GetDimension();
    int dimfluxvec = bli.DimFlux(); 

    const BilinearFormIntegrator & fluxbli =
      bound ? (*fesflux.GetBoundaryEvaluator()) : (*fesflux.GetEvaluator());

    Array<int> cnti(fesflux.GetNDof());
    cnti = 0;

    flux.GetVector() = 0.0;

    int cnt = 0;
    clock_t prevtime = clock();

#pragma omp parallel 
    {
      LocalHeap lh = clh.Split();
      
      // ElementTransformation eltrans;
      Array<int> dnums, dnumsflux;

#pragma omp for 
      for (int i = 0; i < ne; i++)
	{
	  HeapReset hr(lh);
#pragma omp critical(fluxprojetpercent)
	  {
	    cnt++;
	    if (clock()-prevtime > 0.1 * CLOCKS_PER_SEC)
	      {
		cout << "\rpostprocessing element " << cnt << "/" << ne << flush;
		ma.SetThreadPercentage ( 100.0*cnt / ne );
		prevtime = clock();
	      }
	  }

	  int eldom = 
	    bound ? ma.GetSElIndex(i) : ma.GetElIndex(i);
	  
	  if(!domains[eldom])
	    continue;
	  
	  const FiniteElement & fel = 
	    bound ? fes.GetSFE(i, lh) : fes.GetFE (i, lh);
	  
	  const FiniteElement & felflux = 
	    bound ? fesflux.GetSFE(i, lh) : fesflux.GetFE (i, lh);
	  
	  ElementTransformation eltrans = ma.GetTrafo (i, bound);
	  
	  if (bound)
	    {
	      // ma.GetSurfaceElementTransformation (i, eltrans, lh);
	      fes.GetSDofNrs (i, dnums);
	      fesflux.GetSDofNrs (i, dnumsflux);
	    }
	  else
	    {
	      // ma.GetElementTransformation (i, eltrans, lh);
	      fes.GetDofNrs (i, dnums);
	      fesflux.GetDofNrs (i, dnumsflux);
	    }
	  
	  FlatVector<SCAL> elu(dnums.Size() * dim, lh);
	  FlatVector<SCAL> elflux(dnumsflux.Size() * dimflux, lh);
	  FlatVector<SCAL> elfluxi(dnumsflux.Size() * dimflux, lh);
	  FlatVector<SCAL> fluxi(dimfluxvec, lh);
	  
	  
	  u.GetElementVector (dnums, elu);
	  fes.TransformVec (i, bound, elu, TRANSFORM_SOL);
	  
	  const IntegrationRule & ir = 
	    SelectIntegrationRule(fel.ElementType(), max(fel.Order(),felflux.Order())+felflux.Order());


	  BaseMappedIntegrationRule & mir = eltrans(ir, lh);
	  FlatMatrix<SCAL> mfluxi(ir.GetNIP(), dimfluxvec, lh);
	  
	  bli.CalcFlux (fel, mir, elu, mfluxi, applyd, lh);
	  
	  for (int j = 0; j < ir.GetNIP(); j++)
	    mfluxi.Row(j) *= ir[j].Weight() * mir[j].GetMeasure();
	  
	  elflux = 0;
	  fluxbli.ApplyBTrans (felflux, mir, mfluxi, elflux, lh);


	      

	  if (dimflux > 1)
	    {
	      FlatMatrix<SCAL> elmat(dnumsflux.Size(), lh);
	      const BlockBilinearFormIntegrator & bbli = 
		dynamic_cast<const BlockBilinearFormIntegrator&> (fluxbli);
	      bbli . Block() . CalcElementMatrix (felflux, eltrans, elmat, lh);
	      FlatCholeskyFactors<SCAL> invelmat(elmat, lh);
	    
	      FlatVector<SCAL> hv1(dnumsflux.Size(), lh);
	      FlatVector<SCAL> hv2(dnumsflux.Size(), lh);
	      for (int j = 0; j < dimflux; j++)
		{
		  hv1 = elflux.Slice (j, dimflux);
		  invelmat.Mult (hv1, hv2);
		  elfluxi.Slice(j, dimflux) = hv2;
		}
	    }
	  else
	    {
	      FlatMatrix<SCAL> elmat(dnumsflux.Size(), lh);
	      fluxbli.CalcElementMatrix (felflux, eltrans, elmat, lh);
	      FlatCholeskyFactors<SCAL> invelmat(elmat, lh);
	      invelmat.Mult (elflux, elfluxi);
	    }
	  
	  fesflux.TransformVec (i, bound, elfluxi, TRANSFORM_SOL);
	  
	  
#pragma omp critical(fluxprojetadd)
	  {
	    flux.GetElementVector (dnumsflux, elflux);
	    elfluxi += elflux;
	    flux.SetElementVector (dnumsflux, elfluxi);
	    
	    for (int j = 0; j < dnumsflux.Size(); j++)
	      cnti[dnumsflux[j]]++;
	  }
	}
    }
    
    cout << "\rpostprocessing element " << ne << "/" << ne << endl;


    FlatVector<SCAL> fluxi(dimflux, clh);
    Array<int> dnumsflux(1);
    for (int i = 0; i < cnti.Size(); i++)
      if (cnti[i])
	{
	  dnumsflux[0] = i;
	  flux.GetElementVector (dnumsflux, fluxi);
	  fluxi /= double (cnti[i]);
	  flux.SetElementVector (dnumsflux, fluxi);
	}
    
    ma.PopStatus ();
  }


  
  template <class SCAL>
  void CalcFluxProject (const MeshAccess & ma, 
			const S_GridFunction<SCAL> & u,
			S_GridFunction<SCAL> & flux,
			const BilinearFormIntegrator & bli,
			bool applyd, int domain, LocalHeap & lh)
  {
    BitArray domains(ma.GetNDomains());
    
    if(domain == -1)
      domains.Set();
    else
      {
	domains.Clear();
	domains.Set(domain);
      }
    
    CalcFluxProject(ma,u,flux,bli,applyd,domains,lh);
  }


  void CalcFluxProject (const MeshAccess & ma, 
			const GridFunction & bu,
			GridFunction & bflux,
			const BilinearFormIntegrator & bli,
			bool applyd, int domain, LocalHeap & lh)
  {
    if (bu.GetFESpace().IsComplex())
      {
	CalcFluxProject (ma, 
			 dynamic_cast<const S_GridFunction<Complex>&> (bu),
			 dynamic_cast<S_GridFunction<Complex>&> (bflux),
			 bli, applyd, domain, lh);
      }
    else
      {
	CalcFluxProject (ma, 
			 dynamic_cast<const S_GridFunction<double>&> (bu),
			 dynamic_cast<S_GridFunction<double>&> (bflux),
			 bli, applyd, domain, lh);
      }
  }




  template <class SCAL> 
  int CalcPointFlux (const MeshAccess & ma, 
		     const GridFunction & bu,
		     const FlatVector<double> & point,
		     const Array<int> & domains,
		     FlatVector<SCAL> & flux,
		     const BilinearFormIntegrator & bli,
		     bool applyd,
		     LocalHeap & lh,
		     int component)// = 0)
  {
    HeapReset hr(lh);

    int elnr;
    Array<int> dnums;
    ElementTransformation eltrans;

    IntegrationPoint ip(0,0,0,1);

    bool boundary = bli.BoundaryForm();

    if(boundary)
      {
	if(domains.Size() > 0)
	  elnr = ma.FindSurfaceElementOfPoint(point,ip,false,&domains);
	else
	  elnr = ma.FindSurfaceElementOfPoint(point,ip,false);
      }
    else
      {
      	if(domains.Size() > 0)
	  elnr = ma.FindElementOfPoint(point,ip,false,&domains);
	else
	  elnr = ma.FindElementOfPoint(point,ip,false);
      }
    if (elnr < 0) return 0;

    const S_GridFunction<SCAL> & u = 
      dynamic_cast<const S_GridFunction<SCAL>&> (bu);
	
    const FESpace & fes = u.GetFESpace();
    const FiniteElement & fel = (boundary) ? fes.GetSFE (elnr, lh) : fes.GetFE (elnr, lh);

    if (boundary)
      {
	ma.GetSurfaceElementTransformation (elnr, eltrans, lh);
	fes.GetSDofNrs (elnr, dnums);
      }
    else
      {
	ma.GetElementTransformation (elnr, eltrans, lh);
	fes.GetDofNrs (elnr, dnums);
      }
	
    FlatVector<SCAL> elu(dnums.Size() * fes.GetDimension(), lh);
	
    if(bu.GetCacheBlockSize() == 1)
      {
	    u.GetElementVector (dnums, elu);
      }
    else
      {
	FlatVector<SCAL> elu2(dnums.Size() * fes.GetDimension() * bu.GetCacheBlockSize(), lh);
	u.GetElementVector (dnums,elu2);
	for(int i=0; i<elu.Size(); i++)
	  elu[i] = elu2[i*bu.GetCacheBlockSize()+component];
      }
    
    fes.TransformVec (elnr, boundary, elu, TRANSFORM_SOL);
	
    bli.CalcFlux (fel, eltrans(ip, lh), elu, flux, applyd, lh);
  
    /*
    if(boundary)
      {
	if(domains.Size() > 0)
	  elnr = ma.FindSurfaceElementOfPoint(point,ip,false,&domains);
	else
	  elnr = ma.FindSurfaceElementOfPoint(point,ip,false);
	
	if (elnr < 0) return 0;

	const S_GridFunction<SCAL> & u = 
	  dynamic_cast<const S_GridFunction<SCAL>&> (bu);
	
	const FESpace & fes = u.GetFESpace();
	const FiniteElement & fel = fes.GetSFE (elnr, lh);
	ma.GetSurfaceElementTransformation (elnr, eltrans, lh);
	
	fes.GetSDofNrs (elnr, dnums);
	
	FlatVector<SCAL> elu(dnums.Size() * fes.GetDimension(), lh);
	
	if(bu.GetCacheBlockSize() == 1)
	  {
	    u.GetElementVector (dnums, elu);
	  }
	else
	  {
	    FlatVector<SCAL> elu2(dnums.Size() * fes.GetDimension() * bu.GetCacheBlockSize(), lh);
	    u.GetElementVector (dnums,elu2);
	    for(int i=0; i<elu.Size(); i++)
	      elu[i] = elu2[i*bu.GetCacheBlockSize()+component];
	  }
	
	fes.TransformVec (elnr, true, elu, TRANSFORM_SOL);
	
	bli.CalcFlux (fel, eltrans(ip, lh), elu, flux, applyd, lh);
      }
    else
      {
	if(domains.Size() > 0)
	  elnr = ma.FindElementOfPoint(point,ip,false,&domains);
	else
	  elnr = ma.FindElementOfPoint(point,ip,false);
       
	if (elnr < 0) return 0;

	const S_GridFunction<SCAL> & u = 
	  dynamic_cast<const S_GridFunction<SCAL>&> (bu);
	
	const FESpace & fes = u.GetFESpace();
	const FiniteElement & fel = fes.GetFE (elnr, lh);
	ma.GetElementTransformation (elnr, eltrans, lh);
	
	fes.GetDofNrs (elnr, dnums);
	
	FlatVector<SCAL> elu(dnums.Size() * fes.GetDimension(), lh);
	
	if(bu.GetCacheBlockSize() == 1)
	  {
	    u.GetElementVector (dnums, elu);
	  }
	else
	  {
	    FlatVector<SCAL> elu2(dnums.Size() * fes.GetDimension() * bu.GetCacheBlockSize(), lh);
	    u.GetElementVector (dnums,elu2);
	    for(int i=0; i<elu.Size(); i++)
	      elu[i] = elu2[i*bu.GetCacheBlockSize()+component];
	  }
	
	fes.TransformVec (elnr, false, elu, TRANSFORM_SOL);
	bli.CalcFlux (fel, eltrans(ip, lh), elu, flux, applyd, lh);
      }
    */ 
    return 1;
  }
  

  template NGS_DLL_HEADER int CalcPointFlux<double> (const MeshAccess & ma, 
				      const GridFunction & u,
				      const FlatVector<double> & point,
				      const Array<int> & domains,
				      FlatVector<double> & flux,
				      const BilinearFormIntegrator & bli,
				      bool applyd,
				      LocalHeap & lh,
				      int component);
  
  template NGS_DLL_HEADER int CalcPointFlux<Complex> (const MeshAccess & ma, 
				       const GridFunction & u,
				       const FlatVector<double> & point,
				       const Array<int> & domains,
				       FlatVector<Complex> & flux,
				       const BilinearFormIntegrator & bli,
				       bool applyd,
				       LocalHeap & lh,
				       int component);
    

  
  template <class SCAL>
  int CalcPointFlux (const MeshAccess & ma, 
		     const GridFunction & bu,
		     const FlatVector<double> & point,
		     FlatVector<SCAL> & flux,
		     const BilinearFormIntegrator & bli,
		     bool applyd,
		     LocalHeap & lh,
		     int component)
  {
    Array<int> dummy;
    return CalcPointFlux(ma,bu,point,dummy,flux,bli,applyd,lh,component);
  }



  template NGS_DLL_HEADER int CalcPointFlux<double> (const MeshAccess & ma, 
				      const GridFunction & u,
				      const FlatVector<double> & point,
				      FlatVector<double> & flux,
				      const BilinearFormIntegrator & bli,
				      bool applyd,
				      LocalHeap & lh,
				      int component);
  
  template NGS_DLL_HEADER int CalcPointFlux<Complex> (const MeshAccess & ma, 
				       const GridFunction & u,
				       const FlatVector<double> & point,
				       FlatVector<Complex> & flux,
				       const BilinearFormIntegrator & bli,
				       bool applyd,
				       LocalHeap & lh,
				       int component);
    








  template <class SCAL>
  void SetValues (const MeshAccess & ma, 
		  const CoefficientFunction & coef,
		  GridFunction & bu,
		  bool bound,
		  DifferentialOperator * diffop,
		  LocalHeap & clh)
  {
    S_GridFunction<SCAL> & u = dynamic_cast<S_GridFunction<SCAL> &> (bu);

    ma.PushStatus ("setvalues");

    const FESpace & fes = u.GetFESpace();

    int ne      = bound ? ma.GetNSE() : ma.GetNE();
    int dim     = fes.GetDimension();

    const BilinearFormIntegrator & bli =
      bound ? (*fes.GetBoundaryEvaluator()) : (*fes.GetEvaluator());

    if (&bli == NULL)
      throw Exception ("no evaluator available");

    int dimflux = diffop ? diffop->Dim() : bli.DimFlux(); 
    
    Array<int> cnti(fes.GetNDof());
    cnti = 0;

    u.GetVector() = 0.0;

    int cnt = 0;
    clock_t prevtime = clock();


#pragma omp parallel 
    {
      LocalHeap lh = clh.Split();
      
      // Array<int> dnums;

#pragma omp for 
      for (int i = 0; i < ne; i++)
	{
	  HeapReset hr(lh);
#pragma omp critical(fluxprojetpercent)
	  {
	    cnt++;
	    if (clock()-prevtime > 0.1 * CLOCKS_PER_SEC)
	      {
		// cout << "\rsetvalues element " << cnt << "/" << ne << flush;
		ma.SetThreadPercentage ( 100.0*cnt / ne );
		prevtime = clock();
	      }
	  }


	  if (bound && !fes.IsDirichletBoundary(ma.GetSElIndex(i)))
	    continue;

	  const FiniteElement & fel = bound ? fes.GetSFE(i, lh) : fes.GetFE (i, lh);
	  Array<int> dnums(fel.GetNDof(), lh);
	  
	  ElementTransformation eltrans = ma.GetTrafo (i, bound); 
	  fes.GetDofNrs (i, dnums, bound);

	  FlatVector<SCAL> elflux(dnums.Size() * dim, lh);
	  FlatVector<SCAL> elfluxi(dnums.Size() * dim, lh);
	  FlatVector<SCAL> fluxi(dimflux, lh);
	  
	  const IntegrationRule & ir = 
	    GetIntegrationRules().SelectIntegrationRule(fel.ElementType(), 2*fel.Order());

	  FlatMatrix<SCAL> mfluxi(ir.GetNIP(), dimflux, lh);

	  BaseMappedIntegrationRule & mir = eltrans(ir, lh);
	  coef.Evaluate (mir, mfluxi);
	  

	  for (int j = 0; j < ir.GetNIP(); j++)
	    mfluxi.Row(j) *= ir[j].Weight() * mir[j].GetMeasure();

	  if (diffop)
	    diffop -> ApplyTrans (fel, mir, mfluxi, elflux, lh);
	  else
	    bli.ApplyBTrans (fel, mir, mfluxi, elflux, lh);

	  if (dim > 1)
	    {
	      FlatMatrix<SCAL> elmat(dnums.Size(), lh);
	      const BlockBilinearFormIntegrator & bbli = 
		dynamic_cast<const BlockBilinearFormIntegrator&> (bli);
	      bbli . Block() . CalcElementMatrix (fel, eltrans, elmat, lh);
	      FlatCholeskyFactors<SCAL> invelmat(elmat, lh);
	      
	      FlatVector<SCAL> hv1(dnums.Size(), lh);
	      FlatVector<SCAL> hv2(dnums.Size(), lh);
	      for (int j = 0; j < dim; j++)
		{
		  hv1 = elflux.Slice (j, dim);
		  invelmat.Mult (hv1, hv2);
		  elfluxi.Slice(j, dim) = hv2;
		}
	    }
	  else
	    {
	      FlatMatrix<SCAL> elmat(dnums.Size(), lh);
	      bli.CalcElementMatrix (fel, eltrans, elmat, lh);
	      FlatCholeskyFactors<SCAL> invelmat(elmat, lh);
	      invelmat.Mult (elflux, elfluxi);
	    }
	  
	  fes.TransformVec (i, bound, elfluxi, TRANSFORM_SOL);

#pragma omp critical(fluxprojetadd)
	  {
	    u.GetElementVector (dnums, elflux);
	    elfluxi += elflux;
	    u.SetElementVector (dnums, elfluxi);
	    
	    for (int j = 0; j < dnums.Size(); j++)
	      cnti[dnums[j]]++;
	  }
	}
    }


    // cout << "\rsetvalues element " << ne << "/" << ne << endl;


    FlatVector<SCAL> fluxi(dim, clh);
    Array<int> dnums(1);
    for (int i = 0; i < cnti.Size(); i++)
      if (cnti[i])
	{
	  dnums[0] = i;
	  u.GetElementVector (dnums, fluxi);
	  fluxi /= double (cnti[i]);
	  u.SetElementVector (dnums, fluxi);
	}
    
    ma.PopStatus ();
  }
  
  NGS_DLL_HEADER void SetValues (const MeshAccess & ma, 
				 const CoefficientFunction & coef,
				 GridFunction & u,
				 bool bound,
				 DifferentialOperator * diffop,
				 LocalHeap & clh)
  {
    if (u.GetFESpace().IsComplex())
      SetValues<Complex> (ma, coef, u, bound, diffop, clh);
    else
      SetValues<double> (ma, coef, u, bound, diffop, clh);
  }





  template <class SCAL>
  void CalcError (const MeshAccess & ma, 
		  const S_GridFunction<SCAL> & u,
		  const S_GridFunction<SCAL> & flux,
		  const BilinearFormIntegrator & bli,
		  FlatVector<double> & err,
		  const BitArray & domains, LocalHeap & lh)
  {
    static int timer = NgProfiler::CreateTimer ("CalcError");
    NgProfiler::RegionTimer reg (timer);

    ma.PushStatus ("Error estimator");

    const FESpace & fes = u.GetFESpace();
    const FESpace & fesflux = flux.GetFESpace();

    bool bound = bli.BoundaryForm();

    int ne      = bound ? ma.GetNSE() : ma.GetNE();
    int dim     = fes.GetDimension();
    int dimflux = fesflux.GetDimension();
    int dimfluxvec = bli.DimFlux(); // fesflux.GetDimension();

    const BilinearFormIntegrator & fluxbli =
      bound ? (*fesflux.GetBoundaryEvaluator()) : (*fesflux.GetEvaluator());

    ElementTransformation eltrans;

    Array<int> dnums;
    Array<int> dnumsflux;

    double sum = 0;
    for (int i = 0; i < ne; i++)
      {
	HeapReset hr(lh);
	ma.SetThreadPercentage ( 100.0*i / ne );

	int eldom = 
	  bound ? ma.GetSElIndex(i) : ma.GetElIndex(i);
	
	if (!domains[eldom])
	  continue;

	const FiniteElement & fel = 
	  bound ? fes.GetSFE(i, lh) : fes.GetFE (i, lh);

	const FiniteElement & felflux = 
	  (bound ? fesflux.GetSFE(i, lh) : fesflux.GetFE (i, lh));

	if (bound)
	  {
	    ma.GetSurfaceElementTransformation (i, eltrans, lh);
	    fes.GetSDofNrs (i, dnums);
	    fesflux.GetSDofNrs (i, dnumsflux);
	  }
	else
	  {
	    ma.GetElementTransformation (i, eltrans, lh);
	    fes.GetDofNrs (i, dnums);
	    fesflux.GetDofNrs (i, dnumsflux);
	  }

	FlatVector<SCAL> elu(dnums.Size() * dim, lh);
	FlatVector<SCAL> elflux(dnumsflux.Size() * dimflux, lh);
	FlatVector<SCAL> fluxi(dimfluxvec, lh);
	FlatVector<SCAL> fluxi2(dimfluxvec, lh);


	u.GetElementVector (dnums, elu);
	fes.TransformVec (i, bound, elu, TRANSFORM_SOL);
	flux.GetElementVector (dnumsflux, elflux);
	fesflux.TransformVec (i, bound, elflux, TRANSFORM_SOL);


	const IntegrationRule & ir = 
	  SelectIntegrationRule(felflux.ElementType(), 2*felflux.Order());


	FlatMatrix<SCAL> mfluxi(ir.GetNIP(), dimfluxvec, lh);
	FlatMatrix<SCAL> mfluxi2(ir.GetNIP(), dimfluxvec, lh);
	
	
	BaseMappedIntegrationRule & mir = eltrans(ir, lh);
	bli.CalcFlux (fel, mir, elu, mfluxi, 1, lh);
	fluxbli.CalcFlux (felflux, mir, elflux, mfluxi2, 0, lh);
	
	mfluxi -= mfluxi2;
	
	bli.ApplyDMatInv (fel, mir, mfluxi, mfluxi2, lh);
	
	double elerr = 0;
	for (int j = 0; j < ir.GetNIP(); j++)
	  elerr += ir[j].Weight() * mir[j].GetMeasure() *
	    fabs (InnerProduct (mfluxi.Row(j), mfluxi2.Row(j)));

	err(i) += elerr;
	sum += elerr;
      }
    ma.PopStatus ();
  }
  

  
  template <class SCAL>
  void CalcError (const MeshAccess & ma, 
		  const S_GridFunction<SCAL> & u,
		  const S_GridFunction<SCAL> & flux,
		  const BilinearFormIntegrator & bli,
		  FlatVector<double> & err,
		  int domain, LocalHeap & lh)
  {
    BitArray domains(ma.GetNDomains());
    
    if(domain == -1)
      domains.Set();
    else
      {
	domains.Clear();
	domains.Set(domain);
      }

    CalcError(ma,u,flux,bli,err,domains,lh);    
  }


  void CalcError (const MeshAccess & ma, 
		  const GridFunction & bu,
		  const GridFunction & bflux,
		  const BilinearFormIntegrator & bli,
		  FlatVector<double> & err,
		  int domain, LocalHeap & lh)
  {
    if (bu.GetFESpace().IsComplex())
      {
	CalcError (ma, 
		   dynamic_cast<const S_GridFunction<Complex>&> (bu),
		   dynamic_cast<const S_GridFunction<Complex>&> (bflux),
		   bli, err, domain, lh);
      }
    else
      {
	CalcError (ma, 
		   dynamic_cast<const S_GridFunction<double>&> (bu),
		   dynamic_cast<const S_GridFunction<double>&> (bflux),
		   bli, err, domain, lh);
      }
  }
  

  template <class SCAL>
  void CalcDifference (const MeshAccess & ma, 
		       const S_GridFunction<SCAL> & u1,
		       const S_GridFunction<SCAL> & u2,
		       const BilinearFormIntegrator & bli1,
		       const BilinearFormIntegrator & bli2,
		       FlatVector<double> & diff,
		       int domain, LocalHeap & lh)
  {
    ma.PushStatus ("Calc Difference");

    const FESpace & fes1 = u1.GetFESpace();
    const FESpace & fes2 = u2.GetFESpace();

    bool bound1 = bli1.BoundaryForm();
    bool bound2 = bli2.BoundaryForm();


    if(bound1!=bound2) 
      {
	cout << " ERROR: CalcDifference :: bli1.BoundaryForm != bl2.BoundaryForm there is something wrong?" << endl; 
	diff = 0; 
	return; 
      } 

    int ne      = bound1 ? ma.GetNSE() : ma.GetNE();
    int dim1    = fes1.GetDimension();
    int dim2    = fes2.GetDimension();
    int dimflux1 = bli1.DimFlux();
    int dimflux2 = bli2.DimFlux();

    if(dimflux1 != dimflux2) 
      { 
	cout << " ERROR: CalcDifference :: dimflux1 != dimflux2 !!!!! -> set diff = 0" << endl; 
	diff = 0; 
	return; 	
      } 

    ElementTransformation eltrans;

    bool applyd1 = 0;
    bool applyd2 = 0;

    Array<int> dnums1;
    Array<int> dnums2;

    double sum = 0;
    for (int i = 0; i < ne; i++)
      {
	ma.SetThreadPercentage ( 100.0*i / ne );

	lh.CleanUp();

	int eldom = 
	  bound1 ? ma.GetSElIndex(i) : ma.GetElIndex(i);
	
	if ((domain != -1) && (domain != eldom))
	  continue;

	const FiniteElement & fel1 = 
	  bound1 ? fes1.GetSFE(i, lh) : fes1.GetFE (i, lh);

	const FiniteElement & fel2 = 
	  bound1 ? fes2.GetSFE(i, lh) : fes2.GetFE (i, lh);

	if (bound1)
	  {
	    ma.GetSurfaceElementTransformation (i, eltrans, lh);
	    fes1.GetSDofNrs (i, dnums1);
	    fes2.GetSDofNrs (i, dnums2);
	  }
	else
	  {
	    ma.GetElementTransformation (i, eltrans, lh);
	    fes1.GetDofNrs (i, dnums1);
	    fes2.GetDofNrs (i, dnums2);
	  }

	FlatVector<SCAL> elu1(dnums1.Size() * dim1, lh);
	FlatVector<SCAL> elu2(dnums2.Size() * dim2, lh);
	FlatVector<SCAL> fluxi1(dimflux1, lh);
	FlatVector<SCAL> fluxi2(dimflux2, lh);


	u1.GetElementVector (dnums1, elu1);
	fes1.TransformVec (i, bound1, elu1, TRANSFORM_SOL);
	u2.GetElementVector (dnums2, elu2);
	fes2.TransformVec (i, bound2, elu2, TRANSFORM_SOL);

	double elerr = 0;

	int io = max(fel1.Order(),fel2.Order()); 

	const IntegrationRule & ir = 
	  GetIntegrationRules().SelectIntegrationRule(fel1.ElementType(), 
						      2*io+2);
	double det; 
	
	for (int j = 0; j < ir.GetNIP(); j++)
	  {
	    void * heapp = lh.GetPointer();
	    if (!bound1)
	      {
		if (ma.GetDimension() == 2)
		  {
		    SpecificIntegrationPoint<2,2> sip (ir[j], eltrans, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
		    bli2.CalcFlux (fel2, sip, elu2, fluxi2, applyd2, lh);
		    det = fabs(sip.GetJacobiDet()); 
		  }
		else
		  {
		    SpecificIntegrationPoint<3,3> sip (ir[j], eltrans, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
		    bli2.CalcFlux (fel2, sip, elu2, fluxi2, applyd2, lh);
		    det = fabs(sip.GetJacobiDet());  
		  }
	      }
	    else 
	      {
		if (ma.GetDimension() == 3)
		  {
		    SpecificIntegrationPoint<2,3> sip (ir[j], eltrans, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
		    bli2.CalcFlux (fel2, sip, elu2, fluxi2, applyd2, lh);
		    det = fabs(sip.GetJacobiDet()); 
		  }
		else
		  {
		    SpecificIntegrationPoint<1,2> sip (ir[j], eltrans, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
		    bli2.CalcFlux (fel2, sip, elu2, fluxi2, applyd2, lh);
		    det = fabs(sip.GetJacobiDet()); 
		  }
	      }

	    	  
	    fluxi1 -= fluxi2;
	     
	    double dx = ir[j].Weight() * det; 
	    
	    elerr += dx * L2Norm2 (fluxi1);
	    
	    lh.CleanUp (heapp);
	  }

	diff(i) += elerr;
	sum += elerr;
      }
    cout << "difference = " << sqrt(sum) << endl;
    ma.PopStatus ();
  }
  
  template void CalcDifference<double> (const MeshAccess & ma, 
					const S_GridFunction<double> & bu1,
					const S_GridFunction<double> & bu2,
					const BilinearFormIntegrator & bli1,
					const BilinearFormIntegrator & bli2,
					FlatVector<double> & err,
					int domain, LocalHeap & lh);
  
  template void CalcDifference<Complex> (const MeshAccess & ma, 
					 const S_GridFunction<Complex> & bu1,
					 const S_GridFunction<Complex> & bu2,
					 const BilinearFormIntegrator & bli1,
					 const BilinearFormIntegrator & bli2,
					 FlatVector<double> & err,
					 int domain, LocalHeap & lh);    
  




  template <class SCAL>
  void CalcDifference (const MeshAccess & ma, 
		       const S_GridFunction<SCAL> & u1,
		       const BilinearFormIntegrator & bli1,
		       const CoefficientFunction * coef, 
		       FlatVector<double> & diff,
		       int domain, LocalHeap & lh)
  {
    ma.PushStatus ("Calc Difference");

    const FESpace & fes1 = u1.GetFESpace();

    bool bound1 = bli1.BoundaryForm();


    int ne      = bound1 ? ma.GetNSE() : ma.GetNE();
    int dim1    = fes1.GetDimension();
    int dimflux1 = bli1.DimFlux();

    ElementTransformation eltrans;

    bool applyd1 = 0;

    Array<int> dnums1;

    double sum = 0;
    for (int i = 0; i < ne; i++)
      {
	ma.SetThreadPercentage ( 100.0*i / ne );

	lh.CleanUp();

	int eldom = 
	  bound1 ? ma.GetSElIndex(i) : ma.GetElIndex(i);
	
	if ((domain != -1) && (domain != eldom))
	  continue;

	const FiniteElement & fel1 = 
	  bound1 ? fes1.GetSFE(i, lh) : fes1.GetFE (i, lh);

	if (bound1)
	  {
	    ma.GetSurfaceElementTransformation (i, eltrans, lh);
	    fes1.GetSDofNrs (i, dnums1);
	  }
	else
	  {
	    ma.GetElementTransformation (i, eltrans, lh);
	    fes1.GetDofNrs (i, dnums1);
	  }

	FlatVector<SCAL> elu1(dnums1.Size() * dim1, lh);
	FlatVector<SCAL> fluxi1(dimflux1, lh);
	FlatVector<SCAL> fluxi2(dimflux1, lh);


	u1.GetElementVector (dnums1, elu1);
	fes1.TransformVec (i, bound1, elu1, TRANSFORM_SOL);

	double elerr = 0;

	const IntegrationRule & ir = 
	  GetIntegrationRules().SelectIntegrationRule(fel1.ElementType(), 
						      2*fel1.Order()+3);
	double det = 0;
	
	if (bound1) 
	  throw Exception ("CalcDifference on boundary not supported");

	if (ma.GetDimension() == 2)
	  {
	    MappedIntegrationRule<2,2> mir(ir, eltrans, lh);
	    FlatMatrix<SCAL> mfluxi(ir.GetNIP(), dimflux1, lh);
	    FlatMatrix<SCAL> mfluxi2(ir.GetNIP(), dimflux1, lh);
	    
	    bli1.CalcFlux (fel1, mir, elu1, mfluxi, applyd1, lh);
	    // coef->Evaluate(mir, mfluxi2);
	    
	    for (int j = 0; j < ir.GetNIP(); j++)
	      {
		coef->Evaluate(mir[j],fluxi2);
		mfluxi2.Row(j) = fluxi2;
	      }
	    mfluxi -= mfluxi2;

	    for (int j = 0; j < ir.GetNIP(); j++)
	      {
		double dx = fabs (mir[j].GetJacobiDet()) * ir[j].Weight();
		elerr += dx * L2Norm2 (mfluxi.Row(j));
	      }
	    

	    diff(i) += elerr;
	    sum += elerr;
	  }

	else
	  {
	    for (int j = 0; j < ir.GetNIP(); j++)
	      {
		HeapReset hr(lh);
		if (!bound1)
		  {
		    if (ma.GetDimension() == 2)
		      {
			Vec<2> point; 
			SpecificIntegrationPoint<2,2> sip (ir[j], eltrans, lh);
			eltrans.CalcPoint(sip.IP(), point, lh);
			bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
			const_cast<CoefficientFunction*>(coef)->Evaluate(sip,fluxi2);
			det = fabs(sip.GetJacobiDet()); 
		      }
		    else
		      {
			Vec<3> point;
			SpecificIntegrationPoint<3,3> sip (ir[j], eltrans, lh);
			eltrans.CalcPoint(sip.IP(), point, lh);
			bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
			const_cast<CoefficientFunction*>(coef)->Evaluate(sip,fluxi2);
			det = fabs(sip.GetJacobiDet());  
		      }
		  }

	    	  
		(*testout) << "diff: fluxi = " << fluxi1 << " =?= " << fluxi2 << endl;
	    
		fluxi1 -= fluxi2;
	     
		double dx = ir[j].Weight() * det; 
	    
		elerr += dx * L2Norm2 (fluxi1);
	      }

	    diff(i) += elerr;
	    sum += elerr;
	  }
      }
    cout << "difference = " << sqrt(sum) << endl;
    ma.PopStatus ();
  }

  NGS_DLL_HEADER void CalcDifference (const MeshAccess & ma, 
				      const GridFunction & u1,
				      const BilinearFormIntegrator & bfi1,
				      const CoefficientFunction * coef, 
				      FlatVector<double> & diff,
				      int domain, LocalHeap & lh)
  {
    if (u1.GetFESpace().IsComplex())
      CalcDifference (ma, 
		      dynamic_cast<const S_GridFunction<Complex>&> (u1), 
		      bfi1, coef, diff, domain, lh);
    else
      CalcDifference (ma, 
		      dynamic_cast<const S_GridFunction<double>&> (u1), 
		      bfi1, coef, diff, domain, lh);
  }

















  template <class SCAL>
  void CalcGradient (const MeshAccess & ma,
		     const FESpace & fesh1,
		     const S_BaseVector<SCAL> & vech1,
		     const FESpace & feshcurl,
		     S_BaseVector<SCAL> & vechcurl)
  {
    cout << "CalcGrad" << endl;
    const ScalarFiniteElement<2> * h1fe2d;
    const ScalarFiniteElement<3> * h1fe3d;
    const HCurlFiniteElement<2> * hcurlfe2d;
    const HCurlFiniteElement<3> * hcurlfe3d;

    h1fe2d = dynamic_cast<const ScalarFiniteElement<2>*> (&fesh1.GetFE(ET_TRIG));
    hcurlfe2d = dynamic_cast<const HCurlFiniteElement<2>*> (&feshcurl.GetFE(ET_TRIG));
    Matrix<> gradtrig(hcurlfe2d->GetNDof(), h1fe2d->GetNDof());
    ComputeGradientMatrix<2> (*h1fe2d, *hcurlfe2d, gradtrig);
    (*testout) << "gradtrig = " << gradtrig << endl;

    h1fe3d = dynamic_cast<const ScalarFiniteElement<3>*> (&fesh1.GetFE(ET_TET));
    hcurlfe3d = dynamic_cast<const HCurlFiniteElement<3>*> (&feshcurl.GetFE(ET_TET));
    Matrix<> gradtet(hcurlfe3d->GetNDof(), h1fe3d->GetNDof());
    ComputeGradientMatrix<3> (*h1fe3d, *hcurlfe3d, gradtet);
    (*testout) << "gradtet = " << gradtet << endl;


    int ne = ma.GetNE();
    Array<int> dnumsh1, dnumshcurl;
    LocalHeap lh(100000, "CalcGradient");
    
    for (int i = 0; i < ne; i++)
      {
	lh.CleanUp();
	fesh1.GetDofNrs (i, dnumsh1);
	feshcurl.GetDofNrs (i, dnumshcurl);

	FlatVector<SCAL> elhcurl(dnumshcurl.Size(), lh);
	FlatVector<SCAL> elh1(dnumsh1.Size(), lh);



	vech1.GetIndirect (dnumsh1, elh1);
	fesh1.TransformVec (i, 0, elh1, TRANSFORM_RHS);

	switch (fesh1.GetFE(i, lh).ElementType())
	  {
	  case ET_TRIG:
	    elhcurl = gradtrig * elh1;
	    break;
	  case ET_TET:
	    elhcurl = gradtet * elh1;
	    break;
	  default:
	    throw Exception ("CalcGradient: unsupported element");
	  }

	feshcurl.TransformVec (i, 0, elhcurl, TRANSFORM_RHS);
	vechcurl.SetIndirect (dnumshcurl, elhcurl);
      }
  }
  
  template
  void CalcGradient<double> (const MeshAccess & ma,
			     const FESpace & fesh1,
			     const S_BaseVector<double> & vech1,
			     const FESpace & feshcurl,
			     S_BaseVector<double> & vechcurl);







  
  template <class SCAL>
  void CalcGradientT (const MeshAccess & ma,
		      const FESpace & feshcurl,
		      const S_BaseVector<SCAL> & vechcurl1,
		      const FESpace & fesh1,
		      S_BaseVector<SCAL> & vech1)
  {
    cout << "CalcGrad" << endl;
    const ScalarFiniteElement<2> * h1fe2d;
    const ScalarFiniteElement<3> * h1fe3d;
    const HCurlFiniteElement<2> * hcurlfe2d;
    const HCurlFiniteElement<3> * hcurlfe3d;

    h1fe2d = dynamic_cast<const ScalarFiniteElement<2>*> (&fesh1.GetFE(ET_TRIG));
    hcurlfe2d = dynamic_cast<const HCurlFiniteElement<2>*> (&feshcurl.GetFE(ET_TRIG));
    Matrix<> gradtrig(hcurlfe2d->GetNDof(), h1fe2d->GetNDof());
    ComputeGradientMatrix<2> (*h1fe2d, *hcurlfe2d, gradtrig);
    (*testout) << "gradtrig = " << gradtrig << endl;

    h1fe3d = dynamic_cast<const ScalarFiniteElement<3>*> (&fesh1.GetFE(ET_TET));
    hcurlfe3d = dynamic_cast<const HCurlFiniteElement<3>*> (&feshcurl.GetFE(ET_TET));
    Matrix<> gradtet(hcurlfe3d->GetNDof(), h1fe3d->GetNDof());
    ComputeGradientMatrix<3> (*h1fe3d, *hcurlfe3d, gradtet);
    (*testout) << "gradtet = " << gradtet << endl;


    S_BaseVector<SCAL> & vechcurl =
      dynamic_cast<S_BaseVector<SCAL>&> (*vechcurl1.CreateVector());

    int ne = ma.GetNE();
    Array<int> dnumsh1, dnumshcurl;
    LocalHeap lh(100000, "CalcGradientT");
    
    vechcurl = vechcurl1;
    vech1.SetScalar(0); //  = SCAL(0);
    for (int i = 0; i < ne; i++)
      {
	lh.CleanUp();
	fesh1.GetDofNrs (i, dnumsh1);
	feshcurl.GetDofNrs (i, dnumshcurl);

	FlatVector<SCAL> elhcurl(dnumshcurl.Size(), lh);
	FlatVector<SCAL> elh1(dnumsh1.Size(), lh);

	vechcurl.GetIndirect (dnumshcurl, elhcurl);
	feshcurl.TransformVec (i, 0, elhcurl, TRANSFORM_RHS);

	switch (fesh1.GetFE(i, lh).ElementType())
	  {
	  case ET_TRIG:
	    elh1 = Trans (gradtrig) * elhcurl;
	    break;
	  case ET_TET:
	    elh1 = Trans (gradtet) * elhcurl;
	    break;
	  default:
	    throw Exception ("CalcGradientT: unsupported element");
	  }

	fesh1.TransformVec (i, 0, elh1, TRANSFORM_RHS);
	vech1.AddIndirect (dnumsh1, elh1);

	elhcurl = 0;
	vechcurl.SetIndirect (dnumshcurl, elhcurl);
      }
  }
  
  template
  void CalcGradientT<double> (const MeshAccess & ma,
			      const FESpace & feshcurl,
			      const S_BaseVector<double> & vechcurl,
			      const FESpace & fesh1,
			      S_BaseVector<double> & vech1);


  
}
