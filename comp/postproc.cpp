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
    ma.PushStatus ("Post-processing");

    const FESpace & fes = u.GetFESpace();
    const FESpace & fesflux = flux.GetFESpace();

    bool bound = bli.BoundaryForm();

    int ne      = bound ? ma.GetNSE() : ma.GetNE();
    int dim     = fes.GetDimension();
    int dimflux = fesflux.GetDimension();

    // ElementTransformation eltrans;

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
      
      ElementTransformation eltrans;
      Array<int> dnums, dnumsflux;

      void * heapp1 = lh.GetPointer();
#pragma omp for 
      for (int i = 0; i < ne; i++)
	{

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
	  FlatVector<SCAL> elfluxi(dnumsflux.Size() * dimflux, lh);
	  FlatVector<SCAL> fluxi(dimflux, lh);
	  
	  
	  u.GetElementVector (dnums, elu);
	  fes.TransformVec (i, bound, elu, TRANSFORM_SOL);
	  
	  const IntegrationRule & ir = 
	    GetIntegrationRules().SelectIntegrationRule(fel.ElementType(), max(fel.Order(),felflux.Order())+felflux.Order());
	  elflux = 0;
	  
	  for (int j = 0; j < ir.GetNIP(); j++)
	    {
	      HeapReset hr(lh);
	      // bli.CalcFlux (fel, eltrans, ir.GetIP(j), elu, fluxi, applyd, lh);
	      // fluxbli.ApplyBTrans (felflux, eltrans, ir.GetIP(j), fluxi, elfluxi, lh);
	      double fac;
	      if (!bound)
		{
		  if (fel.SpatialDim() == 2)
		    {
		      SpecificIntegrationPoint<2,2> sip (ir.GetIP(j), eltrans, lh);
		      fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());
		      bli.CalcFlux (fel, sip, elu, fluxi, applyd, lh);
		      fluxbli.ApplyBTrans (felflux, sip, fluxi, elfluxi, lh);
		    }
		  else
		    {
		      SpecificIntegrationPoint<3,3> sip (ir.GetIP(j), eltrans, lh);
		      fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());
		      bli.CalcFlux (fel, sip, elu, fluxi, applyd, lh);
		      fluxbli.ApplyBTrans (felflux, sip, fluxi, elfluxi, lh);
		    }
		}
	      else
		{
		  if (fel.SpatialDim() == 2)
		    {
		      SpecificIntegrationPoint<2,3> sip (ir.GetIP(j), eltrans, lh);
		      fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());
		      bli.CalcFlux (fel, sip, elu, fluxi, applyd, lh);
		      fluxbli.ApplyBTrans (felflux, sip, fluxi, elfluxi, lh);
		    }
		  else
		    {
		      SpecificIntegrationPoint<1,2> sip (ir.GetIP(j), eltrans, lh);
		      fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());
		      bli.CalcFlux (fel, sip, elu, fluxi, applyd, lh);
		      fluxbli.ApplyBTrans (felflux, sip, fluxi, elfluxi, lh);
		    }
		}
	      
	      elflux += fac * elfluxi;
	    }
	  
	  if (dimflux > 1)
	    {
	    FlatMatrix<SCAL> elmat(dnumsflux.Size(), lh);
	    const BlockBilinearFormIntegrator & bbli = 
	      dynamic_cast<const BlockBilinearFormIntegrator&> (fluxbli);
	    bbli . Block() . AssembleElementMatrix (felflux, eltrans, elmat, lh);
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
	      fluxbli.AssembleElementMatrix (felflux, eltrans, elmat, lh);
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

	  lh.CleanUp(heapp1);
	}
    }

    cout << "\rpostprocessing element xx " << ne << "/" << ne << endl;


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




  template void CalcFluxProject<double> (const MeshAccess & ma, 
					 const S_GridFunction<double> & bu,
					 S_GridFunction<double> & bflux,
					 const BilinearFormIntegrator & bli,
					 bool applyd, int domain, LocalHeap & lh);
  
  template void CalcFluxProject<Complex> (const MeshAccess & ma, 
					  const S_GridFunction<Complex> & bu,
					  S_GridFunction<Complex> & bflux,
					  const BilinearFormIntegrator & bli,
					  bool applyd, int domain, LocalHeap & lh);






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
    int elnr;
    //double lami[3];
    Array<int> dnums;
    ElementTransformation eltrans;

    IntegrationPoint ip(0,0,0,1);

    bool boundary = bli.BoundaryForm();

    //cout << " point " << point << endl; 
    if(boundary)
      {
	if(domains.Size() > 0)
	  elnr = ma.FindSurfaceElementOfPoint(point,ip,false,&domains);
	else
	  elnr = ma.FindSurfaceElementOfPoint(point,ip,false);

	
	if (elnr < 0) return 0;

	Array<int> vnums; 
	ma.GetSElVertices(elnr, vnums); 
	
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
	
	
	bli.CalcFlux (fel, eltrans, ip, elu, flux, applyd, lh);
      }
    else
      {
	if(domains.Size() > 0)
	  elnr = ma.FindElementOfPoint(point,ip,false,&domains);
	else
	  elnr = ma.FindElementOfPoint(point,ip,false);
       
	if (elnr < 0) return 0;

	Array<int> vnums; 
	ma.GetElVertices(elnr, vnums); 
	
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
	
	
	bli.CalcFlux (fel, eltrans, ip, elu, flux, applyd, lh);
      }
 
    return 1;
  }
  

  template int CalcPointFlux<double> (const MeshAccess & ma, 
				      const GridFunction & u,
				      const FlatVector<double> & point,
				      const Array<int> & domains,
				      FlatVector<double> & flux,
				      const BilinearFormIntegrator & bli,
				      bool applyd,
				      LocalHeap & lh,
				      int component);
  
  template int CalcPointFlux<Complex> (const MeshAccess & ma, 
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



  template int CalcPointFlux<double> (const MeshAccess & ma, 
				      const GridFunction & u,
				      const FlatVector<double> & point,
				      FlatVector<double> & flux,
				      const BilinearFormIntegrator & bli,
				      bool applyd,
				      LocalHeap & lh,
				      int component);
  
  template int CalcPointFlux<Complex> (const MeshAccess & ma, 
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
		  S_GridFunction<SCAL> & u,
		  bool bound,
		  LocalHeap & clh)
  {
    ma.PushStatus ("setvalues");

    const FESpace & fes = u.GetFESpace();
    // bool bound = bli.BoundaryForm();

    int ne      = bound ? ma.GetNSE() : ma.GetNE();
    int dim     = fes.GetDimension();

    const BilinearFormIntegrator & bli =
      bound ? (*fes.GetBoundaryEvaluator()) : (*fes.GetEvaluator());

    Array<int> cnti(fes.GetNDof());
    cnti = 0;

    u.GetVector() = 0.0;

    int cnt = 0;
    clock_t prevtime = clock();

#pragma omp parallel 
    {
      LocalHeap lh = clh.Split();
      
      ElementTransformation eltrans;
      Array<int> dnums;

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

	  // int eldom = 
	  // bound ? ma.GetSElIndex(i) : ma.GetElIndex(i);
	  
	  const FiniteElement & fel = 
	    bound ? fes.GetSFE(i, lh) : fes.GetFE (i, lh);
	  
	  if (bound)
	    {
	      ma.GetSurfaceElementTransformation (i, eltrans, lh);
	      fes.GetSDofNrs (i, dnums);
	    }
	  else
	    {
	      ma.GetElementTransformation (i, eltrans, lh);
	      fes.GetDofNrs (i, dnums);
	    }
	  
	  FlatVector<SCAL> elflux(dnums.Size() * dim, lh);
	  FlatVector<SCAL> elfluxi(dnums.Size() * dim, lh);
	  FlatVector<SCAL> fluxi(dim, lh);
	  
	  const IntegrationRule & ir = 
	    GetIntegrationRules().SelectIntegrationRule(fel.ElementType(), 2*fel.Order());
	  elflux = 0;
	  
	  for (int j = 0; j < ir.GetNIP(); j++)
	    {
	      HeapReset hr(lh);

	      double fac;
	      if (!bound)
		{
		  if (fel.SpatialDim() == 2)
		    {
		      SpecificIntegrationPoint<2,2> sip (ir.GetIP(j), eltrans, lh);
		      fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());
		      coef.Evaluate (sip, fluxi);
		      bli.ApplyBTrans (fel, sip, fluxi, elfluxi, lh);
		    }
		  else
		    {
		      SpecificIntegrationPoint<3,3> sip (ir.GetIP(j), eltrans, lh);
		      fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());
		      coef.Evaluate (sip, fluxi);
		      // bli.CalcFlux (fel, sip, elu, fluxi, applyd, lh);
		      bli.ApplyBTrans (fel, sip, fluxi, elfluxi, lh);
		    }
		}
	      else
		{
		  if (fel.SpatialDim() == 2)
		    {
		      SpecificIntegrationPoint<2,3> sip (ir.GetIP(j), eltrans, lh);
		      fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());
		      coef.Evaluate (sip, fluxi);
		      // bli.CalcFlux (fel, sip, elu, fluxi, applyd, lh);
		      bli.ApplyBTrans (fel, sip, fluxi, elfluxi, lh);
		    }
		  else
		    {
		      SpecificIntegrationPoint<1,2> sip (ir.GetIP(j), eltrans, lh);
		      fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());
		      coef.Evaluate (sip, fluxi);
		      // bli.CalcFlux (fel, sip, elu, fluxi, applyd, lh);
		      bli.ApplyBTrans (fel, sip, fluxi, elfluxi, lh);
		    }
		}
	      
	      elflux += fac * elfluxi;
	    }
	  
	  /*
	  if (dimflux > 1)
	    {
	      FlatMatrix<SCAL> elmat(dnumsflux.Size(), lh);
	      const BlockBilinearFormIntegrator & bbli = 
		dynamic_cast<const BlockBilinearFormIntegrator&> (fluxbli);
	      bbli . Block() . AssembleElementMatrix (felflux, eltrans, elmat, lh);
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
	  */
	    {
	      FlatMatrix<SCAL> elmat(dnums.Size(), lh);
	      bli.AssembleElementMatrix (fel, eltrans, elmat, lh);
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

    
    cout << "\rpostprocessing element " << ne << "/" << ne << endl;


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



  template void SetValues<double> (const MeshAccess & ma, 
				   const CoefficientFunction & coef,
				   S_GridFunction<double> & u,
				   bool bound,
				   LocalHeap & clh);

  /*
  template void SetValues<Complex> (const MeshAccess & ma, 
				    const CoefficientFunction & coef,
				    S_GridFunction<Complex> & u,
				    bool bound,
				    LocalHeap & clh);
  */


  template <class SCAL>
  void CalcError (const MeshAccess & ma, 
		  const S_GridFunction<SCAL> & u,
		  const S_GridFunction<SCAL> & flux,
		  const BilinearFormIntegrator & bli,
		  FlatVector<double> & err,
		  const BitArray & domains, LocalHeap & lh)
  {
    ma.PushStatus ("Error estimator");

    const FESpace & fes = u.GetFESpace();
    const FESpace & fesflux = flux.GetFESpace();

    bool bound = bli.BoundaryForm();

    int ne      = bound ? ma.GetNSE() : ma.GetNE();
    int dim     = fes.GetDimension();
    int dimflux = fesflux.GetDimension();

    const BilinearFormIntegrator & fluxbli =
      bound ? (*fesflux.GetBoundaryEvaluator()) : (*fesflux.GetEvaluator());

    ElementTransformation eltrans;

    Array<int> dnums;
    Array<int> dnumsflux;

    double sum = 0;
    for (int i = 0; i < ne; i++)
      {
	//	cout << "\rerror element " << i << "/" << ne << flush;
	ma.SetThreadPercentage ( 100.0*i / ne );

	lh.CleanUp();

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
	FlatVector<SCAL> fluxi(dimflux, lh);
	FlatVector<SCAL> fluxi2(dimflux, lh);


	u.GetElementVector (dnums, elu);
	fes.TransformVec (i, bound, elu, TRANSFORM_SOL);
	flux.GetElementVector (dnumsflux, elflux);
	fesflux.TransformVec (i, bound, elflux, TRANSFORM_SOL);

	double elerr = 0;

	const IntegrationRule & ir = 
	  SelectIntegrationRule(felflux.ElementType(), 2*felflux.Order());
	// double vol2 = ma.ElementVolume(i);
	// SZ sum(integration weights) trig = 1/2  

	double vol; 

	void * heapp = lh.GetPointer();

	for (int j = 0; j < ir.GetNIP(); j++)
	  {
	    if (!bound)
	      {
		if (fel.SpatialDim() == 2)
		  {
		    SpecificIntegrationPoint<2,2> sip (ir.GetIP(j), eltrans, lh);
		    bli.CalcFlux (fel, sip, elu, fluxi, 1, lh);
		    fluxbli.CalcFlux (felflux, sip, elflux, fluxi2, 0, lh);
		    vol = fabs(sip.GetJacobiDet()); 
		  }
		else
		  {
		    SpecificIntegrationPoint<3,3> sip (ir.GetIP(j), eltrans, lh);
		    bli.CalcFlux (fel, sip, elu, fluxi, 1, lh);
		    fluxbli.CalcFlux (felflux, sip, elflux, fluxi2, 0, lh);
		    vol = fabs(sip.GetJacobiDet()); 
		  }
	      }
	    else
	      {
		if (fel.SpatialDim() == 2)
		  {
		    SpecificIntegrationPoint<2,3> sip (ir.GetIP(j), eltrans, lh);
		    bli.CalcFlux (fel, sip, elu, fluxi, 1, lh);
		    fluxbli.CalcFlux (felflux, sip, elflux, fluxi2, 0, lh);
		    vol = fabs(sip.GetJacobiDet()); 
		  }
		else
		  {
		    SpecificIntegrationPoint<1,2> sip (ir.GetIP(j), eltrans, lh);
		    bli.CalcFlux (fel, sip, elu, fluxi, 1, lh);
		    fluxbli.CalcFlux (felflux, sip, elflux, fluxi2, 0, lh);
		    vol = fabs(sip.GetJacobiDet()); 
		  }
	      }

	    fluxi -= fluxi2;
	    
	    elerr += ir.GetIP(j).Weight() * vol * L2Norm2 (fluxi);


	    lh.CleanUp (heapp);
	  }



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


  template void CalcError<double> (const MeshAccess & ma, 
				   const S_GridFunction<double> & bu,
				   const S_GridFunction<double> & bflux,
				   const BilinearFormIntegrator & bli,
				   FlatVector<double> & err,
				   int domain, LocalHeap & lh);
  
  template void CalcError<Complex> (const MeshAccess & ma, 
				    const S_GridFunction<Complex> & bu,
				    const S_GridFunction<Complex> & bflux,
				    const BilinearFormIntegrator & bli,
				    FlatVector<double> & err,
				    int domain, LocalHeap & lh);    

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
		if (fel1.SpatialDim() == 2)
		  {
		    SpecificIntegrationPoint<2,2> sip (ir.GetIP(j), eltrans, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
		    bli2.CalcFlux (fel2, sip, elu2, fluxi2, applyd2, lh);
		    det = fabs(sip.GetJacobiDet()); 
		  }
	 	else
		  {
		    SpecificIntegrationPoint<3,3> sip (ir.GetIP(j), eltrans, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
		    bli2.CalcFlux (fel2, sip, elu2, fluxi2, applyd2, lh);
		    det = fabs(sip.GetJacobiDet());  
		  }
	      }
	    else 
	      {
		if (fel1.SpatialDim() == 2)
		  {
		    SpecificIntegrationPoint<2,3> sip (ir.GetIP(j), eltrans, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
		    bli2.CalcFlux (fel2, sip, elu2, fluxi2, applyd2, lh);
		    det = fabs(sip.GetJacobiDet()); 
		  }
		else
		  {
		    SpecificIntegrationPoint<1,2> sip (ir.GetIP(j), eltrans, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
		    bli2.CalcFlux (fel2, sip, elu2, fluxi2, applyd2, lh);
		    det = fabs(sip.GetJacobiDet()); 
		  }
	      }

	    	  
	    //  (*testout) << "diff: fluxi = " << fluxi1 << " =?= " << fluxi2 << endl;
	    
	    fluxi1 -= fluxi2;
	     
	    double dx = ir.GetIP(j).Weight() * det; 
	    
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
		       const CoefficientFunction * coef_real, 
		       const CoefficientFunction * coef_imag,
		       FlatVector<double> & diff,
		       int domain, LocalHeap & lh)
  {
    ;
  }



  template <>
  void CalcDifference (const MeshAccess & ma, 
		       const S_GridFunction<double> & u1,
		       const BilinearFormIntegrator & bli1,
		       const CoefficientFunction * coef_real, 
		       const CoefficientFunction * coef_imag,
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

	FlatVector<double> elu1(dnums1.Size() * dim1, lh);
	FlatVector<double> fluxi1(dimflux1, lh);
	FlatVector<double> fluxi2(dimflux1, lh);


	u1.GetElementVector (dnums1, elu1);
	fes1.TransformVec (i, bound1, elu1, TRANSFORM_SOL);

	double elerr = 0;

	int io = 5 + fel1.Order();

	const IntegrationRule & ir = 
	  GetIntegrationRules().SelectIntegrationRule(fel1.ElementType(), 
						      2*io+2);
	double det = 0;
	
	for (int j = 0; j < ir.GetNIP(); j++)
	  {
	    void * heapp = lh.GetPointer();
	    if (!bound1)
	      {
		if (fel1.SpatialDim() == 2)
		  {
		    Vec<2> point; 
		    SpecificIntegrationPoint<2,2> sip (ir.GetIP(j), eltrans, lh);
		    eltrans.CalcPoint(sip.IP(), point, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
		    fluxi2(0) = const_cast<CoefficientFunction*>(coef_real)->Evaluate(sip);
		    det = fabs(sip.GetJacobiDet()); 
		  }
	 	else
		  {
		    Vec<3> point;
		    SpecificIntegrationPoint<3,3> sip (ir.GetIP(j), eltrans, lh);
		    eltrans.CalcPoint(sip.IP(), point, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);
		    fluxi2(0) = const_cast<CoefficientFunction*>(coef_real)->Evaluate(sip);
		    det = fabs(sip.GetJacobiDet());  
		  }
	      }

	    	  
	    (*testout) << "diff: fluxi = " << fluxi1 << " =?= " << fluxi2 << endl;
	    
	    fluxi1 -= fluxi2;
	     
	    double dx = ir.GetIP(j).Weight() * det; 
	    
	    elerr += dx * L2Norm2 (fluxi1);
	    
	    lh.CleanUp (heapp);
	  }

	diff(i) += elerr;
	sum += elerr;
      }
    cout << "difference = " << sqrt(sum) << endl;
    ma.PopStatus ();
  }
  




  template <>
  void CalcDifference (const MeshAccess & ma, 
		       const S_GridFunction<Complex> & u1,
		       const BilinearFormIntegrator & bli1,
		       const CoefficientFunction * coef_real, 
		       const CoefficientFunction * coef_imag,
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

	FlatVector<Complex> elu1(dnums1.Size() * dim1, lh);
	FlatVector<Complex> fluxi1(dimflux1, lh);
	FlatVector<Complex> fluxi2(dimflux1, lh);


	u1.GetElementVector (dnums1, elu1);
	fes1.TransformVec (i, bound1, elu1, TRANSFORM_SOL);

	double elerr = 0;

	int io = 5 + fel1.Order();

	const IntegrationRule & ir = 
	  GetIntegrationRules().SelectIntegrationRule(fel1.ElementType(), 
						      2*io+2);
	double det = 0; 
	
	for (int j = 0; j < ir.GetNIP(); j++)
	  {
	    void * heapp = lh.GetPointer();
	    if (!bound1)
	      {
		if (fel1.SpatialDim() == 2)
		  {
		    Vec<2> point; 
		    SpecificIntegrationPoint<2,2> sip (ir.GetIP(j), eltrans, lh);
		    eltrans.CalcPoint(sip.IP(), point, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);

                    double real, imag = 0;
		    real = const_cast<CoefficientFunction*>(coef_real)->Evaluate(sip);
		    if ( coef_imag ) imag = const_cast<CoefficientFunction*>(coef_imag)->Evaluate(sip);
                    fluxi2(0) = Complex(real,imag);

		    det = fabs(sip.GetJacobiDet()); 
		  }
	 	else
		  {
		    Vec<3> point;
		    SpecificIntegrationPoint<3,3> sip (ir.GetIP(j), eltrans, lh);
		    eltrans.CalcPoint(sip.IP(), point, lh);
		    bli1.CalcFlux (fel1, sip, elu1, fluxi1, applyd1, lh);

		    // fluxi2(0) = const_cast<CoefficientFunction*>(coef_real)->Evaluate(sip);
		    // if ( coef_imag ) fluxi2(0).imag() = const_cast<CoefficientFunction*>(coef_imag)->Evaluate(sip);

                    double real, imag = 0;
		    real = const_cast<CoefficientFunction*>(coef_real)->Evaluate(sip);
		    if ( coef_imag ) imag = const_cast<CoefficientFunction*>(coef_imag)->Evaluate(sip);
                    fluxi2(0) = Complex(real,imag);
		    det = fabs(sip.GetJacobiDet());  
		  }
	      }

	    	  
	    (*testout) << "diff: fluxi = " << fluxi1 << " =?= " << fluxi2 << endl;
	    
	    fluxi1 -= fluxi2;
	     
	    double dx = ir.GetIP(j).Weight() * det; 
	    
	    elerr += dx * L2Norm2 (fluxi1);
	    
	    lh.CleanUp (heapp);
	  }

	diff(i) += elerr;
	sum += elerr;
      }
    cout << "difference = " << sqrt(sum) << endl;
    ma.PopStatus ();
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
    LocalHeap lh(100000);
    
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
    LocalHeap lh(100000);
    
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
