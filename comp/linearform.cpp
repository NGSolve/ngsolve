#include <comp.hpp>


#include <parallelngs.hpp>

namespace ngcomp
{
  using namespace ngcomp;
  using namespace ngparallel;

  LinearForm ::
  LinearForm (const FESpace & afespace, 
	      const string & aname,
	      const Flags & flags)
    : NGS_Object(afespace.GetMeshAccess(), aname), fespace(afespace)
  {
    independent = false;
    print = flags.GetDefineFlag ("print");
    printelvec = flags.GetDefineFlag ("printelvec");
    assembled = false;
    initialassembling = true;
    cacheblocksize = 1;
  }

  LinearForm :: ~LinearForm ()
  { 
    for (int i = 0; i < parts.Size(); i++)
      if ( parts_deletable[i]) delete parts[i];
  }

  void LinearForm :: AddIntegrator (LinearFormIntegrator * lfi, const bool deletable)
  {
    parts.Append (lfi);
    parts_deletable.Append(deletable);
  }

  void LinearForm :: PrintReport (ostream & ost)
  {
    ost << "on space " << GetFESpace().GetName() << endl
	<< "integrators: " << endl;
    for (int i = 0; i < parts.Size(); i++)
      ost << "  " << parts[i]->Name() << endl;
  }

  void LinearForm :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    if (&GetVector())  
      {
	int olds = mu.Size();
	GetVector().MemoryUsage (mu);
	for (int i = olds; i < mu.Size(); i++)
	  mu[i]->AddName (string(" lf ")+GetName());
      }
  }

  bool LinearForm :: InitialAssembling (void) {return initialassembling; }
  void LinearForm :: SetNoInitialAssembling (void) {initialassembling =  false; }

  bool LinearForm :: IsAssembled (void) { return assembled; }

  template <class SCAL>
  BaseVector & T_LinearForm<SCAL> :: GetVector () const 
  { 
#ifdef PARALLEL
    const FESpace & afespace = LinearForm :: GetFESpace();
    vec -> SetParallelDofs ( &afespace.GetParallelDofs() );
    if ( &(afespace.GetParallelDofs()) )
      vec -> SetStatus (DISTRIBUTED);
    else
      vec -> SetStatus ( NOT_PARALLEL );
#endif
    return *vec; 
  }

  template <class SCAL>
  void S_LinearForm<SCAL> :: Assemble (LocalHeap & clh)
  {
    static int vectimer = NgProfiler::CreateTimer ("Vector assembling");
    NgProfiler::RegionTimer reg (vectimer);

    assembled = true;

    if (independent)
      {
	AssembleIndependent(clh);
	return;
      }
    
    try
      {
	ma.PushStatus ("Assemble Vector");

	AllocateVector();

	
	int ne = ma.GetNE();
	int nse = ma.GetNSE();

	bool hasbound = false;
	bool hasinner = false;
	
	for (int j = 0; j < parts.Size(); j++)
	  {
	    if (parts[j] -> BoundaryForm())
	      hasbound = true;
	    else if (!parts[j] -> IntegrationAlongCurve())
	      hasinner = true;
	  }
	
	clock_t prevtime = clock();
	if (hasinner)
	  {
	    int cnt = 0;

#pragma omp parallel
	    {

	      LocalHeap lh = clh.Split();
	      
	      Array<int> dnums;
	      ElementTransformation eltrans;
	      
#pragma omp for	      
	      for (int i = 0; i < ne; i++)
		{
#pragma omp critical(printvecasstatus)
		  {
		    cnt++;
		    if (clock()-prevtime > 0.1 * CLOCKS_PER_SEC)
		      {
			cout << "\rassemble element " << cnt << "/" << ne << flush;
			ma.SetThreadPercentage ( 100.0*cnt / (ne+nse) );
			prevtime = clock();
		      }
		  }
		  
		  HeapReset hr(lh);

		  if ( ma.IsGhostEl( i ) ) continue;

		  const FiniteElement & fel = fespace.GetFE (i, lh);
		  ma.GetElementTransformation (i, eltrans, lh);
		      
		  fespace.GetDofNrs (i, dnums);
		
		  for (int j = 0; j < parts.Size(); j++)
		    {
		      if (parts[j] -> BoundaryForm()) continue;
		      if (!parts[j] -> DefinedOn (ma.GetElIndex (i))) continue;
		      if (parts[j] -> IntegrationAlongCurve()) continue;

                      int elvec_size = dnums.Size()*fespace.GetDimension();
		      FlatVector<TSCAL> elvec(elvec_size, lh);
		  
		      parts[j] -> AssembleElementVector (fel, eltrans, elvec, lh);
                 
		      if (printelvec)
			{
			  testout->precision(8);

			  (*testout) << "elnum= " << i << endl;
			  (*testout) << "integrator " << parts[j]->Name() << endl;
			  (*testout) << "dnums = " << endl << dnums << endl;
			  (*testout) << "element-index = " << eltrans.GetElementIndex() << endl;
			  (*testout) << "elvec = " << endl << elvec << endl;
			}
		
		      fespace.TransformVec (i, false, elvec, TRANSFORM_RHS);
		
#pragma omp critical (addelvec1) 
		      {
			AddElementVector (dnums, elvec, parts[j]->CacheComp()-1);
		      } 
		    
		    }
		}
	    }
	    cout << "\rassemble element " << ne << "/" << ne << endl;
	    
	  }


	if (hasbound)
	  {

#pragma omp parallel
	    {
	      LocalHeap lh = clh.Split();
	      Array<int> dnums;
	      ElementTransformation eltrans;
	      
#pragma omp for	      


	      for (int i = 0; i < nse; i++)
		{
		  if (i % 10 == 0)
		    cout << "\rassemble surface element " << i << "/" << nse << flush;
		  ma.SetThreadPercentage ( 100.0*(ne+i+1) / (ne+nse) );

		  if ( ma.IsGhostSEl ( i ) ) continue;

		  lh.CleanUp();
	      
		  const FiniteElement & fel = fespace.GetSFE (i, lh);
	      
		  ma.GetSurfaceElementTransformation (i, eltrans, lh);
		  fespace.GetSDofNrs (i, dnums);
	      
		  for (int j = 0; j < parts.Size(); j++)
		    {
		      if (!parts[j] -> BoundaryForm()) continue;
		      if (!parts[j] -> DefinedOn (ma.GetSElIndex (i))) continue;
		      if (parts[j] -> IntegrationAlongCurve()) continue;		    
		  
                      int elvec_size = dnums.Size()*fespace.GetDimension();
		      FlatVector<TSCAL> elvec(elvec_size, lh);
		      parts[j] -> AssembleElementVector (fel, eltrans, elvec, lh);
		      if (printelvec)
			{
			  testout->precision(8);

			  (*testout) << "surface-elnum= " << i << endl;
			  (*testout) << "integrator " << parts[j]->Name() << endl;
			  (*testout) << "dnums = " << endl << dnums << endl;
			  (*testout) << "element-index = " << eltrans.GetElementIndex() << endl;
			  (*testout) << "elvec = " << endl << elvec << endl;
			}



		      fespace.TransformVec (i, true, elvec, TRANSFORM_RHS);
		    
#pragma omp critical (addelvec2)		    
		      {
			AddElementVector (dnums, elvec, parts[j]->CacheComp()-1);
		      }
		    }
		}
	    }
	    cout << "\rassemble surface element " << nse << "/" << nse << endl;	  
	  }


	Array<int> dnums;
	ElementTransformation eltrans;

	for (int j=0; j<parts.Size(); j++ )
	  {
	    if (!(parts[j] -> IntegrationAlongCurve())) continue;
	    
	    const FiniteElement * fel = NULL;
	    int oldelement;
	    int element;
	    double oldlength;

	    
	    Array<int> domains;

	    if(parts[j]->DefinedOnSubdomainsOnly())
	      {
		for(int i=0; i<ma.GetNDomains(); i++)
		  if(parts[j]->DefinedOn(i))
		    domains.Append(i);
	      }
	    
	    
	    for(int nc = 0; nc < parts[j]->GetNumCurveParts(); nc++)
	      {
		oldelement = -1;
		oldlength = 0;

		for(int i = parts[j]->GetStartOfCurve(nc); i < parts[j]->GetEndOfCurve(nc); i++)
		  {
		    
		    if (i%500 == 0)
		      {
			cout << "\rassemble curvepoint " << i << "/" << parts[j]->NumCurvePoints() << flush;
			ma.SetThreadPercentage(100.*i/parts[j]->NumCurvePoints());
		      }
		    
		    FlatVector<TSCAL> elvec;
		    IntegrationPoint ip;
		    if(domains.Size() > 0)
		      element = ma.FindElementOfPoint(parts[j]->CurvePoint(i),ip,true,&domains);
		    else
		      element = ma.FindElementOfPoint(parts[j]->CurvePoint(i),ip,true);
		    if(element < 0)
		      throw Exception("element for curvepoint not found");
		    
		    if(element != oldelement)
		      { 
			clh.CleanUp();
			fel = &fespace.GetFE(element,clh);
			ma.GetElementTransformation(element,eltrans,clh);
			fespace.GetDofNrs(element,dnums);
		      }
		    
		    void * heapp = clh.GetPointer();
		    
		    FlatVector<double> tangent(parts[j]->CurvePoint(0).Size(),clh);
		    tangent = parts[j]->CurvePointTangent(i);
		    double length = L2Norm(tangent);
		    if(length < 1e-15)
		      {
			int i1,i2;
			i1 = i-1; i2=i+1;
			if(i1 < parts[j]->GetStartOfCurve(nc)) i1=parts[j]->GetStartOfCurve(nc);
			if(i2 >= parts[j]->GetEndOfCurve(nc)) i2 = parts[j]->GetEndOfCurve(nc)-1;
			
			for(int k=0; k<tangent.Size(); k++)
			  tangent[k] = (parts[j]->CurvePoint(i2))[k]-(parts[j]->CurvePoint(i1))[k];
			
			length = L2Norm(tangent);
		      }
		    
		    tangent *= 1./length;
		    
		    length = 0;
		    if(i < parts[j]->GetEndOfCurve(nc)-1)
		      for(int k=0; k<tangent.Size(); k++)
			length += pow(parts[j]->CurvePoint(i+1)[k]-parts[j]->CurvePoint(i)[k],2);
		    length = 0.5*sqrt(length);
		    
		    
		    if (eltrans.SpaceDim() == 3)
		      {
			SpecificIntegrationPoint<1,3> s_sip(ip,eltrans,clh);
			SpecificIntegrationPoint<3,3> g_sip(ip,eltrans,clh);
			Vec<3> tv;
			tv(0) = tangent(0); tv(1) = tangent(1); tv(2) = tangent(2);
			s_sip.SetTV(tv);
			parts[j]->AssembleElementVectorIndependent(*fel,
								   s_sip,
								   g_sip,
								   elvec,clh,true);
		      }
		    else if (eltrans.SpaceDim() == 2)
		      {
			SpecificIntegrationPoint<1,2> s_sip(ip,eltrans,clh);
			SpecificIntegrationPoint<2,2> g_sip(ip,eltrans,clh);
			Vec<2> tv;
			tv(0) = tangent(0); tv(1) = tangent(1);
			s_sip.SetTV(tv);
			parts[j]->AssembleElementVectorIndependent(*fel,
								   s_sip,
								   g_sip,
								   elvec,clh,true);
		      }
		    fespace.TransformVec (element, false, elvec, TRANSFORM_RHS);

		    elvec *= (oldlength+length); // Des is richtig
		    
		    
		    AddElementVector (dnums, elvec, parts[j]->CacheComp()-1);
		    
		    
		    
		    
		    clh.CleanUp(heapp);
		    
		    oldlength = length;      
		    
		    oldelement = element;
		  }
	      }



	    cout << "\rassemble curvepoint " << parts[j]->NumCurvePoints() << "/" 
		 << parts[j]->NumCurvePoints() << endl;
	    clh.CleanUp();

	  }
	
	




	if (print)
	  {
	    (*testout) << "Linearform " << GetName() << ": " << endl;
	    (*testout) << GetVector() << endl;
	  }

	ma.PopStatus ();
	// (*testout) << "Linearform, vec = " << endl << GetVector() << endl;


      }
    catch (Exception & e)
      {
	stringstream ost;
	ost << "in Assemble LinearForm" << endl;
	e.Append (ost.str());
	throw;
      }
    catch (exception & e)
      {
	throw (Exception (string(e.what()) +
			  string("\n in Assemble LinearForm\n")));
      }
  }








  template <class SCAL>
  void S_LinearForm<SCAL> :: AssembleIndependent (LocalHeap & lh)
  {
    assembled = true;

    try
      {
	AllocateVector();

	// int ne = ma.GetNE();
	int nse = ma.GetNSE();
	
	
	Array<int> dnums;
	ElementTransformation seltrans, geltrans;

	for (int i = 0; i < nse; i++)
	  {
	    lh.CleanUp();
	    
	    const FiniteElement & sfel = fespace.GetSFE (i, lh);
	    ma.GetSurfaceElementTransformation (i, seltrans, lh);
	      	
	    // (*testout) << "el = " << i << ", ind = " << ma.GetSElIndex(i) << endl;
	    if (!parts[0]->DefinedOn (ma.GetSElIndex(i))) continue;
	    // (*testout) << "integrate surf el " << endl;
	    
	    const IntegrationRule & ir =
	      GetIntegrationRules().SelectIntegrationRule (sfel.ElementType(), 5);
	    
	    for (int j = 0; j < ir.GetNIP(); j++)
	      {
		const IntegrationPoint & ip = ir.GetIP(j);
		SpecificIntegrationPoint<2,3> sip(ip, seltrans, lh);
		
		// (*testout) << "point = " << sip.GetPoint() << endl;
		
		IntegrationPoint gip;
		int elnr;
		// elnr = ma.FindElementOfPoint (FlatVector<>(sip.GetPoint()), gip, 1);
                elnr = ma.FindElementOfPoint (sip.GetPoint(), gip, 1);
		
		// (*testout) << "elnr = " << elnr << endl;
		if (elnr == -1) continue;
		
		const FiniteElement & gfel = fespace.GetFE (elnr, lh);
		ma.GetElementTransformation (elnr, geltrans, lh);
		SpecificIntegrationPoint<3,3> gsip(gip, geltrans, lh);
		
		// (*testout) << " =?= p = " << gsip.GetPoint() << endl;

		fespace.GetDofNrs (elnr, dnums);
		
		for (int k = 0; k < parts.Size(); k++)
		  {
		    FlatVector<TSCAL> elvec;
		    
		    if(geltrans.SpaceDim() == 3)
		      {
			SpecificIntegrationPoint<2,3> s_sip(ip,seltrans,lh);
			SpecificIntegrationPoint<3,3> g_sip(gip,geltrans,lh);
			parts[k] -> AssembleElementVectorIndependent
			  (gfel,s_sip,g_sip,elvec,lh);
		      }
		    else if(geltrans.SpaceDim() == 2)
		      {
			SpecificIntegrationPoint<1,2> s_sip(ip,seltrans,lh);
			SpecificIntegrationPoint<2,2> g_sip(gip,geltrans,lh);
			parts[k] -> AssembleElementVectorIndependent
			  (gfel,s_sip,g_sip,elvec,lh);
		      }

		    // (*testout) << "elvec, 1 = " << elvec << endl;

		    elvec *= fabs (sip.GetJacobiDet()) * ip.Weight();
		    fespace.TransformVec (elnr, 0, elvec, TRANSFORM_RHS);

		    // (*testout) << "final vec = " << elvec << endl;
		    // (*testout) << "dnums = " << dnums << endl;
		    AddElementVector (dnums, elvec);
		  }
	      }
	  }
	// (*testout) << "Assembled vector = " << endl << GetVector() << endl;
      }
    
    catch (Exception & e)
      {
	stringstream ost;
	ost << "in Assemble LinearForm Independent" << endl;
	e.Append (ost.str());
	throw;
      }
    catch (exception & e)
      {
	throw (Exception (string(e.what()) +
			  string("\n in Assemble LinearForm Independent\n")));
      }
  }

















  template <typename TV>
  T_LinearForm<TV> ::
  T_LinearForm (const FESpace & afespace, const string & aname, const Flags & flags)
    : S_LinearForm<TSCAL> (afespace, aname, flags), vec(0) 
  { 
    ; 
  }
  
  template <typename TV>
  T_LinearForm<TV> :: ~T_LinearForm ()
  {
    delete vec;
  }

#ifndef PARALLEL
  template <typename TV>
  void T_LinearForm<TV> :: AllocateVector ()
  {
    delete vec;
    vec = new ngla::VVector<TV> (this->fespace.GetNDof());
    (*vec) = TSCAL(0);
  }
#else
  template <typename TV>
  void T_LinearForm<TV> :: AllocateVector ()
  {
    delete vec;
    const FESpace & afespace = this->fespace;
    if ( &afespace.GetParallelDofs() == 0 )
      vec = new VVector<TV> (afespace.GetNDof());
    else
      vec = new ParallelVVector<TV> ( afespace.GetNDof(),& afespace.GetParallelDofs());
    (*vec) = TSCAL(0);
  }
#endif


  template <typename TV>
  void T_LinearForm<TV> :: CleanUpLevel ()
  {
    delete vec;
    vec = 0;
  }


  template <typename TV>
  void T_LinearForm<TV> ::
  AddElementVector (const Array<int> & dnums,
		    const FlatVector<TSCAL> & elvec,
		    const int cachecomp) 
  {
    FlatVector<TV> fv = vec->FV();

    if(cachecomp < 0)
      {
	for (int k = 0; k < dnums.Size(); k++)
	  if (dnums[k] != -1)
	    for (int j = 0; j < HEIGHT; j++)
	      fv(dnums[k])(j) += elvec(k*HEIGHT+j);
      }
    else
      {
	for (int k = 0; k < dnums.Size(); k++)
	  if (dnums[k] != -1)
	    fv(dnums[k])(cachecomp) += elvec(k);
      }
  }
  
  template <> void T_LinearForm<double>::
  AddElementVector (const Array<int> & dnums,
		    const FlatVector<double> & elvec,
		    const int cachecomp) 
  {
    FlatVector<double> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (dnums[k] != -1)
	fv(dnums[k]) += elvec(k);
  }
  
  template <> void T_LinearForm<Complex>::
  AddElementVector (const Array<int> & dnums,
		    const FlatVector<Complex> & elvec,
		    const int cachecomp) 
  {
    FlatVector<Complex> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (dnums[k] != -1)
	fv(dnums[k]) += elvec(k);
  }
  




  template <typename TV>
  void T_LinearForm<TV> ::
  SetElementVector (const Array<int> & dnums,
		    const FlatVector<TSCAL> & elvec) 
  {
    FlatVector<TV> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (dnums[k] != -1)
	for (int j = 0; j < HEIGHT; j++)
	  fv(dnums[k])(j) = elvec(k*HEIGHT+j);
  }
  
  template <> void T_LinearForm<double>::
  SetElementVector (const Array<int> & dnums,
		    const FlatVector<double> & elvec) 
  {
    FlatVector<double> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (dnums[k] != -1)
	fv(dnums[k]) = elvec(k);
  }
  
  template <> void T_LinearForm<Complex>::
  SetElementVector (const Array<int> & dnums,
		    const FlatVector<Complex> & elvec) 
  {
    FlatVector<Complex> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (dnums[k] != -1)
	fv(dnums[k]) = elvec(k);
  }
  




  template <typename TV>
  void T_LinearForm<TV> ::
  GetElementVector (const Array<int> & dnums,
		    FlatVector<TSCAL> & elvec) const
  {
    FlatVector<TV> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (dnums[k] != -1)
	for (int j = 0; j < HEIGHT; j++)
	  elvec(k*HEIGHT+j) = fv(dnums[k])(j);
  }
  
  template <> void T_LinearForm<double>::
  GetElementVector (const Array<int> & dnums,
		    FlatVector<double> & elvec) const
  {
    FlatVector<double> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (dnums[k] != -1)
	elvec(k)= fv(dnums[k]);
  }
  
  template <> void T_LinearForm<Complex>::
  GetElementVector (const Array<int> & dnums,
		    FlatVector<Complex> & elvec) const
  {
    FlatVector<Complex> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (dnums[k] != -1)
	elvec(k)= fv(dnums[k]);
  }
  







  LinearForm * CreateLinearForm (const FESpace * space,
				 const string & name, const Flags & flags)
  {
    /*
      LinearForm * lf;
      CreateVecObject2 (lf, T_LinearForm, 
      space->GetDimension(), space->IsComplex(),   
      *space, name);

      lf->SetIndependent (flags.GetDefineFlag ("independent"));
      return lf;
    */

    LinearForm * lf = 
      CreateVecObject  <T_LinearForm, LinearForm, const FESpace, const string, const Flags>
      (space->GetDimension() * int(flags.GetNumFlag("cacheblocksize",1)), space->IsComplex(), *space, name, flags);
    
    lf->SetIndependent (flags.GetDefineFlag ("independent"));

    if (flags.GetDefineFlag ( "noinitialassembling" )) lf->SetNoInitialAssembling();
    
    lf->SetCacheBlockSize(int(flags.GetNumFlag("cacheblocksize",1)));

    return lf;
  }




  template class S_LinearForm<double>;
  template class S_LinearForm<Complex>;



  template class T_LinearForm<double>;
  template class T_LinearForm<Complex>;

  template class T_LinearForm< Vec<3> >;
}
