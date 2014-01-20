#include <comp.hpp>
#include <parallelngs.hpp>

//#include <unistd.h>



  /*
namespace ngla
{
  extern template class S_BaseVector<double>;
  // template class S_BaseVector<Complex>;
  
  extern template class VFlatVector<double>;
  
  extern template class S_BaseVectorPtr<double>;
  extern template class S_BaseVectorPtr<Complex>;

  extern template class VVector<double>;
  extern template class VVector<Complex>;

  extern template class VVector<Vec<1>>;
  extern template class VVector<Vec<2>>;
  extern template class VVector<Vec<3>>;
  extern template class VVector<Vec<4>>;
  extern template class VVector<Vec<5>>;
  extern template class VVector<Vec<6>>;
  extern template class VVector<Vec<7>>;
  extern template class VVector<Vec<8>>;
  extern template class VVector<Vec<9>>;
  extern template class VVector<Vec<10>>;
  extern template class VVector<Vec<11>>;
  extern template class VVector<Vec<12>>;

  extern template class VVector<Vec<1,Complex>>;
  extern template class VVector<Vec<2,Complex>>;
  extern template class VVector<Vec<3,Complex>>;
  extern template class VVector<Vec<4,Complex>>;
  extern template class VVector<Vec<5,Complex>>;
  extern template class VVector<Vec<6,Complex>>;
  extern template class VVector<Vec<7,Complex>>;
  extern template class VVector<Vec<8,Complex>>;
  extern template class VVector<Vec<9,Complex>>;
  extern template class VVector<Vec<10,Complex>>;
  extern template class VVector<Vec<11,Complex>>;
  extern template class VVector<Vec<12,Complex>>;
}
*/


namespace ngcomp
{

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
    allocated = false;
    initialassembling = true;
    cacheblocksize = 1;
  }

  LinearForm :: ~LinearForm ()
  { 
    for (int i = 0; i < parts.Size(); i++)
      if (parts_deletable[i]) delete parts[i];
  }

  void LinearForm :: AddIntegrator (LinearFormIntegrator * lfi, bool deletable)
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

  template <class TV>
  BaseVector & T_LinearForm<TV> :: GetVector () const 
  { 
    return *vec; 
  }

  template <class SCAL>
  void S_LinearForm<SCAL> :: Assemble (LocalHeap & clh)
  {
    static Timer timer("Vector assembling");
    static Timer timer1("Vector assembling 1", 2);
    static Timer timer2("Vector assembling 2", 2);
    static Timer timer2a("Vector assembling 2a", 2);
    static Timer timer2b("Vector assembling 2b", 2);
    static Timer timer2c("Vector assembling 2c", 2);
    static Timer timer3("Vector assembling 3", 2);
    RegionTimer reg (timer);

    assembled = true;
    
    if (independent)
      {
	AssembleIndependent(clh);
	return;
      }

    try
      {
	ma.PushStatus ("Assemble Vector");

	timer1.Start();

        // check if integrators fit to space
        for (int i = 0; i < NumIntegrators(); i++)
          if (parts[i] -> BoundaryForm())
            {
              for (int j = 0; j < ma.GetNSE(); j++)
                {
                  HeapReset hr(clh);
		  if (!parts[i] -> SkeletonForm()){
		    if (parts[i] -> DefinedOn (ma.GetSElIndex(j)))
		      parts[i] -> CheckElement (fespace.GetSFE(j, clh));
		  }else{
		    //TODO: get aligned volelement and check that
		  }
                }
            }
          else
            {
              if (parts[i] -> SkeletonForm()) 
                { 
                  throw Exception ("There are no LinearformIntegrator which act on the skeleton so far!");
                }
              else
                {
                  for (int j = 0; j < ma.GetNE(); j++)
                    {
                      HeapReset hr(clh);
                      if (parts[i] -> DefinedOn (ma.GetElIndex(j)))
                        parts[i] -> CheckElement (fespace.GetFE(j, clh));
                    }
                }
            }

	if(!allocated || ( this->GetVector().Size() != this->fespace.GetNDof()))
	  {
	    AllocateVector();
	    allocated=true;
	  }
	else
	  {
	    this->GetVector() = TSCAL(0);
	  }

	int ne = ma.GetNE();
	int nf = ma.GetNFacets();
	int nse = ma.GetNSE();

	bool hasbound = false;
	bool hasinner = false;
	bool hasskeletonbound = false;
	bool hasskeletoninner = false;
	
	//check hasbound, hasinner, hasskeletonbound, hasskeletoninner
	for (int j = 0; j < parts.Size(); j++)
	  {
	    if (!parts[j] -> SkeletonForm()){
	      if (parts[j] -> BoundaryForm())
		hasbound = true;
	      else 
		if (!parts[j] -> IntegrationAlongCurve()) //CL: ist das richtig? 
		  hasinner = true;
	    }else{
	      if (parts[j] -> BoundaryForm())
		hasskeletonbound = true;
	      else 
		if (!parts[j] -> IntegrationAlongCurve()) //CL: ist das richtig? 
		  hasskeletoninner = true;
	    }
	  }

	if (print)
	  {
	    *testout << " LINEARFORM TEST:" << endl;
	    *testout << " hasinner = " << hasinner << endl;
	    *testout << " hasouter = " << hasbound << endl;
	    *testout << " hasskeletoninner = " << hasskeletoninner << endl;
	    *testout << " hasskeletonouter = " << hasskeletonbound << endl;
	  }

	int nrcases = 0;
	int loopsteps = 0;
	if (hasinner) {nrcases++; loopsteps+=ne;}
	if (hasbound) {nrcases++; loopsteps+=nse;}
	if (hasskeletoninner) {nrcases++; loopsteps+=nf;}
	if (hasskeletonbound) {nrcases++; loopsteps+=nse;}
	if (fespace.specialelements.Size()>0) {nrcases++; loopsteps+=fespace.specialelements.Size();}
	// int actcase = 0;
	int gcnt = 0; //global count (for all cases)
	
	timer1.Stop();


	if (hasinner)
	  {

	    ProgressOutput progress (ma, "assemble element", ma.GetNE());
	    gcnt += ne;
            
            
	    IterateElements 
              (fespace, VOL, clh, 
               [&] (ElementId ei, LocalHeap & lh)
               {
                 RegionTimer reg2 (timer2);
                 timer2a.Start();
                 progress.Update ();

                 const FiniteElement & fel = fespace.GetFE(ei, lh);
                 const ElementTransformation & eltrans = ma.GetTrafo (ei, lh);
                 // FlatArray<int> dnums = fespace.GetDofNrs (ei, lh);
                 Array<int> dnums (fel.GetNDof(), lh);
                 fespace.GetDofNrs (ei, dnums);
                 timer2a.Stop();

                 for (int j = 0; j < parts.Size(); j++)
                   {
                     if (!parts[j] -> VolumeForm()) continue;
                     if (!parts[j] -> DefinedOn (eltrans.GetElementIndex())) continue;
                     int elvec_size = fel.GetNDof()*fespace.GetDimension();
                     FlatVector<TSCAL> elvec(elvec_size, lh);

		     timer2b.Start();
                     parts[j] -> CalcElementVector (fel, eltrans, elvec, lh);
                     timer2b.Stop();
                     
                     if (printelvec)
                       {
                         testout->precision(8);
                         *testout << "elnum = " << ei.Nr() << endl
                                  << "integrator " << parts[j]->Name() << endl
                                  << "dnums = " << endl << dnums << endl
                                  << "element-index = " << eltrans.GetElementIndex() << endl
                                  << "elvec = " << endl << elvec << endl;
                       }
                     timer2c.Start();
                     fespace.TransformVec (ei, elvec, TRANSFORM_RHS);
                     AddElementVector (dnums, elvec, parts[j]->CacheComp()-1);
                     timer2c.Stop();
                   }
               });
          }
            

	RegionTimer reg3(timer3);

	if (hasbound)
	  {
	    ProgressOutput progress (ma, "assemble surface element", ma.GetNSE());
	    gcnt += nse;


	    IterateElements 
              (fespace, BND, clh, 
               [&] (ElementId ei, LocalHeap & lh)
               {
                 RegionTimer reg2 (timer2);
                 
                 progress.Update ();

                 const FiniteElement & fel = fespace.GetFE(ei, lh);
                 ElementTransformation & eltrans = ma.GetTrafo (ei, lh);
                 // FlatArray<int> dnums = fespace.GetDofNrs (ei, lh);
                 Array<int> dnums(fel.GetNDof(), lh);
                 fespace.GetDofNrs (ei, dnums);

                 for (int j = 0; j < parts.Size(); j++)
                   {
                     if (!parts[j] -> BoundaryForm()) continue;
                     if (!parts[j] -> DefinedOn (eltrans.GetElementIndex())) continue;
		     
                     int elvec_size = fel.GetNDof()*fespace.GetDimension();
                     FlatVector<TSCAL> elvec(elvec_size, lh);
                     
                     parts[j] -> CalcElementVector (fel, eltrans, elvec, lh);
                     
                     if (printelvec)
                       {
                         testout->precision(8);
                         *testout << "surface-elnum = " << ei.Nr() << endl
                                  << "integrator " << parts[j]->Name() << endl
                                  << "dnums = " << endl << dnums << endl
                                  << "element-index = " << eltrans.GetElementIndex() << endl
                                  << "elvec = " << endl << elvec << endl;
                       }
                     
                     fespace.TransformVec (ei, elvec, TRANSFORM_RHS);
                     AddElementVector (dnums, elvec, parts[j]->CacheComp()-1);
                   }
               });

	  }//end of hasbound
	






	if(hasskeletoninner)
	{
	  cout << "\rInnerFacetIntegrators not known - cannot handle it yet" << endl;	  
	}
	if(hasskeletonbound)
	{
#pragma omp parallel
	  {
	    LocalHeap lh = clh.Split();
	    Array<int> dnums;
	    Array<int> fnums, elnums, vnums;
	    //Schleife fuer Facet-Integrators: 
#pragma omp for	      
	      for (int i = 0; i < nse; i++)
		{
#pragma omp critical (linformsurfneighprint)		    
		  {
		    gcnt++;
		    if (i % 10 == 0)
		      cout << "\rassemble facet surface element " << i << "/" << nse << flush;
		    ma.SetThreadPercentage ( 100.0*(gcnt) / (loopsteps) );
		  }

		  HeapReset hr(lh);
		  
		  
		  ma.GetSElFacets(i,fnums);
		  int fac = fnums[0];
		  ma.GetFacetElements(fac,elnums);
		  int el = elnums[0];
		  ma.GetElFacets(el,fnums);
		  int facnr = 0;
		  for (int k=0; k<fnums.Size(); k++)
		    if(fac==fnums[k]) facnr = k;
		  
		  const FiniteElement & fel = fespace.GetFE (el, lh);
		
		  ElementTransformation & eltrans = ma.GetTrafo (el, false, lh);
		  ElementTransformation & seltrans = ma.GetTrafo (i, true, lh);

		  fespace.GetDofNrs (el, dnums);
		  ma.GetElVertices (el, vnums);		
	      
		  for (int j = 0; j < parts.Size(); j++)
		    {
		      if (!parts[j] -> SkeletonForm()) continue;
		      if (!parts[j] -> BoundaryForm()) continue;
		      if (!parts[j] -> DefinedOn (ma.GetSElIndex (i))) continue;
		      if (parts[j] -> IntegrationAlongCurve()) continue;		    
		  
		      int elvec_size = dnums.Size()*fespace.GetDimension();
		      FlatVector<TSCAL> elvec(elvec_size, lh);
		      dynamic_cast<const FacetLinearFormIntegrator*>(parts[j]) 
			  -> CalcFacetVector (fel,facnr,eltrans,vnums,seltrans, elvec, lh);
		      if (printelvec)
			{
			  testout->precision(8);

			  (*testout) << "surface-elnum= " << i << endl;
			  (*testout) << "integrator " << parts[j]->Name() << endl;
			  (*testout) << "dnums = " << endl << dnums << endl;
			  (*testout) << "(vol)element-index = " << eltrans.GetElementIndex() << endl;
			  (*testout) << "elvec = " << endl << elvec << endl;
			}


		      fespace.TransformVec (i, false, elvec, TRANSFORM_RHS);
		    
#pragma omp critical (addelvec3)		    
		      {
			AddElementVector (dnums, elvec, parts[j]->CacheComp()-1);
		      }
		    }
		}
		
	  }//end of parallel
	  cout << "\rassemble facet surface element  " << nse << "/" << nse << endl;	  
	}//endof hasskeletonbound


	Array<int> dnums;

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
			fespace.GetDofNrs(element,dnums);
		      }
		    ElementTransformation & eltrans = ma.GetTrafo (element, false, clh);

		    
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
			MappedIntegrationPoint<1,3> s_sip(ip,eltrans);
			MappedIntegrationPoint<3,3> g_sip(ip,eltrans);
			Vec<3> tv;
			tv(0) = tangent(0); tv(1) = tangent(1); tv(2) = tangent(2);
			s_sip.SetTV(tv);
			parts[j]->CalcElementVectorIndependent(*fel,
								   s_sip,
								   g_sip,
								   elvec,clh,true);
		      }
		    else if (eltrans.SpaceDim() == 2)
		      {
			MappedIntegrationPoint<1,2> s_sip(ip,eltrans);
			MappedIntegrationPoint<2,2> g_sip(ip,eltrans);
			Vec<2> tv;
			tv(0) = tangent(0); tv(1) = tangent(1);
			s_sip.SetTV(tv);
			parts[j]->CalcElementVectorIndependent(*fel,
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




  void LinearForm :: AddElementVector (FlatArray<int> dnums,
				       FlatVector<double> elvec,
				       int cachecomp)
  {
    throw Exception ("LinearForm::AddElementVector: real elvec for complex li-form");
  }
  void LinearForm :: SetElementVector (FlatArray<int> dnums,
				       FlatVector<double> elvec)
  {
    throw Exception ("LinearForm::SetElementVector: real elvec for complex li-form");
  }

  void LinearForm :: GetElementVector (FlatArray<int> dnums,
                                       FlatVector<double> elvec) const
  {
    throw Exception ("LinearForm::GetElementVector: real elvec for complex li-form");
  }

  
  void LinearForm :: AddElementVector (FlatArray<int> dnums,
				       FlatVector<Complex> elvec,
				       int cachecomp)
  {
    throw Exception ("LinearForm::AddElementVector: complex elvec for real li-form");
  }

  void LinearForm :: SetElementVector (FlatArray<int> dnums,
				       FlatVector<Complex> elvec)
  {
    throw Exception ("LinearForm::SetElementVector: complex elvec for real li-form");
  }

  void LinearForm :: GetElementVector (FlatArray<int> dnums,
				       FlatVector<Complex> elvec) const
  {
    throw Exception ("LinearForm::GetElementVector: complex elvec for real li-form");
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
	// ElementTransformation seltrans, geltrans;

	for (int i = 0; i < nse; i++)
	  {
	    lh.CleanUp();
	    
	    const FiniteElement & sfel = fespace.GetSFE (i, lh);
	    // ma.GetSurfaceElementTransformation (i, seltrans);
	    ElementTransformation & seltrans = ma.GetTrafo (i, true, lh);

	      	
	    // (*testout) << "el = " << i << ", ind = " << ma.GetSElIndex(i) << endl;
	    if (!parts[0]->DefinedOn (ma.GetSElIndex(i))) continue;
	    // (*testout) << "integrate surf el " << endl;
	    
	    const IntegrationRule & ir = SelectIntegrationRule (sfel.ElementType(), 5);
	    
	    for (int j = 0; j < ir.GetNIP(); j++)
	      {
		const IntegrationPoint & ip = ir[j];
		MappedIntegrationPoint<2,3> sip(ip, seltrans);
		
		// (*testout) << "point = " << sip.GetPoint() << endl;
		
		IntegrationPoint gip;
		int elnr;
		// elnr = ma.FindElementOfPoint (FlatVector<>(sip.GetPoint()), gip, 1);
                elnr = ma.FindElementOfPoint (sip.GetPoint(), gip, 1);
		
		// (*testout) << "elnr = " << elnr << endl;
		if (elnr == -1) continue;
		
		const FiniteElement & gfel = fespace.GetFE (elnr, lh);
		// ma.GetElementTransformation (elnr, geltrans);
		ElementTransformation & geltrans = ma.GetTrafo (elnr, false, lh);

		MappedIntegrationPoint<3,3> gsip(gip, geltrans);
		
		// (*testout) << " =?= p = " << gsip.GetPoint() << endl;

		fespace.GetDofNrs (elnr, dnums);
		
		for (int k = 0; k < parts.Size(); k++)
		  {
		    FlatVector<TSCAL> elvec;
		    
		    if(geltrans.SpaceDim() == 3)
		      {
			MappedIntegrationPoint<2,3> s_sip(ip,seltrans);
			MappedIntegrationPoint<3,3> g_sip(gip,geltrans);
			parts[k] -> CalcElementVectorIndependent
			  (gfel,s_sip,g_sip,elvec,lh);
		      }
		    else if(geltrans.SpaceDim() == 2)
		      {
			MappedIntegrationPoint<1,2> s_sip(ip,seltrans);
			MappedIntegrationPoint<2,2> g_sip(gip,geltrans);
			parts[k] -> CalcElementVectorIndependent
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


  template <typename TV>
  void T_LinearForm<TV> :: AllocateVector ()
  {
    delete vec;
    const FESpace & fes = this->fespace;
#ifdef PARALLEL
    if ( &fes.GetParallelDofs() )
      vec = new ParallelVVector<TV> (fes.GetNDof(), &fes.GetParallelDofs(), DISTRIBUTED);
    else
#endif
      vec = new VVector<TV> (fes.GetNDof());
    (*vec) = TSCAL(0);
    vec->SetParallelStatus (DISTRIBUTED);
  }





  template <typename TV>
  void T_LinearForm<TV> :: CleanUpLevel ()
  {
    delete vec;
    vec = 0;
    this -> allocated = false;
  }


  template <typename TV>
  void T_LinearForm<TV> ::
  AddElementVector (FlatArray<int> dnums, FlatVector<TSCAL> elvec,
		    int cachecomp) 
  {
    FlatVector<TV> fv = vec->FV();
    if(cachecomp < 0)
      {
        Scalar2ElemVector<TV, TSCAL> ev(elvec);
	for (int k = 0; k < dnums.Size(); k++)
	  if (dnums[k] != -1)
	    fv(dnums[k]) += ev(k);
      }
    else
      {
	for (int k = 0; k < dnums.Size(); k++)
	  if (dnums[k] != -1)
	    fv(dnums[k])(cachecomp) += elvec(k);
      }
  }
  
  template <> void T_LinearForm<double>::
  AddElementVector (FlatArray<int> dnums,
		    FlatVector<double> elvec,
		    int cachecomp) 
  {
    vec -> AddIndirect (dnums, elvec);
  }
  
  template <> void T_LinearForm<Complex>::
  AddElementVector (FlatArray<int> dnums,
		    FlatVector<Complex> elvec,
		    int cachecomp) 
  {
    vec -> AddIndirect (dnums, elvec);
  }

  template <typename TV>
  void T_LinearForm<TV> ::
  SetElementVector (FlatArray<int> dnums,
		    FlatVector<TSCAL> elvec) 
  {
    FlatVector<TV> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (dnums[k] != -1)
	for (int j = 0; j < HEIGHT; j++)
	  fv(dnums[k])(j) = elvec(k*HEIGHT+j);
  }
  
  template <> void T_LinearForm<double>::
  SetElementVector (FlatArray<int> dnums,
		    FlatVector<double> elvec) 
  {
    vec -> SetIndirect (dnums, elvec);
  }
  
  template <> void T_LinearForm<Complex>::
  SetElementVector (FlatArray<int> dnums,
		    FlatVector<Complex> elvec) 
  {
    vec -> SetIndirect (dnums, elvec);
  }
  


  template <typename TV>
  void T_LinearForm<TV> ::
  GetElementVector (FlatArray<int> dnums,
		    FlatVector<TSCAL> elvec) const
  {
    FlatVector<TV> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (dnums[k] != -1)
	for (int j = 0; j < HEIGHT; j++)
	  elvec(k*HEIGHT+j) = fv(dnums[k])(j);
  }
  
  template <> void T_LinearForm<double>::
  GetElementVector (FlatArray<int> dnums,
		    FlatVector<double> elvec) const
  {
    vec -> GetIndirect (dnums, elvec);
  }
  
  template <> void T_LinearForm<Complex>::
  GetElementVector (FlatArray<int> dnums,
		    FlatVector<Complex> elvec) const
  {
    vec -> GetIndirect (dnums, elvec);
  }
  







  LinearForm * CreateLinearForm (const FESpace * space,
				 const string & name, const Flags & flags)
  {
    LinearForm * lf = 
      CreateVecObject  <T_LinearForm, LinearForm> 
      (space->GetDimension() * int(flags.GetNumFlag("cacheblocksize",1)), 
       space->IsComplex(), *space, name, flags);
    
    lf->SetIndependent (flags.GetDefineFlag ("independent"));
    if (flags.GetDefineFlag ( "noinitialassembling" )) lf->SetNoInitialAssembling();
    lf->SetCacheBlockSize(int(flags.GetNumFlag("cacheblocksize",1)));

    return lf;
  }



  template class S_LinearForm<double>;
  template class S_LinearForm<Complex>;

  template class T_LinearForm<double>;
  template class T_LinearForm<Complex>;

  // template class T_LinearForm< Vec<3> >;
}
