#include <comp.hpp>
#include <parallelngs.hpp>


namespace ngcomp
{

  LinearForm ::
  LinearForm (shared_ptr<FESpace> afespace, 
	      const string & aname,
	      const Flags & flags)
    : NGS_Object(afespace->GetMeshAccess(), flags, aname), fespace(afespace)
  {
    independent = false;
    print = flags.GetDefineFlag ("print");
    printelvec = flags.GetDefineFlag ("printelvec");
    assembled = false;
    allocated = false;
    initialassembling = true;
    checksum = flags.GetDefineFlag ("checksum");
    cacheblocksize = 1;
  }

  LinearForm & LinearForm :: AddIntegrator (shared_ptr<LinearFormIntegrator> lfi)
  {
    auto anydim = dynamic_pointer_cast<LinearFormIntegratorAnyDim> (lfi);
    if (anydim) lfi = anydim->GetLFI (ma->GetDimension());

    parts.Append (lfi);
    if (!lfi->SkeletonForm() && !lfi->IntegrationAlongCurve())
      VB_parts[lfi->VB()].Append(lfi);
    return *this;
  }
  
  void LinearForm :: PrintReport (ostream & ost) const
  {
    ost << "on space " << GetFESpace()->GetName() << endl
	<< "integrators: " << endl;
    for (int i = 0; i < parts.Size(); i++)
      ost << "  " << parts[i]->Name() << endl;
  }

  Array<MemoryUsage> LinearForm :: GetMemoryUsage () const
  {
    if (GetVectorPtr())  
      {
	auto mu = GetVectorPtr()->GetMemoryUsage ();
	for (int i = 0; i < mu.Size(); i++)
	  mu[i].AddName (string(" lf ")+GetName());
        return mu;
      }
    return Array<MemoryUsage>();
  }

  template <class SCAL>
  void S_LinearForm<SCAL> :: Assemble (LocalHeap & clh)
  {
    static Timer timer("Vector assembling");
    static Timer timer1("Vector assembling 1", NoTracing);
    static Timer timer2("Vector assembling 2", NoTracing);
    static Timer timer3("Vector assembling 3", NoTracing);
    static mutex linformsurfneighprint_mutex;
    static mutex addelvec3_mutex;
    RegionTimer reg (timer);

    assembled = true;
    
    if (independent)
      {
	AssembleIndependent(clh);
	return;
      }

    try
      {
	ma->PushStatus ("Assemble Vector");
	timer1.Start();

        // check if integrators fit to space
        /*
        for (int i = 0; i < NumIntegrators(); i++)
          if (parts[i] -> BoundaryForm())
            {
              for (int j = 0; j < ma->GetNSE(); j++)
                {
                  HeapReset hr(clh);
		  if (!parts[i] -> SkeletonForm()){
		    if (parts[i] -> DefinedOn (ma->GetSElIndex(j)))
		      parts[i] -> CheckElement (fespace->GetSFE(j, clh));
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
                  for (int j = 0; j < ma->GetNE(); j++)
                    {
                      HeapReset hr(clh);
                      if (parts[i] -> DefinedOn (ma->GetElIndex(j)))
                        parts[i] -> CheckElement (fespace->GetFE(j, clh));
                    }
                }
            }
        */

	for (VorB vb : {VOL,BND,BBND,BBBND})
	  {
            if (int(vb) > ma->GetDimension()) continue;
	    for (auto lfi : parts)
	      {
		if(lfi->VB()==vb)
		  {
		    if(lfi->SkeletonForm())
		      {
			if(vb==VOL)
			  throw Exception ("There are no LinearformIntegrator which act on the skeleton so far!");
			else {
		      //TODO: get aligned volelement and check that
			}
		      }
		    IterateElements
		      (*fespace,vb,clh, [&](FESpace::Element el, LocalHeap & lh)
		       {
                         if (lfi->DefinedOn(el.GetIndex()) && lfi->DefinedOnElement(el.Nr()))
			   lfi->CheckElement(el.GetFE());
		       });
		    
		  }
	      }
	  }

	if(!allocated || ( this->GetVector().Size() != this->fespace->GetNDof()))
	  {
	    AllocateVector();
	    allocated=true;
	  }
	else
	  {
	    this->GetVector() = TSCAL(0);
	    this->GetVector().SetParallelStatus(DISTRIBUTED);
	  }


	bool hasparts[] = {false,false,false,false};
	bool hasskeletonparts[] = {false,false};

	for(int j = 0; j < parts.Size(); j++)
	  {
	    if(!parts[j]->SkeletonForm())
	      {
		hasparts[parts[j]->VB()] = true;
	      }
	    else
	      {
		if(parts[j]->VB()==BBND)
		  throw Exception("Skeleton not implemented for co dim 2!");
		hasskeletonparts[parts[j]->VB()] = true;
	      }
	  }

	// ended here...
        
	if (print)
	  {
	    *testout << " LINEARFORM TEST:" << endl;
	    *testout << " hasinner = " << hasparts[VOL] << endl;
	    *testout << " hasouter = " << hasparts[BND] << endl;
	    *testout << " hasbbound = " << hasparts[BBND] << endl;
	    *testout << " hasskeletoninner = " << hasskeletonparts[VOL] << endl;
	    *testout << " hasskeletonouter = " << hasskeletonparts[BND] << endl;
	  }

	int loopsteps = 0;
	for(VorB vb : {VOL,BND,BBND,BBBND})
	  {
	    if(hasparts[vb]){
	      int ne = ma->GetNE(vb);
	      loopsteps += ne;
	    }
	  }
	if(hasskeletonparts[VOL])
	  loopsteps += ma->GetNFacets();
	if(hasskeletonparts[BND])
	  loopsteps += ma->GetNSE();
	// if(fespace->specialelements.Size())
        // loopsteps += fespace->specialelements.Size();
        
	// int actcase = 0;
	int gcnt = 0; //global count (for all cases)
	
	timer1.Stop();

	for(VorB vb : {VOL,BND,BBND,BBBND})
	  {
	    if(hasparts[vb])
	      {
		int ne = ma->GetNE(vb);
		// string vb_str = vb==VOL ? "VOL" : (vb==BND ? "BND" : "BBND");
		// ProgressOutput progress (ma, string("assemble ") + vb_str + string(" element"),ne);
                ProgressOutput progress (ma, string("assemble ") + ToString(vb) + string(" element"),ne);
		gcnt += ne;
		IterateElements
		  (*fespace,vb,clh,[&] (FESpace::Element el, LocalHeap &lh)
		   {
		     // RegionTimer reg2(timer2);
		     progress.Update();

		     auto & fel = el.GetFE();
		     auto & eltrans = el.GetTrafo();

		     // for(int j = 0; j<parts.Size(); j++)
                     for (auto & lfip : VB_parts[vb])
		       {
			 // if(parts[j]->VB() != vb) continue;
			 // if(!parts[j]->DefinedOn(el.GetIndex())) continue;
                         // if (parts[j]->SkeletonForm()) continue;
			 // if (parts[j] -> IntegrationAlongCurve()) continue;
                         if(!lfip->DefinedOn(el.GetIndex())) continue;
                         if(!lfip->DefinedOnElement(el.Nr())) continue;
                         
			 int elvec_size = fel.GetNDof()*fespace->GetDimension();
			 FlatVector<TSCAL> elvec(elvec_size, lh);
			 auto & mapped_trafo = eltrans.AddDeformation(lfip->GetDeformation().get(), lh);
			 lfip -> CalcElementVector (fel, mapped_trafo, elvec, lh);
			 
			 if (printelvec)
			   {
			     testout->precision(8);
			     *testout << "elnum = " << el.Nr() << endl
				      << "integrator " << lfip->Name() << endl
				      << "dnums = " << endl << el.GetDofs() << endl
				      << "element-index = " << eltrans.GetElementIndex() << endl
				      << "elvec = " << endl << elvec << endl;
			   }
			 
			 fespace->TransformVec (el, elvec, TRANSFORM_RHS);
			 AddElementVector (el.GetDofs(), elvec, lfip->CacheComp()-1);
		       }
		   });
	      }
	  }

        
        for (auto pe : pnteval)
          {
            LocalHeap & lh = clh;
            HeapReset hr(lh);

            IntegrationPoint ip;
            auto elnr = ma->FindElementOfPoint(pe->point, ip, false);

            ElementId ei(VOL, elnr);
            auto & trafo = ma->GetTrafo(ei, lh);
            auto & fel = fespace->GetFE(ei, lh);
            Array<int> dnums(fel.GetNDof(), lh);
            fespace->GetDofNrs(ei, dnums);

            IntegrationRule ir (1, &ip);
            auto & mir = trafo(ir, lh);

            pe->cf->TraverseTree
              ([&] (CoefficientFunction& cf)
               {
                 if (auto * proxy = dynamic_cast<ProxyFunction*> (&cf))
                   {
                     FlatVector<SCAL> elvec(dnums.Size()*fespace->GetDimension(), lh);
                     FlatMatrix<SCAL> values(1, 1, lh);
                     ProxyUserData ud;
                     const_cast<ElementTransformation&>(trafo).userdata = &ud;
        
                     elvec = 0;

                     FlatMatrix<SCAL> proxyvalues(1, proxy->Dimension(), lh);
                     for (int k = 0; k < proxy->Dimension(); k++)
                       {
                         ud.testfunction = proxy;
                         ud.test_comp = k;
                         
                         pe->cf -> Evaluate (mir, values);
                         proxyvalues(0,k) = values(0,0);
                       }
                     proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec, lh);
                     AddElementVector (dnums, elvec);
                   }
               });
          }

	if(hasskeletonparts[VOL])
	{
	  cout << "\rInnerFacetIntegrators not known - cannot handle it yet" << endl;	  
	}
	if(hasskeletonparts[BND])
	{
          int nse = ma->GetNE(BND);
	  ParallelForRange( IntRange(nse), [&] ( IntRange r )
	  {
	    LocalHeap lh = clh.Split();
	    Array<int> dnums;
	    Array<int> fnums, elnums, vnums;
	    // loop for facet integrators: 
            for (int i : r)
		{
		  {
                    lock_guard<mutex> guard(linformsurfneighprint_mutex);
		    gcnt++;
		    if (i % 10 == 0)
		      cout << IM(3) << "\rassemble facet surface element " << i << "/" << nse << flush;
		    ma->SetThreadPercentage ( 100.0*(gcnt) / (loopsteps) );
		  }

		  HeapReset hr(lh);
		  
                  ElementId sei(BND, i);
		  
		  fnums = ma->GetElFacets(sei);
		  int fac = fnums[0];
		  ma->GetFacetElements(fac,elnums);
		  int el = elnums[0];
                  ElementId ei(VOL, el);
		  fnums = ma->GetElFacets(ei);
		  int facnr = 0;
		  for (int k=0; k<fnums.Size(); k++)
		    if(fac==fnums[k]) facnr = k;

		  const FiniteElement & fel = fespace->GetFE (ei, lh);
		
		  ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
		  ElementTransformation & seltrans = ma->GetTrafo (sei, lh);

		  fespace->GetDofNrs (ei, dnums);
		  vnums = ma->GetElVertices (ei);		
	      
		  for (int j = 0; j < parts.Size(); j++)
		    {
		      if (!parts[j] -> SkeletonForm()) continue;
		      if (parts[j] -> VB()!=BND) continue;
		      if (!parts[j] -> DefinedOn (ma->GetElIndex (sei))) continue;
		      if (!parts[j] -> DefinedOnElement (i)) continue;
		      if (parts[j] -> IntegrationAlongCurve()) continue;		    
		  
		      int elvec_size = dnums.Size()*fespace->GetDimension();
		      FlatVector<TSCAL> elvec(elvec_size, lh);
		      dynamic_cast<const FacetLinearFormIntegrator*>(parts[j].get()) 
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

		      fespace->TransformVec (ei, elvec, TRANSFORM_RHS);
		    
		      {
                        lock_guard<mutex> guard(addelvec3_mutex);
			AddElementVector (dnums, elvec, parts[j]->CacheComp()-1);
		      }
		    }
		}
		
	  });//end of parallel
	  cout << IM(3) << "\rassemble facet surface element  " << nse << "/" << nse << endl;	  
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
		for(int i=0; i<ma->GetNDomains(); i++)
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
			cout << IM(3) << "\rassemble curvepoint " << i << "/" << parts[j]->NumCurvePoints() << flush;
			ma->SetThreadPercentage(100.*i/parts[j]->NumCurvePoints());
		      }
		    
		    FlatVector<TSCAL> elvec;
		    IntegrationPoint ip;
		    if(domains.Size() > 0)
		      element = ma->FindElementOfPoint(parts[j]->CurvePoint(i),ip,true,&domains);
		    else
		      element = ma->FindElementOfPoint(parts[j]->CurvePoint(i),ip,true);
		    if(element < 0)
		      throw Exception("element for curvepoint not found");
		    
		    if(element != oldelement)
		      { 
			clh.CleanUp();
			fel = &fespace->GetFE(ElementId(VOL, element),clh);
			fespace->GetDofNrs(ElementId(VOL,element),dnums);
		      }
		    ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL, element), clh);

		    
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
			MappedIntegrationPoint<1,3> s_sip(ip,eltrans, -42 /* Don't call CalcPointJacobian (eltrans expects DIMR==DIMS ) */);
			MappedIntegrationPoint<3,3> g_sip(ip,eltrans);
			s_sip.Point() = g_sip.Point();
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
			MappedIntegrationPoint<1,2> s_sip(ip,eltrans, -42 /* Don't call CalcPointJacobian (eltrans expects DIMR==DIMS ) */);
			MappedIntegrationPoint<2,2> g_sip(ip,eltrans);
			s_sip.Point() = g_sip.Point();
			Vec<2> tv;
			tv(0) = tangent(0); tv(1) = tangent(1);
			s_sip.SetTV(tv);
			parts[j]->CalcElementVectorIndependent(*fel,
								   s_sip,
								   g_sip,
								   elvec,clh,true);
		      }
		    fespace->TransformVec (ElementId(VOL, element), elvec, TRANSFORM_RHS);

		    elvec *= (oldlength+length); // Des is richtig
		    
		    
		    AddElementVector (dnums, elvec, parts[j]->CacheComp()-1);
		    
		    
		    
		    
		    clh.CleanUp(heapp);
		    
		    oldlength = length;      
		    
		    oldelement = element;
		  }
	      }



	    cout << IM(3) << "\rassemble curvepoint " << parts[j]->NumCurvePoints() << "/" 
		 << parts[j]->NumCurvePoints() << endl;
	    clh.CleanUp();

	  }
	
	
	if (print)
	  {
	    (*testout) << "Linearform " << GetName() << ": " << endl;
	    (*testout) << GetVector() << endl;
	  }


        if (checksum)
          cout << "|vector| = " 
               << setprecision(16) << L2Norm (GetVector()) << endl;


	ma->PopStatus ();
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

	// int ne = ma->GetNE();
	int nse = ma->GetNSE();
	
	
	Array<int> dnums;
	// ElementTransformation seltrans, geltrans;

	for (int i = 0; i < nse; i++)
	  {
	    lh.CleanUp();
	    ElementId sei(BND, i);
	    const FiniteElement & sfel = fespace->GetFE (sei, lh);
	    ElementTransformation & seltrans = ma->GetTrafo (sei, lh);

	      	
	    // (*testout) << "el = " << i << ", ind = " << ma->GetSElIndex(i) << endl;
	    if (!parts[0]->DefinedOn (ma->GetElIndex(sei))) continue;
	    if (!parts[0]->DefinedOnElement (i)) continue;
	    // (*testout) << "integrate surf el " << endl;
	    
	    const IntegrationRule & ir = SelectIntegrationRule (sfel.ElementType(), 5);
	    
	    for (int j = 0; j < ir.GetNIP(); j++)
	      {
		const IntegrationPoint & ip = ir[j];
		MappedIntegrationPoint<2,3> sip(ip, seltrans);
		
		// (*testout) << "point = " << sip.GetPoint() << endl;
		
		IntegrationPoint gip;
		int elnr;
		// elnr = ma->FindElementOfPoint (FlatVector<>(sip.GetPoint()), gip, 1);
                elnr = ma->FindElementOfPoint (sip.GetPoint(), gip, 1);
		ElementId ei(VOL, elnr);
		// (*testout) << "elnr = " << elnr << endl;
		if (elnr == -1) continue;
		
		const FiniteElement & gfel = fespace->GetFE (ei, lh);
		// ma->GetElementTransformation (elnr, geltrans);
		ElementTransformation & geltrans = ma->GetTrafo (ei, lh);

		MappedIntegrationPoint<3,3> gsip(gip, geltrans);
		
		// (*testout) << " =?= p = " << gsip.GetPoint() << endl;

		fespace->GetDofNrs (ElementId(VOL,elnr), dnums);
		
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
		    fespace->TransformVec (ElementId(VOL,elnr), elvec, TRANSFORM_RHS);

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
  T_LinearForm<TV> :: ~T_LinearForm () { ; }


  template <typename TV>
  void T_LinearForm<TV> :: AllocateVector ()
  {
    // delete vec;
    // using this::fespace;
    auto fes = this->fespace;
    if ( fes->IsParallel() )
      vec = make_shared<ParallelVVector<TV>> (fes->GetNDof(), 
                                              fes->GetParallelDofs(), DISTRIBUTED);
    else
      vec = make_shared<VVector<TV>> (fes->GetNDof());
    (*vec) = TSCAL(0);
    vec->SetParallelStatus (DISTRIBUTED);
  }



  template <typename TV>
  void T_LinearForm<TV> :: CleanUpLevel ()
  {
    vec.reset();
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
	  if (IsRegularDof(dnums[k]))
	    fv(dnums[k]) += ev(k);
      }
    else
      {
	for (int k = 0; k < dnums.Size(); k++)
	  if (IsRegularDof(dnums[k]))
	    fv(dnums[k])(cachecomp) += elvec(k);
      }
  }
  
  template <> void T_LinearForm<double>::
  AddElementVector (FlatArray<int> dnums,
		    FlatVector<double> elvec,
		    int cachecomp) 
  {
    vec -> AddIndirect (dnums, elvec, fespace->HasAtomicDofs());
  }
  
  template <> void T_LinearForm<Complex>::
  AddElementVector (FlatArray<int> dnums,
		    FlatVector<Complex> elvec,
		    int cachecomp) 
  {
    vec -> AddIndirect (dnums, elvec, fespace->HasAtomicDofs());
  }

  template <typename TV>
  void T_LinearForm<TV> ::
  SetElementVector (FlatArray<int> dnums,
		    FlatVector<TSCAL> elvec) 
  {
    FlatVector<TV> fv = vec->FV();
    for (int k = 0; k < dnums.Size(); k++)
      if (IsRegularDof(dnums[k]))
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
      if (IsRegularDof(dnums[k]))
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
  


  ComponentLinearForm :: ComponentLinearForm (shared_ptr<LinearForm> abase_lf, int acomp, int ancomp)
    : LinearForm( (*dynamic_pointer_cast<CompoundFESpace> (abase_lf->GetFESpace()))[acomp], "comp-lf", Flags()), 
      base_lf(abase_lf), comp(acomp) // , ncomp(ancomp)
  { 
    ;
  }

  LinearForm & ComponentLinearForm :: AddIntegrator (shared_ptr<LinearFormIntegrator> lfi)
  {
    auto block_lfi = make_shared<CompoundLinearFormIntegrator> (lfi, comp);
    base_lf -> AddIntegrator (block_lfi);
    return *this;
  }




  shared_ptr<LinearForm> CreateLinearForm (shared_ptr<FESpace> space,
                                           const string & name, const Flags & flags)
  {
    shared_ptr<LinearForm> lf = 
      CreateSharedVecObject  <T_LinearForm, LinearForm> 
      (space->GetDimension() * int(flags.GetNumFlag("cacheblocksize",1)), 
       space->IsComplex(), space, name, flags);
  
    // shared_ptr<LinearForm> lf(lfp);
    
    lf->SetIndependent (flags.GetDefineFlag ("independent"));
    if (flags.GetDefineFlag ( "noinitialassembling" )) lf->SetNoInitialAssembling();
    lf->SetCacheBlockSize(int(flags.GetNumFlag("cacheblocksize",1)));

    // cout << "type = " << typeid(*lf).name() << ", addr = " << lfp << endl;

    return lf;
  }



  template class S_LinearForm<double>;
  template class S_LinearForm<Complex>;

  template class T_LinearForm<double>;
  template class T_LinearForm<Complex>;

  // template class T_LinearForm< Vec<3> >;



  
  SparseVector PointEvaluationFunctional :: Assemble() const
  {
    typedef double SCAL;  // only real for the moment ...
    
    LocalHeap lh(10000);
    shared_ptr<FESpace> fespace;
    shared_ptr<MeshAccess> ma;
    
    cf->TraverseTree
      ([&] (CoefficientFunction& cf)
       {
         if (auto * proxy = dynamic_cast<ProxyFunction*> (&cf))
           {
             fespace = proxy->GetFESpace();
             ma = fespace->GetMeshAccess();
           }
       });

    if (!fespace)
      throw Exception("PointFunctional does not have Testfunction");
    
    SparseVector sv(fespace->GetNDof());
    

    
    IntegrationPoint ip;
    auto elnr = ma->FindElementOfPoint(point, ip, false);

    ElementId ei(VOL, elnr);
    auto & trafo = ma->GetTrafo(ei, lh);
    auto & fel = fespace->GetFE(ei, lh);
    Array<int> dnums(fel.GetNDof(), lh);
    fespace->GetDofNrs(ei, dnums);
    
    IntegrationRule ir (1, &ip);
    auto & mir = trafo(ir, lh);
    
    cf->TraverseTree
      ([&] (CoefficientFunction& cfi)
       {
         if (auto * proxy = dynamic_cast<ProxyFunction*> (&cfi))
           {
             FlatVector<SCAL> elvec(dnums.Size()*fespace->GetDimension(), lh);
             FlatMatrix<SCAL> values(1, 1, lh);
             ProxyUserData ud;
             const_cast<ElementTransformation&>(trafo).userdata = &ud;
             
             elvec = 0;
             
             FlatMatrix<SCAL> proxyvalues(1, proxy->Dimension(), lh);
             for (int k = 0; k < proxy->Dimension(); k++)
               {
                 ud.testfunction = proxy;
                 ud.test_comp = k;
                 
                 cf -> Evaluate (mir, values);
                 proxyvalues(0,k) = values(0,0);
               }
             proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec, lh);
             
             // AddElementVector (dnums, elvec);
             for (size_t i = 0; i < dnums.Size(); i++)
               sv[dnums[i]] = elvec[i];
           }
       });
    
    return sv;
  }
  
}
