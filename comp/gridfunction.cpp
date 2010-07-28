#include <comp.hpp>
#include <multigrid.hpp>

#ifdef PARALLEL
#include <parallelngs.hpp>
#endif

namespace ngcomp
{
  using namespace ngcomp;
  using namespace ngmg;

#ifdef PARALLEL
  using namespace ngparallel;
#endif

  GridFunction :: GridFunction (const FESpace & afespace, const string & name,
				const Flags & flags)
    : NGS_Object (afespace.GetMeshAccess(), name), fespace(afespace)
  { 
    nested = flags.GetDefineFlag ("nested");
    visual = !flags.GetDefineFlag ("novisual");
    multidim = int (flags.GetNumFlag ("multidim", 1));
    level_updated = -1;
    cacheblocksize = 1;
    // vis = 0;
  }

  GridFunction :: ~GridFunction()
  {
    for (int i = 0; i < vec.Size(); i++)
      delete vec[i];
  }



  void GridFunction :: Update ()
  {
    throw Exception("GridFunction::Update not overloaded");
  }

  bool GridFunction :: IsUpdated () const
  {
    int ndof = GetFESpace().GetNDof();
    for (int i = 0; i < multidim; i++)
      {
	if (!vec[i]) return false;
	if (ndof != vec[i]->Size()) return false;
      }
    return true;
  }



  void GridFunction :: PrintReport (ostream & ost)
  {
    ost << "on space " << GetFESpace().GetName() << endl
	<< "nested = " << nested << endl;
  }

  void GridFunction :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    //if (&const_cast<GridFunction&> (*this).GetVector())
    if (&(*this).GetVector())
      {
	int olds = mu.Size();
	//const_cast<GridFunction&> (*this).GetVector().MemoryUsage (mu);
	(*this).GetVector().MemoryUsage (mu);
	for (int i = olds; i < mu.Size(); i++)
	  mu[i]->AddName (string(" gf ")+GetName());
      }
  }







  void GridFunction :: Visualize(const string & given_name)
  {
    if (!visual) return;

    // cout << "visualize gridfunction " << name << endl;

    const BilinearFormIntegrator * bfi2d = 0, * bfi3d = 0;

    if (ma.GetDimension() == 2)
      {
	bfi2d = fespace.GetEvaluator();
      }
    else
      {
	bfi3d = fespace.GetEvaluator();
	bfi2d = fespace.GetBoundaryEvaluator();
      }

    if (bfi2d || bfi3d)
      {
	// if (ma.GetNLevels() > 1) return;
        // if (!vis)

	netgen::SolutionData * vis;
	if (!fespace.IsComplex())
	  vis = new VisualizeGridFunction<double> (ma, this, bfi2d, bfi3d, 0);
	else
	  vis = new VisualizeGridFunction<Complex> (ma, this, bfi2d, bfi3d, 0);

	Ng_SolutionData soldata;
	Ng_InitSolutionData (&soldata);
	
	soldata.name = given_name.c_str();
	soldata.data = 0;
	soldata.components = vis->GetComponents();
	soldata.iscomplex = vis->IsComplex();
	soldata.draw_surface = bfi2d != 0;
	soldata.draw_volume  = bfi3d != 0;
	soldata.dist = soldata.components;
	soldata.soltype = NG_SOLUTION_VIRTUAL_FUNCTION;
	soldata.solclass = vis;
	Ng_SetSolutionData (&soldata);    
      }
  }



  template <class SCAL>
  GridFunction * S_GridFunction<SCAL> :: GetComponent (int compound_comp) const
  {
    GridFunction * gf = new S_ComponentGridFunction<SCAL> (*this, compound_comp);
    if (level_updated != -1) gf->Update();
    return gf;
  }


  template <class SCAL>
  S_ComponentGridFunction<SCAL> :: S_ComponentGridFunction (const S_GridFunction<SCAL> & agf, int acomp)
    : S_GridFunction<SCAL> (*dynamic_cast<const CompoundFESpace&> (agf.GetFESpace())[acomp], agf.GetName(), Flags()),
      gf(agf), comp(acomp)
  { 
    ;
  }
  

  template <class SCAL>
  void S_ComponentGridFunction<SCAL> :: Update()
  {
    const CompoundFESpace & cfes = dynamic_cast<const CompoundFESpace&> (gf.GetFESpace());

    int begin = cfes.GetStorageStart(comp);
    int end = cfes.GetStorageEnd(comp);
  
    this -> vec.SetSize (gf.GetMultiDim());
    for (int i = 0; i < gf.GetMultiDim(); i++)
      (this->vec)[i] = gf.GetVector(i).Range (begin, end);
  
    this -> level_updated = this -> ma.GetNLevels();
  }







  template <class TV>
  T_GridFunction<TV> ::
  T_GridFunction (const FESpace & afespace, const string & aname, const Flags & flags)
    : S_GridFunction<TSCAL> (afespace, aname, flags)
  {
    vec.SetSize (this->multidim);
    vec = 0;

    Visualize (this->name);
  }

  /*
  template <class TV>
  T_GridFunction<TV> :: ~T_GridFunction()
  {
    for (int i = 0; i < vec.Size(); i++)
      delete vec[i];
  }
  */

  /*
  template <class TV>
  bool T_GridFunction<TV> :: IsUpdated (void) const
  {
    int ndof = this->GetFESpace().GetNDof();
    bool retval = true;

    for (int i = 0; retval && i < this->multidim; i++)
      retval = (vec[i] && ndof == vec[i]->Size());

    return retval;
  }
  */

  template <class TV>
  void T_GridFunction<TV> :: Update () 
  {
    try
      {
	int ndof = this->GetFESpace().GetNDof();

	for (int i = 0; i < this->multidim; i++)
	  {
	    
	    if (vec[i] && ndof == vec[i]->Size())
	      break;
	    
	    T_BaseVector<TV> * ovec = dynamic_cast<T_BaseVector<TV>*> (vec[i]);
	
#ifdef PARALLEL
	    *testout << &this->GetFESpace().GetParallelDofs() << endl;
	    if ( & this->GetFESpace().GetParallelDofs() )
	      vec[i] = new ParallelVVector<TV> (ndof, &this->GetFESpace().GetParallelDofs() );
	    else
#endif
 	      vec[i] = new VVector<TV> (ndof);


	    if (this->nested && ovec && this->GetFESpace().GetProlongation())
	      {
		(*vec[i]) = TSCAL(0);
		*vec[i]->Range (0, ovec->Size()) += (*ovec);

		const_cast<ngmg::Prolongation&> (*this->GetFESpace().GetProlongation()).Update();
		
		cout << "prolongate gridfunction" << endl;
		this->GetFESpace().GetProlongation()->ProlongateInline
		  (this->GetMeshAccess().GetNLevels()-1, *vec[i]);
	      }
	    else
	      {
		(*vec[i]) = TSCAL(0);
		// cout << "do not prolongate gridfunction" << endl;
	      }

	    //	    if (i == 0)
            // cout << "visualize" << endl;
            // Visualize (this->name);
	    
	    delete ovec;
	  }
	
	this -> level_updated = this -> ma.GetNLevels();


#ifdef PARALLEL
        const FESpace & afespace = GridFunction :: GetFESpace();
        for ( int i = 0; i < vec.Size(); i++ )
          {
            (vec[i]) -> SetParallelDofs ( &afespace.GetParallelDofs() );
            if ( &(afespace.GetParallelDofs() ) )
              (vec[i]) -> SetStatus ( CUMULATED );
            else
              (vec[i]) -> SetStatus ( NOT_PARALLEL );
          }
#endif

      }
    catch (exception & e)
      {
	Exception e2 (e.what());
	e2.Append ("\nIn GridFunction::Update()\n");
	throw e2;
      }
    catch (Exception & e)
      {
	e.Append ("In GridFunction::Update()\n");
	throw e;
      }    
  }





  GridFunction * CreateGridFunction (const FESpace * space,
				     const string & name, const Flags & flags)
  {
    GridFunction * gf = 
      CreateVecObject  <T_GridFunction, GridFunction, const FESpace, const string, const Flags>
      (space->GetDimension() * int(flags.GetNumFlag("cacheblocksize",1)), 
       space->IsComplex(), *space, name, flags);
  
    gf->SetCacheBlockSize(int(flags.GetNumFlag("cacheblocksize",1)));

    return gf;
  }



  GridFunctionCoefficientFunction :: 
  GridFunctionCoefficientFunction (GridFunction & agf, int acomp)
    : gf(dynamic_cast<S_GridFunction<double>&> (agf)),
      diffop (NULL),
      comp (acomp) 
  {
    ;
  }



  GridFunctionCoefficientFunction :: 
  GridFunctionCoefficientFunction (GridFunction & agf, DifferentialOperator * adiffop, int acomp)
    : gf(dynamic_cast<S_GridFunction<double>&> (agf)),
      diffop (adiffop),
      comp (acomp) 
  {
    ;
  }

  GridFunctionCoefficientFunction :: 
  ~GridFunctionCoefficientFunction ()
  {
    ;
  }

  int GridFunctionCoefficientFunction::Dimension() const{ 
    int res = -1;
    if (diffop==NULL){
      res = gf.GetFESpace().GetEvaluator()->DimFlux();
    }
    else{
      res = diffop->Dim();
    }
    return res;
  }

  

  double GridFunctionCoefficientFunction :: Evaluate (const BaseSpecificIntegrationPoint & ip) const
  {
    LocalHeapMem<100000> lh2 ("GridFunctionCoefficientFunction - evaluate");
    
    const int elnr = ip.GetTransformation().GetElementNr();
    bool boundary = ip.GetTransformation().Boundary();

    const FESpace & fes = gf.GetFESpace();
    const FiniteElement & fel = (boundary) ? fes.GetSFE(elnr, lh2) : fes.GetFE (elnr, lh2);
    const int dim     = fes.GetDimension();
    

    ArrayMem<int, 50> dnums;
    if(boundary)
      fes.GetSDofNrs(elnr, dnums);
    else
      fes.GetDofNrs (elnr, dnums);
    
    VectorMem<50> elu(dnums.Size()*dim);

    gf.GetElementVector (comp, dnums, elu);
    fes.TransformVec (elnr, boundary, elu, TRANSFORM_SOL);


    if (diffop)
      {
	VectorMem<10> flux(diffop->Dim());
	diffop->Apply (fel, ip, elu, flux, lh2);
	return flux(0);
      }
    else
      {
	const BilinearFormIntegrator * bfi = boundary ? fes.GetBoundaryEvaluator() : fes.GetEvaluator();
	VectorMem<10> flux(bfi->DimFlux());
	bfi->CalcFlux (fel, ip, elu, flux, false, lh2);
	return flux(0); 
      }
  }

  void GridFunctionCoefficientFunction :: Evaluate (const BaseSpecificIntegrationPoint & ip,
						    FlatVector<> result) const
  {
    LocalHeapMem<100000> lh2 ("GridFunctionCoefficientFunction, Eval 2");
    
    const int elnr = ip.GetTransformation().GetElementNr();
    bool boundary = ip.GetTransformation().Boundary();

    const FESpace & fes = gf.GetFESpace();
    const FiniteElement & fel = (boundary) ? fes.GetSFE(elnr, lh2) : fes.GetFE (elnr, lh2);
    const int dim     = fes.GetDimension();
    

    ArrayMem<int, 50> dnums;
    if(boundary)
      fes.GetSDofNrs(elnr, dnums);
    else
      fes.GetDofNrs (elnr, dnums);
    
    VectorMem<50> elu(dnums.Size()*dim);

    gf.GetElementVector (comp, dnums, elu);
    fes.TransformVec (elnr, boundary, elu, TRANSFORM_SOL);

    if (diffop)
      {
	diffop->Apply (fel, ip, elu, result, lh2);
      }
    else
      {
	const BilinearFormIntegrator * bfi = boundary ? fes.GetBoundaryEvaluator() : fes.GetEvaluator();
	bfi->CalcFlux (fel, ip, elu, result, false, lh2);
      }

  }


  void GridFunctionCoefficientFunction :: 
  Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    LocalHeapMem<100000> lh2("GridFunctionCoefficientFunction - Evalute 3");
    
    const int elnr = ir.GetTransformation().GetElementNr();
    bool boundary = ir.GetTransformation().Boundary();

    const FESpace & fes = gf.GetFESpace();
    const FiniteElement & fel = (boundary) ? fes.GetSFE(elnr, lh2) : fes.GetFE (elnr, lh2);
    const int dim     = fes.GetDimension();
    

    ArrayMem<int, 50> dnums;
    if(boundary)
      fes.GetSDofNrs(elnr, dnums);
    else
      fes.GetDofNrs (elnr, dnums);
    
    VectorMem<50> elu(dnums.Size()*dim);

    gf.GetElementVector (comp, dnums, elu);
    fes.TransformVec (elnr, boundary, elu, TRANSFORM_SOL);

    
    if (diffop)
      {
	diffop->Apply (fel, ir, elu, values, lh2);
      }
    else
      {
	const BilinearFormIntegrator * bfi = boundary ? fes.GetBoundaryEvaluator() : fes.GetEvaluator();
	bfi->CalcFlux (fel, ir, elu, values, false, lh2);
      }
  }





  template <class SCAL>
  VisualizeGridFunction<SCAL> ::
  VisualizeGridFunction (const MeshAccess & ama,
			 const GridFunction * agf,
			 const BilinearFormIntegrator * abfi2d,
			 const BilinearFormIntegrator * abfi3d,
			 bool aapplyd)

    : SolutionData (agf->GetName(), -1, agf->GetFESpace().IsComplex()),
      ma(ama), gf(dynamic_cast<const S_GridFunction<SCAL>*> (agf)), 
      applyd(aapplyd), cache_elnr(-1), lh(10000013, "VisualizedGridFunction 2"), fel(NULL)
  { 
    if(abfi2d)
      bfi2d.Append(abfi2d);
    if(abfi3d)
      bfi3d.Append(abfi3d);

    if (abfi2d) components = abfi2d->DimFlux();
    if (abfi3d) components = abfi3d->DimFlux();
    if (iscomplex) components *= 2;
    multidimcomponent = 0;
  }

  template <class SCAL>
  VisualizeGridFunction<SCAL> ::
  VisualizeGridFunction (const MeshAccess & ama,
			 const GridFunction * agf,
			 const Array<BilinearFormIntegrator *> & abfi2d,
			 const Array<BilinearFormIntegrator *> & abfi3d,
			 bool aapplyd)

    : SolutionData (agf->GetName(), -1, agf->GetFESpace().IsComplex()),
      ma(ama), gf(dynamic_cast<const S_GridFunction<SCAL>*> (agf)), 
      applyd(aapplyd), cache_elnr(-1), lh(10000002, "VisualizeGridFunction"), fel(NULL)
  { 
    for(int i=0; i<abfi2d.Size(); i++)
      bfi2d.Append(abfi2d[i]);
    for(int i=0; i<abfi3d.Size(); i++)
      bfi3d.Append(abfi3d[i]);
    

    if (bfi2d.Size()) components = bfi2d[0]->DimFlux();
    if (bfi3d.Size()) components = bfi3d[0]->DimFlux();

    if (iscomplex) components *= 2;
    multidimcomponent = 0;
  }

  template <class SCAL>
  VisualizeGridFunction<SCAL> :: ~VisualizeGridFunction ()
  {}
  

  template <class SCAL>
  bool VisualizeGridFunction<SCAL> :: GetValue (int elnr, 
						double lam1, double lam2, double lam3,
						double * values) 
  { 
    if (!bfi3d.Size()) return 0;
    if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

    const FESpace & fes = gf->GetFESpace();

    int dim     = fes.GetDimension();

    // NgLock lock(const_cast<S_GridFunction<SCAL>*> (gf) -> Mutex(), 1);

    if ( !fes.DefinedOn(ma.GetElIndex(elnr)) ) 
      return 0;

    if (cache_elnr != elnr || cache_bound)
      {
	lh.CleanUp();
	ma.GetElementTransformation (elnr, eltrans, lh);
	fes.GetDofNrs (elnr, dnums);
	elu.AssignMemory (dnums.Size() * dim, lh);

	if(gf->GetCacheBlockSize() == 1)
	  {
	    gf->GetElementVector (multidimcomponent, dnums, elu);
	  }
	else
	  {
	    FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
	    //gf->GetElementVector (multidimcomponent, dnums, elu2);
	    gf->GetElementVector (0, dnums, elu2);
	    for(int i=0; i<elu.Size(); i++)
	      elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
	  }

	fes.TransformVec (elnr, 0, elu, TRANSFORM_SOL);

	fel = &fes.GetFE (elnr, lh);

	cache_elnr = elnr;
	cache_bound = 0;

      }
    
    IntegrationPoint ip(lam1, lam2, lam3, 0);
    SpecificIntegrationPoint<3,3> sip (ip, eltrans, lh);


    for(int j = 0; j<bfi3d.Size(); j++)
      {
	HeapReset hr(lh);
	FlatVector<SCAL> flux(bfi3d[j] -> DimFlux(), lh);
	bfi3d[j]->CalcFlux (*fel, sip, elu, flux, applyd, lh);

	for (int i = 0; i < components; i++)
	  {
	    if(j == 0)
	      values[i] = ((double*)(void*)&flux(0))[i];
	    else
	      values[i] += ((double*)(void*)&flux(0))[i];
	  }
      }

    return 1; 
  }


  template <class SCAL>
  bool VisualizeGridFunction<SCAL> :: GetValue (int elnr, 
						const double xref[], const double x[], const double dxdxref[],
						double * values) 
  { 
    if (!bfi3d.Size()) return 0;
    if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

    //LocalHeapMem<10000> lh2;

    const FESpace & fes = gf->GetFESpace();
    // const FiniteElement & fel = fes.GetFE (elnr, lh2);

    int dim     = fes.GetDimension();

    // NgLock lock(const_cast<S_GridFunction<SCAL>*> (gf) -> Mutex(), 1);

    // added 07.04.2004 (FB): don't draw, if not defined on this domain
    if ( !fes.DefinedOn(ma.GetElIndex(elnr)) ) return 0;

    if (cache_elnr != elnr || cache_bound)
      {
	lh.CleanUp();
	ma.GetElementTransformation (elnr, eltrans, lh);
	fes.GetDofNrs (elnr, dnums);
    
	elu.AssignMemory (dnums.Size() * dim, lh);
	if(gf->GetCacheBlockSize() == 1)
	  {
	    gf->GetElementVector (multidimcomponent, dnums, elu);
	  }
	else
	  {
	    FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
	    //gf->GetElementVector (multidimcomponent, dnums, elu2);
	    gf->GetElementVector (0, dnums, elu2);
	    for(int i=0; i<elu.Size(); i++)
	      elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
	  }
	fes.TransformVec (elnr, 0, elu, TRANSFORM_SOL);

	fel = &fes.GetFE (elnr, lh);

	cache_elnr = elnr;
	cache_bound = 0;
      }

    void * hp = lh.GetPointer();

    

    Vec<3> vx;
    Mat<3,3> mdxdxref;
    for (int i = 0; i < 3; i++)
      {
	vx(i) = x[i];
	for (int j = 0; j < 3; j++)
	  mdxdxref(i,j) = dxdxref[3*i+j];
      }

    IntegrationPoint ip(xref[0], xref[1], xref[2], 0);
    SpecificIntegrationPoint<3,3> sip (ip, eltrans, vx, mdxdxref); // , lh);

    for(int j = 0; j<bfi3d.Size(); j++)
      {
	FlatVector<SCAL> flux (bfi3d[j]->DimFlux(), lh);
	bfi3d[j]->CalcFlux (*fel, sip, elu, flux, applyd, lh);


	for (int i = 0; i < components; i++)
	  {
	    if(j == 0)
	      values[i] = ((double*)(void*)&flux(0))[i];
	    else
	      values[i] += ((double*)(void*)&flux(0))[i];
	  }
      }

    lh.CleanUp(hp);
    return 1; 
  }




  template <class SCAL>
  bool VisualizeGridFunction<SCAL> ::
  GetMultiValue (int elnr, int npts,
		 const double * xref, int sxref,
		 const double * x, int sx,
		 const double * dxdxref, int sdxdxref,
		 double * values, int svalues)
  {
    try
      {
        if (!bfi3d.Size()) return 0;
        if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

        const FESpace & fes = gf->GetFESpace();
        int dim = fes.GetDimension();
	
        HeapReset hr(lh);
	
	ma.GetElementTransformation (elnr, eltrans, lh);
	fes.GetDofNrs (elnr, dnums);
	fel = &fes.GetFE (elnr, lh);
        
        elu.AssignMemory (dnums.Size() * dim, lh);

        if(gf->GetCacheBlockSize() == 1)
          {
            gf->GetElementVector (multidimcomponent, dnums, elu);
          }
        else
          {
            FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
            gf->GetElementVector (0, dnums, elu2);
            for(int i=0; i<elu.Size(); i++)
              elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
          }
        
        fes.TransformVec (elnr, false, elu, TRANSFORM_SOL);
        
	if (!fes.DefinedOn(eltrans.GetElementIndex()))return 0;

	IntegrationRule ir(npts);
	for (int i = 0; i < npts; i++)
	  ir.Append (IntegrationPoint (xref[i*sxref], xref[i*sxref+1], xref[i*sxref+2]));

	// ElementTransformation eltrans;
	// ma.GetElementTransformation (elnr, eltrans, lh);

	MappedIntegrationRule<3,3> mir(ir, eltrans, 1, lh);

	for (int k = 0; k < npts; k++)
	  {
	    Mat<3,3> & mdxdxref = *new((double*)(dxdxref+k*sdxdxref)) Mat<3,3>;
	    FlatVec<3> vx( (double*)x + k*sx);
	    mir[k] = SpecificIntegrationPoint<3,3> (ir[k], eltrans, vx, mdxdxref);
	  }

	for (int k = 0; k < npts; k++)
	  for (int i = 0; i < components; i++)
	    values[k*svalues+i] = 0.0;

	for(int j = 0; j < bfi3d.Size(); j++)
	  {
	    HeapReset hr(lh);

	    FlatMatrix<SCAL> flux(npts, bfi3d[j]->DimFlux(), lh);
	    bfi3d[j]->CalcFlux (*fel, mir, elu, flux, applyd, lh);
	    
	    for (int k = 0; k < npts; k++)
	      for (int i = 0; i < components; i++)
		values[k*svalues+i] += ((double*)(void*)&flux(k,0))[i];
          }

        return 1; 
      }
    catch (Exception & e)
      {
        cout << "GetMultiValue caught exception" << endl
             << e.What();
        return 0;
      }
  }
































  template <class SCAL>
  bool VisualizeGridFunction<SCAL> :: GetSurfValue (int elnr,
                                                    double lam1, double lam2, 
                                                    double * values) 
  { 
    if (!bfi2d.Size()) return 0;
    if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

    bool bound = (ma.GetDimension() == 3);
    const FESpace & fes = gf->GetFESpace();


    // NgLock lock(const_cast<S_GridFunction<SCAL>*> (gf) -> Mutex(), 1);

    int dim = fes.GetDimension();

    /*
      if ( bound ? 
      !fes.DefinedOnBoundary(ma.GetSElIndex(elnr)) :
      !fes.DefinedOn(ma.GetElIndex(elnr)) ) return 0;
    */

    if (cache_elnr != elnr || !cache_bound)
      {
	lh.CleanUp();

	if (bound)
	  {
	    ma.GetSurfaceElementTransformation (elnr, eltrans, lh);
	    fes.GetSDofNrs (elnr, dnums);
	    fel = &fes.GetSFE (elnr, lh);
	  }
	else
	  {
	    ma.GetElementTransformation (elnr, eltrans, lh);
	    fes.GetDofNrs (elnr, dnums);
	    fel = &fes.GetFE (elnr, lh);
	  }

	elu.AssignMemory (dnums.Size() * dim, lh);

	if(gf->GetCacheBlockSize() == 1)
	  {
	    gf->GetElementVector (multidimcomponent, dnums, elu);
	  }
	else
	  {
	    FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
	    //gf->GetElementVector (multidimcomponent, dnums, elu2);
	    gf->GetElementVector (0, dnums, elu2);
	    for(int i=0; i<elu.Size(); i++)
	      elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
	  }

	fes.TransformVec (elnr, bound, elu, TRANSFORM_SOL);

	cache_elnr = elnr;
	cache_bound = 1;
      }

    if ( bound ? 
	 !fes.DefinedOnBoundary(eltrans.GetElementIndex()) : 
	 !fes.DefinedOn(eltrans.GetElementIndex()) ) return 0;


    HeapReset hr(lh);


    IntegrationPoint ip(lam1, lam2, 0, 0);


    if (bound)
      {
	SpecificIntegrationPoint<2,3> sip (ip, eltrans, lh);
	for(int j = 0; j<bfi2d.Size(); j++)
	  {
            FlatVector<SCAL> flux(bfi2d[j]->DimFlux(), lh);
	    bfi2d[j]->CalcFlux (*fel, sip, elu, flux, applyd, lh);

	    for (int i = 0; i < components; i++)
	      {
		if(j == 0)
		  values[i] = ((double*)(void*)&flux(0))[i];
		else
		  values[i] += ((double*)(void*)&flux(0))[i];
	      }
	  }
      }
    else
      {
	SpecificIntegrationPoint<2,2> sip (ip, eltrans, lh);
	for(int j = 0; j<bfi2d.Size(); j++)
	  {
            FlatVector<SCAL> flux(bfi2d[j]->DimFlux(), lh);
	    bfi2d[j]->CalcFlux (*fel, sip, elu, flux, applyd, lh);

	    for (int i = 0; i < components; i++)
	      {
		if(j == 0)
		  values[i] = ((double*)(void*)&flux(0))[i];
		else
		  values[i] += ((double*)(void*)&flux(0))[i];
	      }
	  }
      }

    return 1; 
  }




  template <class SCAL>
  bool VisualizeGridFunction<SCAL> :: GetSurfValue (int elnr,
						    const double xref[], const double x[], const double dxdxref[],
						    double * values) 
  { 
    try
      {
        if (!bfi2d.Size()) return 0;
        if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

        bool bound = (ma.GetDimension() == 3);
        const FESpace & fes = gf->GetFESpace();


        // NgLock lock(const_cast<S_GridFunction<SCAL>*> (gf) -> Mutex(), 1);

        int dim     = fes.GetDimension();


        if (cache_elnr != elnr || !cache_bound)
          {
            lh.CleanUp();

            if (bound)
              {
                ma.GetSurfaceElementTransformation (elnr, eltrans, lh);
                fes.GetSDofNrs (elnr, dnums);
                fel = &fes.GetSFE (elnr, lh);
              }
            else
              {
                ma.GetElementTransformation (elnr, eltrans, lh);
                fes.GetDofNrs (elnr, dnums);
                fel = &fes.GetFE (elnr, lh);
              }

            elu.AssignMemory (dnums.Size() * dim, lh);

            if(gf->GetCacheBlockSize() == 1)
              {
                gf->GetElementVector (multidimcomponent, dnums, elu);
              }
            else
              {
                FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
                //gf->GetElementVector (multidimcomponent, dnums, elu2);
                gf->GetElementVector (0, dnums, elu2);
                for(int i=0; i<elu.Size(); i++)
                  elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
              }

            fes.TransformVec (elnr, bound, elu, TRANSFORM_SOL);

            cache_elnr = elnr;
            cache_bound = 1;
          }


        /*
          if ( bound ? 
          !fes.DefinedOnBoundary(ma.GetSElIndex(elnr)) :
          !fes.DefinedOn(ma.GetElIndex(elnr)) ) return 0;
        */
        if ( bound ? 
             !fes.DefinedOnBoundary(eltrans.GetElementIndex()) : 
             !fes.DefinedOn(eltrans.GetElementIndex()) ) return 0;



        HeapReset hr(lh);


        IntegrationPoint ip(xref[0], xref[1], 0, 0);
        if (bound)
          {
            Vec<3> vx;
            Mat<3,2> mdxdxref;
            for (int i = 0; i < 3; i++)
              {
                vx(i) = x[i];
                for (int j = 0; j < 2; j++)
                  mdxdxref(i,j) = dxdxref[2*i+j];
              }
            SpecificIntegrationPoint<2,3> sip (ip, eltrans, vx, mdxdxref); 
            for (int i = 0; i < components; i++)
              values[i] = 0.0;
            for(int j = 0; j<bfi2d.Size(); j++)
              {
		FlatVector<SCAL> flux(bfi3d[j]->DimFlux(), lh);
                bfi2d[j]->CalcFlux (*fel, sip, elu, flux, applyd, lh);
                for (int i = 0; i < components; i++)
                  values[i] += ((double*)(void*)&flux(0))[i];
              }
          }
        else
          {
            Vec<2> vx;
            Mat<2,2> mdxdxref;
            for (int i = 0; i < 2; i++)
              {
                vx(i) = x[i];
                for (int j = 0; j < 2; j++)
                  mdxdxref(i,j) = dxdxref[2*i+j];
              }
            SpecificIntegrationPoint<2,2> sip (ip, eltrans, vx, mdxdxref); 

            for (int i = 0; i < components; i++)
              values[i] = 0.0;
            for(int j = 0; j<bfi2d.Size(); j++)
              {
                FlatVector<SCAL> flux(bfi2d[j]->DimFlux(), lh);
                bfi2d[j]->CalcFlux (*fel, sip, elu, flux, applyd, lh);
                for (int i = 0; i < components; i++)
                  values[i] += ((double*)(void*)&flux(0))[i];
              }
          }

        return 1; 
      }
    catch (Exception & e)
      {
        cout << "GetSurfaceValue2 caught exception" << endl
             << e.What();

        return 0;
      }
      
  }






  template <class SCAL>
  bool VisualizeGridFunction<SCAL> ::
  GetMultiSurfValue (int elnr, int npts,
                     const double * xref, int sxref,
                     const double * x, int sx,
                     const double * dxdxref, int sdxdxref,
                     double * values, int svalues)
  {
    try
      {
        if (!bfi2d.Size()) return 0;
        if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

        bool bound = (ma.GetDimension() == 3);
	
        const FESpace & fes = gf->GetFESpace();
        int dim = fes.GetDimension();

        
        HeapReset hr(lh);
        
        if (bound)
          {
            ma.GetSurfaceElementTransformation (elnr, eltrans, lh);
            fes.GetSDofNrs (elnr, dnums);
            fel = &fes.GetSFE (elnr, lh);
          }
        else
          {
            ma.GetElementTransformation (elnr, eltrans, lh);
            fes.GetDofNrs (elnr, dnums);
            fel = &fes.GetFE (elnr, lh);
          }
        
        elu.AssignMemory (dnums.Size() * dim, lh);

        if(gf->GetCacheBlockSize() == 1)
          {
            gf->GetElementVector (multidimcomponent, dnums, elu);
          }
        else
          {
            FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
            gf->GetElementVector (0, dnums, elu2);
            for(int i=0; i<elu.Size(); i++)
              elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
          }
        
        fes.TransformVec (elnr, bound, elu, TRANSFORM_SOL);

        
        if ( bound ? 
             !fes.DefinedOnBoundary(eltrans.GetElementIndex()) : 
             !fes.DefinedOn(eltrans.GetElementIndex()) ) return 0;




        if (bound)
          {

            for (int k = 0; k < npts; k++)
              for (int i = 0; i < components; i++)
                values[k*svalues+i] = 0.0;
            
            for (int k = 0; k < npts; k++)
              {
                HeapReset hr(lh);
                
                IntegrationPoint ip(xref[k*sxref], xref[k*sxref+1], 0, 0);
                Vec<3> vx;
                Mat<3,2> mdxdxref;

                for (int i = 0; i < 3; i++)
                  vx(i) = x[k*sx+i];
                for (int i = 0; i < 3; i++)
                  for (int j = 0; j < 2; j++)
                    mdxdxref(i,j) = dxdxref[k*sdxdxref+2*i+j];

		SpecificIntegrationPoint<2,3> sip (ip, eltrans, vx, mdxdxref); 
		// SpecificIntegrationPoint<2,3> sip (ip, eltrans, x+k*sx, mdxdxref); 
                
                for(int j = 0; j<bfi2d.Size(); j++)
                  {
		    FlatVector<SCAL> flux (bfi2d[j]->DimFlux(), lh);
                    bfi2d[j]->CalcFlux (*fel, sip, elu, flux, applyd, lh);
                    for (int i = 0; i < components; i++)
                      values[k*svalues+i] += ((double*)(void*)&flux(0))[i];
                  }
              }
          }
        else
          {
	    IntegrationRule ir(npts);
	    for (int i = 0; i < npts; i++)
	      ir.Append (IntegrationPoint (xref[i*sxref], xref[i*sxref+1]));

	    ElementTransformation eltrans;
	    ma.GetElementTransformation (elnr, eltrans, lh);

	    MappedIntegrationRule<2,2> mir(ir, eltrans, 1, lh);

	    for (int k = 0; k < npts; k++)
	      {
		Mat<2,2> & mdxdxref = *new((double*)(dxdxref+k*sdxdxref)) Mat<2,2>;
		FlatVec<2> vx( (double*)x + k*sx);
		mir[k] = SpecificIntegrationPoint<2,2> (ir[k], eltrans, vx, mdxdxref);
	      }

            for (int k = 0; k < npts; k++)
              for (int i = 0; i < components; i++)
                values[k*svalues+i] = 0.0;
	    
	    for(int j = 0; j<bfi2d.Size(); j++)
	      {
		FlatMatrix<SCAL> flux(npts, bfi2d[j]->DimFlux(), lh);
		bfi2d[j]->CalcFlux (*fel, mir, elu, flux, applyd, lh);

		for (int k = 0; k < npts; k++)
		  for (int i = 0; i < components; i++)
		    values[k*svalues+i] += ((double*)(void*)&flux(k,0))[i];
	      }
          }

        return 1; 
      }
    catch (Exception & e)
      {
        cout << "GetMultiSurfaceValue caught exception" << endl
             << e.What();

        return 0;
      }
  }














  
  template <class SCAL>
  void VisualizeGridFunction<SCAL> :: 
  Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages, int component)
  {
    int ndomains = 0;

    if (bfi3d.Size()) 
      ndomains = ma.GetNDomains();
    else if(bfi2d.Size()) 
      ndomains = ma.GetNBoundaries();

    Array<double> volumes(ndomains);

    Analyze(minima,maxima,averages,volumes,component);
    
    for(int i=0; i<ndomains; i++)
      {
	if(component == -1)
	  for(int j=0; j<components; j++)
	    {
	      averages[i*components+j] /= volumes[i];
	    }
	else
	  averages[i] /= volumes[i];
      }
  }
  

  template <class SCAL>
  void VisualizeGridFunction<SCAL> :: 
  Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages_times_volumes, Array<double> & volumes, int component)
  {
    const FESpace & fes = gf->GetFESpace();

    int domain;
    double *val;
    int pos;
    double vol;

    int ndomains;


    if(bfi3d.Size()) ndomains = ma.GetNDomains();
    else if(bfi2d.Size()) ndomains = ma.GetNBoundaries();

    Array<double> posx;
    Array<double> posy;
    Array<double> posz;
    ELEMENT_TYPE cache_type = ET_SEGM;
	
    LocalHeapMem<10000> lh2("Gridfunction - Analyze");
	
    val = new double[components];
			

    for(int i=0; i<minima.Size(); i++)
      {
	minima[i] = 1e100;
	maxima[i] = -1e100;
	averages_times_volumes[i] = 0;
      }

    for(int i=0; i<volumes.Size(); i++) volumes[i] = 0;
    
    void * heapp = lh2.GetPointer();
    if(bfi3d.Size())
      {
	for(int i=0; i<ma.GetNE(); i++)
	  {
	    const FiniteElement & fel = fes.GetFE(i,lh2);
	    
	    domain = ma.GetElIndex(i);
	    
	    vol = ma.ElementVolume(i);
	    
	    volumes[domain] += vol;
	    
	    if(fel.ElementType() != cache_type)
	      {
		posx.DeleteAll(); posy.DeleteAll(); posz.DeleteAll(); 
		switch(fel.ElementType())
		  {
		  case ET_TET:
		    posx.Append(0.25); posy.Append(0.25); posz.Append(0.25); 
		    posx.Append(0); posy.Append(0); posz.Append(0); 
		    posx.Append(1); posy.Append(0); posz.Append(0);
		    posx.Append(0); posy.Append(1); posz.Append(0);
		    posx.Append(0); posy.Append(0); posz.Append(1);
		    break;
		  case ET_HEX:
		    posx.Append(0.5); posy.Append(0.5); posz.Append(0.5);  
		    posx.Append(0); posy.Append(0); posz.Append(0); 
		    posx.Append(0); posy.Append(0); posz.Append(1);
		    posx.Append(0); posy.Append(1); posz.Append(0);
		    posx.Append(0); posy.Append(1); posz.Append(1);  
		    posx.Append(1); posy.Append(0); posz.Append(0); 
		    posx.Append(1); posy.Append(0); posz.Append(1);
		    posx.Append(1); posy.Append(1); posz.Append(0);
		    posx.Append(1); posy.Append(1); posz.Append(1);
		    break;
                  default:
                    {
                      bool firsttime = 1;
                      if (firsttime)
                        {
                          cerr << "WARNING::VisGridFunction::Analyze: unsupported element "
                               << ElementTopology::GetElementName(fel.ElementType()) << endl;
                          firsttime = 0;
                        }
                      break;
                    } 
		  }
		cache_type = fel.ElementType();
	      }
	    
	    
	    for(int k=0; k<posx.Size(); k++)
	      {
		GetValue(i,posx[k],posy[k],posz[k],val);


		if(component == -1)
		  {
		    for(int j=0; j<components; j++)
		      {
			pos = domain*components+j;
			if(val[j] > maxima[pos]) maxima[pos] = val[j];
			if(val[j] < minima[pos]) minima[pos] = val[j];
			if(k==0) averages_times_volumes[pos] += val[j]*vol;
		      }
		  }
		else
		  {
		    pos = domain;
		    if(val[component] > maxima[pos]) maxima[pos] = val[component];
		    if(val[component] < minima[pos]) minima[pos] = val[component];
		    if(k==0) averages_times_volumes[pos] += val[component]*vol;
		  }
	      }
	    lh2.CleanUp(heapp);
	  }
      }
    else if (bfi2d.Size())
      {
	for(int i=0; i<ma.GetNSE(); i++)
	  {
	    const FiniteElement & fel = fes.GetSFE(i,lh2);

	    domain = ma.GetSElIndex(i);

	    vol = ma.SurfaceElementVolume(i);

	    volumes[domain] += vol;

	    if(fel.ElementType() != cache_type)
	      {
		posx.DeleteAll(); posy.DeleteAll(); posz.DeleteAll(); 
		switch(fel.ElementType())
		  {
		  case ET_SEGM: 
		  case ET_TET: case ET_HEX: case ET_PRISM: case ET_PYRAMID:
		    break;
		    
		  case ET_TRIG:
		    posx.Append(0.33); posy.Append(0.33);
		    posx.Append(0); posy.Append(0);
		    posx.Append(0); posy.Append(1);
		    posx.Append(1); posy.Append(0);
		    break;
		  case ET_QUAD:
		    posx.Append(0.5); posy.Append(0.5);
		    posx.Append(0); posy.Append(0);
		    posx.Append(0); posy.Append(1); 
		    posx.Append(1); posy.Append(0);
		    posx.Append(1); posy.Append(1);
		    break;
		  }
		cache_type = fel.ElementType();
	      }
	    for(int k=0; k<posx.Size(); k++)
	      {
		GetSurfValue(i,posx[k],posy[k],val);
		if(component == -1)
		  {
		    for(int j=0; j<components; j++)
		      {
			pos = domain*components+j;
			if(val[j] > maxima[pos]) maxima[pos] = val[j];
			if(val[j] < minima[pos]) minima[pos] = val[j];
			if(k==0) averages_times_volumes[pos] += val[j]*vol;
		      }
		  }
		else
		  {
		    pos = domain;
		    if(val[component] > maxima[pos]) maxima[pos] = val[component];
		    if(val[component] < minima[pos]) minima[pos] = val[component];
		    if(k==0) averages_times_volumes[pos] += val[component]*vol;
		  }
	      }
	    lh2.CleanUp(heapp);
	  }
      }

    delete [] val;
  }



  template class T_GridFunction<double>;
  template class T_GridFunction<Vec<2> >;
  template class T_GridFunction<Vec<3> >;
  template class T_GridFunction<Vec<4> >;
  

  template class  VisualizeGridFunction<double>;
  template class  VisualizeGridFunction<Complex>;

}
