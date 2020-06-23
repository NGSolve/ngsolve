/*********************************************************************/
/* File:   integrator.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/
/* 
   Finite Element Integrators
*/

#include <fem.hpp>
  
namespace ngfem
{

  DifferentialOperator :: ~DifferentialOperator ()
  {
    ;
  }

  Integrator :: Integrator() throw ()
  {
    SetHigherIntegrationOrder(20);
    SetIntegrationOrder (-1);
    SetCacheComp (0);
    // fast = false;
    const_coef = false;
    // checkfast = false;
    name = "Integrator";
  }

  ///
  Integrator :: ~Integrator() 
  {
    DeleteCurveIPs();
  }

  ostream & operator << (ostream & ost, const Integrator & igt)
  {
    ost << igt.Name();
    return ost;
  }
  
  
  bool Integrator :: DefinedOn (int mat) const
  {
    if(mat<0)
      return false;
    if (definedon.Size())
      {
        if (mat >= definedon.Size()) return false;
        return definedon.Test(mat);
      }
    return true;
  }
  

  void Integrator :: SetDefinedOn (const BitArray & adefinedon)
  {
    definedon = adefinedon;
  }

  void Integrator :: SetDefinedOn (const Array<int> & regions)
  {
    int maxval = 0;
    for (int val : regions) maxval = max2(maxval, val);

    definedon.SetSize (maxval+1);
    definedon.Clear();

    for (int val : regions)
      if (val >= 0)
        definedon.SetBit (val);
      else
        throw Exception ("SetDefinedOn called with negative index");
  }

  void Integrator :: SetName (const string & aname)
  { 
    name = aname; 
  }

  string Integrator :: Name () const
  {
    return name;   // string("Integrator ") + typeid(*this).name(); 
  }

    

  void Integrator :: SetIntegrationAlongCurve ( const int npoints )
  {
    is_curve_integrator = true;
    curve_ips.SetSize(npoints); 
    curve_ip_tangents.SetSize(npoints);
    for(int i=0; i<npoints; i++)
      {
	curve_ips[i] = new Vector<>(3);
	curve_ip_tangents[i] = new Vector<>(3);
	*(curve_ip_tangents[i]) = 0;
      }
  }

  void Integrator :: DeleteCurveIPs ( void )
  {
    for(int i=0; i<curve_ips.Size(); i++)
      delete curve_ips[i];
    curve_ips.DeleteAll();

    for(int i=0; i<curve_ip_tangents.Size(); i++)
      delete curve_ip_tangents[i];
    curve_ip_tangents.DeleteAll();

    continuous_curveparts.DeleteAll();
  }

  void Integrator :: UnSetIntegrationAlongCurve ( void )
  {
    DeleteCurveIPs();
  }

  

  void Integrator :: AppendCurvePoint(const FlatVector<double> & point, const FlatVector<double> & tangent)
  {
    is_curve_integrator = true;
    if(continuous_curveparts.Size() == 0)
      continuous_curveparts.Append(0);
    Vector<> * vec = new Vector<>(3);
    *vec = point;
    curve_ips.Append(vec);

    vec = new Vector<>(3);
    *vec = tangent;
    curve_ip_tangents.Append(vec);

  }
  void Integrator :: AppendCurvePoint(const FlatVector<double> & point)
  {
    Vector<> tangent(3);
    tangent = 0;
    AppendCurvePoint(point,tangent);    
  }
  
  void Integrator :: SetCurveClearance(void)
  {
    continuous_curveparts.Append(curve_ips.Size());
  }

  int Integrator :: GetNumCurveParts(void) const
  {
    return continuous_curveparts.Size();
  }
      
  int Integrator :: GetStartOfCurve(const int i) const
  {
    return continuous_curveparts[i];
  }
  int Integrator :: GetEndOfCurve(const int i) const
  {
    if(i+1 < continuous_curveparts.Size())
      return continuous_curveparts[i+1];
    else
      return curve_ips.Size();
  }
  
  void Integrator :: SetFileName(const string & filename)
  {
    cerr << "SetFileName not defined for Integrator base class" << endl;
  }
  

  BilinearFormIntegrator :: ~BilinearFormIntegrator ()
  {
    ;
  }


  void BilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<Complex> elmat,
		     LocalHeap & lh) const
  {
    FlatMatrix<double> rmat (elmat.Height(), elmat.Width(), lh);
    CalcElementMatrix (fel, eltrans, rmat, lh);
    elmat = rmat;
  }

  
  void BilinearFormIntegrator ::
  CalcElementMatrixAdd (const FiniteElement & fel,
                        const ElementTransformation & eltrans, 
                        FlatMatrix<double> elmat,
                        bool & symmetric_so_far,
                        LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<double> helmat(elmat.Height(), elmat.Width(), lh);
    CalcElementMatrix(fel, eltrans, helmat, lh);
    elmat += helmat;
    if (!IsSymmetric().IsTrue()) symmetric_so_far = false;
  }

  void BilinearFormIntegrator ::    
  CalcElementMatrixAdd (const FiniteElement & fel,
                        const ElementTransformation & eltrans, 
                        FlatMatrix<Complex> elmat,
                        bool & symmetric_so_far,                        
                        LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<Complex> helmat(elmat.Height(), elmat.Width(), lh);
    CalcElementMatrix(fel, eltrans, helmat, lh);
    elmat += helmat;
    if (!IsSymmetric().IsTrue()) symmetric_so_far = false;    
  }
  


  
  void BilinearFormIntegrator ::
  CalcElementMatrixDiag (const FiniteElement & fel,
			     const ElementTransformation & eltrans, 
			     FlatVector<double> diag,
			     LocalHeap & lh) const
  {
    static bool first_time = true;
    if (first_time)
      {
        cout << "WARNING: base class, assemble diag is highly inefficient" << endl;
        first_time = false;
      }

    FlatMatrix<> elmat(diag.Size(), lh);
    CalcElementMatrix (fel, eltrans, elmat, lh);

    // diag.AssignMemory (elmat.Height(), lh);
    for (int i = 0; i < diag.Size(); i++)
      diag(i) = elmat(i,i);
  }


  void BilinearFormIntegrator ::
  CalcLinearizedElementMatrix (const FiniteElement & fel, 
			       const ElementTransformation & eltrans, 
			       FlatVector<double> elveclin,
			       FlatMatrix<double> elmat,
			       LocalHeap & lh) const
  {
    CalcElementMatrix (fel, eltrans, elmat, lh);
  }

  void BilinearFormIntegrator ::
  CalcLinearizedElementMatrix (const FiniteElement & fel, 
			       const ElementTransformation & eltrans, 
			       FlatVector<Complex> elveclin,
			       FlatMatrix<Complex> elmat,
			       LocalHeap & lh) const
  {
    CalcElementMatrix (fel, eltrans, elmat, lh);
  }



  void BilinearFormIntegrator ::
  ApplyElementMatrix (const FiniteElement & fel, 
		      const ElementTransformation & eltrans, 
                      FlatVector<double> elx, 
		      FlatVector<double> ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    static atomic<int> cnt(0);
    static mutex m;
    if (cnt < 3)
      {
        lock_guard<mutex> guard(m);
        if (cnt < 3) 
          cout << "WARNING: call baseclass ApplyElementMatrix, type = " << typeid(*this).name() << endl;
	if (cnt == 2) cout << "(further warnings suppressed)" << endl;
	cnt++;
      }
    FlatMatrix<double> mat(elx.Size(), lh);
    CalcElementMatrix (fel, eltrans, mat, lh);
    ely = mat * elx;
  }

  void BilinearFormIntegrator ::
  ApplyElementMatrix (const FiniteElement & fel, 
		      const ElementTransformation & eltrans, 
                      FlatVector<Complex> elx, 
		      FlatVector<Complex> ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    static atomic<int> cnt(0);
    static mutex m;
    if (cnt < 3)
      {
        lock_guard<mutex> guard(m);
        if (cnt < 3)
          cout << "WARNING: call baseclass ApplyElementMatrix complex, type = " << typeid(*this).name() << endl;
        if (cnt == 2) cout << "(further warnings suppressed)" << endl;
	cnt++;
      }

    FlatMatrix<Complex> mat(elx.Size(), lh);
    CalcElementMatrix (fel, eltrans, mat, lh);
    ely = mat * elx;
  }


  void BilinearFormIntegrator ::
  ApplyLinearizedElementMatrix (const FiniteElement & fel, 
				const ElementTransformation & eltrans, 
                                FlatVector<double> ellin,
                                FlatVector<double> elx, 
				FlatVector<double> ely,
				LocalHeap & lh) const
  {
    ApplyElementMatrix (fel, eltrans, elx, ely, 0, lh);
  }

  void BilinearFormIntegrator :: 
  ApplyLinearizedElementMatrix (const FiniteElement & fel, 
				const ElementTransformation & eltrans, 
                                FlatVector<Complex> ellin, 
                                FlatVector<Complex> elx, 
				FlatVector<Complex> ely,
				LocalHeap & lh) const
  {
    ApplyElementMatrix (fel, eltrans, elx, ely, 0, lh);
  }


  double BilinearFormIntegrator ::
  Energy (const FiniteElement & fel, 
	  const ElementTransformation & eltrans, 
          FlatVector<double> elx, 
	  LocalHeap & lh) const
  {
    FlatVector<double> ely (elx.Size(), lh);
    ApplyElementMatrix (fel, eltrans, elx, ely, 0, lh);
    return 0.5 * InnerProduct (elx, ely);
  }

  
  double BilinearFormIntegrator :: 
  Energy (const FiniteElement & fel, 
	  const ElementTransformation & eltrans, 
          FlatVector<Complex> elx, 
	  LocalHeap & lh) const
  {
    cout << "error: Energy for Complex vector called" << endl;
    return 0;
  }


  void BilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationPoint & bmip,
            BareSliceVector<double> elx, 
	    FlatVector<double> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    /*
    cerr << "calcflux<double> called for base class, should be overloaded in " 
	 << typeid(*this).name()
	 << endl;
    */
    throw Exception(string("calcflux<double> called for base class, should be overloaded in ")+
                    typeid(*this).name());
  }


  void BilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationPoint & bmip,
            BareSliceVector<Complex> elx, 
	    FlatVector<Complex> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    cerr << "calcflux<Complex> called for base class, should be overloaded in " 
	 << typeid(*this).name()
	 << endl;
  }



  void BilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & fel,
	    const FiniteElement & felflux,
	    const ElementTransformation & eltrans,
            BareSliceVector<> elx, 
	    FlatVector<> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    cerr << "calcflux<double> fel-fel called for base class " 
	 << typeid(*this).name()
	 << endl;
  }


  void BilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationRule & mir,
            BareSliceVector<double> elx, 
	    BareSliceMatrix<double> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    for (int l = 0; l < mir.Size(); l++)
      {
	FlatVector<> res = flux.Row(l).Range(0,DimFlux());
	CalcFlux (fel, mir[l], elx, res, applyd, lh);
      }
  }
  
  
  void BilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationRule & mir,
            BareSliceVector<Complex> elx, 
	    BareSliceMatrix<Complex> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    for (int l = 0; l < mir.Size(); l++)
      {
	FlatVector<Complex> res = flux.Row(l).Range(0,DimFlux());
	CalcFlux (fel, mir[l], elx, res, applyd, lh);
      }
  }





  void BilinearFormIntegrator :: 
  CalcFluxMulti (const FiniteElement & fel,
		 const BaseMappedIntegrationPoint & bmip,
		 int m,
                 FlatVector<double> elx, 
		 FlatVector<double> flux,
		 bool applyd,
		 LocalHeap & lh) const
  {
    FlatVector<double> selx(elx.Size()/m, lh);
    FlatVector<double> sflux(DimFlux(), lh);
    // flux.AssignMemory (m*DimFlux(), lh);
    for (int j = 0; j < m; j++)
      {
	for (int i = 0; i < selx.Size(); i++)
	  selx(i) = elx(m*i+j);
	CalcFlux (fel, bmip, selx, sflux, applyd, lh);
	for (int i = 0; i < sflux.Size(); i++)
	  flux(m*i+j) = sflux(i);
      }
  }

  void BilinearFormIntegrator :: 
  CalcFluxMulti (const FiniteElement & fel,
		 const BaseMappedIntegrationPoint & bmip,
		 int n,
                 FlatVector<Complex> elx, 
		 FlatVector<Complex> flux,
		 bool applyd,
		 LocalHeap & lh) const
  {
    cerr << "calcflux<Complex> for Specific called for base class " 
	 << typeid(*this).name()
	 << endl;
  }














  void BilinearFormIntegrator :: 
  ApplyBTrans (const FiniteElement & fel,
	       // const ElementTransformation & eltrans,
	       // const IntegrationPoint & ip,
	       const BaseMappedIntegrationPoint & bmip,
               FlatVector<double> elx, 
	       FlatVector<double> ely,
	       LocalHeap & lh) const
  {
    cerr << "ApplyBTrans<double> called for class " 
	 << typeid(*this).name()
	 << endl;
  }

  void BilinearFormIntegrator :: 
  ApplyBTrans (const FiniteElement & fel,
	       const BaseMappedIntegrationPoint & bmip,
               FlatVector<Complex> elx, 
	       FlatVector<Complex> ely,
	       LocalHeap & lh) const
  {
    cerr << "ApplyBTrans<Complex> called for class " 
	 << typeid(*this).name()
	 << endl;
  }

  void BilinearFormIntegrator :: 
  ApplyBTrans (const FiniteElement & fel,
	       const BaseMappedIntegrationRule & mir,
               FlatMatrix<double> elx, 
	       FlatVector<double> ely,
	       LocalHeap & lh) const
  {
    FlatVector<> hv(ely.Size(), lh);
    ely = 0.0;
    for (int l = 0; l < mir.Size(); l++)
      {
	FlatVector<> elxi = elx.Row(l);
	ApplyBTrans (fel, mir[l], elxi, hv, lh);
	ely += hv;
      }
  }
  
  void BilinearFormIntegrator :: 
  ApplyBTrans (const FiniteElement & fel,
	       const BaseMappedIntegrationRule & mir,
               FlatMatrix<Complex> elx, 
	       FlatVector<Complex> ely,
	       LocalHeap & lh) const
  {
    FlatVector<Complex> hv(ely.Size(), lh);
    ely = 0.0;
    for (int l = 0; l < mir.Size(); l++)
      {
	FlatVector<Complex> elxi = elx.Row(l);
	ApplyBTrans (fel, mir[l], elxi, hv, lh);
	ely += hv;
      }
  }
  
  



  
  void BilinearFormIntegrator :: 
  ApplyDMat (const FiniteElement & bfel,
	     const BaseMappedIntegrationPoint & bmip,
             FlatVector<double> elx, 
	     FlatVector<double> eldx,
	     LocalHeap & lh) const
  {
    cerr << "ApplyDMat<double> called for class " 
	 << typeid(*this).name()
	 << endl;
  }

  void BilinearFormIntegrator :: 
  ApplyDMat (const FiniteElement & bfel,
	     const BaseMappedIntegrationPoint & bmip,
             FlatVector<Complex> elx, 
	     FlatVector<Complex> eldx,
	     LocalHeap & lh) const
  {
    cerr << "ApplyDMat<Complex> called for class " 
	 << typeid(*this).name()
	 << endl;

  }

  
  void BilinearFormIntegrator :: 
  ApplyDMat (const FiniteElement & bfel,
	     const BaseMappedIntegrationRule & mir,
             FlatMatrix<double> elx, 
	     FlatMatrix<double> eldx,
	     LocalHeap & lh) const
  {
    cerr << "ApplyDMat<double>, MappedIR called for class " 
	 << typeid(*this).name()
	 << endl;

  }

  void BilinearFormIntegrator :: 
  ApplyDMat (const FiniteElement & bfel,
	     const BaseMappedIntegrationRule & mir,
             FlatMatrix<Complex> elx, 
	     FlatMatrix<Complex> eldx,
	     LocalHeap & lh) const
  {
    cerr << "ApplyDMat<Complex>, MappedIR called for class " 
	 << typeid(*this).name()
	 << endl;
  }

  
  void BilinearFormIntegrator :: 
  ApplyDMatInv (const FiniteElement & bfel,
                const BaseMappedIntegrationPoint & bmip,
                FlatVector<double> elx, 
                FlatVector<double> eldx,
                LocalHeap & lh) const
  {
    cerr << "ApplyDMatInv<double> called for class " 
	 << typeid(*this).name()
	 << endl;
  }

  void BilinearFormIntegrator :: 
  ApplyDMatInv (const FiniteElement & bfel,
                const BaseMappedIntegrationPoint & bmip,
                FlatVector<Complex> elx, 
                FlatVector<Complex> eldx,
                LocalHeap & lh) const
  {
    cerr << "ApplyDMatInv<Complex> called for class " 
	 << typeid(*this).name()
	 << endl;

  }

  
  void BilinearFormIntegrator :: 
  ApplyDMatInv (const FiniteElement & bfel,
		const BaseMappedIntegrationRule & mir,
                FlatMatrix<double> elx, 
		FlatMatrix<double> eldx,
		LocalHeap & lh) const
  {
    cerr << "ApplyDMatInv<double>, MappedIR called for class " 
	 << typeid(*this).name()
	 << endl;

  }

  void BilinearFormIntegrator :: 
  ApplyDMatInv (const FiniteElement & bfel,
		const BaseMappedIntegrationRule & mir,
                FlatMatrix<Complex> elx, 
		FlatMatrix<Complex> eldx,
		LocalHeap & lh) const
  {
    cerr << "ApplyDMatInv<Complex>, MappedIR called for class " 
	 << typeid(*this).name()
	 << endl;
  }

  









  /*
  const IntegrationRule &  BilinearFormIntegrator :: 
  GetIntegrationRule (const FiniteElement & fel,
                      const bool use_higher_integration_order) const
  {
    cerr << "GetIntegrationRule  called for class " 
	 << typeid(*this).name()
	 << endl;
    static IntegrationRule dummy; return dummy;
  }
  */







  BlockBilinearFormIntegrator :: 
  BlockBilinearFormIntegrator (shared_ptr<BilinearFormIntegrator> abfi, int adim, int acomp)
    : bfi(abfi), dim(adim), comp(acomp) 
  { 
    ;
  }

  BlockBilinearFormIntegrator ::
  BlockBilinearFormIntegrator (shared_ptr<BilinearFormIntegrator> abfi, int adim)
    : bfi(abfi), dim(adim), comp(-1)
  {
    ; 
  }

  BlockBilinearFormIntegrator :: 
  ~BlockBilinearFormIntegrator ()
  {
    ;
  }

  void BlockBilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & bfel,
		     const ElementTransformation & eltrans,
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    FlatMatrix<double> mat1(bfel.GetNDof(), lh);
    bfi->CalcElementMatrix (bfel, eltrans, mat1, lh);
    
    elmat = 0;

    if (comp == -1)
      for (int i = 0; i < mat1.Height(); i++)
	for (int j = 0; j < mat1.Width(); j++)
	  for (int k = 0; k < dim; k++)
	    elmat(dim*i+k, dim*j+k) = mat1(i,j);
    else
      for (int i = 0; i < mat1.Height(); i++)
	for (int j = 0; j < mat1.Width(); j++)
	  elmat(dim*i+comp, dim*j+comp) = mat1(i,j);
  }  




  void BlockBilinearFormIntegrator :: 
  CalcElementMatrix (const FiniteElement & bfel, 
		     const ElementTransformation & eltrans, 
		     FlatMatrix<Complex> elmat,
		     LocalHeap & lh) const
  {
    FlatMatrix<Complex> mat1(bfel.GetNDof(), lh);
    bfi->CalcElementMatrix (bfel, eltrans, mat1, lh);
    
    elmat = 0;

    if (comp == -1)
      for (int i = 0; i < mat1.Height(); i++)
	for (int j = 0; j < mat1.Width(); j++)
	  for (int k = 0; k < dim; k++)
	    elmat(dim*i+k, dim*j+k) = mat1(i,j);
    else
      for (int i = 0; i < mat1.Height(); i++)
	for (int j = 0; j < mat1.Width(); j++)
	  elmat(dim*i+comp, dim*j+comp) = mat1(i,j);
     }  


  
  void BlockBilinearFormIntegrator ::
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
                      FlatVector<double> elx, 
		      FlatVector<double> ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    int smallsizex = elx.Size()/dim;
    int smallsizey = ely.Size()/dim;

    Vector<double> small_elx(smallsizex);
    Vector<double> small_ely(smallsizey);

    ely = 0;
    if(comp == -1)
      {
	for(int d = 0; d < dim; d++)
	  {
            small_elx = elx.Slice(d, dim);
	    // for(int i=0; i<smallsizex; i++)
            // small_elx(i) = elx(i*dim+d);

	    bfi->ApplyElementMatrix(bfel,eltrans,small_elx,small_ely, 0, lh);
            ely.Slice(d, dim) = small_ely;

	    // for(int i=0; i<smallsizey; i++)
            // ely(i*dim+d) = small_ely(i);
	  }
      }
    else
      {
	for(int i=0; i<smallsizex; i++)
	  {
	    small_elx(i) = elx(i*dim+comp);
	  }
	bfi->ApplyElementMatrix(bfel,eltrans,small_elx,small_ely,0, lh);
	for(int i=0; i<smallsizey; i++)
	  {
	    ely(i*dim+comp) = small_ely(i);
	  }
      }
  }


  void BlockBilinearFormIntegrator ::
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
                      FlatVector<Complex> elx, 
		      FlatVector<Complex> ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    const int smallsizex = elx.Size()/dim;
    const int smallsizey = ely.Size()/dim;

    Vector<Complex> small_elx(smallsizex);
    Vector<Complex> small_ely(smallsizey);

    ely = 0;
    if(comp == -1)
      {
	for(int d=0; d<dim; d++)
	  {
	    for(int i=0; i<smallsizex; i++)
	      {
		small_elx(i) = elx(i*dim+d);
	      }
	    bfi->ApplyElementMatrix(bfel,eltrans,small_elx,small_ely,0, lh);
	    for(int i=0; i<smallsizey; i++)
	      {
		ely(i*dim+d) = small_ely(i);
	      }
	  }
      }
    else
      {
	for(int i=0; i<smallsizex; i++)
	  {
	    small_elx(i) = elx(i*dim+comp);
	  }
	bfi->ApplyElementMatrix(bfel,eltrans,small_elx,small_ely,0, lh);
	for(int i=0; i<smallsizey; i++)
	  {
	    ely(i*dim+comp) = small_ely(i);
	  }
      }
  }


  void BlockBilinearFormIntegrator ::
  CalcLinearizedElementMatrix (const FiniteElement & bfel,
			       const ElementTransformation & eltrans,
			       FlatVector<double> elveclin,
			       FlatMatrix<double> elmat,
			       LocalHeap & lh) const
  {
    FlatMatrix<double> mat1(bfel.GetNDof(), lh);
    const int smallsizelin = elveclin.Size()/dim;

    FlatVector<double> small_elveclin(smallsizelin, lh);

    if (comp == -1)
      {
	for(int d=0; d<dim; d++)
	  {
	    for(int i=0; i<smallsizelin; i++)
              small_elveclin(i) = elveclin(i*dim+d);

	    bfi->CalcLinearizedElementMatrix (bfel, eltrans, small_elveclin, mat1, lh);
            elmat = 0;

	    for (int i = 0; i < mat1.Height(); i++)
	      for (int j = 0; j < mat1.Width(); j++)
		elmat(i*dim+d, j*dim+d) = mat1(i,j);
	  }
      }
    else
      {
	for(int i=0; i<smallsizelin; i++)
          small_elveclin(i) = elveclin(i*dim+comp);

	bfi->CalcLinearizedElementMatrix (bfel, eltrans, small_elveclin, mat1, lh);
	elmat = 0;
	for (int i = 0; i < mat1.Height(); i++)
	  for (int j = 0; j < mat1.Width(); j++)
	    elmat(dim*i+comp, dim*j+comp) = mat1(i,j);
      }  
  }




  void BlockBilinearFormIntegrator :: 
  CalcLinearizedElementMatrix (const FiniteElement & bfel, 
			       const ElementTransformation & eltrans, 
			       FlatVector<Complex> elveclin,
			       FlatMatrix<Complex> elmat,
			       LocalHeap & lh) const
  {
    throw Exception ("old style matrix allocation in BlockBilinearFormIntegrator :: CalcLinearizedElementMatrix");

    const int smallsizelin = elveclin.Size()/dim;

    Vector<Complex> small_elveclin(smallsizelin);

    FlatMatrix<Complex> mat1;

    if (comp == -1)
      {
	bool allocated(false);

	for(int d=0; d<dim; d++)
	  {
	    for(int i=0; i<smallsizelin; i++)
	      {
		small_elveclin(i) = elveclin(i*dim+d);
	      }
	    bfi->CalcLinearizedElementMatrix (bfel, eltrans, small_elveclin, mat1, lh);
	    if(!allocated)
	      {
		elmat.AssignMemory (mat1.Height()*dim, mat1.Width()*dim, lh);
		elmat = 0;
		allocated = true;
	      }
	    for (int i = 0; i < mat1.Height(); i++)
	      for (int j = 0; j < mat1.Width(); j++)
		elmat(i*dim+d, j*dim+d) = mat1(i,j);
	  }
      }
    else
      {
	for(int i=0; i<smallsizelin; i++)
	  {
	    small_elveclin(i) = elveclin(i*dim+comp);
	  }
	bfi->CalcLinearizedElementMatrix (bfel, eltrans, small_elveclin, mat1, lh);
	elmat.AssignMemory (mat1.Height()*dim, mat1.Width()*dim, lh);
	elmat = 0;
	for (int i = 0; i < mat1.Height(); i++)
	  for (int j = 0; j < mat1.Width(); j++)
	    elmat(dim*i+comp, dim*j+comp) = mat1(i,j);
      }  
  }  


#ifdef OLD_REMOVED_FOR_CLANG
  void  BlockBilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & fel,
	    const ElementTransformation & eltrans,
	    const IntegrationPoint & ip,
            FlatVector<double> elx, 
	    FlatVector<double> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    if (comp >= 0)
      {
        FlatVector<double> selx = flux.Slice(comp, dim) | lh; 
        /*
	FlatVector<double> selx(elx.Size()/dim, lh);
	for (int i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+comp);
        */
	bfi->CalcFlux (fel, eltrans(ip, lh), selx, flux, applyd, lh);
      }
    else
      {
	FlatVector<double> selx(elx.Size()/dim, lh);
	FlatVector<double> sflux(bfi->DimFlux(), lh);
	// flux.AssignMemory (DimFlux(), lh);
	for (int j = 0; j < dim; j++)
	  {
	    for (int i = 0; i < selx.Size(); i++)
	      selx(i) = elx(dim*i+j);
	    bfi->CalcFlux (fel, eltrans(ip, lh), selx, sflux, applyd, lh);
	    for (int i = 0; i < sflux.Size(); i++)
	      flux(dim*i+j) = sflux(i);
	  }
      }
    /*
    int i, j;
    int mincomp = 0;
    int maxcomp = dim-1;
    if (comp >= 0)
      mincomp = maxcomp = comp;

    FlatVector<double> selx(elx.Size()/dim, lh);
    FlatVector<double> sflux(bfi.DimFlux(), lh);
    flux.AssignMemory (DimFlux(), lh);
    for (j = mincomp; j <= maxcomp; j++)
      {
	for (i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+j);
	bfi.CalcFlux (fel, eltrans, ip, selx, sflux, applyd, lh);
	for (i = 0; i < sflux.Size(); i++)
	  flux(dim*i+j) = sflux(i);
      }
    */
  }

  void  BlockBilinearFormIntegrator ::
  CalcFlux (const FiniteElement & fel,
	    const ElementTransformation & eltrans,
	    const IntegrationPoint & ip,
            FlatVector<Complex> elx, 
	    FlatVector<Complex> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    if (comp >= 0)
      {
	FlatVector<Complex> selx(elx.Size()/dim, lh);
	for (int i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+comp);
	bfi->CalcFlux (fel, eltrans(ip, lh), selx, flux, applyd, lh);
      }
    else
      {
	FlatVector<Complex> selx(elx.Size()/dim, lh);
	FlatVector<Complex> sflux(bfi->DimFlux(), lh);
	// flux.AssignMemory (DimFlux(), lh);
	for (int j = 0; j < dim; j++)
	  {
	    for (int i = 0; i < selx.Size(); i++)
	      selx(i) = elx(dim*i+j);
	    bfi->CalcFlux (fel, eltrans(ip, lh), selx, sflux, applyd, lh);
	    for (int i = 0; i < sflux.Size(); i++)
	      flux(dim*i+j) = sflux(i);
	  }
      }

    /*
    int i, j;
    int mincomp = 0;
    int maxcomp = dim-1;
    if (comp >= 0)
      mincomp = maxcomp = comp;
    
    FlatVector<Complex> selx(elx.Size()/dim, lh);
    FlatVector<Complex> sflux(bfi.DimFlux(), lh);
    flux.AssignMemory (DimFlux(), lh);
    for (j = mincomp; j <= maxcomp; j++)
      {
	for (i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+j);
	bfi.CalcFlux (fel, eltrans, ip, selx, sflux, applyd, lh);
	for (i = 0; i < sflux.Size(); i++)
	  flux(dim*i+j) = sflux(i);
      }
    */
  }
#endif


  void  BlockBilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationPoint & bmip,
            BareSliceVector<double> elx, 
	    FlatVector<double> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    if (comp >= 0)
      {
        /*
	FlatVector<double> selx(elx.Size()/dim, lh);
	for (int i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+comp);
	bfi->CalcFlux (fel, bmip, selx, flux, applyd, lh);
        */
        bfi->CalcFlux (fel, bmip, elx.Slice(comp, dim), flux, applyd, lh);
      }
    else
      {
        // only possible for elx.Dist() = 1
        if (elx.Dist() == 1)
          {
            FlatVector<> felx(dim*fel.GetNDof(), &elx(0));
            bfi->CalcFluxMulti (fel, bmip, dim, felx, flux, applyd, lh);
            return;
          }
        
	// FlatVector<double> selx(elx.Size()/dim, lh);
	FlatVector<double> sflux(bfi->DimFlux(), lh);
	// flux.AssignMemory (DimFlux(), lh);
	for (int j = 0; j < dim; j++)
	  {
            /*
	    for (int i = 0; i < selx.Size(); i++)
	      selx(i) = elx(dim*i+j);
	    bfi.CalcFlux (fel, bmip, selx, sflux, applyd, lh);
            */
	    bfi->CalcFlux (fel, bmip, elx.Slice(j,dim), sflux, applyd, lh);            
	    for (int i = 0; i < sflux.Size(); i++)
	      flux(dim*i+j) = sflux(i);
	  }
      }


    /*
    int i, j;
    int mincomp = 0;
    int maxcomp = dim-1;
    if (comp >= 0)
      mincomp = maxcomp = comp;


    FlatVector<double> selx(elx.Size()/dim, lh);
    FlatVector<double> sflux(bfi.DimFlux(), lh);
    flux.AssignMemory (DimFlux(), lh);
    for (j = mincomp; j <= maxcomp; j++)
      {
	for (i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+j);
	bfi.CalcFlux (fel, bmip, selx, sflux, applyd, lh);
	for (i = 0; i < sflux.Size(); i++)
	  flux(dim*i+j) = sflux(i);
      }
    */
  }

  void  BlockBilinearFormIntegrator ::
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationPoint & bmip,
            BareSliceVector<Complex> elx, 
	    FlatVector<Complex> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    if (comp >= 0)
      {
        /*
	FlatVector<Complex> selx(elx.Size()/dim, lh);
	for (int i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+comp);
	bfi->CalcFlux (fel, bmip, selx, flux, applyd, lh);
        */
	bfi->CalcFlux (fel, bmip, elx.Slice(comp,dim), flux, applyd, lh);        
      }
    else
      {
	// FlatVector<Complex> selx(elx.Size()/dim, lh);
	FlatVector<Complex> sflux(bfi->DimFlux(), lh);
	// flux.AssignMemory (DimFlux(), lh);
	for (int j = 0; j < dim; j++)
	  {
            /*
	    for (int i = 0; i < selx.Size(); i++)
	      selx(i) = elx(dim*i+j);
	    bfi->CalcFlux (fel, bmip, selx, sflux, applyd, lh);
            */
	    bfi->CalcFlux (fel, bmip, elx.Slice(j,dim), sflux, applyd, lh);            
	    for (int i = 0; i < sflux.Size(); i++)
	      flux(dim*i+j) = sflux(i);
	  }
      }

    /*
    int i, j;
    int mincomp = 0;
    int maxcomp = dim-1;
    if (comp >= 0)
      mincomp = maxcomp = comp;
    
    FlatVector<Complex> selx(elx.Size()/dim, lh);
    FlatVector<Complex> sflux(bfi.DimFlux(), lh);
    flux.AssignMemory (DimFlux(), lh);
    for (j = mincomp; j <= maxcomp; j++)
      {
	for (i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+j);
	bfi.CalcFlux (fel, bmip, selx, sflux, applyd, lh);
	for (i = 0; i < sflux.Size(); i++)
	  flux(dim*i+j) = sflux(i);
      }
    */
  }


  

  void  BlockBilinearFormIntegrator ::
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationRule & mir,
            BareSliceVector<double> elx, 
	    BareSliceMatrix<double> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    int mincomp = 0;
    int maxcomp = dim-1;
    if (comp >= 0)
      mincomp = maxcomp = comp;
    
    // FlatVector<> selx(elx.Size()/dim, lh);
    FlatMatrix<> sflux(mir.Size(), bfi->DimFlux(), lh);

    for (int j = mincomp; j <= maxcomp; j++)
      {
        /*
	for (int i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+j);
	bfi->CalcFlux (fel, mir, selx, sflux, applyd, lh);
        */
	bfi->CalcFlux (fel, mir, elx.Slice(j,dim), sflux, applyd, lh);        
	for (int k = 0; k < mir.Size(); k++)
	  for (int i = 0; i < sflux.Width(); i++)
	    flux.Row(k)(dim*i+j) = sflux.Row(k)(i);
      }
  }






  double BlockBilinearFormIntegrator :: 
  Energy (const FiniteElement & fel, 
	  const ElementTransformation & eltrans, 
          FlatVector<double> elx, 
	  LocalHeap & lh) const
  {
    int i, j;
    int mincomp = 0;
    int maxcomp = dim-1;
    if (comp >= 0)
      mincomp = maxcomp = comp;

    double energy = 0;
    FlatVector<double> selx(elx.Size()/dim, lh);
    for (j = mincomp; j <= maxcomp; j++)
      {
	for (i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+j);
	energy += bfi->Energy (fel, eltrans, selx, lh);
      }
    return energy;
  }




  void  BlockBilinearFormIntegrator ::
  ApplyBTrans (const FiniteElement & fel,
	       const BaseMappedIntegrationPoint & bmip,
               FlatVector<double> elx, 
	       FlatVector<double> ely,
	       LocalHeap & lh) const
  {
    int mincomp = 0;
    int maxcomp = dim-1;
    if (comp >= 0)
      mincomp = maxcomp = comp;

    FlatVector<double> selx(elx.Size()/dim, lh);
    FlatVector<double> sely(ely.Size()/dim, lh);

    for (int j = mincomp; j <= maxcomp; j++)
      {
	/*
	for (int i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+j);
	*/
	selx = elx.Slice(j, dim);
	bfi->ApplyBTrans (fel, bmip, selx, sely, lh);
	ely.Slice(j, dim) = sely;
	/*
	for (int i = 0; i < sely.Size(); i++)
	  ely(dim*i+j) = sely(i);
	*/
      }
  }




  void  BlockBilinearFormIntegrator ::
  ApplyBTrans (const FiniteElement & fel,
	       const BaseMappedIntegrationPoint & bmip,
               FlatVector<Complex> elx, 
	       FlatVector<Complex> ely,
	       LocalHeap & lh) const
  {
    int i, j;
    int mincomp = 0;
    int maxcomp = dim-1;
    if (comp >= 0)
      mincomp = maxcomp = comp;
    
    FlatVector<Complex> selx(elx.Size()/dim, lh);
    FlatVector<Complex> sely(ely.Size()/dim, lh);
    for (j = mincomp; j <= maxcomp; j++)
      {
	for (i = 0; i < selx.Size(); i++)
	  selx(i) = elx(dim*i+j);
	bfi->ApplyBTrans (fel, bmip, selx, sely, lh);
	for (i = 0; i < sely.Size(); i++)
	  ely(dim*i+j) = sely(i);
      }
  }


  void  BlockBilinearFormIntegrator ::
  ApplyBTrans (const FiniteElement & fel,
	       const BaseMappedIntegrationRule & mir,
               FlatMatrix<double> elx, 
	       FlatVector<double> ely,
	       LocalHeap & lh) const
  {
    int mincomp = 0;
    int maxcomp = dim-1;
    if (comp >= 0) mincomp = maxcomp = comp;

    FlatMatrix<double> selx(mir.Size(), elx.Width()/dim, lh);
    FlatVector<double> sely(ely.Size()/dim, lh);

    if (comp >= 0) ely = 0.0;

    for (int j = mincomp; j <= maxcomp; j++)
      {
	for (int i = 0; i < selx.Width(); i++)
	  selx.Col(i) = elx.Col(dim*i+j);

	bfi->ApplyBTrans (fel, mir, selx, sely, lh);

        ely.Slice(j, dim) = sely;
      }
  }
  

  string BlockBilinearFormIntegrator :: Name () const
  {
    return
      string ("BlockIntegrator (") + bfi->Name() +
      string (")");
  }






  
  ComplexBilinearFormIntegrator :: 
  ComplexBilinearFormIntegrator (shared_ptr<BilinearFormIntegrator> abfi, 
				 Complex afactor)
    : bfi(abfi), factor(afactor)
  {
    ;
  }
  

  void ComplexBilinearFormIntegrator :: 
  CalcElementMatrix (const FiniteElement & fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    throw Exception ("ComplexBilinearFormIntegrator: cannot assemble double matrix");
  }



  void ComplexBilinearFormIntegrator :: 
  CalcElementMatrix (const FiniteElement & fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<Complex> elmat,
		     LocalHeap & lh) const
  {
    FlatMatrix<double> rmat(elmat.Height(), lh);
    bfi->CalcElementMatrix (fel, eltrans, rmat, lh);
    elmat = factor * rmat; 
  }  


  void ComplexBilinearFormIntegrator :: 
  ApplyElementMatrix (const FiniteElement & fel, 
		      const ElementTransformation & eltrans, 
                      FlatVector<Complex> elx, 
		      FlatVector<Complex> ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    bfi->ApplyElementMatrix (fel, eltrans, elx, ely, 0, lh);
    ely *= factor;
  }


  
  void ComplexBilinearFormIntegrator :: 
  CalcElementMatrixIndependent (const FiniteElement & bfel_master,
				    const FiniteElement & bfel_master_element,				    
				    const FiniteElement & bfel_other,
				    const ElementTransformation & eltrans_master, 
				    const ElementTransformation & eltrans_master_element, 
				    const ElementTransformation & eltrans_other,
				    const IntegrationPoint & ip_master,
				    const IntegrationPoint & ip_master_element,
				    const IntegrationPoint & ip_other,
				    FlatMatrix<double> & elmat,
				    LocalHeap & lh) const
  {
    throw Exception ("ComplexBilinearFormIntegrator: cannot assemble double matrix");
  }

  

  void ComplexBilinearFormIntegrator :: 
  CalcElementMatrixIndependent (const FiniteElement & bfel_master,
				    const FiniteElement & bfel_master_element,				    
				    const FiniteElement & bfel_other,
				    const ElementTransformation & eltrans_master, 
				    const ElementTransformation & eltrans_master_element, 
				    const ElementTransformation & eltrans_other,
				    const IntegrationPoint & ip_master,
				    const IntegrationPoint & ip_master_element,
				    const IntegrationPoint & ip_other,
				    FlatMatrix<Complex> & elmat,
				    LocalHeap & lh) const
  {
    FlatMatrix<double> rmat;
    bfi->CalcElementMatrixIndependent(bfel_master,bfel_master_element,bfel_other,
                                     eltrans_master, eltrans_master_element, eltrans_other,
                                     ip_master, ip_master_element, ip_other,
                                     rmat, lh);
    elmat.AssignMemory(rmat.Height(), rmat.Width(), lh);
    elmat = factor * rmat;
  }

  

  void ComplexBilinearFormIntegrator :: 
  CalcElementMatrixIndependent (const FiniteElement & bfel_master,
                                const FiniteElement & bfel_other,
                                const ElementTransformation & eltrans_master, 
                                const ElementTransformation & eltrans_other,
                                const IntegrationPoint & ip_master,
                                const IntegrationPoint & ip_other,
                                FlatMatrix<double> & elmat,
                                LocalHeap & lh) const
  {
    throw Exception ("ComplexBilinearFormIntegrator: cannot assemble double matrix");
  }



  void ComplexBilinearFormIntegrator :: 
  CalcElementMatrixIndependent (const FiniteElement & bfel_master,
				    const FiniteElement & bfel_other,
                                const ElementTransformation & eltrans_master, 
                                const ElementTransformation & eltrans_other,
                                const IntegrationPoint & ip_master,
                                const IntegrationPoint & ip_other,
                                FlatMatrix<Complex> & elmat,
                                LocalHeap & lh) const
  {
    FlatMatrix<double> rmat;
    bfi->CalcElementMatrixIndependent(bfel_master,bfel_other,
					 eltrans_master, eltrans_other,
					 ip_master, ip_other,
					 rmat, lh);
    elmat.AssignMemory(rmat.Height(), rmat.Width(), lh);
    elmat = factor * rmat;
  }



  
  /*
  void ComplexBilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & fel,
	    const ElementTransformation & eltrans,
	    const IntegrationPoint & ip,
            FlatVector<Complex> elx, 
	    FlatVector<Complex> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    bfi->CalcFlux (fel, eltrans(ip, lh), elx, flux, applyd, lh);
    flux *= factor;
  }
  */

  void ComplexBilinearFormIntegrator ::
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationPoint & bmip,
            BareSliceVector<Complex> elx, 
	    FlatVector<Complex> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    bfi->CalcFlux (fel, bmip, elx, flux, applyd, lh);
    flux *= factor;
  }

  string ComplexBilinearFormIntegrator :: Name () const
  {
    return
      string ("ComplexIntegrator (") + bfi->Name() +
      string (")");
  }



  /*
  const IntegrationRule & ComplexBilinearFormIntegrator :: GetIntegrationRule (const FiniteElement & fel,
									       const bool use_higher_integration_order) const
  {
    return bfi->GetIntegrationRule(fel,use_higher_integration_order);
  }
  */

  void TransposeBilinearFormIntegrator :: 
  CalcElementMatrix (const FiniteElement & bfel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<double> temp (bfel.GetNDof(), lh);
    bfi->CalcElementMatrix (bfel, eltrans, temp, lh);
    elmat = Trans (temp);
  }



  void BilinearFormIntegratorAnyDim ::  
  CheckElement (const FiniteElement & el) const
  {
    ;
  }

  void BilinearFormIntegratorAnyDim ::  
  CalcElementMatrix (const FiniteElement & bfel, 
                     const ElementTransformation & eltrans, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    int dim = eltrans.SpaceDim();
    if (bfi[dim])
      bfi[dim] -> CalcElementMatrix(bfel, eltrans, elmat, lh);
    else
      throw Exception (ToString("Integrator-Anydim not available for dimension ")+
                       ToString(dim));
  }

  void BilinearFormIntegratorAnyDim ::  
  CalcElementMatrix (const FiniteElement & bfel, 
                     const ElementTransformation & eltrans, 
                     FlatMatrix<Complex> elmat,
                     LocalHeap & lh) const
  {
    int dim = eltrans.SpaceDim();
    if (bfi[dim])
      bfi[dim] -> CalcElementMatrix(bfel, eltrans, elmat, lh);
    else
      throw Exception (ToString("Integrator-Anydim not available for dimension ")+
                       ToString(dim));
  }

  



  CompoundBilinearFormIntegrator :: 
  CompoundBilinearFormIntegrator (shared_ptr<BilinearFormIntegrator> abfi, int acomp)
    : bfi(abfi), comp(acomp)
  {
    auto & comp_evaluators = bfi->GetEvaluators();
    for (size_t i = 0; i < comp_evaluators.Size(); i++)
      evaluators.Set(comp_evaluators.GetName(i),
                     make_shared<CompoundDifferentialOperator>(comp_evaluators[i], comp));
  }
  

  void CompoundBilinearFormIntegrator :: 
  CalcElementMatrix (const FiniteElement & bfel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    FlatMatrix<double> mat1(fel[comp].GetNDof(), lh);
    bfi->CalcElementMatrix (fel[comp], eltrans, mat1, lh);

    elmat = 0;

    IntRange range = fel.GetRange (comp);
    elmat.Rows (range).Cols(range) = mat1;
  }  


  void CompoundBilinearFormIntegrator :: 
  CalcElementMatrix (const FiniteElement & bfel,
			 const ElementTransformation & eltrans, 
			 FlatMatrix<Complex> elmat,
			 LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    FlatMatrix<Complex> mat1(fel[comp].GetNDof(), lh);
    bfi->CalcElementMatrix (fel[comp], eltrans, mat1, lh);

    elmat = 0;

    IntRange range = fel.GetRange (comp);
    elmat.Rows (range).Cols(range) = mat1;
    /*
    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();

    elmat.Rows (base, base+mat1.Height()).Cols (base, base+mat1.Width()) = mat1;
    */
  }  




  void CompoundBilinearFormIntegrator :: 
  CalcLinearizedElementMatrix (const FiniteElement & bfel, 
			       const ElementTransformation & eltrans, 
			       FlatVector<double> ellin,
			       FlatMatrix<double> elmat,
			       LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    int nd = fel[comp].GetNDof();
    FlatMatrix<double> mat1(nd, lh);
    FlatVector<double> ellin1(nd, lh);

    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();
    
    for (int j = 0; j < nd; j++)
      ellin1(j) = ellin(base+j);

    bfi->CalcLinearizedElementMatrix (fel[comp], eltrans, ellin1, mat1, lh);
    
    elmat = 0;
    for (int i = 0; i < mat1.Height(); i++)
      for (int j = 0; j < mat1.Width(); j++)
	elmat(base+i, base+j) = mat1(i,j);
  }  


  void CompoundBilinearFormIntegrator :: 
  CalcLinearizedElementMatrix (const FiniteElement & bfel, 
			       const ElementTransformation & eltrans, 
			       FlatVector<Complex> ellin,
			       FlatMatrix<Complex> elmat,
			       LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      dynamic_cast<const CompoundFiniteElement&> (bfel);


    int nd = fel[comp].GetNDof();

    FlatMatrix<Complex> mat1(nd, lh);
    FlatVector<Complex> ellin1(nd, lh);

    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();
    
    for (int j = 0; j < nd; j++)
      ellin1(j) = ellin(base+j);

    bfi->CalcLinearizedElementMatrix (fel[comp], eltrans, ellin1, mat1, lh);
    
    elmat = 0;
    
    for (int i = 0; i < mat1.Height(); i++)
      for (int j = 0; j < mat1.Width(); j++)
	elmat(base+i, base+j) = mat1(i,j);
  }  




  void CompoundBilinearFormIntegrator :: 
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
                      FlatVector<double> elx,
		      FlatVector<double> ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    int nd = fel[comp].GetNDof();
    FlatVector<double> elx1(nd, lh);
    FlatVector<double> ely1(nd, lh);

    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();
    
    for (int j = 0; j < nd; j++)
      elx1(j) = elx(base+j);

    bfi->ApplyElementMatrix (fel[comp], eltrans, elx1, ely1, 0, lh);
    
    ely = 0;
    for (int j = 0; j < nd; j++)
      ely(base+j) = ely1(j);

  }  





  




  void CompoundBilinearFormIntegrator :: 
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
                      FlatVector<Complex> elx,
		      FlatVector<Complex> ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      static_cast<const CompoundFiniteElement&> (bfel);

    int nd = fel[comp].GetNDof();
    FlatVector<Complex> elx1(nd, lh);
    FlatVector<Complex> ely1(nd, lh);

    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();
    

    for (int j = 0; j < nd; j++)
      elx1(j) = elx(base+j);

    bfi->ApplyElementMatrix (fel[comp], eltrans, elx1, ely1, 0, lh);
    
    ely = 0;
    for (int j = 0; j < nd; j++)
      ely(base+j) = ely1(j);
  }  






  void CompoundBilinearFormIntegrator :: 
  ApplyLinearizedElementMatrix (const FiniteElement & bfel, 
				const ElementTransformation & eltrans, 
                                FlatVector<double> ellin,
                                FlatVector<double> elx,
				FlatVector<double> ely,
				LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      static_cast<const CompoundFiniteElement&> (bfel);

    int nd = fel[comp].GetNDof();
    FlatVector<double> ellin1(nd, lh);
    FlatVector<double> elx1(nd, lh);
    FlatVector<double> ely1(nd, lh);

    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();
    
    for (int j = 0; j < nd; j++)
      {
	ellin1(j) = ellin(base+j);
	elx1(j) = elx(base+j);
      }

    bfi->ApplyLinearizedElementMatrix (fel[comp], eltrans, ellin1, 
				      elx1, ely1, lh);
    
    ely = 0;
    for (int j = 0; j < nd; j++)
      ely(base+j) = ely1(j);
  }  


  void CompoundBilinearFormIntegrator :: 
  ApplyLinearizedElementMatrix (const FiniteElement & bfel, 
				const ElementTransformation & eltrans, 
                                FlatVector<Complex> ellin,
                                FlatVector<Complex> elx,
				FlatVector<Complex> ely,
				LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      static_cast<const CompoundFiniteElement&> (bfel);

    int nd = fel[comp].GetNDof();
    FlatVector<Complex> ellin1(nd, lh);
    FlatVector<Complex> elx1(nd, lh);
    FlatVector<Complex> ely1(nd, lh);

    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();
    
    for (int j = 0; j < nd; j++)
      {
	ellin1(j) = ellin(base+j);
	elx1(j) = elx(base+j);
      }

    bfi->ApplyLinearizedElementMatrix (fel[comp], eltrans, ellin1,
				      elx1, ely1, lh);
    
    ely = 0;
    for (int j = 0; j < nd; j++)
      ely(base+j) = ely1(j);
  }  



  void CompoundBilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & bfel,
	    const BaseMappedIntegrationPoint & ip,
            BareSliceVector<double> elx, 
	    FlatVector<double> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      static_cast<const CompoundFiniteElement&> (bfel);
    
    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();

    bfi->CalcFlux (fel[comp], ip, elx.Range (base, base+fel[comp].GetNDof()), 
		  flux, applyd, lh);
  }


  void CompoundBilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & bfel,
	    const BaseMappedIntegrationPoint & ip,
            BareSliceVector<Complex> elx, 
	    FlatVector<Complex> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      static_cast<const CompoundFiniteElement&> (bfel);

    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();

    bfi->CalcFlux (fel[comp], ip, elx.Range (base, base+fel[comp].GetNDof()), 
		  flux, applyd, lh);
  }


  void CompoundBilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & bfel,
	    const BaseMappedIntegrationRule & mir,
            BareSliceVector<double> elx, 
	    BareSliceMatrix<double> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      static_cast<const CompoundFiniteElement&> (bfel);
    
    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();

    bfi->CalcFlux (fel[comp], mir, elx.Range (base, base+fel[comp].GetNDof()), 
		  flux, applyd, lh);
  }


  void CompoundBilinearFormIntegrator :: 
  CalcFlux (const FiniteElement & bfel,
	    const BaseMappedIntegrationRule & mir,
            BareSliceVector<Complex> elx, 
	    BareSliceMatrix<Complex> flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      static_cast<const CompoundFiniteElement&> (bfel);
    
    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();

    bfi->CalcFlux (fel[comp], mir, elx.Range (base, base+fel[comp].GetNDof()), 
		  flux, applyd, lh);
  }



  void CompoundBilinearFormIntegrator :: 
  ApplyBTrans (const FiniteElement & bfel,
	       const BaseMappedIntegrationPoint & mip,
               FlatVector<double> elx, 
	       FlatVector<double> ely,
	       LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      static_cast<const CompoundFiniteElement&> (bfel);

    ely = 0;
    FlatVector<double> hy = ely.Range (fel.GetRange(comp)); 
    bfi->ApplyBTrans (fel[comp], mip, elx, hy, lh);
  }
  
  void CompoundBilinearFormIntegrator :: 
  ApplyBTrans (const FiniteElement & bfel,
	       const BaseMappedIntegrationPoint & mip,
               FlatVector<Complex> elx, 
	       FlatVector<Complex> ely,
	       LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      static_cast<const CompoundFiniteElement&> (bfel);

    ely = 0;
    FlatVector<Complex> hy = ely.Range (fel.GetRange(comp)); 
    bfi->ApplyBTrans (fel[comp], mip, elx, hy, lh);
  }



  string CompoundBilinearFormIntegrator :: Name () const
  {
    return
      string ("CompoundIntegrator (") + bfi->Name() +
      string (")");
  }




  void LinearFormIntegrator ::
  CalcElementVector (const FiniteElement & fel,
                     const ElementTransformation & eltrans, 
                     FlatVector<double> elvec,
                     LocalHeap & lh) const
  { 
    cerr << "LinearFormIntegrator::CalcElementVector: base class called" << endl;
  }
  

  void LinearFormIntegrator ::
  CalcElementVector (const FiniteElement & fel,
                     const ElementTransformation & eltrans, 
                     FlatVector<Complex> elvec,
                     LocalHeap & lh) const
  {
    FlatVector<double> rvec(elvec.Size(), lh);
    CalcElementVector (fel, eltrans, rvec, lh);
    elvec = rvec;
  }




  BlockLinearFormIntegrator :: 
  BlockLinearFormIntegrator (shared_ptr<LinearFormIntegrator> alfi, int adim, int acomp)
    : lfi(alfi), dim(adim), comp(acomp)
  {
    ;
  }

  void BlockLinearFormIntegrator :: 
  CalcElementVector (const FiniteElement & bfel, 
		     const ElementTransformation & eltrans, 
		     FlatVector<double> elvec,
		     LocalHeap & lh) const
  {
    FlatVector<double> vec1(bfel.GetNDof(), lh);
    lfi->CalcElementVector (bfel, eltrans, vec1, lh);
    elvec = 0;
    if (comp == -1)
      for (int i = 0; i < vec1.Size(); i++)
        for (int k=0; k<dim; k++)
          elvec(dim*i+k) = vec1(i);
    else
      for (int i = 0; i < vec1.Size(); i++)
        elvec(dim*i+comp) = vec1(i);
  }  










  void CompoundLinearFormIntegrator ::  
  CalcElementVector (const FiniteElement & bfel, 
		     const ElementTransformation & eltrans, 
		     FlatVector<double> elvec,
		     LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    FlatVector<double> vec1(fel[comp].GetNDof(), lh);
    lfi->CalcElementVector (fel[comp], eltrans, vec1, lh);
    
    elvec = 0;

    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();

    for (int i = 0; i < vec1.Size(); i++)
      elvec(base+i) = vec1(i);
  }  



  void CompoundLinearFormIntegrator ::  
  CalcElementVector (const FiniteElement & bfel, 
		     const ElementTransformation & eltrans, 
		     FlatVector<Complex> elvec,
		     LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel =
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    FlatVector<Complex> vec1(fel[comp].GetNDof(), lh);
    lfi->CalcElementVector (fel[comp], eltrans, vec1, lh);
    
    elvec = 0;

    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();

    for (int i = 0; i < vec1.Size(); i++)
      elvec(base+i) = vec1(i);
  }  


  void CompoundLinearFormIntegrator ::  
  CalcElementVectorIndependent (const FiniteElement & gfel,
                                const BaseMappedIntegrationPoint & s_mip,
                                const BaseMappedIntegrationPoint & g_mip,
                                FlatVector<double> & elvec,
                                LocalHeap & lh,
                                const bool curveint) const
  {
    const CompoundFiniteElement & fel =
      dynamic_cast<const CompoundFiniteElement&> (gfel);

    int i;
    FlatVector<double> vec1;
    lfi->CalcElementVectorIndependent (fel[comp], s_mip, g_mip, vec1, lh, curveint);
    
    elvec.AssignMemory (fel.GetNDof(), lh);
    elvec = 0;

    int base = 0;
    for (i = 0; i < comp; i++)
      base += fel[i].GetNDof();

    for (i = 0; i < vec1.Size(); i++)
      elvec(base+i) = vec1(i);
  }

  void CompoundLinearFormIntegrator ::  
  CalcElementVectorIndependent (const FiniteElement & gfel,
                                const BaseMappedIntegrationPoint & s_mip,
                                const BaseMappedIntegrationPoint & g_mip,
                                FlatVector<Complex> & elvec,
                                LocalHeap & lh,
                                const bool curveint) const
  {
    const CompoundFiniteElement & fel =
      dynamic_cast<const CompoundFiniteElement&> (gfel);

    FlatVector<Complex> vec1;
    lfi->CalcElementVectorIndependent (fel[comp], s_mip, g_mip, vec1, lh, curveint);
    
    elvec.AssignMemory (fel.GetNDof(), lh);
    elvec = 0;

    int base = 0;
    for (int i = 0; i < comp; i++)
      base += fel[i].GetNDof();

    for (int i = 0; i < vec1.Size(); i++)
      elvec(base+i) = vec1(i);
  }



  void LinearFormIntegratorAnyDim ::  
  CheckElement (const FiniteElement & el) const
  {
    ;
  }

  void LinearFormIntegratorAnyDim ::  
  CalcElementVector (const FiniteElement & bfel, 
                     const ElementTransformation & eltrans, 
                     FlatVector<double> elvec,
                     LocalHeap & lh) const
  {
    int dim = eltrans.SpaceDim();
    if (lfi[dim])
      lfi[dim] -> CalcElementVector(bfel, eltrans, elvec, lh);
    else
      throw Exception (ToString("Integrator-Anydim not available for dimension ")+
                       ToString(dim));
  }

  void LinearFormIntegratorAnyDim ::  
  CalcElementVector (const FiniteElement & bfel, 
                     const ElementTransformation & eltrans, 
                     FlatVector<Complex> elvec,
                     LocalHeap & lh) const
  {
    int dim = eltrans.SpaceDim();
    if (lfi[dim])
      lfi[dim] -> CalcElementVector(bfel, eltrans, elvec, lh);
    else
      throw Exception (ToString("Integrator-Anydim not available for dimension ")+
                       ToString(dim));
  }






  template <typename T>
  Integrators::IntegratorInfo<T>::
  IntegratorInfo (const string & aname,
		   int aspacedim, int anumcoeffs,
		   shared_ptr<T> (*acreator)(const Array<shared_ptr<CoefficientFunction>>&))
    : name(aname), spacedim(aspacedim), numcoeffs(anumcoeffs), creator(acreator)
  {
    ;
  }
  
  
  Integrators :: Integrators ()
  {
    ;
  }

  Integrators :: ~Integrators()
  {
    for (int i = 0; i < bfis.Size(); i++)
      delete bfis[i];
    for (int i = 0; i < lfis.Size(); i++)
      delete lfis[i];
  }
  
  void Integrators :: 
  AddBFIntegrator (const string & aname, int aspacedim, int anumcoeffs,
                   shared_ptr<BilinearFormIntegrator>(*creator)(const Array<shared_ptr<CoefficientFunction>> &))
  {
    bfis.Append (new IntegratorInfo<BilinearFormIntegrator>(aname, aspacedim, anumcoeffs, creator));
  }

  const Integrators::IntegratorInfo<BilinearFormIntegrator> * 
  Integrators::GetBFI(const string & name, int spacedim) const
  {
    for (int i = 0; i < bfis.Size(); i++)
      if (name == bfis[i]->name && spacedim == bfis[i]->spacedim)
	return bfis[i];

    throw Exception (string ("GetBFI: Unknown integrator ") + name + "\n");
  }

  shared_ptr<BilinearFormIntegrator>
  Integrators::CreateBFI(const string & name, int dim, 
			 const Array<shared_ptr<CoefficientFunction>> & coeffs) const
  {
    if (dim == -1)
      {
        shared_ptr<BilinearFormIntegrator> abfi[4];
        for (int i = 0; i < bfis.Size(); i++)
          if (name == bfis[i]->name)
            {
              int spacedim = bfis[i]->spacedim;
              abfi[spacedim] = bfis[i] -> creator(coeffs);
            }
        return make_shared<BilinearFormIntegratorAnyDim> (abfi);
      }
    auto bfi = GetBFI(name, dim)->creator(coeffs);
    bfi -> SetName (name);
    return bfi;
  }

  shared_ptr<BilinearFormIntegrator>
  Integrators::CreateBFI(const string & name, int dim, 
			 shared_ptr<CoefficientFunction> coef) const
  {
    Array<shared_ptr<CoefficientFunction>> coeffs(1);
    coeffs[0] = coef;
    return CreateBFI (name, dim, coeffs);
  }

  shared_ptr<BilinearFormIntegrator>
  Integrators::CreateBFI(const string & name, int dim, 
			 const CoefficientFunction * coef) const
  {
    Array<shared_ptr<CoefficientFunction>> coeffs(1);
    coeffs[0] = shared_ptr<CoefficientFunction> 
      (const_cast<CoefficientFunction*>(coef), NOOP_Deleter);
    return CreateBFI (name, dim, coeffs);
  }


  void Integrators :: 
  AddLFIntegrator (const string & aname, int aspacedim, int anumcoeffs,
		   shared_ptr<LinearFormIntegrator>(*creator)(const Array<shared_ptr<CoefficientFunction>>&))
  {
    lfis.Append (new IntegratorInfo<LinearFormIntegrator>(aname, aspacedim, anumcoeffs, creator));
  }

  const Integrators::IntegratorInfo<LinearFormIntegrator> * 
  Integrators::GetLFI(const string & name, int spacedim) const
  {
    for (int i = 0; i < lfis.Size(); i++)
      if (name == lfis[i]->name && spacedim == lfis[i]->spacedim)
	return lfis[i];

    throw Exception (string ("GetLFI: Unknown integrator ") + name + "\n");
  }


  shared_ptr<LinearFormIntegrator>
  Integrators::CreateLFI(const string & name, int dim, 
			 const Array<shared_ptr<CoefficientFunction>> & coeffs) const
  {
    // LinearFormIntegrator * lfi =
    // dynamic_cast<LinearFormIntegrator*> (GetLFI(name, dim)->creator(coeffs));

    if (dim == -1)
      {
        shared_ptr<LinearFormIntegrator> alfi[4];
        for (int i = 0; i < lfis.Size(); i++)
          if (name == lfis[i]->name)
            {
              int spacedim = lfis[i]->spacedim;
              try
                {
                  alfi[spacedim] = lfis[i] -> creator(coeffs);
                }
              catch (Exception e)
                {
                  // creates lfi only if dimension of coef-function fits
                }
            }
        return make_shared<LinearFormIntegratorAnyDim> (alfi);
      }

    auto lfi = GetLFI(name, dim)->creator(coeffs);
    lfi -> SetName (name);
    return lfi;
  }

  shared_ptr<LinearFormIntegrator>
  Integrators::CreateLFI(const string & name, int dim, 
			 shared_ptr<CoefficientFunction> coef) const
  {
    Array<shared_ptr<CoefficientFunction>> coeffs(1);
    coeffs[0] = coef;
    return CreateLFI (name, dim, coeffs);
  }




  void Integrators :: Print (ostream & ost) const
  {
    int i;
    ost << endl << "Bilinear-form integrators:" << endl;
    ost <<         "--------------------------" << endl;
    ost << setw(20) << "Name" 
	<< setw(4) << "dim"
	<< setw(4) << "nco" << endl;
    for (i = 0; i < bfis.Size(); i++)
      ost << setw(20) << bfis[i]->name
	  << setw(4) << bfis[i]->spacedim
	  << setw(4) << bfis[i]->numcoeffs
	  << endl;

    ost << endl << "Linear-form integrators:" << endl;
    ost <<         "------------------------" << endl;
    ost << setw(20) << "Name" 
	<< setw(4) << "dim"
	<< setw(4) << "nco" << endl;
    for (i = 0; i < lfis.Size(); i++)
      ost << setw(20) << lfis[i]->name
	  << setw(4) << lfis[i]->spacedim
	  << setw(4) << lfis[i]->numcoeffs
	  << endl;
  }




 
  Integrators & GetIntegrators ()
  {
    static Integrators itgs;
    return itgs;
  }


}
  
