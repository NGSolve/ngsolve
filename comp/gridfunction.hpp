#ifndef FILE_GRIDFUNCTION
#define FILE_GRIDFUNCTION

/*********************************************************************/
/* File:   gridfunction.hpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Jul. 2000                                             */
/*********************************************************************/

namespace ngcomp
{

  /** 
      Grid-functions
  */
  class NGS_DLL_HEADER GridFunction : public NGS_Object
  {
  protected:
    /// the finite element space
    const FESpace & fespace;
    /// should we do a prolongation from one multigrid-level to the next ?
    bool nested;
    /// should we visualize the gridfunction ?
    bool visual;
    /// how many functions
    int multidim;
    /// highest multigrid-level for which Update was called (memory allocation)
    int level_updated;

    /// used for many right-hand-sides
    int cacheblocksize;
    /// the actual data, array for multi-dim 
    Array<BaseVector*> vec;
    /// component GridFunctions if fespace is a CompoundFESpace
    Array<GridFunction*> compgfs;

  public:
    /// 
    GridFunction (const FESpace & afespace, const string & name, const Flags & flags);
    ///
    virtual ~GridFunction ();
    ///
    virtual void Update ();
    ///  
    virtual BaseVector & GetVector (int comp = 0) { return *(vec[comp]); }
    ///  
    virtual const BaseVector & GetVector (int comp = 0) const  { return *(vec[comp]); }


    /*
    operator BaseVector& () { return GetVector(); }
    template <typename T> 
    GridFunction & operator= (const VVecExpr<T> & v) { GetVector() = v; return *this; }
    GridFunction & operator= (const BaseVector & v) { GetVector() = v; return *this; }
    template <typename T> 
    GridFunction & operator+= (const VVecExpr<T> & v) { GetVector() += v; return *this; }
    GridFunction & operator+= (const BaseVector & v) { GetVector() += v; return *this; }
    template <typename T> 
    GridFunction & operator-= (const VVecExpr<T> & v) { GetVector() -= v; return *this; }
    GridFunction & operator-= (const BaseVector & v) { GetVector() -= v; return *this; }
    */

    ///
    void SetNested (int anested = 1) { nested = anested; }
    ///
    void SetVisual (bool avisual = 1) { visual = avisual; }

    int GetMultiDim () const { return multidim; }
  
    int GetLevelUpdated() const { return level_updated; }
    ///
    const FESpace & GetFESpace() const
    { return fespace; }

    ///
    virtual string GetClassName () const
    {
      return "GridFunction";
    }

    ///
    virtual void PrintReport (ostream & ost);
    ///
    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;

    ///
    void Visualize(const string & name);

    ///
    virtual void SetCacheBlockSize (const int size) {cacheblocksize = size;}
    ///
    virtual int GetCacheBlockSize (void) const { return cacheblocksize;}
    ///
    virtual bool IsUpdated () const; 

    
    int GetNComponents () const { return compgfs.Size(); }
    GridFunction * GetComponent (int compound_comp) const;


    ///
    virtual void GetElementVector (const FlatArray<int> & dnums,
				   FlatVector<double> & elvec) const
    { vec[0] -> GetIndirect (dnums, elvec); }
      
    ///
    virtual void SetElementVector (const FlatArray<int> & dnums,
				   const FlatVector<double> & elvec) 
    { vec[0] -> SetIndirect (dnums, elvec); }


    ///
    virtual void GetElementVector (int comp,
				   const FlatArray<int> & dnums,
				   FlatVector<double> & elvec) const 
    { vec[comp] -> GetIndirect (dnums, elvec); }


    ///
    virtual void SetElementVector (int comp,
				   const FlatArray<int> & dnums,
				   const FlatVector<double> & elvec) 
    { vec[comp] -> SetIndirect (dnums, elvec); }

    ///
    virtual void GetElementVector (const FlatArray<int> & dnums,
				   FlatVector<Complex> & elvec) const
    { vec[0] -> GetIndirect (dnums, elvec); }
      
    ///
    virtual void SetElementVector (const FlatArray<int> & dnums,
				   const FlatVector<Complex> & elvec) 
    { vec[0] -> SetIndirect (dnums, elvec); }


    ///
    virtual void GetElementVector (int comp,
				   const FlatArray<int> & dnums,
				   FlatVector<Complex> & elvec) const 
    { vec[comp] -> GetIndirect (dnums, elvec); }


    ///
    virtual void SetElementVector (int comp,
				   const FlatArray<int> & dnums,
				   const FlatVector<Complex> & elvec) 
    { vec[comp] -> SetIndirect (dnums, elvec); }


  };



  template <class SCAL>
  class NGS_DLL_HEADER S_GridFunction : public GridFunction
  {
  public:
    S_GridFunction (const FESpace & afespace, const string & aname, const Flags & flags)
      : GridFunction (afespace, aname, flags) { ; }
  

    // virtual GridFunction * GetComponent (int compound_comp) const;
  };



  template <class TV>
  class NGS_DLL_HEADER T_GridFunction : public S_GridFunction<typename mat_traits<TV>::TSCAL>
  {
    using S_GridFunction<typename mat_traits<TV>::TSCAL>::vec;
    using S_GridFunction<typename mat_traits<TV>::TSCAL>::compgfs;

  public:
    typedef typename mat_traits<TV>::TSCAL TSCAL;
    enum { VDIM = mat_traits<TV>::HEIGHT };

    T_GridFunction (const FESpace & afespace, const string & aname, const Flags & flags);
    virtual ~T_GridFunction ();

    virtual void Update ();
  };


  extern NGS_DLL_HEADER GridFunction * CreateGridFunction (const FESpace * space,
					    const string & name, const Flags & flags);




  template <class SCAL>
  class NGS_DLL_HEADER S_ComponentGridFunction : public S_GridFunction<SCAL>
  {
    const S_GridFunction<SCAL> & gf_parent;
    int comp;
  public:
    S_ComponentGridFunction (const S_GridFunction<SCAL> & agf_parent, int acomp);
    virtual ~S_ComponentGridFunction ();
    virtual void Update ();
  };
  



  class NGS_DLL_HEADER GridFunctionCoefficientFunction : public CoefficientFunction
  {
  protected:
    GridFunction & gf;
    const DifferentialOperator * diffop;
    int comp;
  public:
    GridFunctionCoefficientFunction (GridFunction & agf, int acomp = 0);
    GridFunctionCoefficientFunction (GridFunction & agf, DifferentialOperator * adiffop, int acomp = 0);
    
    virtual ~GridFunctionCoefficientFunction ();
    /// scalar valued or vector valued
    virtual bool IsComplex() const;
    virtual int Dimension() const;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const;
    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const;
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
			   FlatMatrix<double> values) const;
  };



  template <class SCAL>
  class NGS_DLL_HEADER VisualizeGridFunction : public netgen::SolutionData
  {
    const MeshAccess & ma;
    const S_GridFunction<SCAL> * gf;
    Array<const BilinearFormIntegrator *> bfi2d;
    Array<const BilinearFormIntegrator *> bfi3d;
    bool applyd;
    //
    int cache_elnr;
    bool cache_bound;
    LocalHeap lh;
    // ElementTransformation eltrans;
    const FiniteElement * fel;
    Array<int> dnums;
    FlatVector<SCAL> elu;

  public:
    VisualizeGridFunction (const MeshAccess & ama,
			   const GridFunction * agf,
			   const BilinearFormIntegrator * abfi2d,
			   const BilinearFormIntegrator * abfi3d,
			   bool aapplyd);
    VisualizeGridFunction (const MeshAccess & ama,
			   const GridFunction * agf,
			   const Array<BilinearFormIntegrator *> & abfi2d,
			   const Array<BilinearFormIntegrator *> & abfi3d,
			   bool aapplyd);

    virtual ~VisualizeGridFunction ();
  
    virtual bool GetValue (int elnr, 
			   double lam1, double lam2, double lam3,
			   double * values) ;

    virtual bool GetValue (int elnr, 
			   const double xref[], const double x[], const double dxdxref[],
			   double * values) ;

    virtual bool GetMultiValue (int elnr, int npts,
				const double * xref, int sxref,
				const double * x, int sx,
				const double * dxdxref, int sdxdxref,
				double * values, int svalues);

    virtual bool GetSurfValue (int elnr, int facetnr,
			       double lam1, double lam2, 
			       double * values) ;

    virtual bool GetSurfValue (int selnr, int facetnr, 
			       const double xref[], const double x[], const double dxdxref[],
			       double * values);

    virtual bool GetMultiSurfValue (int selnr, int facetnr, int npts,
                                    const double * xref, int sxref,
                                    const double * x, int sx,
                                    const double * dxdxref, int sdxdxref,
                                    double * values, int svalues);


    virtual int GetNumMultiDimComponents ()
    {
      return gf->GetMultiDim();
    }


    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages, int component = -1);
    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages_times_volumes, Array<double> & volumes, int component = -1);
    

  };













  class NGS_DLL_HEADER VisualizeCoefficientFunction : public netgen::SolutionData
  {
    const MeshAccess & ma;
    const CoefficientFunction * cf;
    LocalHeap lh;
  public:
    VisualizeCoefficientFunction (const MeshAccess & ama,
				  const CoefficientFunction * acf);
    
    virtual ~VisualizeCoefficientFunction ();
  
    virtual bool GetValue (int elnr, 
			   double lam1, double lam2, double lam3,
			   double * values) ;

    virtual bool GetValue (int elnr, 
			   const double xref[], const double x[], const double dxdxref[],
			   double * values) ;

    virtual bool GetMultiValue (int elnr, int npts,
				const double * xref, int sxref,
				const double * x, int sx,
				const double * dxdxref, int sdxdxref,
				double * values, int svalues);

    virtual bool GetSurfValue (int elnr, int facetnr,
			       double lam1, double lam2, 
			       double * values) ;

    virtual bool GetSurfValue (int selnr, int facetnr, 
			       const double xref[], const double x[], const double dxdxref[],
			       double * values);

    virtual bool GetMultiSurfValue (int selnr, int facetnr, int npts,
                                    const double * xref, int sxref,
                                    const double * x, int sx,
                                    const double * dxdxref, int sdxdxref,
                                    double * values, int svalues);


    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages, int component = -1);
    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages_times_volumes, Array<double> & volumes, int component = -1);
  };



}




#endif
