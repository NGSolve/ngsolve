#ifndef FILE_GRIDFUNCTION
#define FILE_GRIDFUNCTION

/*********************************************************************/
/* File:   gridfunction.hpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Jul. 2000                                             */
/*********************************************************************/

namespace ngcomp
{


  class GridFunction;

  class NGS_DLL_HEADER GridFunctionCoefficientFunction : public CoefficientFunctionNoDerivative
  {
  protected:
    shared_ptr<GridFunction> gf_shared_ptr;
    GridFunction* gf;
    shared_ptr<FESpace> fes;
    shared_ptr<DifferentialOperator> diffop[4];
    int comp;
    GridFunctionCoefficientFunction (shared_ptr<DifferentialOperator> adiffop,
                                     shared_ptr<DifferentialOperator> atrace_diffop = nullptr,
				     shared_ptr<DifferentialOperator> attrace_diffop = nullptr,
                                     int acomp = 0);
  public:
    GridFunctionCoefficientFunction (shared_ptr<GridFunction> agf, int acomp = 0);
    GridFunctionCoefficientFunction (shared_ptr<GridFunction> agf, 
                                     shared_ptr<DifferentialOperator> adiffop,
                                     shared_ptr<DifferentialOperator> atrace_diffop = nullptr,
				     shared_ptr<DifferentialOperator> attrace_diffop = nullptr,
                                     int acomp = 0);
    GridFunctionCoefficientFunction (shared_ptr<GridFunction> agf, 
                                     shared_ptr<BilinearFormIntegrator> abfi, int acomp = 0);
    
    virtual ~GridFunctionCoefficientFunction ();
    /// scalar valued or vector valued
    virtual bool IsComplex() const;
    virtual int Dimension() const;
    virtual Array<int> Dimensions() const;
    virtual bool DefinedOn (const ElementTransformation & trafo) override;
    void SelectComponent (int acomp) { comp = acomp; }
    const GridFunction & GetGridFunction() const { return *gf; }
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override;
    virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const override; 

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const override;
    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const override;
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
			   BareSliceMatrix<double> values) const override;
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
			   FlatMatrix<Complex> values) const override;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const override;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const override
    { Evaluate (ir, values); }
    /*
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    { Evaluate (ir, result); deriv = 0.0; }
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input,
                                FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    { Evaluate (ir, result); deriv = 0.0; }
    */
    virtual bool StoreUserData() const override { return true; }

    virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                                 FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
    {
      nonzero = true;
      nonzero_deriv = false;
      nonzero_dderiv = false;
    }

    // generation information for pickling:
    bool generated_from_deriv = false;
    string generated_from_operator;
  };




 

  /** 
      Grid-functions
  */
  class NGS_DLL_HEADER GridFunction 
    : public NGS_Object, public GridFunctionCoefficientFunction
  {
  protected:
    /// the finite element space
    shared_ptr<FESpace> fespace;
    /// should we do a prolongation from one multigrid-level to the next ?
    bool nested;
    /// should we visualize the gridfunction ?
    bool visual;
    /// how many functions
    int multidim;
    /// highest multigrid-level for which Update was called (memory allocation)
    int level_updated = -1;

    /// used for many right-hand-sides
    int cacheblocksize = 1;
    /// the actual data, array for multi-dim 
    Array<shared_ptr<BaseVector>> vec;
    /// component GridFunctions if fespace is a CompoundFESpace
    Array<shared_ptr<GridFunction>> compgfs;
  public:
    /// 
    GridFunction (shared_ptr<FESpace> afespace, 
		  const string & name = "gfu", 
		  const Flags & flags = Flags());
    ///
    virtual ~GridFunction ();
    ///
    virtual void Update ();
    ///
    virtual void DoArchive (Archive & archive);
    ///  
    virtual BaseVector & GetVector (int comp = 0) { return *vec[comp]; }
    virtual const BaseVector & GetVector (int comp = 0) const  { return *vec[comp]; }
    ///  
    virtual shared_ptr<BaseVector> GetVectorPtr (int comp = 0) const  { return vec[comp]; }
    ///
    void SetNested (int anested = 1) { nested = anested; }
    ///
    bool GetVisual () const { return visual; }
    ///
    void SetVisual (bool avisual = 1) { visual = avisual; }

    int GetMultiDim () const { return multidim; }

    /// increase multidim and copy vec to new component
    void AddMultiDimComponent (BaseVector & vec);
  
    int GetLevelUpdated() const { return level_updated; }
    ///

    // const FESpace & GetFESpace() const { return *fespace; }
    ///
    shared_ptr<FESpace> GetFESpace() const { return fespace; }

    ///
    virtual string GetClassName () const
    {
      return "GridFunction";
    }

    virtual string GetName () const
    {
      return name;
    }
    

    ///
    virtual void PrintReport (ostream & ost) const;
    ///
    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;

    ///
    // void Visualize(const string & name);

    ///
    virtual void SetCacheBlockSize (const int size) {cacheblocksize = size;}
    ///
    virtual int GetCacheBlockSize (void) const { return cacheblocksize;}
    ///
    virtual bool IsUpdated () const; 

    
    int GetNComponents () const { return compgfs.Size(); }
    shared_ptr<GridFunction> GetComponent (int compound_comp) const;


    ///
    virtual void GetElementVector (FlatArray<int> dnums,
				   FlatVector<double> elvec) const
    { vec[0] -> GetIndirect (dnums, elvec); }
      
    ///
    virtual void SetElementVector (FlatArray<int> dnums,
				   FlatVector<double> elvec) 
    { vec[0] -> SetIndirect (dnums, elvec); }


    ///
    virtual void GetElementVector (int comp,
				   FlatArray<int> dnums,
				   FlatVector<double> elvec) const 
    { vec[comp] -> GetIndirect (dnums, elvec); }


    ///
    virtual void SetElementVector (int comp,
				   FlatArray<int> dnums,
				   FlatVector<double> elvec) 
    { vec[comp] -> SetIndirect (dnums, elvec); }

    ///
    virtual void GetElementVector (FlatArray<int> dnums,
				   FlatVector<Complex> elvec) const
    { vec[0] -> GetIndirect (dnums, elvec); }
      
    ///
    virtual void SetElementVector (FlatArray<int> dnums,
				   FlatVector<Complex> elvec) 
    { vec[0] -> SetIndirect (dnums, elvec); }


    ///
    virtual void GetElementVector (int comp,
				   FlatArray<int> dnums,
				   FlatVector<Complex> elvec) const 
    { vec[comp] -> GetIndirect (dnums, elvec); }


    ///
    virtual void SetElementVector (int comp,
				   FlatArray<int> dnums,
				   FlatVector<Complex> elvec) 
    { vec[comp] -> SetIndirect (dnums, elvec); }



    virtual void Load (istream & ist) = 0;
    virtual void Save (ostream & ost) const = 0;
  };


  inline ostream & operator<< (ostream & ost, const GridFunction & gf)
  {
    gf.PrintReport (ost);
    return ost;
  }

  
  extern NGS_DLL_HEADER void Visualize(shared_ptr<GridFunction> gf, const string & name);
  
  template <class SCAL>
  class NGS_DLL_HEADER S_GridFunction : public GridFunction
  {
  public:
    /*
    S_GridFunction (const FESpace & afespace, const string & aname, const Flags & flags)
      : GridFunction (afespace, aname, flags) { ; }
    S_GridFunction (shared_ptr<FESpace> afespace, const string & aname, const Flags & flags)
      : GridFunction (afespace, aname, flags) { ; }
    */
    using GridFunction::GridFunction;


    // parallel Load/Save by Martin Huber and Lothar Nannen 
    virtual void Load (istream & ist);
    virtual void Save (ostream & ost) const;

  private:
    template <int N, NODE_TYPE NT> void LoadNodeType (istream & ist);

    template <int N, NODE_TYPE NT> void SaveNodeType (ostream & ost) const;
  };



  template <class TV>
  class NGS_DLL_HEADER T_GridFunction : public S_GridFunction<typename mat_traits<TV>::TSCAL>
  {
    using S_GridFunction<typename mat_traits<TV>::TSCAL>::vec;
    using S_GridFunction<typename mat_traits<TV>::TSCAL>::compgfs;

  public:
    typedef typename mat_traits<TV>::TSCAL TSCAL;
    enum { VDIM = mat_traits<TV>::HEIGHT };

    T_GridFunction (const FESpace & afespace, 
		    const string & aname = "gfu", 
		    const Flags & flags = Flags());
    T_GridFunction (shared_ptr<FESpace> afespace, 
		    const string & aname = "gfu", 
		    const Flags & flags = Flags());

    virtual ~T_GridFunction ();

    virtual void Update ();
  };




  extern NGS_DLL_HEADER 
  shared_ptr<GridFunction> CreateGridFunction (shared_ptr<FESpace> space,
                                               const string & name, const Flags & flags);

  /// compatibility with old codes
  inline 
  shared_ptr<GridFunction> CreateGridFunction (const FESpace * space,
                                               const string & name, const Flags & flags)
  {
    return 
      CreateGridFunction (shared_ptr<FESpace> (const_cast<FESpace*>(space), NOOP_Deleter), 
                          name, flags);
  }



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
  





  template <class SCAL>
  class NGS_DLL_HEADER VisualizeGridFunction : public netgen::SolutionData
  {
    shared_ptr<MeshAccess> ma;
    shared_ptr<S_GridFunction<SCAL>> gf;
    Array<shared_ptr<BilinearFormIntegrator>> bfi2d;
    Array<shared_ptr<BilinearFormIntegrator>> bfi3d;
    bool applyd;
    //
    // int cache_elnr;
    // bool cache_bound;
    // LocalHeap lh;
    // ElementTransformation eltrans;
    // const FiniteElement * fel;
    // Array<int> dnums;
    // FlatVector<SCAL> elu;

  public:
    VisualizeGridFunction (shared_ptr<MeshAccess> ama,
			   shared_ptr<GridFunction> agf,
			   const shared_ptr<BilinearFormIntegrator> abfi2d,
			   const shared_ptr<BilinearFormIntegrator> abfi3d,
			   bool aapplyd);
    VisualizeGridFunction (shared_ptr<MeshAccess> ama,
			   shared_ptr<GridFunction> agf,
			   const Array<shared_ptr<BilinearFormIntegrator>> & abfi2d,
			   const Array<shared_ptr<BilinearFormIntegrator>> & abfi3d,
			   bool aapplyd);

    virtual ~VisualizeGridFunction ();
  
    virtual bool GetValue (int elnr, 
			   double lam1, double lam2, double lam3,
			   double * values) ;

    virtual bool GetValue (int elnr, 
			   const double xref[], const double x[], const double dxdxref[],
			   double * values) ;

    virtual bool GetMultiValue (int elnr, int facetnr, int npts,
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

#ifdef __AVX__
    virtual bool GetMultiSurfValue (size_t selnr, size_t facetnr, size_t npts,
                                    const tAVXd * xref, 
                                    const tAVXd * x, 
                                    const tAVXd * dxdxref, 
                                    tAVXd * values);
#endif

    virtual bool GetSegmentValue (int segnr, double xref, double * values);

    virtual int GetNumMultiDimComponents ()
    {
      return gf->GetMultiDim();
    }


    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages, int component = -1);
    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages_times_volumes, Array<double> & volumes, int component = -1);
    

  };













  class NGS_DLL_HEADER VisualizeCoefficientFunction : public netgen::SolutionData
  {
    shared_ptr<MeshAccess> ma;
    shared_ptr<CoefficientFunction> cf;
    // LocalHeap lh;
  public:
    VisualizeCoefficientFunction (shared_ptr<MeshAccess> ama,
				  shared_ptr<CoefficientFunction> acf);
    
    virtual ~VisualizeCoefficientFunction ();
  
    virtual bool GetValue (int elnr, 
			   double lam1, double lam2, double lam3,
			   double * values) ;

    virtual bool GetValue (int elnr, 
			   const double xref[], const double x[], const double dxdxref[],
			   double * values) ;

    virtual bool GetMultiValue (int elnr, int facetnr, int npts,
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

#ifdef __AVX__
    virtual bool GetMultiSurfValue (size_t selnr, size_t facetnr, size_t npts,
                                    const tAVXd * xref, 
                                    const tAVXd * x, 
                                    const tAVXd * dxdxref, 
                                    tAVXd * values);
#endif

    virtual bool GetSegmentValue (int segnr, double xref, double * values);

    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages, int component = -1);
    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages_times_volumes, Array<double> & volumes, int component = -1);


    virtual int GetNumMultiDimComponents ();
    virtual void SetMultiDimComponent (int mc);
  };



}



#endif
