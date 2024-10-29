#ifndef FILE_GRIDFUNCTION
#define FILE_GRIDFUNCTION

/*********************************************************************/
/* File:   gridfunction.hpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Jul. 2000                                             */
/*********************************************************************/


#include <meshing/soldata.hpp>   // netgen visualization
#include "fespace.hpp"

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
				     shared_ptr<DifferentialOperator> atttrace_diffop = nullptr,                                     
                                     int acomp = 0);
  public:
    GridFunctionCoefficientFunction () = default;
    GridFunctionCoefficientFunction (shared_ptr<GridFunction> agf, int acomp = 0);
    GridFunctionCoefficientFunction (shared_ptr<GridFunction> agf, 
                                     shared_ptr<DifferentialOperator> adiffop,
                                     shared_ptr<DifferentialOperator> atrace_diffop = nullptr,
				     shared_ptr<DifferentialOperator> attrace_diffop = nullptr,
				     shared_ptr<DifferentialOperator> atttrace_diffop = nullptr,
                                     int acomp = 0);
    GridFunctionCoefficientFunction (shared_ptr<GridFunction> agf, 
                                     shared_ptr<BilinearFormIntegrator> abfi, int acomp = 0);
    GridFunctionCoefficientFunction (shared_ptr<GridFunction> agf, shared_ptr<ProxyFunction> proxy);
    
    virtual ~GridFunctionCoefficientFunction ();
    void DoArchive(Archive& ar) override;
    /// scalar valued or vector valued
    virtual bool IsComplex() const;
    virtual int Dimension() const;
    virtual Array<int> Dimensions() const;
    virtual bool DefinedOn (const ElementTransformation & trafo) override;
    void SelectComponent (int acomp) { comp = acomp; }
    const GridFunction & GetGridFunction() const { return *gf; }
      using CoefficientFunction::Evaluate;
      using CoefficientFunctionNoDerivative::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override;
    virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const override; 

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const override;
    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const override;
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
			   BareSliceMatrix<double> values) const override;
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
			   BareSliceMatrix<Complex> values) const override;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const override;

    virtual bool StoreUserData() const override { return true; }

    virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<AutoDiffDiff<1,NonZero>> nonzero) const override
    {
      nonzero = AutoDiffDiff<1,NonZero> (true);
    }

    // generation information for pickling:
    bool generated_from_deriv = false;
    string generated_from_operator;
    shared_ptr<GridFunction> GetGridFunctionPtr() const { return gf_shared_ptr; }
    const auto & GetDifferentialOperator (VorB vb) const { return diffop[vb]; }

    virtual shared_ptr<CoefficientFunction>
      Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override;

    shared_ptr<CoefficientFunction> Primary() const override;
    shared_ptr<GridFunctionCoefficientFunction> GetTrace() const;
  };




 

  /** 
      Grid-functions
  */
  class NGS_DLL_HEADER GridFunction 
    : public GridFunctionCoefficientFunction // , public NGS_Object
  {
  protected:
    string name;
    Flags flags;
    Flags flaglist;
    /// the finite element space
    shared_ptr<FESpace> fespace;
    /// should we do a prolongation from one multigrid-level to the next ?
    bool nested;
    /// should the GF be updated automatically on an update of the underliying space?
    bool autoupdate;
    /// should we visualize the gridfunction ?
    bool visual;
    /// how many functions
    int multidim;
    /// highest multigrid-level for which Update was called (memory allocation)
    int level_updated = -1;

    /// used for many right-hand-sides
    int cacheblocksize = 1;
    /// Component gridfunctions living somewhere - need to call update on them
    Array<weak_ptr<GridFunction>> compgfs;
    /// the actual data, array for multi-dim 
    Array<shared_ptr<BaseVector>> vec;
    /// component GridFunctions if fespace is a CompoundFESpace
    weak_ptr<GridFunctionCoefficientFunction> derivcf;
  public:
    /// 
    explicit GridFunction () = default;
    GridFunction (shared_ptr<FESpace> afespace, 
		  const string & name = "gfu", 
		  const Flags & flags = Flags());
    ///
    virtual ~GridFunction ();
    ///
    virtual void Update ();
    ///
    bool DoesAutoUpdate () const { return autoupdate; }
    void ConnectAutoUpdate();
    
    ///
    virtual void DoArchive (Archive & archive) override;
    ///  
    virtual BaseVector & GetVector (int comp = 0) { return *vec[comp]; }
    virtual const BaseVector & GetVector (int comp = 0) const  { return *vec[comp]; }
    ///  
    virtual shared_ptr<BaseVector> GetVectorPtr (int comp = 0) const  { return vec[comp]; }
    auto GetVectors() const { return MultiVector(vec); }
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
    shared_ptr<MeshAccess> GetMeshAccess() const { return fespace->GetMeshAccess(); }
    ///
    virtual string GetClassName () const
    {
      return "GridFunction";
    }

    void SetName(const string & aname) { name = aname; }
    string GetName () const { return name; }

    const Flags& GetFlags() const { return flags; }
    Flags& GetFlags() { return flags; }

    virtual void Interpolate (const CoefficientFunction & cf,
                              const Region * reg, int mdcomp, LocalHeap & lh);

    ///
    virtual void PrintReport (ostream & ost) const override;
    ///
    virtual Array<MemoryUsage> GetMemoryUsage () const;

    // void Visualize(const string & name);

    ///
    virtual void SetCacheBlockSize (const int size) {cacheblocksize = size;}
    ///
    virtual int GetCacheBlockSize (void) const { return cacheblocksize;}
    ///
    virtual bool IsUpdated () const; 

    
    int GetNComponents () const
    {
      auto compfes = dynamic_pointer_cast<CompoundFESpace>(fespace);
      if(compfes)
        return compfes->GetNSpaces();
      return 0;
    }
    shared_ptr<GridFunction> GetComponent (int compound_comp);

    [[deprecated("Use GridFunction.Deriv instead!")]]
    shared_ptr<GridFunctionCoefficientFunction> GetDeriv()
    {
      return Deriv();
    }

    shared_ptr<GridFunctionCoefficientFunction> Deriv()
    {
      auto res = derivcf.lock();
      if(res) return res;

      // derivcf is not set -> initialize it
      res =
        make_shared<GridFunctionCoefficientFunction> (dynamic_pointer_cast<GridFunction> (shared_from_this()),
                                                      GetFESpace()->GetFluxEvaluator(),
                                                      GetFESpace()->GetFluxEvaluator(BND),
                                                      GetFESpace()->GetFluxEvaluator(BBND));
      res -> generated_from_deriv = true;
      derivcf = res;
      return res;
    }

    shared_ptr<CoefficientFunction> Operator (shared_ptr<DifferentialOperator> diffop) const override;
    shared_ptr<CoefficientFunction> Operator (const string& name) const override;
    
    
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


    // multidim component, if -1 then all components are loaded/saved
    virtual void Load (istream & ist, int mdcomp = -1) = 0;
    virtual void Save (ostream & ost, int mdcomp = -1) const = 0;
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
    S_GridFunction () = default;
    S_GridFunction (shared_ptr<FESpace> afespace, 
		    const string & aname = "gfu", 
		    const Flags & flags = Flags());
      
    // parallel Load/Save by Martin Huber and Lothar Nannen 
    virtual void Load (istream & ist, int mdcomp);
    virtual void Save (ostream & ost, int mdcomp) const;

    virtual void Update ();

  private:
    template <int N, NODE_TYPE NT> void LoadNodeType (istream & ist, int mdcomp);
    template <int N, NODE_TYPE NT> void SaveNodeType (ostream & ost, int mdcomp) const;
  };



  template <class TV>
  class NGS_DLL_HEADER T_GridFunction : public S_GridFunction<typename mat_traits<TV>::TSCAL>
  {
  public:
    typedef typename mat_traits<TV>::TSCAL TSCAL;

    [[deprecated("Use T_GridFunction(shared_ptr<FESpace> ... instead!")]] 
    T_GridFunction (const FESpace & afespace, 
		    const string & aname = "gfu", 
		    const Flags & flags = Flags())
      : T_GridFunction(shared_ptr<FESpace> (const_cast<FESpace*>(&afespace),NOOP_Deleter), aname, flags)
      { ; } 

    T_GridFunction (shared_ptr<FESpace> afespace, 
		    const string & aname = "gfu", 
		    const Flags & flags = Flags())
      : S_GridFunction<TSCAL> (afespace, aname, flags)
      { ; }

    virtual ~T_GridFunction () { ; } 
  };




  extern NGS_DLL_HEADER 
  shared_ptr<GridFunction> CreateGridFunction (shared_ptr<FESpace> space,
                                               const string & name, const Flags & flags);

  [[deprecated("Use CreateGridFunction(shared_ptr<FESpace> ... instead!")]]   
  inline 
  shared_ptr<GridFunction> CreateGridFunction (const FESpace * space,
                                               const string & name, const Flags & flags)
  {
    return 
      CreateGridFunction (shared_ptr<FESpace> (const_cast<FESpace*>(space), NOOP_Deleter), 
                          name, flags);
  }



  class NGS_DLL_HEADER ComponentGridFunction : public GridFunction
  {
    shared_ptr<GridFunction> gf_parent;
    int comp;
  public:
    ComponentGridFunction (shared_ptr<GridFunction> agf_parent, int acomp);
    ~ComponentGridFunction () override;
    void Update () override;
    void Load(istream& ist, int mdcomp) override { throw Exception("Load not implemented for ComponentGF"); }
    void Save(ostream& ost, int mdcomp) const override { throw Exception("Save not implemented for ComponentGF"); }
    shared_ptr<GridFunction> GetParent() const { return gf_parent; }
    int GetComponent() const { return comp; }
  };
  




  
  template <class SCAL>
  class NGS_DLL_HEADER VisualizeGridFunction : public netgen::SolutionData
  {
    shared_ptr<MeshAccess> ma;
    shared_ptr<GridFunction> gf;
    Array<shared_ptr<BilinearFormIntegrator>> bfi2d;
    Array<shared_ptr<BilinearFormIntegrator>> bfi3d;
    bool applyd;

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

    virtual bool GetMultiSurfValue (size_t selnr, size_t facetnr, size_t npts,
                                    const SIMD<double> * xref,
                                    const SIMD<double> * x,
                                    const SIMD<double> * dxdxref,
                                    SIMD<double> * values);

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

    virtual bool GetMultiSurfValue (size_t selnr, size_t facetnr, size_t npts,
                                    const SIMD<double> * xref,
                                    const SIMD<double> * x,
                                    const SIMD<double> * dxdxref,
                                    SIMD<double> * values);

    virtual bool GetSegmentValue (int segnr, double xref, double * values);

    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages, int component = -1);
    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages_times_volumes, Array<double> & volumes, int component = -1);


    virtual int GetNumMultiDimComponents ();
    virtual void SetMultiDimComponent (int mc);
  };



}



#endif
