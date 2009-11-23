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
  class GridFunction : public NGS_Object
  {
  protected:
    const FESpace & fespace;
    bool nested;
    bool visual;
    int multidim;
    int level_updated;

    // netgen::SolutionData * vis;
    int cacheblocksize;


    Array<BaseVector*> vec;

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

    ///
    void SetNested (int anested = 1) { nested = anested; }
    ///
    void SetVisual (bool avisual = 1) { visual = avisual; }

    int GetMultiDim () const { return multidim; }
  
    int GetLevelUpdated() const { return level_updated; }
    ///
    const FESpace & GetFESpace() const
    { return fespace; }

  
    // NgMutex & Mutex () { return mutex; }

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


    virtual GridFunction * GetComponent (int compound_comp) const = 0;
  };



  template <class SCAL>
  class S_GridFunction : public GridFunction
  {
  public:
    S_GridFunction (const FESpace & afespace, const string & aname, const Flags & flags)
      : GridFunction (afespace, aname, flags) { ; }
  
    ///
    virtual void GetElementVector (const Array<int> & dnums,
				   FlatVector<SCAL> & elvec) const
    { 
      vec[0] -> GetIndirect (dnums, elvec); 
    }
      
    ///
    virtual void SetElementVector (const Array<int> & dnums,
				   const FlatVector<SCAL> & elvec) 
    {
      vec[0] -> SetIndirect (dnums, elvec);
    }


    ///
    virtual void GetElementVector (int comp,
				   const Array<int> & dnums,
				   FlatVector<SCAL> & elvec) const 
    {
      vec[comp] -> GetIndirect (dnums, elvec);
    }


    virtual void SetElementVector (int comp,
				   const Array<int> & dnums,
				   const FlatVector<SCAL> & elvec) 
    {
      vec[comp] -> SetIndirect (dnums, elvec);
    }

    virtual GridFunction * GetComponent (int compound_comp) const;
  };



  template <class TV>
  class T_GridFunction : public S_GridFunction<typename mat_traits<TV>::TSCAL>
  {
    using S_GridFunction<typename mat_traits<TV>::TSCAL>::vec;

  public:
    typedef typename mat_traits<TV>::TSCAL TSCAL;
    enum { VDIM = mat_traits<TV>::HEIGHT };

    T_GridFunction (const FESpace & afespace, const string & aname, const Flags & flags);
    // virtual ~T_GridFunction ();

    virtual void Update ();

    // virtual bool IsUpdated () const;

    virtual T_BaseVector<TV> & GetVector (int comp = 0)
    { return dynamic_cast<T_BaseVector<TV>&>(*(vec[comp])); }
    // virtual const T_BaseVector<TV> & GetVector (int comp = 0) const  { return *vec[comp]; }
    virtual const T_BaseVector<TV> & GetVector (int comp = 0) const
    { return dynamic_cast<T_BaseVector<TV>&>(*(vec[comp])); }


    friend class GridFunctionCoefficientFunction;
  };


  extern GridFunction * CreateGridFunction (const FESpace * space,
					    const string & name, const Flags & flags);




  template <class SCAL>
  class S_ComponentGridFunction : public S_GridFunction<SCAL>
  {
    const S_GridFunction<SCAL> & gf;
    int comp;
  public:
    S_ComponentGridFunction (const S_GridFunction<SCAL> & agf, int acomp);
    virtual void Update ();
  };
  



  class GridFunctionCoefficientFunction : public CoefficientFunction
  {
  protected:
    ///
    S_GridFunction<double> & gf;
    /*
    mutable LocalHeap lh;
    mutable FlatVector<double> elu;
    mutable Array<int> dnums;
    mutable int cache_elnr;
    */
    int comp;
  public:
    ///
    GridFunctionCoefficientFunction (GridFunction & agf)
      : gf(dynamic_cast<S_GridFunction<double>&> (agf)),comp(0) // ,lh(1000000)
    { 
      // cache_elnr = -1;
    }			
			
    GridFunctionCoefficientFunction (GridFunction & agf, const int acomp)
      : gf(dynamic_cast<S_GridFunction<double>&> (agf)),comp(acomp) // ,lh(1000000)
    { 
      // cache_elnr = -1;
    }
		
    ///
    virtual ~GridFunctionCoefficientFunction () {}

    virtual double Evaluate (const BaseSpecificIntegrationPoint & ip) const;

    virtual void Evaluate(const BaseSpecificIntegrationPoint & ip,
			  FlatVector<> result) const;
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;

  };



  template <class SCAL>
  class VisualizeGridFunction : public netgen::SolutionData
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
    ElementTransformation eltrans;
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

    virtual bool GetSurfValue (int elnr,
			       double lam1, double lam2, 
			       double * values) ;

    virtual bool GetSurfValue (int selnr,
			       const double xref[], const double x[], const double dxdxref[],
			       double * values);

    virtual bool GetMultiSurfValue (int selnr, int npts,
                                    const double * xref, int sxref,
                                    const double * x, int sx,
                                    const double * dxdxref, int sdxdxref,
                                    double * values, int svalues);


    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages, int component = -1);
    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages_times_volumes, Array<double> & volumes, int component = -1);
    

  };

}




#endif
