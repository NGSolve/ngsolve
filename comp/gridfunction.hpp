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
    NgMutex mutex;
    netgen::SolutionData * vis;

    int cacheblocksize;
  public:
    ///
    GridFunction (const FESpace & afespace, const string & name, const Flags & flags);
    ///
    virtual ~GridFunction () { ; }
  
    ///
    virtual void Update () = 0;
    ///  
    virtual BaseVector & GetVector (int comp = 0) = 0;
    ///  
    virtual const BaseVector & GetVector (int comp = 0) const = 0;

    ///
    void SetNested (int anested = 1)
    { nested = anested; }

    ///
    void SetVisual (bool avisual = 1)
    { visual = avisual; }

    int GetMultiDim () const { return multidim; }
  
    int GetLevelUpdated() const { return level_updated; }
    ///
    const FESpace & GetFESpace() const
    { return fespace; }

  
    NgMutex & Mutex () { return mutex; }

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
    ///
    void Visualize(const string & name);

    ///
    virtual void SetCacheBlockSize (const int size)
    {
      cacheblocksize = size;
    }

    virtual int GetCacheBlockSize (void) const
    {
      return cacheblocksize;
    }

    virtual bool IsUpdated (void) const
    {
      return false;
    }
  
  };



  template <class SCAL>
  class S_GridFunction : public GridFunction
  {
  public:
    S_GridFunction (const FESpace & afespace, const string & aname, const Flags & flags)
      : GridFunction (afespace, aname, flags) { ; }
  
    ///
    virtual void GetElementVector (const Array<int> & dnums,
				   FlatVector<SCAL> & elvec) const = 0;
    ///
    virtual void SetElementVector (const Array<int> & dnums,
				   const FlatVector<SCAL> & elvec) = 0;


    ///
    virtual void GetElementVector (int comp,
				   const Array<int> & dnums,
				   FlatVector<SCAL> & elvec) const = 0;
    ///
    virtual void SetElementVector (int comp,
				   const Array<int> & dnums,
				   const FlatVector<SCAL> & elvec) = 0;
  };



  template <class TV>
  class T_GridFunction : public S_GridFunction<typename mat_traits<TV>::TSCAL>
  {
  protected:
    Array<VVector<TV>*> vec;
  
  public:
    typedef typename mat_traits<TV>::TSCAL TSCAL;
    enum { VDIM = mat_traits<TV>::HEIGHT };

    T_GridFunction (const FESpace & afespace, const string & aname, const Flags & flags);
    virtual ~T_GridFunction ();

    virtual void Update ();

    virtual bool IsUpdated (void) const;

    virtual BaseVector & GetVector (int comp = 0);
    virtual const BaseVector & GetVector (int comp = 0) const;

    ///
    virtual void GetElementVector (const Array<int> & dnums,
				   FlatVector<TSCAL> & elvec) const;

    ///
    virtual void SetElementVector (const Array<int> & dnums,
				   const FlatVector<TSCAL> & elvec);


    ///
    virtual void GetElementVector (int comp,
				   const Array<int> & dnums,
				   FlatVector<TSCAL> & elvec) const;

    ///
    virtual void SetElementVector (int comp,
				   const Array<int> & dnums,
				   const FlatVector<TSCAL> & elvec);


    /*
    // only implemented double element-const gf
    template <int S, int R>
    double Evaluate (const SpecificIntegrationPoint<S,R> & ip, int comp = 1)
    {
    
    const FESpace & fes = GetFESpace();
    
    const ElementTransformation & eltrans = ip.GetElementTransformation();
    
    int elnr = eltrans.GetElementNr();
    
    const FiniteElement & el = fes.GetFE(elnr);
    int nd = el.GetNDof();
    
    Array<int> dnums (nd); 
    fes.GetDofNrs (elnr, dnums);
    
    Vector<double> elemu (nd);
    
    GetVector().GetIndirect (dnums, elemu);
    fes.TransformVec (elnr, false, elemu, TRANSFORM_SOL);
    
    return  InnerProduct (el.GetShape(ip), elemu);

    return 0;
    }
    */
    friend class GridFunctionCoefficientFunction;
  };


  extern GridFunction * CreateGridFunction (const FESpace * space,
					    const string & name, const Flags & flags);





  class GridFunctionCoefficientFunction : public CoefficientFunction
  {
  protected:
    ///
    S_GridFunction<double> & gf;
    LocalHeap lh;
    FlatVector<double> elu;
    Array<int> dnums;
    int cache_elnr;
    int comp;

  public:
    ///
    GridFunctionCoefficientFunction (GridFunction & agf)
      : gf(dynamic_cast<S_GridFunction<double>&> (agf)),lh(1000000),comp(0)
    { 
      cache_elnr = -1;
      //cout << "Created GridFunctionCoefficientFunction with gf = " << &gf << endl;
    }			
			
    GridFunctionCoefficientFunction (GridFunction & agf, const int acomp)
      : gf(dynamic_cast<S_GridFunction<double>&> (agf)),lh(1000000),comp(acomp)
    { 
      cache_elnr = -1;
      //cout << "Created GridFunctionCoefficientFunction with gf = " << &gf << endl;
    }
		
    ///
    virtual ~GridFunctionCoefficientFunction () {}
    ///

    //   template <int S, int R>
    //   double Evaluate (const SpecificIntegrationPoint<S,R> & ip)
    //   {
    //     int elnr = ip.GetTransformation().GetElementNr();
    
    //     if (elnr < 0 || elnr >= gf.vec[0]->Size())
    //       {
    // 	ostringstream ost;
    // 	ost << "GridFunctionCoefficientFunction: Element nr. "
    // 	    << elnr << " out of range 0 - " << gf.vec[0]->Size()-1 << endl;
    // 	throw Exception (ost.str());
    //       }

    //     return (gf.vec[0]->FV())(elnr);
    //   }
  
    virtual double Evaluate (const BaseSpecificIntegrationPoint & ip);


    //   {
    //     if (elnr < 0 || elnr >= gf.vec[0]->Size())
    //       {
    // 	ostringstream ost;
    // 	ost << "GridFunctionCoefficientFunction: Element nr. "
    // 	    << elnr << " out of range 0 - " << gf.vec[0]->Size()-1 << endl;
    // 	throw Exception (ost.str());
    //       }

    //     return (gf.vec[0]->FV())(elnr);
    //  }
  
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

    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages, int component = -1);
    void Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages_times_volumes, Array<double> & volumes, int component = -1);
    

  };

}




#endif
