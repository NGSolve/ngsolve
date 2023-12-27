#ifndef FILE_LINEARFORM
#define FILE_LINEARFORM

/*********************************************************************/
/* File:   linearform.hh                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#include "fespace.hpp"

namespace ngcomp
{

  /** 
      Linearform
  */

  
  class NGS_DLL_HEADER LinearForm : public NGS_Object
  {
  protected:
    ///
    shared_ptr<FESpace> fespace;
    ///
    Array<shared_ptr<LinearFormIntegrator>> parts;
    Array<shared_ptr<LinearFormIntegrator>> VB_parts[4];
    Array<shared_ptr<class PointEvaluationFunctional>> pnteval;
    /// do the integration on independent meshes
    bool independent;
    /// print the assembled vector to testout
    bool print;
    /// print element vectos to testout
    bool printelvec;

    bool allocated;
    bool assembled;
    bool initialassembling;

    int cacheblocksize;
    /// output of norm of matrix entries
    bool checksum;

  public:
    ///
    LinearForm (shared_ptr<FESpace> afespace, 
		const string & aname, const Flags & flags);

    ///
    virtual ~LinearForm () { ; }
  
    ///
    shared_ptr<FESpace> GetFESpace() const { return fespace; }

    ///
    virtual LinearForm & AddIntegrator (shared_ptr<LinearFormIntegrator> lfi);

    LinearForm & operator+= (shared_ptr<LinearFormIntegrator> lfi)
    {
      return AddIntegrator(lfi);
    }

    ///
    shared_ptr<LinearFormIntegrator> GetIntegrator (int i) const
    {
      return parts[i]; 
    }

    const Array<shared_ptr<LinearFormIntegrator>> & Integrators() const
    {
      return parts;
    }

    ///
    int NumIntegrators () const
    {
      return parts.Size(); 
    }

    LinearForm & operator+= (shared_ptr<PointEvaluationFunctional> _pnteval)
      {
        pnteval += _pnteval;
        return *this;
      }
    
    void SetIndependent (bool aindependent = true)
    { 
      independent = aindependent; 
    }

    ///
    virtual void Assemble (LocalHeap & lh) = 0;
    ///
    virtual void AllocateVector () = 0;

    virtual void CleanUpLevel() { ; }

    virtual bool IsAssembled () { return assembled; }
    bool InitialAssembling () { return initialassembling; }
    void SetNoInitialAssembling () { initialassembling = false; }


    ///
    virtual shared_ptr<BaseVector> GetVectorPtr() const = 0;
    virtual BaseVector & GetVector () const = 0;

    ///
    virtual string GetClassName () const
    {
      return "LinearForm";
    }

    void SetPrint (bool ap);
    void SetPrintElmat (bool ap);

    ///
    virtual void PrintReport (ostream & ost) const;
    ///
    virtual Array<MemoryUsage> GetMemoryUsage () const;
      
    virtual void AddElementVector (FlatArray<int> dnums,
                                   FlatVector<double> elvec,
				   int cachecomp = -1);
    virtual void SetElementVector (FlatArray<int> dnums,
				   FlatVector<double> elvec);
    virtual void GetElementVector (FlatArray<int> dnums,
				   FlatVector<double> elvec) const;

    virtual void AddElementVector (FlatArray<int> dnums,
				   FlatVector<Complex> elvec,
				   int cachecomp = -1);
    virtual void SetElementVector (FlatArray<int> dnums,
                                   FlatVector<Complex> elvec);
    virtual void GetElementVector (FlatArray<int> dnums,
				   FlatVector<Complex> elvec) const;



    virtual void SetCacheBlockSize (const int size)
    {
      cacheblocksize = size;
    }
  };


  template <class SCAL>
  class NGS_DLL_HEADER S_LinearForm : public LinearForm
  {
  protected:
    shared_ptr<BaseVector> vec;    
    
  public:
    typedef SCAL TSCAL;
    using LinearForm::LinearForm;

    virtual BaseVector & GetVector () const override { return *vec; }
    virtual shared_ptr<BaseVector> GetVectorPtr () const override { return vec; }

    virtual void AllocateVector () override;
    virtual void CleanUpLevel() override;    
    ///
    virtual void AddElementVector (FlatArray<int> dnums,
                                   FlatVector<SCAL> elvec,
                                   int cachecomp = -1) override = 0;
    virtual void SetElementVector (FlatArray<int> dnums,
                                   FlatVector<SCAL> elvec) override;
    virtual void GetElementVector (FlatArray<int> dnums,
				   FlatVector<SCAL> elvec) const override;

    ///
    virtual void Assemble (LocalHeap & lh) override;
    void AssembleIndependent (LocalHeap & lh);
  };






  /// Template argument specifies vector type
  template <class TV>
  class NGS_DLL_HEADER T_LinearForm : public S_LinearForm<typename mat_traits<TV>::TSCAL>
  {
    ///
    // shared_ptr<VVector<TV>> vec;
    typedef S_LinearForm<typename mat_traits<TV>::TSCAL> BASE;
    using BASE::vec;
  public:

    typedef typename mat_traits<TV>::TSCAL TSCAL;
    // enum { HEIGHT = mat_traits<TV>::HEIGHT };

    ///
    using S_LinearForm<TSCAL>::S_LinearForm;

    ///
    virtual ~T_LinearForm () { };

    ///
    // virtual void AllocateVector () override;
    ///
    // virtual void CleanUpLevel() override;

    ///
    virtual void AddElementVector (FlatArray<int> dnums,
				   FlatVector<TSCAL> elvec,
                                   int cachecomp = -1) override;
    /*
    virtual void SetElementVector (FlatArray<int> dnums,
                                   FlatVector<TSCAL> elvec) override;
    virtual void GetElementVector (FlatArray<int> dnums,
				   FlatVector<TSCAL> elvec) const override;
    */
  };





  class ComponentLinearForm : public LinearForm
  {
    shared_ptr<LinearForm> base_lf;
    int comp;
  public:
    ComponentLinearForm (shared_ptr<LinearForm> abase_lf, int acomp, int ancomp);
    virtual LinearForm & AddIntegrator (shared_ptr<LinearFormIntegrator> lfi);

    virtual void AllocateVector () { cerr << "comp - allocate is illegal" << endl; }
    virtual void Assemble (LocalHeap & lh) { cerr << "comp - assemble is illegal" << endl; }
    virtual shared_ptr<BaseVector> GetVectorPtr() const
    {
      auto fes = dynamic_pointer_cast<CompoundFESpace> (base_lf->GetFESpace());
      return base_lf->GetVectorPtr()->Range(fes->GetRange(comp));
      // throw Exception ("comp - GetVectorPtr is illegal");
    }
    virtual BaseVector & GetVector () const 
    { throw Exception ("comp - GetVector is illegal"); }
  };



  extern NGS_DLL_HEADER shared_ptr<LinearForm> CreateLinearForm (shared_ptr<FESpace> space,
                                                                 const string & name,
                                                                 const Flags & flags);


  class PointEvaluationFunctional
  {
  public:
    shared_ptr<CoefficientFunction> cf;
    Vector<> point;
  public:
    PointEvaluationFunctional (shared_ptr<CoefficientFunction> acf,
                               Vector<> apoint)
      : cf(acf), point(apoint) { }
    
    SparseVector<double> Assemble() const;
  };
  

}

#endif


