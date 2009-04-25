#ifndef FILE_LINEARFORM
#define FILE_LINEARFORM

/*********************************************************************/
/* File:   linearform.hh                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngcomp
{

  /** 
      Linearform
  */

  class LinearForm : public NGS_Object
  {
  protected:
    ///
    const FESpace & fespace;
    ///
    Array<LinearFormIntegrator*> parts;
    ///
    Array<bool> parts_deletable;
    /// do the integration on independent meshes
    bool independent;

    bool print;
    bool printelvec;

    bool assembled;
    bool initialassembling;

    int cacheblocksize;
  public:
    ///
    LinearForm (const FESpace & afespace, 
		const string & aname, const Flags & flags);

    ///
    virtual ~LinearForm ();
  
    ///
    const FESpace & GetFESpace() const
    { return fespace; }


    ///
    virtual void AddIntegrator (LinearFormIntegrator * lfi, const bool deletable = true);

    ///
    virtual const LinearFormIntegrator * GetIntegrator (int i) const
    {
      return parts[i]; 
    }

    virtual LinearFormIntegrator * GetIntegrator (int i) 
    {
      return parts[i]; 
    }

    ///
    virtual int NumIntegrators () const
    {
      return parts.Size(); 
    }


    void SetIndependent (int aindependent = true)
    { independent = aindependent; }


    ///
    virtual void Assemble (LocalHeap & lh) = 0;

    virtual void CleanUpLevel() { ; }

    virtual bool IsAssembled (void);
    bool InitialAssembling (void);
    void SetNoInitialAssembling (void);


    ///
    virtual BaseVector & GetVector () const = 0;
    ///
    virtual string GetClassName () const
    {
      return "LinearForm";
    }

    void SetPrint (bool ap);
    void SetPrintElmat (bool ap);

    ///
    virtual void PrintReport (ostream & ost);
    ///
    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;

    virtual void SetCacheBlockSize (const int size)
    {
      cacheblocksize = size;
    }
  };


  template <class SCAL>
  class S_LinearForm : public LinearForm
  {
  public:
    typedef SCAL TSCAL;

    ///
    S_LinearForm (const FESpace & afespace, 
		  const string & aname,
		  const Flags & flags)
      : LinearForm (afespace, aname, flags)
    {
      ;
    }
  
    ///
    virtual void AllocateVector () = 0;

    ///
    virtual void AddElementVector (const Array<int> & dnums,
				   const FlatVector<SCAL> & elvec,
				   const int cachecomp = -1) = 0;
    virtual void SetElementVector (const Array<int> & dnums,
				   const FlatVector<SCAL> & elvec) = 0;
    virtual void GetElementVector (const Array<int> & dnums,
				   FlatVector<SCAL> & elvec) const = 0;

    ///
    virtual void Assemble (LocalHeap & lh);
    void AssembleIndependent (LocalHeap & lh);
  };






  /// Template argument specifies vector type
  template <class TV>
  class T_LinearForm : public S_LinearForm<typename mat_traits<TV>::TSCAL>
  {
    ///
    ngla::VVector<TV> * vec;

  public:

    typedef typename mat_traits<TV>::TSCAL TSCAL;
    enum { HEIGHT = mat_traits<TV>::HEIGHT };

    ///
    T_LinearForm (const FESpace & afespace, const string & aname,
		  const Flags & flags);

    ///
    virtual ~T_LinearForm ();

    ///
    virtual BaseVector & GetVector () const ;

    ///
    virtual void AllocateVector ();
    ///
    virtual void CleanUpLevel();

    ///
    virtual void AddElementVector (const Array<int> & dnums,
				   const FlatVector<TSCAL> & elvec,
				   const int cachecomp = -1);
    virtual void SetElementVector (const Array<int> & dnums,
				   const FlatVector<TSCAL> & elvec);
    virtual void GetElementVector (const Array<int> & dnums,
				   FlatVector<TSCAL> & elvec) const;
  };


  extern LinearForm * CreateLinearForm (const FESpace * space,
					const string & name,
					const Flags & flags);

}

#endif


