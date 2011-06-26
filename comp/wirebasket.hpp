#ifndef FILE_HO_BIFORM
#define FILE_HO_BIFORM

wer braucht das ???


/*********************************************************************/
/* File:   ho_biform.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   06. Dec. 2004                                             */
/*********************************************************************/

namespace ngcomp
{

template<class SCAL>
class HO_BilinearForm : public S_BilinearForm<SCAL>
{
  SparseMatrix<SCAL> * inexact_schur;
public:
  HO_BilinearForm (const FESpace & afespace, const string & aname,
		   const Flags & flags);
  virtual ~HO_BilinearForm ();
  
  virtual void AllocateMatrix ();
  virtual BaseVector * CreateVector() const;
  
  virtual void AddElementMatrix (const Array<int> & dnums1,
				 const Array<int> & dnums2,
				 const FlatMatrix<SCAL> & elmat,
				 bool inner_element, int elnr,
				 LocalHeap & lh);    

  const SparseMatrix<SCAL> & InexactSchur() const { return *inexact_schur; }

#ifdef PARALLEL
  virtual void AllocateConsistentMatrix ()
  { 
    ; 
  }

  virtual void BuildConsistentMatrix (LocalHeap & lh) 
  { ; }
#endif

};


  template <class SCAL>
  class WireBasketPreconditioner : public Preconditioner
  {
    const HO_BilinearForm<SCAL> * bfa;
    BaseMatrix * pre;
    string inversetype;
    int smoothingtype;
  public:
    WireBasketPreconditioner (const PDE * pde, const Flags & aflags, 
			      const std::string aname = "wirebasketprecond");

    virtual ~WireBasketPreconditioner() 
    { ; }

    static Preconditioner * Create (const PDE & pde, const Flags & flags)
    {
      return new WireBasketPreconditioner<SCAL> (&pde, flags);
    }

    virtual void Update ();

    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa->GetMatrix();
    }

    virtual const BaseMatrix & GetMatrix() const
    {
      return *pre;
    }

    virtual const char * ClassName() const
    { return "Wire-basket Preconditioner"; }
  };

}

#endif
