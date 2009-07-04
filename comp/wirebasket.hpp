#ifndef FILE_HO_BIFORM
#define FILE_HO_BIFORM

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
  class ElementByElementMatrix : public BaseMatrix
  {
    Array<FlatMatrix<SCAL> > elmats;
    Array<FlatArray<int> > dnums;
    int height;
  public:
    ElementByElementMatrix (int h) 
    { height = h; }

    virtual int VHeight() const { return height; }
    virtual int VWidth() const { return height; }

    virtual BaseVector * CreateVector () const
    {
      return new VVector<double> (height);
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    void AddElementMatrix (const Array<int> & dnums1,
			   const Array<int> & dnums2,
			   const FlatMatrix<SCAL> & elmat);
    

    virtual BaseVector & AsVector() 
    {
      return *new VVector<double> (1);
    }

    const FlatMatrix<SCAL> GetElementMatrix( int elnum ) const
    {
      return elmats[elnum];
    }

    const FlatArray<int> GetElementDNums ( int elnum ) const
    {
      return dnums[elnum]; 
    }

    virtual ostream & Print (ostream & ost) const
    {
      ost << "Element-by-Element Matrix:" << endl;
      return ost;
    }

    virtual BaseBlockJacobiPrecond * 
    CreateBlockJacobiPrecond (Table<int> & blocks,
			      const BaseVector * constraint = 0, int * paralleloptions = 0) const;
//     { 
//       return new BlockJacobiPrecond_ElByEl<SCAL> (*this, blocks );
//     }


    virtual BaseMatrix * 
    InverseMatrix (BitArray * subset = 0) const;

  };


template<class SCAL>
class ElementByElement_BilinearForm : public S_BilinearForm<SCAL>
{
public:
  ElementByElement_BilinearForm (const FESpace & afespace, const string & aname,
		   const Flags & flags);
  virtual ~ElementByElement_BilinearForm ();
  
  virtual void AllocateMatrix ();
  virtual BaseVector * CreateVector() const;
  
  virtual void AddElementMatrix (const Array<int> & dnums1,
				 const Array<int> & dnums2,
				 const FlatMatrix<SCAL> & elmat,
				 bool inner_element, int elnr,
				 LocalHeap & lh);    


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
  };





  template <class SCAL>
  class BDDCPreconditioner : public Preconditioner
  {
    const HO_BilinearForm<SCAL> * bfa;
    BaseMatrix * pre;
    string inversetype;
  public:
    BDDCPreconditioner (const PDE * pde, const Flags & aflags,
			const std::string aname = "bddcprecond");

    virtual ~BDDCPreconditioner() 
    { ; }

    static Preconditioner * Create (const PDE & pde, const Flags & flags)
    {
      return new BDDCPreconditioner<SCAL> (&pde, flags);
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
  };







}

#endif
