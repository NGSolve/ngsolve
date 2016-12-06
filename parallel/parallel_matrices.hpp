
#ifndef FILE_NGS_PARALLEL_MATRICES
#define FILE_NGS_PARALLEL_MATRICES

/* ************************************************************************/
/* File:   parallelmatrices.hpp                                           */
/* Author: Astrid Sinwel, Joachim Schoeberl                               */
/* Date:   2007,2011                                                      */
/* ************************************************************************/

namespace ngla
{

#ifdef PARALLEL


  template <typename TM>
  class MasterInverse : public BaseMatrix
  {
    shared_ptr<BaseMatrix> inv;
    shared_ptr<BitArray> subset;
    DynamicTable<int> loc2glob;
    Array<int> select;
    string invtype;
    const ParallelDofs * pardofs;
  public:
    MasterInverse (const SparseMatrixTM<TM> & mat, shared_ptr<BitArray> asubset, 
		   const ParallelDofs * apardofs);
    virtual ~MasterInverse ();
    virtual bool IsComplex() const { return inv->IsComplex(); } 
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    virtual int VHeight() const { return pardofs->GetNDofLocal(); }
    virtual int VWidth() const { return pardofs->GetNDofLocal(); }
  };


  class ParallelMatrix : public BaseMatrix
  {
    shared_ptr<BaseMatrix> mat;
    // const ParallelDofs & pardofs;
  public:
    ParallelMatrix (shared_ptr<BaseMatrix> amat, const ParallelDofs * apardofs);
    // : mat(*amat), pardofs(*apardofs) 
    // {const_cast<BaseMatrix&>(mat).SetParallelDofs (apardofs);}

    virtual ~ParallelMatrix ();
    virtual bool IsComplex() const { return mat->IsComplex(); } 
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;

    virtual BaseVector & AsVector() { return mat->AsVector(); }
    virtual const BaseVector & AsVector() const { return mat->AsVector(); }

    BaseMatrix & GetMatrix() const { return const_cast<BaseMatrix&> (*mat); }
    virtual shared_ptr<BaseMatrix> CreateMatrix () const;
    virtual AutoVector CreateVector () const;

    virtual ostream & Print (ostream & ost) const;

    virtual int VHeight() const;
    virtual int VWidth() const;

    // virtual const ParallelDofs * GetParallelDofs () const {return &pardofs;}


    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = 0) const;
    template <typename TM>
    shared_ptr<BaseMatrix> InverseMatrixTM (shared_ptr<BitArray> subset = 0) const;

    virtual shared_ptr<BaseMatrix> InverseMatrix (const Array<int> * clusters) const;
    virtual INVERSETYPE SetInverseType ( INVERSETYPE ainversetype ) const;
    virtual INVERSETYPE SetInverseType ( string ainversetype ) const;
    virtual INVERSETYPE  GetInverseType () const;
  };

  

#endif
}

#endif
