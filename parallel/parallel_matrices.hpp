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
    shared_ptr<ParallelDofs> pardofs;
  public:
    MasterInverse (const SparseMatrixTM<TM> & mat, shared_ptr<BitArray> asubset, 
		   shared_ptr<ParallelDofs> apardofs);
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
    ParallelMatrix (shared_ptr<BaseMatrix> amat, shared_ptr<ParallelDofs> apardofs);
    // : mat(*amat), pardofs(*apardofs) 
    // {const_cast<BaseMatrix&>(mat).SetParallelDofs (apardofs);}

    virtual ~ParallelMatrix () override;
    virtual bool IsComplex() const override { return mat->IsComplex(); } 
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual BaseVector & AsVector() override { return mat->AsVector(); }
    virtual const BaseVector & AsVector() const override { return mat->AsVector(); }

    shared_ptr<BaseMatrix> GetMatrix() const { return mat; }
    virtual shared_ptr<BaseMatrix> CreateMatrix () const override;
    virtual AutoVector CreateVector () const override;

    virtual ostream & Print (ostream & ost) const override;

    virtual int VHeight() const override;
    virtual int VWidth() const override;

    // virtual const ParallelDofs * GetParallelDofs () const {return &pardofs;}


    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = 0) const override;
    template <typename TM>
    shared_ptr<BaseMatrix> InverseMatrixTM (shared_ptr<BitArray> subset = 0) const;

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<const Array<int>> clusters) const override;
    virtual INVERSETYPE SetInverseType ( INVERSETYPE ainversetype ) const override;
    virtual INVERSETYPE SetInverseType ( string ainversetype ) const override;
    virtual INVERSETYPE  GetInverseType () const override;
  };

  class FETI_Jump_Matrix : public BaseMatrix
  {
  public:
    FETI_Jump_Matrix (shared_ptr<ParallelDofs> pardofs)
      : BaseMatrix(pardofs) { ; }

    virtual bool IsComplex() const { return false; }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;
  };

#endif
}

#endif
