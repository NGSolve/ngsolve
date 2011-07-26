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
    BaseMatrix * inv;
    const BitArray * subset;
    DynamicTable<int> loc2glob;
    Array<int> select;
    string invtype;
  public:
    MasterInverse (const SparseMatrixTM<TM> & mat, const BitArray * asubset, const ParallelDofs * pardofs);
    virtual ~MasterInverse ();
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
  };


  class ParallelMatrix : public BaseMatrix
  {
    const BaseMatrix & mat;
    const ParallelDofs & pardofs;
  public:
    ParallelMatrix (const BaseMatrix & amat, const ParallelDofs & apardofs)
      : mat(amat), pardofs(apardofs) { ; }

    virtual ~ParallelMatrix ();
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;

    virtual BaseVector & AsVector() { return const_cast<BaseMatrix&> (mat).AsVector(); }
    virtual const BaseVector & AsVector() const { return mat.AsVector(); }

    BaseMatrix & GetMatrix() { return const_cast<BaseMatrix&> (mat); }
    virtual BaseVector * CreateVector () const;

    virtual ostream & Print (ostream & ost) const;

    virtual int VHeight() const;
    virtual int VWidth() const;

    virtual BaseMatrix * InverseMatrix (const BitArray * subset = 0) const;
    virtual BaseMatrix * InverseMatrix (const Array<int> * clusters) const;
    virtual INVERSETYPE SetInverseType ( INVERSETYPE ainversetype ) const;
    virtual INVERSETYPE SetInverseType ( string ainversetype ) const;
    virtual INVERSETYPE  GetInverseType () const;
  };

  

#endif
}

#endif
