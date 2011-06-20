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


  class ParallelBaseMatrix : virtual public BaseMatrix
  {
  protected:
    const ngparallel::ParallelDofs * paralleldofs;

  public:
    ///
    ParallelBaseMatrix () : paralleldofs(NULL) { ; }

    ///
    ParallelBaseMatrix ( const ParallelBaseMatrix & amat )
      : paralleldofs ( amat.GetParallelDofs() ) { ; }

    ///
    ParallelBaseMatrix ( const ngparallel::ParallelDofs * aparalleldofs )
      : paralleldofs ( aparalleldofs ) { ; }

    /// destructor
    virtual ~ParallelBaseMatrix ();


    virtual BaseMatrix * ConsistentMat () 
    { cerr << "ERROR -- ParallelBaseMatrix::ConsistentMat() called" << endl; return 0; }

    virtual const BaseMatrix * ConsistentMat () const  
    { cerr << "ERROR -- ParallelBaseMatrix::ConsistentMat() called" << endl; return 0; }

    virtual void SetConsistentMat ( BaseMatrix * aconsistentmat )
    { cerr << "ERROR -- ParallelBaseMatrix::SetConsistentMat called" << endl; }

    virtual void AllocateConsistentMat ()
    { cerr << "ERROR -- ParallelBaseMatrix::AllocateConsistentMat called" << endl; }

    virtual void  AllocateConsistentMat ( const class MatrixGraph & graph )
    { cerr << "ERROR -- ParallelBaseMatrix::AllocateConsistentMat(const MatrixGraph&) called" << endl; }

    virtual void CalcConsistentMat (LocalHeap & lh) 
    { cerr << "ERROR -- ParallelBaseMatrix::CalcConsistentMat called" << endl; }

    virtual void SetParallelDofs ( ngparallel::ParallelDofs * aparalleldofs )
    { paralleldofs = aparalleldofs; }


    virtual const ngparallel::ParallelDofs * GetParallelDofs ( ) const { return paralleldofs; }

    virtual bool IsParallelMatrix() const
    {
      (*testout) << "PARALLELDOFS " <<  paralleldofs << endl;
      if ( paralleldofs ) return true;
      else return false;
    }

  };








  template <typename TM>
  class MasterInverse : public ParallelBaseMatrix
  {
    BaseMatrix * inv;
    const BitArray * subset;
    DynamicTable<int> loc2glob;
    Array<int> select;
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

  };

  

#endif
}

#endif
