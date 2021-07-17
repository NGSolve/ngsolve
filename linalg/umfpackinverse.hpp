#ifndef FILE_UMFPACKINVERSE
#define FILE_UMFPACKINVERSE

#ifdef USE_UMFPACK

/* *************************************************************************/
/* File:   umfpackinverse.hpp                                              */
/* Author: Matthias Hochsteger                                             */
/* Date:   Nov. 15                                                         */
/* *************************************************************************/

/*
  interface to the UMFPACK - package
*/



////////////////////////////////////////////////////////////////////////////////

#include <umfpack.h>

namespace ngla
{

  template<class TM>
  class UmfpackInverseTM : public SparseFactorization
  {
  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;

  protected:
    typedef SuiteSparse_long suite_long;

    int height;             // matrix size in scalars
    int compressed_height;  // matrix size after compression in scalars
    suite_long nze, entrysize;
    bool print;

    void *Symbolic = nullptr;
    void *Numeric = nullptr;

    Array<suite_long> rowstart, indices;
    Array<TSCAL> values;

    bool symmetric, is_complex;

    void SetMatrixType();

    bool compressed;
    Array<int> compress;

  public:

    ///
    UmfpackInverseTM (shared_ptr<const SparseMatrixTM<TM>> a,
		      shared_ptr<BitArray> ainner = nullptr,
		      shared_ptr<const Array<int>> acluster = nullptr,
		      int symmetric = 0);
    ///

    template <typename TSUBSET>
    void GetUmfpackMatrix (const SparseMatrixTM<TM> & a, TSUBSET subset);

    virtual ~UmfpackInverseTM ();
    ///
    int VHeight() const { return height/entrysize; }
    ///
    int VWidth() const { return height/entrysize; }
    ///
    virtual ostream & Print (ostream & ost) const;

    virtual bool SupportsUpdate() const { return true; }     
    virtual void Update();

    virtual Array<MemoryUsage> GetMemoryUsage () const
    {
      return { MemoryUsage ("Umfpack", nze*sizeof(TM), 1) };
    }
  };


  template<class TM,
	   class TV_ROW = typename mat_traits<TM>::TV_ROW,
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class UmfpackInverse : public UmfpackInverseTM<TM>
  {
    using UmfpackInverseTM<TM>::height;
    using UmfpackInverseTM<TM>::is_complex;
    using UmfpackInverseTM<TM>::compressed_height;
    using UmfpackInverseTM<TM>::entrysize;
    using UmfpackInverseTM<TM>::rowstart;
    using UmfpackInverseTM<TM>::indices;
    using UmfpackInverseTM<TM>::compressed;
    using UmfpackInverseTM<TM>::compress;

  public:
    typedef TV_COL TV;
    typedef TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    ///
    UmfpackInverse (shared_ptr<const SparseMatrix<TM,TV_ROW,TV_COL>> a,
		    shared_ptr<BitArray> ainner = nullptr,
		    shared_ptr<const Array<int>> acluster = nullptr,
		    int symmetric = 0)
      : UmfpackInverseTM<TM> (a, ainner, acluster, symmetric) { ; }

    virtual ~UmfpackInverse () { ; }
    ///
    void Mult (const BaseVector & x, BaseVector & y) const override;
    void MultTrans (const BaseVector & x, BaseVector & y) const override;
    ///
    AutoVector CreateRowVector () const override { return make_unique<VVector<TV>> (height/entrysize); }
    AutoVector CreateColVector () const override { return make_unique<VVector<TV>> (height/entrysize); }
  };

}

#endif
#endif
