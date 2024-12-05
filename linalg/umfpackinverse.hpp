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
#include "sparsecholesky.hpp"   // for SparseFactorization
namespace ngla
{

  template<class SCAL>
  class S_UmfpackInverse : public SparseFactorization
  {
  public:
    // typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef SCAL TSCAL;

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

    bool symmetric;
    static constexpr bool is_complex = ngbla::IsComplex<TSCAL>();

    void SetMatrixType();

    bool compressed;
    Array<int> compress;

  public:

    ///
    S_UmfpackInverse (shared_ptr<const S_BaseSparseMatrix<TSCAL>> a,
                      shared_ptr<BitArray> ainner = nullptr,
		      shared_ptr<const Array<int>> acluster = nullptr,
		      int symmetric = 0);
    virtual ~S_UmfpackInverse();
    ///

    template <typename TSUBSET>
    // void GetUmfpackMatrix (const SparseMatrixTM<TM> & a, TSUBSET subset);
    void GetUmfpackMatrix (const S_BaseSparseMatrix<TSCAL> & a, TSUBSET subset);

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
      return { MemoryUsage ("Umfpack", nze*sizeof(TSCAL)*entrysize*entrysize, 1) };
    }
  };


  
  template<class SCAL, class SCAL_VEC>
  class S_UmfpackInverse_SVec : public S_UmfpackInverse<SCAL>
  {
    typedef S_UmfpackInverse<SCAL> BASE;
  protected:
    using typename BASE::TSCAL;    
    using BASE::height;
    using BASE::is_complex;
    using BASE::compressed_height;
    using BASE::entrysize;
    using BASE::rowstart;
    using BASE::indices;
    using BASE::compressed;
    using BASE::compress;
    
    using S_UmfpackInverse<SCAL>::S_UmfpackInverse;
    static constexpr bool is_vector_complex = ngbla::IsComplex<SCAL_VEC>();
    
    AutoVector CreateRowVector () const override { return CreateBaseVector(height/entrysize, is_vector_complex, entrysize); }
    AutoVector CreateColVector () const override { return CreateBaseVector(height/entrysize, is_vector_complex, entrysize); }

    void Mult (const BaseVector & x, BaseVector & y) const override;
    void MultTrans (const BaseVector & x, BaseVector & y) const override;
  };
  

  template<class TM,
	   class TV_ROW = typename mat_traits<TM>::TV_ROW,
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class UmfpackInverse : public S_UmfpackInverse_SVec<typename mat_traits<TM>::TSCAL, typename mat_traits<TV_COL>::TSCAL>
  {
  public:
    typedef S_UmfpackInverse_SVec<typename mat_traits<TM>::TSCAL, typename mat_traits<TV_COL>::TSCAL> BASE;
    
    using typename BASE::TSCAL;
    typedef TV_COL TV;
    typedef TV_ROW TVX;
    
  private:
    using BASE::height;
    using BASE::is_complex;
    using BASE::is_vector_complex;    
    using BASE::compressed_height;
    using BASE::entrysize;
    using BASE::rowstart;
    using BASE::indices;
    using BASE::compressed;
    using BASE::compress;

  public:
    UmfpackInverse (shared_ptr<const SparseMatrix<TM,TV_ROW,TV_COL>> a,
		    shared_ptr<BitArray> ainner = nullptr,
		    shared_ptr<const Array<int>> acluster = nullptr,
		    int symmetric = 0)
      : BASE (a, ainner, acluster, symmetric) { ; }

    // virtual ~UmfpackInverse () { ; }
  };

}

#endif
#endif
