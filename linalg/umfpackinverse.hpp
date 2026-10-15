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
    // void GetUmfpackMatrix (const SparseMatrix<TM> & a, TSUBSET subset);
    void GetUmfpackMatrix (const S_BaseSparseMatrix<TSCAL> & a, TSUBSET subset);

    ///
    int VHeight() const override { return height/entrysize; }
    ///
    int VWidth() const override { return height/entrysize; }
    ///
    virtual ostream & Print (ostream & ost) const override;

    virtual bool SupportsUpdate() const override { return true; }     
    virtual void Update() override;

    virtual Array<MemoryUsage> GetMemoryUsage () const
 override {
      return { MemoryUsage ("Umfpack", nze*sizeof(TSCAL)*entrysize*entrysize, 1) };
    }
  };


  
  // the vector scalar is chosen at run time, a real matrix may solve for complex vectors
  template<class SCAL>
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
    Scalar vecscal = SCAL(0);
    
    using S_UmfpackInverse<SCAL>::S_UmfpackInverse;
    
    VecFormat RowFormat () const override
    { return VecFormat (height/entrysize, vecscal, entrysize); }
    VecFormat ColFormat () const override { return RowFormat(); }

    void Mult (const BaseVector & x, BaseVector & y) const override;
    void MultTrans (const BaseVector & x, BaseVector & y) const override;
  private:
    template <typename SCAL_VEC>
    void T_Mult (const BaseVector & x, BaseVector & y, bool trans) const;
  };
  

  template<class TM>
  class UmfpackInverse : public S_UmfpackInverse_SVec<typename mat_traits<TM>::TSCAL>
  {
  public:
    typedef S_UmfpackInverse_SVec<typename mat_traits<TM>::TSCAL> BASE;
    
    using typename BASE::TSCAL;
    
  public:
    UmfpackInverse (shared_ptr<const SparseMatrix<TM>> a,
		    shared_ptr<BitArray> ainner = nullptr,
		    shared_ptr<const Array<int>> acluster = nullptr,
		    int symmetric = 0)
      : BASE (a, ainner, acluster, symmetric)
    { this->vecscal = a->GetVectorScalar(); }

    // virtual ~UmfpackInverse () { ; }
  };

}

#endif
#endif
