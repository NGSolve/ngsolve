#ifndef FILE_PARDISOINVERSE
#define FILE_PARDISOINVERSE

/* *************************************************************************/
/* File:   pardisoinverse.hpp                                              */
/* Author: Florian Bachinger                                               */
/* Date:   Feb. 04                                                         */
/* *************************************************************************/

/*
  interface to the PARDISO - package
*/



////////////////////////////////////////////////////////////////////////////////

#include "sparsecholesky.hpp"  // for SparseFactorization

namespace ngla
{
  using ngbla::integer;
  extern bool is_pardiso_available;


  /*
  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class PardisoInverse : public SparseFactorization
  {
    integer height;             // matrix size in scalars
    integer compressed_height;  // matrix size after compression in scalars
    integer nze, entrysize;
    bool print;

    integer pt[128];
    integer hparams[64];

    Array<integer,size_t> rowstart, indices; 
    Array<typename mat_traits<TM>::TSCAL,size_t> matrix;

    integer matrixtype;
    bool symmetric, spd;

    void SetMatrixType(); 

    bool compressed;
    Array<int> compress;
  
  public:
    typedef TV_COL TV;
    typedef TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    ///
    PardisoInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
		    const BitArray * ainner = NULL,
		    const Array<int> * acluster = NULL,
		    int symmetric = 0);
    ///

    template <typename TSUBSET>
    void GetPardisoMatrix (const SparseMatrix<TM,TV_ROW,TV_COL> & a, TSUBSET subset);

    virtual ~PardisoInverse ();
    ///
    int VHeight() const { return height/entrysize; }
    ///
    int VWidth() const { return height/entrysize; }
    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const;


    ///
    virtual ostream & Print (ostream & ost) const;

    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const
    {
      mu.Append (new MemoryUsageStruct ("SparseChol", nze*sizeof(TM), 1));
    }
    virtual BaseVector * CreateVector () const
    {
      return new VVector<TV> (height/entrysize);
    }
  };
  */









  template<class TM>
  class PardisoInverseTM : public SparseFactorization
  {
  protected:
    integer height;             // matrix size in scalars
    integer compressed_height;  // matrix size after compression in scalars
    integer nze, entrysize;
    bool print;

    integer pt[128];
    integer hparams[64];

    Array<integer> rowstart, indices; 
    Array<typename mat_traits<TM>::TSCAL> matrix;

    integer matrixtype;
    bool symmetric, spd;

    void SetMatrixType(); 

    bool compressed;
    Array<int> compress;
    size_t memory_allocated_in_pardiso_lib = 0;
  
  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    ///
    PardisoInverseTM (shared_ptr<const SparseMatrixTM<TM>> a,
		      shared_ptr<BitArray> ainner = nullptr,
		      shared_ptr<const Array<int>> acluster = nullptr,
		      int symmetric = 0);
    ///

    template <typename TSUBSET>
    void GetPardisoMatrix (const SparseMatrixTM<TM> & a, TSUBSET subset);

    virtual ~PardisoInverseTM ();
    ///
    int VHeight() const { return height/entrysize; }
    ///
    int VWidth() const { return height/entrysize; }
    ///
    virtual ostream & Print (ostream & ost) const;

    virtual Array<MemoryUsage> GetMemoryUsage () const
    {
      return { MemoryUsage ("Pardiso", nze*sizeof(TM), 1) };
    }
  };


  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class PardisoInverse : public PardisoInverseTM<TM>
  {
    using PardisoInverseTM<TM>::height;
    using PardisoInverseTM<TM>::compressed_height;
    using PardisoInverseTM<TM>::entrysize;
    using PardisoInverseTM<TM>::pt;
    using PardisoInverseTM<TM>::hparams;
    using PardisoInverseTM<TM>::matrixtype;
    using PardisoInverseTM<TM>::rowstart;
    using PardisoInverseTM<TM>::indices;
    using PardisoInverseTM<TM>::matrix;
    using PardisoInverseTM<TM>::compressed;
    using PardisoInverseTM<TM>::compress;

  public:
    typedef TV_COL TV;
    typedef TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    ///
    PardisoInverse (shared_ptr<const SparseMatrix<TM,TV_ROW,TV_COL>> a,
		    shared_ptr<BitArray> ainner = nullptr,
		    shared_ptr<const Array<int>> acluster = nullptr,
		    int symmetric = 0)
      : PardisoInverseTM<TM> (a, ainner, acluster, symmetric) { ; }

    virtual ~PardisoInverse () { ; }
    ///
    void Mult (const BaseVector & x, BaseVector & y) const override;
    void MultTrans (const BaseVector & x, BaseVector & y) const override;
    ///

    AutoVector CreateRowVector() const override
    {
      return make_unique<VVector<TV>> (height/entrysize);
    }
    AutoVector CreateColVector () const override
    {
      return make_unique<VVector<TV>> (height/entrysize);
    }
  };





}

#endif
