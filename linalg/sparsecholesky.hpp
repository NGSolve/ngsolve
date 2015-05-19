#ifndef FILE_SPARSECHOLESKY
#define FILE_SPARSECHOLESKY

/* *************************************************************************/
/* File:   sparsecholesky.hpp                                              */
/* Author: Joachim Schoeberl                                               */
/* Date:   18. Jun. 97                                                     */
/* *************************************************************************/

/*
  sparse cholesky factorization
*/

namespace ngla
{

  class SparseFactorization : public BaseMatrix
  { 
  protected:
    const BaseSparseMatrix & matrix;
    const BitArray * inner;
    const Array<int> * cluster;
    bool smooth_is_projection;

  public:
    SparseFactorization (const BaseSparseMatrix & amatrix,
			 const BitArray * ainner,
			 const Array<int> * acluster);
 
    virtual void Smooth (BaseVector & u, const BaseVector & f, BaseVector & y) const;

    int VHeight(){ return matrix.VWidth();}
    int VWidth(){ return matrix.VHeight();}

    bool SmoothIsProjection () const { return smooth_is_projection; }
  };



#ifdef USE_NUMA

template <typename T>
class NumaInterleavedArray : public Array<T,size_t>
{
  T * numa_ptr;
  size_t numa_size;
public:
  NumaInterleavedArray () { numa_size = 0; numa_ptr = nullptr; }
  NumaInterleavedArray (size_t s)
    : Array<T,size_t> (s, (T*)numa_alloc_interleaved(s*sizeof(T)))
  {
    numa_ptr = this->data;
    numa_size = s;

    /*
    int avail = numa_available();
    int num_nodes = numa_num_configured_nodes();
    size_t pagesize = numa_pagesize();
    
    int npages = ceil ( double(s)*sizeof(T) / pagesize );

    cout << "size = " << numa_size << endl;
    cout << "npages = " << npages << endl;

    for (int i = 0; i < num_nodes; i++)
      {
        int beg = (i * npages) / num_nodes;
        int end = ( (i+1) * npages) / num_nodes;
        cout << "node " << i << " : [" << beg << "-" << end << ")" << endl;
        numa_tonode_memory(numa_ptr+beg*pagesize/sizeof(T), (end-beg)*pagesize, i);
      }
    */
  }

  ~NumaInterleavedArray ()
  {
    numa_free (numa_ptr, numa_size*sizeof(T));
  }

  NumaInterleavedArray & operator= (T val)
  {
    Array<T,size_t>::operator= (val);      
    return *this;
  }

  NumaInterleavedArray & operator= (NumaInterleavedArray && a2)
  {
    Array<T,size_t>::operator= ((Array<T,size_t>&&)a2);  
    ngstd::Swap (numa_ptr, a2.numa_ptr);
    ngstd::Swap (numa_size, a2.numa_size);
    return *this;
  }

  void Swap (NumaInterleavedArray & b)
  {
    Array<T,size_t>::Swap(b);    
    ngstd::Swap (numa_ptr, b.numa_ptr);
    ngstd::Swap (numa_size, b.numa_size);
  }

  void SetSize (size_t size)
  {
    cerr << "************************* NumaDistArray::SetSize not overloaded" << endl;
    Array<T,size_t>::SetSize(size);
  }
};
#else

  template <typename T>
  using NumaInterleavedArray = Array<T>;
  
#endif








  /**
     A sparse cholesky factorization.
     The unknowns are reordered by the minimum degree
     ordering algorithm

     computs A = L D L^t
     L is stored column-wise
  */

  template<class TM>
	   // class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   // class TV_COL = typename mat_traits<TM>::TV_COL>
  class SparseCholeskyTM : public SparseFactorization
  {
  protected:
    // height of the matrix
    int height;

    // number of non-zero entries in the L-factor
    size_t nze;

    // the reordering (original dofnr i -> order[i])
    Array<int, int> order;
    
    // L-factor in compressed storage
    // Array<TM, size_t> lfact;
    NumaInterleavedArray<TM> lfact;

    // index-array to lfact
    Array<size_t, size_t> firstinrow;

    // diagonal 
    Array<TM, size_t> diag;


    // row-indices of non-zero entries
    // all row-indices within one block are identic, and stored just once
    Array<int, size_t> rowindex2;
    // index-array to rowindex
    Array<size_t, size_t> firstinrow_ri;
    
    // blocknr of dof
    Array<int> blocknrs;

    // block i has dofs  [blocks[i], blocks[i+1])
    Array<int> blocks; 

    // dependency graph for elimination
    Table<int> block_dependency; 


    class MicroTask
    {
    public:
      int blocknr;
      bool solveL;
      int bblock;
      int nbblocks;
    };
    
    Array<MicroTask> microtasks;
    Table<int> micro_dependency;     
    Table<int> micro_dependency_trans;     


    //
    MinimumDegreeOrdering * mdo;

    // maximal non-zero entries in a column
    int maxrow;

    // the original matrix
    const SparseMatrixTM<TM> & mat;

  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL_MAT;

    ///
    SparseCholeskyTM (const SparseMatrixTM<TM> & a, 
                      const BitArray * ainner = NULL,
                      const Array<int> * acluster = NULL,
                      bool allow_refactor = 0);
    ///
    virtual ~SparseCholeskyTM ();
    ///
    int VHeight() const { return height; }
    ///
    int VWidth() const { return height; }
    ///
    void Allocate (const Array<int> & aorder, 
		   const Array<MDOVertex> & vertices,
		   const int * blocknr);
    ///
    void Factor (); 
#ifdef LAPACK
    void FactorSPD (); 
#endif
    ///
    void FactorNew (const SparseMatrix<TM> & a);

    /**
       A = L+D+L^T
       y = f - (L+D)^T u
       w = C^{-1} (y - L u)
       u += w
       y -= (L+D)^T w
    **/
    // virtual void Smooth (BaseVector & u, const BaseVector & f, BaseVector & y) const;
    ///
    virtual ostream & Print (ostream & ost) const;

    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const
    {
      mu.Append (new MemoryUsageStruct ("SparseChol", nze*sizeof(TM), 1));
    }


    ///
    void Set (int i, int j, const TM & val);
    ///
    const TM & Get (int i, int j) const;
    ///
    void SetOrig (int i, int j, const TM & val)
    { Set (order[i], order[j], val); }


    // the dofs of block bnr
    IntRange BlockDofs (int bnr) const { return Range(blocks[bnr], blocks[bnr+1]); }

    // the external dofs of block bnr
    FlatArray<int> BlockExtDofs (int bnr) const
    {
      auto range = BlockDofs (bnr);
      auto base = firstinrow_ri[range.begin()] + range.Size()-1;
      auto ext_size =  firstinrow[range.begin()+1]-firstinrow[range.begin()] - range.Size()+1;
      return rowindex2.Range(base, base+ext_size);
    }
  };









  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class SparseCholesky : public SparseCholeskyTM<TM>
  {
    typedef SparseCholeskyTM<TM> BASE;
    using BASE::height;
    using BASE::Height;
    using BASE::inner;
    using BASE::cluster;

    using BASE::lfact;
    using BASE::diag;
    using BASE::order;
    using BASE::firstinrow;

    using BASE::blocks;
    using BASE::MicroTask;
    using BASE::microtasks;
    using BASE::micro_dependency;
    using BASE::micro_dependency_trans;
    using BASE::block_dependency;
    using BASE::BlockDofs;
    using BASE::BlockExtDofs;
  public:
    typedef TV_COL TV;
    typedef TV_ROW TVX;
    typedef typename mat_traits<TV_ROW>::TSCAL TSCAL_VEC;

    
    SparseCholesky (const SparseMatrixTM<TM> & a, 
		    const BitArray * ainner = NULL,
		    const Array<int> * acluster = NULL,
		    bool allow_refactor = 0)
      : SparseCholeskyTM<TM> (a, ainner, acluster, allow_refactor) { ; }

    ///
    virtual ~SparseCholesky () { ; }
    
    virtual void Mult (const BaseVector & x, BaseVector & y) const;

    virtual void MultAdd (TSCAL_VEC s, const BaseVector & x, BaseVector & y) const;

    virtual AutoVector CreateVector () const
    {
      return make_shared<VVector<TV>> (height);
    }


    void SolveBlock (int i, FlatVector<TV> hy) const;
    void SolveBlockT (int i, FlatVector<TV> hy) const;
  };


}

#endif
