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

  class NGS_DLL_HEADER SparseFactorization : public BaseMatrix
  { 
  protected:
    weak_ptr<const BaseSparseMatrix> matrix;
    shared_ptr<BitArray> inner;
    shared_ptr<const Array<int>> cluster;
    bool smooth_is_projection;

  public:
    SparseFactorization (shared_ptr<const BaseSparseMatrix> amatrix,
			 shared_ptr<BitArray> ainner,
			 shared_ptr<const Array<int>> acluster);
    SparseFactorization() {}

    void DoArchive(Archive& ar) override;

    virtual bool IsComplex() const override { return matrix.lock()->IsComplex(); }

    virtual void Smooth (BaseVector & u, const BaseVector & f, BaseVector & y) const;

    int VHeight() const override { return matrix.lock()->VWidth();}
    int VWidth() const override { return matrix.lock()->VHeight();}

    bool SmoothIsProjection () const { return smooth_is_projection; }
    
    auto GetAMatrix() const { return matrix.lock(); }
    virtual bool SupportsUpdate() const { return false; } 
  };








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
  class NGS_DLL_HEADER SparseCholeskyTM : public SparseFactorization
  {
  protected:
    // height of the matrix
    int height;
    // number of real unknowns
    int nused;
    // number of non-zero entries in the L-factor
    size_t nze;
    //
    bool hermitian = false;   // Hermitian or complex-symmetric ? 
    // the reordering (original dofnr i -> order[i])
    Array<int> order;
    Array<int> inv_order;
    
    // L-factor in compressed storage
    // Array<TM, size_t> lfact;
    NumaInterleavedArray<TM> lfact;

    // index-array to lfact
    Array<size_t> firstinrow;

    // diagonal 
    Array<TM> diag;


    // row-indices of non-zero entries
    // all row-indices within one block are identic, and stored just once
    Array<int> rowindex2;
    // index-array to rowindex
    Array<size_t> firstinrow_ri;
    
    // blocknr of dof
    Array<int> blocknrs;

    // block i has dofs  [blocks[i], blocks[i+1])
    Array<int> blocks; 

    // dependency graph for elimination
    Table<int> block_dependency; 

  public:      // needed for gcc 4.9, why  ??? 
    class MicroTask
    {
    public:
      int blocknr;
      enum BT { L_BLOCK, B_BLOCK, LB_BLOCK };
      BT type;
      int bblock;
      int nbblocks;
      template <typename ARCHIVE>
      void DoArchive(ARCHIVE& ar)
      {
        ar & blocknr & type & bblock & nbblocks;
      }
    };
  protected:
    
    Array<MicroTask> microtasks;
    Table<int> micro_dependency;     
    Table<int> micro_dependency_trans;     


    //
    MinimumDegreeOrdering * mdo;

    // maximal non-zero entries in a column
    int maxrow;

    // the original matrix
    // const SparseMatrixTM<TM> & mat;

    int max_bs = 1024;  
    int max_micro_bs = 256; // not yet used
    
  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL_MAT;

    ///
    SparseCholeskyTM (shared_ptr<const SparseMatrixTM<TM>> a,
                      shared_ptr<BitArray> ainner = nullptr,
                      shared_ptr<const Array<int>> acluster = nullptr,
                      bool allow_refactor = 0);
    SparseCholeskyTM() {}
    ///
    virtual ~SparseCholeskyTM ();
    ///
    int VHeight() const override { return height; }
    ///
    int VWidth() const override { return height; }
    ///
    BaseMatrix::OperatorInfo GetOperatorInfo () const override
    { return { string("SparseCholesky-")+typeid(TM).name(), size_t(height), size_t(height) }; }
    
    void SetHermitian (bool herm = true) { hermitian = herm; } 
    ///
    void Allocate (const Array<int> & aorder, 
		   const Array<MDOVertex> & vertices,
		   const int * blocknr);

    void DoArchive(Archive& ar) override;
    ///
    void Factor (); 
#ifdef LAPACK
    void FactorSPD (); 
    template <typename T>
    void FactorSPD1 (T dummy); 
#endif

    virtual bool SupportsUpdate() const override { return true; }
    virtual void Update() override
    {
      // FactorNew (dynamic_cast<const SparseMatrix<TM>&> (*matrix.lock().get()));
      auto castmatrix = dynamic_pointer_cast<const SparseMatrix<TM>>(matrix.lock());
      FactorNew (*castmatrix);
    }
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
    virtual ostream & Print (ostream & ost) const override;

    virtual Array<MemoryUsage> GetMemoryUsage () const override
    {
      return { MemoryUsage ("SparseChol", nze*sizeof(TM), 1) };
    }

    virtual size_t NZE () const override { return nze; }
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
      auto base = firstinrow_ri[range.First()] + range.Size()-1;
      auto ext_size =  firstinrow[range.First()+1]-firstinrow[range.First()] - range.Size()+1;
      return rowindex2.Range(base, base+ext_size);
    }


    FlatArray<MicroTask> GetMicroTasks() const { return microtasks; }
    FlatTable<int> GetMicroDependency() const { return micro_dependency; }
    FlatTable<int> GetMicroDependencyTranspose() const { return micro_dependency_trans; }
    FlatArray<int> GetBlocks() const { return blocks; }

    FlatArray<int> GetOrder() const { return order; }
    FlatArray<int> GetInvOrder() const { return inv_order; }

    FlatArray<int> GetRowIndex2() const { return rowindex2; }
    FlatArray<size_t> GetFirstInRowRI() const { return firstinrow_ri; }
    FlatArray<size_t> GetFirstInRow() const { return firstinrow; }

    FlatArray<TM> GetLFact() const { return lfact; }
    FlatArray<TM> GetDiag() const { return diag; }

    auto GetNUsed() const { return nused; }
  };





  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class NGS_DLL_HEADER SparseCholesky : public SparseCholeskyTM<TM>
  {
    typedef SparseCholeskyTM<TM> BASE;
    using BASE::height;
    using BASE::Height;
    using BASE::inner;
    using BASE::cluster;

    using BASE::lfact;
    using BASE::diag;
    using BASE::order;
    using BASE::inv_order;
    using BASE::firstinrow;

    using BASE::blocks;
    using typename BASE::MicroTask;
    using BASE::microtasks;
    using BASE::micro_dependency;
    using BASE::micro_dependency_trans;
    using BASE::block_dependency;
    using BASE::BlockDofs;
    using BASE::BlockExtDofs;
    using BASE::hermitian;
  public:
    typedef TV_COL TV;
    typedef TV_ROW TVX;
    typedef typename mat_traits<TV_ROW>::TSCAL TSCAL_VEC;

    
    SparseCholesky (shared_ptr<const SparseMatrixTM<TM>> a,
		    shared_ptr<BitArray> ainner = nullptr,
		    shared_ptr<const Array<int>> acluster = nullptr,
		    bool allow_refactor = 0)
      : SparseCholeskyTM<TM> (a, ainner, acluster, allow_refactor) { ; }
    SparseCholesky() {}

    ///
    virtual ~SparseCholesky () { ; }
    
    void Mult (const BaseVector & x, BaseVector & y) const override;

    void MultAdd (TSCAL_VEC s, const BaseVector & x, BaseVector & y) const override;
    void MultTransAdd (TSCAL_VEC s, const BaseVector & x, BaseVector & y) const override
    {
      MultAdd (s, x, y);
    }

    AutoVector CreateRowVector () const override { return make_unique<VVector<TV>> (height); }
    AutoVector CreateColVector () const override { return make_unique<VVector<TV>> (height); }

    void Smooth (BaseVector & u, const BaseVector & f, BaseVector & y) const override;

    void SolveBlock (int i, FlatVector<TV> hy) const;
    void SolveBlockT (int i, FlatVector<TV> hy) const;
  private:
    void SolveReordered(FlatVector<TVX> hy) const;
  };


}

#endif
