#ifndef FILE_BLOCK_JACOBI
#define FILE_BLOCK_JACOBI

/* *************************************************************************/
/* File:   blockjacobi.hpp                                                  */
/* Author: Joachim Schoeberl                                               */
/* Date:   06. Oct. 96                                                     */
/* *************************************************************************/

namespace ngla
{

  /**
     Base class for Block - Jacobi and Block Gauss Seidel smoother.
  */
  class NGS_DLL_HEADER BaseBlockJacobiPrecond : virtual public BaseMatrix
  {
  protected:
    /// the table defining the blocks
    shared_ptr<Table<int>> blocktable;
    /// maximal block size
    int maxbs;

    /// block coloring 
    Table<int> block_coloring;

    /// balancing for each color
    Array<Partitioning> color_balance;

    size_t nze;
  public:
    /// the blocktable define the blocks. ATTENTION: entries will be reordered !
    BaseBlockJacobiPrecond (shared_ptr<Table<int>> ablocktable);

    /// deletes the table
    virtual ~BaseBlockJacobiPrecond ();

    /// performs steps Gauss-Seidel steps for the equation A x = b
    virtual void GSSmooth (BaseVector & x, const BaseVector & b,
			   int steps = 1) const = 0;

    /// performs steps Gauss-Seidel steps for the equation A x = b with partial residual y
    virtual void GSSmoothPartial (BaseVector & x, const BaseVector & b, BaseVector & y) const 
    {
      GSSmooth (x, b, 1);
    }


    /// does smoothing. The vector res contains the residuum (b-Ax) before and after the smoothing
    virtual void GSSmoothResiduum (BaseVector & x, const BaseVector & b,
				   BaseVector & res, int steps = 1) const = 0;

    /// does smoothing in reversed order
    virtual void GSSmoothBack (BaseVector & x, const BaseVector & b,
			       int steps = 1) const = 0;

    virtual void GSSmoothBackPartial (BaseVector & x, const BaseVector & b, BaseVector & y) const 
    {
      GSSmoothBack (x, b, 1);
    }


    /// reorders block entries for band-width minimization
    int Reorder (FlatArray<int> block, const MatrixGraph & graph,
		 FlatArray<int> usedflags,        // in and out: array of -1, size = graph.size
		 LocalHeap & lh);

    auto GetBlockTable() const { return blocktable; } 
  };


  class SymmetricBlockGaussSeidelPrecond : virtual public BaseMatrix
  {
    shared_ptr<BaseBlockJacobiPrecond> jac;
  public:
    SymmetricBlockGaussSeidelPrecond (shared_ptr<const BaseSparseMatrix> mat, shared_ptr<Table<int>> blocktable)
      : jac(mat->CreateBlockJacobiPrecond(blocktable)) { }

    int VHeight() const override { return jac->VHeight(); }
    int VWidth() const override { return jac->VHeight(); }

    void Mult (const BaseVector & x, BaseVector & y) const override
    {
      y = 0;
      jac->GSSmooth (y, x);
      jac->GSSmoothBack (y, x);
    }
    
    AutoVector CreateRowVector() const override { return jac->CreateRowVector(); }
    AutoVector CreateColVector() const override { return jac->CreateColVector(); }
  };

  


  /**
     A block-Jacobi preconditioner.
     The blocks are specified by a table container
  */
  template <class TM, class TV_ROW, class TV_COL>
  class  NGS_DLL_HEADER BlockJacobiPrecond : virtual public BaseBlockJacobiPrecond,
                                         virtual public S_BaseMatrix<typename mat_traits<TM>::TSCAL>
  {
  protected:
    /// a reference to the matrix
    // const SparseMatrix<TM,TV_ROW,TV_COL> & mat;
    shared_ptr<const SparseMatrix<TM,TV_ROW,TV_COL>> mat;
    /// inverses of the small blocks
    Array<FlatMatrix<TM>> invdiag;
    /// the data for the inverses
    Array<TM> bigmem;

  public:
    // typedef typename mat_traits<TM>::TV_ROW TVX;
    typedef TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    ///
    BlockJacobiPrecond (shared_ptr<const SparseMatrix<TM,TV_ROW,TV_COL>> amat, 
			shared_ptr<Table<int>> ablocktable, bool cumulate_block_diags = false);
    ///
    virtual ~BlockJacobiPrecond ();

    size_t Height() const { return mat->Height(); }
    int VHeight() const override { return mat->Height(); }
    size_t Width() const { return mat->Width(); }
    int VWidth() const override { return mat->Width(); }

    AutoVector CreateRowVector() const override { return mat->CreateColVector(); }
    AutoVector CreateColVector() const override { return mat->CreateRowVector(); }

    BaseMatrix::OperatorInfo GetOperatorInfo () const override
    { return { string("BlockJacobi-")+typeid(TM).name(), this->Height(), this->Width() }; }
    
    ///
    void MultAdd (TSCAL s, const BaseVector & x, BaseVector & y) const override;

    ///
    void MultTransAdd (TSCAL s, const BaseVector & x, BaseVector & y) const override;


    ///
    void GSSmooth (BaseVector & x, const BaseVector & b,
                   int steps = 1) const override;

    void GSSmoothBack (BaseVector & x, const BaseVector & b,
                       int steps = 1) const override;
  
    void GSSmoothResiduum (BaseVector & x, const BaseVector & b,
                           BaseVector & res, int steps = 1) const  override
    {
      GSSmooth (x, b, steps);
      res = b - (*mat) * x;
    }

    ///
    virtual void GSSmoothNumbering (BaseVector & x, const BaseVector & b,
				    const Array<int> & numbering, 
				    int forward = 1) const
    {
      ;
    }

    Array<MemoryUsage> GetMemoryUsage () const override
    {
      int nels = 0;
      for (size_t i = 0; i < blocktable->Size(); i++)
	{
	  int bs = (*blocktable)[i].Size();
	  nels += bs*bs;
	}
      return { MemoryUsage ("BlockJac", nels*sizeof(TM), blocktable->Size()) };
    }

    const Array<FlatMatrix<TM>> & GetInverses() const { return invdiag; }
    const Array<TM> & MatrixData() const { return bigmem; } 
  };





  /* **************** SYMMETRIC ****************** */


  ///
  template <class TM, class TV>
  class BlockJacobiPrecondSymmetric : 
    virtual public BaseBlockJacobiPrecond,
    virtual public S_BaseMatrix<typename mat_traits<TM>::TSCAL>
  {
  protected:
    shared_ptr<const SparseMatrixSymmetric<TM,TV>> mat;

    // 
    enum { NBLOCKS = 20 };
    // DynamicMem<int> blockstart, blocksize, blockbw;
    // DynamicMem<TM> data[NBLOCKS];

    Array<int> blockstart, blocksize, blockbw;
    Array<TM> data[NBLOCKS];


    bool lowmem;
  public:
  
    typedef TV TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;
  
    ///
    NGS_DLL_HEADER BlockJacobiPrecondSymmetric (shared_ptr<const SparseMatrixSymmetric<TM,TV>> amat, 
                                                shared_ptr<Table<int>> ablocktable);
    
    /*
    BlockJacobiPrecondSymmetric (const SparseMatrixSymmetric<TM,TV> & amat, 
				 const FlatVector<TVX> & constraint,
				 Table<int> & ablocktable);
    */
    ///
    virtual ~BlockJacobiPrecondSymmetric ();
  
    FlatBandCholeskyFactors<TM> InvDiag (int i) const
    {
      return FlatBandCholeskyFactors<TM> (blocksize[i], 
					  blockbw[i], 
					  const_cast<TM*>(&data[i%NBLOCKS][blockstart[i]]));
    }

    void ComputeBlockFactor (FlatArray<int> block, int bw, FlatBandCholeskyFactors<TM> & inv) const;
  
    ///
    void MultAdd (TSCAL s, const BaseVector & x, BaseVector & y) const override;

    ///
    void MultTransAdd (TSCAL s, const BaseVector & x, BaseVector & y) const override;

    ///
    AutoVector CreateRowVector () const override { return mat->CreateColVector(); }
    AutoVector CreateColVector () const override { return mat->CreateRowVector(); }


    int Height() const { return mat->Height(); }
    int VHeight() const override { return mat->Height(); }
    int Width() const { return mat->Width(); }
    int VWidth() const override { return mat->Width(); }


    ///
    void GSSmooth (BaseVector & x, const BaseVector & b,
                   int steps = 1) const override;
  
    void GSSmoothPartial (BaseVector & x, const BaseVector & b,
                          BaseVector & y) const override;
  


    ///
    void GSSmoothResiduum (BaseVector & x, const BaseVector & b,
                           BaseVector & res, int steps = 1) const override;

  
    ///
    void GSSmoothBack (BaseVector & x, const BaseVector & b,
                       int steps = 1) const override;
 

    void GSSmoothBackPartial (BaseVector & x, const BaseVector & b,
                              BaseVector & y) const override;
 




    void SmoothBlock (int i, 
		      FlatVector<TVX> & x,
		      // const FlatVector<TVX> & b,
		      FlatVector<TVX> & y) const;
 

    ///
    virtual void GSSmoothNumbering (BaseVector & x, const BaseVector & b,
				    const Array<int> & numbering, 
				    int forward = 1) const
    {
      ;
    }
  
    /*  
	virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const
	{
	// Array<BandCholeskyFactors<TM>*> invdiag;


	int nels = 0;
	for (int i = 0; i < blocktable.Size(); i++)
	if (invdiag[i])
	{
	int s = invdiag[i]->Size();
	int bw = invdiag[i]->BandWidth();
	nels += s*bw - (bw * (bw-1)) / 2;
	nels += s;
	}
	mu.Append (new MemoryUsageStruct ("BlockJacSym Table", blocktable.NElements()*sizeof(int), 1));
	mu.Append (new MemoryUsageStruct ("BlockJacSym", nels*sizeof(TM), 2*blocktable.Size()));
	}
    */
  };


}

#endif
