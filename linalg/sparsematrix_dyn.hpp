#ifndef FILE_NGS_SPARSEMATRIX_DYN
#define FILE_NGS_SPARSEMATRIX_DYN

/**************************************************************************/
/* File:   sparsematrix_dyn.hpp                                           */
/* Author: Joachim Schoeberl                                              */
/* Date:   July 2019                                                      */
/**************************************************************************/


namespace ngla
{

/// A general, sparse matrix
  template<class TSCAL>
  class  NGS_DLL_HEADER SparseMatrixDynamic : public BaseSparseMatrix, 
                                              public S_BaseMatrix<TSCAL>
  {
  protected:
    size_t bh, bw, bs;
    Array<TSCAL> data;
    // NumaDistributedArray<TM> data;
    // VFlatVector<typename mat_traits<TM>::TSCAL> asvec;
    TSCAL nul;
    
  public:
    template <typename TM>
      SparseMatrixDynamic (const SparseMatrixTM<TM> & mat)
      : BaseSparseMatrix (mat, false)
    {
      width = mat.Width();
      bh = mat_traits<TM>::HEIGHT;
      bw = mat_traits<TM>::WIDTH;
      bs = bh*bw;
      nze = mat.NZE();
      data.SetSize(nze*bs);
      auto matvec = mat.AsVector().template FV<TM>();
      for (size_t i = 0; i < nze; i++)
        {
          FlatMatrix<TSCAL> fm(bh, bw, &data[i*bs]);
          fm = matvec(i);
        }
    }

    virtual int VHeight() const override { return size; }
    virtual int VWidth() const override { return width; }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
  };




  /*
    SparseMatrixTM (int as, int max_elsperrow)
      : BaseSparseMatrix (as, max_elsperrow),
	data(nze), nul(TSCAL(0))
    {
      ; 
    }

    SparseMatrixTM (const Array<int> & elsperrow, int awidth)
      : BaseSparseMatrix (elsperrow, awidth), 
	data(nze), nul(TSCAL(0))
    {
      ; 
    }

    SparseMatrixTM (int size, int width, const Table<int> & rowelements, 
		    const Table<int> & colelements, bool symmetric)
      : BaseSparseMatrix (size, width, rowelements, colelements, symmetric), 
	data(nze), nul(TSCAL(0))
    { 
      ; 
    }

    SparseMatrixTM (const MatrixGraph & agraph, bool stealgraph)
      : BaseSparseMatrix (agraph, stealgraph), 
	data(nze), nul(TSCAL(0))
    { 
      FindSameNZE();
    }

    SparseMatrixTM (const SparseMatrixTM & amat)
    : BaseSparseMatrix (amat), 
      data(nze), nul(TSCAL(0))
    { 
      AsVector() = amat.AsVector(); 
    }

    static shared_ptr<SparseMatrixTM> CreateFromCOO (FlatArray<int> i, FlatArray<int> j,
                                                     FlatArray<TSCAL> val, size_t h, size_t w);
      
    virtual ~SparseMatrixTM ();

    int Height() const { return size; }
    int Width() const { return width; }
    virtual int VHeight() const { return size; }
    virtual int VWidth() const { return width; }

    TM & operator[] (int i)  { return data[i]; }
    const TM & operator[] (int i) const { return data[i]; }

    TM & operator() (int row, int col)
    {
      return data[CreatePosition(row, col)];
    }

    const TM & operator() (int row, int col) const
    {
      size_t pos = GetPositionTest (row,col);
      if (pos != numeric_limits<size_t>::max())
	return data[pos];
      else
	return nul;
    }

    void PrefetchRow (int rownr) const;
    
    FlatVector<TM> GetRowValues(int i) const
      // { return FlatVector<TM> (firsti[i+1]-firsti[i], &data[firsti[i]]); }
    { return FlatVector<TM> (firsti[i+1]-firsti[i], data+firsti[i]); }

    static bool IsRegularIndex (int index) { return index >= 0; }
    virtual void AddElementMatrix(FlatArray<int> dnums1, 
                                  FlatArray<int> dnums2, 
                                  BareSliceMatrix<TSCAL> elmat,
                                  bool use_atomic = false);

    virtual void AddElementMatrixSymmetric(FlatArray<int> dnums,
                                           BareSliceMatrix<TSCAL> elmat,
                                           bool use_atomic = false);
    
    virtual BaseVector & AsVector() 
    {
      // asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)&data[0]);
      asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)data.Addr(0));
      return asvec; 
    }

    virtual const BaseVector & AsVector() const
    { 
      const_cast<VFlatVector<TSCAL>&> (asvec).
	AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)&data[0]);
      return asvec; 
    }

    virtual void SetZero();


    ///
    virtual ostream & Print (ostream & ost) const;

    ///
    virtual Array<MemoryUsage> GetMemoryUsage () const;    
  };
  



  template<class TM, class TV_ROW, class TV_COL>
  class NGS_DLL_HEADER SparseMatrix : public SparseMatrixTM<TM>
  {
  public:
    using SparseMatrixTM<TM>::firsti;
    using SparseMatrixTM<TM>::colnr;
    using SparseMatrixTM<TM>::data;
    using SparseMatrixTM<TM>::balance;


    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV_ROW TVX;
    typedef TV_COL TVY;

    ///
    SparseMatrix (int as, int max_elsperrow)
      : SparseMatrixTM<TM> (as, max_elsperrow) { ; }

    SparseMatrix (const Array<int> & aelsperrow)
      : SparseMatrixTM<TM> (aelsperrow, aelsperrow.Size()) { ; }

    SparseMatrix (const Array<int> & aelsperrow, int awidth)
      : SparseMatrixTM<TM> (aelsperrow, awidth) { ; }

    SparseMatrix (int height, int width, const Table<int> & rowelements, 
		  const Table<int> & colelements, bool symmetric)
      : SparseMatrixTM<TM> (height, width, rowelements, colelements, symmetric) { ; }

    SparseMatrix (const MatrixGraph & agraph, bool stealgraph);
    // : SparseMatrixTM<TM> (agraph, stealgraph) { ; }

    SparseMatrix (const SparseMatrix & amat)
      : SparseMatrixTM<TM> (amat) { ; }

    SparseMatrix (const SparseMatrixTM<TM> & amat)
      : SparseMatrixTM<TM> (amat) { ; }

    virtual shared_ptr<BaseMatrix> CreateMatrix () const override;
    // virtual BaseMatrix * CreateMatrix (const Array<int> & elsperrow) const;
    ///
    virtual AutoVector CreateVector () const override;
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
    
    virtual shared_ptr<BaseJacobiPrecond>
      CreateJacobiPrecond (shared_ptr<BitArray> inner) const override
    { 
      return make_shared<JacobiPrecond<TM,TV_ROW,TV_COL>> (*this, inner);
    }
    
    virtual shared_ptr<BaseBlockJacobiPrecond>
      CreateBlockJacobiPrecond (shared_ptr<Table<int>> blocks,
                                const BaseVector * constraint = 0,
                                bool parallel = 1,
                                shared_ptr<BitArray> freedofs = NULL) const override
    { 
      return make_shared<BlockJacobiPrecond<TM,TV_ROW,TV_COL>> (*this, blocks );
    }

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;
    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<const Array<int>> clusters) const override;

    virtual shared_ptr<BaseSparseMatrix> Restrict (const SparseMatrixTM<double> & prol,
					 shared_ptr<BaseSparseMatrix> cmat = nullptr) const override;
  
    ///
    inline TVY RowTimesVector (int row, const FlatVector<TVX> vec) const
    {
      typedef typename mat_traits<TVY>::TSCAL TTSCAL;
      TVY sum = TTSCAL(0);
      for (size_t j = firsti[row]; j < firsti[row+1]; j++)
	sum += data[j] * vec(colnr[j]);
      return sum;
    }

    ///
    void AddRowTransToVector (int row, TVY el, FlatVector<TVX> vec) const
    {
      size_t first = firsti[row];
      size_t last = firsti[row+1];

      const int * colpi = colnr.Addr(0);
      const TM * datap = data.Addr(0);

      for (size_t j = first; j < last; j++)
        vec[colpi[j]] += Trans(datap[j]) * el; 
    }
    
    ///
    void AddRowConjTransToVector (int row, TVY el, FlatVector<TVX> vec) const
    {
      size_t first = firsti[row];
      size_t last = firsti[row+1];

      const int * colpi = colnr.Addr(0);
      const TM * datap = data.Addr(0);

      for (size_t j = first; j < last; j++)
        vec[colpi[j]] += Conj(Trans(datap[j])) * el; 
    }


    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultConjTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd1 (double s, const BaseVector & x, BaseVector & y,
			   const BitArray * ainner = NULL,
			   const Array<int> * acluster = NULL) const override;
    
    virtual void DoArchive (Archive & ar) override;
  };
*/

}
#endif
  
