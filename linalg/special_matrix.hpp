#ifndef FILE_NGS_SPECIALMATRIX
#define FILE_NGS_SPECIALMATRIX

/* ************************************************************************/
/* File:   special_matrix.hpp                                             */
/* Author: Joachim Schoeberl                                              */
/* Date:   14 Mar. 02                                                     */
/* ************************************************************************/

namespace ngla
{


  class NGS_DLL_HEADER Projector : public BaseMatrix
  {
    shared_ptr<BitArray> bits;
    bool keep_values;
  public:
    // projector on true / false bits
    Projector (shared_ptr<BitArray> abits, bool akeep_values = true)
      : bits(abits), keep_values(akeep_values) { ; }
    
    virtual bool IsComplex() const override { return false; } 

    virtual int VHeight() const override { return bits->Size(); }
    virtual int VWidth() const override { return bits->Size(); }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;    
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;    
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void Project (BaseVector & x) const;
    virtual void SetValues (BaseVector & x, double val) const;

    bool KeepValues() const { return keep_values; }
    shared_ptr<BitArray> Mask() const { return bits; }
    
    AutoVector CreateRowVector() const override
    { throw Exception("CreateRowVector not implemented for Projector!"); }
    AutoVector CreateColVector() const override
    { throw Exception("CreateColVector not implemented for Projector!"); }

    AutoVector Evaluate(BaseVector & v) const override
    {
      auto res = v.CreateVector();
      Mult (v, res);
      return res;
    }
    
  };


  template <typename TM=double>
  class DiagonalMatrix : public BaseMatrix
  {
    shared_ptr<VVector<TM>> diag;
  public:
    // typedef typename mat_traits<TM>::TV_ROW TV_ROW;
    // typedef typename mat_traits<TM>::TV_COL TV_COL;
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    
    DiagonalMatrix(size_t h)
      : diag(make_shared<VVector<TM>>(h)) { }
    DiagonalMatrix(const VVector<TM> & diag_)
      : diag(make_shared<VVector<TM>>(diag_)) { } 
    DiagonalMatrix(shared_ptr<VVector<TM>> diag_)
      : diag(diag_) { } 
    
    bool IsComplex() const override { return false; } 
    TM & operator() (size_t i) { return (*diag)(i); }
    const TM & operator() (size_t i) const { return (*diag)(i); }
    int VHeight() const override { return diag->Size(); }
    int VWidth() const override { return diag->Size(); }

    BaseVector & AsVector() override { return *diag; }
    const BaseVector & AsVector() const override { return *diag; }
    ostream & Print (ostream & ost) const override;
    
    AutoVector CreateRowVector () const override;
    AutoVector CreateColVector () const override;

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;    
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;
  };



  class BlockDiagonalMatrix : public BaseMatrix
  {
    Tensor<3> blockdiag;  
    int blocks, dimy, dimx;
  public:
    typedef double TSCAL;
    
    BlockDiagonalMatrix(Tensor<3> _blockdiag);
    bool IsComplex() const override { return false; } 

    int VHeight() const override { return blocks*dimy; }
    int VWidth() const override { return blocks*dimx; }

    ostream & Print (ostream & ost) const override;
    
    AutoVector CreateRowVector () const override;
    AutoVector CreateColVector () const override;

    void Mult (const BaseVector & x, BaseVector & y) const override;    
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;    
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;
  };

  // blocks is inner-most dimension of tensor and vectors
  class BlockDiagonalMatrixSoA : public BaseMatrix
  {
    Tensor<3> blockdiag;  
    int blocks, dimy, dimx;
    Matrix<double> nonzero;
  public:
    typedef double TSCAL;
    
    BlockDiagonalMatrixSoA(Tensor<3> _blockdiag);
    bool IsComplex() const override { return false; } 

    int VHeight() const override { return blocks*dimy; }
    int VWidth() const override { return blocks*dimx; }

    ostream & Print (ostream & ost) const override;
    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    
    AutoVector CreateRowVector () const override;
    AutoVector CreateColVector () const override;

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;    
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    // shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;

    FlatTensor<3> GetBlockDiag () const { return blockdiag; }
    FlatMatrix<double> GetNonZeroPattern() const { return nonzero; }
  };


  

  // Convert RowMajor to ColMajor matrix (stored as vector)
  class TransposeVector : public BaseMatrix
  {
    int h, w; // result matrix
  public:
    typedef double TSCAL;
    
    TransposeVector (int ah, int aw);
    bool IsComplex() const override { return false; } 

    int VHeight() const override { return h*w; }
    int VWidth() const override { return h*w; }

    ostream & Print (ostream & ost) const override;
    
    AutoVector CreateRowVector () const override;
    AutoVector CreateColVector () const override;

    void Mult (const BaseVector & x, BaseVector & y) const override;    
    void MultTrans (const BaseVector & x, BaseVector & y) const override;    
  };

  
  
  class PermutationMatrix : public BaseMatrix
  {
    size_t width;
    Array<size_t> ind;
  public:
    PermutationMatrix (size_t awidth, Array<size_t> aind)
      : width(awidth), ind(aind) { ; } 

    virtual bool IsComplex() const override { return false; } 

    virtual int VHeight() const override { return ind.Size(); }
    virtual int VWidth() const override { return width; }

    virtual AutoVector CreateRowVector () const override
    {
      return CreateBaseVector(width, false, 1);
    }
    
    virtual AutoVector CreateColVector () const override
    {
      return CreateBaseVector(ind.Size(), false, 1);
    }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
  };


  class Embedding : public BaseMatrix
  {
    size_t height;
    IntRange range;
    bool is_complex;
  public:
    Embedding (size_t aheight, IntRange arange, bool ais_complex = false)
      : height(aheight), range(arange), is_complex(ais_complex) { ; }

    virtual bool IsComplex() const override { return is_complex; } 

    virtual int VHeight() const override { return height; }
    virtual int VWidth() const override { return range.Size(); }

    virtual AutoVector CreateRowVector () const override
    {
      return CreateBaseVector(range.Size(), is_complex, 1);
    }
    
    virtual AutoVector CreateColVector () const override
    {
      return CreateBaseVector(height, is_complex, 1);
    }

    auto GetRange() const { return range; }
    
    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
  };


  
  class EmbeddedMatrix : public BaseMatrix
  {
    size_t height;
    IntRange range;
    shared_ptr<BaseMatrix> mat;
  public:
    EmbeddedMatrix (size_t aheight, IntRange arange, shared_ptr<BaseMatrix> amat)
      : height(aheight), range(arange), mat(amat) { ; }

    virtual bool IsComplex() const override { return mat->IsComplex(); } 

    virtual int VHeight() const override { return height; }
    virtual int VWidth() const override { return mat->VWidth(); }

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;

    virtual AutoVector CreateRowVector () const override
    {
      return mat->CreateRowVector();      
    }
    
    virtual AutoVector CreateColVector () const override
    {
      return CreateBaseVector(height, IsComplex(), 1);      
    }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    shared_ptr<BaseMatrix> GetMatrix() const { return mat; }
    IntRange GetRange() const { return range; }

    /*
    virtual shared_ptr<BaseMatrix> CreateDeviceMatrix() const override
    {
      return make_shared<EmbeddedMatrix>(height, range, mat->CreateDeviceMatrix());
    }
    */
  };


  class EmbeddingTranspose : public BaseMatrix
  {
    size_t width;
    IntRange range;
    bool is_complex;
  public:
    EmbeddingTranspose (size_t awidth, IntRange arange, bool ais_complex = false)
      : width(awidth), range(arange), is_complex(ais_complex) { ; }

    virtual bool IsComplex() const override { return is_complex; } 
    
    virtual int VHeight() const override { return range.Size(); }
    virtual int VWidth() const override { return width; }

    virtual AutoVector CreateRowVector () const override
    {
      return CreateBaseVector(width, is_complex, 1);
    }
    
    virtual AutoVector CreateColVector () const override
    {
      return CreateBaseVector(range.Size(), is_complex, 1);
    }

    auto GetRange() const { return range; }
    
    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
  };


  class EmbeddedTransposeMatrix : public BaseMatrix
  {
    size_t width;
    IntRange range;
    shared_ptr<BaseMatrix> mat;
  public:
    EmbeddedTransposeMatrix (size_t awidth, IntRange arange, shared_ptr<BaseMatrix> amat)
      : width(awidth), range(arange), mat(amat) { ; }

    virtual bool IsComplex() const override { return mat->IsComplex(); } 

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    
    virtual int VHeight() const override { return mat->Height(); }
    virtual int VWidth() const override { return width; }

    virtual AutoVector CreateRowVector () const override
    {
      return CreateBaseVector(width, mat->IsComplex(), 1);
    }
    
    virtual AutoVector CreateColVector () const override
    {
      return mat->CreateColVector();
    }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    shared_ptr<BaseMatrix> GetMatrix() const { return mat; }
    IntRange GetRange() const { return range; }

    /*
    virtual shared_ptr<BaseMatrix> CreateDeviceMatrix() const override
    {
      return make_shared<EmbeddedTransposeMatrix>(width, range, mat->CreateDeviceMatrix());
    }
    */
  };



  
  
  template <class TVR, class TVC>
  class Real2ComplexMatrix : public BaseMatrix
  {
    shared_ptr<BaseMatrix> realmatrix;
    VVector<TVR> hx, hy;
  public:
    NGS_DLL_HEADER Real2ComplexMatrix (shared_ptr<BaseMatrix> arealmatrix = nullptr);
    bool IsComplex() const override { return true; }     
    void SetMatrix (shared_ptr<BaseMatrix> arealmatrix);
    const BaseMatrix & GetMatrix () const { return *realmatrix; }
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override;

    AutoVector CreateRowVector() const override;
    AutoVector CreateColVector() const override;
  };




  ////////////////////////////////////////////////////////////////////////////////
  // added 08/19/2003
  template <class TVR>
  class NGS_DLL_HEADER Sym2NonSymMatrix : public BaseMatrix
  {
    const BaseMatrix * base;
    VVector<TVR> hx, hy;
  public:
    Sym2NonSymMatrix (const BaseMatrix * abasematrix = 0);
    void SetMatrix (const BaseMatrix * abasematrix);
    const BaseMatrix & GetMatrix () const { return *base; }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    //  virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const;
  };



  ////////////////////////////////////////////////////////////////////////////////
  // added 09/02/2003
  template <class TVSMALL, class TVBIG>
  class NGS_DLL_HEADER Small2BigNonSymMatrix : public BaseMatrix
  {
    const BaseMatrix * base;
    VVector<TVSMALL> hx1, hx2, hy1, hy2;
  public:
    Small2BigNonSymMatrix (const BaseMatrix * abasematrix = 0);
    void SetMatrix (const BaseMatrix * abasematrix);
    bool IsComplex() const override { return base->IsComplex(); }
    const BaseMatrix & GetMatrix () const { return *base; }
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    //  virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const;
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    AutoVector CreateRowVector() const override
    { throw Exception("CreateRowVector not implemented for Small2BigNonSymMatrix!"); }
    AutoVector CreateColVector() const override
    { throw Exception("CreateColVector not implemented for Small2BigNonSymMatrix!"); }
  };


  class BlockMatrix : public BaseMatrix
  {
    Array<Array<shared_ptr<BaseMatrix>>> mats;
    size_t h, w;

    // one matrix per row/col that can be used to create vectors etc.
    Array<shared_ptr<BaseMatrix>> row_reps;
    Array<shared_ptr<BaseMatrix>> col_reps;

  public:
    BlockMatrix (const Array<Array<shared_ptr<BaseMatrix>>> & amats);

    bool IsComplex() const override { return row_reps[0]->IsComplex(); }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    const shared_ptr<BaseMatrix> & operator()(size_t i, size_t j)
    {
      if (i >= h) throw Exception("Tried to access BlockMatrix row that is out of range");
      if (j >= w) throw Exception("Tried to access BlockMatrix col that is out of range");
      return mats[i][j];
    }

    size_t BlockRows() const { return h; }
    size_t BlockCols() const { return w; }

    virtual int VHeight() const override { throw Exception("VHeight does not make sense for BlockMatrix");}
    virtual int VWidth() const override { throw Exception("VWidth does not make sense for BlockMatrix");}

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
  };


  
  class BaseMatrixFromVector : public BaseMatrix
  {
    shared_ptr<BaseVector> vec;

  public:
    BaseMatrixFromVector (shared_ptr<BaseVector> avec);

    bool IsComplex() const override { return vec->IsComplex(); }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual int VHeight() const override { return vec->Size(); }
    virtual int VWidth() const override { return 1; }

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
  };


  class BaseMatrixFromMultiVector : public BaseMatrix
  {
    shared_ptr<MultiVector> vec;

  public:
    BaseMatrixFromMultiVector (shared_ptr<MultiVector> avec);

    bool IsComplex() const override { return vec->IsComplex(); }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual int VHeight() const override { return vec->RefVec()->Size(); }
    virtual int VWidth() const override { return vec->Size(); }

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
  };


  template <typename T = double>
  class BaseMatrixFromMatrix : public BaseMatrix
  {
    Matrix<T> mat;

  public:
    NGS_DLL_HEADER BaseMatrixFromMatrix (Matrix<T> amat);

    bool IsComplex() const override { return typeid(T) == typeid(Complex); }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override;

    virtual int VHeight() const override { return mat.Height(); }
    virtual int VWidth() const override { return mat.Width(); }
    virtual size_t NZE() const override { return mat.Height()*mat.Width(); }
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
  };




  class LoggingMatrix : public BaseMatrix
  {
    shared_ptr<BaseMatrix> mat;
    string label;
    unique_ptr<ostream> out;
    optional<NgMPI_Comm> comm;
  public:
    LoggingMatrix (shared_ptr<BaseMatrix> amat, string alabel, string filename,
                   optional<NgMPI_Comm> acomm);
    ~LoggingMatrix ();
    int VHeight() const override { return mat->VHeight(); }
    int VWidth() const override { return mat->VWidth(); }
    bool IsComplex() const override { return mat->IsComplex(); }
    
    BaseVector & AsVector() override;
    const BaseVector & AsVector() const override;
    void SetZero() override;

    ostream & Print (ostream & ost) const override { return mat->Print(ost); }
    Array<MemoryUsage> GetMemoryUsage () const override { return mat->GetMemoryUsage(); }
    size_t NZE () const override { return mat->NZE(); }

    void Update() override { mat->Update(); }
    shared_ptr<BaseMatrix> CreateMatrix () const override { return mat->CreateMatrix(); }
    AutoVector CreateRowVector () const override;
    AutoVector CreateColVector () const override;

    void Mult (const BaseVector & x, BaseVector & y) const override;
    void MultTrans (const BaseVector & x, BaseVector & y) const override;
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override;
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override;
    void MultConjTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override;
    void MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const override;
    
    void MultAdd1 (double s, const BaseVector & x, BaseVector & y,
                   const BitArray * ainner = NULL,
                   const Array<int> * acluster = NULL) const override
    { mat->MultAdd1 (s, x, y, ainner, acluster); }
    
    void MultAdd2 (double s, const BaseVector & x, BaseVector & y,
                   const BitArray * ainner = NULL,
                   const Array<int> * acluster = NULL) const override
    { mat->MultAdd2 (s, x, y, ainner, acluster); }

    shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override
    { return mat->InverseMatrix(subset); }
    
    shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<const Array<int>> clusters) const override
    { return mat->InverseMatrix(clusters); }
    
    INVERSETYPE SetInverseType ( INVERSETYPE ainversetype ) const override
    { return mat->SetInverseType(ainversetype); }
    INVERSETYPE SetInverseType ( string ainversetype ) const override
    { return mat->SetInverseType(ainversetype); }
    INVERSETYPE  GetInverseType () const override
    { return mat->GetInverseType(); }
    void DoArchive (Archive & ar) override
    { mat->DoArchive(ar); }
  };
  



  
}


#endif
