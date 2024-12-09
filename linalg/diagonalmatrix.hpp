#ifndef FILE_NGS_DIAGONALMATRIX
#define FILE_NGS_DIAGONALMATRIX



#include "basematrix.hpp"
#include "vvector.hpp"


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

    virtual shared_ptr<BaseSparseMatrix> CreateSparseMatrix() const override;
    
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
  class NGS_DLL_HEADER DiagonalMatrix : public BaseMatrix
  {
    shared_ptr<VVector<TM>> diag;
  public:
    // typedef typename mat_traits<TM>::TV_ROW TV_ROW;
    // typedef typename mat_traits<TM>::TV_COL TV_COL;
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    
    DiagonalMatrix(size_t h)
      : diag(make_shared<VVector<TM>>(h)) { }
    DiagonalMatrix(const VVector<TM> & diag_);
    DiagonalMatrix(shared_ptr<VVector<TM>> diag_);
    virtual ~DiagonalMatrix();
    
    bool IsComplex() const override { return false; } 
    TM & operator() (size_t i) { return (*diag)(i); }
    const TM & operator() (size_t i) const { return (*diag)(i); }
    int VHeight() const override { return diag->Size(); }
    int VWidth() const override { return diag->Size(); }

    BaseVector & AsVector() override { return *diag; }
    const BaseVector & AsVector() const override { return *diag; }
    ostream & Print (ostream & ost) const override;

    virtual shared_ptr<BaseSparseMatrix> CreateSparseMatrix() const override;
    
    AutoVector CreateRowVector () const override;
    AutoVector CreateColVector () const override;

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;    
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;
  };


  template <typename TM=double>  
  class BlockDiagonalMatrix : public BaseMatrix
  {
    Tensor<3,TM> blockdiag;  
    int blocks, dimy, dimx;
  public:
    // typedef double TSCAL;
    
    BlockDiagonalMatrix(Tensor<3,TM> _blockdiag);
    bool IsComplex() const override { return ngbla::IsComplex<TM>(); } 

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
    Matrix<bool> nonzero;
    Table<int> sparse, sparseT; // sparse non-zero pattern
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

    void Mult (const BaseVector & x, BaseVector & y) const override;    
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    void MultTrans (const BaseVector & x, BaseVector & y) const override;    
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    // shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;

    FlatTensor<3> GetBlockDiag () const { return blockdiag; }
    FlatMatrix<bool> GetNonZeroPattern() const { return nonzero; }
  };

}



#endif

