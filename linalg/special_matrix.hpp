#ifndef FILE_NGS_SPECIALMATRIX
#define FILE_NGS_SPECIALMATRIX

/* ************************************************************************/
/* File:   special_matrix.hpp                                             */
/* Author: Joachim Schoeberl                                              */
/* Date:   14 Mar. 02                                                     */
/* ************************************************************************/

#include "basematrix.hpp"
#include "vvector.hpp"

namespace ngla
{

  // Convert RowMajor to ColMajor matrix (stored as vector)
  class TransposeVector : public BaseMatrix
  {
    int h, w; // result matrix
  public:
    typedef double TSCAL;
    
    TransposeVector (int ah, int aw);


    ostream & Print (ostream & ost) const override;
    
    VecFormat RowFormat () const override { return VVectorFormat<double> (h*w); }
    VecFormat ColFormat () const override { return VVectorFormat<double> (h*w); }

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



    VecFormat RowFormat () const override { return VecFormat(width); }
    VecFormat ColFormat () const override { return VecFormat(ind.Size()); }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual shared_ptr<BaseSparseMatrix> CreateSparseMatrix() const override;
    FlatArray<size_t> GetIndices() const { return ind; }
  };


  class Embedding : public BaseMatrix
  {
    size_t height;
    IntRange range;
    bool is_complex;
  public:
    Embedding (size_t aheight, IntRange arange, bool ais_complex = false)
      : height(aheight), range(arange), is_complex(ais_complex) { ; }



    // complex if requested, otherwise the scalar stays open
    VecFormat RowFormat () const override
    { return is_complex ? VecFormat(range.Size(), Complex(0)) : VecFormat(range.Size()); }
    VecFormat ColFormat () const override
    { return is_complex ? VecFormat(height, Complex(0)) : VecFormat(height); }

    auto GetRange() const { return range; }
    
    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual shared_ptr<BaseSparseMatrix> CreateSparseMatrix() const override;    

    NGS_DLL_HEADER shared_ptr<BaseMatrix> CreateDeviceMatrix () const override;
  };


  
  class EmbeddedMatrix : public BaseMatrix
  {
    size_t height;
    IntRange range;
    shared_ptr<BaseMatrix> mat;
  public:
    EmbeddedMatrix (size_t aheight, IntRange arange, shared_ptr<BaseMatrix> amat)
      : height(aheight), range(arange), mat(amat) { ; }



    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;

    VecFormat RowFormat () const override { return mat->RowFormat(); }
    VecFormat ColFormat () const override
    { return VecFormat::Merge (VecFormat(height), mat->ColFormat().ValueAxes()); }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    shared_ptr<BaseMatrix> GetMatrix() const { return mat; }
    IntRange GetRange() const { return range; }

    virtual shared_ptr<BaseSparseMatrix> CreateSparseMatrix() const override;
    
    
    NGS_DLL_HEADER shared_ptr<BaseMatrix> CreateDeviceMatrix() const override;
  };


  class EmbeddingTranspose : public BaseMatrix
  {
    size_t width;
    IntRange range;
    bool is_complex;
  public:
    EmbeddingTranspose (size_t awidth, IntRange arange, bool ais_complex = false)
      : width(awidth), range(arange), is_complex(ais_complex) { ; }

    

    VecFormat RowFormat () const override
    { return is_complex ? VecFormat(width, Complex(0)) : VecFormat(width); }
    VecFormat ColFormat () const override
    { return is_complex ? VecFormat(range.Size(), Complex(0)) : VecFormat(range.Size()); }

    auto GetRange() const { return range; }
    
    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual shared_ptr<BaseSparseMatrix> CreateSparseMatrix() const override;    
    NGS_DLL_HEADER shared_ptr<BaseMatrix> CreateDeviceMatrix () const override;
  };


  class EmbeddedTransposeMatrix : public BaseMatrix
  {
    size_t width;
    IntRange range;
    shared_ptr<BaseMatrix> mat;
  public:
    EmbeddedTransposeMatrix (size_t awidth, IntRange arange, shared_ptr<BaseMatrix> amat)
      : width(awidth), range(arange), mat(amat) { ; }


    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    

    VecFormat RowFormat () const override
    { return VecFormat::Merge (VecFormat(width), mat->RowFormat().ValueAxes()); }
    VecFormat ColFormat () const override { return mat->ColFormat(); }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    shared_ptr<BaseMatrix> GetMatrix() const { return mat; }
    virtual shared_ptr<BaseSparseMatrix> CreateSparseMatrix() const override;
    IntRange GetRange() const { return range; }

    NGS_DLL_HEADER shared_ptr<BaseMatrix> CreateDeviceMatrix() const override;
  };



  
  
  template <class TVR, class TVC>
  class Real2ComplexMatrix : public BaseMatrix
  {
    shared_ptr<BaseMatrix> realmatrix;
    VVector<TVR> hx, hy;
  public:
    NGS_DLL_HEADER Real2ComplexMatrix (shared_ptr<BaseMatrix> arealmatrix = nullptr);
    void SetMatrix (shared_ptr<BaseMatrix> arealmatrix);
    const BaseMatrix & GetMatrix () const { return *realmatrix; }
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override;

    VecFormat RowFormat () const override { return VVectorFormat<TVC> (realmatrix->Width()); }
    VecFormat ColFormat () const override { return VVectorFormat<TVC> (realmatrix->Width()); }
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
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
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

    VecFormat RowFormat () const override { return VecFormat(); }
    VecFormat ColFormat () const override { return VecFormat(); }
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

    VecFormat RowFormat () const override;
    VecFormat ColFormat () const override;
  };


  
  class BaseMatrixFromVector : public BaseMatrix
  {
    shared_ptr<BaseVector> vec;

  public:
    BaseMatrixFromVector (shared_ptr<BaseVector> avec);

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;


    // missing parallel: 1 dof for all
    VecFormat RowFormat () const override { return VVectorFormat<double> (1); }
    VecFormat ColFormat () const override { return vec->GetFormat(); }
  };


  class BaseMatrixFromMultiVector : public BaseMatrix
  {
    shared_ptr<MultiVector> vec;

  public:
    BaseMatrixFromMultiVector (shared_ptr<MultiVector> avec);

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;


    // missing parallel: 1 dof for all
    VecFormat RowFormat () const override { return VVectorFormat<double> (vec->Size()); }
    VecFormat ColFormat () const override { return vec->RefVec()->GetFormat(); }
  };


  template <typename T = double>
  class BaseMatrixFromMatrix : public BaseMatrix
  {
    Matrix<T> mat;

  public:
    NGS_DLL_HEADER BaseMatrixFromMatrix (Matrix<T> amat);

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override;

    virtual size_t NZE() const override { return mat.Height()*mat.Width(); }
    VecFormat RowFormat () const override { return VVectorFormat<T> (mat.Width()); }
    VecFormat ColFormat () const override { return VVectorFormat<T> (mat.Height()); }
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
    
    BaseVector & AsVector() override;
    const BaseVector & AsVector() const override;
    void SetZero() override;

    ostream & Print (ostream & ost) const override { return mat->Print(ost); }
    Array<MemoryUsage> GetMemoryUsage () const override { return mat->GetMemoryUsage(); }
    size_t NZE () const override { return mat->NZE(); }

    void Update() override { mat->Update(); }
    shared_ptr<BaseMatrix> CreateMatrix () const override { return mat->CreateMatrix(); }
    VecFormat RowFormat () const override { return mat->RowFormat(); }
    VecFormat ColFormat () const override { return mat->ColFormat(); }
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
    
    void DoArchive (Archive & ar) override
    { mat->DoArchive(ar); }
  };
  



  
}


#endif
