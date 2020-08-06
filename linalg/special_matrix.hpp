#ifndef FILE_NGS_SPECIALMATRIX
#define FILE_NGS_SPECIALMATRIX

/* ************************************************************************/
/* File:   special_matrix.hpp                                             */
/* Author: Joachim Schoeberl                                              */
/* Date:   14 Mar. 02                                                     */
/* ************************************************************************/

namespace ngla
{


  class Projector : public BaseMatrix
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

    AutoVector CreateRowVector() const override
    { throw Exception("CreateRowVector not implemented for Projector!"); }
    AutoVector CreateColVector() const override
    { throw Exception("CreateColVector not implemented for Projector!"); }

  };


  template <typename TM=double>
  class DiagonalMatrix : public BaseMatrix
  {
    VVector<TM> diag;
  public:
    // typedef typename mat_traits<TM>::TV_ROW TV_ROW;
    // typedef typename mat_traits<TM>::TV_COL TV_COL;
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    
    DiagonalMatrix(size_t h)
      : diag(h) { }
    DiagonalMatrix(const VVector<TM> & diag_)
      : diag(diag_) { } 

    bool IsComplex() const override { return false; } 
    TM & operator() (size_t i) { return diag(i); }
    const TM & operator() (size_t i) const { return diag(i); }
    int VHeight() const override { return diag.Size(); }
    int VWidth() const override { return diag.Size(); }

    BaseVector & AsVector() override { return diag; }
    ostream & Print (ostream & ost) const override;
    
    AutoVector CreateRowVector () const override;
    AutoVector CreateColVector () const override;

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;    
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;
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

    AutoVector CreateRowVector() const override
    { return realmatrix->CreateRowVector(); }
    AutoVector CreateColVector() const override
    { throw Exception("CreateColVector not implemented for Real2ComplexMatrix!"); }
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
}


#endif
