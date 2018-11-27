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

    virtual bool IsComplex() const { return false; } 

    virtual int VHeight() const { return bits->Size(); }
    virtual int VWidth() const { return bits->Size(); }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;    
    virtual void Project (BaseVector & x) const;    
  };

  class Embedding : public BaseMatrix
  {
    size_t height;
    IntRange range;
  public:
    Embedding (size_t aheight, IntRange arange)
      : height(aheight), range(arange) { ; }

    virtual bool IsComplex() const override { return false; } 

    virtual int VHeight() const override { return height; }
    virtual int VWidth() const override { return range.Size(); }

    virtual AutoVector CreateRowVector () const override
    {
      return CreateBaseVector(range.Size(), false, 1);
    }
    
    virtual AutoVector CreateColVector () const override
    {
      return CreateBaseVector(height, false, 1);
    }
    
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
  };


  
  template <class TVR, class TVC>
  class Real2ComplexMatrix : public BaseMatrix
  {
    const BaseMatrix * realmatrix;
    VVector<TVR> hx, hy;
  public:
    NGS_DLL_HEADER Real2ComplexMatrix (const BaseMatrix * arealmatrix = 0);
    virtual bool IsComplex() const { return true; }     
    void SetMatrix (const BaseMatrix * arealmatrix);
    const BaseMatrix & GetMatrix () const { return *realmatrix; }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const;
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
    virtual bool IsComplex() const { return base->IsComplex(); }     
    const BaseMatrix & GetMatrix () const { return *base; }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    //  virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;
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
