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
    const BitArray & bits;
    bool keep_values;
  public:
    // projector on true / false bits
    Projector (const BitArray & abits, bool akeep_values = true)
      : bits(abits), keep_values(akeep_values) { ; }

    virtual bool IsComplex() const { return false; } 

    virtual int VHeight() const { return bits.Size(); }
    virtual int VWidth() const { return bits.Size(); }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;    
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


}


#endif
