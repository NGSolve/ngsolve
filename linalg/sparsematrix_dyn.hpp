#ifndef FILE_NGS_SPARSEMATRIX_DYN
#define FILE_NGS_SPARSEMATRIX_DYN

/**************************************************************************/
/* File:   sparsematrix_dyn.hpp                                           */
/* Author: Joachim Schoeberl                                              */
/* Date:   July 2019                                                      */
/**************************************************************************/


#include "sparsematrix.hpp"

namespace ngla
{

/// Sparse matrix with dynamic block size (still experimental)
  template<class TSCAL>
  class  NGS_DLL_HEADER SparseMatrixDynamic : public BaseSparseMatrix
  // public S_BaseMatrix<TSCAL>
  {
  protected:
    size_t bh, bw, bs;
    Array<TSCAL> data;
    TSCAL nul;
    
  public:
    template <typename TM>
      SparseMatrixDynamic (const SparseMatrixTM<TM> & mat)
      : BaseSparseMatrix (mat) // , false)
    {
      width = mat.Width();
      bh = ngbla::Height<TM>(); // mat_traits<TM>::HEIGHT;
      bw = ngbla::Width<TM>();  // mat_traits<TM>::WIDTH;
      bs = bh*bw;
      is_complex = ngbla::IsComplex<TM>();
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

    AutoVector CreateRowVector() const override
    { throw Exception("CreateRowVector not implemented for SparseMatrixDynamic!"); }
    AutoVector CreateColVector() const override
    { throw  Exception("CreateColVector not implemented for SparseMatrixDynamic!"); }

    virtual tuple<int,int> EntrySizes() const override { return { bh, bw }; }

  };
  



  template <class TSCAL>
  class  NGS_DLL_HEADER SparseMatrixVariableBlocks : public S_BaseMatrix<TSCAL>
  {
  protected:
    size_t height, width, nblocks;
    Array<int> colnr;
    Array<TSCAL> data;
    Array<size_t> firsti_colnr, firsti_data;
    Array<int> cum_block_size;
    TSCAL nul;
    
  public:
    SparseMatrixVariableBlocks (const SparseMatrixTM<TSCAL> & mat);

    int VHeight() const override { return height; }
    int VWidth() const override { return width; }

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;

    AutoVector CreateRowVector () const override;
    AutoVector CreateColVector () const override;
  };



}
#endif
  
