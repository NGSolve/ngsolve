#ifndef FILE_DEVICE_DIAGONALMATRIX_HPP
#define FILE_DEVICE_DIAGONALMATRIX_HPP

/*********************************************************************/
/* File:   device_diagonalmatrix.hpp                                 */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   4. Sep. 2026                                              */
/*********************************************************************/

/*
  Backend-independent diagonal matrix on the gpu, the device counterpart
  of DiagonalMatrix and of JacobiPrecond (whose inverse diagonal is
  computed on the host). The diagonal is a DeviceVector<T>.
*/

#include "devicevector.hpp"
#include "diagonalmatrix.hpp"
#include "jacobi.hpp"

namespace ngla
{

  template <typename T>
  class NGS_DLL_HEADER DeviceDiagonalMatrix : public BaseMatrix
  {
  protected:
    DeviceVector<T> diag;

  public:
    // entries are converted to T
    template <typename TS>
    DeviceDiagonalMatrix (FlatVector<TS> adiag);
    virtual ~DeviceDiagonalMatrix () { }

    virtual int VHeight() const override { return diag.Size(); }
    virtual int VWidth() const override { return diag.Size(); }
    virtual bool IsComplex() const override { return false; }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override
    { Mult (x, y); }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    { MultAdd (s, x, y); }
    void Launch (const BaseVector & x, BaseVector & y, T s, T beta) const;

    virtual BaseVector & AsVector() override { return diag; }
    virtual const BaseVector & AsVector() const override { return diag; }
    const DeviceVector<T> & Diag() const { return diag; }

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    virtual ostream & Print (ostream & ost) const override;
  };


  /*
    Device counterpart of BlockDiagonalMatrixSoA: dimy x dimx blocks of
    'blocks' diagonal entries each, block index fastest. Only the non-zero
    (i,j) block positions are visited, as compressed row lists for A and A^T.
  */
  template <typename T>
  class NGS_DLL_HEADER DeviceBlockDiagonalMatrixSoA : public BaseMatrix
  {
  protected:
    int blocks, dimy, dimx;
    MemType memtype;
    ngs_gpu::TypedBuffer<T> dev_data;                 // (dimy*dimx) x blocks
    ngs_gpu::TypedBuffer<int> first, aind, xind;      // per row j: entries k, a-row aind[k], x-row xind[k]
    ngs_gpu::TypedBuffer<int> firstT, aindT, xindT;   // the same for the transpose

    void Launch (const BaseVector & x, BaseVector & y, T s, T beta, bool trans) const;

  public:
    DeviceBlockDiagonalMatrixSoA (const BlockDiagonalMatrixSoA & mat);
    virtual ~DeviceBlockDiagonalMatrixSoA () { }

    virtual int VHeight() const override { return blocks*dimy; }
    virtual int VWidth() const override { return blocks*dimx; }
    virtual bool IsComplex() const override { return false; }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    virtual ostream & Print (ostream & ost) const override;
  };


#if !defined(FILE_DEVICE_DIAGONALMATRIX_CPP)
  extern template class DeviceDiagonalMatrix<double>;
  extern template class DeviceDiagonalMatrix<float>;
  extern template class DeviceBlockDiagonalMatrixSoA<double>;
  extern template class DeviceBlockDiagonalMatrixSoA<float>;
#endif
}

#endif
