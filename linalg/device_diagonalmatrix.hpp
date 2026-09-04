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

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    { MultAdd (s, x, y); }

    virtual BaseVector & AsVector() override { return diag; }
    virtual const BaseVector & AsVector() const override { return diag; }
    const DeviceVector<T> & Diag() const { return diag; }

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    virtual ostream & Print (ostream & ost) const override;
  };


#if !defined(FILE_DEVICE_DIAGONALMATRIX_CPP)
  extern template class DeviceDiagonalMatrix<double>;
  extern template class DeviceDiagonalMatrix<float>;
#endif
}

#endif
