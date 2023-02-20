#ifndef CUDA_LINALG_HPP
#define CUDA_LINALG_HPP

#include <la.hpp>

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>

#include "cuda_ngstd.hpp"
using namespace ngs_cuda;

namespace ngla
{
  cublasHandle_t Get_CuBlas_Handle ();
}

  
#include "cuda_ngbla.hpp"


#include "linalg_kernels.hpp"

#include "unifiedvector.hpp"


// from CUDA C++ Programming Guide:
// https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
#if __CUDA_ARCH__ < 600
__device__ inline double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif



namespace ngla
{
  using namespace ngs_cuda;
  
  void InitCuLinalg();

  class DevSparseMatrix;


  /* AutoVector CreateUnifiedVector(size_t size); */

  class DevMatrix : public BaseMatrix
  {
  public:
    DevMatrix() { }

    virtual AutoVector CreateRowVector() const { return make_unique<UnifiedVector>(Width()); }
    virtual AutoVector CreateColVector() const { return make_unique<UnifiedVector>(Height()); }
  };

  shared_ptr<BaseMatrix> CreateDevMatrix (BaseMatrix &mat);
  shared_ptr<BaseMatrix> CreateDevMatrix (Matrix<> &mat);


  class DevSparseMatrix : public DevMatrix
  {
  protected:
    //cusparseMatDescr_t * descr;
    cusparseSpMatDescr_t descr;
    int * dev_ind;
    int * dev_col;
    double * dev_val;
    int height, width, nze;
  public:
    DevSparseMatrix () { }
    DevSparseMatrix (const SparseMatrix<double> & mat);
    virtual ~DevSparseMatrix ();

    virtual void Mult (const BaseVector & x, BaseVector & y) const;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;

    virtual int VHeight() const { return height; }
    virtual int VWidth() const { return width; }
  };


  class DevDiagonalMatrix : public DevMatrix
  {
  protected:
    const UnifiedVector diag;

  public:
    DevDiagonalMatrix (const UnifiedVector _diag) : diag(_diag) { }

    virtual void Mult (const BaseVector & x, BaseVector & y) const;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    virtual int VHeight() const { return diag.Size(); }
    virtual int VWidth() const { return diag.Size(); }
  };


  
  class DevConstantElementByElementMatrix : public DevMatrix
  {
    size_t h, w; // big matrix shape

    Matrix<Dev<double>> devmat;
    
    DevTable<int> rowdnums, coldnums;
    DevDataTable<int> row_coloring, col_coloring;

    bool disjoint_rows, disjoint_cols;
    size_t numblocks;
  public:
    DevConstantElementByElementMatrix (const ConstantElementByElementMatrix & mat);
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    int VHeight() const override { return h; }
    int VWidth() const override { return w; }
  };
  


  class DevBlockDiagonalMatrixSoA : public DevMatrix
  {
    double * dev_data; // Tensor<3> blockdiag;  
    int blocks, dimy, dimx;
    Matrix<double> nonzero;
    Array<Dev<int>> indices, indices_trans;
 public:
    DevBlockDiagonalMatrixSoA (const BlockDiagonalMatrixSoA & mat);
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    int VHeight() const override { return dimy*blocks; }
    int VWidth() const override { return dimx*blocks; }
  };

  class DevBlockJacobiMatrix : public DevMatrix
  {
    double h, w;
    Array<Dev<int>> indices;
    Array<Dev<double>> matrices;
    Array<Dev<BlockJacobiCtr>> ctrstructs;
 public:
    DevBlockJacobiMatrix (const BlockJacobiPrecond<double> & mat);
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    // void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    int VHeight() const override { return h; }
    int VWidth() const override { return w; }
  };

  
  class DevEmbeddedMatrix : public EmbeddedMatrix
  {
  public:
    using EmbeddedMatrix::EmbeddedMatrix;
    AutoVector CreateColVector() const override { return make_unique<UnifiedVector>(Height()); }      
  };
  
  class DevEmbeddedTransposeMatrix : public EmbeddedTransposeMatrix
  {
  public:
    using EmbeddedTransposeMatrix::EmbeddedTransposeMatrix;
    AutoVector CreateRowVector() const override { return make_unique<UnifiedVector>(Width()); }      
  };


  class DevProjector : public DevMatrix
  {
  private:
    shared_ptr<DevBitArray> bits;
    bool keep_values;
  public:
    DevProjector (const Projector & proj)
      : bits(make_shared<DevBitArray>(*proj.Mask())), 
        keep_values(proj.KeepValues()) { ; }

    void Mult (const BaseVector & x, BaseVector & y) const;
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    void Project (BaseVector & x) const;

    virtual int VHeight() const override { return bits->Size(); }
    virtual int VWidth() const override { return bits->Size(); }

    AutoVector CreateRowVector() const override
    { throw Exception("CreateRowVector not implemented for DevProjector!"); }
    AutoVector CreateColVector() const override
    { throw Exception("CreateColVector not implemented for DevProjector!"); }

  };
  



  
  
  shared_ptr<DevSparseMatrix> MatMult (const DevSparseMatrix& mata, const DevSparseMatrix& matb);

  // dense device matrix
  class DevDMatrix : public DevMatrix
  {
  private:
    size_t height, width;
    double* dev_data;
    /* cusparseDnMatDescr_t descr; */

  public:
    DevDMatrix ();
    DevDMatrix (size_t height, size_t width);
    DevDMatrix (const Matrix<>& mat);
    DevDMatrix (const DevDMatrix& mat);
    ~DevDMatrix ();

    /* const DevDMatrix & operator= (const Expr<TBxx> & mat) const; */
    const DevDMatrix & operator= (double d) const;
    const DevDMatrix & operator= (const DevDMatrix & mat) const;

    int VHeight() const { return height; }
    int VWidth() const { return width; }

    virtual AutoVector CreateRowVector () const;
    virtual AutoVector CreateColVector () const;

    virtual void Add (const BaseMatrix& b);
    virtual void Scale (double d);
    virtual void Mult (const BaseVector& x, BaseVector& y) const;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    void SetZero ();

    double* DevData () const;

    virtual ostream & Print (ostream & ost) const;    
  };

  shared_ptr<DevDMatrix> MatMult (const DevDMatrix& mata, const DevDMatrix& matb);

  class DevEBEMatrix : public DevMatrix
  {
  private:
    size_t height, width;
    DevDMatrix devmat;
    Table<int> col_dnums, row_dnums;
    bool disjointrows, disjointcols;
  public:
    DevEBEMatrix (const ConstantElementByElementMatrix& mat);
    ~DevEBEMatrix ();

    virtual AutoVector CreateRowVector () const;
    virtual AutoVector CreateColVector () const;

    void MultAdd (double s, const UnifiedVector & x, UnifiedVector & y);
    /* void Scale (double d); */
  };

  // currently uses Mult and MultAdd of DevSparseMatrix
  class DevJacobiPrecond : public DevSparseMatrix
  {
  private:
    /* cusparseSpMatDescr_t descr; */
    shared_ptr<BitArray> inner;
    double* dev_invdiag;
  public:
    DevJacobiPrecond (const SparseMatrix<double> & amat, 
      shared_ptr<BitArray> ainner=nullptr, bool use_par=true);

    /* DevJacobiPrecond (const JacobiPrecond<double> & amat); */

    virtual ~DevJacobiPrecond ();

    /* void Mult (const BaseVector & x, BaseVector & y) const; */
    /* void MultAdd (double s, const BaseVector & x, BaseVector & y) const; */

    /* int VHeight() const override { return height; } */
    /* int VWidth() const override { return height; } */

  };


  // old version
  /* class DevJacobiPreconditioner : public BaseMatrix */
  /* { */
  /*   // should be like this: */
  /*   // double * dev_diag; */
  /*   // int size; */

  /*   // stored as sparse matrix ... */
  /*   /1* cusparseMatDescr_t * descr; *1/ */

    /* // used to be cusparseMatSpDescr_t... not sure about the difference */
    /* cusparseSpMatDescr_t descr; */
  /*   int * dev_ind; */
  /*   int * dev_col; */
  /*   double * dev_val; */
  /*   int height, width, nze; */
    

  /* public: */
  /*   DevJacobiPreconditioner (const SparseMatrix<double> & mat, const BitArray & freedofs); */
  /*   virtual void Mult (const BaseVector & x, BaseVector & y) const; */
  /*   virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const; */
  /* }; */

}


#endif
