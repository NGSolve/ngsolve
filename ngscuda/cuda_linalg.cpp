/*********************************************************************/
/* File:   cuda_linalg.cpp                                           */
/* Author: Joachim Schoeberl, Matthias Hochsteger                    */
/* Date:   11. Aug. 2014                                             */
/*********************************************************************/


#include <la.hpp>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cuda_runtime_api.h>

#include "cuda_linalg.hpp"

// TODO: why use own kernel instead of cublas?
//extern void SetScalar (double val, int n, double * dev_ptr);

namespace ngla
{

  /*
  UnifiedVector * dynamic_cast_UnifiedVector (BaseVector * x)
  {
    // cout << "my dynamic * cast" << endl;
    AutoVector * ax = dynamic_cast<AutoVector*> (x);
    if (ax)
      return dynamic_cast<UnifiedVector*> (&**ax);
    return dynamic_cast<UnifiedVector*> (x);
  }

  const UnifiedVector * dynamic_cast_UnifiedVector (const BaseVector * x)
  {
    // cout << "my dynamic const * cast" << endl;
    const AutoVector * ax = dynamic_cast<const AutoVector*> (x);
    if (ax)
      { 
        // cout << "is an autovector" << endl; 
        return dynamic_cast<const UnifiedVector*> (&**ax);
      }
    return dynamic_cast<const UnifiedVector*> (x);
  }
  
  UnifiedVector & dynamic_cast_UnifiedVector (BaseVector & x)
  {
    // cout << "my dynamic cast" << endl;
    AutoVector * ax = dynamic_cast<AutoVector*> (&x);
    if (ax)
      return dynamic_cast<UnifiedVector&> (**ax);
    return dynamic_cast<UnifiedVector&> (x);
  }
  const UnifiedVector & dynamic_cast_UnifiedVector (const BaseVector & x)
  {
    // cout << "my dynamic cast" << endl;
    const AutoVector * ax = dynamic_cast<const AutoVector*> (&x);
    if (ax)
      return dynamic_cast<const UnifiedVector&> (**ax);
    return dynamic_cast<const UnifiedVector&> (x);
  }
  */

  /*
  cublasHandle_t handle;
  cusparseHandle_t cusparseHandle;
  */

  cublasHandle_t Get_CuBlas_Handle ()
  {
    static Timer tblashandle("CUDA create cublas handle");
    RegionTimer reg(tblashandle);

    static cublasHandle_t handle;
    static bool first_call = true;

    if (first_call)
      {
        first_call = false;
        cublasCreate (&handle);
      }
    return handle;
  }
  cusparseHandle_t Get_CuSparse_Handle ()
  {
    static Timer tsparsehandle("CUDA create cusparse handle");
    RegionTimer reg(tsparsehandle);

    static cusparseHandle_t handle;
    static bool first_call = true;

    if (first_call)
      {
        first_call = false;
        cusparseCreate (&handle);
      }
    return handle;
  }


  // cublasHandle_t handle;
  // cusparseHandle_t cusparseHandle;


  UnifiedVector :: UnifiedVector (int asize)
  {
    initialize_unified(asize);
    (*this) = 0.0;
  }

  UnifiedVector :: UnifiedVector (const BaseVector& vec)
  {
    initialize_unified(vec.Size());
    (*this) = vec;
  }

  UnifiedVector :: UnifiedVector (const UnifiedVector & uvec)
  {
    initialize_unified(uvec.Size());
    (*this) = uvec;
  }

  UnifiedVector & UnifiedVector :: operator= (const UnifiedVector & uvec)
  {
    if (uvec.dev_uptodate)
      {
        cudaMemcpy (dev_data, uvec.dev_data, sizeof(double)*size, cudaMemcpyDeviceToDevice);    
        dev_uptodate = true;
        host_uptodate = false;
      }
    else if (uvec.host_uptodate)
      {
        VFlatVector<double> fv(size, host_data);
        VFlatVector<double> fv2 = uvec.FVDouble();

        fv = fv2;

        host_uptodate = true;
        dev_uptodate = false;
      }
    else
      {
        cerr << "operator= (BaseVector) : undefined vector" << endl;
      }
    return *this;
  }

  void UnifiedVector :: initialize_unified (size_t size)
  {
    this->size = size;

    host_data = new double[size];
    cudaMalloc((void**)&dev_data, size*sizeof(double));

    cusparseCreateDnVec (&descr, size, dev_data, CUDA_R_64F);

    host_uptodate = false;
    dev_uptodate = false;
  }
  
  UnifiedVector :: ~UnifiedVector ()
  {
    cerr << "dtor UnifiedVector" << endl;

    cusparseDestroyDnVec(descr);
    cudaFree(dev_data);
    delete[] host_data;
  }

  BaseVector & UnifiedVector :: operator= (double d)
  {
    for (int i = 0; i < size; i++) host_data[i] = d;
    /* ::SetScalar (d, size, dev_data); */

  
    cublasDscal(Get_CuBlas_Handle(), size, &d, dev_data, 1); 

    host_uptodate = true;
    dev_uptodate = true;
    
    return *this;
    /*
    host_uptodate = true;
    dev_uptodate = false;
    UpdateDevice();
    return *this;
    */
  }

  BaseVector & UnifiedVector :: operator= (const BaseVector & v2)
  {
    // moved to operator=(UnifiedVector &v2) (ran into problems while working on python bindings)
    /* UnifiedVector * uv2 = dynamic_cast<UnifiedVector*> (&v2); */
    /* if (uv2) */
    /*   { */
    /*     if (uv2->dev_uptodate) */
    /*       { */
    /*         cudaMemcpy (dev_data, uv2->dev_data, sizeof(double)*size, cudaMemcpyDeviceToDevice); */    
    /*         dev_uptodate = true; */
    /*         host_uptodate = false; */
    /*       } */
    /*     else if (uv2->host_uptodate) */
    /*       { */
    /*         VFlatVector<double> fv(size, host_data); */
    /*         fv = 1.0*v2; */
    /*         host_uptodate = true; */
    /*         dev_uptodate = false; */
    /*       } */
    /*     else */
    /*       { */
    /*         cerr << "operator= (BaseVector) : undefined vector" << endl; */
    /*       } */
    /*     return *this; */
    /*   } */

    VFlatVector<double> fv(size, host_data);
    fv = 1.0*v2;

    host_uptodate = true;
    dev_uptodate = false;
    return *this;
  }

  /* UnifiedVector & UnifiedVector :: operator= (UnifiedVector & uv2) */
  /* { */
  /*   if (uv2.dev_uptodate) */
  /*     { */
  /*       cudaMemcpy (dev_data, uv2.dev_data, sizeof(double)*size, cudaMemcpyDeviceToDevice); */    
  /*       dev_uptodate = true; */
  /*       host_uptodate = false; */
  /*     } */
  /*   else if (uv2.host_uptodate) */
  /*     { */
  /*       VFlatVector<double> fv(size, host_data); */
  /*       /1* fv = uv2.FVDouble(); *1/ */
  /*       cerr << "here should be a cpy..." << endl; */
  /*       host_uptodate = true; */
  /*       dev_uptodate = false; */
  /*     } */
  /*   else */
  /*     { */
  /*       cerr << "operator= (BaseVector) : undefined vector" << endl; */
  /*     } */
  /*   return *this; */
  /* } */

  size_t UnifiedVector :: Size () const throw()
  {
    return size;
  }

  /*
   *  TODO:
   *    avoid cpy between device and host
   *    maybe leave the update to the user?
   * */
  const double & UnifiedVector :: operator [] (const int ind) const
  {
    UpdateHost(); 
    return host_data[ind];
  }

  double & UnifiedVector :: operator [] (const int ind)
  {
    UpdateHost(); 
    dev_uptodate = false; // TODO: not sure, that this is the best approach
    return host_data[ind];
  }

  const cusparseDnVecDescr_t& UnifiedVector :: GetDescr() const
  {
    return descr;
  }

  cusparseDnVecDescr_t& UnifiedVector :: GetDescr()
  {
    return descr;
  }


  
  BaseVector & UnifiedVector :: Scale (double scal)
  {
    UpdateDevice();
    cublasDscal (Get_CuBlas_Handle(), size, &scal, dev_data, 1);
    host_uptodate = false;
    return *this;
  }

  BaseVector & UnifiedVector :: SetScalar (double scal)
  {
    (*this) = scal;
    return *this;
  }
  
  BaseVector & UnifiedVector :: Set (double scal, const BaseVector & v)
  {
    (*this) = 0.0;
    Add (scal, v);
    return *this;
  }
  
  
  BaseVector & UnifiedVector :: Add (double scal, const BaseVector & v)
  {
    const UnifiedVector * v2 = dynamic_cast<const UnifiedVector*> (&v);

    if (v2)
      {
        UpdateDevice();
        v2->UpdateDevice();

        cublasDaxpy (Get_CuBlas_Handle(), 
                           size, &scal, v2->dev_data, 1, dev_data, 1);
        host_uptodate = false;
      }
    else
      {
        UpdateHost();
        VFlatVector<> (size, host_data) += scal * v;
        dev_uptodate = false;
      }

    return *this;
  }

  // TODO: not optimal.
  BaseVector& UnifiedVector :: operator- () const
  {
    shared_ptr<UnifiedVector> uvecptr = make_shared<UnifiedVector>(*this);
    return uvecptr->Scale(-1);
  }
  
  double UnifiedVector :: InnerProduct (const BaseVector & v2) const
  {
    static Timer tdot("CUDA InnerProduct");
    RegionTimer reg(tdot);
    // cout << "Inner Prod" << endl;

    const UnifiedVector * uv2 = dynamic_cast<const UnifiedVector*> (&v2);
    if (uv2)
    {
      UpdateDevice();
      uv2->UpdateDevice();
      
      double res;
      cublasDdot (Get_CuBlas_Handle(), 
                        size, dev_data, 1, uv2->dev_data, 1, &res);
      return res;
    }

    FlatVector<> fv = FVDouble();
    FlatVector<> fv2 = v2.FVDouble();
    return ngbla::InnerProduct (fv, fv2);
  }


  ostream & UnifiedVector :: Print (ostream & ost) const
  {
    cout << "output unified vector of size " << size;
    cout << ", host = " << host_uptodate << ", dev = " << dev_uptodate << endl;
    if (!host_uptodate)
      cout << "host not up-to-date!" << endl;
    ost << FVDouble();
    return ost;
  }
  
  /* void UnifiedVector :: PrintDevice () const */
  /* { */
  /*   cerr << " i am the device print " << endl; */
  /*   int DSIZE = size * sizeof(double); */
  /*   double *tmp = (double*) malloc(DSIZE); */
  /*   cudaMemcpy(tmp, dev_data, DSIZE, cudaMemcpyDeviceToHost); */
  /*   cout << "device up-to-date: " << dev_uptodate << endl; */
  /*   for (int i=0; i<size; i++) */
  /*     cout << tmp[i] << endl; */
  /* } */

  void UnifiedVector :: UpdateHost () const
  {
    if (host_uptodate) return;
    if (!dev_uptodate) cout << "ERROR UnifiedVector::UpdateHost non is uptodate" << endl;
    cudaMemcpy (host_data, dev_data, sizeof(double)*size, cudaMemcpyDeviceToHost);    
    host_uptodate = true;
  }

  void UnifiedVector :: UpdateDevice () const
  {
    if (dev_uptodate) return;
    if (!host_uptodate) cout << "ERROR UnifiedVector::UpdateDevice non is uptodate" << endl;
    cout << "Host2Device copy !!!!!!!!!!!!!!" << endl;
    cudaMemcpy (dev_data, host_data, sizeof(double)*size, cudaMemcpyHostToDevice);
    dev_uptodate = true;
  }
  
  FlatVector<double> UnifiedVector :: FVDouble () const
  {
    /* UpdateHost(); */
    return FlatVector<> (size, host_data);
  }
  
  FlatVector<Complex> UnifiedVector :: FVComplex () const
  {
    throw Exception ("unified complex not yet supported");
  }
    
  void * UnifiedVector :: Memory() const throw()
  { 
    UpdateHost(); 
    return host_data;
  }

  AutoVector UnifiedVector :: CreateVector () const
  {
    return make_unique<UnifiedVector> (size);
  }

  
  void UnifiedVector :: GetIndirect (const FlatArray<int> & ind, 
            const FlatVector<double> & v) const
  {
    cout << "UnifiedVector :: GetIndirect not supported" << endl;
  }
  void UnifiedVector :: GetIndirect (const FlatArray<int> & ind, 
            const FlatVector<Complex> & v) const
  {
    cout << "UnifiedVector :: GetIndirect not supported" << endl;
  }


  /* BaseVector& operator+ (const UnifiedVector& v1, const BaseVector& v2) */
  /* { */
  /*   shared_ptr<UnifiedVector> res = make_shared<UnifiedVector>(v1); */
  /*   return res->Add(1, v2); */
  /* } */
  /* BaseVector& operator+ (const BaseVector& v1, const UnifiedVector& v2) */
  /* { */
  /*   return v2 + v1; */
  /* } */
  /* BaseVector& operator+ (const UnifiedVector& v1, const UnifiedVector& v2) */
  /* { */
  /*   shared_ptr<UnifiedVector> res = make_shared<UnifiedVector>(v1); */
  /*   return res->Add(1, v2); */
  /* } */

  /* BaseVector& operator- (const UnifiedVector& v1, const BaseVector& v2) */
  /* { */
  /*   shared_ptr<UnifiedVector> res = make_shared<UnifiedVector>(v1); */
  /*   return res->Add(-1, v2); */
  /* } */
  /* BaseVector& operator- (const BaseVector& v1, const UnifiedVector& v2) */
  /* { */
  /*   shared_ptr<UnifiedVector> res = make_shared<UnifiedVector>(v1); */
  /*   return res->Add(-1, v2); */
  /* } */
  /* BaseVector& operator- (const UnifiedVector& v1, const UnifiedVector& v2) */
  /* { */
  /*   shared_ptr<UnifiedVector> res = make_shared<UnifiedVector>(v1); */
  /*   return res->Add(-1, v2); */
  /* } */

  /* BaseVector& operator* (double d, const BaseVector& v) */
  /* { */
  /*   shared_ptr<UnifiedVector> res = make_shared<UnifiedVector>(v); */
  /*   return res->Scale(d); */
  /* } */

  /* BaseVector& operator* (double d, const UnifiedVector& v) */
  /* { */
  /*   shared_ptr<UnifiedVector> res = make_shared<UnifiedVector>(v); */
  /*   return res->Scale(d); */
  /* } */





  DevSparseMatrix :: DevSparseMatrix (const SparseMatrix<double> & mat)
  {
    height = mat.Height();
    width = mat.Width();
    nze = mat.NZE();

    /*
    descr = new cusparseMatDescr_t;
    cusparseCreateMatDescr (descr);

    cusparseSetMatType(*descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(*descr, CUSPARSE_INDEX_BASE_ZERO);
    */

    cout << "create device sparse matrix, n = " << height << ", nze = " << nze << endl;
    
    Array<int> temp_ind (mat.Height()+1); 
    for (int i = 0; i <= mat.Height(); i++) temp_ind[i] = mat.First(i); // conversion to 32-bit integer

    cudaMalloc ((void**)&dev_ind, (mat.Height()+1) * sizeof(int));
    cudaMalloc ((void**)&dev_col, (mat.NZE()) * sizeof(int));
    cudaMalloc ((void**)&dev_val, (mat.NZE()) * sizeof(double));

    cudaMemcpy (dev_ind, &temp_ind[0], (mat.Height()+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy (dev_col, &mat.GetRowIndices(0)[0], mat.NZE()*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy (dev_val, &mat.GetRowValues(0)[0], mat.NZE()*sizeof(double), cudaMemcpyHostToDevice);

    /* cerr << "CREATING!" << endl; */
    /* int *host_row = (int*) malloc((mat.Height()+1) * sizeof(int)); */
    /* int *host_col = (int*) malloc(mat.Width() * sizeof(int)); */
    /* double *host_val = (double*) malloc(mat.NZE() * sizeof(double)); */
    /* cudaMemcpy(host_row, dev_ind, (mat.Height() + 1) * sizeof(int), cudaMemcpyDeviceToHost); */
    /* cudaMemcpy(host_col, dev_col, mat.Width() * sizeof(int), cudaMemcpyDeviceToHost); */
    /* cudaMemcpy(host_val, dev_val, mat.NZE() * sizeof(double), cudaMemcpyDeviceToHost); */
    /* for (int i=0; i<mat.Height() +1 ; i++) */
    /* { */
    /*   cerr << host_row[i] << " "; */
    /* } */
    /* cerr << endl; */
    /* for (int i=0; i<mat.Width(); i++) */
    /* { */
    /*   cerr << host_col[i] << " "; */
    /* } */
    /* cerr << endl; */
    /* for (int i=0; i<mat.NZE(); i++) */
    /* { */
    /*   cerr << host_val[i] << " "; */
    /* } */
    /* cerr << endl; */


    /* for (int i=0; i<mat.Height()+1; i++) */
    /* { */
    /*   cerr << temp_ind[i] << " "; */
    /* } */
    /* cerr << endl; */
    /* for (int i=0; i<mat.NZE()+1; i++) */
    /* { */
    /*   cerr << mat.GetRowIndices(0)[i] << " "; */
    /* } */
    /* cerr << endl; */

    cusparseCreateCsr(&descr, height, width, nze,
                      dev_ind, dev_col, dev_val,
                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
                      CUDA_R_64F);
  }

  DevSparseMatrix :: ~DevSparseMatrix ()
  {
    cusparseDestroySpMat(descr);
  }
  
  void DevSparseMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer tprel("CUDA PREL");
    static Timer tbuffer("CUDA BUFFER");
    static Timer tmv("CUDA MV");
    static Timer tbufferallocate("CUDA BUFFER ALLOCATE");

    tprel.Start();
    // cout << "device mult sparse" << endl;
    const UnifiedVector & ux = dynamic_cast<const UnifiedVector&> (x);
    UnifiedVector & uy = dynamic_cast<UnifiedVector&> (y);

    ux.UpdateDevice();
    uy = 0.0;
    uy.UpdateDevice();

    double alpha= 1;
    double beta = 0;
    tprel.Stop();

    // deprecated
    /*
    cusparseDcsrmv (Get_CuSparse_Handle(), 
                    CUSPARSE_OPERATION_NON_TRANSPOSE, height, width, nze, 
        &alpha, *descr, 
        dev_val, dev_ind, dev_col, 
        ux.dev_data, &beta, uy.dev_data);
    */

    tbuffer.Start();
    size_t bufferSize = 0;
    void* dBuffer = NULL;

    cusparseSpMV_bufferSize(Get_CuSparse_Handle(), CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha, descr, ux.descr, &beta, uy.descr, CUDA_R_64F,
                            CUSPARSE_MV_ALG_DEFAULT, &bufferSize);
    tbuffer.Stop();
    tbufferallocate.Start();
    cudaMalloc(&dBuffer, bufferSize);
    tbufferallocate.Stop();

    cusparseStatus_t status;
    tmv.Start();
    status = cusparseSpMV(Get_CuSparse_Handle(), 
                 CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, descr,
                 ux.descr, &beta, uy.descr, CUDA_R_64F,
                 CUSPARSE_SPMV_ALG_DEFAULT, dBuffer);
    tmv.Stop();
    /* cerr << "Status: " << cusparseGetErrorString(status) << endl; */
    
    uy.host_uptodate = false;
    // cout << "mult complete" << endl;
  }


  void DevSparseMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    // cout << "device multadd sparse matrix" << endl;
    const UnifiedVector & ux = dynamic_cast<const UnifiedVector&> (x);
    UnifiedVector & uy = dynamic_cast<UnifiedVector&> (y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    double alpha= 1;
    //double beta = 1;
    double beta = s; // s was not in use. I guess this is where it's supposed to be

    // deprecated
    /*
    cusparseDcsrmv (Get_CuSparse_Handle(), 
                    CUSPARSE_OPERATION_NON_TRANSPOSE, height, width, nze, 
        &alpha, *descr, 
        dev_val, dev_ind, dev_col, 
        ux.dev_data, &beta, uy.dev_data);
    */

    cusparseSpMatDescr_t matA;
    size_t bufferSize = 0;
    void* dBuffer = NULL;

    cusparseSpMV_bufferSize(Get_CuSparse_Handle(), CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha, matA, ux.descr, &beta, uy.descr, CUDA_R_64F,
                            CUSPARSE_MV_ALG_DEFAULT, &bufferSize);
    cudaMalloc(&dBuffer, bufferSize);

    cusparseSpMV(Get_CuSparse_Handle(), 
                 CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA,
                 ux.descr, &beta, uy.descr, CUDA_R_64F,
                 CUSPARSE_MV_ALG_DEFAULT, &bufferSize);

    uy.host_uptodate = false;

    // cout << "mult complete" << endl;
  }
  
  void DevSparseMatrix :: Scale (double d)
  {
    throw Exception("not implemented yet.");
  }
  void DevSparseMatrix :: Mult (const DevSparseMatrix& a)
  {
    // cusparseSp
    // cusparseSpGEMMreuse
    throw Exception("not implemented yet.");
  }

  /* DevSparseMatrix& operator+ (const DevSparseMatrix& a, const DevSparseMatrix& b) */
  /* { */
  /*   throw Exception("TODO"); */
  /* } */

  /* DevSparseMatrix& operator- (const DevSparseMatrix& a, const DevSparseMatrix& b) */
  /* { */
  /*   throw Exception("TODO"); */
  /* } */

  /* DevSparseMatrix& operator* (const DevSparseMatrix& a, const DevSparseMatrix& b) */
  /* { */
  /*   throw Exception("TODO"); */
  /* } */
  
  /* DevSparseMatrix& operator* (const DevSparseMatrix& a, double d) */
  /* { */
  /*   throw Exception("TODO"); */
  /* } */


  // TODO:
  //   redo. similiar to DevSparseMatrix
  /* DevJacobiPreconditioner :: DevJacobiPreconditioner (const SparseMatrix<double> & mat, */
                  /* const BitArray & freedofs) */
  /* { */
  /*   height = mat.Height(); */
  /*   width = mat.Height(); */
  /*   nze = mat.Height(); */

  /*   descr = new cusparseMatDescr_t; */
  /*   cusparseCreateMatDescr (descr); */

  /*   cusparseSetMatType(*descr, CUSPARSE_MATRIX_TYPE_GENERAL); */
  /*   cusparseSetMatIndexBase(*descr, CUSPARSE_INDEX_BASE_ZERO); */

  /*   cout << "create Jacobi preconditioner" << endl; */
    
  /*   Array<int> temp_ind (height+1); */ 
  /*   Array<int> temp_cols (height); */
  /*   Array<double> temp_vals (height); */

  /*   for (int i = 0; i <= height; i++) temp_ind[i] = i; */ 

  /*   for (int i = 0; i < height; i++) */
  /*     { */
  /* temp_cols[i] = i; */
  /* if (freedofs.Test(i)) */
    /* temp_vals[i] = 1.0 / mat(i,i); */
  /* else */
    /* temp_vals[i] = 0.0; */
  /*     } */

  /*   cudaMalloc ((void**)&dev_ind, (height+1) * sizeof(int)); */
  /*   cudaMalloc ((void**)&dev_col, (height) * sizeof(int)); */
  /*   cudaMalloc ((void**)&dev_val, (height) * sizeof(double)); */

  /*   cudaMemcpy (dev_ind, &temp_ind[0], (height+1)*sizeof(int), cudaMemcpyHostToDevice); */
  /*   cudaMemcpy (dev_col, &temp_cols[0], height*sizeof(int), cudaMemcpyHostToDevice); */
  /*   cudaMemcpy (dev_val, &temp_vals[0], height*sizeof(double), cudaMemcpyHostToDevice); */
  /* } */
  
  /* void DevJacobiPreconditioner :: Mult (const BaseVector & x, BaseVector & y) const */
  /* { */
  /*   // cout << "device mult precond" << endl; */

  /*   const UnifiedVector & ux = dynamic_cast<UnifiedVector*> (x); */
  /*   UnifiedVector & uy = dynamic_cast<UnifiedVector*> (y); */

  /*   ux.UpdateDevice(); */
  /*   uy = 0.0; */
  /*   uy.UpdateDevice(); */

  /*   double alpha= 1; */
  /*   double beta = 0; */

    /* // TODO: fix this */
  /*   /1* cusparseDcsrmv (Get_CuSparse_Handle (), *1/ */ 
  /*   /1*                 CUSPARSE_OPERATION_NON_TRANSPOSE, height, width, nze, *1/ */ 
        /* /1* &alpha, *descr, *1/ */ 
        /* /1* dev_val, dev_ind, dev_col, *1/ */ 
        /* /1* ux.dev_data, &beta, uy.dev_data); *1/ */

  /*   uy.host_uptodate = false; */

  /*   // cout << "mult complete" << endl; */
  /* } */


  /* void DevJacobiPreconditioner :: MultAdd (double s, const BaseVector & x, BaseVector & y) const */
  /* { */
  /*   cout << "device multadd precond" << endl; */

  /*   const UnifiedVector & ux = dynamic_cast<UnifiedVector*>(x); */
  /*   UnifiedVector & uy = dynamic_cast<UnifiedVector*>(y); */

  /*   ux.UpdateDevice(); */
  /*   uy.UpdateDevice(); */

  /*   double alpha= 1; */
  /*   double beta = 1; */


    /* // TODO: fix this */
  /*   /1* cusparseDcsrmv (Get_CuSparse_Handle (), *1/ */ 
  /*   /1*                 CUSPARSE_OPERATION_NON_TRANSPOSE, height, width, nze, *1/ */ 
        /* /1* &alpha, *descr, *1/ */ 
        /* /1* dev_val, dev_ind, dev_col, *1/ */ 
        /* /1* ux.dev_data, &beta, uy.dev_data); *1/ */

  /*   uy.host_uptodate = false; */
  /*   cout << "mult complete" << endl; */
  /* } */
  















  /*
  class InitCuBlasHandle
  {
  public:
    InitCuBlasHandle ()
    {
      cout << "create cublas handle" << endl;
      cublasCreate (&handle);
      cout << "creaet cusparse handle" << endl;
      cusparseCreate(&cusparseHandle);
    }
  };
  InitCuBlasHandle init;
  */

  DevDMatrix :: DevDMatrix ()
  { }

  DevDMatrix :: DevDMatrix (const Matrix<>& mat)
  {
    height = mat.Height();
    width = mat.Width();

    host_data = mat.Data();
    cudaMalloc((void**) &dev_data, height * width * sizeof(double));
    cudaMemcpy(dev_data, host_data, height * width * sizeof(double), cudaMemcpyHostToDevice);

    /* double* tmp = new double [height * width]; */
    /* cudaMemcpy(tmp, dev_data, height * width * sizeof(double), cudaMemcpyDeviceToHost); */
    /* for (int i=0; i<height*width; i++) */
    /*   cerr << tmp[i] << endl; */
    /* delete[] tmp; */

    // in case we need cusparse operations, we could add an cusparse descriptor
    /* cusparseCreateDnMat (&descr, mat.Height(), mat.Width(), mat.Height(), */
    /*                      dev_data, CUDA_R_64F, CUSPARSE_ORDER_ROW); */
  }

  DevDMatrix :: ~DevDMatrix ()
  {
    /* cusparseDestroyDnMat (descr); */
    cudaFree(dev_data);
  }

  AutoVector DevDMatrix :: CreateRowVector () const
  {
    return UnifiedVector(width).CreateVector();
  }

  AutoVector DevDMatrix :: CreateColVector () const
  {
    return UnifiedVector(height).CreateVector();
  }

  // TODO: currently in-place. change?
  void DevDMatrix :: Add (const DevDMatrix& b)
  {
    double alpha = 1;
    /* double beta = 1; */

    cublasAxpyEx (Get_CuBlas_Handle(), height*width, &alpha, CUDA_R_64F,
                  b.DevData(), CUDA_R_64F, 1, dev_data, CUDA_R_64F, 1, CUDA_R_64F);

    /* cublasDgeam(Get_CuBlas_Handle(), CUBLAS_OP_N, CUBLAS_OP_N, height, width, */
    /*             &alpha, dev_data, width, &beta, b.DevData(), width, dev_data, width); */
  }

  // TODO: fix
  void DevDMatrix :: Mult (const DevDMatrix& b, DevDMatrix& c)
  {
    int m = height;
    int n = width;
    int k = b.Width();
    
    double alpha = 1;
    double beta = 0;

    throw Exception ("causes errors (prbly due to col-major order)")

    cerr << "dense device matrix prod " << m << "x" << n << " and " << n << "x" << k << endl;
    cerr << b.Height() << " " << b.Width() << endl;
    cublasDgemm (Get_CuBlas_Handle(), CUBLAS_OP_N, CUBLAS_OP_N, m, n, k,
                 &alpha, dev_data, m, b.DevData(), k, &beta, c.DevData(), m);
  }

  void DevDMatrix :: Scale (double d)
  {
    cublasScalEx(Get_CuBlas_Handle(), height*width, &d, CUDA_R_64F, 
                 dev_data, CUDA_R_64F, 1, CUDA_R_64F);
  }

  void DevDMatrix :: SetZero ()
  {
    double alpha = 0;
    double beta = 0;

    // special case
    // see cublas documentation
    cublasDgeam (Get_CuBlas_Handle(), CUBLAS_OP_N, CUBLAS_OP_N, height, width,
                 &alpha, nullptr, height, &beta, nullptr, height, dev_data, height);
  }

  double* DevDMatrix :: DevData () const
  {
    return dev_data;
  }

  ostream & DevDMatrix :: Print (ostream & ost) const
  {
    cout << "output dense device Matrix of size " << height << "x" << width << endl;
    Matrix<> tmp (height, width);
    cudaMemcpy(tmp.Data(), dev_data, height * width * sizeof(double), cudaMemcpyDeviceToHost);
    ost << tmp << endl;
    return ost;
  }

  DevEBEMatrix :: DevEBEMatrix (const ConstantElementByElementMatrix& ebemat)
    : devmat(ebemat.GetMatrix()), col_dnums(ebemat.GetColDNums()), row_dnums(ebemat.GetRowDNums())
  { }

  DevEBEMatrix :: ~DevEBEMatrix ()
  { }

  AutoVector DevEBEMatrix :: CreateRowVector () const
  {
    return UnifiedVector(width).CreateVector();
  }

  AutoVector DevEBEMatrix :: CreateColVector () const
  {
    return UnifiedVector(height).CreateVector();
  }

  void DevEBEMatrix :: MultAdd (double s, const UnifiedVector& x, UnifiedVector& y)
  {
    static Timer timer("Dev-EBE-Matrix::MultAdd");
    RegionTimer reg(timer);

    size_t maxs = 0;
    for (size_t i=0; i<col_dnums.Size(); i++)
      maxs = max2 (maxs, col_dnums[i].Size());

    /* ebe_multadd_kernel(); */

    throw Exception("not implemented yet.");
  }

  void DevEBEMatrix :: Scale (double d)
  {
    throw Exception("not implemented yet.");
  }


  shared_ptr<DevMatrix> CreateDevMatrix (BaseMatrix & mat)
  /* DevMatrix :: DevMatrix (BaseMatrix & mat) */
  {
    cerr << "creating dev matrix." << endl;
    if (typeid(mat) == typeid(SparseMatrix<double>))
    {
      SparseMatrix<double>& sparse_mat = dynamic_cast<SparseMatrix<double>&>(mat);
      cerr << "creating sparse matrix" << endl;
      return make_shared<DevSparseMatrix>(sparse_mat);
      /* *this = DevSparseMatrix(sparse_mat); */
      /* return *this; */
    }
    else if (typeid(mat) == typeid(ConstantElementByElementMatrix))
    {
      ConstantElementByElementMatrix& ebe_mat = dynamic_cast<ConstantElementByElementMatrix&>(mat);
      cerr << "creating ebe matrix" << endl;
      return make_shared<DevEBEMatrix>(ebe_mat);
    }
    else
    {
      throw Exception(string("matrix type not supported: ") + typeid(mat).name());
    }
  }

  // TODO: also complex
  shared_ptr<DevMatrix> CreateDevMatrix (Matrix<> & mat)
  {
    Matrix<double>& dmat = dynamic_cast<Matrix<double>&>(mat);
    cerr << "creating dense device matrix" << endl;
    return make_shared<DevDMatrix>(dmat);
  }

  BaseVector& operator* (const UnifiedVector& v, const DevSparseMatrix& mat)
  {
    shared_ptr<UnifiedVector> res = make_shared<UnifiedVector>(v.Size());
    mat.Mult(v, *res);
    return *res;
  }


}
