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
    cerr << "UnifiedVector (int) ctor of size " << asize << endl;
    this->size = asize;

    host_data = new double[size];
    cudaMalloc((void**)&dev_data, size*sizeof(double));

    cusparseCreateDnVec (&descr, size, dev_data, CUDA_R_64F);

    /* host_uptodate = false; */
    // not the best solution, since the values in host_data are random host_uptodate and dev_uptodate false causes problems
    // alternative: check for this case in FVDouble()
    host_uptodate = true;
    dev_uptodate = true;
  }

  UnifiedVector :: UnifiedVector (const BaseVector& vec) : UnifiedVector(vec.Size())
  {
    (*this) = vec;
  }

  // TODO: remove
  /* UnifiedVector :: UnifiedVector (const UnifiedVector & uvec) */
  /* { */
  /*   initialize_unified(uvec.Size()); */
  /*   (*this) = uvec; */
  /* } */

  
  UnifiedVector :: ~UnifiedVector ()
  {
    cerr << "dtor UnifiedVector" << endl;

    cusparseDestroyDnVec(descr);
    cudaFree(dev_data);
    delete[] host_data;
  }

  BaseVector & UnifiedVector :: operator= (double d)
  {
    /* for (int i = 0; i < size; i++) host_data[i] = d; */
    host_uptodate = false;

    /* ::SetScalar (d, size, dev_data); */
    cublasDscal(Get_CuBlas_Handle(), size, &d, dev_data, 1); 
    dev_uptodate = true;
    
    return *this;
    /*
    host_uptodate = true;
    dev_uptodate = false;
    UpdateDevice();
    return *this;
    */
  }

  /* UnifiedVector & UnifiedVector :: operator= (const UnifiedVector & uvec) */
  /* { */
  /*   if (uvec.dev_uptodate) */
  /*     { */
  /*       cudaMemcpy (dev_data, uvec.dev_data, sizeof(double)*size, cudaMemcpyDeviceToDevice); */    
  /*       dev_uptodate = true; */
  /*       host_uptodate = false; */
  /*     } */
  /*   else if (uvec.host_uptodate) */
  /*     { */
  /*       VFlatVector<double> fv(size, host_data); */
  /*       VFlatVector<double> fv2 = uvec.FVDouble(); */

  /*       fv = fv2; */

  /*       host_uptodate = true; */
  /*       dev_uptodate = false; */
  /*     } */
  /*   else */
  /*     { */
  /*       cerr << "operator= (BaseVector) : undefined vector" << endl; */
  /*     } */
  /*   return *this; */
  /* } */


  BaseVector & UnifiedVector :: operator= (const BaseVector & v2)
  {
    // moved to operator=(UnifiedVector &v2) (ran into problems while working on python bindings)
    const UnifiedVector * uv2 = dynamic_cast<const UnifiedVector*> (&v2);
    if (uv2)
      {
        if (uv2->dev_uptodate)
          {
            cudaMemcpy (dev_data, uv2->dev_data, sizeof(double)*size, cudaMemcpyDeviceToDevice);    
            dev_uptodate = true;
            host_uptodate = false;
          }
        else if (uv2->host_uptodate)
          {
            /* VFlatVector<double> fv(size, host_data); */
            /* fv = 1.0*v2; */
            FVDouble() = uv2->FVDouble();
            host_uptodate = true;
            dev_uptodate = false;
          }
        else
          {
            cerr << "operator= (BaseVector) : undefined vector" << endl;
          }
        return *this;
      }

    /* VFlatVector<double> fv(size, host_data); */
    /* fv = 1.0*v2; */
    FVDouble() = v2.FVDouble();

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

  /* size_t UnifiedVector :: Size () const throw() */
  /* { */
  /*   return size; */
  /* } */

  /*
   *  TODO:
   *    avoid cpy between device and host
   *    maybe leave the update to the user?
   * */
  const double & UnifiedVector :: operator [] (const int ind) const
  {
    cerr << "UnifiedVector operator[]" << endl;
    UpdateHost(); 
    return host_data[ind];
  }

  double & UnifiedVector :: operator [] (const int ind)
  {
    cerr << "UnifiedVector operator[] const" << endl;
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
    cerr << "Unified::Scale" << endl;
    UpdateDevice();
    cublasDscal (Get_CuBlas_Handle(), size, &scal, dev_data, 1);
    host_uptodate = false;
    return *this;
  }

  BaseVector & UnifiedVector :: SetScalar (double scal)
  {
    cerr << "uvec setscalar" << endl;
    (*this) = scal;
    return *this;
  }
  
  BaseVector & UnifiedVector :: Set (double scal, const BaseVector & v)
  {
    cerr << "uvec set" << endl;
    (*this) = 0.0;
    Add (scal, v);
    return *this;
  }
  
  
  BaseVector & UnifiedVector :: Add (double scal, const BaseVector & v)
  {
    const UnifiedVector * v2 = dynamic_cast<const UnifiedVector*> (&v);

    if (v2)
      {
        cerr << "uvec add uvec" << endl;
        UpdateDevice();
        v2->UpdateDevice();

        cublasDaxpy (Get_CuBlas_Handle(), 
                           size, &scal, v2->dev_data, 1, dev_data, 1);
        host_uptodate = false;
      }
    else
      {
        cerr << "uvec add basevec" << endl;
        UpdateHost();
        VFlatVector<> (size, host_data) += scal * v;

        // does not work, since operator+=(FlatVector<double>, VVecExpr<...>) not defined
        /* FVDouble() += scal * v; */
        dev_uptodate = false;
      }

    return *this;
  }
  
  double UnifiedVector :: InnerProduct (const BaseVector & v2) const
  {
    static Timer tdot("CUDA InnerProduct");
    RegionTimer reg(tdot);
    // cout << "Inner Prod" << endl;

    const UnifiedVector * uv2 = dynamic_cast<const UnifiedVector*> (&v2);
    if (uv2)
    {
      cerr << "uvec ip uvec" << endl;
      UpdateDevice();
      uv2->UpdateDevice();
      
      double res;
      cublasDdot (Get_CuBlas_Handle(), 
                        size, dev_data, 1, uv2->dev_data, 1, &res);
      return res;
    }
    cerr << "uvec ip basevec" << endl;

    FlatVector<> fv = FVDouble();
    FlatVector<> fv2 = v2.FVDouble();
    return ngbla::InnerProduct (fv, fv2);
  }


  ostream & UnifiedVector :: Print (ostream & ost) const
  {
    cout << "output unified vector of size " << size;
    cout << ", host = " << host_uptodate << ", dev = " << dev_uptodate << endl;
    if (!host_uptodate)
    {
      if (dev_uptodate)
      {
        cout << "host not up-to-data. printing device data" << endl;
        Vector<double> tmp(size);
        cudaMemcpy(tmp.Data(), dev_data, size * sizeof(double), cudaMemcpyDeviceToHost);
        ost << tmp << endl;
      }
      else
      {
        cout << "undefined vector" << endl;
      }
    }
    else
    {
      ost << FVDouble();
    }
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
    if (host_uptodate || !dev_uptodate) return;
    /* if (!dev_uptodate) cout << "ERROR UnifiedVector::UpdateHost non is uptodate" << endl; */
    cudaMemcpy (host_data, dev_data, sizeof(double)*size, cudaMemcpyDeviceToHost);    
    host_uptodate = true;
  }

  void UnifiedVector :: UpdateDevice () const
  {
    if (dev_uptodate || !host_uptodate) return;
    /* if (!host_uptodate) cout << "ERROR UnifiedVector::UpdateDevice non is uptodate" << endl; */
    cout << "Host2Device copy !!!!!!!!!!!!!!" << endl;
    cudaMemcpy (dev_data, host_data, sizeof(double)*size, cudaMemcpyHostToDevice);
    dev_uptodate = true;
  }
  
  FlatVector<double> UnifiedVector :: FVDouble () const
  {
    UpdateHost();
    dev_uptodate = false;
    return FlatVector<> (size, host_data);
  }
  
  FlatVector<Complex> UnifiedVector :: FVComplex () const
  {
    throw Exception ("unified complex not yet supported");
  }
    
  void * UnifiedVector :: Memory() const throw()
  { 
    cerr << "access2memory" << endl;
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

    cusparseCreateCsr(&descr, height, width, nze,
                      dev_ind, dev_col, dev_val,
                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
                      CUDA_R_64F);
  }

  DevSparseMatrix :: ~DevSparseMatrix ()
  {
    cusparseDestroySpMat(descr);
    cudaFree(dev_ind);
    cudaFree(dev_col);
    cudaFree(dev_val);
  }
  
  void DevSparseMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer tprel("CUDA PREL");
    static Timer tbuffer("CUDA BUFFER");
    static Timer tmv("CUDA MV");
    static Timer tbufferallocate("CUDA BUFFER ALLOCATE");

    tprel.Start();
    cout << "device mult sparse" << endl;
    cout << "vec0: " << typeid(x).name() << endl;
    cout << "vec1: " << typeid(y).name() << endl;
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
    cout << "device multadd sparse matrix" << endl;
    cout << "vec1: " << typeid(x).name() << endl;
    cout << "vec2: " << typeid(y).name() << endl;
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

  /* // y += s L * x */
  /* virtual void MultAdd1 (const BaseVector & x, BaseVector & y, */ 
  /*                        const BitArray * ainner = NULL, */
  /*                        const Array<int> * acluster = NULL) const */
  /* { */
  /*   throw Exception("not implemented yet."); */
  /* } */

  /* // y += s (D + L^T) * x) */
  /* void MultAdd2 (double s, const BaseVector & x, BaseVector & y, */
  /*                        const BitArray * ainner = NULL, */
  /*                        const Array<int> * acluster = NULL) const */
  /* { */
  /*   throw Exception("not implemented yet."); */
  /* } */



  /* DevJacobiPreconditioner :: DevJacobiPreconditioner (const SparseMatrix<double> & mat, */
  /*                 const BitArray & freedofs) */
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
  /*       temp_cols[i] = i; */
  /*       if (freedofs.Test(i)) */
  /*         temp_vals[i] = 1.0 / mat(i,i); */
  /*       else */
  /*         temp_vals[i] = 0.0; */
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
    int k = width;
    int n = b.Width();
    
    double alpha = 1;
    double beta = 0;

    /* throw Exception ("causes errors (prbly due to col-major order)") */

    cerr << "dense device matrix prod " << m << "x" << n << " and " << n << "x" << k << endl;
    cerr << b.Height() << " " << b.Width() << endl;
    cublasDgemm (Get_CuBlas_Handle(), CUBLAS_OP_N, CUBLAS_OP_N, m, n, k,
                 &alpha, dev_data, m, b.DevData(), k, &beta, c.DevData(), m);
  }

  void DevDMatrix :: MultAdd (const DevDMatrix& b, DevDMatrix& c)
  {
    throw Exception("not implemented yet.");
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


  shared_ptr<BaseMatrix> CreateDevMatrix (BaseMatrix & mat)
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
  shared_ptr<BaseMatrix> CreateDevMatrix (Matrix<> & mat)
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

  DevJacobiPrecond :: DevJacobiPrecond (const SparseMatrix<double> & amat, 
    shared_ptr<BitArray> ainner, bool use_par)
      : inner(ainner)
  {
    cerr << "Creating Device Jacobi Smoother of size " << height << endl;

    height = amat.Height();
    width = amat.Height();
    nze = 0;
    Array<int> tmp_ind(height+1);
    Array<int> tmp_col(height);
    Array<double> tmp_val(height);

    tmp_ind[0] = 0;
    for (int i=0; i<height; i++)
    {
      if (!inner || inner->Test(i))
      {
        tmp_col[nze] = i;
        tmp_ind[i+1] = tmp_ind[i]+1;
        tmp_val[nze] = 1 / amat(i,i);
        nze++;
      }
      else
      {
        tmp_ind[i+1] = tmp_ind[i];
      }
    }
    
    cerr << "malloc." << endl;
    cudaMalloc((void**) &dev_ind, (height+1) * sizeof(int));
    cudaMalloc((void**) &dev_col, nze * sizeof(int));
    cudaMalloc((void**) &dev_val, nze * sizeof(double));

    cerr << "copying" << endl;
    cudaMemcpy(dev_ind, tmp_ind.begin(), (height + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_col, tmp_col.begin(), nze * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_val, tmp_val.begin(), nze * sizeof(double), cudaMemcpyHostToDevice);

    cerr << "creating matrix" << endl;
    cusparseCreateCsr(&descr, height, height, nze,
                      dev_ind, dev_col, dev_val,
                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
                      CUDA_R_64F);
    cerr << "finished." << endl;
  }

  DevJacobiPrecond :: ~DevJacobiPrecond ()
  {
    cerr << "DevJacobiPrecond dtor" << endl;
  }

}
