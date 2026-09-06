/*********************************************************************/
/* File:   cuda_linalg.cpp                                           */
/* Author: Joachim Schoeberl, Matthias Hochsteger                    */
/* Date:   11. Aug. 2014                                             */
/*********************************************************************/

#include <la.hpp>
#include "cuda_linalg.hpp"

namespace ngla
{
  extern void InitSparseCholesky();

  
  cublasHandle_t Get_CuBlas_Handle ()
  {
    // static Timer tblashandle("CUDA create cublas handle");
    // RegionTimer reg(tblashandle);

    static cublasHandle_t handle;
    static bool first_call = true;

    if (first_call)
      {
        first_call = false;
        cublasCreate_v2 (&handle);

      }
    return handle;
  }

  // Workspace pointer — only allocated when graphs are used
  static void* cublas_workspace = nullptr;
  static const size_t cublas_workspace_size = 32 * 1024 * 1024; // 32 MiB (H100)

  // Call before graph capture to ensure workspace is allocated
  void EnsureCuBlasWorkspace()
  {
    if (!cublas_workspace)
      cudaMalloc(&cublas_workspace, cublas_workspace_size);
    cublasSetWorkspace(Get_CuBlas_Handle(), cublas_workspace, cublas_workspace_size);
  }

  // Wrapper: set stream AND re-apply workspace if already allocated
  // (cublasSetStream resets workspace to default pool)
  void SetCuBlasStream(cudaStream_t stream)
  {
    cublasSetStream(Get_CuBlas_Handle(), stream);
    if (cublas_workspace)
      cublasSetWorkspace(Get_CuBlas_Handle(), cublas_workspace, cublas_workspace_size);
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
        std::cerr << "[cusparse] handle created" << std::endl;
      }
    return handle;
  }

  void InitCuLinalg()
  {
    cout << "Initializing cublas and cusparse." << endl;

    Get_CuBlas_Handle();
    Get_CuSparse_Handle();

    cusparseSetStream(Get_CuSparse_Handle(), ngs_cuda::ngs_cuda_stream);
    std::cerr << "[InitCuLinalg] cusparseSetStream bound to ngs_cuda_stream" << std::endl;

    ngs_cuda::CudaGraph::stream_change_callback = [](cudaStream_t s) {
      cusparseSetStream(Get_CuSparse_Handle(), s);
    };

    // with lazy module loading the first kernel launch of this library
    // finalizes the whole fatbin (~300ms) - pay that at import, not
    // inside the user's first operator application
    DeviceParallelFor (1, [] DEVICE_LAMBDA (size_t) { });
    cudaDeviceSynchronize();

    std::cerr << "[InitCuLinalg] callback wired, registering creators..." << std::endl;
    BaseVector::RegisterDeviceVectorCreator(typeid(S_BaseVectorPtr<double>),
                                            [] (const BaseVector & vec, bool unified) -> shared_ptr<BaseVector>
                                            {
                                              return make_shared<UnifiedVector>(vec);
                                            });
    BaseVector::RegisterDeviceVectorCreator(typeid(VVector<double>),
                                            [] (const BaseVector & vec, bool unified) -> shared_ptr<BaseVector>
                                            {
                                              return make_shared<UnifiedVector>(vec);
                                            });
    
    BaseMatrix::RegisterDeviceMatrixCreator(typeid(SparseMatrix<double>),
                                            [] (const BaseMatrix & mat) -> shared_ptr<BaseMatrix>
                                            {
                                              auto & sparse_mat = dynamic_cast<const SparseMatrix<double>&>(mat);
                                              return make_shared<DevSparseMatrix>(sparse_mat);
                                            });



    
    
  }


  /******************** DevMatrix ********************/

  shared_ptr<BaseMatrix> CreateDevMatrix (BaseMatrix & mat)
  {
    if (auto res = mat.CreateDeviceMatrix())
      return res;
    else
      throw Exception(string("matrix type not supported: ") + typeid(mat).name());
  }


  /******************** DevSparseMatrix ********************/

  DevSparseMatrix :: DevSparseMatrix (const SparseMatrix<double> & mat)
  {
    height = mat.Height();
    width = mat.Width();
    nze = mat.NZE();

    cout << IM(7) << "DevSparseMatrix" << endl
         << " height = " << height << ", width = " << width << ", nze = " << nze << endl;
    
    Array<int> temp_ind (height+1); 
    for (int i = 0; i <= height; i++) temp_ind[i] = mat.First(i); // conversion to 32-bit integer

    cudaMalloc ((void**)&dev_ind, (mat.Height()+1) * sizeof(int));
    cudaMalloc ((void**)&dev_col, (mat.NZE()) * sizeof(int));
    cudaMalloc ((void**)&dev_val, (mat.NZE()) * sizeof(double));
    
    cudaMemcpy (dev_ind, temp_ind.Data(), (mat.Height()+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy (dev_col, mat.GetRowIndices(0).Data(), mat.NZE()*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy (dev_val, mat.GetRowValues(0).Data(), mat.NZE()*sizeof(double), cudaMemcpyHostToDevice);

    cusparseCreateCsr(&descr, height, width, nze,
                      dev_ind, dev_col, dev_val,
                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
                      CUDA_R_64F);

    // pre-compute buffer size, preprocess for graph capture
    { double alpha=1, beta=0;
//      std::cerr << "[preprocess] entering block, handle=" << Get_CuSparse_Handle() << std::endl;
      double *raw_x = nullptr, *raw_y = nullptr;
      auto err1 = cudaMalloc(&raw_x, width  * sizeof(double));
      auto err2 = cudaMalloc(&raw_y, height * sizeof(double));
//      std::cerr << "[preprocess] raw_x=" << raw_x << " err1=" << err1 << " raw_y=" << raw_y << " err2=" << err2 << std::endl;
      cusparseDnVecDescr_t dx, dy;
      cusparseCreateDnVec(&dx, width,  raw_x, CUDA_R_64F);
      cusparseCreateDnVec(&dy, height, raw_y, CUDA_R_64F);
//      std::cerr << "[preprocess] calling bufferSize" << std::endl;
      cusparseSpMV_bufferSize(Get_CuSparse_Handle(), CUSPARSE_OPERATION_NON_TRANSPOSE,
          &alpha, descr, dx, &beta, dy, CUDA_R_64F,
          CUSPARSE_SPMV_ALG_DEFAULT, &spmv_bufferSize);
//      std::cerr << "[preprocess] bufferSize=" << spmv_bufferSize << " calling cudaMalloc" << std::endl;
      cudaMalloc(&spmv_buffer, spmv_bufferSize);
//      std::cerr << "[preprocess] buffer allocated OK, skipping preprocess" << std::endl;
      cusparseDestroyDnVec(dx); cusparseDestroyDnVec(dy);
      cudaFree(raw_x); cudaFree(raw_y); }
  }


  DevSparseMatrix :: ~DevSparseMatrix ()
  {
    cusparseDestroySpMat(descr);
    cudaFree(spmv_buffer);
    cudaFree(dev_ind);
    cudaFree(dev_col);
    cudaFree(dev_val);
  }


  void DevSparseMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer tmv("DevSparseMatrix :: MultAdd");
    CudaRegionTimer rt(tmv);
    // RegionTimer reg(tmv);

    DeviceVectorWrapper<double> ux(x);
    DeviceVectorWrapper<double> uy(y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    double alpha= s;
    double beta = 1;

    cusparseDnVecDescr_t descr_x, descr_y;
    cusparseCreateDnVec (&descr_x, ux.Size(), DevPtr(ux), CUDA_R_64F);
    cusparseCreateDnVec (&descr_y, uy.Size(), DevPtr(uy), CUDA_R_64F);

    { cudaStreamCaptureStatus cap_status;
      cudaStreamIsCapturing(ngs_cuda::ngs_cuda_stream, &cap_status);
//      std::cerr << "[SpMV] capture status before SpMV: " << cap_status << " (1=capturing)" << std::endl; }
      }
    cusparseSetStream(Get_CuSparse_Handle(), ngs_cuda::ngs_cuda_stream);
    cusparseSpMV(Get_CuSparse_Handle(),
                 CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, descr,
                 descr_x, &beta, descr_y, CUDA_R_64F,
                 CUSPARSE_SPMV_ALG_DEFAULT, spmv_buffer);

    cusparseDestroyDnVec(descr_x);
    cusparseDestroyDnVec(descr_y);
    uy.InvalidateHost();
  }


  void DevSparseMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer tmv("DevSparseMatrix :: MultTransAdd");
    CudaRegionTimer reg(tmv);

    DeviceVectorWrapper<double> ux(x);
    DeviceVectorWrapper<double> uy(y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    double alpha= s;
    double beta = 1;

    size_t bufferSize = 0;
    void* dBuffer = NULL;

    cusparseDnVecDescr_t descr_x, descr_y;
    cusparseCreateDnVec (&descr_x, ux.Size(), DevPtr(ux), CUDA_R_64F);
    cusparseCreateDnVec (&descr_y, uy.Size(), DevPtr(uy), CUDA_R_64F);

    cusparseSpMV_bufferSize(Get_CuSparse_Handle(), CUSPARSE_OPERATION_TRANSPOSE,
                            &alpha, descr, descr_x, &beta, descr_y, CUDA_R_64F,
                            CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);
    cudaMalloc(&dBuffer, bufferSize);

    cusparseSpMV(Get_CuSparse_Handle(), 
                 CUSPARSE_OPERATION_TRANSPOSE, &alpha, descr,
                 descr_x, &beta, descr_y, CUDA_R_64F,
                 CUSPARSE_SPMV_ALG_DEFAULT, dBuffer);

    cudaFree(dBuffer);

    cusparseDestroyDnVec(descr_x);
    cusparseDestroyDnVec(descr_y);

    uy.InvalidateHost();
  }


  

  void DevDiagonalMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    DeviceVectorWrapper<double> ux(x);
    DeviceVectorWrapper<double> uy(y);
    
    /*
    ux.UpdateDevice();
    uy.UpdateDevice();

    // MultDiagonal (diag.Size(), diag.DevData(), DevPtr(ux), DevPtr(uy));
    DeviceParallelFor
      (diag.Size(),
       [ddiag=diag.DevData(), dx=DevPtr(ux), dy=DevPtr(uy)] DEVICE_LAMBDA (auto tid)
           {
             dy[tid] = ddiag[tid]*dx[tid];
           });

    uy.InvalidateHost();
    */

    DeviceParallelFor
      (diag.Size(),
       [ddiag=diag.DevData(), dx=FVDevRO(ux), dy=FVDev(uy)] DEVICE_LAMBDA (auto tid)
           {
             dy(tid) = ddiag[tid]*dx(tid);
           });
  }
  
  void DevDiagonalMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    DeviceVectorWrapper<double> ux(x);
    DeviceVectorWrapper<double> uy(y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    // MultAddDiagonal (diag.Size(), s, diag.DevData(), DevPtr(ux), DevPtr(uy));
    DeviceParallelFor
      (diag.Size(),
       [ddiag=diag.DevData(), dx=DevPtr(ux), dy=DevPtr(uy), s] DEVICE_LAMBDA (auto tid)
           {
             dy[tid] += s*ddiag[tid]*dx[tid];
           });

    uy.InvalidateHost();    
  }
  
  bool synckernels = true;
}
