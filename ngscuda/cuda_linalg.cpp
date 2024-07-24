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

  void InitCuLinalg()
  {
    cout << "Initializing cublas and cusparse." << endl;

    Get_CuBlas_Handle();
    Get_CuSparse_Handle();

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

    BaseMatrix::RegisterDeviceMatrixCreator(typeid(JacobiPrecond<double>),
                                            [] (const BaseMatrix & mat) -> shared_ptr<BaseMatrix>
                                            {
                                              auto & Jacobimat = dynamic_cast<const JacobiPrecond<double>&>(mat);
                                              auto diagarray = Jacobimat.GetInverse();

                                              VVector<double> diag(diagarray.Size());
                                              auto fv = diag.FVDouble();
                                              for (size_t i = 0; i < fv.Size(); i++)
                                                fv[i] = diagarray[i];
                                              
                                              return make_shared<DevDiagonalMatrix>(diag);
                                            });

    BaseMatrix::RegisterDeviceMatrixCreator(typeid(DiagonalMatrix<double>),
                                            [] (const BaseMatrix & mat) -> shared_ptr<BaseMatrix>
                                            {
                                              auto & diagmat = dynamic_cast<const DiagonalMatrix<double>&>(mat);
                                              
                                              return make_shared<DevDiagonalMatrix>(diagmat.AsVector());
                                            });
    
    BaseMatrix::RegisterDeviceMatrixCreator(typeid(ConstantElementByElementMatrix<>),
                                            [] (const BaseMatrix & mat) -> shared_ptr<BaseMatrix>
                                            {
                                              auto & ebe_mat = dynamic_cast<const ConstantElementByElementMatrix<>&>(mat);
                                              return make_shared<DevConstantElementByElementMatrix>(ebe_mat);
                                            });

    BaseMatrix::RegisterDeviceMatrixCreator(typeid(BlockDiagonalMatrixSoA),
                                            [] (const BaseMatrix & mat) -> shared_ptr<BaseMatrix>
                                            {
                                              auto & bdm_mat = dynamic_cast<const BlockDiagonalMatrixSoA&>(mat);
                                              return make_shared<DevBlockDiagonalMatrixSoA>(bdm_mat);
                                            });

    
    BaseMatrix::RegisterDeviceMatrixCreator(typeid(EmbeddedMatrix),
                                            [] (const BaseMatrix & bmat) -> shared_ptr<BaseMatrix>
                                            {
                                              auto & mat = dynamic_cast<const EmbeddedMatrix&>(bmat);
                                              return make_shared<DevEmbeddedMatrix>(mat.Height(), mat.GetRange(),
                                                                                    mat.GetMatrix()->CreateDeviceMatrix());
                                            });
    
    BaseMatrix::RegisterDeviceMatrixCreator(typeid(EmbeddedTransposeMatrix),
                                            [] (const BaseMatrix & bmat) -> shared_ptr<BaseMatrix>
                                            {
                                              auto & mat = dynamic_cast<const EmbeddedTransposeMatrix&>(bmat);
                                              return make_shared<DevEmbeddedTransposeMatrix>(mat.Width(), mat.GetRange(),
                                                                                             mat.GetMatrix()->CreateDeviceMatrix());
                                            });

    BaseMatrix::RegisterDeviceMatrixCreator(typeid(Projector),
                                            [] (const BaseMatrix & bmat) -> shared_ptr<BaseMatrix>
                                            {
                                              auto & proj = dynamic_cast<const Projector&>(bmat);
                                              return make_shared<DevProjector>(proj);
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
  }


  DevSparseMatrix :: ~DevSparseMatrix ()
  {
    cusparseDestroySpMat(descr);
    cudaFree(dev_ind);
    cudaFree(dev_col);
    cudaFree(dev_val);
  }


  void DevSparseMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer tmv("DevSparseMatrix :: MultAdd");
    CudaRegionTimer rt(tmv);
    // RegionTimer reg(tmv);

    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    double alpha= s;
    double beta = 1;

    size_t bufferSize = 0;
    void* dBuffer = NULL;

    cusparseDnVecDescr_t descr_x, descr_y;
    cusparseCreateDnVec (&descr_x, ux.Size(), ux.DevData(), CUDA_R_64F);
    cusparseCreateDnVec (&descr_y, uy.Size(), uy.DevData(), CUDA_R_64F);

    cusparseSpMV_bufferSize(Get_CuSparse_Handle(), CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha, descr, descr_x, &beta, descr_y, CUDA_R_64F,
                            CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);
    cudaMalloc(&dBuffer, bufferSize);

    cusparseSpMV(Get_CuSparse_Handle(), 
                 CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, descr,
                 descr_x, &beta, descr_y, CUDA_R_64F,
                 CUSPARSE_SPMV_ALG_DEFAULT, dBuffer);

    cudaFree(dBuffer);

    cusparseDestroyDnVec(descr_x);
    cusparseDestroyDnVec(descr_y);
    uy.InvalidateHost();
  }


  void DevSparseMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer tmv("DevSparseMatrix :: MultTransAdd");
    CudaRegionTimer reg(tmv);

    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    double alpha= s;
    double beta = 1;

    size_t bufferSize = 0;
    void* dBuffer = NULL;

    cusparseDnVecDescr_t descr_x, descr_y;
    cusparseCreateDnVec (&descr_x, ux.Size(), ux.DevData(), CUDA_R_64F);
    cusparseCreateDnVec (&descr_y, uy.Size(), uy.DevData(), CUDA_R_64F);

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
    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    // MultDiagonal (diag.Size(), diag.DevData(), ux.DevData(), uy.DevData());
    DeviceParallelFor
      (diag.Size(),
       [ddiag=diag.DevData(), dx=ux.DevData(), dy=uy.DevData()] DEVICE_LAMBDA (auto tid)
           {
             dy[tid] = ddiag[tid]*dx[tid];
           });

    uy.InvalidateHost();    
  }
  
  void DevDiagonalMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    // MultAddDiagonal (diag.Size(), s, diag.DevData(), ux.DevData(), uy.DevData());
    DeviceParallelFor
      (diag.Size(),
       [ddiag=diag.DevData(), dx=ux.DevData(), dy=uy.DevData(), s] DEVICE_LAMBDA (auto tid)
           {
             dy[tid] += s*ddiag[tid]*dx[tid];
           });

    uy.InvalidateHost();    
  }
  
  /******************** DevConstantEBEMatrix ********************/

  bool synckernels = true;
  
  DevConstantElementByElementMatrix ::
  DevConstantElementByElementMatrix (const ConstantElementByElementMatrix<> & mat)
    : h(mat.Height()), w(mat.Width()),
      devmat(mat.GetMatrix()),
      rowdnums(mat.GetRowDNums()), coldnums(mat.GetColDNums()),
      row_coloring(mat.GetRowColoring()), col_coloring(mat.GetColColoring()),
      numblocks(mat.GetRowDNums().Size())
  {
    disjoint_rows = (mat.GetRowColoring().Size() == 0);
    disjoint_cols = (mat.GetColColoring().Size() == 0);
  }

  
  void DevConstantElementByElementMatrix ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DevConstantEBEMatrix::MultAdd"); RegionTimer reg(t);
    static Timer tmult("DevConstantEBEMatrix::MultAdd - mult");
    static Timer tcopyin("DevConstantEBEMatrix::MultAdd - copyin");
    static Timer tcopyout("DevConstantEBEMatrix::MultAdd - copyout");
    
    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    if (synckernels)
     cudaDeviceSynchronize();
    
    // if (disjoint_cols)
    if (true)
      {
        DevStackArray<double> dev_hx(numblocks*devmat.Width());
        DevStackArray<double> dev_hy(numblocks*devmat.Height());

        tcopyin.Start();
        // ConstEBEKernelCopyIn (numblocks, devmat.Width(), rowdnums.DevData(), ux.DevData(), dev_hx.DevData());
	DeviceParallelFor
          (numblocks*devmat.Width(),
           [locx=dev_hx.DevData(), globx=ux.DevData(), idx=rowdnums.DevData()] DEVICE_LAMBDA (auto tid)
           {
             locx[tid] = globx[idx[tid]];
           });
        if (synckernels) cudaDeviceSynchronize();
        tcopyin.Stop();
        
        // dev_hy = dev_hx * Trans(mat)
        tmult.Start();
        FlatMatrix<Dev<double>> matx(numblocks, devmat.Width(), dev_hx.Data());
        FlatMatrix<Dev<double>> maty(numblocks, devmat.Height(), dev_hy.Data());
        // MultMatMat (matx, Trans(devmat), maty, s, 0);
        maty = s * matx * Trans(devmat);
        if (synckernels) cudaDeviceSynchronize();
        tmult.Stop();
        
        tcopyout.Start();        
        // ConstEBEKernelCopyOut (numblocks, devmat.Height(), coldnums.DevData(), dev_hy.DevData(), uy.DevData());
        DeviceParallelFor
          (numblocks*devmat.Height(),
           [globy=uy.DevData(), locy=dev_hy.DevData(), idx=coldnums.DevData() ] DEVICE_LAMBDA (auto tid)
           {
             atomicAdd((double*)globy+idx[tid], locy[tid]);
           });
 
        if (synckernels) cudaDeviceSynchronize();
        tcopyout.Stop();
      }
    else
      {
        for (auto c : col_coloring)
          {
            DevStackArray<double> dev_hx(c.Size()*devmat.Width());
            DevStackArray<double> dev_hy(c.Size()*devmat.Height());

            tcopyin.Start();            
            ConstEBEKernelCopyInIdx (c.Size(), (int*)c.Data(), devmat.Width(), rowdnums.DevData(), (double*)ux.DevData(), dev_hx.DevData());
            if (synckernels) cudaDeviceSynchronize();            
            tcopyin.Stop();
            
            // dev_hy = dev_hx * Trans(mat)
            tmult.Start();
            FlatMatrix<Dev<double>> matx(c.Size(), devmat.Width(), dev_hx.Data());
            FlatMatrix<Dev<double>> maty(c.Size(), devmat.Height(), dev_hy.Data());
            // MultMatMat (matx, Trans(devmat), maty, s, 0);
            maty = s * matx * Trans(devmat); 
            if (synckernels) cudaDeviceSynchronize();            
            tmult.Stop();

            tcopyout.Start();
            ConstEBEKernelCopyOutIdx (c.Size(), (int*)c.Data(), devmat.Height(), coldnums.DevData(), dev_hy.DevData(), (double*)uy.DevData());
            if (synckernels) cudaDeviceSynchronize();            
            tcopyout.Stop();
          }
      }

    uy.InvalidateHost();    
  }
  

  void DevConstantElementByElementMatrix ::
  MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DevConstantEBEMatrix::MultTransAdd"); RegionTimer reg(t);
    
    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);
    
    ux.UpdateDevice();
    uy.UpdateDevice();
    
    auto hm = devmat.Height();
    auto wm = devmat.Width();
    // if (disjoint_rows)
    if (true)      
      {
        DevStackArray<double> dev_hx(numblocks*hm);
        DevStackArray<double> dev_hy(numblocks*wm);

        // ConstEBEKernelCopyIn (numblocks, hm, coldnums.DevData(), (double*)ux.DevData(), dev_hx.DevData());
        DeviceParallelFor
          (numblocks*hm,
           [locx=dev_hx.DevData(), globx=ux.DevData(), idx=coldnums.DevData()] DEVICE_LAMBDA (auto tid)
           {
             locx[tid] = globx[idx[tid]];
           });
        
        // dev_hy = dev_hx * mat
        
        FlatMatrix<Dev<double>> matx(numblocks, hm, dev_hx.Data());
        FlatMatrix<Dev<double>> maty(numblocks, wm, dev_hy.Data());
        MultMatMat (matx, devmat, maty, s, 0);
        
        // ConstEBEKernelCopyOut (numblocks, wm, rowdnums.DevData(), dev_hy.DevData(), (double*)uy.DevData());
        DeviceParallelFor
          (numblocks*wm,
           [globy=uy.DevData(), locy=dev_hy.DevData(), idx=rowdnums.DevData()] DEVICE_LAMBDA (auto tid)
           {
             atomicAdd((double*)globy+idx[tid], locy[tid]);
           });
        
      }
    else
      {
        for (auto c : row_coloring)
          {
            DevStackArray<double> dev_hx(c.Size()*hm);
            DevStackArray<double> dev_hy(c.Size()*wm);

            ConstEBEKernelCopyInIdx (c.Size(), (int*)c.Data(), hm, coldnums.DevData(), (double*)ux.DevData(), dev_hx.DevData());
            // dev_hy = dev_hx * mat

            FlatMatrix<Dev<double>> matx(c.Size(), hm, dev_hx.Data());
            FlatMatrix<Dev<double>> maty(c.Size(), wm, dev_hy.Data());
            MultMatMat (matx, devmat, maty, s, 0);
           
            ConstEBEKernelCopyOutIdx (c.Size(), (int*)c.Data(), wm, rowdnums.DevData(), dev_hy.DevData(), (double*)uy.DevData());
          }
      }
    
    if (synckernels) cudaDeviceSynchronize();
    uy.InvalidateHost();    
  }
  



  
  
  DevBlockDiagonalMatrixSoA ::
  DevBlockDiagonalMatrixSoA (const BlockDiagonalMatrixSoA & mat)
    : nonzero(mat.GetNonZeroPattern())
  {
    FlatTensor<3> blockdiag = mat.GetBlockDiag ();
    // tie(dimy, dimx, blocks) = blockdiag.Shape();

    dimy = blockdiag.GetSize();
    dimx = blockdiag.GetSubTensor().GetSize();
    blocks = blockdiag.GetSubTensor().GetSubTensor().GetSize();
    
    auto err = cudaMalloc((void**)&dev_data, dimx*dimy*blocks*sizeof(double));
    if (err != 0)
      throw Exception("DevBlockDiagonalMatrixSoA allocation error, ec="+ToString(err));

    cudaMemcpy (dev_data, blockdiag.Data(), sizeof(double)*dimx*dimy*blocks, cudaMemcpyHostToDevice);


    Array<int> nonzeroinds, nonzeroinds_trans;
    for (int i = 0; i < dimy; i++)
      for (int j = 0; j < dimx; j++)
        if (nonzero(i,j) != 0)
          {
            nonzeroinds.Append(i*dimx+j);
            nonzeroinds.Append(j);
            nonzeroinds.Append(i);

            nonzeroinds_trans.Append(i*dimx+j);
            nonzeroinds_trans.Append(i);
            nonzeroinds_trans.Append(j);
        
        }

    indices = nonzeroinds;
    indices_trans = nonzeroinds_trans;
}
  
  void DevBlockDiagonalMatrixSoA :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DevBlockDiagonalMatrixSoA::MultAdd"); RegionTimer reg(t);
    
    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);
    ux.UpdateDevice();
    uy.UpdateDevice();
/*
    for (int i = 0; i < dimy; i++)
      for (int j = 0; j < dimx; j++)
        if (nonzero(i,j) != 0)
          DevBlockDiagonalMatrixSoAMultAddVecs (s, blocks, dev_data + blocks*(i*dimx+j), ux.DevData()+blocks*j, uy.DevData()+blocks*i);
*/

    FlatMatrix<Dev<double>> a(dimx*dimy, blocks, (Dev<double>*)dev_data);
    FlatMatrix<Dev<double>> b(dimx, blocks,  (Dev<double>*)ux.DevData());
    FlatMatrix<Dev<double>> res(dimy, blocks,  (Dev<double>*)uy.DevData());
    // DevBlockDiagonalMatrixSoAMultAddVecs (s, indices, a, b, res);

    DeviceParallelFor
      (res.Width(),
       [a,b,res,inds=FlatArray(indices),s] DEVICE_LAMBDA (auto i)
       {
         for (int j = 0; j < inds.Size(); j+=3)
           {
             int rowa = inds[j];
             int rowb = inds[j+1];
             int rowres = inds[j+2];
             res(rowres,i) += s * a(rowa,i) * b(rowb,i);
           }
       });
    
    if (synckernels) cudaDeviceSynchronize();

    uy.InvalidateHost();
  }

  void DevBlockDiagonalMatrixSoA :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DevBlockDiagonalMatrixSoA::MultTransAdd"); RegionTimer reg(t);
    
    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);
    ux.UpdateDevice();
    uy.UpdateDevice();
    
      /*
    for (int i = 0; i < dimy; i++)
      for (int j = 0; j < dimx; j++)
        if (nonzero(i,j) != 0)
          DevBlockDiagonalMatrixSoAMultAddVecs (s, blocks, dev_data + blocks*(i*dimx+j), ux.DevData()+blocks*i, uy.DevData()+blocks*j);
*/
      
    FlatMatrix<Dev<double>> a(dimx*dimy, blocks, (Dev<double>*)dev_data);
    FlatMatrix<Dev<double>> b(dimy, blocks,  (Dev<double>*)ux.DevData());
    FlatMatrix<Dev<double>> res(dimx, blocks,  (Dev<double>*)uy.DevData());
    // DevBlockDiagonalMatrixSoAMultAddVecs (s, indices_trans, a, b, res);

    DeviceParallelFor
      (res.Width(),
       [a,b,res,inds=FlatArray(indices),s] DEVICE_LAMBDA (auto i)
       {
         for (int j = 0; j < inds.Size(); j+=3)
           {
             int rowa = inds[j];
             int rowb = inds[j+2];
             int rowres = inds[j+1];
             res(rowres,i) += s * a(rowa,i) * b(rowb,i);
           }
       });

    
    if (synckernels) cudaDeviceSynchronize();
      
    uy.InvalidateHost();
  }


  

  

  void DevProjector :: Mult (const BaseVector & x, BaseVector & y) const
  {
    y = x;
    Project (y);
  }

  void DevProjector :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DevProjector::MultAdd"); RegionTimer reg(t);
    if (x.EntrySize() != 1)
      throw Exception("DevProjector :: MultAdd not implemented for EntrySize > 1");

    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);
    ux.UpdateDevice();
    uy.UpdateDevice();

    // DevProjectorMultAdd (s, bits->Size(), ux.DevData(), uy.DevData(), bits->Data(), keep_values);
    if (keep_values)
      DeviceParallelFor
        (bits->Size(),
         [s, x=ux.DevData(), y=uy.DevData(), bits=bits->Data()] DEVICE_LAMBDA (auto i)
         {
           unsigned char mask = (char(1) << (i % CHAR_BIT));
           unsigned int addr = i / CHAR_BIT;
           
           if (bits[addr] & mask)
             y[i] += s * x[i];
         });
    else
      DeviceParallelFor
        (bits->Size(),
         [s, x=ux.DevData(), y=uy.DevData(), bits=bits->Data()] DEVICE_LAMBDA (auto i)
         {
           unsigned char mask = (char(1) << (i % CHAR_BIT));
           unsigned int addr = i / CHAR_BIT;
           
           if (! (bits[addr] & mask) )
             y[i] += s * x[i];
         });

    uy.InvalidateHost();
  }

  void DevProjector :: Project (BaseVector & x) const
  {
    static Timer t("DevProjector::Project"); RegionTimer reg(t);
    if (x.EntrySize() != 1)
      throw Exception("DevProjector :: Project not implemented for EntrySize > 1");

    UnifiedVectorWrapper ux(x);

    ux.UpdateDevice();

    // DevProjectorProject (bits->Size(), (double*)ux.DevData(), bits->Data(), keep_values);
    if (keep_values)
      DeviceParallelFor
        (bits->Size(),
         [x=ux.DevData(), bits=bits->Data()] DEVICE_LAMBDA (auto i)
         {
           unsigned char mask = (char(1) << (i % CHAR_BIT));
           unsigned int addr = i / CHAR_BIT;
           
           if (! (bits[addr] & mask) )
             x[i] = 0.0;
         });
    else
      DeviceParallelFor
        (bits->Size(),
         [x=ux.DevData(), bits=bits->Data()] DEVICE_LAMBDA (auto i)
         {
           unsigned char mask = (char(1) << (i % CHAR_BIT));
           unsigned int addr = i / CHAR_BIT;
           
           if (bits[addr] & mask)
             x[i] = 0;
         });
      
      
    ux.InvalidateHost();
  }

    
}
