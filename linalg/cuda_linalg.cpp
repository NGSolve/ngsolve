#ifdef CUDA 

/*********************************************************************/
/* File:   cuda_linalg.cpp                                           */
/* Author: Joachim Schoeberl, Matthias Hochsteger                    */
/* Date:   11. Aug. 2014                                             */
/*********************************************************************/


#include <la.hpp>
#include <cublas_v2.h>
#include <cusparse.h>

extern void SetScalar (double val, int n, double * dev_ptr);



namespace ngla
{

  /*
  cublasHandle_t handle;
  cusparseHandle_t cusparseHandle;
  */

  cublasHandle_t Get_CuBlas_Handle ()
  {
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
    size = asize;
    host_data = new double[asize];
    cudaMalloc(&dev_data, size*sizeof(double));
    host_uptodate = false;
    dev_uptodate = false;

    (*this) = 0.0;
  }

  BaseVector & UnifiedVector :: operator= (double d)
  {
    for (int i = 0; i < size; i++) host_data[i] = d;
    ::SetScalar (d, size, dev_data);
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

  BaseVector & UnifiedVector :: operator= (BaseVector & v2)
  {
    UnifiedVector * uv2 = dynamic_cast<UnifiedVector*> (&v2);
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
            VFlatVector<double> fv(size, host_data);
            fv = 1.0*v2;
            host_uptodate = true;
            dev_uptodate = false;
          }
        else
          {
            cerr << "operator= (BaseVector) : undefined vector" << endl;
          }
        return *this;
      }

    VFlatVector<double> fv(size, host_data);
    fv = 1.0*v2;

    host_uptodate = true;
    dev_uptodate = false;
    return *this;
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
  
  double UnifiedVector :: InnerProduct (const BaseVector & v2) const
  {
    cout << "Inner Prod" << endl;

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
    cout << "output unified vector, host = " << host_uptodate << ", dev = " << dev_uptodate << endl;
    ost << FVDouble();
    return ost;
  }
  
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
    UpdateHost();
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
    return make_shared<UnifiedVector> (size);
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

    descr = new cusparseMatDescr_t;
    cusparseCreateMatDescr (descr);

    cusparseSetMatType(*descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(*descr, CUSPARSE_INDEX_BASE_ZERO);

    cout << "create device sparse matrix, n = " << height << ", nze = " << nze << endl;
    
    Array<int> temp_ind (mat.Height()+1); 
    for (int i = 0; i <= mat.Height(); i++) temp_ind[i] = mat.First(i); // conversion to 32-bit integer

    cudaMalloc ((void**)&dev_ind, (mat.Height()+1) * sizeof(int));
    cudaMalloc ((void**)&dev_col, (mat.NZE()) * sizeof(int));
    cudaMalloc ((void**)&dev_val, (mat.NZE()) * sizeof(double));

    cudaMemcpy (dev_ind, &temp_ind[0], (mat.Height()+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy (dev_col, &mat.GetRowIndices(0)[0], mat.NZE()*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy (dev_val, &mat.GetRowValues(0)[0], mat.NZE()*sizeof(double), cudaMemcpyHostToDevice);
  }
  
  void DevSparseMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    cout << "device mult sparse" << endl;
    const UnifiedVector & ux = dynamic_cast<const UnifiedVector&> (x);
    UnifiedVector & uy = dynamic_cast<UnifiedVector&> (y);

    ux.UpdateDevice();
    uy = 0.0;
    uy.UpdateDevice();

    double alpha= 1;
    double beta = 0;
    cusparseDcsrmv (Get_CuSparse_Handle(), 
                    CUSPARSE_OPERATION_NON_TRANSPOSE, height, width, nze, 
		    &alpha, *descr, 
		    dev_val, dev_ind, dev_col, 
		    ux.dev_data, &beta, uy.dev_data);

    uy.host_uptodate = false;
    cout << "mult complete" << endl;
  }


  void DevSparseMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    cout << "device multadd sparse matrix" << endl;
    const UnifiedVector & ux = dynamic_cast<const UnifiedVector&> (x);
    UnifiedVector & uy = dynamic_cast<UnifiedVector&> (y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    double alpha= 1;
    double beta = 1;
    cusparseDcsrmv (Get_CuSparse_Handle(), 
                    CUSPARSE_OPERATION_NON_TRANSPOSE, height, width, nze, 
		    &alpha, *descr, 
		    dev_val, dev_ind, dev_col, 
		    ux.dev_data, &beta, uy.dev_data);

    uy.host_uptodate = false;

    cout << "mult complete" << endl;
  }
  







  DevJacobiPreconditioner :: DevJacobiPreconditioner (const SparseMatrix<double> & mat,
						      const BitArray & freedofs)
  {
    height = mat.Height();
    width = mat.Height();
    nze = mat.Height();

    descr = new cusparseMatDescr_t;
    cusparseCreateMatDescr (descr);

    cusparseSetMatType(*descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(*descr, CUSPARSE_INDEX_BASE_ZERO);

    cout << "create Jacobi preconditioner" << endl;
    
    Array<int> temp_ind (height+1); 
    Array<int> temp_cols (height);
    Array<double> temp_vals (height);

    for (int i = 0; i <= height; i++) temp_ind[i] = i; 

    for (int i = 0; i < height; i++)
      {
	temp_cols[i] = i;
	if (freedofs.Test(i))
	  temp_vals[i] = 1.0 / mat(i,i);
	else
	  temp_vals[i] = 0.0;
      }

    cudaMalloc ((void**)&dev_ind, (height+1) * sizeof(int));
    cudaMalloc ((void**)&dev_col, (height) * sizeof(int));
    cudaMalloc ((void**)&dev_val, (height) * sizeof(double));

    cudaMemcpy (dev_ind, &temp_ind[0], (height+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy (dev_col, &temp_cols[0], height*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy (dev_val, &temp_vals[0], height*sizeof(double), cudaMemcpyHostToDevice);
  }
  
  void DevJacobiPreconditioner :: Mult (const BaseVector & x, BaseVector & y) const
  {
    cout << "device mult precond" << endl;

    const UnifiedVector & ux = dynamic_cast<const UnifiedVector&> (x);
    UnifiedVector & uy = dynamic_cast<UnifiedVector&> (y);

    ux.UpdateDevice();
    uy = 0.0;
    uy.UpdateDevice();

    double alpha= 1;
    double beta = 0;
    cusparseDcsrmv (Get_CuSparse_Handle (), 
                    CUSPARSE_OPERATION_NON_TRANSPOSE, height, width, nze, 
		    &alpha, *descr, 
		    dev_val, dev_ind, dev_col, 
		    ux.dev_data, &beta, uy.dev_data);

    uy.host_uptodate = false;

    cout << "mult complete" << endl;
  }


  void DevJacobiPreconditioner :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    cout << "device multadd precond" << endl;

    const UnifiedVector & ux = dynamic_cast<const UnifiedVector&> (x);
    UnifiedVector & uy = dynamic_cast<UnifiedVector&> (y);

    ux.UpdateDevice();
    uy.UpdateDevice();

    double alpha= 1;
    double beta = 1;
    cusparseDcsrmv (Get_CuSparse_Handle (), 
                    CUSPARSE_OPERATION_NON_TRANSPOSE, height, width, nze, 
		    &alpha, *descr, 
		    dev_val, dev_ind, dev_col, 
		    ux.dev_data, &beta, uy.dev_data);

    uy.host_uptodate = false;
    cout << "mult complete" << endl;
  }
  















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

}




#endif
