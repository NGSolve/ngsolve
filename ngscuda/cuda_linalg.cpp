/*********************************************************************/
/* File:   cuda_linalg.cpp                                           */
/* Author: Joachim Schoeberl, Matthias Hochsteger                    */
/* Date:   11. Aug. 2014                                             */
/*********************************************************************/


/*
 *	TODO:
 *		-) operator+
 *		-) float / template
 *		-) complex
 * */

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
    size = asize;
    host_data = new double[asize];

		cudaMalloc((void**)&dev_data, size*sizeof(double));

		cusparseCreateDnVec (&descr, size, dev_data, CUDA_R_64F);


    host_uptodate = false;
    dev_uptodate = false;

    (*this) = 0.0;
  }

	UnifiedVector :: UnifiedVector (BaseVector& vec)
	{
		// TODO: check, fix and pybind
		cerr << "(BaseVector) ctor UnifiedVector" << endl;
    host_data = new double[size];
    cudaMalloc((void**)&dev_data, size*sizeof(double));

		cusparseCreateDnVec (&descr, size, dev_data, CUDA_R_64F);

    host_uptodate = false;
    dev_uptodate = false;

		(*this) = vec;
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
    // cout << "Inner Prod" << endl;

    const UnifiedVector * uv2 = dynamic_cast<const UnifiedVector*> (&v2);
    if (uv2)
		{
			cerr << "cublas dot" << endl;
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
  
	void UnifiedVector :: PrintDevice () const
	{
		int DSIZE = size * sizeof(double);
		double *tmp = (double*) malloc(DSIZE);
		cudaMemcpy(tmp, dev_data, DSIZE, cudaMemcpyDeviceToHost);
		cerr << "device up-to-date: " << dev_uptodate << endl;
		for (int i=0; i<size; i++)
		{
			cerr << tmp[i] << endl;
		}
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

		/* cerr << "CREATING!" << endl; */
		/* int *host_row = (int*) malloc((mat.Height()+1) * sizeof(int)); */
		/* int *host_col = (int*) malloc(mat.Width() * sizeof(int)); */
		/* double *host_val = (double*) malloc(mat.NZE() * sizeof(double)); */
		/* cudaMemcpy(host_row, dev_ind, (mat.Height() + 1) * sizeof(int), cudaMemcpyDeviceToHost); */
		/* cudaMemcpy(host_col, dev_col, mat.Width() * sizeof(int), cudaMemcpyDeviceToHost); */
		/* cudaMemcpy(host_val, dev_val, mat.NZE() * sizeof(double), cudaMemcpyDeviceToHost); */
		/* for (int i=0; i<mat.Height() +1 ; i++) */
		/* { */
		/* 	cerr << host_row[i] << " "; */
		/* } */
		/* cerr << endl; */
		/* for (int i=0; i<mat.Width(); i++) */
		/* { */
		/* 	cerr << host_col[i] << " "; */
		/* } */
		/* cerr << endl; */
		/* for (int i=0; i<mat.NZE(); i++) */
		/* { */
		/* 	cerr << host_val[i] << " "; */
		/* } */
		/* cerr << endl; */


		/* for (int i=0; i<mat.Height()+1; i++) */
		/* { */
		/* 	cerr << temp_ind[i] << " "; */
		/* } */
		/* cerr << endl; */
		/* for (int i=0; i<mat.NZE()+1; i++) */
		/* { */
		/* 	cerr << mat.GetRowIndices(0)[i] << " "; */
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
		cerr << "buffersize: " << bufferSize << endl;
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
  





	// TODO:
	// 	redo. similiar to DevSparseMatrix
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
		else
		{
			throw Exception(string("matrix type not supported: ") + typeid(mat).name());
		}
	}




}
