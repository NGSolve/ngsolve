#include <la.hpp>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cuda_runtime_api.h>

#include "cuda_linalg.hpp"


namespace ngla
{
  
  UnifiedVector :: UnifiedVector (int asize)
  {
    this->size = asize;
    cout << IM(7) << "Create unified vector, size = " << size << endl;

    host_data = new double[size];
    auto err = cudaMalloc((void**)&dev_data, size*sizeof(double));
    if (err != 0)
      throw Exception("UnifiedVector allocation error, ec="+ToString(err));
    
    cusparseCreateDnVec (&descr, size, dev_data, CUDA_R_64F);

    host_uptodate = false;
    dev_uptodate = false;
  }

  UnifiedVector :: UnifiedVector (const BaseVector& vec) : UnifiedVector(vec.Size())
  {
    (*this) = vec;
  }

  UnifiedVector :: ~UnifiedVector ()
  {
    /* cerr << "dtor UnifiedVector" << endl; */

    cusparseDestroyDnVec(descr);
    cudaFree(dev_data);
    delete[] host_data;
  }

  BaseVector & UnifiedVector :: operator= (double d)
  {
    /* for (int i = 0; i < size; i++) host_data[i] = d; */
    host_uptodate = false;

    ::SetScalar (d, size, dev_data); 
    // cublasDscal(Get_CuBlas_Handle(), size, &d, dev_data, 1); 
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

  const double & UnifiedVector :: operator [] (const int ind) const
  {
    /* cerr << "UnifiedVector operator[]" << endl; */
    UpdateHost(); 
    return host_data[ind];
  }

  double & UnifiedVector :: operator [] (const int ind)
  {
    UpdateHost(); 
    dev_uptodate = false;
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
        FVDouble() += scal * v.FVDouble();
      }

    return *this;
  }
  
  double UnifiedVector :: InnerProduct (const BaseVector & v2, bool conjugate) const
  {
    if (conjugate)
      throw Exception("conjugate in innerproduct not implemented yet.");

    static Timer tdot("CUDA InnerProduct");
    RegionTimer reg(tdot);

    const UnifiedVector * uv2 = dynamic_cast<const UnifiedVector*> (&v2);
    if (uv2)
    {
      static Timer tdot("CUDA InnerProduct devdev");
      RegionTimer reg(tdot);
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

  // TODO: maybe remove. mainly for testing
  ostream & UnifiedVector :: PrintStatus (ostream & ost) const
  {
    cout << "output unified vector of size " << size;
    cout << ", host = " << host_uptodate << ", dev = " << dev_uptodate << endl;
    return ost;
  }
  
  /* void UnifiedVector :: PrintDevice () const */
  /* { */
  /*   int DSIZE = size * sizeof(double); */
  /*   double *tmp = (double*) malloc(DSIZE); */
  /*   cudaMemcpy(tmp, dev_data, DSIZE, cudaMemcpyDeviceToHost); */
  /*   cout << "device up-to-date: " << dev_uptodate << endl; */
  /*   for (int i=0; i<size; i++) */
  /*     cout << tmp[i] << endl; */
  /* } */

  AutoVector UnifiedVector :: CreateVector () const
  {
    return make_unique<UnifiedVector> (size);
  }

  void UnifiedVector :: UpdateHost () const
  {
    if (host_uptodate)
      return;

    if (dev_uptodate)
      cudaMemcpy (host_data, dev_data, sizeof(double)*size, cudaMemcpyDeviceToHost);    
    /* else */
    /*   cout << "ERROR UnifiedVector::UpdateHost non is uptodate" << endl; */

    host_uptodate = true;
  }

  void UnifiedVector :: UpdateDevice () const
  {
    if (dev_uptodate)
      return;

    if (host_uptodate)
      cudaMemcpy (dev_data, host_data, sizeof(double)*size, cudaMemcpyHostToDevice);
    /* else */
    /*   cout << "ERROR UnifiedVector::UpdateDevice non is uptodate" << endl; */

    cout << "Host2Device copy!" << endl;
    
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
    UpdateHost(); 
    return host_data;
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

}
