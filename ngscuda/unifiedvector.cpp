#include <la.hpp>
#include "cuda_linalg.hpp"

namespace ngla
{
  
  UnifiedVector :: UnifiedVector (int asize)
  {
    this->size = asize;

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
    cusparseDestroyDnVec(descr);
    cudaFree(dev_data);
    delete[] host_data;
  }

  BaseVector & UnifiedVector :: operator= (double d)
  {
    ::SetScalar (d, size, dev_data); 

    host_uptodate = false;
    dev_uptodate = true;
    
    return *this;
  }

  BaseVector & UnifiedVector :: operator= (const BaseVector & v2)
  {
    if (auto uv2 = dynamic_cast<const UnifiedVector*> (&v2))
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

    FVDouble() = v2.FVDouble();

    host_uptodate = true;
    dev_uptodate = false;
    return *this;
  }

  const double & UnifiedVector :: operator [] (const int ind) const
  {
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
    if (auto v2 = dynamic_cast<const UnifiedVector*> (&v))
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

    if (auto uv2 = dynamic_cast<const UnifiedVector*> (&v2))
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

  
  AutoVector UnifiedVector :: CreateVector () const
  {
    return make_unique<UnifiedVector> (size);
  }

  void UnifiedVector :: UpdateHost () const
  {
    if (host_uptodate)
      return;

    if (dev_uptodate)
      {
        cudaMemcpy (host_data, dev_data, sizeof(double)*size, cudaMemcpyDeviceToHost);
        cout << IM(5) << "Device2Host copy!" << endl;        
      }
    
    host_uptodate = true;
  }

  void UnifiedVector :: UpdateDevice () const
  {
    if (dev_uptodate)
      return;

    if (host_uptodate)
      {
        cudaMemcpy (dev_data, host_data, sizeof(double)*size, cudaMemcpyHostToDevice);
        cout << IM(5) << "Host2Device copy!" << endl;
      }
    
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
}
