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
    
    host_uptodate = false;
    dev_uptodate = false;
  }

  UnifiedVector :: UnifiedVector (const BaseVector& vec) : UnifiedVector(vec.Size())
  {
    (*this) = vec;
    UpdateDevice();    
  }

  UnifiedVector :: UnifiedVector (const UnifiedVector& vec) : UnifiedVector(vec.Size())
  {
    (*this) = vec;
    UpdateDevice();    
  }
  
  // to be improved
  UnifiedVector :: UnifiedVector (UnifiedVector && vec)
    : UnifiedVector (vec.Size())
  {
    (*this) = vec;
    UpdateDevice();
  }
  
  UnifiedVector :: ~UnifiedVector ()
  {
    cudaFree(dev_data);
    delete[] host_data;
  }

  BaseVector & UnifiedVector :: operator= (double d)
  {
    // ::SetScalar (d, size, dev_data); 
    // ::SetScalar (d, FlatVector<Dev<double>> (size, dev_data));
    ::SetScalar (d, FVDev()); 
      
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
            UpdateDevice();            
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
    UpdateDevice();
    return *this;
  }

  UnifiedVector & UnifiedVector :: operator= (const UnifiedVector & v2)
  {
    if (v2.dev_uptodate)
      {
        cudaMemcpy (dev_data, v2.dev_data, sizeof(double)*size, cudaMemcpyDeviceToDevice);    
        dev_uptodate = true;
        host_uptodate = false;
      }
    else if (v2.host_uptodate)
      {
        FVDouble() = v2.FVDouble();
        host_uptodate = true;
        dev_uptodate = false;
        UpdateDevice();            
      }
    else
      {
        cerr << "operator=UnifiedVector - not up to date" << endl;
      }
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

  AutoVector UnifiedVector :: Range (T_Range<size_t> range) const
  {
    return make_unique<UnifiedVectorWrapper>(*this, range);
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
    UnifiedVectorWrapper uv(v);
    uv.UpdateDevice();
    SetVector (scal, Size(), uv.DevData(), DevData());
    host_uptodate = false;
    return *this;
  }
  
  
  BaseVector & UnifiedVector :: Add (double scal, const BaseVector & v)
  {
    if (auto v2 = dynamic_cast<const UnifiedVector*> (&v))
      {
        UpdateDevice();
        v2->UpdateDevice();
        /*
        cublasDaxpy (Get_CuBlas_Handle(), 
                           size, &scal, v2->dev_data, 1, dev_data, 1);
        */
        // MyDaxpy (scal, size, v2->dev_data, dev_data);

        DeviceParallelFor
          (size, [scal, x=v2->dev_data, y=dev_data] DEVICE_LAMBDA (auto tid) 
           {
             y[tid] += scal*x[tid];
           }); 
        
        host_uptodate = false;
      }
    else
      {
        FVDouble() += scal * v.FVDouble();
      }

    return *this;
  }
  
  double UnifiedVector :: InnerProductD (const BaseVector & v2) const
  {
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

  double UnifiedVector :: L2Norm() const
  {
    UpdateDevice();
    double res;
    cublasDnrm2(Get_CuBlas_Handle(), size, dev_data, 1, &res);
    return res;
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

  FlatVector<Dev<double>> UnifiedVector :: FVDev() const
  {
    UpdateDevice();
    InvalidateHost();
    return { Size(), (Dev<double>*)dev_data };
  }

  FlatVector<Dev<double>> UnifiedVector :: FVDevRO() const
  {
    UpdateDevice();
    return { Size(), (Dev<double>*)dev_data };
  }
  
  void * UnifiedVector :: Memory() const throw()
  { 
    UpdateHost(); 
    return host_data;
  }

  UnifiedVectorWrapper :: UnifiedVectorWrapper(const BaseVector & vec_, optional<IntRange> opt_range)
    : vec(vec_)
  {
    IntRange range = {0, vec.Size()};
    if(opt_range)
      range = *opt_range;
    this->size = range.Size();

    auto uptr = dynamic_cast<const UnifiedVector*>(&vec_);
    if(uptr)
    {
      host_data = uptr->HostData() + range.First();
      dev_data = uptr->DevData() + range.First();
      uptr->UpdateDevice();
      uptr->InvalidateHost();
      initial_host_uptodate =  uptr->IsHostUptodate();
      initial_dev_uptodate =  uptr->IsDevUptodate();
    }
    else
    {
      auto err = cudaMalloc((void**)&dev_data, size*sizeof(double));
      if (err != 0)
        throw Exception("UnifiedVector allocation error, ec="+ToString(err));
      initial_host_uptodate = true;
      initial_dev_uptodate = false;
      host_data = vec.FVDouble().Data() + range.First();
    }
    host_uptodate = initial_host_uptodate;
    dev_uptodate = initial_dev_uptodate;
  }

  UnifiedVectorWrapper :: ~UnifiedVectorWrapper()
  {
    if(initial_host_uptodate && !host_uptodate)
      UpdateHost();
    if(initial_dev_uptodate && !dev_uptodate)
      UpdateDevice();

    host_data = nullptr;
    auto uptr = dynamic_cast<const UnifiedVector*>(&vec);
    if(uptr)
      dev_data = nullptr;
  }

}
