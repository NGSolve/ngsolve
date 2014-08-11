#ifdef CUDA 

/*********************************************************************/
/* File:   cuda_linalg.cpp                                           */
/* Author: Joachim Schoeberl, Matthias Hochsteger                    */
/* Date:   11. Aug. 2014                                             */
/*********************************************************************/


#include <la.hpp>


namespace ngla
{
  
  UnifiedVector :: UnifiedVector (int asize)
  {
    size = asize;
    host_data = new double[asize];
    // cudaMalloc((double**)&dev_data, size*sizeof(double));
    cudaMalloc(&dev_data, size*sizeof(double));
    host_uptodate = false;
    dev_uptodate = false;
  }

  BaseVector & UnifiedVector :: operator= (double d)
  {
    for (int i = 0; i < size; i++) host_data[i] = d;
    host_uptodate = true;
    UpdateDevice();
    return *this;

    /*
    if (dev_uptodate || !host_uptodate)
      {
        // dev_dat = d  // howto ????
        dev_uptodate = true;
        host_uptodate = false;
      }
    else
      {
        for (int i = 0; i < size; i++)
          host_data[i] = d;

        dev_uptodate = false;
        host_uptodate = true;
      }
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

  
  void UnifiedVector :: UpdateHost ()
  {
    if (!dev_uptodate) cout << "device not uptodate" << endl;
    cudaMemcpy (host_data, dev_data, sizeof(double)*size, cudaMemcpyDeviceToHost);    
    host_uptodate = true;
  }

  void UnifiedVector :: UpdateDevice ()
  {
    if (!host_uptodate) cout << "host not uptodate" << endl;
    cudaMemcpy (dev_data, host_data, sizeof(double)*size, cudaMemcpyHostToDevice);
    dev_uptodate = true;
  }

  
  
}

#endif
