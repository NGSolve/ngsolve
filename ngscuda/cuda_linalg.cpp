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

  
  void InitCuLinalg()
  {
    // with lazy module loading the first kernel launch of this library
    // finalizes the whole fatbin (~300ms) - pay that at import, not
    // inside the user's first operator application
    DeviceParallelFor (1, [] DEVICE_LAMBDA (size_t) { });
    cudaDeviceSynchronize();

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
    



    
    
  }


  /******************** DevMatrix ********************/

  shared_ptr<BaseMatrix> CreateDevMatrix (BaseMatrix & mat)
  {
    if (auto res = mat.CreateDeviceMatrix())
      return res;
    else
      throw Exception(string("matrix type not supported: ") + typeid(mat).name());
  }

  bool synckernels = true;
}
