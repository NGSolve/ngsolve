/*********************************************************************/
/* File:   cuda_linalg.cpp                                           */
/* Author: Joachim Schoeberl, Matthias Hochsteger                    */
/* Date:   11. Aug. 2014                                             */
/*********************************************************************/

#include <la.hpp>
#include "cuda_linalg.hpp"
#include <gpuwrapper.hpp>
#include <gpukernel.hpp>

namespace ngla
{
  void InitCuLinalg()
  {
    // pay the first nvrtc compile and module load at import, not
    // inside the user's first operator application
    auto dev = ngs_gpu::GetDevice();
    if (!dev) return;
    auto kernel = dev->CompileSource (ngs_gpu::code_gpukernel + "KERNEL(noop, VALUE(int,n)) { }")->GetKernel ("noop");
    dev->DefaultQueue()->Launch (*kernel, ngs_gpu::Dim3(1), ngs_gpu::Dim3(1), { ngs_gpu::KernelArg(int(0)) });
    cudaDeviceSynchronize();
  }


  shared_ptr<BaseMatrix> CreateDevMatrix (BaseMatrix & mat)
  {
    if (auto res = mat.CreateDeviceMatrix())
      return res;
    else
      throw Exception(string("matrix type not supported: ") + typeid(mat).name());
  }
}
