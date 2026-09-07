#include <core/python_ngcore.hpp>

#include "../ngstd/python_ngstd.hpp"
#include "cuda_linalg.hpp"
#include "cuda_profiler.hpp"
#include "cuda_device.hpp"

// TODO: always use ngs_cuda?
using namespace ngbla;
using namespace ngla;
using namespace ngs_cuda;

namespace ngla {
  extern bool synckernels;
}

PYBIND11_MODULE(_ngscuda, m) {

  InitCUDA(1);
  InitCudaDevice();      // register as ngs_gpu backend
  InitCuLinalg();

  m.def("InitCuLinalg", &InitCuLinalg, "Initializing cuda linalg.");
  
  py::class_<DevMatrix, BaseMatrix, shared_ptr<DevMatrix>>
    (m, "DevBaseMatrix", "device matrix for CUDA applications");


  m.def("CreateDevMatrix", [] (BaseMatrix &mat)
          {
            return CreateDevMatrix(mat);
          });
    

  m.def("__time_tracer__", TimeProfiler);
  m.def("SetCudaTimer", CudaRegionTimer::SetCudaTimer);
  
  m.def("SetSyncKernels", [](bool sync) { synckernels = sync; });
  


  // ExportDemo(m);
}

