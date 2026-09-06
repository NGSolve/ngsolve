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
  extern void InitApplyIntegrationPoints ();
  extern void InitBTDTB ();
  extern bool synckernels;
  
}

PYBIND11_MODULE(_ngscuda, m) {

  InitCUDA(1);
  InitCudaDevice();      // register as ngs_gpu backend
  InitCuLinalg();
  InitApplyIntegrationPoints();
  InitBTDTB();

  m.def("InitCuLinalg", &InitCuLinalg, "Initializing cuda linalg.");
  
  py::class_<UnifiedVector, BaseVector, shared_ptr<UnifiedVector>>
    (m, "UnifiedVector", "UnifiedVector for CUDA applications", py::multiple_inheritance())
    
    .def(py::init([] (int size)
                  { 
                    return make_shared<UnifiedVector>(size); 
                  }))
    .def(py::init([] (const BaseVector &vec) 
                  {
                    return make_shared<UnifiedVector>(vec);
                  }))
    .def(py::init([] (py::array_t<double> bvec)
                  {
                    auto vec = bvec.template unchecked<1>();
                    shared_ptr<UnifiedVector> uv = make_shared<UnifiedVector>(vec.size());
                    FlatVector<double> fv = uv->FV<double>();
                    for (size_t i = 0; i < vec.size(); i++)
                      {
                        fv(i) = vec(i);
                      }
                    return uv;
                  }))

    .def("UpdateHost", &UnifiedVector::UpdateHost)
    .def("UpdateDevice", &UnifiedVector::UpdateDevice)
    .def_property_readonly("__cuda_array_interface__", [](UnifiedVector& self)
    {
        self.UpdateDevice();
        auto ptr = reinterpret_cast<uintptr_t>(self.DevData());
        py::dict cai;
        cai["version"] = 2;
        cai["shape"]   = py::make_tuple(self.Size());
        // "<f8" = little-endian float64
        cai["typestr"] = "<f8";
        // data: (ptr, readonly_flag)
        cai["data"] = py::make_tuple(ptr, false);
        // contiguous 1D, so no strides
        cai["strides"] = py::none();
        return cai;
    })
    .def_property_readonly("dev_ptr", [](UnifiedVector& self)
    {
        return reinterpret_cast<uintptr_t>(self.DevData());
    })
    ;


  // the scalar of a UnifiedVector, with the operators of ngsolve.la.DeviceScalarD
  py::class_<UnifiedScalar, DeviceScalar<double>, shared_ptr<UnifiedScalar>>
    (m, "UnifiedScalar", "scalar on the cuda device")
    .def(py::init<double>(), py::arg("value")=0.0);

  py::class_<DevMatrix, BaseMatrix, shared_ptr<DevMatrix>>
    (m, "DevBaseMatrix", "device matrix for CUDA applications");



  m.def("CreateDevMatrix", [] (BaseMatrix &mat)
          {
            return CreateDevMatrix(mat);
          });
    

  m.def("__time_tracer__", TimeProfiler);
  m.def("SetCudaTimer", CudaRegionTimer::SetCudaTimer);
  
  m.def("SetSyncKernels", [](bool sync) { synckernels = sync; });
  


  py::class_<CudaGraph> (m, "CudaGraph")
    .def(py::init<>())
    .def("BeginCapture", &CudaGraph::BeginCapture)
    .def("EndCapture", &CudaGraph::EndCapture)
    .def("Launch", &CudaGraph::Launch)
    ;
  // DevCGSolver — preconditioned CG with CUDA graph capture of iteration body
  m.def("DevCGSolver",
  [](shared_ptr<BaseMatrix> mat,
     shared_ptr<BaseMatrix> pre,
     shared_ptr<BaseMatrix> adev_raw,
     shared_ptr<BaseMatrix> cdev_raw,
     double precision,
     int    maxsteps,
     bool   printrates)
  {
    auto solver = make_shared<DevCGSolver>(
        mat, pre, adev_raw, cdev_raw);
    solver->SetPrecision(precision);
    solver->SetMaxSteps(maxsteps);
    solver->SetPrintRates(printrates);
    return shared_ptr<KrylovSpaceSolver>(solver);
  },
  py::arg("mat"),
  py::arg("pre"),
  py::arg("adev_raw")   = nullptr,
  py::arg("cdev_raw")   = nullptr,
  py::arg("precision")  = 1e-12,
  py::arg("maxsteps")   = 400,
  py::arg("printrates") = false,
  "Preconditioned CG solver with CUDA graph capture of iteration body (convergence check via DtoH remains outside graph).");

  // DevTFQMRSolver — preconditioned TFQMR for non-symmetric systems
  m.def("DevTFQMRSolver",
  [](shared_ptr<BaseMatrix> mat,
     shared_ptr<BaseMatrix> pre,
     shared_ptr<BaseMatrix> adev_raw,
     shared_ptr<BaseMatrix> cdev_raw,
     double precision,
     int    maxsteps,
     bool   printrates)
  {
    auto solver = make_shared<DevTFQMRSolver>(
        mat, pre, adev_raw, cdev_raw);
    solver->SetPrecision(precision);
    solver->SetMaxSteps(maxsteps);
    solver->SetPrintRates(printrates);
    return shared_ptr<KrylovSpaceSolver>(solver);
  },
  py::arg("mat"),
  py::arg("pre"),
  py::arg("adev_raw")   = nullptr,
  py::arg("cdev_raw")   = nullptr,
  py::arg("precision")  = 1e-8,
  py::arg("maxsteps")   = 400,
  py::arg("printrates") = false,
  "Preconditioned TFQMR solver for non-symmetric systems on GPU.");
  // ExportDemo(m);
}

