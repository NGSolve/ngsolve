#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>

#include <comp.hpp>
#include <python_comp.hpp>

#include "ngsmetal.hpp"
#include "metal_vector.hpp"
#include "metal_btdtb.hpp"

#include "metal_benchmarks.hpp"
#include "metal_device.hpp"

using namespace ngsmetal;

namespace ngsmetal
{
  extern double Metal_Shared_Benchmark ();
}

namespace
{
  class MetalCaptureContext
  {
    std::string filename;
    bool capturing = false;

  public:
    explicit MetalCaptureContext (std::string afilename)
      : filename(std::move(afilename))
    {
      if (filename.empty())
        throw Exception("MetalCaptureManager: filename must not be empty");
      if (!filename.ends_with(".gputrace"))
        throw Exception("MetalCaptureManager: filename must end with '.gputrace'");
    }

    MetalCaptureContext (const MetalCaptureContext &) = delete;
    MetalCaptureContext & operator= (const MetalCaptureContext &) = delete;

    ~MetalCaptureContext()
    {
      Stop();
    }

    MetalCaptureContext & Enter()
    {
      auto capture_manager = MTL::CaptureManager::sharedCaptureManager();
      if (capturing || capture_manager->isCapturing())
        throw Exception("MetalCaptureManager: a Metal capture is already active");

      auto command_queue = GetCommandQueue();
      if (!command_queue)
        throw Exception("MetalCaptureManager: no Metal command queue is available");

      if (!capture_manager->supportsDestination(MTL::CaptureDestinationGPUTraceDocument))
        throw Exception("MetalCaptureManager: GPU trace documents are not supported\n"
                        "Ensure MTL_CAPTURE_ENABLED=1 was set before importing ngsolve.ngsmetal");

      auto descriptor = MTL::CaptureDescriptor::alloc()->init();
      descriptor->setDestination(MTL::CaptureDestinationGPUTraceDocument);
      descriptor->setCaptureObject(command_queue);

      auto path = NS::String::string(filename.c_str(), NS::UTF8StringEncoding);
      descriptor->setOutputURL(NS::URL::fileURLWithPath(path));

      NS::Error * error = nullptr;
      bool started = capture_manager->startCapture(descriptor, &error);
      descriptor->release();

      if (!started)
        {
          std::string message = "MetalCaptureManager: failed to start capture";
          if (error && error->localizedDescription())
            message += std::string(": ") + error->localizedDescription()->utf8String();
          message += ". Ensure MTL_CAPTURE_ENABLED=1 was set before importing ngsolve.ngsmetal";
          throw Exception(message);
        }

      capturing = true;
      return *this;
    }

    void Stop()
    {
      if (capturing)
        {
          MTL::CaptureManager::sharedCaptureManager()->stopCapture();
          capturing = false;
        }
    }
  };
}


PYBIND11_MODULE(libngsmetal, m)
{
  InitMetalDevice();      // register as ngs_gpu backend

  py::class_<MetalCaptureContext> (m, "MetalCaptureManager", R"doc(
Capture Metal commands submitted inside a Python ``with`` block.

Set ``MTL_CAPTURE_ENABLED=1`` before importing ``ngsolve.ngsmetal``. If the
module is already loaded, restart the Python process or notebook kernel. The
filename must end with ``.gputrace``; the resulting document can be opened in
Xcode.
)doc")
    .def(py::init<std::string>(), py::arg("filename"))
    .def("__enter__", &MetalCaptureContext::Enter,
         py::return_value_policy::reference_internal)
    .def("__exit__", [] (MetalCaptureContext & self,
                         py::object, py::object, py::object)
         {
           self.Stop();
           return false;
         })
    ;
  
  py::class_<MetalVector, BaseVector, shared_ptr<MetalVector>> (m, "MetalVector")
    .def(py::init<size_t>(), py::arg("size"))
    .def(py::init<const BaseVector&>(), py::arg("size"))
    .def("D2H", [](MetalVector& mv) {
      auto tmp = make_shared<MetalVector>(mv.Size(), MemType::Shared);
      tmp -> Set(1.0, mv);
      return tmp;
    })
    .def("WaitUntilCompleted", [](MetalVector& mv) { mv.GetQueue()->Finish(); })
    ;


  BaseVector::RegisterDeviceVectorCreator(typeid(S_BaseVectorPtr<double>),
                                          [] (const BaseVector & vec, bool unified) -> shared_ptr<BaseVector>
                                          {
                                            return make_shared<MetalVector>
                                              (vec, unified ? MemType::Shared : MemType::Device);
                                          });
  BaseVector::RegisterDeviceVectorCreator(typeid(VVector<double>),
                                          [] (const BaseVector & vec, bool unified) -> shared_ptr<BaseVector>
                                          {
                                            return make_shared<MetalVector>
                                              (vec, unified ? MemType::Shared : MemType::Device);
                                          });
  
  


  py::class_<MetalBTDTBMatrix, BaseMatrix, shared_ptr<MetalBTDTBMatrix>> (m, "MetalBTDTBMatrix")
    .def(py::init<const BaseMatrix&>(), py::arg("mat"))
    ;

  BaseMatrix::RegisterDeviceMatrixCreator(typeid(MatrixFreeBTDTB),
                                          [] (const BaseMatrix & bmat) -> shared_ptr<BaseMatrix>
                                          {
                                            return make_shared<MetalBTDTBMatrix>(bmat);
                                          });

  py::class_<Metal_MM_Benchmark> (m, "Metal_MM_Benchmark")
    .def(py::init<int,int,int,int,int>(), py::arg("n"), py::arg("m"), py::arg("k"), py::arg("lda"), py::arg("ldb"))
    .def("Run", &Metal_MM_Benchmark::Run, py::arg("timing")=true)
    ;
                                  
  
  m.def("Metal_Shared_Benchmark", Metal_Shared_Benchmark);
}
