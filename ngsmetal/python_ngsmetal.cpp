#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>

#include <comp.hpp>
#include <python_comp.hpp>

#include "ngsmetal.hpp"
#include "metal_vector.hpp"
#include "metal_btdtb.hpp"

using namespace ngsmetal;


PYBIND11_MODULE(libngsmetal, m)
{
  cout << "Loading ngs-metal library" << endl;

  py::class_<MetalVector, BaseVector, shared_ptr<MetalVector>> (m, "MetalVector")
    .def(py::init<size_t>(), py::arg("size"))
    .def(py::init<const BaseVector&>(), py::arg("size"))    
    ;


  BaseVector::RegisterDeviceVectorCreator(typeid(S_BaseVectorPtr<double>),
                                          [] (const BaseVector & vec, bool unified) -> shared_ptr<BaseVector>
                                          {
                                            return make_shared<MetalVector>(vec);
                                          });
  BaseVector::RegisterDeviceVectorCreator(typeid(VVector<double>),
                                          [] (const BaseVector & vec, bool unified) -> shared_ptr<BaseVector>
                                          {
                                            return make_shared<MetalVector>(vec);
                                          });
  
  


  py::class_<MetalBTDTBMatrix, BaseMatrix, shared_ptr<MetalBTDTBMatrix>> (m, "MetalBTDTBMatrix")
    .def(py::init<const BaseMatrix&>(), py::arg("mat"))
    ;

  BaseMatrix::RegisterDeviceMatrixCreator(typeid(MatrixFreeBTDTB),
                                          [] (const BaseMatrix & bmat) -> shared_ptr<BaseMatrix>
                                          {
                                            // auto &mat  = dynamic_cast<const MatrixFreeBTDTB&>(bmat);
                                            return make_shared<MetalBTDTBMatrix>(bmat);
                                          });
  
  
}
