#ifndef NGSOLVE_PYTHON_COMP_HPP
#define NGSOLVE_PYTHON_COMP_HPP

#include <python_ngstd.hpp>

#include "comp.hpp"

namespace ngcomp
{
  // TODO: use Archive structure to pickle fespaces
  inline py::tuple fesPickle (const FESpace& fes)
  {
    auto flags = fes.GetFlags();
    auto mesh = fes.GetMeshAccess();
    auto type = fes.type;
    return py::make_tuple(type,mesh,flags);
  }
  
  template<typename FESPACE>
  shared_ptr<FESPACE> fesUnpickle(py::tuple state)
  {
    auto fes = CreateFESpace(state[0].cast<string>(),
                             state[1].cast<shared_ptr<MeshAccess>>(),
                             state[2].cast<Flags>());

    LocalHeap glh(10000000, "Unpickl-lh");
    fes->Update(glh);
    fes->FinalizeUpdate(glh);
    return dynamic_pointer_cast<FESPACE>(fes);
  };

  template <typename FES, typename BASE=FESpace>
  auto ExportFESpace (py::module & m, string pyname)
  {
    auto docu = FES::GetDocu();
    string docuboth = docu.short_docu + "\n\n" + docu.long_docu;
    auto pyspace = py::class_<FES, shared_ptr<FES>,BASE> (m, pyname.c_str(), docuboth.c_str());

    pyspace
      .def(py::init([pyspace](shared_ptr<MeshAccess> ma, bool autoupdate, py::kwargs kwargs)
                    {
                      py::list info;
                      info.append(ma);
                      auto flags = CreateFlagsFromKwArgs(kwargs, pyspace, info);
                      auto fes = make_shared<FES>(ma,flags);
                      LocalHeap glh(10000000, "init-fes-lh");                    
                      fes->Update(glh);
                      fes->FinalizeUpdate(glh);
                      if(autoupdate)
                        {
                          auto fesptr = fes.get();
                          ma->updateSignal.Connect(fesptr, [fesptr]()
                                     {
                                       LocalHeap lh(1000000);
                                       fesptr->Update(lh);
                                       fesptr->FinalizeUpdate(lh);
                                     });
                        }
                      return fes;
                    }),py::arg("mesh"), py::arg("autoupdate")=false)
    
      .def(py::pickle(&fesPickle,
                      (shared_ptr<FES>(*)(py::tuple)) fesUnpickle<FES>))
      ;

    pyspace.def_static("__flags_doc__", [docu]()
                                        {
                                          auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                                              attr("FESpace").
                                                                              attr("__flags_doc__")());
                                          for (auto & flagdoc : docu.arguments)
                                            flags_doc[get<0> (flagdoc).c_str()] = get<1> (flagdoc);
                                          return flags_doc;
                                        });
      
    return pyspace;
  }
}
#endif // NGSOLVE_PYTHON_COMP_HPP
