#ifndef NGSOLVE_PYTHON_COMP_HPP
#define NGSOLVE_PYTHON_COMP_HPP

#include <python_ngstd.hpp>

#include "comp.hpp"

namespace pybind11
{
  // specialized implementation of polymorphic_type_hook to allow
  // unregistered derived types convert to registered derived types
  // by trying to dynamic cast
  template<>
  struct polymorphic_type_hook<ngcomp::FESpace>
  {
    static const void* get(const ngcomp::FESpace* src, const type_info*& type)
    {
      // for example this could be a QuasiPeriodicFESpace<double>, ...
      if(auto cast = dynamic_cast<const ngcomp::PeriodicFESpace*>(src))
        {
          type = &typeid(ngcomp::PeriodicFESpace);
          return cast;
        }
      // if it's not one of these spaces use the default implementation
      // for dynamic types
      type = src ? &typeid(*src) : nullptr;
      return dynamic_cast<const void*>(src);
    }
  };
}

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
    fes->Update();
    fes->FinalizeUpdate();
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
                      fes->Update();
                      fes->FinalizeUpdate();
                      if(autoupdate)
                        {
                          auto fesptr = fes.get();
                          ma->updateSignal.Connect(fesptr, [fesptr]()
                                     {
                                       fesptr->Update();
                                       fesptr->FinalizeUpdate();
                                     });
                        }
                      return fes;
                    }),py::arg("mesh"), py::arg("autoupdate")=false)
    
      .def(py::pickle(&fesPickle,
                      (shared_ptr<FES>(*)(py::tuple)) fesUnpickle<FES>))
      ;

    pyspace.def_static("__flags_doc__", [docu]()
                                        {
                                          py::dict flags_doc;
                                          for (auto & flagdoc : FES::GetDocu().arguments)
                                            flags_doc[get<0> (flagdoc).c_str()] = get<1> (flagdoc);
                                          return flags_doc;
                                         });
      
    return pyspace;
  }
}
#endif // NGSOLVE_PYTHON_COMP_HPP
