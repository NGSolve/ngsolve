#ifndef NGSOLVE_PYTHON_COMP_HPP
#define NGSOLVE_PYTHON_COMP_HPP

#include <python_ngstd.hpp>
#include <python_linalg.hpp>


#include "gridfunction.hpp"
#include "periodic.hpp"

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
  inline void PyCheckId (const MeshAccess & ma, ElementId ei)
  {
    PyCheckIndex (ei.Nr(), ma.GetNE(ei.VB()), "element number");
  }

  inline void PyCheckBonusOrder (int order)
  {
    if (order < -1000 || order > 1000)
      throw py::value_error("bonus_intorder " + std::to_string(order) + " out of range [-1000, 1000]");
  }

  inline void PyCheckId (const MeshAccess & ma, NodeId ni)
  {
    PyCheckIndex (ni.GetNr(), ma.GetNNodes(ni.GetType()), "node number");
  }

  // TODO: use Archive structure to pickle fespaces
  inline py::tuple fesPickle (const FESpace& fes)
  {
    auto flags = fes.GetFlags();
    auto mesh = fes.GetMeshAccess();
    auto type = fes.type;

    // orders of regions, elements and single nodes are not in the flags.
    // a linked space has the base orders in its flag orderlink
    auto & spec = fes.GetOrderSpecification();
    auto & base = *spec.base;
    bool nodes = spec.node_order_level == mesh->GetNLevels();
    bool withbase = !flags.AnyFlagDefined ("orderlink") &&
      (base.region_orders.Size() || base.element_order.Size() || base.has_links);
    if (!withbase && !nodes)
      return py::make_tuple(type,mesh,flags);

    py::list nodeorders;
    for (int i = 0; i < 3; i++)
      nodeorders.append (nodes ? py::bytes ((const char*)spec.node_order[i].Data(),
                                            spec.node_order[i].Size()*sizeof(short)) : py::bytes());
    return py::make_tuple(type,mesh,flags,
                          py::make_tuple (withbase ? py::cast(spec.base) : py::object(py::none()),
                                          nodes, nodeorders));
  }
  
  template<typename FESPACE>
  shared_ptr<FESPACE> fesUnpickle(py::tuple state)
  {
    auto mesh = state[1].cast<shared_ptr<MeshAccess>>();
    auto fes = CreateFESpace(state[0].cast<string>(), mesh,
                             state[2].cast<Flags>());
    if (state.size() > 3)
      {
        FESpace::OrderSpecification spec;
        py::tuple extra = state[3];
        bool creator = !extra[0].is_none();
        spec.base = creator ? extra[0].cast<shared_ptr<FESpace::BaseOrders>>()
          : fes->GetOrderSpecification().base;
        if (extra[1].cast<bool>())
          {
            auto nodeorders = extra[2].cast<py::list>();
            for (int i = 0; i < 3; i++)
              {
                std::string s = nodeorders[i].cast<std::string>();
                spec.node_order[i].SetSize (s.size()/sizeof(short));
                memcpy ((void*)spec.node_order[i].Data(), s.data(), s.size());
              }
            spec.node_order_level = mesh->GetNLevels();
          }
        fes->SetOrderSpecification (spec, creator);
      }
    fes->Update();
    fes->FinalizeUpdate();
    // MR: connect_auto_update?
    return dynamic_pointer_cast<FESPACE>(fes);
  };


  inline void connect_auto_update(FESpace* fes) {
    fes->ConnectAutoUpdate();
  }

  /*
  inline void connect_auto_update(GridFunction* gf) {
    gf->ConnectAutoUpdate();
  }
  */
  
  template <typename FES, typename BASE=FESpace>
  auto ExportFESpace (py::module & m, string pyname, bool module_local = false)
  {
    auto docu = FES::GetDocu();
    string docstring = docu.GetPythonDocString();
    // string docuboth = docu.short_docu + "\n\n" + docu.long_docu;
    auto pyspace = py::class_<FES, shared_ptr<FES>,BASE> (m, pyname.c_str(), docstring.c_str(), py::module_local(module_local));

    pyspace
      .def(py::init([pyspace](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                    {
                      py::list info;
                      info.append(ma);
                      auto flags = CreateFlagsFromKwArgs(kwargs, pyspace, info);
                      auto fes = make_shared<FES>(ma,flags);
                      if (flags.AnyFlagDefined ("orderlink") && !fes->IsOrderLinked())
                        throw Exception ("orderlink is not available for space " + fes->GetClassName());
                      fes->Update();
                      fes->FinalizeUpdate();
                      // connect_auto_update(fes.get());
                      fes -> ConnectAutoUpdate();
                      return fes;
                    }),py::arg("mesh"))
    
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
