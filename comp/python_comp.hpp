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

    // orders of regions, elements and single nodes are not in the flags
    auto & spec = fes.GetOrderSpecification();
    bool nodes = spec.node_order_level == mesh->GetNLevels();
    if (!spec.region_orders.Size() && !spec.element_order.Size() && !nodes)
      return py::make_tuple(type,mesh,flags);

    auto bytes = [] (const auto & a)
    { return py::bytes ((const char*)a.Data(), a.Size()*sizeof(a[0])); };
    py::list regions;
    for (auto & ro : spec.region_orders)
      regions.append (py::make_tuple (int(ro.vb), make_shared<BitArray>(ro.mask),
                                      ro.order, ro.order_dir, ro.direction));
    py::list nodeorders;
    for (int i = 0; i < 3; i++)
      nodeorders.append (nodes ? bytes(spec.node_order[i]) : py::bytes());
    return py::make_tuple(type,mesh,flags,
                          py::make_tuple (regions, bytes(spec.element_order), nodes, nodeorders));
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
        py::tuple orders = state[3];
        auto frombytes = [] (py::handle h, auto & a)
        {
          std::string s = h.cast<std::string>();
          a.SetSize (s.size()/sizeof(a[0]));
          memcpy ((void*)a.Data(), s.data(), a.Size()*sizeof(a[0]));
        };
        for (auto item : orders[0].cast<py::list>())
          {
            auto t = item.cast<py::tuple>();
            spec.region_orders.Append
              (FESpace::OrderSpecification::RegionOrder
               { VorB(t[0].cast<int>()), BitArray(*t[1].cast<shared_ptr<BitArray>>()),
                 t[2].cast<int>(), t[3].cast<int>(),
                 t[4].is_none() ? nullptr : t[4].cast<shared_ptr<CoefficientFunction>>() });
          }
        frombytes (orders[1], spec.element_order);
        if (orders[2].cast<bool>())
          {
            auto nodeorders = orders[3].cast<py::list>();
            for (int i = 0; i < 3; i++)
              frombytes (nodeorders[i], spec.node_order[i]);
            spec.node_order_level = mesh->GetNLevels();
          }
        fes->SetOrderSpecification (spec);
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
