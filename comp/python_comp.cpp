#ifdef NGS_PYTHON
#include <regex>

#include "../ngstd/python_ngstd.hpp"
#include "python_comp.hpp"
#include <comp.hpp>
#include <multigrid.hpp> 
#include <pybind11/functional.h>

#include "hdivdivfespace.hpp"
#include "hcurldivfespace.hpp"
#include "hcurlcurlfespace.hpp"
#include "normalfacetfespace.hpp"
#include "../fem/hdivdivfe.hpp"
#include "hdivdivsurfacespace.hpp"
#include "numberfespace.hpp"
#include "irspace.hpp"
#include "compressedfespace.hpp"
#include "../fem/integratorcf.hpp"
#include "../fem/h1lofe.hpp"
#include "contact.hpp"
#include "globalinterfacespace.hpp"
#include "globalspace.hpp"
using namespace ngcomp;

using ngfem::ELEMENT_TYPE;

typedef GridFunction GF;


namespace ngcomp
{
  void PatchwiseSolve (shared_ptr<SumOfIntegrals> bf,
                       shared_ptr<SumOfIntegrals> lf,
                       shared_ptr<GridFunction> gf,
                       LocalHeap & lh);
}

namespace ngfem
{
  extern bool symbolic_integrator_uses_diff;
}


class PyNumProc : public NumProc
{
public:
  using NumProc::NumProc;
  virtual void Do (LocalHeap &lh) override
  {
      auto pylh = py::cast(lh, py::return_value_policy::reference);
      try{
          PYBIND11_OVERLOAD_PURE(
                                 void,       /* Return type */
                                 PyNumProc,  /* Parent class */
                                 Do,         /* Name of function */
                                 pylh
          );
      }
      catch (py::error_already_set const &e) {
          cerr << e.what() << endl;
          PyErr_Print();
      }
  }
};



py::object ProxyNode2Py (const ProxyNode & node)
{
  if (auto proxy = *node)
    return py::cast(proxy);
  
  py::list l;
  for (auto & sub : node.list)
    l.append (ProxyNode2Py(sub));
  return py::object(l);
}

py::object MakeProxyFunction (shared_ptr<FESpace> fes,
                              bool testfunction) 
{
  return ProxyNode2Py(fes->GetProxyFunction(testfunction));
}


  shared_ptr<SumOfIntegrals> DualProxyFunction :: operator() (shared_ptr<CoefficientFunction> u) const
  {
    VorB vb = evaluator -> VB();
    Array<VorB> node_types { fes->GetDualShapeNodes(vb) };
    
    auto sum = make_shared<SumOfIntegrals>();
    for (auto nt : node_types)
      {
        DifferentialSymbol dx(vb, nt, false, 0);
        if (this -> Dimension() == 1)
          sum->icfs += make_shared<Integral> (const_cast<DualProxyFunction*>(this)->shared_from_this()  * u, dx);
        else
          sum->icfs += make_shared<Integral> (InnerProduct(const_cast<DualProxyFunction*>(this)->shared_from_this(), u), dx);
          
      }
    return sum;
  }
  


class GlobalDummyVariables 
{
public:
  int GetMsgLevel() { return printmessage_importance; }
  void SetMsgLevel(int msg_level) 
  {
    // cout << "set printmessage_importance to " << msg_level << endl;
    if(msg_level == 0)
      Logger::SetGlobalLoggingLevel(level::off);
    else if(msg_level > 0)
      Logger::SetGlobalLoggingLevel(level::err);
    else if(msg_level > 3)
      Logger::SetGlobalLoggingLevel(level::info);
    else if(msg_level > 6)
      Logger::SetGlobalLoggingLevel(level::debug);
    printmessage_importance = msg_level; 
    netgen::printmessage_importance = msg_level; 
  }
  string GetTestoutFile () const
  {
    return dynamic_cast<ofstream*>(testout) ? "testout set" : "no testout set";
  }
  void SetTestoutFile(string filename) 
  {
    // cout << "set testout-file to " << filename << endl;
    delete testout;
    testout = new ofstream(filename);
  }
  
};
static GlobalDummyVariables globvar;


void ExportNgcompMesh (py::module &m);

void NGS_DLL_HEADER ExportNgcomp(py::module &m)
{

  ExportNgcompMesh(m);
  //////////////////////////////////////////////////////////////////////////////////////////

  static size_t global_heapsize =
    (sizeof(size_t) == 8) ? 100000000 : 10000000;
  static LocalHeap glh(global_heapsize, "python-comp lh", true);

  class LocalHeapProvider
  {
    Array<LocalHeap*> heaps;
    std::mutex m;
    
    class BorrowedLocalHeap 
    {
      LocalHeap * lh;
      LocalHeapProvider * provider;
    public:
      BorrowedLocalHeap (LocalHeap * alh, LocalHeapProvider * aprovider)
        : lh(alh), provider(aprovider) { }
    
      ~BorrowedLocalHeap()
      {
        // cout << "returns lh to provider" << endl;
        provider -> ReturnLH(lh);
      }

      LocalHeap & LH() { return *lh; }
      operator LocalHeap& () { return *lh; }
    };
    
  
  public:
    BorrowedLocalHeap GetLH()
    {
      std::lock_guard lock(m);
      if (heaps.Size())
        {
          auto tmp = heaps.Last();
          heaps.SetSize(heaps.Size()-1);
          // cout << "reuse existing lh" << endl;        
          return BorrowedLocalHeap (tmp, this);
        }
      
      // cout << "create new lh" << endl;
      return BorrowedLocalHeap(new LocalHeap(global_heapsize, "python-comp lh", true), this);
    }
    
    void ReturnLH (LocalHeap * lh)
    {
      std::lock_guard lock(m);
      heaps.Append (lh);
    }
    void Clear()
    {
        for (auto lh : heaps)
            delete lh;
        heaps.SetSize(0);
    }
  };
  static LocalHeapProvider lhp;

  
  //////////////////////////////////////////////////////////////////////////////////////////

  py::enum_<COUPLING_TYPE> (m, "COUPLING_TYPE", docu_string(R"raw_string(
Enum specifying the coupling type of a degree of freedom, each dof is
either UNUSED_DOF, LOCAL_DOF, INTERFACE_DOF or WIREBASKET_DOF, other values
are provided as combinations of these:

UNUSED_DOF: Dof is not used, i.e the minion dofs in a :any:`Periodic` finite
    element space.

LOCAL_DOF: Inner degree of freedom, will be eliminated by static
    condensation and reconstructed afterwards.

HIDDEN_DOF: Inner degree of freedom, that will be eliminated by static
    condensation and *not* reconstruced afterwards(spares some entries).
    Note: 
     * without static condensation a HIDDEN_DOF is treated as any other
       DOF, e.g. as a LOCAL_DOF
     * To a HIDDEN_DOF the r.h.s. vector must have zero entries.
     * When static condensation is applied (eliminate_hidden/
       eliminate_internal) the block corresponding to HIDDEN_DOFs
       has to be invertible.

CONDENSABLE_DOF: Inner degree of freedom, that will be eliminated by static
    condensation (LOCAL_DOF or HIDDEN_DOF)

INTERFACE_DOF: Degree of freedom between two elements, these will not be
    eliminated by static condensation, but not be put into the wirebasket
    system for i.e. a bddc :any:`Preconditioner`.

NONWIREBASKET_DOF: Either a LOCAL_DOF or an INTERFACE_DOF

WIREBASKET_DOF: Degree of freedom coupling with many elements (more than
    one). These will be put into the system for a bddc preconditioner.
    The HCurl space also treats degrees of freedom of badly shaped
    elements as WIREBASKET_DOFs.

EXTERNAL_DOF: Either INTERFACE_DOF or WIREBASKET_DOF

VISIBLE_DOF: not UNUSED_DOF or HIDDEN_DOF

ANY_DOF: Any used dof (LOCAL_DOF or INTERFACE_DOF or WIREBASKET_DOF)

)raw_string"))
    .value("UNUSED_DOF", UNUSED_DOF)
    .value("HIDDEN_DOF", HIDDEN_DOF)
    .value("LOCAL_DOF", LOCAL_DOF)
    .value("CONDENSABLE_DOF", CONDENSABLE_DOF)
    .value("INTERFACE_DOF", INTERFACE_DOF)
    .value("NONWIREBASKET_DOF", NONWIREBASKET_DOF)
    .value("WIREBASKET_DOF", WIREBASKET_DOF)
    .value("EXTERNAL_DOF", EXTERNAL_DOF)
    .value("VISIBLE_DOF", VISIBLE_DOF)
    .value("ANY_DOF", ANY_DOF)
    // .export_values()
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  py::class_<ElementRange, IntRange> (m, "ElementRange")
    .def(py::init<const MeshAccess&,VorB,IntRange>(), py::arg("mesh"), py::arg("vb"), py::arg("range"))
    .def("__iter__", [] (ElementRange &er)
      { return py::make_iterator(er.begin(), er.end()); },
      py::keep_alive<0,1>()
    );

  py::class_<FESpace::ElementRange, IntRange> (m, "FESpaceElementRange")
    .def("__iter__", [] (FESpace::ElementRange &er)
      { return py::make_iterator(er.begin(), er.end()); },
      py::keep_alive<0,1>()
    );

  //////////////////////////////////////////////////////////////////////////////////////////

  
  py::enum_<ORDER_POLICY>(m, "ORDER_POLICY", "Enumeration of all supported order policies")
    .value("CONSTANT", CONSTANT_ORDER)
    .value("NODETYPE", NODE_TYPE_ORDER)
    .value("VARIABLE", VARIABLE_ORDER)
    .value("OLDSTYLE", OLDSTYLE_ORDER)
    ;


  //////////////////////////////////////////////////////////////////////////////////////////

  py::class_<FESpace::Element,Ngs_Element>(m, "FESpaceElement")
    .def_property_readonly("dofs",
                           [](FESpace::Element & el) 
                           {
                             // don't use the cached dofs, not cheap anyway
                             // this allows to keep the elements without the iterator
                             // return MakePyList (el.GetDofs()); 
                             Array<DofId> dofs;
                             el.GetFESpace().GetDofNrs(el, dofs);
                             return MakePyList (dofs);
                           },
                           "degrees of freedom of element"
                           )

    
    .def("GetFE",[](FESpace::Element & el)
         {
           // return shared_ptr<FiniteElement>(const_cast<FiniteElement*>(&el.GetFE()), NOOP_Deleter);
           return shared_ptr<FiniteElement> (&el.GetFESpace().GetFE(el, global_alloc));
         },
         // py::return_value_policy::reference,
         "the finite element containing shape functions"
         )

    .def("GetTrafo",[](FESpace::Element & el)
         {
           // return shared_ptr<ElementTransformation>(const_cast<ElementTransformation*>(&el.GetTrafo()), NOOP_Deleter);
           return shared_ptr<ElementTransformation> (&el.GetFESpace().GetMeshAccess()->GetTrafo(el, global_alloc));
         },
         // py::return_value_policy::reference,
         "the transformation from reference element to physical element"
         )

    ;
  //////////////////////////////////////////////////////////////////////////////////////////


  py::class_<GlobalDummyVariables> (m, "GlobalVariables")
    .def_property("msg_level", 
                 &GlobalDummyVariables::GetMsgLevel,
                  &GlobalDummyVariables::SetMsgLevel, "message level")
    .def_property("testout", 
                 &GlobalDummyVariables::GetTestoutFile,
                  &GlobalDummyVariables::SetTestoutFile, "testout file")

    .def_property("pajetrace",
		  &GlobalDummyVariables::GetTestoutFile,
		 [] (GlobalDummyVariables&, int size)
				  {
                                    TaskManager::SetPajeTrace(size > 0);
                                    PajeTrace::SetMaxTracefileSize(size);
                                  })
    .def_property("numthreads",
		 [] (GlobalDummyVariables&)
				  {
                                    return TaskManager::GetMaxThreads ();
                                  },
		 [] (GlobalDummyVariables&, int numthreads)
				  {
                                    TaskManager::SetNumThreads (numthreads);
                                  })

    .def_property("symbolic_integrator_uses_diff",
                  [] (GlobalDummyVariables&)
                  {
                    return symbolic_integrator_uses_diff;
                  },
                  [] (GlobalDummyVariables&, bool val)
                  {
                    symbolic_integrator_uses_diff = val;
                  }, "New treatment of symobolic forms using differentiation by proxies")
                  
    ;

  m.attr("ngsglobals") = py::cast(&globvar);





  //////////////////////////////////////////////////////////////////////////////////////////
  
  py::class_<NGS_Object, shared_ptr<NGS_Object>>(m, "NGS_Object")
    .def_property("name", &NGS_Object::GetName, &NGS_Object::SetName)
    .def_property_readonly("__memory__",
                           [] (const NGS_Object & self)
                           {
                             std::vector<tuple<string,size_t, size_t>> ret;
                             for (auto mui : self.GetMemoryUsage())
                               ret.push_back ( make_tuple(mui.Name(), mui.NBytes(), mui.NBlocks()));
                             return ret;
                           })
    .def_property_readonly("flags", &NGS_Object::GetFlags)
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  typedef shared_ptr<CoefficientFunction> spCF;
  typedef shared_ptr<ProxyFunction> spProxy;

  py::class_<ProxyFunction, spProxy, CoefficientFunction> (m, "ProxyFunction", docu_string(R"raw_string(
Either FESpace.TrialFunction or FESpace.TestFunction. Is a
placeholder coefficient function for Symbolic Integrators. The
integrators will replace it with the basis functions of the finite element space
when building the system matrices.

)raw_string"))
    .def(py::init([](spProxy self) { return self; }))
    .def_property_readonly("space", [](ProxyFunction & self) { return self.GetFESpace(); },
                           "the finite element space")
    .def("Deriv", 
         [](const spProxy self)
         { return self->Deriv(); },
         "take canonical derivative (grad, curl, div)")
    .def("Trace", 
         [](const spProxy self)
         { return self->Trace(); },
         "take canonical boundary trace")
    .def("Other", 
         [](const spProxy self, py::object bnd)
         {
           if (py::extract<double> (bnd).check())
             return self->Other(make_shared<ConstantCoefficientFunction>(py::extract<double> (bnd)()));
           if (py::extract<spCF> (bnd).check())
             return self->Other(py::extract<spCF> (bnd)());
           else
             return self->Other(nullptr);
         },
         "take value from neighbour element (DG)",
         py::arg("bnd") = DummyArgument()
         )
    .def_property_readonly("derivname",
                  [](const spProxy self) -> string
                   {
                     if (!self->Deriv() || !self->DerivEvaluator()) return "";
                     return self->DerivEvaluator()->Name();
                   }, "name of the canonical derivative")
    .def("Operator",
         [] (const spProxy self, string name)
          {
            auto op = self->GetAdditionalProxy(name);
            if (!op)
              throw Exception(string("Operator \"") + name + string("\" does not exist for ") + self->GetFESpace()->GetClassName() + string("!"));

            if (name == "dual")
              op = make_shared<DualProxyFunction> (*op);
            return op;
	  }, py::arg("name"), "Use an additional operator of the finite element space")
    .def("Operators",
         [] (const spProxy self)
         {
           py::list l;
           auto ops = self->GetAdditionalEvaluators();
           for (size_t i = 0; i < ops.Size(); i++)
             l.append (ops.GetName(i));
           return l;
         },"returns list of available differential operators")
    .def("__diffop__", &ProxyFunction::Evaluator)
    ;

  m.def("SetHeapSize",
        [](size_t heapsize)
        {
          if (heapsize > global_heapsize)
            {
              global_heapsize = heapsize;
              glh = LocalHeap (heapsize, "python-comp lh", true);
              lhp.Clear();
            }
        }, py::arg("size"), docu_string(R"raw_string(
Set a new heapsize.

Parameters:

size : int
  input heap size

)raw_string"));
  
  m.def("SetTestoutFile",
        [](string filename)
        {
          delete testout;
          testout = new ofstream (filename);
        }, py::arg("file"), docu_string(R"raw_string(
Enable some logging into file with given filename

Parameters:

file : string
  input file name

)raw_string")
        );

  py::class_<DualProxyFunction, shared_ptr<DualProxyFunction>, ProxyFunction> (m, "DualProxyFunction")
    .def("__call__", [](shared_ptr<DualProxyFunction> self, shared_ptr<CoefficientFunction> u)
         {
           return (*self)(u);
         });
  ;
    

  //////////////////////////////////////////////////////////////////////////////////////////

  ExportArray<COUPLING_TYPE> (m);
  
  auto fes_class = py::class_<FESpace, shared_ptr<FESpace>, NGS_Object>(m, "FESpace",
		    docu_string(R"raw_string(Finite Element Space

Provides the functionality for finite element calculations.

Some available FESpaces are:

H1, HCurl, HDiv, L2, FacetFESpace, HDivDiv

2 __init__ overloads:
  1) To create a registered FESpace
  2) To create a compound FESpace from multiple created FESpaces

1)

Parameters:

type : string
  Type of the finite element space. This parameter is automatically
  set if the space is constructed with a generator function.

mesh : ngsolve.Mesh
  Mesh on which the finite element space is defined on.

kwargs : kwargs
  For a description of the possible kwargs have a look a bit further down.

2)

Parameters:

spaces : list of ngsolve.FESpace
  List of the spaces for the compound finite element space

kwargs : kwargs
  For a description of the possible kwargs have a look a bit further down.

)raw_string"), py::dynamic_attr());
  fes_class
    .def(py::init([fes_class] (py::list lspaces, py::kwargs kwargs)
                  {
		    py::list info;
		    auto flags = CreateFlagsFromKwArgs(kwargs, fes_class, info);
                    Array<shared_ptr<FESpace>> spaces;
                    bool dgjumps = flags.GetDefineFlag("dgjumps");
                    for (auto fes : lspaces )
                    {
                      auto fespace = py::extract<shared_ptr<FESpace>>(fes)(); 
                      if (fespace->UsesDGCoupling())
                        dgjumps = true;
                      spaces.Append(fespace);
                    }
                    flags.SetFlag("dgjumps", dgjumps);
                    if (spaces.Size() == 0)
                      throw Exception("Compound space must have at least one space");

                    bool autoupdate = spaces[0]->DoesAutoUpdate();
                    for (auto space : spaces)
                      {
                        if (space->DoesAutoUpdate() != autoupdate)
                          throw Exception("All spaces must have the same autoupdate setting.");
                      }
                    flags.SetFlag("autoupdate", autoupdate || flags.GetDefineFlag("autoupdate"));

                    int dim = spaces[0]->GetDimension();
                    for (auto space : spaces)
                      if (space->GetDimension() != dim)
                        throw Exception("Compound space of spaces with different dimensions is not allowed");
                    flags.SetFlag ("dim", dim);
                    bool is_complex = spaces[0]->IsComplex() || flags.GetDefineFlag("complex");
                    for (auto space : spaces)
                      if (space->IsComplex() != is_complex)
                        throw Exception("Compound space of spaces with complex and real spaces is not allowed");
                    if (is_complex)
                      flags.SetFlag ("complex");
                    shared_ptr<CompoundFESpace>
                            fes = make_shared<CompoundFESpace> (spaces[0]->GetMeshAccess(), spaces, flags);
                    fes->SetDoSubspaceUpdate(false);
                    fes->Update();
                    fes->FinalizeUpdate();
                    if (!spaces[0]->DoesAutoUpdate())
                        fes->SetDoSubspaceUpdate(true);
                    connect_auto_update(fes.get());
                    return shared_ptr<FESpace>{fes};
                    //                              py::cast(*instance).attr("flags") = bpflags;
                  }),
                  py::arg("spaces"),
                  "construct product space (compound-space) from list of component spaces"
                  )
    .def(py::init([fes_class] (const string & type, shared_ptr<MeshAccess> ma,
                      py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(kwargs, fes_class, info);
                    auto fes = CreateFESpace (type, ma, flags);
                    fes->Update();
                    fes->FinalizeUpdate();
                    connect_auto_update(fes.get());
                    return fes;
                  }),
                  py::arg("type"), py::arg("mesh"),
                  "allowed types are: 'h1ho', 'l2ho', 'hcurlho', 'hdivho' etc."
                  )


    .def_static("__flags_doc__", [] ()
         {
           py::dict flags_doc;
           for (auto & flagdoc : FESpace::GetDocu().arguments)
             flags_doc[get<0> (flagdoc).c_str()] = get<1> (flagdoc);
           return flags_doc;
         })
    .def_static("__special_treated_flags__", [] ()
                {
                  // CAUTION: Do not use references as arguments for cpp function lambdas
                  // pybind11 somehow removes the references and creates copies!
                  // pass arguments either by value or by pointer!
                  py::dict special
                    (
                     py::arg("dirichlet") = py::cpp_function
                     ([] (py::object dirichlet, Flags* flags, py::list info)
                      {
                        auto ma = py::cast<shared_ptr<MeshAccess>>(info[0]);
                        if(py::isinstance<py::list>(dirichlet))
                          {
                            flags->SetFlag("dirichlet",
                                           makeCArray<double>(py::list(dirichlet)));
                            return;
                          }
                        else if(py::isinstance<Region>(dirichlet))
                          {
                            Array<double> dir_indices;
                            auto dir_region = py::cast<Region>(dirichlet);
                            flags->SetFlag("dirichlet", dir_region);
                          }
                        else if (py::isinstance<py::str>(dirichlet))
                          {
                            std::regex pattern(dirichlet.cast<string>());
                            Array<double> dirlist;
                            for (int i = 0; i < ma->GetNBoundaries(); i++)
                              if (std::regex_match (ma->GetMaterial(BND, i), pattern))
                                {
                                  dirlist.Append (i+1);
                                }
                            flags->SetFlag("dirichlet", dirlist);
                          }
                        else
                          throw py::type_error("dirichlet parameter has wrong type!");
                      }),
                     py::arg("dirichlet_bbnd") = py::cpp_function
                     ([](py::object dirichlet_bbnd, Flags* flags, py::list info)
                     {
                       if(py::isinstance<py::str>(dirichlet_bbnd))
                         flags->SetFlag("dirichlet_bbnd", py::cast<std::string>(dirichlet_bbnd));
                       else if(py::isinstance<Region>(dirichlet_bbnd))
                         flags->SetFlag("dirichlet_bbnd", py::cast<Region>(dirichlet_bbnd));
                       else
                         throw py::type_error("dirichlet_bbnd has wrong type!");
                     }),
                     py::arg("dirichlet_bbbnd") = py::cpp_function
                     ([](py::object dirichlet_bbbnd, Flags* flags, py::list info)
                     {
                       if(py::isinstance<py::str>(dirichlet_bbbnd))
                         flags->SetFlag("dirichlet_bbbnd", py::cast<std::string>(dirichlet_bbbnd));
                       else if(py::isinstance<Region>(dirichlet_bbbnd))
                         flags->SetFlag("dirichlet_bbbnd", py::cast<Region>(dirichlet_bbbnd));
                       else
                         throw py::type_error("dirichlet_bbbnd has wrong type!");
                     }),
                     py::arg("definedon") = py::cpp_function
                     ([] (py::object definedon, Flags* flags, py::list info)
                      {
                        auto ma = py::cast<shared_ptr<MeshAccess>>(info[0]);
                        if (py::isinstance<py::str>(definedon))
                          {
                            std::regex pattern(definedon.cast<string>());
                            Array<double> defonlist;
                            for (int i = 0; i < ma->GetNDomains(); i++)
                              if (regex_match(ma->GetMaterial(VOL,i), pattern))
                                defonlist.Append(i+1);
                            flags->SetFlag ("definedon", defonlist);
                          }

                        if (py::isinstance<py::list> (definedon))
                          flags->SetFlag ("definedon", makeCArray<double> (definedon));

                        try
                          {
                            auto reg = py::cast<Region>(definedon);
                            flags->SetFlag("definedon", reg);
                          }
                        catch(py::cast_error)
                          {}
                        try
                          {
                            auto map = py::cast<std::map<VorB, Region>>(definedon);
                            flags->SetFlag("definedon", map);
                          }
                        catch(py::cast_error)
                          {}
                        // py::extract<Region> definedon_reg(definedon);
                        // if (definedon_reg.check() && definedon_reg().IsVolume())
                        //   {
                        //     Array<double> defonlist;
                        //     for (int i = 0; i < definedon_reg().Mask().Size(); i++)
                        //       if (definedon_reg().Mask().Test(i))
                        //         defonlist.Append(i+1);
                        //     flags->SetFlag ("definedon", defonlist);
                        //   }
                        // if (definedon_reg.check() && definedon_reg().VB()==BND)
                        //   {
                        //     Array<double> defonlist;
                        //     flags->SetFlag("definedon", defonlist); //empty
                        //     for (auto i : Range(definedon_reg().Mask().Size()))
                        //       if(definedon_reg().Mask().Test(i))
                        //         defonlist.Append(i+1);
                        //     flags->SetFlag("definedonbound", defonlist);
                        //   }
                      }),
                     py::arg("order_policy") = py::cpp_function
                     ([] (ORDER_POLICY op, Flags* flags, py::list info)
                      {
                        flags->SetFlag("order_policy", int(op));
                      })
                     );
                     return special;
                     })
    .def(py::pickle(&fesPickle, (shared_ptr<FESpace>(*)(py::tuple)) fesUnpickle<FESpace>))
    .def("Update", [](shared_ptr<FESpace> self)
         { 
           self->Update();
           self->FinalizeUpdate();
         },
         "update space after mesh-refinement")
     .def("UpdateDofTables", [](shared_ptr<FESpace> self)
         {
           self->UpdateDofTables();
           self->UpdateCouplingDofArray();
           self->FinalizeUpdate();
         },
         "update dof-tables after changing polynomial order distribution")
     .def("FinalizeUpdate", [](shared_ptr<FESpace> self)
         { 
           self->FinalizeUpdate();
         },
         "finalize update")
     .def("HideAllDofs", [](shared_ptr<FESpace> self, py::object acomp)
         {
           shared_ptr<FESpace> space = self;
           if (! py::extract<DummyArgument> (acomp).check())
           {
             auto comp = py::extract<int>(acomp)();
             auto compspace = dynamic_pointer_cast<CompoundFESpace> (self);
             if (!compspace)
               throw py::type_error("'components' is available only for product spaces");
             space = (*compspace)[comp];
             IntRange range = compspace->GetRange(comp);
             for (auto d : range)
             {
               auto doftype = compspace->GetDofCouplingType(d);
               if (doftype != UNUSED_DOF)
                 compspace->SetDofCouplingType(d,HIDDEN_DOF);
             }
           }
           for (DofId d : Range(space->GetNDof()))
             {
               auto doftype = space->GetDofCouplingType(d);
               if (doftype != UNUSED_DOF)
                 space->SetDofCouplingType(d,HIDDEN_DOF);
             }
           self->FinalizeUpdate(); //Update FreeDofs
         }, py::arg("component")=DummyArgument(), 
         "set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())")

    .def_property_readonly("components", 
                  [](shared_ptr<FESpace> self)-> py::tuple
                   { 
                     if (auto compspace = dynamic_pointer_cast<CompoundFESpace>(self))
                       {
                         py::tuple vecs(compspace->GetNSpaces());
                         for (int i = 0; i < compspace -> GetNSpaces(); i++) 
                           vecs[i]= py::cast((*compspace)[i]);
                         return vecs;
                       }
                     throw Exception("components only available for ProductSpace");                     
                   }, "deprecated, will be only available for ProductSpace")
    .def("Range",
         [] (shared_ptr<FESpace> self, int comp)
         {
           if (auto compspace = dynamic_pointer_cast<CompoundFESpace>(self))
             return compspace->GetRange(comp);
           throw Exception("Range only available for ProductSpace");
         }, "deprecated, will be only available for ProductSpace")

    
    .def_property_readonly ("ndof", [](shared_ptr<FESpace> self) { return self->GetNDof(); },
                            "number of degrees of freedom")

    .def_property_readonly ("ndofglobal",
                            [](shared_ptr<FESpace> self) { return self->GetNDofGlobal(); },
                            "global number of dofs on MPI-distributed mesh")
    .def_property_readonly ("dim",
                            [](shared_ptr<FESpace> self) { return self->GetDimension(); },
                            "multi-dim of FESpace")
    .def("__str__", [] (shared_ptr<FESpace> self) { return ToString(*self); } )
    .def("__timing__", [] (shared_ptr<FESpace> self) { return py::cast(self->Timing()); })
    .def_property_readonly("lospace", [](shared_ptr<FESpace> self) -> shared_ptr<FESpace>
			   { return self->LowOrderFESpacePtr(); })
    .def_property_readonly("loembedding", [](shared_ptr<FESpace> self) -> shared_ptr<BaseMatrix>
			   { return self->LowOrderEmbedding(); })
    .def_property_readonly("mesh",
                           [](shared_ptr<FESpace> self) -> shared_ptr<MeshAccess>
                           { return self->GetMeshAccess(); }, "mesh on which the FESpace is created")

    .def_property_readonly("globalorder", [] (shared_ptr<FESpace> self) { return self->GetOrder(); },
                  "query global order of space")    
    .def_property_readonly("type", [] (shared_ptr<FESpace> self) { return self->type; },
                  "type of finite element space")

    .def_property_readonly("is_complex", &FESpace::IsComplex)

    .def("SetDefinedOn", [] (FESpace& self, Region& reg)
         {
           self.SetDefinedOn(reg.VB(),reg.Mask());
         }, py::arg("region"), docu_string(R"raw_string(
Set the regions on which the FESpace is defined.

Parameters:

region : ngsolve.comp.Region
  input region

)raw_string"))

    .def("SetOrder",
         [](shared_ptr<FESpace> self, ELEMENT_TYPE et, int order)
         {
           self->SetOrder (et, order);
         },
         py::arg("element_type"),
         py::arg("order"), docu_string(R"raw_string(

Parameters:

element_type : ngsolve.fem.ET
  input element type

order : object
  input polynomial order
)raw_string")
         )

    .def("SetOrder",
         [](shared_ptr<FESpace> self, NodeId ni, int order)
         {
           self->SetOrder(ni, order);
         },
         py::arg("nodeid"),
         py::arg("order"), docu_string(R"raw_string(

Parameters:

nodeid : ngsolve.comp.NodeId
  input node id

order : int
  input polynomial order

)raw_string")
         )

    .def("GetOrder",
         [](shared_ptr<FESpace> self, NodeId ni) -> int
         {
           return self->GetOrder(ni);
         },
         py::arg("nodeid"),
         "return order of node.\n"
         "by now, only isotropic order is supported here\n")
    
    .def("Elements", 
         [](shared_ptr<FESpace> self, VorB vb)
         { return FESpace::ElementRange(self->Elements(vb, glh)); },
         py::arg("VOL_or_BND")=VOL, docu_string(R"raw_string(
Returns an iterable range of elements.

Parameters:

VOL_or_BND : ngsolve.comp.VorB
  input VOL, BND, BBND,...

)raw_string"))

    .def("GetDofNrs", [](shared_ptr<FESpace> self, ElementId ei)
         {
           Array<DofId> tmp; self->GetDofNrs(ei,tmp);
           return MakePyTuple(tmp);
         }, py::arg("ei"), docu_string(R"raw_string(

Parameters:

ei : ngsolve.comp.ElementId
  input element id

)raw_string"))

    .def("GetDofNrs", [](shared_ptr<FESpace> self, NodeId ni)
         {
           Array<DofId> tmp; self->GetDofNrs(ni,tmp);
           return MakePyTuple(tmp);
         }, py::arg("ni"), docu_string(R"raw_string(

Parameters:

ni : ngsolve.comp.NodeId
  input node id

)raw_string"))

    .def ("GetDofs", [](shared_ptr<FESpace> self, Region reg)
          {
            return self->GetDofs(reg);
          }, py::arg("region"), docu_string(R"raw_string(
Returns all degrees of freedom in given region.

Parameters:

region : ngsolve.comp.Region
  input region

)raw_string"))
    
    .def("CouplingType", [](shared_ptr<FESpace> self, DofId dofnr) -> COUPLING_TYPE
         { return self->GetDofCouplingType(dofnr); },
         py::arg("dofnr"), docu_string(R"raw_string(
         Get coupling type of a degree of freedom.

Parameters:

dofnr : int
  input dof number

)raw_string")
         )
    .def("SetCouplingType", [](shared_ptr<FESpace> self, DofId dofnr, COUPLING_TYPE ct)
         { self->SetDofCouplingType(dofnr,ct); },
         py::arg("dofnr"), py::arg("coupling_type"), docu_string(R"raw_string(
         Set coupling type of a degree of freedom.

Parameters:

dofnr : int
  input dof number

coupling_type : ngsolve.comp.COUPLING_TYPE
  input coupling type

)raw_string")
         )

    .def("SetCouplingType", [](shared_ptr<FESpace> self, IntRange dofnrs, COUPLING_TYPE ct)
         {
           for (auto d : dofnrs)
             self->SetDofCouplingType(DofId(d),ct);
         },
         py::arg("dofnrs"), py::arg("coupling_type"), docu_string(R"raw_string(
         Set coupling type for interval of dofs.

Parameters:

dofnrs : Range
  range of dofs

coupling_type : ngsolve.comp.COUPLING_TYPE
  input coupling type

)raw_string")
         )

    .def_property_readonly("couplingtype", [] (shared_ptr<FESpace> self)
                           { return FlatArray<COUPLING_TYPE>(self->CouplingTypes()); })
    .def_property_readonly("autoupdate", [] (shared_ptr<FESpace> self) {return self->DoesAutoUpdate();})
    .def ("GetFE", [](shared_ptr<FESpace> self, ElementId ei) -> py::object
          {
            auto fe = shared_ptr<FiniteElement> (&self->GetFE(ei, global_alloc));
            
            auto scalfe = dynamic_pointer_cast<BaseScalarFiniteElement> (fe);
            if (scalfe) return py::cast(scalfe);

            auto hcurlfe = dynamic_pointer_cast<BaseHCurlFiniteElement> (fe);
            if (hcurlfe) return py::cast(hcurlfe);

            auto hdivfe = dynamic_pointer_cast<BaseHDivFiniteElement> (fe);
            if (hdivfe) return py::cast(hdivfe);

            auto hdivdivfe = dynamic_pointer_cast<BaseHDivDivFiniteElement> (fe);
            if (hdivdivfe) return py::cast(hdivdivfe);

            return py::cast(fe);
          }, py::arg("ei"), docu_string(R"raw_string(
Get the finite element to corresponding element id.

Parameters:

ei : ngsolve.comp.ElementId
   input element id

)raw_string"))
    
    .def("FreeDofs",
         [] (const shared_ptr<FESpace>self, bool coupling)
         { return self->GetFreeDofs(coupling); },
         py::arg("coupling")=false,docu_string(R"raw_string(

Return BitArray of free (non-Dirichlet) dofs\n
coupling=False ... all free dofs including local dofs\n
coupling=True .... only element-boundary free dofs

Parameters:

coupling : bool
  input coupling

)raw_string")
         )

    .def("ParallelDofs",
         [] (const shared_ptr<FESpace> self)
         { return self->GetParallelDofs(); },
         "Return dof-identification for MPI-distributed meshes")

    .def("Prolongation",
         [] (const shared_ptr<FESpace> self)
         { return self->GetProlongation(); },
         "Return prolongation operator for use in multi-grid")

    .def("TrialFunction",
         [] (const shared_ptr<FESpace> self)
         {
           return MakeProxyFunction (self, false);
         },
         docu_string("Return a proxy to be used as a trialfunction in :any:`Symbolic Integrators<symbolic-integrators>`"))
    
    .def("TestFunction",
         [] (const shared_ptr<FESpace> self)
           {
             return MakeProxyFunction (self, true);
           },
         docu_string("Return a proxy to be used as a testfunction for :any:`Symbolic Integrators<symbolic-integrators>`"))

    .def("TnT",
         [] (const shared_ptr<FESpace> self)
         {
           return std::make_tuple(MakeProxyFunction (self, false), MakeProxyFunction (self, true));
         },
         docu_string("Return a tuple of trial and testfunction"))

    .def("InvM",
         [] (const shared_ptr<FESpace> self,
             shared_ptr<CoefficientFunction> rho) -> shared_ptr<BaseMatrix>
         {
           return make_shared<ApplyMass> (self, rho, true, nullptr, glh); 
         }, py::arg("rho") = nullptr)
    .def("Mass",
         [] (const shared_ptr<FESpace> self,
             shared_ptr<CoefficientFunction> rho,
             optional<Region> definedon) -> shared_ptr<BaseMatrix>
         {
           shared_ptr<Region> spdefon;
           if (definedon) spdefon = make_shared<Region> (*definedon);
           // return make_shared<ApplyMass> (self, rho, false, spdefon, glh);
           return self->GetMassOperator(rho, spdefon, glh);
         }, py::arg("rho") = nullptr, py::arg("definedon") = nullptr)
    
    .def("SolveM",
         [] (const shared_ptr<FESpace> self,
             BaseVector& vec, spCF rho, Region * definedon) 
         { self->SolveM(rho.get(), vec, definedon, glh); },
         py::arg("vec"), py::arg("rho")=nullptr, py::arg("definedon")=nullptr, docu_string(R"raw_string(
         Solve with the mass-matrix. Available only for L2-like spaces.

Parameters:

vec : ngsolve.la.BaseVector
  input right hand side vector

rho : ngsolve.fem.CoefficientFunction
  input CF

)raw_string"))
    .def("ApplyM",
         [] (const shared_ptr<FESpace> self,
             BaseVector& vec, spCF rho, Region * definedon)
         { self->ApplyM(rho.get(), vec, definedon, glh); },
         py::arg("vec"), py::arg("rho")=nullptr, py::arg("definedon")=nullptr,
         "Apply mass-matrix. Available only for L2-like spaces")
    .def ("TraceOperator", [] (shared_ptr<FESpace> self, shared_ptr<FESpace> tracespace,
                               bool avg) -> shared_ptr<BaseMatrix>
          {
            return self->GetTraceOperator(tracespace, avg);
            // return make_shared<ApplyTrace> (self, tracespace, avg, glh);             
          }, py::arg("tracespace"), py::arg("average"))
    .def ("ConvertL2Operator", [] (shared_ptr<FESpace> self, shared_ptr<FESpace> l2space)
          {
            return self->ConvertL2Operator(l2space);
          }, py::arg("l2space"))
    .def ("GetTrace", [] (shared_ptr<FESpace> self, const FESpace & tracespace,
                          BaseVector & in, BaseVector & out, bool avg)
          {
            self->GetTrace(tracespace, in, out, avg, glh);
          })
    .def ("GetTraceTrans", [] (shared_ptr<FESpace> self, const FESpace & tracespace,
                               BaseVector & in, BaseVector & out, bool avg)
          {
            self->GetTraceTrans(tracespace, in, out, avg, glh);
          })
    
    .def("__eq__",
         [] (shared_ptr<FESpace> self, shared_ptr<FESpace> other)
         {
           return self == other;
         }, py::arg("space"))
    
    .def("__mul__", [] (shared_ptr<FESpace> space1, shared_ptr<FESpace> space2) {
        bool is_complex = space1->IsComplex();
        if (is_complex != space2->IsComplex())
          throw Exception("Product space must consist of all real or all complex spaces");
        int dim = space1->GetDimension();
        if (dim != space2->GetDimension())
          throw Exception("Product space needs same dimension for all components");
        Flags flags;
        if (is_complex) flags.SetFlag("complex");
        flags.SetFlag ("dim", dim);
        flags.SetFlag ("dgjumps", space1->UsesDGCoupling() || space2->UsesDGCoupling());

        if (space1->DoesAutoUpdate() != space2->DoesAutoUpdate())
          throw Exception("Both spaces need same autoupdate setting.");
        flags.SetFlag ("autoupdate", space1->DoesAutoUpdate());

        if(space1->LowOrderFESpacePtr() && space2->LowOrderFESpacePtr())
          flags.SetFlag("low_order_space");
        auto productspace = make_shared<CompoundFESpace> (space1->GetMeshAccess(), flags);

        for (auto s : { space1, space2 })
          {
            const FESpace & fes = *s; // avoid side effect warning
            if (typeid(fes) == typeid(CompoundFESpace)) // exactly Compound and not more
              for (auto spc : dynamic_pointer_cast<CompoundFESpace> (s) -> Spaces())
                productspace->AddSpace(spc);
            else
              productspace->AddSpace (s);
          }
        productspace->SetDoSubspaceUpdate(false);
        productspace->Update();
        productspace->FinalizeUpdate();
        if (!space1->DoesAutoUpdate())
          productspace->SetDoSubspaceUpdate(true);
        connect_auto_update(productspace.get());
        return productspace;
      })

    .def("__pow__", [] (shared_ptr<FESpace> space, int p) {
        bool is_complex = space->IsComplex();
        int dim = space->GetDimension();
        Flags flags;
        if (is_complex) flags.SetFlag("complex");
        flags.SetFlag ("dim", dim);
        flags.SetFlag ("autoupdate", space->DoesAutoUpdate());
        auto productspace = make_shared<CompoundFESpaceAllSame> (space, p, flags);
        productspace->SetDoSubspaceUpdate(false);
        productspace->Update();
        productspace->FinalizeUpdate();
        if (!space->DoesAutoUpdate())
          productspace->SetDoSubspaceUpdate(true);
        connect_auto_update(productspace.get());
        return productspace;
      })

    // not working because shared_ptr<Array<int>> cannot be pybind return type?
    // TODO CL: Find out why...
    // .def("CreateDirectSolverCluster", &FESpace::CreateDirectSolverClusters)
    .def("CreateDirectSolverCluster", [](FESpace& self, py::kwargs kwargs)
    {
      auto flags = CreateFlagsFromKwArgs(kwargs);
      auto cluster = self.CreateDirectSolverClusters(flags);
      py::list pycluster(self.GetNDof());
      if(cluster)
        for(auto i : Range(self.GetNDof()))
          pycluster[i] = (*cluster)[i];
      else
        for(auto i : Range(self.GetNDof()))
          pycluster[i] = 0;
      return pycluster;
    })
    ;

  py::class_<CompoundFESpace, shared_ptr<CompoundFESpace>, FESpace>
    (m,"ProductSpace")
    .def(py::init([] (py::args pyspaces) {

          /*
          Array<shared_ptr<FESpace>> spaces;
          for (auto fes : pyspaces )
            spaces.Append(py::extract<shared_ptr<FESpace>>(fes)());
          */
          auto spaces = makeCArray<shared_ptr<FESpace>> (pyspaces);
            
          Flags flags;
          if (spaces.Size() == 0)
            throw Exception("Product space must have at least one space");
          int dim = spaces[0]->GetDimension();
          for (auto space : spaces)
            if (space->GetDimension() != dim)
              throw Exception("Product space of spaces with different dimensions is not allowed");
          flags.SetFlag ("dim", dim);
          bool is_complex = spaces[0]->IsComplex();
          for (auto space : spaces)
            if (space->IsComplex() != is_complex)
              throw Exception("Product space of spaces with complex and real spaces is not allowed");
          if (is_complex)
            flags.SetFlag ("complex");

          bool autoupdate = spaces[0]->DoesAutoUpdate();
          for (auto space : spaces)
            if (space->DoesAutoUpdate() != autoupdate)
                throw Exception("All spaces must have the same autoupdate setting.");

          auto fes = make_shared<CompoundFESpace> (spaces[0]->GetMeshAccess(), spaces, flags);
          fes->SetDoSubspaceUpdate(false);
          fes->Update();
          fes->FinalizeUpdate();
          if (!autoupdate)
            fes->SetDoSubspaceUpdate(true);
          connect_auto_update(fes.get());
          return fes;
        }))
    .def(py::pickle([] (py::object pyfes)
                    {
                      auto fes = py::cast<shared_ptr<CompoundFESpace>>(pyfes);
                      auto flags = fes->GetFlags();
                      py::list lst;
                      for(auto i : Range(fes->GetNSpaces()))
                        lst.append((*fes)[i]);
                      return py::make_tuple(lst,flags,pyfes.attr("__dict__"));
                    },
                    [] (py::tuple state)
                    {
                      Array<shared_ptr<FESpace>> spaces;
                      for (auto pyfes : state[0].cast<py::list>())
                        spaces.Append(pyfes.cast<shared_ptr<FESpace>>());
                      auto fes = make_shared<CompoundFESpace>
                        (spaces[0]->GetMeshAccess(), spaces, state[1].cast<Flags>());
                      LocalHeap lh (1000000, "FESpace::Update-heap");
                      fes->Update();
                      fes->FinalizeUpdate();
                      py::cast(fes).attr("__dict__") = state[2];
                      return fes;
                    }))

    .def("Range",
         [] (shared_ptr<CompoundFESpace> self, int comp)
         {
           return self->GetRange(comp);
         },
         py::arg("component"), docu_string(R"raw_string(
         Return interval of dofs of a component of a product space.

Parameters:

component : int
  input component

)raw_string"))

    .def("Embedding", [] (shared_ptr<CompoundFESpace> self, int comp)
         {
           return self->EmbeddingOperator(comp);
           // return make_shared<Embedding> (self->GetNDof(), self->GetRange(comp), self->IsComplex());
         },
         py::arg("component"), "create embedding operator for this component")

    .def("Restriction", [] (shared_ptr<CompoundFESpace> self, int comp)
         {
           return self->RestrictionOperator(comp);
           // return make_shared<Embedding> (self->GetNDof(), self->GetRange(comp), self->IsComplex());
         },
         py::arg("component"), "create restriction operator onto this component")

    
    .def_property_readonly("embeddings", 
                  [](shared_ptr<CompoundFESpace> self)-> py::list
                   { 
                     py::list embeddings(self->GetNSpaces());
                     for (int i = 0; i < self -> GetNSpaces(); i++)
                       embeddings[i] = self->EmbeddingOperator(i);                       
                       /*
                       embeddings[i]= py::cast(make_shared<Embedding> (self->GetNDof(), self->GetRange(i),
                                                                       self->IsComplex()));
                       */
                     return embeddings;
                   },
                  "returns a list of embeddings for the component spaces")

    .def_property_readonly("restrictions", 
                  [](shared_ptr<CompoundFESpace> self)-> py::list
                   { 
                     py::list restrictions(self->GetNSpaces());
                     for (int i = 0; i < self -> GetNSpaces(); i++)
                       restrictions[i] = self->RestrictionOperator(i);                       
                     return restrictions;
                   },
                  "returns a list of restrictions onto the component spaces")
    
    .def_property_readonly("components", 
                           [](shared_ptr<CompoundFESpace> self)
                           {
                             return Array<shared_ptr<FESpace>>(self->Spaces());
                           },
                  "Return a list of the components of a product space")
    .def("SetDoSubspaceUpdate", &CompoundFESpace::SetDoSubspaceUpdate)
    ;

  py::class_<CompoundFESpaceAllSame, shared_ptr<CompoundFESpaceAllSame>, CompoundFESpace>
    (m,"VectorValued")
    .def(py::init([] (shared_ptr<FESpace> space, optional<int> optdim, bool autoupdate) {
          Flags flags;
          flags.SetFlag("autoupdate", autoupdate || space->DoesAutoUpdate());
          int sdim = optdim.value_or (space->GetSpatialDimension());
          auto vecspace = make_shared<CompoundFESpaceAllSame> (space, sdim, flags);
          vecspace->SetDoSubspaceUpdate(false);
          vecspace->Update();
          vecspace->FinalizeUpdate();
          if (!space->DoesAutoUpdate())
            vecspace->SetDoSubspaceUpdate(true);
          connect_auto_update(vecspace.get());
          return vecspace;
        }),py::arg("space"), py::arg("dim")=nullopt, py::arg("autoupdate")=false)
    .def(py::pickle([] (py::object pyfes)
                    {
                      auto fes = py::cast<shared_ptr<CompoundFESpaceAllSame>>(pyfes);
                      auto flags = fes->GetFlags();
                      return py::make_tuple((*fes)[0], fes->GetNSpaces(), flags, pyfes.attr("__dict__"));
                    },
                    [] (py::tuple state)
                    {
                      shared_ptr<FESpace> spc1 = state[0].cast<shared_ptr<FESpace>>();
                      int dim = state[1].cast<int>();
                      auto fes = make_shared<CompoundFESpaceAllSame>(spc1, dim, state[2].cast<Flags>());
                      LocalHeap lh (1000000, "FESpace::Update-heap");
                      fes->Update();
                      fes->FinalizeUpdate();
                      py::cast(fes).attr("__dict__") = state[3];
                      return fes;
                    }))
    
    ;

  py::class_<MatrixFESpace, shared_ptr<MatrixFESpace>, CompoundFESpace>
    (m,"MatrixValued")
    .def(py::init([] (shared_ptr<FESpace> space, optional<int> optdim, bool symmetric, bool deviatoric,
            bool autoupdate) {
          Flags flags;
          if (symmetric) flags.SetFlag("symmetric");
          if (deviatoric) flags.SetFlag("deviatoric");
          flags.SetFlag("autoupdate", autoupdate || space->DoesAutoUpdate());
          int sdim = optdim.value_or (space->GetSpatialDimension());
          auto matspace = make_shared<MatrixFESpace> (space, sdim, flags);
          matspace->SetDoSubspaceUpdate(false);
          matspace->Update();
          matspace->FinalizeUpdate();
          if (!space->DoesAutoUpdate())
            matspace->SetDoSubspaceUpdate(true);
          connect_auto_update(matspace.get());
          return matspace;
        }), py::arg("space"), py::arg("dim")=nullopt, py::arg("symmetric")=false, py::arg("deviatoric")=false,
        py::arg("autoupdate")=false)
    ;
  
  ExportFESpace<HCurlHighOrderFESpace> (m, "HCurl")
    .def("CreateGradient", [](shared_ptr<HCurlHighOrderFESpace> self) {
        auto fesh1 = self->CreateGradientSpace();
        shared_ptr<BaseMatrix> grad = self->CreateGradient(*fesh1);
        return py::make_tuple(grad, shared_ptr<FESpace>(fesh1));
      })
    ;
  
  ExportFESpace<HDivHighOrderFESpace> (m, "HDiv")
    .def("Average", &HDivHighOrderFESpace::Average,
         py::arg("vector"))
    ;
  
  ExportFESpace<H1HighOrderFESpace> (m, "H1");

  ExportFESpace<VectorH1FESpace, CompoundFESpace> (m, "VectorH1");
 
  ExportFESpace<VectorL2FESpace, CompoundFESpace> (m, "VectorL2");

  ExportFESpace<L2SurfaceHighOrderFESpace> (m, "SurfaceL2");

  ExportFESpace<TangentialSurfaceL2FESpace> (m, "TangentialSurfaceL2");

  ExportFESpace<NumberFESpace> (m, "NumberSpace");
  ExportFESpace<IntegrationRuleSpace> (m, "IntegrationRuleSpace")
    .def("GetIntegrationRules", &IntegrationRuleSpace::GetIntegrationRules)
    ;
  ExportFESpace<IntegrationRuleSpaceSurface> (m, "IntegrationRuleSpaceSurface")
    .def("GetIntegrationRules", &IntegrationRuleSpaceSurface::GetIntegrationRules)
    ;
  ExportFESpace<H1LumpingFESpace> (m, "H1LumpingFESpace")
    .def("GetIntegrationRules", &H1LumpingFESpace::GetIntegrationRules)
    ;
    

  ExportFESpace<L2HighOrderFESpace> (m, "L2");

  ExportFESpace<HDivDivFESpace> (m, "HDivDiv");
  
  ExportFESpace<HCurlDivFESpace> (m, "HCurlDiv");

  ExportFESpace<HCurlCurlFESpace> (m, "HCurlCurl");
  
  ExportFESpace<HDivDivSurfaceSpace> (m, "HDivDivSurface");
  
  // ExportFESpace<VectorFacetFESpace> (m, "VectorFacet");
  ExportFESpace<VectorFacetFESpace> (m, "TangentialFacetFESpace");
  ExportFESpace<NormalFacetFESpace> (m, "NormalFacetFESpace");

  ExportFESpace<FacetFESpace> (m, "FacetFESpace");
  
  ExportFESpace<FacetSurfaceFESpace> (m, "FacetSurface");
  ExportFESpace<NormalFacetSurfaceFESpace> (m, "NormalFacetSurface");
  
  ExportFESpace<HDivHighOrderSurfaceFESpace> (m, "HDivSurface")
    .def("Average", &HDivHighOrderSurfaceFESpace::Average,
         py::arg("vector"))
    ;

  ExportFESpace<VectorFESpace<L2SurfaceHighOrderFESpace>> (m, "VectorSurfaceL2");
  ExportFESpace<VectorFESpace<FacetFESpace>> (m, "VectorFacetFESpace");
  ExportFESpace<VectorFESpace<FacetSurfaceFESpace>> (m, "VectorFacetSurface");


  ExportFESpace<NodalFESpace> (m, "NodalFESpace");  
  ExportFESpace<VectorFESpace<NodalFESpace>> (m, "VectorNodalFESpace");
  
  // py::class_<CompoundFESpace, shared_ptr<CompoundFESpace>, FESpace>
  //   (m, "CompoundFESpace")
  //   .def("Range", &CompoundFESpace::GetRange)
  //   ;

  auto export_global = ExportFESpace<GlobalSpace> (m, "GlobalSpace");
  export_global.def("AddOperator", [](shared_ptr<GlobalSpace> space, string name,
                                      VorB vb, shared_ptr<CoefficientFunction> dbasis)
                    {
                      space->AddOperator(name, vb, dbasis);
                    });
  export_global.def_static("__special_treated_flags__", [fes_class] ()
                           {
                             py::dict special = fes_class.attr("__special_treated_flags__")();
                             special["basis"] = py::cpp_function
                               ([] (py::object pybasis, Flags* flags, py::list info)
                                {
                                  auto cppbasis = py::cast<shared_ptr<CoefficientFunction>>(pybasis);
                                  flags -> SetFlag("basis", std::any(cppbasis));
                                });
                             return special;
                           });


  
  py::class_<PeriodicFESpace, shared_ptr<PeriodicFESpace>, FESpace>(m, "Periodic",
	docu_string(R"delimiter(Periodic or quasi-periodic Finite Element Spaces.
The periodic fespace is a wrapper around a standard fespace with an 
additional dof mapping for the periodic degrees of freedom. All dofs 
on minion boundaries are mapped to their master dofs. Because of this, 
the mesh needs to be periodic. Low order fespaces are currently not
supported, so methods using them will not work.

Parameters:

fespace : ngsolve.comp.FESpace
    finite element space

phase : list of Complex = None
    phase shift for quasi-periodic finite element space. The basis
    functions on the minion boundary are multiplied by the factor
    given in this list. If None (default) is given, a periodic
    fespace is created. The order of the list must match the order
    of the definition of the periodic boundaries in the mesh.

used_idnrs : list of int = None
    identification numbers to be made periodic if you don't want to
    use all periodic identifications defined in the mesh, if None
    (default) all available periodic identifications are used.

)delimiter"))
    .def(py::init([] (shared_ptr<FESpace> & fes,
                      optional<py::list> phase, py::object use_idnrs, bool autoupdate)
                  {
                    Flags flags = fes->GetFlags();
                    flags.SetFlag("autoupdate", autoupdate || fes->DoesAutoUpdate());
                    shared_ptr<Array<int>> a_used_idnrs;
                    if(py::extract<py::list>(use_idnrs).check())
                      a_used_idnrs = make_shared<Array<int>>(makeCArray<int>(py::extract<py::list>(use_idnrs)()));
                    else
                      throw Exception("Argument for use_idnrs in Periodic must be list of identification numbers (int)");
                    shared_ptr<PeriodicFESpace> perfes;
                    if(phase.has_value() && py::len(*phase) > 0)
                      {
                        auto lphase = *phase;
                        py::extract<double> ed(lphase[0]);
                        if(ed.check())
                          {
                            auto a_phase = make_shared<Array<double>>(py::len(*phase));
                            for (auto i : Range(a_phase->Size()))
                              (*a_phase)[i] = py::cast<double>(lphase[i]);
                            perfes = make_shared<QuasiPeriodicFESpace<double>>(fes,flags,a_used_idnrs,a_phase);
                          }
                        else
                          {
                            auto a_phase = make_shared<Array<Complex>>(py::len(*phase));
                            for (auto i : Range(a_phase->Size()))
                              (*a_phase)[i] = py::cast<Complex>(lphase[i]);
                            perfes = make_shared<QuasiPeriodicFESpace<Complex>>(fes,flags,a_used_idnrs,a_phase);
                          }
                      }
                    else
                      perfes = make_shared<PeriodicFESpace>(fes,flags,a_used_idnrs);
                    // MR: Note that PeriodicFESpace always updates the wrapped space
                    perfes->Update();
                    perfes->FinalizeUpdate();
                    connect_auto_update(perfes.get());
                    return perfes;
                  }), py::arg("fespace"), py::arg("phase")=nullopt,
                  py::arg("use_idnrs")=py::list(), py::arg("autoupdate")=false)
    .def(py::pickle([](const PeriodicFESpace* per_fes)
                    {
                      py::list idnrs;
                      for (auto idnr : *per_fes->GetUsedIdnrs())
                        idnrs.append(idnr);
                      auto quasiper_fes_d = dynamic_cast<const QuasiPeriodicFESpace<double>*>(per_fes);
                      if(quasiper_fes_d)
                        {
                          py::list fac;
                          for(auto factor : *quasiper_fes_d->GetFactors())
                            fac.append(factor);
                          return py::make_tuple(per_fes->GetBaseSpace(),idnrs,fac);
                        }
                      auto quasiper_fes_c = dynamic_cast<const QuasiPeriodicFESpace<Complex>*>(per_fes);
                      if(quasiper_fes_c)
                        {
                          py::list fac;
                          for(auto factor : *quasiper_fes_c->GetFactors())
                            fac.append(factor);
                          return py::make_tuple(per_fes->GetBaseSpace(),idnrs,fac);
                        }
                      return py::make_tuple(per_fes->GetBaseSpace(),idnrs);
                    },
                    [] (py::tuple state) -> shared_ptr<PeriodicFESpace>
                    {
                      shared_ptr<PeriodicFESpace> fes;
                      auto idnrs = make_shared<Array<int>>();
                      for (auto id : state[1].cast<py::list>())
                        idnrs->Append(id.cast<int>());
                      if(py::len(state)==3)
                        {
                          auto pyfacs = state[2].cast<py::list>();
                          if(py::extract<double>(pyfacs[0]).check())
                            {
                              auto facs = make_shared<Array<double>>();
                              for(auto fac : pyfacs)
                                facs->Append(fac.cast<double>());
                              fes = make_shared<QuasiPeriodicFESpace<double>>
                                (state[0].cast<shared_ptr<FESpace>>(), Flags(),idnrs,facs);
                            }
                          else
                            {
                              auto facs = make_shared<Array<Complex>>();
                              for (auto fac : pyfacs)
                                facs->Append(fac.cast<Complex>());
                              fes = make_shared<QuasiPeriodicFESpace<Complex>>
                                (state[0].cast<shared_ptr<FESpace>>(), Flags(),idnrs,facs);
                            }
                        }
                      else
                        fes = make_shared<PeriodicFESpace>(state[0].cast<shared_ptr<FESpace>>(),
                                                           Flags(),idnrs);
                      fes->Update();
                      fes->FinalizeUpdate();
                      return fes;
                    }))
    ;



  auto disc_class = py::class_<DiscontinuousFESpace, shared_ptr<DiscontinuousFESpace>, FESpace, NGS_Object>(m, "Discontinuous",
	docu_string(R"delimiter(Discontinuous Finite Element Spaces.
FESpace that splits up all dofs that are shared by several (volume or surface) elements. Every element gets a single copy of that dof. Basis functions become element-local.

Parameters:

fespace : ngsolve.comp.FESpace
    finite element space

BND : boolean or None
    separate across surface elements instead of volume elements (for surface FESpaces)
)delimiter"), py::dynamic_attr());
  disc_class
    .def(py::init([disc_class] (shared_ptr<FESpace> & fes, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(kwargs, disc_class);
                    flags.SetFlag("autoupdate", flags.GetDefineFlag("autoupdate") || fes->DoesAutoUpdate());
                    auto dcfes = make_shared<DiscontinuousFESpace>(fes, flags);
                    // MR: Note that dcfes updates the wrapped space!
                    dcfes->Update();
                    dcfes->FinalizeUpdate();
                    connect_auto_update(dcfes.get());
                    return dcfes;
                  }), py::arg("fespace"))
    /*
    .def(py::pickle([](const PeriodicFESpace* per_fes)
                    {
                      py::list idnrs;
                      for (auto idnr : *per_fes->GetUsedIdnrs())
                        idnrs.append(idnr);
                      auto quasiper_fes = dynamic_cast<const QuasiPeriodicFESpace*>(per_fes);
                      if(quasiper_fes)
                        {
                          py::list fac;
                          for(auto factor : *quasiper_fes->GetFactors())
                            fac.append(factor);
                          return py::make_tuple(per_fes->GetBaseSpace(),idnrs,fac);
                        }
                      return py::make_tuple(per_fes->GetBaseSpace(),idnrs);
                    },
                    [] (py::tuple state) -> shared_ptr<PeriodicFESpace>
                    {
                      auto idnrs = make_shared<Array<int>>();
                      for (auto id : state[1].cast<py::list>())
                        idnrs->Append(id.cast<int>());
                      if(py::len(state)==3)
                        {
                          auto facs = make_shared<Array<Complex>>();
                          for (auto fac : state[2].cast<py::list>())
                            facs->Append(fac.cast<Complex>());
                          auto fes = make_shared<QuasiPeriodicFESpace>
                            (state[0].cast<shared_ptr<FESpace>>(), Flags(),idnrs,facs);
                          fes->Update(glh);
                          fes->FinalizeUpdate(glh);
                          return fes;
                        }
                      auto fes = make_shared<PeriodicFESpace>(state[0].cast<shared_ptr<FESpace>>(),
                                                              Flags(),idnrs);
                      fes->Update(glh);
                      fes->FinalizeUpdate(glh);
                      return fes;
                    }))
    */
    ;



  auto hiddenfes_class = py::class_<HiddenFESpace, shared_ptr<HiddenFESpace>, FESpace, NGS_Object>(m, "Hidden",
	docu_string(R"delimiter(Hidden Finite Element Spaces.
FESpace has elements, but no gobally enumerated dofs, i.e. all dofs are hidden.

Parameters:

fespace : ngsolve.comp.FESpace
    finite element space
)delimiter"), py::dynamic_attr());
  hiddenfes_class
    .def(py::init([disc_class] (shared_ptr<FESpace> & fes, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(kwargs, disc_class);          
                    // MR: Note that HiddenFESpace always updates the wrapped space!
                    flags.SetFlag("autoupdate", flags.GetDefineFlag("autoupdate") || fes->DoesAutoUpdate());
                    auto hiddenfes = make_shared<HiddenFESpace>(fes, flags);
                    hiddenfes->Update();
                    hiddenfes->FinalizeUpdate();
                    connect_auto_update(hiddenfes.get());
                    return hiddenfes;
                  }), py::arg("fespace"))
    ;


  
  py::class_<ReorderedFESpace, shared_ptr<ReorderedFESpace>, FESpace>(m, "Reorder",
	docu_string(R"delimiter(Reordered Finite Element Spaces.
...
)delimiter"))
    .def(py::init([] (shared_ptr<FESpace> & fes, bool autoupdate)
                  {
                    Flags flags = fes->GetFlags();
                    flags.SetFlag("autoupdate", autoupdate || fes->DoesAutoUpdate());
                    auto refes = make_shared<ReorderedFESpace>(fes, flags);
                    // MR: Update() always updates wrapped space
                    refes->Update();
                    refes->FinalizeUpdate();
                    connect_auto_update(refes.get());
                    return refes;
                  }), py::arg("fespace"), py::arg("autoupdate")=false)
    .def("GetClusters", &ReorderedFESpace::GetClusters)
    /*
    .def(py::pickle([](const PeriodicFESpace* per_fes)
                    {
                      py::list idnrs;
                      for (auto idnr : *per_fes->GetUsedIdnrs())
                        idnrs.append(idnr);
                      auto quasiper_fes = dynamic_cast<const QuasiPeriodicFESpace*>(per_fes);
                      if(quasiper_fes)
                        {
                          py::list fac;
                          for(auto factor : *quasiper_fes->GetFactors())
                            fac.append(factor);
                          return py::make_tuple(per_fes->GetBaseSpace(),idnrs,fac);
                        }
                      return py::make_tuple(per_fes->GetBaseSpace(),idnrs);
                    },
                    [] (py::tuple state) -> shared_ptr<PeriodicFESpace>
                    {
                      auto idnrs = make_shared<Array<int>>();
                      for (auto id : state[1].cast<py::list>())
                        idnrs->Append(id.cast<int>());
                      if(py::len(state)==3)
                        {
                          auto facs = make_shared<Array<Complex>>();
                          for (auto fac : state[2].cast<py::list>())
                            facs->Append(fac.cast<Complex>());
                          auto fes = make_shared<QuasiPeriodicFESpace>
                            (state[0].cast<shared_ptr<FESpace>>(), Flags(),idnrs,facs);
                          fes->Update(glh);
                          fes->FinalizeUpdate(glh);
                          return fes;
                        }
                      auto fes = make_shared<PeriodicFESpace>(state[0].cast<shared_ptr<FESpace>>(),
                                                              Flags(),idnrs);
                      fes->Update(glh);
                      fes->FinalizeUpdate(glh);
                      return fes;
                    }))
    */
    ;



  
  py::class_<CompressedFESpace, shared_ptr<CompressedFESpace>, FESpace>(m, "Compress",
	docu_string(R"delimiter(Wrapper Finite Element Spaces.
The compressed fespace is a wrapper around a standard fespace which removes
certain dofs (e.g. UNUSED_DOFs).

Parameters:

fespace : ngsolve.comp.FESpace
    finite element space

active_dofs : BitArray or None
    don't use the COUPLING_TYPEs of dofs to compress the FESpace, 
    but use a BitArray directly to compress the FESpace
)delimiter"))
    .def(py::init([] (shared_ptr<FESpace> & fes,
                      py::object active_dofs)
                  {
                    // shared_ptr<CompoundFESpace> compspace = dynamic_pointer_cast<CompoundFESpace> (fes);
                    // if (compspace)
                      // cout << "yes, we can also compress a CompoundFESpace" << endl;
                    // throw py::type_error("cannot apply compression on CompoundFESpace - Use CompressCompound(..)");
                    auto ret = make_shared<CompressedFESpace> (fes);
                    shared_ptr<BitArray> actdofs = nullptr;
                    if (! py::extract<DummyArgument> (active_dofs).check())
                      dynamic_pointer_cast<CompressedFESpace>(ret)->SetActiveDofs(py::extract<shared_ptr<BitArray>>(active_dofs)());
                    ret->Update();
                    ret->FinalizeUpdate();
                    connect_auto_update(ret.get());
                    return ret;                    
                  }), py::arg("fespace"), py::arg("active_dofs")=DummyArgument())
    .def("SetActiveDofs", [](CompressedFESpace & self, shared_ptr<BitArray> active_dofs)
         {
           self.SetActiveDofs(active_dofs);
         },
         py::arg("dofs"))
    .def("GetActiveDofs", [](CompressedFESpace& self)
    {
      return self.GetActiveDofs();
    })
    .def("GetBaseSpace", [](CompressedFESpace & self)
         {
           return self.GetBaseSpace();
         })
    .def(py::pickle([](const CompressedFESpace* compr_fes)
                    {
                      return py::make_tuple(compr_fes->GetBaseSpace(),compr_fes->GetActiveDofs());
                    },
                    [] (py::tuple state) -> shared_ptr<CompressedFESpace>
                    {
                      auto fes = make_shared<CompressedFESpace>(state[0].cast<shared_ptr<FESpace>>());
                      if (state[1].cast<shared_ptr<BitArray>>())
                        fes->SetActiveDofs(state[1].cast<shared_ptr<BitArray>>());
                      fes->Update();
                      fes->FinalizeUpdate();
                      return fes;
                    }))
    ;


   m.def("CompressCompound", [](shared_ptr<FESpace> & fes, py::object active_dofs) -> shared_ptr<FESpace>
            {
              shared_ptr<CompoundFESpace> compspace = dynamic_pointer_cast<CompoundFESpace> (fes);
              if (!compspace)
                throw py::type_error("Not a CompoundFESpace!");
              else
              {
                if (! py::extract<DummyArgument> (active_dofs).check())
                  throw py::type_error("cannot apply compression on CompoundFESpace with active_dofs");
                Array<shared_ptr<FESpace>> spaces(compspace->GetNSpaces());
                for (int i = 0; i < compspace->GetNSpaces(); i++)
                  spaces[i] = make_shared<CompressedFESpace> ((*compspace)[i]);
                auto ret = make_shared<CompoundFESpace>(compspace->GetMeshAccess(),spaces, compspace->GetFlags());
                // do not deactivate SubSpaceUpdate here
                ret->Update();
                ret->FinalizeUpdate();
                if (spaces[0]->DoesAutoUpdate())
                  ret->SetDoSubspaceUpdate(false);
                connect_auto_update(ret.get());
                return ret;
              }
            }, py::arg("fespace"), py::arg("active_dofs")=DummyArgument());

   py::class_<GlobalInterfaceSpace, shared_ptr<GlobalInterfaceSpace>,
              FESpace>
     (m, "GlobalInterfaceSpace")
     .def(py::init([](shared_ptr<MeshAccess> ma,
                      shared_ptr<CoefficientFunction> mapping,
                      optional<Region> definedon,
                      bool periodic, bool periodicu, bool periodicv,
                      int order, bool complex, bool polar, bool autoupdate)
     {
       auto fes = CreateGlobalInterfaceSpace(ma, mapping, definedon,
                                             periodic, periodicu,
                                             periodicv, order,
                                             complex, polar, autoupdate);
       fes->Update();
       fes->FinalizeUpdate();
       connect_auto_update(fes.get());
       return fes;
     }), "mesh"_a, "mapping"_a, "definedon"_a = nullopt,
          "periodic"_a = false, "periodicu"_a = false,
          "periodicv"_a = false, "order"_a = 3, "complex"_a = false,
          "polar"_a = false, "autoupdate"_a = false)
     ;


  /////////////////////////////// GridFunctionCoefficientFunction /////////////

  py::class_<GridFunctionCoefficientFunction, shared_ptr<GridFunctionCoefficientFunction>, CoefficientFunction>
    (m, "GridFunctionCoefficientFunction")
    .def(py::pickle([] (const GridFunctionCoefficientFunction & gfcf)
                    {
                      return py::make_tuple(gfcf.GetGridFunctionPtr(),
                                            gfcf.generated_from_deriv,
                                            gfcf.generated_from_operator);
                    },
                    [] (py::tuple state) -> shared_ptr<GridFunctionCoefficientFunction>
                    {
                      auto gf = state[0].cast<shared_ptr<GridFunction>>();
                      auto fes = gf->GetFESpace();
                      bool generated_from_deriv = state[1].cast<bool>();
                      string generated_from_operator = state[2].cast<string>();
                      if (generated_from_deriv)
                        return make_shared<GridFunctionCoefficientFunction> (gf,
                                                                             fes->GetFluxEvaluator(),
                                                                             fes->GetFluxEvaluator(BND),
                                                                             fes->GetFluxEvaluator(BBND));
                      
                      if (fes->GetAdditionalEvaluators().Used(generated_from_operator))
                        {
                          auto diffop = fes->GetAdditionalEvaluators()[generated_from_operator];
                          shared_ptr<GridFunctionCoefficientFunction> coef;
                          switch(diffop->VB())
                            {
                            case VOL:
                              return make_shared<GridFunctionCoefficientFunction> (gf, diffop);
                            case BND:
                              return make_shared<GridFunctionCoefficientFunction> (gf, nullptr,diffop);
                            case BBND:
                              return make_shared<GridFunctionCoefficientFunction> (gf, nullptr,nullptr,diffop);
                              break;
                            case BBBND:
                              throw Exception ("there are no Operators with BBBND");
                            }
                        }

                      throw Exception("cannot unpickle GridFunctionCoefficientFunction");
                    }))
    .def("Trace",  [](shared_ptr<GridFunctionCoefficientFunction> self)
         { return self->GetTrace(); },
         "take canonical boundary trace.")    
    
    ;
    

  ////////////////////////////////////// GridFunction //////////////////////////
  
  auto gf_class = py::class_<GF,shared_ptr<GF>, CoefficientFunction, NGS_Object>
    (m, "GridFunction",  "a field approximated in some finite element space", py::dynamic_attr());
  gf_class
    .def(py::init([gf_class](shared_ptr<FESpace> fes, string & name,
                                 py::kwargs kwargs)
    {
      auto flags = CreateFlagsFromKwArgs(kwargs, gf_class);
      flags.SetFlag("novisual");
      auto gf = CreateGridFunction(fes, name, flags);
      gf->Update();
      connect_auto_update(gf.get());
      return gf;
    }), py::arg("space"), py::arg("name")="gfu",
         "creates a gridfunction in finite element space")
    .def_static("__flags_doc__", [] ()
                {
                  return py::dict
                    (
                     py::arg("multidim") = "\n"
                     " Multidimensional GridFunction",
                     py::arg("nested") = "bool = False\n"
		             " Generates prolongation matrices for each mesh level and prolongates\n"
		             " the solution onto the finer grid after a refinement.",
                     py::arg("autoupdate") = "\n"
                     " Automatically update on FE space update"
                     );
                })
    .def(py::pickle([] (const GridFunction& gf)
                    {
                      if (ngcore::parallel_pickling && gf.GetMeshAccess()->GetCommunicator().Size() > 1)
                        {
                          ostringstream str; // (ios::binary);
                          gf.Save(str);
                          string s = str.str();
                          VFlatVector<double> vec(str.str().length()/8, (double*)(void*)&s[0]);
                          shared_ptr<BaseVector> v2 = make_shared<VVector<double>> (vec.Size());
                          *v2 = vec;
                          
                          if (gf.GetMeshAccess()->GetCommunicator().Rank() == 0)
                            return
                              py::make_tuple(gf.GetFESpace(),
                                             gf.GetName(),
                                             Flags(gf.GetFlags()).SetFlag("parallel"),
                                             v2);
                        }

                      py::list state;
                      state.append (gf.GetFESpace());
                      state.append (gf.GetName());
                      state.append (gf.GetFlags());
                      for (int i = 0; i < gf.GetMultiDim(); i++)
                        state.append (gf.GetVectorPtr(i));
                      return py::tuple(state);
                      /*
                      return py::make_tuple(gf.GetFESpace(),
                                            gf.GetName(),
                                            gf.GetFlags(),
                                            gf.GetVectorPtr());
                      */
                    },
                    [] (py::tuple state)
                    {
                      auto gf = CreateGridFunction(state[0].cast<shared_ptr<FESpace>>(),
                                                   state[1].cast<string>(),
                                                   state[2].cast<Flags>());
                      gf->Update();
                      if (!state[2].cast<Flags>().GetDefineFlag("parallel"))
                        {
                          for (int i = 0; i < gf->GetMultiDim(); i++)                          
                            gf->GetVector(i) = *py::cast<shared_ptr<BaseVector>>(state[3+i]);
                        }
                      else
                        {
                          auto vec = py::cast<shared_ptr<BaseVector>>(state[3]);
                          // cout << "unpickle cf, vec = " << *vec << endl;
                          string str((char*)(void*)vec->FVDouble().Data(), 8*vec->Size());
                          istringstream in(str);
                          gf->Load(in);
                        }
                      return gf;
                    }
                    ))
    .def("__str__", [] (GF & self) { return ToString(self); } )
    .def_property_readonly("space", [](GF & self) { return self.GetFESpace(); },
                           "the finite element space")
    .def("Update", [](GF& self) { self.Update(); },
         "update vector size to finite element space dimension after mesh refinement")
    .def_property_readonly("autoupdate", [] (shared_ptr<GridFunction> self) {return self->DoesAutoUpdate();})
    
    .def("Save", [](GF& self, string filename, bool parallel)
         {
           ofstream out(filename, ios::binary);
           if (parallel)
             self.Save(out);
           else
             for (auto d : self.GetVector().FVDouble())
               SaveBin(out, d);
         },
         py::arg("filename"), py::arg("parallel")=false, docu_string(R"raw_string(
Saves the gridfunction into a file.

Parameters:

filename : string
  input file name

parallel : bool
  input parallel

)raw_string"))
    .def("Load", [](GF& self, string filename, bool parallel)
         {
           ifstream in(filename, ios::binary);
           if(in.fail())
             throw Exception("File " + filename + " does not exist!");
           if (parallel)
             self.Load(in);
           else
             for (auto & d : self.GetVector().FVDouble())
               LoadBin(in, d);
         },
         py::arg("filename"), py::arg("parallel")=false, docu_string(R"raw_string(       
Loads a gridfunction from a file.

Parameters:

filename : string
  input file name

parallel : bool
  input parallel

)raw_string"))
    .def("Set", 
         [](shared_ptr<GF> self, spCF cf,
            VorB vb, py::object definedon, bool dualdiffop, bool use_simd, int mdcomp, optional<shared_ptr<BitArray>> definedonelements, int bonus_intorder)
         {
           shared_ptr<TPHighOrderFESpace> tpspace = dynamic_pointer_cast<TPHighOrderFESpace>(self->GetFESpace());          
            Region * reg = nullptr;
            if (py::extract<Region&> (definedon).check())
              reg = &py::extract<Region&>(definedon)();
            
            py::gil_scoped_release release;

            if(tpspace)
            {
              Transfer2TPMesh(cf.get(),self.get(),glh);
              return;
            }            
            if (reg)
              SetValues (cf, *self, *reg, NULL, glh, dualdiffop, use_simd, mdcomp, definedonelements, bonus_intorder);
            else
              SetValues (cf, *self, vb, NULL, glh, dualdiffop, use_simd, mdcomp, definedonelements, bonus_intorder);
         },
         py::arg("coefficient"),
         py::arg("VOL_or_BND")=VOL,
         py::arg("definedon")=DummyArgument(),
	 py::arg("dual")=false,
         py::arg("use_simd")=true,
         py::arg("mdcomp")=0, 
         py::arg("definedonelements")=nullopt,
         py::arg("bonus_intorder")=0,
         docu_string(R"raw_string(
Set values

Parameters:

coefficient : ngsolve.fem.CoefficientFunction
  input CF to set

VOL_or_BND : ngsolve.comp.VorB
  input VOL, BND, BBND, ...

definedon : object
  input definedon region

dual : bool
  If set to true dual shapes are used, otherwise local L2-projection is used.
  Default is False.

use_simd : bool
  If set to false does not use SIMD (for debugging).

mdcomp : int
  .

definedonelements : nullopt
  .

bonus_intorder : int
  Increase numerical integration order.

)raw_string"))

    .def("Interpolate", 
         [](shared_ptr<GF> self, spCF cf,
            py::object definedon, int mdcomp)
         {
           Region * reg = nullptr;
           if (py::extract<Region&> (definedon).check())
             reg = &py::extract<Region&>(definedon)();
           
           py::gil_scoped_release release;
           self->Interpolate (*cf, reg, mdcomp, glh);
         },
         py::arg("coefficient"),
         py::arg("definedon")=DummyArgument(),
         py::arg("mdcomp")=0)
    
    .def_property_readonly("name", &GridFunction::GetName, "Name of the Gridfunction")

    .def_property_readonly("components",
                           [](shared_ptr<GF> self)-> py::tuple
                   { 
                     py::tuple vecs(self->GetNComponents());
                     for (int i = 0; i < self->GetNComponents(); i++)
                       vecs[i] = self->GetComponent(i);
                     return vecs;
                   },
                  "list of gridfunctions for compound gridfunction")

    .def_property_readonly("vec",
                           [](shared_ptr<GF> self)
                           { return self->GetVectorPtr(); },
                           "coefficient vector")

    .def_property_readonly("vecs", 
                           [](shared_ptr<GF> self)-> py::list
                   { 
                     py::list vecs(self->GetMultiDim());
                     for (int i = 0; i < self->GetMultiDim(); i++)
                       vecs[i] = py::cast(self->GetVectorPtr(i));
                     return vecs;
                   },
                  "list of coefficient vectors for multi-dim gridfunction")

    .def("MDComponent", [] (shared_ptr<GF> self, int mdcomp)
         {
           return make_shared<GridFunctionCoefficientFunction> (self, mdcomp); 
         }, py::arg("mdcomp"), "select component of multidim GridFunction")
    
    .def("Deriv",
         [](shared_ptr<GF> self) -> spCF
          {
            return self->Deriv();
          }, "Returns the canonical derivative of the space behind the GridFunction if possible.")

    .def("Trace",  [](shared_ptr<GF> self)
         { return self->GetTrace(); },
         "take canonical boundary trace. This function is optional, added for consistency with proxies")

    .def("Operators", [] (shared_ptr<GF> self)
         {
           py::list l;
           auto ops = self->GetFESpace()->GetAdditionalEvaluators();
           for (size_t i = 0; i < ops.Size(); i++)
             l.append (ops.GetName(i));
           return l;
         },
         "returns list of available differential operators")
    
    .def("Operator",
         [](shared_ptr<GF> self, string name, VorB vb)
          {
            if (!self->GetFESpace()->GetAdditionalEvaluators().Used(name))
              throw Exception(string("Operator \"") + name + string("\" does not exist for ") + self->GetFESpace()->GetClassName() + string("!"));
            auto diffop = self->GetFESpace()->GetAdditionalEvaluators()[name];

            if (!diffop->SupportsVB(vb))
              throw Exception(string("Operator \"") + name + string("\" does not support vb = ") + ToString(vb) + string("!"));

            shared_ptr<GridFunctionCoefficientFunction> coef;
            switch(vb)
              {
              case VOL:
                coef = make_shared<GridFunctionCoefficientFunction> (self, diffop);
                break;
              case BND:
                coef = make_shared<GridFunctionCoefficientFunction> (self, nullptr,diffop);
                break;
              case BBND:
                coef = make_shared<GridFunctionCoefficientFunction> (self, nullptr,nullptr,diffop);
                break;
              case BBBND:
                throw Exception ("there are no Operators with BBBND");
              }
            coef->SetDimensions(diffop->Dimensions());
            coef->generated_from_operator = name;
            return coef;
          }, py::arg("name"), py::arg("VOL_or_BND")=VOL, docu_string(R"raw_string(
Get access to an operator depending on the FESpace.

Parameters:

name : string
  input name of the requested operator

VOL_or_BND : ngsolve.comp.VorB
  input VOL, BND, BBND, ...

)raw_string"))

    
    .def_property_readonly("derivname", 
                           [](shared_ptr<GF> self) -> string
                   {
                     /*
                     auto deriv = self->GetFESpace()->GetFluxEvaluator();
                     if (!deriv) return "";
                     return deriv->Name();
                     */
                     for (auto vb : { VOL, BND, BBND })
                       if (auto deriv = self->GetFESpace()->GetFluxEvaluator(vb))
                         return deriv->Name();
                     return "";
                   }, "Name of canonical derivative of the space behind the GridFunction.")

    .def("__call__", 
         [](shared_ptr<GF> self, double x, double y, double z, VorB vb)
          {
            HeapReset hr(glh);
            auto space = self->GetFESpace();
            auto evaluator = space->GetEvaluator();
            IntegrationPoint ip;
            int elnr = -1;
            if (vb == VOL)
              elnr = space->GetMeshAccess()->FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            else
              elnr = space->GetMeshAccess()->FindSurfaceElementOfPoint(Vec<3>(x, y, z), ip, true);
            if (elnr < 0) throw Exception ("point out of domain");
            ElementId ei(vb, elnr);
            
            const FiniteElement & fel = space->GetFE(ei, glh);

            Array<int> dnums(fel.GetNDof(), glh);
            space->GetDofNrs(ei, dnums);
            auto & trafo = space->GetMeshAccess()->GetTrafo(ei, glh);

            if (space->IsComplex())
              {
                Vector<Complex> elvec(fel.GetNDof()*space->GetDimension());
                Vector<Complex> values(evaluator->Dim());
                self->GetElementVector(dnums, elvec);

                evaluator->Apply(fel, trafo(ip, glh), elvec, values, glh);
                return (values.Size() > 1) ? py::cast(values) : py::cast(values(0));
              }
            else
              {
                Vector<> elvec(fel.GetNDof()*space->GetDimension());
                Vector<> values(evaluator->Dim());
                self->GetElementVector(dnums, elvec);

                evaluator->Apply(fel, trafo(ip, glh), elvec, values, glh);
                return (values.Size() > 1) ? py::cast(values) : py::cast(values(0));
              }
          },
         py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0, py::arg("VOL_or_BND") = VOL)

    // expose CF __call__ to GF, because pybind11 doesn't do that if a function gets overloaded
    .def("__call__", [](shared_ptr<GF> self, py::args args, py::kwargs kwargs)
         {
           return py::module::import("ngsolve").attr("CoefficientFunction").attr("__call__")(self, *args, **kwargs);
         })
    
    .def("CF", 
         [](shared_ptr<GF> self, shared_ptr<DifferentialOperator> diffop) -> spCF
          {
            if (!diffop->Boundary())
              return make_shared<GridFunctionCoefficientFunction> (self, diffop);
            else
              return make_shared<GridFunctionCoefficientFunction> (self, nullptr, diffop);
          }, py::arg("diffop"), docu_string(R"raw_string(

Parameters:

diffop : ngsolve.fem.DifferentialOperator
  input differential operator

)raw_string"))

    .def("AddMultiDimComponent", &GridFunction::AddMultiDimComponent)
    ;

  py::class_<S_GridFunction<double>, shared_ptr<S_GridFunction<double>>, GridFunction>
    (m, "GridFunctionD")
    .def(py::pickle([](const S_GridFunction<double> gf)
                    {
                      return py::make_tuple(gf.GetFESpace(),
                                            gf.GetName(),
                                            gf.GetFlags(),
                                            gf.GetVectorPtr());
                    },
                    [](py::tuple state)
                    {
                      auto gf = CreateGridFunction(state[0].cast<shared_ptr<FESpace>>(),
                                                   state[1].cast<string>(),
                                                   state[2].cast<Flags>());
                      gf->Update();
                      gf->GetVector() = *py::cast<shared_ptr<BaseVector>>(state[3]);
                      return dynamic_pointer_cast<S_GridFunction<double>>(gf);
                    }))
    ;
  py::class_<S_GridFunction<Complex>, shared_ptr<S_GridFunction<Complex>>, GridFunction>
    (m, "GridFunctionC")
    .def(py::pickle([](const S_GridFunction<Complex> gf)
                    {
                      return py::make_tuple(gf.GetFESpace(),
                                            gf.GetName(),
                                            gf.GetFlags(),
                                            gf.GetVectorPtr());
                    },
                    [](py::tuple state)
                    {
                      auto gf = CreateGridFunction(state[0].cast<shared_ptr<FESpace>>(),
                                                   state[1].cast<string>(),
                                                   state[2].cast<Flags>());
                      gf->Update();
                      gf->GetVector() = *py::cast<shared_ptr<BaseVector>>(state[3]);
                      return dynamic_pointer_cast<S_GridFunction<Complex>>(gf);
                    }))
    ;


  py::class_<ComponentGridFunction, shared_ptr<ComponentGridFunction>, GridFunction>
    (m, "ComponentGridFunction")
    .def(py::pickle(
                    [](ComponentGridFunction& cgf)
                    {
                      return py::make_tuple(cgf.GetParent(), cgf.GetComponent());
                    },
                    [](py::tuple state)
                    {
                      auto self = make_shared<ComponentGridFunction>(py::cast<shared_ptr<GridFunction>>(state[0]),
                                                                     py::cast<int>(state[1]));
                      self->Update();
                      return self;
                    }))
    ;

  ///////////////////////////// BilinearForm   ////////////////////////////////////////

  py::class_<DifferentialSymbol>(m, "DifferentialSymbol")
    .def(py::init<VorB>())
    .def("__call__", [](DifferentialSymbol & self,
                        optional<variant<Region,string>> definedon,
                        bool element_boundary,
                        VorB element_vb, bool skeleton,
                        int bonus_intorder,
                        std::map<ELEMENT_TYPE,IntegrationRule> intrules,
                        shared_ptr<GridFunction> deformation,
                        shared_ptr<BitArray> definedonelements)
         {
           if (element_boundary) element_vb = BND;
           auto dx = DifferentialSymbol(self.vb, element_vb, skeleton, /* defon, */ bonus_intorder);
           if (definedon)
             {
               if (auto definedon_region = get_if<Region>(&*definedon); definedon_region)
                 {
                   dx.definedon = definedon_region->Mask();
                   dx.vb = VorB(*definedon_region);
                 }
               if (auto definedon_string = get_if<string>(&*definedon); definedon_string)
                 dx.definedon = *definedon_string;
             }
           dx.deformation = deformation;
           dx.definedonelements = definedonelements;
           for (auto both : intrules)
             dx.userdefined_intrules[both.first] =
               make_shared<IntegrationRule> (both.second.Copy());
           return dx;
         },
         py::arg("definedon")=nullptr,
         py::arg("element_boundary")=false,
         py::arg("element_vb")=VOL,
         py::arg("skeleton")=false,
         py::arg("bonus_intorder")=0,
         py::arg("intrules")=std::map<ELEMENT_TYPE,IntegrationRule>{},
         py::arg("deformation")=nullptr,
         py::arg("definedonelements")=nullptr)
    ;


  py::class_<Integral, shared_ptr<Integral>> (m, "Integral")
    .def_property_readonly("coef", [] (shared_ptr<Integral> igl) { return igl->cf; })
    .def_property_readonly("symbol", [] (shared_ptr<Integral> igl) { return igl->dx; })
    .def("MakeBFI", [](const Integral & igl) { return igl.MakeBilinearFormIntegrator(); })
    .def("MakeLFI", [](const Integral & igl) { return igl.MakeLinearFormIntegrator(); })
    .def("__radd__", [](shared_ptr<Integral> igl, int i) {
        if (i != 0) throw Exception("can only add integer 0 to Integral (for Python sum(list))");
        return make_shared<SumOfIntegrals>(igl); })
    ;
     
  py::class_<SumOfIntegrals, shared_ptr<SumOfIntegrals>>(m, "SumOfIntegrals")
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(float() * py::self)
    .def_property("linearization",
                  [](const SumOfIntegrals& ints)
                  {
                    if(ints.icfs.Size() > 1)
                      throw Exception("linearization property only availaible for single integrals!");
                  },
                  [](SumOfIntegrals& ints, shared_ptr<SumOfIntegrals> linearization)
                  {
                    if(ints.icfs.Size() > 1)
                      throw Exception("linearization property only availaible for single integrals!");
                    if(linearization->icfs.Size() > 1)
                      throw Exception("linearization property only availaible for single integrals!");
                    ints.icfs[0]->linearization = linearization->icfs[0];
                  })
    .def("__len__", [](shared_ptr<SumOfIntegrals> igls)
         { return igls->icfs.Size(); })
    .def ("__getitem__", [](shared_ptr<SumOfIntegrals> igls, int nr)
          {
            if (nr < 0) nr += igls->icfs.Size();
            if (nr < 0 || nr >= igls->icfs.Size())
              throw py::index_error();
            return igls->icfs[nr];
          })
    .def ("Diff", &SumOfIntegrals::Diff)
    .def ("DiffShape", &SumOfIntegrals::DiffShape)
    .def ("Derive", &SumOfIntegrals::Diff, "depricated: use 'Diff' instead")
    .def ("Compile", &SumOfIntegrals::Compile, py::arg("realcompile")=false, py::arg("wait")=false)
    .def("__str__",  [](shared_ptr<SumOfIntegrals> igls) { return ToString(*igls); } )
    .def("__radd__", [](shared_ptr<SumOfIntegrals> igls, int i) {
        if (i != 0) throw Exception("can only add integer 0 to SumOfIntegrals (for Python sum(list))");
        return igls; })
    .def("SetDefinedOnElements", &SumOfIntegrals::SetDefinedOnElements)
    ;

  py::class_<Variation> (m, "Variation")
    .def(py::init<shared_ptr<SumOfIntegrals>>())
    .def ("Compile", &Variation::Compile, py::arg("realcompile")=false, py::arg("wait")=false)
    ;
  
  
  typedef BilinearForm BF;
  auto bf_class = py::class_<BF, shared_ptr<BilinearForm>, NGS_Object>(m, "BilinearForm",
                                             docu_string(R"raw_string(
Used to store the left hand side of a PDE. integrators (ngsolve.BFI)
to it to implement your PDE. If the left hand side is linear
you can use BilinearForm.Assemble to assemble it after adding
your integrators. For nonlinear usage use BilinearForm.Apply or
BilinearForm.AssembleLinearization instead of Bilinearform.Assemble.

Parameters:

space : ngsolve.FESpace
  The finite element space the bilinearform is defined on. This
  can be a compound FESpace for a mixed formulation.

)raw_string"));
  bf_class
    .def(py::init([bf_class] (shared_ptr<FESpace> fespace, const string& name, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(kwargs, bf_class);
                    auto biform = CreateBilinearForm (fespace, name, flags);
                    return biform;
                  }),
      py::arg("space"), "name"_a = "biform_from_py")
    .def(py::init([bf_class](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space, const string& name,
                             py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(kwargs, bf_class);
                    auto biform = CreateBilinearForm (trial_space, test_space, name, flags);
                    return biform;
                  }),
         py::arg("trialspace"),
         py::arg("testspace"), "name"_a = "biform_from_py")

    .def(py::init([bf_class](shared_ptr<SumOfIntegrals> igls, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(kwargs, bf_class);
                    shared_ptr<FESpace> trial_space, test_space;
                    bool found_trial=false, found_test=false;
                    for (auto igl : *igls)
                      igl->cf -> TraverseTree ([&] (CoefficientFunction& cf) {
                          if (auto * proxy = dynamic_cast<ProxyFunction*>(&cf))
                            {
                              if (proxy->IsTrialFunction())
                                {
                                  found_trial = true;
                                  trial_space = proxy->GetFESpace();
                                }
                              else
                                {
                                  found_test = true;
                                  test_space = proxy->GetFESpace();
                                }
                            }
                        });
                    if ( !found_trial || !found_test)
                      throw Exception("BilinearForm must have Trial- and TestFunction");
                    auto biform = (trial_space == test_space) ?
                      CreateBilinearForm (trial_space, "biform_from_py", flags)
                      :
                      CreateBilinearForm (trial_space, test_space, "biform_from_py", flags);
                    py::cast(biform) += py::cast(igls);
                    return biform;
                  }))
    
    .def_static("__flags_doc__", [] ()
                {
                  return py::dict
                    (
                     py::arg("condense") = "bool = False\n"
                     "  (formerly known as 'eliminate_internal')\n"
                     "  Set up BilinearForm for static condensation of internal\n"
                     "  bubbles. Static condensation has to be done by user,\n"
                     "  this enables only the use of the members harmonic_extension,\n"
                     "  harmonic_extension_trans and inner_solve. Have a look at the\n"
                     "  documentation for further information.",
                     py::arg("eliminate_internal") = "bool = False\n"
                     "  deprecated for static condensation, replaced by 'condense'\n",
                     py::arg("eliminate_hidden") = "bool = False\n"
                     "  Set up BilinearForm for static condensation of hidden\n"
                     "  dofs. May be overruled by eliminate_internal.",
                     py::arg("print") = "bool = False\n"
                     "  Write additional information to testout file. \n"
                     "  This file must be set by ngsolve.SetTestoutFile. Use \n"
                     "  ngsolve.SetNumThreads(1) for serial output",
                     py::arg("printelmat") = "bool = False\n"
                     "  Write element matrices to testout file",
                     py::arg("symmetric") = "bool = False\n"
                     "  BilinearForm is symmetric.\n"
                     "  does not imply symmetric_storage, as used to be earlier\n",
                     py::arg("symmetric_storage") = "bool = False\n"
                     "  Store only lower triangular part of sparse matrix.",
                     py::arg("nonassemble") = "bool = False\n"
                     "  BilinearForm will not allocate memory for assembling.\n"
                     "  optimization feature for (nonlinear) problems where the\n"
                     "  form is only applied but never assembled.",
                     py::arg("project") = "bool = False\n"
                     "  When calling bf.Assemble, all saved coarse matrices from\n"
                     "  mesh refinements are updated as well using a Galerkin projection\n"
                     "  of the matrix on the finest grid. This is needed to use the multigrid\n"
                     "  preconditioner with a changing bilinearform.",
		     py::arg("nonsym_storage") = "bool = False\n"
		     "  (deprecated) The full matrix is stored, even if the symmetric flag is set.",
                     py::arg("diagonal") = "bool = False\n"
                     "  Stores only the diagonal of the matrix.",
                     py::arg("geom_free") = "bool = False\n"
                     "  when element matrices are independent of geometry, we store them \n"
                     "  only for the referecne elements",
                     py::arg("check_unused") = "bool = True\n"
		     "  If set prints warnings if not UNUSED_DOFS are not used."
                     );
                })

    .def("__str__",  []( BF & self ) { return ToString<BilinearForm>(self); } )

    .def("Add", [](BF& self, shared_ptr<BilinearFormIntegrator> bfi) -> BF&
                                 { self.AddIntegrator (bfi); return self; },
         py::return_value_policy::reference, py::arg("integrator"), docu_string(R"raw_string(
         Add integrator to bilinear form.

Parameters:

integrator : ngsolve.fem.BFI
  input bilinear form integrator

)raw_string"))
    
    .def("Add", [](py::object self, shared_ptr<SumOfIntegrals> sum) 
         {
           self += py::cast(sum);
           return self;
         })

         .def("__iadd__",[](BF& self, shared_ptr<BilinearFormIntegrator> other) -> BilinearForm& { self += other; return self; }, py::arg("other") )
    .def("__iadd__", [](BF & self, shared_ptr<SumOfIntegrals> sum) -> BilinearForm& 
         {
           for (auto icf : sum->icfs)
             {
               /*
               auto & dx = icf->dx;

               // check for DG terms
               bool has_other = false;
               icf->cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                                      {
                                        if (dynamic_cast<ProxyFunction*> (&cf))
                                          if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                                            has_other = true;
                                      });
               if (has_other && (dx.element_vb != BND) && !dx.skeleton)
                 throw Exception("DG-facet terms need either skeleton=True or element_boundary=True");

               shared_ptr<BilinearFormIntegrator> bfi;
               if (!has_other && !dx.skeleton)
                 bfi = make_shared<SymbolicBilinearFormIntegrator> (icf->cf, dx.vb, dx.element_vb);
               else
                 bfi = make_shared<SymbolicFacetBilinearFormIntegrator> (icf->cf, dx.vb, !dx.skeleton);
               if (dx.definedon)
                 {
                   if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon); definedon_bitarray)
                     bfi->SetDefinedOn(*definedon_bitarray);
                   if (auto definedon_string = get_if<string> (&*dx.definedon); definedon_string)
                     {
                       Region reg(self.GetFESpace()->GetMeshAccess(), dx.vb, *definedon_string);
                       bfi->SetDefinedOn(reg.Mask());
                     }
                 }
               bfi->SetDeformation(dx.deformation);               
               bfi->SetBonusIntegrationOrder(dx.bonus_intorder);
               if(dx.definedonelements)
                 bfi->SetDefinedOnElements(dx.definedonelements);
               for (auto both : dx.userdefined_intrules)
                 bfi->SetIntegrationRule(both.first, *both.second);
               */

               
               shared_ptr<BilinearFormIntegrator> bfi = icf->MakeBilinearFormIntegrator();
               auto & dx = icf->dx;
               if (dx.definedon)
                 {
                   if (auto definedon_string = get_if<string> (&*dx.definedon); definedon_string)
                     {
                       Region reg(self.GetFESpace()->GetMeshAccess(), dx.vb, *definedon_string);
                       bfi->SetDefinedOn(reg.Mask());
                     }
                 }
               
               self += bfi;
             }
           return self;
         })
    
    .def("__iadd__", [](BF & self, Variation variation) -> BilinearForm& 
         {
           for (auto icf : variation.igls->icfs)
             {
               auto & dx = icf->dx;

               auto bfi = make_shared<SymbolicEnergy> (icf->cf, dx.vb, dx.element_vb);
               if (dx.definedon)
                 {
                   if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon); definedon_bitarray)
                     bfi->SetDefinedOn(*definedon_bitarray);
                   if (auto definedon_string = get_if<string> (&*dx.definedon); definedon_string)
                     {
                       Region reg(self.GetFESpace()->GetMeshAccess(), dx.vb, *definedon_string);
                       bfi->SetDefinedOn(reg.Mask());
                     }
                 }
               bfi->SetDeformation(dx.deformation);               
               bfi->SetBonusIntegrationOrder(dx.bonus_intorder);
               for (auto both : dx.userdefined_intrules)
                 bfi->SetIntegrationRule(both.first, *both.second);
               self += bfi;
             }
           return self;
         })
         
    .def_property_readonly("space", [](BF& self) { return self.GetFESpace(); }, "fespace on which the bilinear form is defined on")

    .def_property_readonly("integrators", [](BF & self)
                           { return MakePyTuple (self.Integrators()); }, "integrators of the bilinear form")

    .def_property_readonly("loform", [](shared_ptr<BilinearForm> self) 
			   { return self->GetLowOrderBilinearForm(); })
    
    .def("Assemble", [](shared_ptr<BilinearForm> self, bool reallocate)
         {
           // self->ReAssemble(glh,reallocate);
           self->ReAssemble(lhp.GetLH(),reallocate);           
           return self;
         }, py::call_guard<py::gil_scoped_release>(),
         py::arg("reallocate")=false, docu_string(R"raw_string(
Assemble the bilinear form.

Parameters:

reallocate : bool
  input reallocate

)raw_string"))

    .def_property_readonly("mat", [](shared_ptr<BF> self) -> shared_ptr<BaseMatrix>
                                         {
                                           if (self->NonAssemble())
                                             {
                                               auto app = make_shared<BilinearFormApplication> (self, glh);
                                               if (self->GetTrialSpace()->IsParallel())
                                                 return make_shared<ParallelMatrix>(app,
                                                                                    self->GetTrialSpace()->GetParallelDofs(),
                                                                                    self->GetTestSpace()->GetParallelDofs(), C2D);
                                               return app;
                                             }
                                           auto mat = self->GetMatrixPtr();
                                           if (!mat)
                                             throw py::type_error("matrix not ready - assemble bilinearform first");
                                           return mat;
                                         }, "matrix of the assembled bilinear form")

    .def_property_readonly("components", [](shared_ptr<BilinearForm> self)-> py::list
                   { 
                     py::list bfs;
                     auto fes = dynamic_pointer_cast<CompoundFESpace> (self->GetFESpace());
                     if (!fes)
                       throw py::type_error("not a compound-fespace\n");
                       
                     int ncomp = fes->GetNSpaces();
                     for (int i = 0; i < ncomp; i++)
                       bfs.append(shared_ptr<BilinearForm>(make_shared<ComponentBilinearForm>(self, i, ncomp)));
                     return bfs;
                   },
                  "list of components for bilinearforms on compound-space")
    
    .def_property_readonly("condense", [](shared_ptr<BilinearForm> self)
                           { return self->UsesEliminateInternal(); },
                           "use static condensation ?"
                           )
                           
    .def("__call__", [](BF & self, const GridFunction & u, const GridFunction & v)
          {
            auto au = self.GetMatrix().CreateColVector();
            au = self.GetMatrix() * u.GetVector();
            return InnerProduct (au, v.GetVector());
          }, py::arg("gfu"), py::arg("gfv"))

    .def("Energy",[](BF & self, shared_ptr<BaseVector> x)
          {
            return self.Energy(*x, glh);
          }, py::call_guard<py::gil_scoped_release>(), py::arg("x"), docu_string(R"raw_string(
Computes the energy of EnergyIntegrators like SymbolicEnergy for given input vector.

Parameters:

x : ngsolve.la.BaseVector
  input vector

)raw_string"))
    
    .def("Apply", [](BF & self, BaseVector& x, BaseVector & y)
	  {
	    self.ApplyMatrix (x, y, glh);
	  }, py::call_guard<py::gil_scoped_release>(),
         py::arg("x"),py::arg("y"), docu_string(R"raw_string(
Applies a (non-)linear variational formulation to x and stores the result in y.

Parameters:

x : ngsolve.BaseVector
  input vector

y : ngsolve.BaseVector
  output vector

)raw_string"))

    .def("ComputeInternal", [](BF & self, BaseVector & u, BaseVector & f)
	  {
	    self.ComputeInternal (u, f, glh );
	  }, py::call_guard<py::gil_scoped_release>(),
         py::arg("u"),py::arg("f"), docu_string(R"raw_string(

Parameters:

u : ngsolve.la.BaseVector
  input vector

f : ngsolve.la.BaseVector
  input right hand side

)raw_string"))
    .def("DeleteSpecialElements", &BF::DeleteSpecialElements)

    .def("DeleteMatrix", &BF::DeleteMatrix)
    .def("AssembleLinearization", [](BF & self, BaseVector & ulin,
                                     bool reallocate)
	  {
	    self.AssembleLinearization (ulin, glh, reallocate);
	  }, py::call_guard<py::gil_scoped_release>(),
         py::arg("ulin"), py::arg("reallocate")=false,
         docu_string(R"raw_string(
Computes linearization of the bilinear form at given vecor.

Parameters:

ulin : ngsolve.la.BaseVector
  input vector

)raw_string"))

    .def("Flux", [](BF & self, shared_ptr<GridFunction> gf) -> spCF
          {
            return make_shared<GridFunctionCoefficientFunction> (gf, self.GetIntegrator(0));
          }, py::arg("gf"), docu_string(R"raw_string(

Parameters:

gf : ngsolve.comp.GridFunction
  input GridFunction

)raw_string"))
    
    .def_property_readonly("harmonic_extension", [](BF & self)
                   {
                     return self.GetHarmonicExtension();
                   }, "harmonic_extension used for static condensaition"
                  )
    .def_property_readonly("harmonic_extension_trans", [](BF & self)
                   {
                     return self.GetHarmonicExtensionTrans();
                   }, "harmonic_extension_trans used for static condensation"
                  )
    .def_property_readonly("inner_solve", [](BF & self)
                   {
                     return self.GetInnerSolve();
                   }, "inner_solve used for static condensation"
                  )
    .def_property_readonly("inner_matrix", [](BF & self)
                   {
                     return self.GetInnerMatrix();
                   }, "inner_matrix of the bilinear form"
                  )
    .def("SetPreconditioner", &BF::SetPreconditioner)
    .def("UnsetPreconditioner", &BF::UnsetPreconditioner)
    ;

  ///////////////////////////////// LinearForm //////////////////////////////////////////

  typedef LinearForm LF;
  auto lf_class = py::class_<LF, shared_ptr<LF>, NGS_Object>(m, "LinearForm", docu_string(R"raw_string(
Used to store the left hand side of a PDE. Add integrators
(ngsolve.LFI) to it to implement your PDE.

Parameters:

space : ngsolve.FESpace
  The space the linearform is defined on. Can be a compound
  FESpace for a mixed formulation.

flags : dict
  Additional options for the linearform, for example:

    print : bool
      Write additional debug information to testout file. This
      file must be set by ngsolve.SetTestoutFile. Use
      ngsolve.SetNumThreads(1) for serial output.

)raw_string"));
  lf_class
    .def(py::init([lf_class] (shared_ptr<FESpace> fespace, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(kwargs, lf_class);
                    auto f = CreateLinearForm (fespace, "lff_from_py", flags);
                    f->AllocateVector();
                    return f;
                  }),
         py::arg("space"))

    .def(py::init([lf_class](shared_ptr<SumOfIntegrals> igls, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(kwargs, lf_class);
                    shared_ptr<FESpace> test_space;
                    bool found = false;
                    for (auto igl : *igls)
                      igl->cf -> TraverseTree ([&] (CoefficientFunction& cf) {
                          if (auto * proxy = dynamic_cast<ProxyFunction*>(&cf))
                            {
                              if (proxy->IsTrialFunction())
                                throw Exception("Linearform must not have TrialFunction");
                              else
                                {
                                  found = true;
                                  test_space = proxy->GetFESpace();
                                }
                            }
                        });
                    if (!found) throw Exception("Linearform must have TestFunction");
                    
                    auto liform = CreateLinearForm (test_space, "liform_from_py", flags);
                    py::cast(liform) += py::cast(igls);
                    liform->AllocateVector();
                    return liform;
                  }))


    
    .def_static("__flags_doc__", [] ()
                {
                  return py::dict
                    (
                     py::arg("print") = "bool\n"
                     "  Write additional debug information to testout file.\n"
                     "  This file must be set by ngsolve.SetTestoutFile. Use\n"
                     "  ngsolve.SetNumThreads(1) for serial output.",
                     py::arg("printelvec") = "bool\n"
                     "  print element vectors to testout file"
                     );
                })
    .def("__str__",  [](LF & self ) { return ToString<LinearForm>(self); } )

    .def_property_readonly("vec", [] (shared_ptr<LF> self)
                           { return self->GetVectorPtr();}, "vector of the assembled linear form")
    .def_property_readonly("space", [](LF& self)
                                    { return self.GetFESpace(); })

    .def("Add", [](shared_ptr<LF> self, shared_ptr<LinearFormIntegrator> lfi)
          { 
            self->AddIntegrator (lfi);
            return self; 
          },
         py::arg("integrator"), docu_string(R"raw_string(
Add integrator to linear form.

Parameters:

integrator : ngsolve.fem.LFI
  input linear form integrator

)raw_string"))

    .def("Add", [](py::object self, shared_ptr<SumOfIntegrals> sum) 
         {
           self += py::cast(sum);
           return self;
         })
    
    .def("__iadd__",[](shared_ptr<LF> self, shared_ptr<LinearFormIntegrator> lfi)
         { (*self)+=lfi; return self; }, py::arg("lfi"))

    .def("__iadd__",[](shared_ptr<LF> self, shared_ptr<PointEvaluationFunctional> lfi)
         { (*self)+=lfi; return self; }, py::arg("lfi"))

    .def("__iadd__", [](shared_ptr<LF> self, shared_ptr<SumOfIntegrals> sum) 
         {


           for (auto icf : (*sum))
             {
               shared_ptr<LinearFormIntegrator> lfi = icf->MakeLinearFormIntegrator();
               auto & dx = icf->dx;
               if (dx.definedon)
                 {
                   if (auto definedon_string = get_if<string> (&*dx.definedon); definedon_string)
                     {
                       Region reg(self->GetFESpace()->GetMeshAccess(), dx.vb, *definedon_string);
                       lfi->SetDefinedOn(reg.Mask());
                     }
                 }
               
               *self += lfi;
             }
           return self;

         })

    .def_property_readonly("integrators", [](shared_ptr<LF> self)
                           { return MakePyTuple (self->Integrators()); }, "returns tuple of integrators of the linear form")

    .def("Assemble", [](shared_ptr<LF> self)
         {
           // self->Assemble(glh);
           self->Assemble(lhp.GetLH());
           return self;
         },
         py::call_guard<py::gil_scoped_release>(), "Assemble linear form")
    
    .def_property_readonly("components", [](shared_ptr<LF> self)
                   { 
                     py::list lfs;
                     auto fes = dynamic_pointer_cast<CompoundFESpace> (self->GetFESpace());
                     if (!fes)
                       throw py::type_error("not a compound-fespace\n");
                       
                     int ncomp = fes->GetNSpaces();
                     for (int i = 0; i < ncomp; i++)
                       lfs.append(shared_ptr<LinearForm>(make_shared<ComponentLinearForm>(self, i, ncomp)));
                     return lfs;
                   }, "list of components for linearforms on compound-space")
    
    .def("__call__", [](shared_ptr<LF> self, const GridFunction & v)
          {
            return InnerProduct (self->GetVector(), v.GetVector());
          }, py::arg("gf"))

    ;

  ////////////////////////////// Prolongation ///////////////////////////////

  py::class_<Prolongation, shared_ptr<Prolongation>> (m, "Prolongation")
    .def ("Prolongate", &Prolongation::ProlongateInline, py::arg("finelevel"), py::arg("vec"))
    .def ("Restrict", &Prolongation::RestrictInline, py::arg("finelevel"), py::arg("vec"))
    .def ("CreateMatrix", &Prolongation::CreateProlongationMatrix, py::arg("finelevel"))
    .def ("LevelDofs", &Prolongation::LevelDofs, py::arg("level"))
    .def ("Operator", [](shared_ptr<Prolongation> prol, int level) -> shared_ptr<BaseMatrix>
          {
            return make_shared<ngmg::ProlongationOperator>(prol, level);
          }, py::arg("finelevel"))
    ;
  
  /////////////////////////////// Preconditioner /////////////////////////////////////////////

  auto prec_class = py::class_<Preconditioner, shared_ptr<Preconditioner>, BaseMatrix, NGS_Object>(m, "Preconditioner", py::dynamic_attr());
  prec_class
    .def(py::init([prec_class](shared_ptr<BilinearForm> bfa, const string & type, py::kwargs kwargs)
         {
           auto flags = CreateFlagsFromKwArgs(kwargs, prec_class);

           if (kwargs.contains("blockcreator"))
             {
               auto bc = kwargs["blockcreator"];
               py::print("createor: ", bc);
                 // if it is a compiled C++ user-function we should stay within C++
               // if (auto func = py::cast<function<shared_ptr<Table<DofId>>(const FESpace&)>>(bc))
               if (py::function pyf=bc; pyf.is_cpp_function())
                 {
                   // cout << "it's a C++ function" << endl;
                   auto func = py::cast<function<shared_ptr<Table<DofId>>(const FESpace&)>>(pyf.cpp_function());
                   // cout << "have func object, type(func) = " << typeid(func).name() << endl;
                   // cout << "type cppfunc = " << typeid(pyf.cpp_function()).name() << endl;
                   // typedef shared_ptr<Table<DofId>>(*callbackfunc)(const FESpace &);
                   // cout << "func-ptr = " << func.target<callbackfunc>() << endl;
                   // cout << "function pointer " << (void*)(*func.target<callbackfunc>()) << endl;

                   flags.SetFlag("blockcreator", func);
                 }
               else
                 {
                   cout << "could not extract C++ function" << endl;                   
                   // cout << "create a wrapper" << endl;
                   function<shared_ptr<Table<DofId>>(const FESpace&)> lam =
                     [bc](const FESpace & fes) -> shared_ptr<Table<DofId>>
                     {
                       py::gil_scoped_acquire aq;
                       py::object blocks = bc(py::cast(fes));

                       if (py::isinstance<Table<DofId>>(blocks))
                         return py::cast<shared_ptr<Table<DofId>>>(blocks);

                       // not yet working:
                       // print ("blocks:", blocks);
                       // py::cast<shared_ptr<Table<DofId>>> (blocks);
                       // cout << "cast did work" << endl;
                       // return py::cast<shared_ptr<Table<DofId>>>(blocks);  

                       size_t size = py::len(blocks);
                       Array<int> cnt(size);
                       size_t i = 0;
                       for (auto block : blocks)
                         cnt[i++] = py::len(block);
                       
                       i = 0;
                       auto blocktable = make_shared<Table<DofId>>(cnt);
                       for (auto block : blocks)
                         {
                           auto row = (*blocktable)[i++];
                           size_t j = 0;
                           for (auto val : block)
                             row[j++] = val.cast<int>();
                         }
                       // cout << "blocktable = " << *blocktable << endl;
                       return blocktable;
                     };
                   flags.SetFlag("blockcreator", lam);
                 }
             }
           
           auto creator = GetPreconditionerClasses().GetPreconditioner(type);
           if (creator == nullptr)
             throw Exception(string("nothing known about preconditioner '") + type + "'");
           return creator->creatorbf(bfa, flags, type);
         }),
         py::arg("bf"), py::arg("type"))

    .def_static("__flags_doc__", [] ()
                {
                  return py::dict
                    (
                     py::arg("inverse") = "\n"
                     "  Inverse type used in Preconditioner.",
                     py::arg("test") = "bool = False\n"
                     "  Computes condition number for preconditioner, if testout file\n"
                     "  is set, prints eigenvalues to file."
                     );
                })
    .def ("Test", [](Preconditioner &pre) { pre.Test();}, py::call_guard<py::gil_scoped_release>())
    .def ("Update", [](Preconditioner &pre) { pre.Update();}, py::call_guard<py::gil_scoped_release>(), "Update preconditioner")
    .def_property_readonly("mat", [](Preconditioner &self)
                   {
                     return self.GetMatrixPtr();
                   }, "matrix of the preconditioner")
    ;

  auto prec_multigrid = py::class_<MGPreconditioner, shared_ptr<MGPreconditioner>, Preconditioner>
    (m,"MultiGridPreconditioner");
  prec_multigrid
    .def(py::init([prec_multigrid](shared_ptr<BilinearForm> bfa, const string& name, optional<shared_ptr<Preconditioner>> lo_precond, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(kwargs, prec_multigrid);
                    auto mgpre = make_shared<MGPreconditioner>(bfa,flags, name);
                    if(lo_precond.has_value())
                      mgpre->SetCoarsePreconditioner(lo_precond.value());
                    return mgpre;
                  }), py::arg("bf"), "name"_a = "multigrid", "lo_preconditioner"_a = nullopt)
    .def_static("__flags_doc__", [prec_class] ()
                {
                  auto mg_flags = py::cast<py::dict>(prec_class.attr("__flags_doc__")());
                  mg_flags["updateall"] = "bool = False\n"
                    "  Update all smoothing levels when calling Update";
                  mg_flags["smoother"] = "string = 'point'\n"
                    "  Smoother between multigrid levels, available options are:\n"
                    "    'point': Gauss-Seidel-Smoother\n"
                    "    'line':  Anisotropic smoother\n"
                    "    'block': Block smoother";
                  mg_flags["coarsetype"] = "string = direct\n"
                    "  How to solve coarse problem.";
                  mg_flags["coarsesmoothingsteps"] = "int = 1\n"
                    "  If coarsetype is smoothing, then how many smoothingsteps will be done.";
                  mg_flags["updatealways"] = "bool = False\n";
                  mg_flags["blocktype"] = "str = vertexpatch\n"
                    "  Blocktype used in compound FESpace for smoothing\n"
                    "  blocks. Options: vertexpatch, edgepatch";
                  return mg_flags;
                })

    // not working because shared_ptr<Array<int>> cannot be pybind arg type?
    // TODO CL: Find out why...
    // .def("SetDirectSolverCluster", &MGPreconditioner::SetDirectSolverCluster)
    .def("SetDirectSolverCluster", [](MGPreconditioner& self, py::list pycluster)
    {
      auto cluster = make_shared<Array<int>>(makeCArray<int>(pycluster));
      self.SetDirectSolverCluster(cluster);
    })
    ;


  class PythonPreconditioner : public Preconditioner
  {
  protected:
    shared_ptr<BitArray> freedofs;
    py::object creator;
    shared_ptr<const BaseMatrix>  mat;
    shared_ptr<BaseMatrix> premat;
  public:
    PythonPreconditioner (shared_ptr<BilinearForm> bfa, const Flags & flags,
                          py::object acreator)
      : Preconditioner(bfa, flags), creator(acreator)
    {
      if (bfa->GetMatrixPtr())
        Update();
    }
    void InitLevel (shared_ptr<BitArray> afreedofs) override
    { freedofs = afreedofs; }
    void FinalizeLevel (const ngla::BaseMatrix * amat) override
    {
      // mat = amat;
      // shared_ptr<BaseMatrix> dummy_sp(const_cast<BaseMatrix*>(amat), NOOP_Deleter);
      mat = amat->shared_from_this();
      
      py::gil_scoped_acquire agil;
      // premat = py::cast<shared_ptr<BaseMatrix>> (creator(dummy_sp, freedofs, flags));
      premat = py::cast<shared_ptr<BaseMatrix>> (creator(mat, freedofs, flags));
    }

    void Update() override
    {
      // cout << "update pre" << endl;
      auto bfa = GetBilinearForm();
      freedofs = bfa->GetFESpace()->GetFreeDofs(bfa->UsesEliminateInternal());
      mat = bfa->GetMatrixPtr();

      py::gil_scoped_acquire agil;
      premat = py::cast<shared_ptr<BaseMatrix>> (creator(mat, freedofs, flags));
    }

    const BaseMatrix & GetAMatrix() const override
    {
      return *mat;
    }

    const BaseMatrix& GetMatrix() const override
    {
      return *premat;
    }
    
    shared_ptr<BaseMatrix> GetMatrixPtr() override
    {
      return premat;
    }

    void Mult (const BaseVector & x, BaseVector & y) const override
    {
      premat->Mult(x, y);
    }

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      premat->MultAdd(s, x, y);
    }
  };


  m.def("RegisterPreconditioner", [] (string name, py::object makepre, py::dict docflags)
    {
      // cout << "register Python preconditioner " << name << endl;
      
      auto creator_function = [makepre]
        (shared_ptr<BilinearForm> bfa, const Flags & flags, const string & name)
        -> shared_ptr<Preconditioner> 
        {
          py::gil_scoped_acquire aq;
          return make_shared<PythonPreconditioner> (bfa, flags, makepre);
        };

      DocInfo docinfo;
      for (auto [key,value] : docflags)
        docinfo.Arg(py::cast<string>(key)) = py::cast<string>(value);
      GetPreconditionerClasses().AddPreconditioner (name, nullptr, creator_function);
    }, py::arg("name"), py::arg("makepre"), py::arg("docflags") = py::dict(), "register creator-function makepre(BaseMatrix,FreeDofs)->BaseMatrix");
  
  
  //////////////////////////////////////////////////////////////////////////////////////////

  py::class_<NumProc, NGS_Object, shared_ptr<NumProc>> (m, "NumProc")
    .def("Do", [](NumProc & self)
         {
           self.Do(glh);
         }, py::call_guard<py::gil_scoped_release>())
    ;

  py::class_<PyNumProc, NumProc, shared_ptr<PyNumProc>> (m, "PyNumProc")
    .def(py::init<> ([](shared_ptr<PDE> pde, Flags & flags)
                     { return new PyNumProc(pde, flags); }), py::arg("pde"), py::arg("flags"))
    .def_property_readonly("pde", [](NumProc &self) { return self.GetPDE(); }, "PDE of the NumProc")
    .def("Do", [](NumProc & self, LocalHeap & lh)
                               {
                                 self.Do(lh);
                               }, py::arg("lh"), py::call_guard<py::gil_scoped_release>())
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  PyExportSymbolTable<shared_ptr<FESpace>> (m);
  PyExportSymbolTable<shared_ptr<CoefficientFunction>> (m);
  PyExportSymbolTable<shared_ptr<GridFunction>> (m);
  PyExportSymbolTable<shared_ptr<BilinearForm>>(m);
  PyExportSymbolTable<shared_ptr<LinearForm>>(m);
  PyExportSymbolTable<shared_ptr<Preconditioner>> (m);
  PyExportSymbolTable<shared_ptr<NumProc>> (m);
  PyExportSymbolTable<double> (m);
  PyExportSymbolTable<shared_ptr<double>> (m);

    py::class_<PDE, shared_ptr<PDE>> (m, "PDE")

#ifndef PARALLEL
      .def(py::init([] (const string & filename)
                    {
                      return LoadPDE (filename);
                    }), py::arg("filename"))
#else
      .def(py::init([](const string & filename)
                           { 
                             // ngs_comm = MPI_COMM_WORLD;

                             //cout << "Rank = " << MyMPI_GetId(ngs_comm) << "/"
                             //     << MyMPI_GetNTasks(ngs_comm) << endl;

                             // NGSOStream::SetGlobalActive (MyMPI_GetId(MPI_COMM_WORLD)==0);
                             return LoadPDE (filename);
                           }), py::arg("filename"))
#endif

    .def(py::init<>())

    .def("LoadSolution", []( shared_ptr<PDE> self, string filename, bool ascii )
        {
          return self->LoadSolution(filename, ascii);
        },
         py::arg("filename"), py::arg("ascii")=false
      )

    .def("__str__", [] (shared_ptr<PDE> self) { return ToString(*self); } )

    .def("Mesh",  [](shared_ptr<PDE> self, int meshnr)
        {
          return self->GetMeshAccess(meshnr);
        },
       py::arg("meshnr")=0
       )

    .def("Solve", [](shared_ptr<PDE> self) { self->Solve(); } )


    .def("Add", [](shared_ptr<PDE> self, shared_ptr<MeshAccess> mesh)
                                {
                                  self->AddMeshAccess (mesh);
                                }, py::arg("mesh"))

    .def("Add", [](shared_ptr<PDE> self, const string & name, double val)
                                {
                                  self->AddConstant (name, val);
                                }, py::arg("name"), py::arg("value"))

    .def("Add", [](shared_ptr<PDE> self, shared_ptr<FESpace> space)
                                {
                                  self->AddFESpace (space->GetName(), space);
                                }, py::arg("space"))

    .def("Add", [](shared_ptr<PDE> self, shared_ptr<GridFunction> gf)
                                {
                                  self->AddGridFunction (gf->GetName(), gf);
                                }, py::arg("gf"))

    .def("Add", [](shared_ptr<PDE> self, shared_ptr<BilinearForm> bf)
                                {
                                  self->AddBilinearForm (bf->GetName(), bf);
                                }, py::arg("bf"))

    .def("Add", [](shared_ptr<PDE> self, shared_ptr<LinearForm> lf)
                                {
                                  self->AddLinearForm (lf->GetName(), lf);
                                }, py::arg("lf"))

    .def("Add", [](shared_ptr<PDE> self, shared_ptr<Preconditioner> pre)
                                {
                                  self->AddPreconditioner (pre->GetName(), pre);
                                }, py::arg("pre"))

// TODO
//     .def("Add", [](PyPDE self, shared_ptr<NumProcWrap> np)
//                                 {
//                                   cout << "add pynumproc" << endl;
//                                   self->AddNumProc ("pynumproc", np);
//                                 })
    
    .def("Add", [](shared_ptr<PDE> self, shared_ptr<NumProc> np)
                                {
				  static int cnt = 0;
				  cnt++;
				  string name = "np_from_py" + ToString(cnt);
                                  self->AddNumProc (name, np);
                                }, py::arg("np"))

      .def("Add", [](shared_ptr<PDE> self, shared_ptr<CoefficientFunction> cf, const string& name)
                  {
                    self->AddCoefficientFunction(name, cf);
                  }, py::arg("cf"), py::arg("name"))

      .def("AddPDE_File", [](shared_ptr<PDE> self, const string& other)
                  {
                    LoadPDE(self, other);
                  }, py::arg("filename"), "Adds definitions of other PDE file into existing one")

    .def("Add", [](shared_ptr<PDE> self, const py::list &l)
                                {
                                  for (int i=0; i<py::len(l); i++)
                                    {
                                      py::extract<shared_ptr<PyNumProc>> np(l[i]);
                                      if(np.check())
                                        {
                                          self->AddNumProc (np()->GetName(), np());
                                          continue;
                                        }
                                      
                                      py::extract<shared_ptr<NumProc>> pnp(l[i]);
                                      if(np.check())
                                        {
                                          self->AddNumProc (pnp()->GetName(), pnp());
                                          continue;
                                        }
                                      
                                      py::extract<shared_ptr<GridFunction>> gf(l[i]);
                                      if(gf.check())
                                        {
                                          self->AddGridFunction (gf()->GetName(), gf());
                                          continue;
                                        }
                                      
                                      py::extract<shared_ptr<BilinearForm>> bf(l[i]);
                                      if(gf.check())
                                        {
                                          self->AddBilinearForm (bf()->GetName(), bf());
                                          continue;
                                        }
                                      
                                      py::extract<shared_ptr<LinearForm>> lf(l[i]);
                                      if(gf.check())
                                        {
                                          self->AddLinearForm (lf()->GetName(), lf());
                                          continue;
                                        }
                                      
                                      py::extract<shared_ptr<Preconditioner>> pre(l[i]);
                                      if(gf.check())
                                        {
                                          self->AddPreconditioner (pre()->GetName(), pre());
                                          continue;
                                        }
                                      
                                      cout << "warning: unknown object at position " << i << endl;
                                    }
                                }, py::arg("list"))

    .def("SetCurveIntegrator", [](shared_ptr<PDE> self, const string & filename, shared_ptr<LinearFormIntegrator> lfi)
          {
            self->SetLineIntegratorCurvePointInfo(filename, lfi.get());
          }, py::arg("filename"), py::arg("lfi"))

    .def_property_readonly ("constants", [](shared_ptr<PDE> self) { return py::cast(self->GetConstantTable()); })
    .def_property_readonly ("variables", [](shared_ptr<PDE> self) { return py::cast(self->GetVariableTable()); })
    .def_property_readonly ("coefficients", [](shared_ptr<PDE> self) { return py::cast(self->GetCoefficientTable()); })
    .def_property_readonly ("spaces", [](shared_ptr<PDE> self) {
          auto table = self->GetSpaceTable();
          SymbolTable<shared_ptr<FESpace>> pytable;
          for ( auto i : Range(table.Size() ))
                pytable.Set(table.GetName(i), shared_ptr<FESpace>(table[i]));
          return py::cast(pytable);
          })
    .def_property_readonly ("gridfunctions", [](shared_ptr<PDE> self) {
          auto table = self->GetGridFunctionTable();
          SymbolTable<shared_ptr<GridFunction>> pytable;
          for ( auto i : Range(table.Size() ))
                pytable.Set(table.GetName(i), shared_ptr<GridFunction>(table[i]));
          return py::cast(pytable);
          })
    .def_property_readonly ("bilinearforms", [](shared_ptr<PDE> self) {
          auto table = self->GetBilinearFormTable();
          SymbolTable<shared_ptr<BilinearForm>> pytable;
          for ( auto i : Range(table.Size() ))
                pytable.Set(table.GetName(i), shared_ptr<BilinearForm>(table[i]));
          return py::cast(pytable);
          })
    .def_property_readonly ("linearforms", [](shared_ptr<PDE> self) {
          auto table = self->GetLinearFormTable();
          SymbolTable<shared_ptr<LinearForm>> pytable;
          for ( auto i : Range(table.Size() ))
                pytable.Set(table.GetName(i), shared_ptr<LinearForm>(table[i]));
          return py::cast(pytable);
          })
    .def_property_readonly ("preconditioners", [](shared_ptr<PDE> self) { return py::cast(self->GetPreconditionerTable()); })
    .def_property_readonly ("numprocs", [](shared_ptr<PDE> self) { return py::cast(self->GetNumProcTable()); })
    ;
  
  m.def("Integrate", 
        [](spCF cf,
           variant<shared_ptr<MeshAccess>,Region> mesh_or_reg, 
           VorB vb, int order,
           // std::optional<Region> definedon,
           Region * definedon,
	   bool region_wise, bool element_wise)
        {
          static Timer t("Integrate CF"); RegionTimer reg(t);
          // static mutex addcomplex_mutex;

          shared_ptr<MeshAccess> ma;
          if(auto is_region = get_if<Region>(&mesh_or_reg)) // ; is_region)
            {
              // cout << "is region" << endl;
              definedon = is_region;
              ma = is_region->Mesh();
              vb = is_region->VB();
            }
          if(auto is_mesh = get_if<shared_ptr<MeshAccess>>(&mesh_or_reg)) // ; is_mesh)
            {
              // cout << "is mesh" << endl;
              ma = *is_mesh;
            }
          
          
          BitArray mask;
          if (definedon)
            {
              vb = VorB(*definedon);
              mask = BitArray((*definedon).Mask());
            }
          if(!mask.Size()){
            mask = BitArray(ma->GetNRegions(vb));
            mask.Set();
          }
 
          int dim = cf->Dimension();
          if((region_wise || element_wise) && dim != 1)
            throw Exception("region_wise and element_wise only implemented for 1 dimensional coefficientfunctions");

          cf -> TraverseTree
            ([&] (CoefficientFunction & stepcf)
             {
               if (dynamic_cast<ProxyFunction*>(&stepcf))
                 throw Exception("Cannot integrate ProxFunction!");
             });
                   
          if (!cf->IsComplex())
            {
              Vector<> sum(dim);
              sum = 0.0;
              Vector<> region_sum(region_wise ? ma->GetNRegions(vb) : 0);
              Vector<> element_sum(element_wise ? ma->GetNE(vb) : 0);
              region_sum = 0;
              element_sum = 0;
              bool use_simd = true;
              
              ma->IterateElements
                (vb, glh, [&] (Ngs_Element el, LocalHeap & lh)
                 {
                   if(!mask.Test(el.GetIndex())) return;
                   auto & trafo = ma->GetTrafo (el, lh);
                   FlatVector<> hsum(dim, lh);
                   hsum = 0.0;
                   bool this_simd = use_simd;
                   
                   if (this_simd)
                     {
                       try
                         {
                           SIMD_IntegrationRule ir(trafo.GetElementType(), order);
                           auto & mir = trafo(ir, lh);
                           FlatMatrix<SIMD<double>> values(dim,ir.Size(), lh);
                           cf -> Evaluate (mir, values);
                           FlatVector<SIMD<double>> vsum(dim, lh);
                           vsum = 0;
                           for (size_t j = 0; j < dim; j++)
                             for (size_t i = 0; i < values.Width(); i++)
                               vsum(j) += mir[i].GetWeight() * values(j,i);
                           for(int i = 0; i< dim; i++)
                             hsum[i] = HSum(vsum[i]);
                         }
                       catch (ExceptionNOSIMD e)
                         {
                           this_simd = false;
                           use_simd = false;
                           hsum = 0.0;
                         }
                     }
                   if (!this_simd)
                     {
                       IntegrationRule ir(trafo.GetElementType(), order);
                       BaseMappedIntegrationRule & mir = trafo(ir, lh);
                       FlatMatrix<> values(ir.Size(), dim, lh);
                       cf -> Evaluate (mir, values);
                       for (int i = 0; i < values.Height(); i++)
                         hsum += mir[i].GetWeight() * values.Row(i);
                     }
                   for(size_t i = 0; i<dim;i++)
                     AtomicAdd(sum(i), hsum(i));
                   if(region_wise)
                     AtomicAdd(region_sum(el.GetIndex()), hsum(0));
                   if (element_wise)
                     element_sum(el.Nr()) = hsum(0);
                 });
              py::gil_scoped_acquire aq;
              py::object result;
              if (region_wise) {
#ifdef PARALLEL
                if (ma->GetCommunicator().Size() > 1)
                  MPI_Allreduce(MPI_IN_PLACE, &region_sum(0), ma->GetNRegions(vb), MPI_DOUBLE, MPI_SUM, ma->GetCommunicator());                  
                  /*
                  {
                    Vector<> rs2(ma->GetNRegions(vb));
                    MPI_Allreduce(&region_sum(0), &rs2(0), ma->GetNRegions(vb), MPI_DOUBLE, MPI_SUM, ma->GetCommunicator());
                    region_sum = rs2;
                  }
                  */
#endif
                // result = py::list(py::cast(region_sum));  // crashes ?!?!
                result = py::cast(region_sum);
              }
              else if (element_wise)
                result = py::cast(element_sum);
              else if(dim==1) {
                sum(0) = ma->GetCommunicator().AllReduce(sum(0), MPI_SUM);
                result = py::cast(sum(0));
              }
              else {
#ifdef PARALLEL
                if (ma->GetCommunicator().Size() > 1)
                  MPI_Allreduce(MPI_IN_PLACE, &sum(0), dim, MPI_DOUBLE, MPI_SUM, ma->GetCommunicator());
                /*                  
                  {
                    Vector<> gsum(dim);
                    MPI_Allreduce(&sum(0), &gsum(0), dim, MPI_DOUBLE, MPI_SUM, ma->GetCommunicator());
                    sum = gsum;
                  }
                */
#endif
                result = py::cast(sum);
              }
              return result;
            }
          else
            {
              Vector<Complex> sum(dim);
              sum = 0.0;
              Vector<Complex> region_sum(region_wise ? ma->GetNRegions(vb) : 0);
              Vector<Complex> element_sum(element_wise ? ma->GetNE(vb) : 0);
              region_sum = 0;
              element_sum = 0;
              
              bool use_simd = true;
              
              ma->IterateElements
                (vb, glh, [&] (Ngs_Element el, LocalHeap & lh)
                 {
                   if(!mask.Test(el.GetIndex())) return;
                   auto & trafo = ma->GetTrafo (el, lh);
                   FlatVector<Complex> hsum(dim, lh);
                   hsum = 0.0;
                   
                   bool this_simd = use_simd;
                   
                   if (this_simd)
                     {
                       try
                         {
                           SIMD_IntegrationRule ir(trafo.GetElementType(), order);
                           auto & mir = trafo(ir, lh);
                           FlatMatrix<SIMD<Complex>> values(dim, ir.Size(), lh);
                           cf -> Evaluate (mir, values);
                           FlatVector<SIMD<Complex>> vsum(dim,lh);
                           vsum = Complex(0.0);
                           for (size_t j = 0; j < dim; j++)
                             for (size_t i = 0; i < values.Width(); i++)
                               vsum(j) += mir[i].GetWeight() * values(j,i);
                           for(size_t i = 0; i < dim; i++)
                             hsum[i] = HSum(vsum[i]);
                         }
                       catch (ExceptionNOSIMD e)
                         {
                           this_simd = false;
                           use_simd = false;
                           hsum = 0.0;
                         }
                     }
                   if (!this_simd)
                     {
                       IntegrationRule ir(trafo.GetElementType(), order);
                       BaseMappedIntegrationRule & mir = trafo(ir, lh);
                       FlatMatrix<Complex> values(ir.Size(), dim, lh);
                       cf -> Evaluate (mir, values);
                       for (int i = 0; i < values.Height(); i++)
                         hsum += mir[i].GetWeight() * values.Row(i);
                     }
                   for(size_t i = 0; i<dim; i++)
                     AtomicAdd (sum(i), hsum(i));
                   if(region_wise)
                     AtomicAdd (region_sum(el.GetIndex()), hsum(0));
                   if (element_wise)
                     element_sum(el.Nr()) = hsum(0);
                 });
              
              py::gil_scoped_acquire aq;
              py::object result;
              if (region_wise) {
#ifdef PARALLEL
                if (ma->GetCommunicator().Size() > 1)
                  MPI_Allreduce(MPI_IN_PLACE, &region_sum(0), ma->GetNRegions(vb),
                                MPI_typetrait<Complex>::MPIType(), MPI_SUM, ma->GetCommunicator());
                
                /*
                Vector<Complex> rs2(ma->GetNRegions(vb));
                if (ma->GetCommunicator().Size() > 1)
                  MPI_Allreduce(&region_sum(0), &rs2(0), ma->GetNRegions(vb), MPI_typetrait<Complex>::MPIType(), MPI_SUM, ma->GetCommunicator());
                region_sum = rs2;
                */
#endif
                // result = py::list(py::cast(region_sum));
                result = py::cast(region_sum);
              }
              else if (element_wise)
                result = py::cast(element_sum);
              else if(dim==1) {
                sum(0) = ma->GetCommunicator().AllReduce(sum(0), MPI_SUM);
                result = py::cast(sum(0));
              }
              else {
#ifdef PARALLEL
                if (ma->GetCommunicator().Size() > 1)
                  MPI_Allreduce(MPI_IN_PLACE, &sum(0), dim, MPI_typetrait<Complex>::MPIType(), MPI_SUM, ma->GetCommunicator());
                /*
                Vector<Complex> gsum(dim);
                if (ma->GetCommunicator().Size() > 1) {
                  MPI_Allreduce(&sum(0), &gsum(0), dim, MPI_typetrait<Complex>::MPIType(), MPI_SUM, ma->GetCommunicator());
		  sum = gsum;
		}
                */
#endif
                result = py::cast(sum);
              }
              return result;
            }
        },
	py::arg("cf"), py::arg("mesh"), py::arg("VOL_or_BND")=VOL, 
	py::arg("order")=5,
	py::arg("definedon") = nullptr, // =DummyArgument(),
        py::arg("region_wise")=false,
	py::arg("element_wise")=false,
        R"raw(
Parameters
----------

cf: ngsolve.CoefficientFunction
  Function to be integrated. Can be vector valued, then the result is an array. If you want to integrate
  a lot of functions on the same domain, it will be faster to put them into a vector valued function,
  NGSolve will then be able to use parallelization and SIMD vectorization more efficiently.

mesh: ngsolve.Mesh
  The mesh to be integrated on.

VOL_or_BND: ngsolve.VorB = VOL
  Co-dimension to be integrated on. Historically this could be volume (VOL) or boundary (BND). If your mesh
  contains co-dim 2 elements this can now be BBND (edges in 3d) as well.

order: int = 5
  Integration order, polynomials up to this order will be integrated exactly.

definedon: ngsolve.Region
  Region to be integrated on. Such region can be created with mesh.Boundaries('bcname') or mesh.Materials('matname')
  it will overwrite the VOL_or_BND argument if given.

region_wise: bool = False
  Integrates region wise on the co-dimension given by VOL_or_BND. Returns results as an array, matching the array
  returned by mesh.GetMaterials() or mesh.GetBoundaries(). Does not support vector valued CoefficientFunctions.

element_wise: bool = False
  Integrates element wise and returns result in a list. This is typically used for local error estimators.
  Does not support vector valued CoefficientFunctions
)raw",
        py::call_guard<py::gil_scoped_release>())
    ;


  m.def ("Integrate",
         [] (const SumOfIntegrals & igls, const MeshAccess & ma, bool element_wise) -> py::object
         {
           bool iscomplex = false;
           for (auto & ci : igls)
             {
               iscomplex |= ci->cf->IsComplex();
               if (ci->cf->Dimension() > 1)
                 throw Exception("Integrate(cf*dx) needs scalar-valued CoefficientFunction");
             }

           auto integrate = [&] (auto tscal) 
           {
             typedef decltype(tscal) TSCAL;
             
             TSCAL sum = 0;
             Vector<TSCAL> elvals(element_wise ? ma.GetNE() : 0);
             elvals = TSCAL(0.0);
             
             for (auto & ci : igls)
               sum += ci->Integrate(ma, elvals);
             if (element_wise) return py::cast(elvals);
             return py::cast(sum);
           };

           if (iscomplex)
             return integrate(Complex(0.0));
           else
             return integrate(double(0.0));
           /*
           if (iscomplex)
             {
               Complex sum = 0;
               for (auto & ci : igls.icfs)
                 sum += ci->Integrate<Complex>(ma);
               return py::cast(sum);
             }
           else
             {
               double sum = 0;
               for (auto & ci : igls.icfs)
                 sum += ci->Integrate<double>(ma);
               return py::cast(sum);
             }
           */
         }, py::arg("igls"), py::arg("mesh"), py::arg("element_wise")=false);

  
  m.def("SymbolicLFI",
          [](spCF cf, VorB vb, bool element_boundary,
             bool skeleton, optional<variant<Region, py::list>> definedon,
             IntegrationRule ir, int bonus_intorder, shared_ptr<BitArray> definedonelem,
             bool simd_evaluate, VorB element_vb,
             shared_ptr<GridFunction> deformation) 
           {
             if(definedon.has_value())
               if(auto defregion = get_if<Region>(&*definedon); defregion)
                 vb = VorB(*defregion);

             if (element_boundary) element_vb = BND;

             shared_ptr<LinearFormIntegrator> lfi;
             if (!skeleton)
               lfi = make_shared<SymbolicLinearFormIntegrator> (cf, vb, element_vb);
             else
               lfi = make_shared<SymbolicFacetLinearFormIntegrator> (cf, vb /* , element_boundary */);
             
             if(definedon.has_value())
               {
                 if(auto defpylist = get_if<py::list>(&*definedon); defpylist)
                   {
                     Array<int> defon = makeCArray<int> (*defpylist);
                     for (int & d : defon) d--;
                     lfi -> SetDefinedOn (defon);
                   }
                 if(auto defregion = get_if<Region>(&*definedon); defregion)
                   lfi->SetDefinedOn(defregion->Mask());
               }

             lfi->SetSimdEvaluate (simd_evaluate);
             lfi->SetDeformation (deformation);

             lfi -> SetBonusIntegrationOrder(bonus_intorder);
	     if (ir.Size())
               {
                 cout << IM(1) << "WARNING: Setting the integration rule for all element types is deprecated, use LFI.SetIntegrationRule(ELEMENT_TYPE, IntegrationRule) instead!" << endl;
                 dynamic_pointer_cast<SymbolicLinearFormIntegrator>
		   (lfi)->SetIntegrationRule(ir);                   
               }

             if (definedonelem)
               lfi -> SetDefinedOnElements (definedonelem);

             return shared_ptr<LinearFormIntegrator>(lfi);
           },
           py::arg("form"),
           py::arg("VOL_or_BND")=VOL,
           py::arg("element_boundary")=false,
           py::arg("skeleton")=false,           
           py::arg("definedon")=nullptr,
	   py::arg("intrule")=IntegrationRule(),
           py::arg("bonus_intorder")=0,
           py::arg("definedonelements")=nullptr,
           py::arg("simd_evaluate")=true,
           py::arg("element_vb")=VOL,
           py::arg("deformation")=shared_ptr<GridFunction>(),
        docu_string(R"raw_string(
A symbolic linear form integrator, where test and trial functions, CoefficientFunctions, etc. can be used to formulate right hand sides in a symbolic way.

Parameters:

form : ngsolve.fem.CoefficientFunction
  input the symbolic right hand side form

VOL_or_BND : ngsolve.comp.VorB
  input VOL, BND, BBND, ...

element_boundary : bool
  input element_boundary. True -> iterates over all element boundaries, but uses volume transformations

skeleton : bool
  input skeleton. True -> iterates over all faces, but uses volume transformations

definedon : object
  input definedon region

intrule : ngsolve.fem.IntegrationRule
  input integration rule

bonus_intorder : int
  input additional integration order

definedonelements : object
  input BitArray that marks all elements or facets (for skeleton-integrators) that the integrator is applied on

simd_evaluate : bool
  input simd_evaluate. True -> tries to use SIMD for faster evaluation

element_vb : ngsolve.fem.VorB
  input element VorB

deformation : ngsolve.comp.GridFunction
  input GridFunction to transform/deform the linear form with

)raw_string")
          );

  m.def("SymbolicBFI",
          [](spCF cf, VorB vb, bool element_boundary,
             bool skeleton, optional<variant<Region, py::list>> definedon,
             IntegrationRule ir, int bonus_intorder, shared_ptr<BitArray> definedonelem,
             bool simd_evaluate, VorB element_vb, bool geom_free,
             shared_ptr<GridFunction> deformation)
           {
             if(definedon.has_value())
               if(auto defregion = get_if<Region>(&*definedon); defregion)
                 vb = VorB(*defregion);

             if (element_boundary) element_vb = BND;
             // check for DG terms
             bool has_other = false;

             cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                               {
                                 if (dynamic_cast<ProxyFunction*> (&cf))
                                   if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                                     has_other = true;
                               });
             // if (has_other && !element_boundary && !skeleton)
             if (has_other && (element_vb != BND) && !skeleton)
               throw Exception("DG-facet terms need either skeleton=True or element_boundary=True");
             
             shared_ptr<BilinearFormIntegrator> bfi;
             if (!has_other && !skeleton)
               bfi = make_shared<SymbolicBilinearFormIntegrator> (cf, vb, element_vb);
             else
               bfi = make_shared<SymbolicFacetBilinearFormIntegrator> (cf, vb, element_boundary);
             bfi->geom_free = geom_free;
             if(definedon.has_value())
               {
                 if(auto defpylist = get_if<py::list>(&*definedon); defpylist)
                   {
                     Array<int> defon = makeCArray<int> (*defpylist);
                     for (int & d : defon) d--;
                     bfi -> SetDefinedOn (defon);
                   }
                 if(auto defregion = get_if<Region>(&*definedon); defregion)
                   bfi->SetDefinedOn(defregion->Mask());
               }
             bfi -> SetBonusIntegrationOrder(bonus_intorder);
             if (ir.Size())
               {
                 cout << IM(1) << "WARNING: Setting the integration rule for all element types is deprecated, use BFI.SetIntegrationRule(ELEMENT_TYPE, IntegrationRule) instead!" << endl;
                 /*
                 dynamic_pointer_cast<SymbolicBilinearFormIntegrator> (bfi)
                   ->SetIntegrationRule(ir);
                 */
                 bfi->SetIntegrationRule(ir);
               }

             bfi->SetSimdEvaluate (simd_evaluate);
             bfi->SetDeformation (deformation);
             if (definedonelem)
               bfi -> SetDefinedOnElements (definedonelem);
             return shared_ptr<BilinearFormIntegrator>(bfi);
           },
        py::arg("form"), py::arg("VOL_or_BND")=VOL,
        py::arg("element_boundary")=false,
        py::arg("skeleton")=false,
        py::arg("definedon")=nullptr,
        py::arg("intrule")=IntegrationRule(),
        py::arg("bonus_intorder")=0,
        py::arg("definedonelements")=nullptr,
        py::arg("simd_evaluate")=true,
        py::arg("element_vb")=VOL,
        py::arg("geom_free")=false,        
        py::arg("deformation")=shared_ptr<GridFunction>(),
        docu_string(R"raw_string(
A symbolic bilinear form integrator, where test and trial functions, CoefficientFunctions, etc. can be used to formulate PDEs in a symbolic way.

Parameters:

form : ngsolve.fem.CoefficientFunction
  input the symbolic right hand side form

VOL_or_BND : ngsolve.comp.VorB
  input VOL, BND, BBND, ...

element_boundary : bool
  input element_boundary. True -> iterates over all element boundaries, but uses volume transformations

skeleton : bool
  input skeleton. True -> iterates over all faces, but uses volume transformations

definedon : object
  input definedon region

intrule : ngsolve.fem.IntegrationRule
  input integration rule

bonus_intorder : int
  input additional integration order

definedonelements : object
  input definedonelements

simd_evaluate : bool
  input simd_evaluate. True -> tries to use SIMD for faster evaluation

element_vb : ngsolve.comp.VorB
  input element_vb. Used for skeleton formulation. VOL -> interior faces, BND -> boundary faces

deformation : ngsolve.comp.GridFunction
  input GridFunction to transform/deform the bilinear form with

)raw_string")
        );
          
  m.def("SymbolicTPBFI",
          [](spCF cf, VorB vb, bool element_boundary,
              bool skeleton, py::object definedon)
           {
             py::extract<Region> defon_region(definedon);
             if (defon_region.check())
               vb = VorB(defon_region());

             // check for DG terms
             bool has_other = false;

             cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                               {
                                 if (dynamic_cast<ProxyFunction*> (&cf))
                                   if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                                     has_other = true;
                               });
             if (has_other && !element_boundary && !skeleton)
               throw Exception("DG-facet terms need either skeleton=True or element_boundary=True");
             
             shared_ptr<BilinearFormIntegrator> bfi;
             if (!has_other && !skeleton)
               bfi = make_shared<TensorProductBilinearFormIntegrator> (cf, vb, element_boundary);
             else
               bfi = make_shared<TensorProductFacetBilinearFormIntegrator> (cf, vb, element_boundary);
             
             if (py::extract<py::list> (definedon).check())
               bfi -> SetDefinedOn (makeCArray<int> (definedon));

             if (defon_region.check())
               {
                 cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
                 bfi->SetDefinedOn(defon_region().Mask());
               }
             
             return shared_ptr<BilinearFormIntegrator>(bfi);
           },
           py::arg("form"), py::arg("VOL_or_BND")=VOL,
           py::arg("element_boundary")=false,
           py::arg("skeleton")=false,
           py::arg("definedon")=DummyArgument()
          );
          
  m.def("SymbolicEnergy",
        [](spCF cf, VorB vb, py::object definedon, bool element_boundary,
           int bonus_intorder, py::object definedonelem, bool simd_evaluate,
           VorB element_vb, shared_ptr<GridFunction> deformation)
        -> shared_ptr<BilinearFormIntegrator>
           {
             py::extract<Region> defon_region(definedon);
             if (defon_region.check())
               vb = VorB(defon_region());

             if (element_boundary) element_vb = BND;


             auto bfi = make_shared<SymbolicEnergy> (cf, vb, element_vb);
             bfi -> SetBonusIntegrationOrder(bonus_intorder);
             if (defon_region.check())
               {
                 cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
                 bfi->SetDefinedOn(defon_region().Mask());
               }
             if (! py::extract<DummyArgument> (definedonelem).check())
               bfi -> SetDefinedOnElements (py::extract<shared_ptr<BitArray>>(definedonelem)());
             bfi->SetSimdEvaluate (simd_evaluate);
             bfi->SetDeformation (deformation);
             return bfi;
           },
        py::arg("form"), py::arg("VOL_or_BND")=VOL,
        py::arg("definedon")=DummyArgument(), py::arg("element_boundary")=false,
        py::arg("bonus_intorder")=0,
        py::arg("definedonelements")=DummyArgument(),
        py::arg("simd_evaluate")=true,
        py::arg("element_vb")=VOL,
        py::arg("deformation")=shared_ptr<GridFunction>(),
        docu_string(R"raw_string(
A symbolic energy form integrator, where test and trial functions, CoefficientFunctions, etc. can be used to formulate PDEs in a symbolic way.

Parameters:

form : ngsolve.fem.CoefficientFunction
  input the symbolic right hand side form

VOL_or_BND : ngsolve.comp.VorB
  input VOL, BND, BBND, ...

definedon : object
  input definedon region

element_boundary : bool
  input element_boundary. True -> iterates over all element boundaries, but uses volume transformations

bonus_intorder : int
  input additional integration order

definedonelements : object
  input definedonelements

simd_evaluate : bool
  input simd_evaluate. True -> tries to use SIMD for faster evaluation

element_vb : ngsolve.fem.VorB
  input eleemnt VorB

deformation : ngsolve.comp.GridFunction
  input GridFunction to transform/deform the bilinear form with

)raw_string")
          );

   m.def("KSpaceCoeffs", [](shared_ptr<GF> gf_tp,shared_ptr<GF> gf_k, double x, double y)
           {
			 HeapReset hr(glh);
             auto tpfes = dynamic_pointer_cast<TPHighOrderFESpace>(gf_tp->GetFESpace());
             IntegrationPoint ip;
             int elnr = tpfes->Spaces(0)[0]->GetMeshAccess()->FindElementOfPoint(Vec<2>(x,y),ip,false);
             int ndofx = tpfes->Spaces(0)[0]->GetFE(ElementId(VOL,elnr),glh).GetNDof();
             int ndofy = tpfes->Spaces(0)[1]->GetFE(ElementId(VOL,0),glh).GetNDof();
             Array<int> indices(2);
             Array<int> dofs(ndofx*ndofy);
             tpfes->GetDofNrs(tpfes->GetIndex(elnr,0),dofs);
             FlatMatrix<> elmat(ndofx,ndofy,glh);
             gf_tp->GetVector().GetIndirect(dofs,elmat.AsVector());
             FlatVector<> shape(ndofx,glh);
             dynamic_cast<const BaseScalarFiniteElement& >(tpfes->Spaces(0)[0]->GetFE(ElementId(VOL,elnr),glh)).CalcShape(ip,shape);
             FlatVector<> result(ndofy,glh);
             result = Trans(elmat)*shape;
             gf_k->GetVector().FVDouble() = result;
             
           });
  
   m.def("TensorProductFESpace", [](py::list spaces_list, const Flags & flags ) -> shared_ptr<FESpace>
            {
              auto spaces = makeCArray<shared_ptr<FESpace>> (spaces_list);
              if(spaces.Size() == 2)
              {
                shared_ptr<FESpace> space( new TPHighOrderFESpace( spaces, flags ) );
                return space;
              }
              else
              {
                Array<shared_ptr<FESpace>> spaces_y(spaces.Size()-1);
                for(int i=1;i<spaces.Size();i++)
                  spaces_y[i-1] = spaces[i];
                shared_ptr<FESpace> space( new TPHighOrderFESpace( spaces[0],spaces_y, flags ) );
                return space;             
              }
            },
            py::arg("spaces"),
            py::arg("flags")=Flags()
           );

   m.def("TensorProductIntegrate", [](shared_ptr<GF> gf_tp, py::list ax0, spCF coef) -> double
           {
             static Timer tall("comp.TensorProductIntegrate - single point"); RegionTimer rall(tall);
             Array<double> x0_help = makeCArray<double> (ax0);
             LocalHeap lh(10000000,"TensorProductIntegrate");
             shared_ptr<TPHighOrderFESpace> tpfes = dynamic_pointer_cast<TPHighOrderFESpace>(gf_tp->GetFESpace());
             const Array<shared_ptr<FESpace> > & spaces = tpfes->Spaces(0);
             FlatVector<> x0(spaces[0]->GetSpatialDimension(),&x0_help[0]);
             IntegrationPoint ip;
             int elnr = spaces[0]->GetMeshAccess()->FindElementOfPoint(x0,ip,true);
             auto & felx = spaces[0]->GetFE(ElementId(elnr),lh);
             FlatVector<> shapex(felx.GetNDof(),lh);
             dynamic_cast<const BaseScalarFiniteElement &>(felx).CalcShape(ip,shapex);
             FlatVector<> val(tpfes->GetDimension(),lh);
             val = 0.0;
             int index = tpfes->GetIndex(elnr,0);
             Array<int> dnums;
             for(int i=index;i<index+spaces[1]->GetMeshAccess()->GetNE();i++)
             {
               auto & fely = spaces[1]->GetFE(ElementId(i-index),lh);
               tpfes->GetDofNrs(i,dnums);
               int tpndof = felx.GetNDof()*fely.GetNDof();
               FlatVector<> elvec(tpndof*tpfes->GetDimension(),lh);
               gf_tp->GetElementVector(dnums,elvec);
               FlatMatrix<> coefmat(felx.GetNDof(),fely.GetNDof()*tpfes->GetDimension(), &elvec(0));
               FlatMatrix<> coefyasmat(fely.GetNDof(),tpfes->GetDimension(),lh);
               // FlatVector<> coefy(fely.GetNDof()*tpfes->GetDimension(),lh);
               coefyasmat.AsVector() = Trans(coefmat)*shapex;
               const IntegrationRule & ir = SelectIntegrationRule(fely.ElementType(),2*fely.Order());
               BaseMappedIntegrationRule & mir = spaces[1]->GetMeshAccess()->GetTrafo(ElementId(i-index),lh)(ir,lh);
               FlatMatrix<> coefvals(ir.Size(), tpfes->GetDimension(),lh);
               coef->Evaluate(mir,coefvals);
               FlatMatrix<> shapesy(fely.GetNDof(),ir.Size(),lh);
               dynamic_cast<const BaseScalarFiniteElement & >(fely).CalcShape(ir,shapesy);
               FlatMatrix<> helpermat(ir.Size(),tpfes->GetDimension(),lh);
               helpermat = Trans(shapesy)*coefyasmat;
               for(int ip=0;ip<ir.Size();ip++)
                 for(int k=0;k<tpfes->GetDimension();k++)
                   val(k)+=helpermat(ip,k)*mir[ip].GetWeight()*coefvals(ip,k); // This still uses only the first coefficient!!!
             }
             double return_val = 0.0;
             for(int j: Range(tpfes->GetDimension()))
               return_val+=val(j);
             return return_val;
           },
        py::call_guard<py::gil_scoped_release>());
   m.def("TensorProductIntegrate",[](shared_ptr<GF> gf_tp, shared_ptr<GF> gf_x, spCF coef)
           {
             static Timer tall("comp.TensorProductIntegrate - total domain integral"); RegionTimer rall(tall);
             BaseVector & vec_in = gf_tp->GetVector();
             BaseVector & vec_out = gf_x->GetVector();
             LocalHeap clh(100000000,"TensorProductIntegrate");
             shared_ptr<TPHighOrderFESpace> tpfes = dynamic_pointer_cast<TPHighOrderFESpace>(gf_tp->GetFESpace());
             const Array<shared_ptr<FESpace> > & spaces = tpfes->Spaces(0);
             int ndofxspace = spaces[0]->GetNDof();
             auto & meshy = spaces[1]->GetMeshAccess();
             FlatVector<> elvec_out(ndofxspace*tpfes->GetDimension(),clh);
             FlatMatrix<> elvec_outmat(meshy->GetNE(),ndofxspace*tpfes->GetDimension(),clh);
             elvec_outmat = 0.0;
             elvec_out = 0.0;
            auto & element_coloring1 = spaces[1]->ElementColoring(VOL);
            for (FlatArray<int> els_of_col : element_coloring1)
            {
              SharedLoop sl(els_of_col.Range());
              task_manager -> CreateJob
              ( [&] (const TaskInfo & ti) 
              {
                LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
                for (int mynr : sl)
                {
                  HeapReset hr(lh);
                  int i = els_of_col[mynr];
                  auto & fely = spaces[1]->GetFE(ElementId(i),lh);
                  int ndofy = fely.GetNDof();
                  FlatMatrix<> elvec_slicemat(ndofy,ndofxspace*tpfes->GetDimension(),lh);
                  Array<int> dnumsslice(ndofy*ndofxspace, lh);
                  tpfes->GetSliceDofNrs(ElementId(i),0,dnumsslice,lh);
                  vec_in.GetIndirect(dnumsslice, elvec_slicemat.AsVector());
                  ElementTransformation & trafo = spaces[1]->GetMeshAccess()->GetTrafo(ElementId(i),lh);
                  const IntegrationRule & ir = SelectIntegrationRule(fely.ElementType(),2*fely.Order());
                  FlatMatrix<> shape(fely.GetNDof(),ir.Size(),lh);
                  dynamic_cast<const BaseScalarFiniteElement &>(fely).CalcShape(ir,shape);
                  BaseMappedIntegrationRule & mir = trafo(ir,lh);
                  FlatMatrix<> vals(mir.Size(), tpfes->GetDimension(),lh);
                  if(coef)
                    coef->Evaluate(mir, vals);
                  else
                    vals = 1.0;
                  for(int s=0;s<ir.Size();s++)
                    vals.Row(s)*=mir[s].GetWeight();
                  
                  int firstxdof = 0;
                  for(int j=0;j<spaces[0]->GetMeshAccess()->GetNE();j++)
                  {
                    int ndofx = spaces[0]->GetFE(ElementId(j),lh).GetNDof();
                    IntRange dnumsx(firstxdof, firstxdof+ndofx*tpfes->GetDimension());
                    FlatMatrix<> coefmat(ndofy,ndofx*tpfes->GetDimension(),lh);
                    coefmat = elvec_slicemat.Cols(dnumsx);
                    FlatMatrix<> tempmat(ndofx*tpfes->GetDimension(),ir.Size(),lh);
                    tempmat = Trans(coefmat)*shape;
                    for(int s=0;s<ir.Size();s++)
                    {
                      for( int dof : Range(ndofx) )
                        for(int d: Range(tpfes->GetDimension()) )
                          tempmat(dof*tpfes->GetDimension()+d,s)*=vals(s,d);
                      elvec_outmat.Cols(dnumsx).Row(i)+=(tempmat.Col(s));
                    }
                    firstxdof+=ndofx*tpfes->GetDimension();
                  }

                }
              }
              );
            }
            for(int i=0;i<elvec_outmat.Height();i++)
              elvec_out+=elvec_outmat.Row(i);
            int firstxdof = 0;
            // In case the x gridfunction and the Tensor space have equal number
            // of components write the integral of each component into the component
            // x gridfunction
            if(tpfes->GetDimension() == gf_x->GetFESpace()->GetDimension())
              for(int i=0;i<spaces[0]->GetMeshAccess()->GetNE();i++)
              {
                int ndofx = spaces[0]->GetFE(ElementId(i),clh).GetNDof();
                IntRange dnumsx(firstxdof, firstxdof+ndofx*tpfes->GetDimension());
                firstxdof+=ndofx*tpfes->GetDimension();
                Array<int> dofsx;
                spaces[0]->GetDofNrs(ElementId(i),dofsx);
                vec_out.SetIndirect(dofsx,elvec_out.Range(dnumsx));
              }
            // In case the x gridfunction has 1 component and the Tensor space has more than
            // one component, write the sum of the integral of each component into the x gridfunction
            else if(tpfes->GetDimension() > 1 && gf_x->GetFESpace()->GetDimension() == 1)
            {
              FlatVector<> elvec_sum(gf_x->GetFESpace()->GetNDof(),clh);
              elvec_sum = 0.0;
              for(int i: Range(tpfes->GetDimension()) )
              {
                SliceVector<> elvec_comp(gf_x->GetFESpace()->GetNDof(), tpfes->GetDimension(), &elvec_out(i));
                elvec_sum+=elvec_comp;
              }
              for(int i=0;i<spaces[0]->GetMeshAccess()->GetNE();i++)
              {
                int ndofx = spaces[0]->GetFE(ElementId(i),clh).GetNDof();
                IntRange dnumsx(firstxdof, firstxdof+ndofx);
                firstxdof+=ndofx;
                Array<int> dofsx;
                spaces[0]->GetDofNrs(ElementId(i),dofsx);
                vec_out.SetIndirect(dofsx,elvec_sum.Range(dnumsx));
              }
            }
   },
   py::arg("gftp"),
   py::arg("gfx"),
   py::arg("weight")=nullptr,
        py::call_guard<py::gil_scoped_release>()
   );

   m.def("ProlongateCoefficientFunction", [](spCF cf_x, int prolongateto, shared_ptr<FESpace> tpfes) -> shared_ptr<CoefficientFunction>
           {
             int dimx = dynamic_pointer_cast<TPHighOrderFESpace>(tpfes)->Spaces(0)[0]->GetMeshAccess()->GetDimension();
             int dimy = dynamic_pointer_cast<TPHighOrderFESpace>(tpfes)->Spaces(0)[1]->GetMeshAccess()->GetDimension();
             auto pcf = make_shared<ProlongateCoefficientFunction>(cf_x,prolongateto,cf_x->Dimension(),dimx,dimy,false);
             pcf->SetDimension(pcf->Dimension());
             return pcf;
           },
        py::call_guard<py::gil_scoped_release>());
   m.def("Prolongate", [](shared_ptr<GF> gf_x, shared_ptr<GF> gf_tp )
            {
              static Timer tall("comp.Prolongate"); RegionTimer rall(tall);
              shared_ptr<TPHighOrderFESpace> tpfes = dynamic_pointer_cast<TPHighOrderFESpace>(gf_tp->GetFESpace());
              LocalHeap lh(100000,"ProlongateFromXSpace");
              if(gf_x->GetFESpace() == tpfes->Space(-1) )
                tpfes->ProlongateFromXSpace(gf_x,gf_tp,lh);
              else
                cout << "GridFunction gf_x is not defined on first space"<<endl;
            },
         py::call_guard<py::gil_scoped_release>());
   m.def("Transfer2StdMesh", [](const shared_ptr<GF> gfutp,
                                shared_ptr<GF> gfustd)
            {
              static Timer tall("comp.Transfer2StdMesh"); RegionTimer rall(tall);
              Transfer2StdMesh(gfutp.get(),gfustd.get(),glh);
              return;
             },
             py::arg("gftp"),
             py::arg("gfstd"),
         py::call_guard<py::gil_scoped_release>()
        );
   
   m.def("Transfer2StdMesh", [](const spCF cftp, shared_ptr<GF> gfustd )
            {
              cout << cftp << endl;
              static Timer tall("comp.Transfer2StdMesh"); RegionTimer rall(tall);
              return;
             });

   py::class_<BaseVTKOutput, shared_ptr<BaseVTKOutput>>(m, "VTKOutput")
    .def(py::init([] (shared_ptr<MeshAccess> ma, py::list coefs_list,
                      py::list names_list, string filename, int subdivision, 
                      int only_element, string floatsize, bool legacy)
         -> shared_ptr<BaseVTKOutput>
         {
           Array<shared_ptr<CoefficientFunction> > coefs
             = makeCArray<shared_ptr<CoefficientFunction>> (coefs_list);
           Array<string > names
             = makeCArray<string> (names_list);
           shared_ptr<BaseVTKOutput> ret;
           if (ma->GetDimension() == 2)
             ret = make_shared<VTKOutput<2>> (ma, coefs, names, filename, subdivision, only_element, floatsize, legacy);
           else
             ret = make_shared<VTKOutput<3>> (ma, coefs, names, filename, subdivision, only_element, floatsize, legacy);
           return ret;
         }),
         py::arg("ma"),
         py::arg("coefs")= py::list(),
         py::arg("names") = py::list(),
         py::arg("filename") = "vtkout",
         py::arg("subdivision") = 0,
         py::arg("only_element") = -1,
         py::arg("floatsize") = "double",
         py::arg("legacy") = false,
         docu_string(R"raw_string(
VTK output class. Allows to put mesh and field information of several CoefficientFunctions into a VTK file.
(Can be used by independent visualization software, e.g. ParaView).

When run in parallel, rank 0 stores no vtk output, but writes the pvd-file that links all parallel
output together.

Parameters:

ma : ngsolve mesh
  mesh (Note: if a deformation is set, the output will be w.r.t. the deformed state of the mesh)

coefs: list of CoefficientFunctions
  list of CFs that are stored as fields in the Paraview output

names : list of strings
  labels for the fields that are put in the output file

filename : string (default: \"output\")
  name of the output file ( .vtu file ending is added or .vtk file ending is added (legacy mode) ).
  If run in parallel, the suffix \"_procxyz\" is added (xyz a number). 
  If output is written several times, the ending \"_stepxyz\" is added (xyz a counter). 
  If run in parallel or the output is called several times a meta file with ending .pvd is also generated for convenience.

subdivision : int
  Number of subdivision (bisections in each direction) that are applied
  (Note that only vertex values are stored otherwise rendering the output information piecewise linear only)

only_element : int
  only work on one specific element (default: -1 which means `draw all elements`)

floatsize : string in {\"single\", \"double\" }   object
  defines the precision of the output data (default is \"double\", \"single\" can be used to reduce output)

legacy : bool (default: False)
  defines if legacy-VTK output shall be used 
            .)raw_string")
         )
     .def("Do", [](shared_ptr<BaseVTKOutput> self, double time, VorB vb)
          { 
            self->Do(glh,time, vb);
            return self->lastoutputname;
          },
          py::arg("time")=-1,
          py::arg("vb")=VOL,
          py::call_guard<py::gil_scoped_release>(),
         docu_string(R"raw_string(
Write mesh and fields to file. When called several times on the same object
an index is added to the output file name. A meta file (.pvd) is written 
(unless in legacy mode).

Returns string of the output filename.

Parameters:

time : 
  associate a time to the current output

vb: VOL_or_BND (default VOL)
  defines if output is done on the volume (VOL) or surface mesh (BND).
            .)raw_string")          
          )
     .def("Do", [](shared_ptr<BaseVTKOutput> self, double time, VorB vb, const BitArray * drawelems)
          { 
            self->Do(glh,time, vb, drawelems);
            return self->lastoutputname;
          },
          py::arg("time")=-1,
          py::arg("vb")=VOL,
          py::arg("drawelems"),
          py::call_guard<py::gil_scoped_release>(),
         docu_string(R"raw_string(
Write mesh and fields to file. When called several times on the same object
an index is added to the output file name. A meta file (.pvd) is written 
(unless in legacy mode).

Returns string of the output filename.

Parameters:

time : 
  associate a time to the current output (default: output counter)

vb: VOL_or_BND (default VOL)
  defines if output is done on the volume (VOL) or surface mesh (BND).

drawelems: BitArray
  defines the submesh (set of elements) that are (only) used for drawing. 
            .)raw_string")          
          )
     ;
   
   m.def("PatchwiseSolve",
         [&] (shared_ptr<SumOfIntegrals> bf,
              shared_ptr<SumOfIntegrals> lf,
              shared_ptr<GridFunction> gf)
         { PatchwiseSolve(bf, lf, gf, glh); },
         py::arg("bf"), py::arg("lf"), py::arg("gf"));

   py::class_<InterpolateProxy, shared_ptr<InterpolateProxy>, ProxyFunction> (m, "InterpolateProxy");
   m.def("Interpolate", 
         [] (shared_ptr<CoefficientFunction> cf, shared_ptr<FESpace> fes, int bonus_intorder)
         {
           if (!fes)
             throw Exception("In Interpolate: invalid space");
           return InterpolateCF(cf, fes, bonus_intorder);
         }, py::arg("cf"), py::arg("space"), py::arg("bonus_intorder")=0,
         docu_string(R"raw_string(Interpolate a CoefficientFunction into the finite element space.
The interpolation is canonical interpolation using dual shapes.
The result is a CoefficientFunction.
Interpolation is done on the fly for each element, no global GridFunction is allocated.)raw_string")
          );
    
   
   m.def("ConvertOperator", [&](shared_ptr<FESpace> spacea, shared_ptr<FESpace> spaceb,
				shared_ptr<ProxyFunction> trial_proxy, shared_ptr<CoefficientFunction> trial_cf,
				optional<Region> definedon, VorB vb, shared_ptr<BitArray> range_dofs, bool localop, bool parmat, bool use_simd,
				int bonus_io_ab, int bonus_io_bb, bool geom_free) -> shared_ptr<BaseMatrix> {

	   const Region* reg = NULL;
	   if( definedon.has_value() ) {
	     reg = &(*definedon);
	     vb = VorB(*definedon);
	   }

	   shared_ptr<BaseMatrix> op;

	   if (trial_proxy != nullptr) {
	     if ( !trial_proxy->IsTrialFunction() )
	       { throw Exception("Need a trial-proxy, but got a test-proxy!"); }
	     shared_ptr<DifferentialOperator> eval;
	     if ( vb == VOL )
	       { eval = trial_proxy->Evaluator(); }
	     else if ( vb == BND )
	       { eval = trial_proxy->TraceEvaluator(); }
	     else if ( vb == BBND )
	       { eval = trial_proxy->TTraceEvaluator(); }
	     else
	       { throw Exception("ProxyFunction has no BBBND evaluator!"); }
	     if ( eval == nullptr )
	       { throw Exception(string("trial-proxy has no evaluator vor vb = ") + to_string(vb) + string("!")); }
	     op = ConvertOperator(spacea, spaceb, vb, glh, eval, trial_cf, reg, range_dofs, localop, parmat, use_simd, bonus_io_ab, bonus_io_bb, geom_free);
	   }
	   else
	     { op = ConvertOperator(spacea, spaceb, vb, glh, nullptr, trial_cf, reg, range_dofs, localop, parmat, use_simd, bonus_io_ab, bonus_io_bb, geom_free); }

	   return op;
	 },
	 py::arg("spacea"), py::arg("spaceb"),
	 py::arg("trial_proxy") = nullptr,
	 py::arg("trial_cf") = nullptr,
	 py::arg("definedon") = nullptr,
	 py::arg("vb") = VOL,
	 py::arg("range_dofs") = nullptr,
	 py::arg("localop") = false,
	 py::arg("parmat") = true,
	 py::arg("use_simd") = true,
	 py::arg("bonus_intorder_ab") = 0,
	 py::arg("bonus_intorder_bb") = 0,
	 py::arg("geom_free") = false,
     docu_string(R"raw_string(
A conversion operator between FESpaces. Embedding if spacea is a subspace of spaceb, otherwise an interpolation operator defined by element-wise application of dual shapes (and averaging between elements).

Parameters:

spacea: ngsolve.comp.FESpace
  the origin space

spaceb: ngsolve.comp.FESpace
  the goal space

trial_proxy: ngsolve.comp.ProxyFunction
  (optional) Must be a trial-proxy on spacea. If given, instead of a FE-function funca from spacea, the operator converts trial_proxy(funca) to spaceb.

trial_proxy: ngsolve.comp.CoefficientFunction
  (optional) Same as trial_proxy, but takes any CoefficientFunction. Use at your own peril.

definedon: object
  what part of the domain to restrict the operator to

vb: ngsolve.comp.VorB
  what kind of co-dimension elements to convert on VOL, BND, BBND, ...

range_dofs: ngsolve.ngstd.BitArray
  Projects out DOFs in the range where range_dofs are not set

localop: bool
  True -> do not average across MPI boundaries. No effect for non MPI-paralell space. Use carefully!!

parmat: bool
  If True, returns a ParallelMatrix for MPI-parallel spaces. If False, or for non MPI-parallel spaces, returns a local BaseMatrix.

use_simd:
  False -> Do not use SIMD for setting up the Matrix. (for debugging purposes).

bonus_intorder_ab/bb: int
  Bonus integration order for spacea/spaceb and spaceb/spaceb integrals. Can be useful for curved elements. Should only be necessary for
spacea/spaceb integrals.

geom_free:
  If True, assembles a matrix-free operator.
)raw_string")
	 );

   m.def("MPI_Init", [&]()
	 {
	   const char * progname = "ngslib";
	   typedef const char * pchar;
	   pchar ptrs[2] = { progname, nullptr };
	   pchar * pptr = &ptrs[0];
          
	   static MyMPI mympi(1, (char**)pptr);
	   return NgMPI_Comm(MPI_COMM_WORLD);
	 });

   py::class_<ContactBoundary, shared_ptr<ContactBoundary>>
     (m, "ContactBoundary")
     .def(py::init([](shared_ptr<FESpace> fes, Region master, Region minion,
                      bool draw_pairs)
     {
       cout << "WARNING: ContactBoundary constructor with FESpace is deprecated, fes will be set correctly in Update!" << endl;
       return make_shared<ContactBoundary>(master, minion, draw_pairs);
     }), "fes"_a, "master"_a, "minion"_a, "draw_pairs"_a = false)
     .def(py::init<Region, Region, bool>(),
          R"delimiter(
Class for managing contact interfaces.
The created object must be kept alive in python as long as
operations of it are used!
)delimiter", "master"_a, "minion"_a, "draw_pairs"_a=false)
     .def("AddEnergy", &ContactBoundary::AddEnergy,
          "form"_a, "deformed"_a = false)
     .def("AddIntegrator", &ContactBoundary::AddIntegrator,
          "form"_a, "deformed"_a = false)
     .def("Update", &ContactBoundary::Update,
          py::arg("gf") = nullptr, py::arg("bf") = nullptr,
          py::arg("intorder") = 4, py::arg("maxdist") = 0.,
          R"delimiter(
Update searchtree for gap function.
If bf is given add specialelements corresponding to
integrationrules of order 'intorder' on each master
element to BilinearForm bf.
`maxdist` is the maximum distance where this function is accurate.
If `maxdist` == 0. then 2*meshsize is used.
)delimiter")
     .def_property_readonly("gap", &ContactBoundary::Gap)
     .def_property_readonly("normal", &ContactBoundary::Normal)
     .def("_GetWebguiData", [] (shared_ptr<ContactBoundary> contact) {
             auto [primary_points, secondary_points] = contact->GetDrawingPairs();
             std::vector<double> p;
             p.reserve(primary_points.Size()*6);

             for(auto i : Range(primary_points))
             {
                 for(auto d : Range(3))
                     p.push_back(primary_points[i][d]);
                 for(auto d : Range(3))
                     p.push_back(secondary_points[i][d]);
             }

             py::dict gap_data;
             gap_data["type"] = "lines";
             gap_data["color"] = "black";
             gap_data["name"] = "Contact Pairs";
             gap_data["position"] = py::cast(p);
             return gap_data;
     })
     ;

  /////////////////////////////////////////////////////////////////////////////////////
}

#endif // NGS_PYTHON
