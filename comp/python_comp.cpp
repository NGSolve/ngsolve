#ifdef NGS_PYTHON

#include <regex>

#include "../ngstd/python_ngstd.hpp"
#include <comp.hpp>

#include "hdivdivfespace.hpp"
#include "hdivdivsurfacespace.hpp"
#include "hcurlcurlfespace.hpp"
#include "numberfespace.hpp"
#include "compressedfespace.hpp"
using namespace ngcomp;

using ngfem::ELEMENT_TYPE;

typedef GridFunction GF;




/*
static size_t global_heapsize = 1000000;
static LocalHeap glh(global_heapsize, "python-comp lh", true);
*/

/*
template <> class cl_NonElement<ElementId>
{
public:
  static ElementId Val() { return ElementId(VOL,-1); }
};
*/



/*
namespace pybind11
{
  // would like to replace MakePyTuple by this cast, but doesn't work
  template <typename T>
  py::tuple cast (const BaseArrayObject<T> & ao)
  {
    size_t s = ao.Size();
    py::tuple tup(s);
    for (size_t i = 0; i < s; i++)
      tup[i] = ao[i];
    return tup;
  }
}
*/


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


py::object MakeProxyFunction2 (shared_ptr<FESpace> fes,
                               bool testfunction,
                               const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock)
{
  auto compspace = dynamic_pointer_cast<CompoundFESpace> (fes);
  if (compspace && !fes->GetEvaluator())
    {
      py::list l;
      int nspace = compspace->GetNSpaces();
      for (int i = 0; i < nspace; i++)
        {
          l.append (MakeProxyFunction2 ((*compspace)[i], testfunction,
                                         [&] (shared_ptr<ProxyFunction> proxy)
                                         {
                                           auto block_eval = make_shared<CompoundDifferentialOperator> (proxy->Evaluator(), i);
                                           shared_ptr <CompoundDifferentialOperator> block_deriv_eval = nullptr;
                                           if (proxy->DerivEvaluator() != nullptr)
                                             block_deriv_eval = make_shared<CompoundDifferentialOperator> (proxy->DerivEvaluator(), i);
                                           shared_ptr <CompoundDifferentialOperator> block_trace_eval = nullptr;
                                           if (proxy->TraceEvaluator() != nullptr)
                                             block_trace_eval = make_shared<CompoundDifferentialOperator> (proxy->TraceEvaluator(), i);
					   shared_ptr <CompoundDifferentialOperator> block_ttrace_eval = nullptr;
					   if (proxy->TTraceEvaluator() != nullptr)
					     block_ttrace_eval = make_shared<CompoundDifferentialOperator> (proxy->TTraceEvaluator(),i);
                                           shared_ptr <CompoundDifferentialOperator> block_trace_deriv_eval = nullptr;
                                           if (proxy->TraceDerivEvaluator() != nullptr)
                                             block_trace_deriv_eval = make_shared<CompoundDifferentialOperator> (proxy->TraceDerivEvaluator(), i);
					   shared_ptr <CompoundDifferentialOperator> block_ttrace_deriv_eval = nullptr;
					   if (proxy->TTraceDerivEvaluator() != nullptr)
					     block_ttrace_deriv_eval = make_shared<CompoundDifferentialOperator> (proxy->TTraceDerivEvaluator(),i);
                                           auto block_proxy = make_shared<ProxyFunction> (fes, testfunction, fes->IsComplex(),                                                                                          block_eval, block_deriv_eval, block_trace_eval, block_trace_deriv_eval,
					  block_ttrace_eval, block_ttrace_deriv_eval);

                                           SymbolTable<shared_ptr<DifferentialOperator>> add = proxy->GetAdditionalEvaluators();
                                           for (int j = 0; j < add.Size(); j++)
                                             block_proxy->SetAdditionalEvaluator(add.GetName(j),
                                                                                 make_shared<CompoundDifferentialOperator> (add[j], i));
                                           
                                           block_proxy = addblock(block_proxy);
                                           return block_proxy;
                                         }));
        }
      return l;
    }

  auto proxy = make_shared<ProxyFunction>  (fes, testfunction, fes->IsComplex(),
                                            fes->GetEvaluator(),
                                            fes->GetFluxEvaluator(),
                                            fes->GetEvaluator(BND),
                                            fes->GetFluxEvaluator(BND),
					    fes->GetEvaluator(BBND),
					    fes->GetFluxEvaluator(BBND));
  auto add_diffops = fes->GetAdditionalEvaluators();
  for (int i = 0; i < add_diffops.Size(); i++)
    proxy->SetAdditionalEvaluator (add_diffops.GetName(i), add_diffops[i]);

  proxy = addblock(proxy);
  return py::cast(proxy);
}

py::object MakeProxyFunction (shared_ptr<FESpace> fes,
                              bool testfunction) 
{
  return 
    MakeProxyFunction2 (fes, testfunction, 
                        [&] (shared_ptr<ProxyFunction> proxy) { return proxy; });
}



class GlobalDummyVariables 
{
public:
  int GetMsgLevel() { return printmessage_importance; }
  void SetMsgLevel(int msg_level) 
  {
    // cout << "set printmessage_importance to " << msg_level << endl;
    printmessage_importance = msg_level; 
    netgen::printmessage_importance = msg_level; 
  }
  string GetTestoutFile () const
  {
    return "no-filename-here";
  }
  void SetTestoutFile(string filename) 
  {
    // cout << "set testout-file to " << filename << endl;
    testout = new ofstream(filename);
  }
  
};
static GlobalDummyVariables globvar;


auto fesPickle = [](const FESpace& fes)
{
  auto flags = fes.GetFlags();
  auto mesh = fes.GetMeshAccess();
  auto type = fes.type;
  // TODO: pickle order policies
  return py::make_tuple(type,mesh,flags);
};


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
  auto pyspace = py::class_<FES, shared_ptr<FES>,BASE> (m, pyname.c_str());
  pyspace
    .def(py::init([pyspace](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(pyspace, kwargs, info);
                    auto fes = make_shared<FES>(ma,flags);
                    LocalHeap glh(10000000, "init-fes-lh");                    
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    
    .def(py::pickle(fesPickle,
                    (shared_ptr<FES>(*)(py::tuple)) fesUnpickle<FES>))
    ;
  
  return pyspace;
}




void ExportNgcompMesh (py::module &m);

void NGS_DLL_HEADER ExportNgcomp(py::module &m)
{

  ExportNgcompMesh(m);
  //////////////////////////////////////////////////////////////////////////////////////////

  static size_t global_heapsize = 1000000;
  static LocalHeap glh(global_heapsize, "python-comp lh", true);

  

  //////////////////////////////////////////////////////////////////////////////////////////

  py::enum_<COUPLING_TYPE> (m, "COUPLING_TYPE", docu_string(R"raw_string(
Enum specifying the coupling type of a degree of freedom, each dof is
either UNUSED_DOF, LOCAL_DOF, INTERFACE_DOF or WIREBASKET_DOF, other values
are provided as combinations of these:

UNUSED_DOF: Dof is not used, i.e the slave dofs in a :any:`Periodic` finite
    element space.

LOCAL_DOF: Inner degree of freedom, will be eliminated by static
    condensation and reconstructed afterwards.

HIDDEN_DOF: Inner degree of freedom, that will be eliminated by static
    condensation and *not* reconstruced afterwards(spares some entries).
    Note: 
     * without static condensation a HIDDEN_DOF is treated as any other
       DOF, e.g. as a LOCAL_DOF
     * To a HIDDEN_DOF the r.h.s. vector must have zero entries.

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
    .def(py::init<const MeshAccess&,VorB,IntRange>())
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

  
  py::enum_<ORDER_POLICY>(m, "ORDER_POLICY")
    .value("CONSTANT", CONSTANT_ORDER)
    .value("NODETYPE", NODE_TYPE_ORDER)
    .value("VARIABLE", VARIABLE_ORDER)
    .value("OLDSTYLE", OLDSTYLE_ORDER)
    ;


  //////////////////////////////////////////////////////////////////////////////////////////

  /*
    // where do we need that ? 
  py::class_<FlatArray<NodeId> > class_flatarrayNI (m, "FlatArrayNI");
  PyDefVector<FlatArray<NodeId>, NodeId>(m, class_flatarrayNI);
  PyDefToString<FlatArray<NodeId> >(m, class_flatarrayNI);
  class_flatarrayNI.def(py::init<size_t, NodeId *>());

  py::class_<Array<NodeId>, FlatArray<NodeId> >(m, "ArrayNI")
    .def(py::init<size_t>())
    ;
  */
  

  py::class_<FESpace::Element,Ngs_Element>(m, "FESpaceElement")
    .def_property_readonly("dofs",
                           [](FESpace::Element & el) 
                           { return MakePyList (el.GetDofs()); },
                           "degrees of freedom of element"
                           )

    /*
      // don't give LocalHeap to Python !
    .def("GetLH",[](FESpace::Element & el) -> LocalHeap & 
                                  {
                                    return el.GetLH();
                                  },
         py::return_value_policy::reference
         )
    */
    
    .def("GetFE",[](FESpace::Element & el)
         {
           return shared_ptr<FiniteElement>(const_cast<FiniteElement*>(&el.GetFE()), NOOP_Deleter);
         },
         py::return_value_policy::reference,
         "the finite element containing shape functions"
         )

    .def("GetTrafo",[](FESpace::Element & el)
         {
           return shared_ptr<ElementTransformation>(const_cast<ElementTransformation*>(&el.GetTrafo()), NOOP_Deleter);
         },
         py::return_value_policy::reference,
         "the transformation from reference element to physical element"
         )

    ;
  //////////////////////////////////////////////////////////////////////////////////////////


  py::class_<GlobalDummyVariables> (m, "GlobalVariables")
    .def_property("msg_level", 
                 &GlobalDummyVariables::GetMsgLevel,
                 &GlobalDummyVariables::SetMsgLevel)
    .def_property("testout", 
                 &GlobalDummyVariables::GetTestoutFile,
                 &GlobalDummyVariables::SetTestoutFile)
    /*
    .def_property_readonly("pajetrace",
		  &GlobalDummyVariables::GetTestoutFile,
		 [] (GlobalDummyVariables&, bool use)
				  { TaskManager::SetPajeTrace(use); });
    */
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
    // &GlobalDummyVariables::SetTestoutFile)
    ;

  m.attr("ngsglobals") = py::cast(&globvar);





  //////////////////////////////////////////////////////////////////////////////////////////
  
  py::class_<NGS_Object, shared_ptr<NGS_Object>>(m, "NGS_Object")
    // .def_property_readonly("name", [](const NGS_Object & self)->string { return self.GetName();})
    .def_property("name", &NGS_Object::GetName, &NGS_Object::SetName)
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
                     if (!self->Deriv()) return "";
                     return self->DerivEvaluator()->Name();
                   })
    .def("Operator",
         [] (const spProxy self, string name) -> py::object
          {
            auto op = self->GetAdditionalProxy(name);
            if (op)
              return py::cast(op);
            return py::none();
	  }, "Use an additional operator of the finite element space")
    ;



  /*
  struct OrderProxy 
  {
    FESpace & fes;
    OrderProxy (FESpace & afes) : fes(afes) { ; }
  };


  py::class_<OrderProxy> (m, "OrderProxy")
    .def("__setitem__",
         [] (OrderProxy & self, ElementId ei, int o) 
          {
            cout << "set order of el " << ei << " to order " << o << endl;
            cout << "(not implemented)" << endl;
          })

    .def("__setitem__",
         [] (OrderProxy & self, ELEMENT_TYPE et, int o) 
          {
            cout << "set order of eltype " << et << " to order " << o << endl;
            self.fes.SetBonusOrder (et, o - self.fes.GetOrder());

            LocalHeap lh (100000, "FESpace::Update-heap", true);
            self.fes.Update(lh);
            self.fes.FinalizeUpdate(lh);
          })
    
    .def("__setitem__",
         [] (OrderProxy & self, NODE_TYPE nt, int o) 
          {
            cout << "set order of nodetype " << int(nt) << " to order " << o << endl;
            nt = StdNodeType (nt, self.fes.GetMeshAccess()->GetDimension());
            cout << "canonical nt = " << int(nt) << endl;
            int bonus = o-self.fes.GetOrder();
            switch (nt)
              {
              case 1: 
                self.fes.SetBonusOrder(ET_SEGM, bonus); break;
              case 2: 
                self.fes.SetBonusOrder(ET_QUAD, bonus); 
                self.fes.SetBonusOrder(ET_TRIG, bonus); break;
              case 3: 
                self.fes.SetBonusOrder(ET_TET, bonus); 
                self.fes.SetBonusOrder(ET_PRISM, bonus);
                self.fes.SetBonusOrder(ET_PYRAMID, bonus);
                self.fes.SetBonusOrder(ET_HEX, bonus); break;
              default: ;
              }

            LocalHeap lh (100000, "FESpace::Update-heap", true);
            self.fes.Update(lh);
            self.fes.FinalizeUpdate(lh);
          })

    .def("__setitem__",
         [] (OrderProxy & self, NODE_TYPE nt, int nr, int o) 
          {
            cout << "set order of " << nt << " " << nr << " to " << o << endl;
            cout << "(not implemented)" << endl;
          })

    .def("__setitem__",
         [] (OrderProxy & self, py::tuple tup, int o) 
          {
            NODE_TYPE nt = py::extract<NODE_TYPE>(tup[0])();
            int nr = py::extract<int>(tup[1])();
            cout << "set order of " << nt << " " << nr << " to " << o << endl;
            cout << "(not implemented)" << endl;
          })
    ;
  */

  m.def("SetHeapSize",
        [](size_t heapsize)
        {
          if (heapsize > global_heapsize)
            {
              global_heapsize = heapsize;
              glh = LocalHeap (heapsize, "python-comp lh", true);
            }
        });
  
  m.def("SetTestoutFile",
        [](string filename)
        {
          testout = new ofstream (filename);
        }, "Enable some logging into file with given filename");


  //////////////////////////////////////////////////////////////////////////////////////////


  auto fes_class = py::class_<FESpace, shared_ptr<FESpace>>(m, "FESpace",
		    docu_string(R"raw_string(Finite Element Space

Provides the functionality for finite element calculations.

Some available FESpaces are:

H1
HCurl
HDiv
L2
FacetFESpace
HDivDiv

2 __init__ overloads:
  1) To create a registered FESpace
  2) To create a compound FESpace from multiple created FESpaces

1)

Parameters

type : string
  Type of the finite element space. This parameter is automatically
  set if the space is constructed with a generator function.

mesh : ngsolve.Mesh
  Mesh on which the finite element space is defined on.

kwargs : For a description of the possible kwargs have a look a bit further down.

2)

Parameters:

spaces : list of ngsolve.FESpace
  List of the spaces for the compound finite element space

kwargs : For a description of the possible kwargs have a look a bit further down.

)raw_string"), py::dynamic_attr());
  fes_class
    .def(py::init([fes_class] (py::list lspaces, py::kwargs kwargs)
                  {
		    py::list info;
		    auto flags = CreateFlagsFromKwArgs(fes_class, kwargs, info);
                    Array<shared_ptr<FESpace>> spaces;
                    for (auto fes : lspaces )
                      spaces.Append(py::extract<shared_ptr<FESpace>>(fes)());
                    if (spaces.Size() == 0)
                      throw Exception("Compound space must have at least one space");
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
                    shared_ptr<FESpace>
                      fes = make_shared<CompoundFESpace> (spaces[0]->GetMeshAccess(), spaces, flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
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
                    auto flags = CreateFlagsFromKwArgs(fes_class, kwargs, info);
                    auto fes = CreateFESpace (type, ma, flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),
                  py::arg("type"), py::arg("mesh"),
                  "allowed types are: 'h1ho', 'l2ho', 'hcurlho', 'hdivho' etc."
                  )

    .def_static("__flags_doc__", [] ()
         {
           return py::dict
             (
              py::arg("order") = "int = 1\n"
              "  order of finite element space",
              py::arg("complex") = "bool = False",
              py::arg("dirichlet") = "regexpr\n"
              "  Regular expression string defining the dirichlet boundary.\n"
              "  More than one boundary can be combined by the | operator,\n"
              "  i.e.: dirichlet = 'top|right'",
              py::arg("definedon") = "Region or regexpr\n"
              "  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n"
              "  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n"
              "  mesh.Materials('regexpr').",
              py::arg("dim") = "int = 1\n"
              "  Create multi dimensional FESpace (i.e. [H1]^3)",
              py::arg("dgjumps") = "bool = False\n"
              "  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n"
              "  since the dofs have a different coupling then and this changes the sparsity\n"
              "  pattern of matrices."
              );
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
                        if (py::isinstance<py::str>(dirichlet))
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

                        py::extract<Region> definedon_reg(definedon);
                        if (definedon_reg.check() && definedon_reg().IsVolume())
                          {
                            Array<double> defonlist;
                            for (int i = 0; i < definedon_reg().Mask().Size(); i++)
                              if (definedon_reg().Mask().Test(i))
                                defonlist.Append(i+1);
                            flags->SetFlag ("definedon", defonlist);
                          }
                        if (definedon_reg.check() && definedon_reg().VB()==BND)
                          {
                            Array<double> defonlist;
                            for (auto i : Range(definedon_reg().Mask().Size()))
                              if(definedon_reg().Mask().Test(i))
                                defonlist.Append(i+1);
                            flags->SetFlag("definedonbound", defonlist);
                          }
                      }),
                     py::arg("order_policy") = py::cpp_function
                     ([] (ORDER_POLICY op, Flags* flags, py::list info)
                      {
                        flags->SetFlag("order_policy", int(op));
                      })
                     );
                     return special;
                     })
    .def(py::pickle(fesPickle, (shared_ptr<FESpace>(*)(py::tuple)) fesUnpickle<FESpace>))
    .def("Update", [](shared_ptr<FESpace> self)
         { 
           self->Update(glh);
           self->FinalizeUpdate(glh);
         },
         "update space after mesh-refinement")
     .def("UpdateDofTables", [](shared_ptr<FESpace> self)
         {
           self->UpdateDofTables();
           self->UpdateCouplingDofArray();
           self->FinalizeUpdate(glh);
         },
         "update dof-tables after changing polynomial order distribution")
     .def("FinalizeUpdate", [](shared_ptr<FESpace> self)
         { 
           self->FinalizeUpdate(glh);
         },
         "finalize update")
    .def_property_readonly ("ndof", [](shared_ptr<FESpace> self) { return self->GetNDof(); },
                            "number of degrees of freedom")

    .def_property_readonly ("ndofglobal",
                            [](shared_ptr<FESpace> self) { return self->GetNDofGlobal(); },
                            "global number of dofs on MPI-distributed mesh")
    .def("__str__", [] (shared_ptr<FESpace> self) { return ToString(*self); } )
    .def("__timing__", [] (shared_ptr<FESpace> self) { return py::cast(self->Timing()); })
    .def_property_readonly("lospace", [](shared_ptr<FESpace> self) -> shared_ptr<FESpace>
			   { return self->LowOrderFESpacePtr(); })
    .def_property_readonly("mesh",
                           [](shared_ptr<FESpace> self) -> shared_ptr<MeshAccess>
                           { return self->GetMeshAccess(); })

    // .def_property_readonly("order", [] (shared_ptr<FESpace> self) { return OrderProxy(*self); },
    // "proxy to set order for individual nodes")
    .def_property_readonly("globalorder", [] (shared_ptr<FESpace> self) { return self->GetOrder(); },
                  "query global order of space")    
    .def_property_readonly("type", [] (shared_ptr<FESpace> self) { return self->type; },
                  "type of finite element space")

    .def_property_readonly("is_complex", &FESpace::IsComplex)

    .def("SetDefinedOn", [] (FESpace& self, Region& reg)
         {
           self.SetDefinedOn(reg.VB(),reg.Mask());
         }, py::arg("Region"))

    .def("SetOrder",
         [](shared_ptr<FESpace> self, ELEMENT_TYPE et, int order /*, py::object order_left, py::object order_right*/)
         {
           self->SetOrder (et, order);
           /*
           if (py::isinstance<py::int_> (order))
             {
               self->SetOrderLeft (et, order.cast<py::int_>());
               self->SetOrderRight (et, order.cast<py::int_>());
             }
           */
           /*
           if (py::isinstance<py::int_> (order_left))
             self->SetOrderLeft (et, order_left.cast<py::int_>());
           if (py::isinstance<py::int_> (order_right))
             self->SetOrderRight (et, order_right.cast<int>());
           */
         },
         py::arg("element_type"),
         py::arg("order")
         // py::arg("order_left")=DummyArgument(),
         // py::arg("order_right")=DummyArgument()
         )

    .def("SetOrder",
         [](shared_ptr<FESpace> self, NodeId ni, int order)
         {
           self->SetOrder(ni, order);
         },
         py::arg("nodeid"),
         py::arg("order")
         )
    
    .def("Elements", 
         [](shared_ptr<FESpace> self, VorB vb)
         { return FESpace::ElementRange(self->Elements(vb, glh)); },
         py::arg("VOL_or_BND")=VOL)

    .def("GetDofNrs", [](shared_ptr<FESpace> self, ElementId ei)
         {
           Array<DofId> tmp; self->GetDofNrs(ei,tmp);
           return MakePyTuple(tmp);           
         })

    .def("GetDofNrs", [](shared_ptr<FESpace> self, NodeId ni)
         {
           Array<DofId> tmp; self->GetDofNrs(ni,tmp);
           return MakePyTuple(tmp);
         })

    .def ("GetDofs", [](shared_ptr<FESpace> self, Region reg)
          {
            return self->GetDofs(reg);
          })
    
    .def("CouplingType", [](shared_ptr<FESpace> self, DofId dofnr) -> COUPLING_TYPE
         { return self->GetDofCouplingType(dofnr); },
         py::arg("dofnr"),
         "get coupling type of a degree of freedom"
         )
    .def("SetCouplingType", [](shared_ptr<FESpace> self, DofId dofnr, COUPLING_TYPE ct)
         { self->SetDofCouplingType(dofnr,ct); },
         py::arg("dofnr"), py::arg("coupling_type"),
         "set coupling type of a degree of freedom"
         )

    .def ("GetFE", [](shared_ptr<FESpace> self, ElementId ei) -> py::object
          {
            Allocator alloc;
            
            auto fe = shared_ptr<FiniteElement> (&self->GetFE(ei, alloc)); 
            
            auto scalfe = dynamic_pointer_cast<BaseScalarFiniteElement> (fe);
            if (scalfe) return py::cast(scalfe);

            auto hcurlfe = dynamic_pointer_cast<BaseHCurlFiniteElement> (fe);
            if (hcurlfe) return py::cast(hcurlfe);
            
            return py::cast(fe);
          })

    /*
    .def ("GetFE", [](shared_ptr<FESpace> self, ElementId ei, LocalHeap & lh)
          {
            return shared_ptr<FiniteElement>(&self->GetFE(ei, lh), NOOP_Deleter);
          },
          py::return_value_policy::reference)
    */
    
    .def("FreeDofs",
         [] (const shared_ptr<FESpace>self, bool coupling)
         { return self->GetFreeDofs(coupling); },
         py::arg("coupling")=false,
         "Return BitArray of free (non-Dirichlet) dofs\n"
         "coupling=False ... all free dofs including local dofs\n"
         "coupling=True .... only element-boundary free dofs"
         )

    .def("ParallelDofs",
         [] (const shared_ptr<FESpace>self)
         { return self->GetParallelDofs(); },
         "Return dof-identification for MPI-distributed meshes")

    .def("Range",
         [] (const shared_ptr<FESpace> self, int comp) -> py::slice
         {
           auto compspace = dynamic_pointer_cast<CompoundFESpace> (self);
           if (!compspace)
             throw py::type_error("'Range' is available only for product spaces");
           IntRange r = compspace->GetRange(comp);
           return py::slice(py::int_(r.First()), py::int_(r.Next()),1);
         },
         "component"_a,
         "Return interval of dofs of a component of a product space")
    
    .def_property_readonly("components", 
                  [](shared_ptr<FESpace> self)-> py::tuple
                   { 
                     auto compspace = dynamic_pointer_cast<CompoundFESpace> (self);
                     if (!compspace)
                       throw py::type_error("'components' is available only for product spaces");
                     py::tuple vecs(compspace->GetNSpaces());
                     for (int i = 0; i < compspace -> GetNSpaces(); i++) 
                       vecs[i]= py::cast((*compspace)[i]);
                     return vecs;
                   },
                  "Return a list of the components of a product space")

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

    .def("SolveM",
         [] (const shared_ptr<FESpace> self,
             BaseVector& vec, spCF rho)
         { self->SolveM(rho.get(), vec, glh); },
         py::arg("vec"), py::arg("rho")=nullptr,
         "Solve with the mass-matrix. Available only for L2-like spaces")
    .def("ApplyM",
         [] (const shared_ptr<FESpace> self,
             BaseVector& vec, spCF rho)
         { self->ApplyM(rho.get(), vec, glh); },
         py::arg("vec"), py::arg("rho")=nullptr,
         "Apply mass-matrix. Available only for L2-like spaces")
        
    .def("__eq__",
         [] (shared_ptr<FESpace> self, shared_ptr<FESpace> other)
         {
           return self == other;
         })
    ;

  py::class_<CompoundFESpace, shared_ptr<CompoundFESpace>, FESpace>
    (m,"CompoundFESpace")
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
                      fes->Update(lh);
                      fes->FinalizeUpdate(lh);
                      py::cast(fes).attr("__dict__") = state[2];
                      return fes;
                    }))
    ;

  
  ExportFESpace<HCurlHighOrderFESpace> (m, "HCurl")
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                      attr("FESpace").
                                                      attr("__flags_doc__")());
                  flags_doc["nograds"] = "bool = False\n"
                    "  Remove higher order gradients of H1 basis functions from HCurl FESpace";
                  flags_doc["type1"] = "bool = False\n"
                    "  Use type 1 Nedelec elements";
                  flags_doc["discontinuous"] = "bool = False\n"
                    "  Create discontinuous HCurl space";
                  return flags_doc;
                })
    .def("CreateGradient", [](shared_ptr<HCurlHighOrderFESpace> self) {
        auto fesh1 = self->CreateGradientSpace();
        shared_ptr<BaseMatrix> grad = self->CreateGradient(*fesh1);
        return py::make_tuple(grad, shared_ptr<FESpace>(fesh1));
      })
    ;
  
  ExportFESpace<HDivHighOrderFESpace> (m, "HDiv")
    .def_static("__flags_doc__", [] ()
              {
                auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                    attr("FESpace").
                                                    attr("__flags_doc__")());
                flags_doc["discontinuous"] = "bool = False\n"
                  "  Create discontinuous HDiv space";
                flags_doc["hodivfree"] = "bool = False\n"
                  "  Remove high order element bubbles with non zero divergence";
                flags_doc["highest_order_dc"] = "bool = False\n"
                  "  Activates relaxed H(div)-conformity. Allows normal discontinuity of highest order facet basis functions";
                return flags_doc;
              })
    .def("Average", &HDivHighOrderFESpace::Average,
          py::arg("vector"))
     ;
  
  auto h1 = ExportFESpace<H1HighOrderFESpace> (m, "H1");

  auto vectorh1 = ExportFESpace<VectorH1FESpace, CompoundFESpace> (m, "VectorH1");
 
  auto l2 = ExportFESpace<L2HighOrderFESpace> (m, "L2");

  auto vectorl2 = ExportFESpace<VectorL2FESpace, CompoundFESpace> (m, "VectorL2");

  auto l2surface = ExportFESpace<L2SurfaceHighOrderFESpace> (m, "SurfaceL2");

  auto numberfes = ExportFESpace<NumberFESpace> (m, "NumberSpace");

  ExportFESpace<HDivDivFESpace> (m, "HDivDiv")
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
		  flags_doc["discontinuous"] = "bool = False\n"
                    "  Create discontinuous HDivDiv space";
		  flags_doc["plus"] = "bool = False\n"
                    "  Add additional internal element bubble";

                  return flags_doc;
                })
    ;

  ExportFESpace<HDivDivSurfaceSpace> (m, "HDivDivSurface")    
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
		  flags_doc["discontinuous"] = "bool = False\n"
                    "  Create discontinuous HDivDivSurface space";
                  return flags_doc;
                })
    ;

  ExportFESpace<VectorFacetFESpace> (m, "VectorFacet")
    .def_static("__flags_doc__", [] ()
              {
                auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                    attr("FESpace").
                                                    attr("__flags_doc__")());
                flags_doc["highest_order_dc"] = "bool = False\n"
                  "  Splits highest order facet functions into two which are associated with\n  the corresponding neighbors and are local dofs on the corresponding element\n (used to realize projected jumps)";
                flags_doc["hide_highest_order_dc"] = "bool = False\n"
                  "  if highest_order_dc is used this flag marks the corresponding local dofs\n  as hidden dofs (reduces number of non-zero entries in a matrix). These dofs\n  can also be compressed.";
                return flags_doc;
              })
    ;

  ExportFESpace<FacetFESpace> (m, "FacetFESpace")
    .def_static("__flags_doc__", [] ()
              {
                auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                    attr("FESpace").
                                                    attr("__flags_doc__")());
                flags_doc["highest_order_dc"] = "bool = False\n"
                  "  Splits highest order facet functions into two which are associated with\n  the corresponding neighbors and are local dofs on the corresponding element\n (used to realize projected jumps)";
                flags_doc["hide_highest_order_dc"] = "bool = False\n"
                  "  if highest_order_dc is used this flag marks the corresponding local dofs\n  as hidden dofs (reduces number of non-zero entries in a matrix). These dofs\n  can also be compressed.";
                return flags_doc;
              })
    ;

  ExportFESpace<FacetSurfaceFESpace> (m, "FacetSurface");

  ExportFESpace<HDivHighOrderSurfaceFESpace> (m, "HDivSurface")
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  flags_doc["discontinuous"] = "bool = False\n"
                    "  Create discontinuous HDivSurface space";
                  return flags_doc;
                })
    .def("Average", &HDivHighOrderSurfaceFESpace::Average,
         py::arg("vector"))
    ;

  
  
  // py::class_<CompoundFESpace, shared_ptr<CompoundFESpace>, FESpace>
  //   (m, "CompoundFESpace")
  //   .def("Range", &CompoundFESpace::GetRange)
  //   ;

  py::class_<PeriodicFESpace, shared_ptr<PeriodicFESpace>, FESpace>(m, "Periodic",
	docu_string(R"delimiter(Periodic or quasi-periodic Finite Element Spaces.
The periodic fespace is a wrapper around a standard fespace with an 
additional dof mapping for the periodic degrees of freedom. All dofs 
on slave boundaries are mapped to their master dofs. Because of this, 
the mesh needs to be periodic. Low order fespaces are currently not
supported, so methods using them will not work.

Parameters:

fespace : ngsolve.comp.FESpace
    finite element space

phase : list of Complex = None
    phase shift for quasi-periodic finite element space. The basis
    functions on the slave boundary are multiplied by the factor
    given in this list. If None (default) is given, a periodic
    fespace is created. The order of the list must match the order
    of the definition of the periodic boundaries in the mesh.

used_idnrs : list of int = None
    identification numbers to be made periodic if you don't want to
    use all periodic identifications defined in the mesh, if None
    (default) all available periodic identifications are used.

)delimiter"))
    .def(py::init([] (shared_ptr<FESpace> & fes,
                      py::object phase, py::object use_idnrs )
                  {
                    Flags flags = fes->GetFlags();
                    shared_ptr<Array<int>> a_used_idnrs;
                    if(py::extract<py::list>(use_idnrs).check())
                      a_used_idnrs = make_shared<Array<int>>(makeCArray<int>(py::extract<py::list>(use_idnrs)()));
                    else
                      throw Exception("Argument for use_idnrs in Periodic must be list of identification numbers (int)");
                    shared_ptr<PeriodicFESpace> perfes;
                    auto ext = py::extract<py::list>(phase);
                    if(ext.check())
                      {
                        auto a_phase = make_shared<Array<Complex>>(py::len(ext()));
                        for (auto i : Range(a_phase->Size()))
                          {
                            auto ext_value = py::extract<Complex>(ext()[i]);
                            if(ext_value.check())
                              (*a_phase)[i] = ext_value();
                            else
                              throw Exception("Periodic FESpace needs a list of complex castable values as parameter phase");
                          }
                        perfes = make_shared<QuasiPeriodicFESpace>(fes,flags,a_used_idnrs,a_phase);
                      }
                    else if (py::isinstance<DummyArgument>(phase) || phase.is_none())
                      {
                        perfes = make_shared<PeriodicFESpace>(fes,flags,a_used_idnrs);
                      }
                    else
                      throw Exception("Periodic FESpace needs a list of complex castable values as parameter 'phase'");
                    perfes->Update(glh);
                    perfes->FinalizeUpdate(glh);
                    return perfes;
                  }), py::arg("fespace"), py::arg("phase")=DummyArgument(),
                  py::arg("use_idnrs")=py::list())
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
                    auto ret = make_shared<CompressedFESpace> (fes);
                    shared_ptr<BitArray> actdofs = nullptr;
                    if (! py::extract<DummyArgument> (active_dofs).check())
                      ret->SetActiveDofs(py::extract<shared_ptr<BitArray>>(active_dofs)());
                    ret->Update(glh);
                    return ret;                    
                  }), py::arg("fespace"), py::arg("active_dofs")=DummyArgument())
    .def("SetActiveDofs", [](CompressedFESpace & self, shared_ptr<BitArray> active_dofs)
         {
           self.SetActiveDofs(active_dofs);
         },
         py::arg("dofs"))
    ;




  /////////////////////////////// GridFunctionCoefficientFunction /////////////

  py::class_<GridFunctionCoefficientFunction, shared_ptr<GridFunctionCoefficientFunction>, CoefficientFunction>
    (m, "CoefficientFunction")
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
    ;
    

  ////////////////////////////////////// GridFunction //////////////////////////
  
  auto gf_class = py::class_<GF,shared_ptr<GF>, CoefficientFunction, NGS_Object>
    (m, "GridFunction",  "a field approximated in some finite element space", py::dynamic_attr());
  gf_class
    .def(py::init([gf_class](shared_ptr<FESpace> fes, string & name,
                                 py::kwargs kwargs)
    {
      auto flags = CreateFlagsFromKwArgs(gf_class, kwargs);
      flags.SetFlag("novisual");
      auto gf = CreateGridFunction(fes, name, flags);
      gf->Update();
      return gf;
    }), py::arg("space"), py::arg("name")="gfu",
         "creates a gridfunction in finite element space")
    .def_static("__flags_doc__", [] ()
                {
                  return py::dict
                    (
                     py::arg("multidim") = "Multidimensional GridFunction",
                     py::arg("nested") = "bool = False\n"
		     " Generates prolongation matrices for each mesh level and prolongates\n"
		     " the solution onto the finer grid after a refinement."
                     );
                })
    .def(py::pickle([] (const GridFunction& gf)
                    {
                      return py::make_tuple(gf.GetFESpace(),
                                            gf.GetName(),
                                            gf.GetFlags(),
                                            gf.GetVectorPtr());
                    },
                    [] (py::tuple state)
                    {
                      auto gf = CreateGridFunction(state[0].cast<shared_ptr<FESpace>>(),
                                                   state[1].cast<string>(),
                                                   state[2].cast<Flags>());
                      gf->Update();
                      gf->GetVector() = *py::cast<shared_ptr<BaseVector>>(state[3]);
                      return gf;
                    }
                    ))
    .def("__str__", [] (GF & self) { return ToString(self); } )
    .def_property_readonly("space", [](GF & self) { return self.GetFESpace(); },
                           "the finite element space")
    .def("Update", [](GF& self) { self.Update(); },
         "update vector size to finite element space dimension after mesh refinement")
    
    .def("Save", [](GF& self, string filename, bool parallel)
         {
           ofstream out(filename, ios::binary);
           if (parallel)
             self.Save(out);
           else
             for (auto d : self.GetVector().FVDouble())
               SaveBin(out, d);
         },
         py::arg("filename"), py::arg("parallel")=false)
    .def("Load", [](GF& self, string filename, bool parallel)
         {
           ifstream in(filename, ios::binary);
           if (parallel)
             self.Load(in);
           else
             for (auto & d : self.GetVector().FVDouble())
               LoadBin(in, d);
         },
         py::arg("filename"), py::arg("parallel")=false)         
         
    .def("Set", 
         [](shared_ptr<GF> self, spCF cf,
            VorB vb, py::object definedon)
         {
           shared_ptr<TPHighOrderFESpace> tpspace = dynamic_pointer_cast<TPHighOrderFESpace>(self->GetFESpace());          
            Region * reg = nullptr;
            if (py::extract<Region&> (definedon).check())
              reg = &py::extract<Region&>(definedon)();
            
            if(tpspace)
            {
              Transfer2TPMesh(cf.get(),self.get(),glh);
              return;
            }            
            if (reg)
              SetValues (cf, *self, *reg, NULL, glh);
            else
              SetValues (cf, *self, vb, NULL, glh);
         },
          py::arg("coefficient"),
          py::arg("VOL_or_BND")=VOL,
         py::arg("definedon")=DummyArgument(),
         "Set values"
      )
    .def_property_readonly("name", &GridFunction::GetName)

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

    .def("Deriv",
         [](shared_ptr<GF> self) -> spCF
          {
            return self->GetDeriv();
          })
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
         [](shared_ptr<GF> self, string name, VorB vb) -> py::object // shared_ptr<CoefficientFunction>
          {
            if (self->GetFESpace()->GetAdditionalEvaluators().Used(name))
              {
                auto diffop = self->GetFESpace()->GetAdditionalEvaluators()[name];
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
                return py::cast(shared_ptr<CoefficientFunction>(coef));
              }
            return py::none(); 
          }, py::arg("name"), py::arg("VOL_or_BND")=VOL)

    
    .def_property_readonly("derivname", 
                           [](shared_ptr<GF> self) -> string
                   {
                     auto deriv = self->GetFESpace()->GetFluxEvaluator();
                     if (!deriv) return "";
                     return deriv->Name();
                   })

    .def("__call__", 
         [](shared_ptr<GF> self, double x, double y, double z)
          {
            HeapReset hr(glh);
            auto space = self->GetFESpace();
            auto evaluator = space->GetEvaluator();
            IntegrationPoint ip;
            int elnr = space->GetMeshAccess()->FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            if (elnr < 0) throw Exception ("point out of domain");
            ElementId ei(VOL, elnr);
            
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
         py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0)


   .def("__call__", 
        [](shared_ptr<GF> self, const BaseMappedIntegrationPoint & mip)
          {
            HeapReset hr(glh);
            auto space = self->GetFESpace();

            ElementId ei = mip.GetTransformation().GetElementId();
            auto evaluator = space->GetEvaluator(VorB(ei));
            const FiniteElement & fel = space->GetFE(ei, glh);

            Array<int> dnums(fel.GetNDof());
            space->GetDofNrs(ei, dnums);

            if (space->IsComplex())
              {
                Vector<Complex> elvec(fel.GetNDof()*space->GetDimension());
                Vector<Complex> values(evaluator->Dim());
                self->GetElementVector(dnums, elvec);

                evaluator->Apply(fel, mip, elvec, values, glh);
                return (values.Size() > 1) ? py::cast(values) : py::cast(values(0));
              }
            else
              {
                Vector<> elvec(fel.GetNDof()*space->GetDimension());
                Vector<> values(evaluator->Dim());
                self->GetElementVector(dnums, elvec);
                evaluator->Apply(fel, mip, elvec, values, glh);
                return (values.Size() > 1) ? py::cast(values) : py::cast(values(0));
              }
          }, 
        py::arg("mip"))
    
    /*
    .def("D", 
         [](shared_ptr<GF> self, const double &x, const double &y, const double &z)
          {
            HeapReset hr(glh);
            const FESpace & space = *self->GetFESpace();
            IntegrationPoint ip;
            int dim_mesh = space.GetMeshAccess()->GetDimension();
            auto evaluator = space.GetFluxEvaluator();
            cout << evaluator->Name() << endl;
            int dim = evaluator->Dim();
            int elnr = space.GetMeshAccess()->FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            ElementId ei(VOL, elnr);
            Array<int> dnums;
            space.GetDofNrs(ei, dnums);
            const FiniteElement & fel = space.GetFE(ei, glh);
            if (space.IsComplex())
              {
                Vector<Complex> elvec;
                Vector<Complex> values(dim);
                elvec.SetSize(fel.GetNDof());
                self->GetElementVector(dnums, elvec);
                if (dim_mesh == 2)
                  {
                    MappedIntegrationPoint<2, 2> mip(ip, space.GetMeshAccess()->GetTrafo(ei, glh));
                    evaluator->Apply(fel, mip, elvec, values, glh);
                  }
                else if (dim_mesh == 3)
                  {
                    MappedIntegrationPoint<3, 3> mip(ip, space.GetMeshAccess()->GetTrafo(ei, glh));
                    evaluator->Apply(fel, mip, elvec, values, glh);
                  }
                if (dim > 1)
                  return py::cast(values);
                else
                  return py::cast(values(0));
              }
            else
              {
                Vector<> elvec;
                Vector<> values(dim);
                elvec.SetSize(fel.GetNDof());
                self->GetElementVector(dnums, elvec);
                ElementId ei(VOL, elnr);
                if (dim_mesh == 2)
                  {
                    MappedIntegrationPoint<2, 2> mip(ip, space.GetMeshAccess()->GetTrafo(ei, glh));
                    evaluator->Apply(fel, mip, elvec, values, glh);
                  }
                else if (dim_mesh == 3)
                  {
                    MappedIntegrationPoint<3, 3> mip(ip, space.GetMeshAccess()->GetTrafo(ei, glh));
                    evaluator->Apply(fel, mip, elvec, values, glh);
                  }
                if (dim > 1)
                  return py::cast(values);
                else
                  return py::cast(values(0));
              }
          },
         py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0)
    */

    .def("CF", 
         [](shared_ptr<GF> self, shared_ptr<DifferentialOperator> diffop) -> spCF
          {
            if (!diffop->Boundary())
              return make_shared<GridFunctionCoefficientFunction> (self, diffop);
            else
              return make_shared<GridFunctionCoefficientFunction> (self, nullptr, diffop);
          })
    ;



  ///////////////////////////// BilinearForm   ////////////////////////////////////////


  typedef BilinearForm BF;
  auto bf_class = py::class_<BF, shared_ptr<BilinearForm>>(m, "BilinearForm",
                                             docu_string(R"raw_string(
Used to store the left hand side of a PDE. integrators (ngsolve.BFI)
to it to implement your PDE. If the left hand side is linear
you can use BilinearForm.Assemble to assemble it after adding
your integrators. For nonlinear usage use BilinearForm.Apply or
BilinearForm.AssembleLinearization instead of Bilinearform.Assemble.

Parameters

space : ngsolve.FESpace
  The finite element space the bilinearform is defined on. This
  can be a compound FESpace for a mixed formulation.

)raw_string"));
  bf_class
    .def(py::init([bf_class] (shared_ptr<FESpace> fespace, /* bool check_unused, */
                              py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(bf_class,kwargs);
                    auto biform = CreateBilinearForm (fespace, "biform_from_py", flags);
                    // biform -> SetCheckUnused (check_unused);
                    return biform;
                  }),
         py::arg("space"))
    .def(py::init([bf_class](shared_ptr<FESpace> trial_space,
                             shared_ptr<FESpace> test_space, 
                             /* bool check_unused, */ py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(bf_class,kwargs);
                    auto biform = CreateBilinearForm (trial_space, test_space, "biform_from_py", flags);
                    // biform -> SetCheckUnused (check_unused);
                    return biform;
                  }),
         py::arg("trialspace"),
         py::arg("testspace"))
    // py::arg("check_unused")=true

    .def_static("__flags_doc__", [] ()
                {
                  return py::dict
                    (
                     py::arg("eliminate_internal") = "bool = False\n"
                     "  Set up BilinearForm for static condensation of internal\n"
                     "  bubbles. Static condensation has to be done by user,\n"
                     "  this enables only the use of the members harmonic_extension,\n"
                     "  harmonic_extension_trans and inner_solve. Have a look at the\n"
                     "  documentation for further information.",
                     py::arg("print") = "bool = False\n"
                     "  Write additional information to testout file. \n"
                     "  This file must be set by ngsolve.SetTestoutFile. Use \n"
                     "  ngsolve.SetNumThreads(1) for serial output",
                     py::arg("printelmat") = "bool = False\n"
                     "  Write element matrices to testout file",
                     py::arg("symmetric") = "bool = False\n"
                     "  If set true, only half the matrix is stored",
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
		     " The full matrix is stored, even if the symmetric flag is set.",
                     py::arg("check_unused") = "bool = True\n"
		     " If set prints warnings if not UNUSED_DOFS are not used."
                     );
                })

    .def("__str__",  []( BF & self ) { return ToString<BilinearForm>(self); } )

    .def("Add", [](BF& self, shared_ptr<BilinearFormIntegrator> bfi) -> BF&
                                 { self.AddIntegrator (bfi); return self; },
         py::return_value_policy::reference,
         "add integrator to bilinear-form")
    
    .def("__iadd__",[](BF& self, shared_ptr<BilinearFormIntegrator> other) -> BilinearForm& { self += other; return self; } )
    .def_property_readonly("space", [](BF& self) { return self.GetFESpace(); })

    .def_property_readonly("integrators", [](BF & self)
                           { return MakePyTuple (self.Integrators()); })
    
    .def("Assemble", [](BF & self, bool reallocate)
         {
           self.ReAssemble(glh,reallocate);
         }, py::call_guard<py::gil_scoped_release>(),
         py::arg("reallocate")=false)

    .def_property_readonly("mat", [](BF & self)
                                         {
                                           auto mat = self.GetMatrixPtr();
                                           if (!mat)
                                             throw py::type_error("matrix not ready - assemble bilinearform first");
                                           return mat;
                                         })

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

    .def("__call__", [](BF & self, const GridFunction & u, const GridFunction & v)
          {
            auto au = self.GetMatrix().CreateVector();
            au = self.GetMatrix() * u.GetVector();
            return InnerProduct (au, v.GetVector());
          })

    .def("Energy",[](BF & self, shared_ptr<BaseVector> x)
          {
            return self.Energy(*x, glh);
          }, py::call_guard<py::gil_scoped_release>())
    
    .def("Apply", [](BF & self, BaseVector& x, BaseVector & y)
	  {
	    self.ApplyMatrix (x, y, glh);
	  }, py::call_guard<py::gil_scoped_release>(),
         py::arg("x"),py::arg("y"), docu_string(R"raw_string(
Applies a (non-)linear variational formulation to x and stores the result in y.

Parameters

x : ngsolve.BaseVector
  input vector

y : ngsolve.BaseVector
  output vector

)raw_string"))

    .def("ComputeInternal", [](BF & self, BaseVector & u, BaseVector & f)
	  {
	    self.ComputeInternal (u, f, glh );
	  }, py::call_guard<py::gil_scoped_release>(),
         py::arg("u"),py::arg("f"))

    .def("AssembleLinearization", [](BF & self, BaseVector & ulin)
	  {
	    self.AssembleLinearization (ulin, glh);
	  }, py::call_guard<py::gil_scoped_release>(),
         py::arg("ulin"))

    .def("Flux", [](BF & self, shared_ptr<GridFunction> gf) -> spCF
          {
            return make_shared<GridFunctionCoefficientFunction> (gf, self.GetIntegrator(0));
          })
    
    .def_property_readonly("harmonic_extension", [](BF & self)
                   {
                     return self.GetHarmonicExtension();
                   }
                  )
    .def_property_readonly("harmonic_extension_trans", [](BF & self)
                   {
                     return self.GetHarmonicExtensionTrans();
                   }
                  )
    .def_property_readonly("inner_solve", [](BF & self)
                   {
                     return self.GetInnerSolve();
                   }
                  )
    .def_property_readonly("inner_matrix", [](BF & self)
                   {
                     return self.GetInnerMatrix();
                   }
                  )
    ;

  ///////////////////////////////// LinearForm //////////////////////////////////////////

  typedef LinearForm LF;
  auto lf_class = py::class_<LF, shared_ptr<LF>, NGS_Object>(m, "LinearForm", docu_string(R"raw_string(
Used to store the left hand side of a PDE. Add integrators
(ngsolve.LFI) to it to implement your PDE.

Parameters

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
                    auto flags = CreateFlagsFromKwArgs(lf_class,kwargs);
                    auto f = CreateLinearForm (fespace, "lff_from_py", flags);
                    f->AllocateVector();
                    return f;
                  }),
         py::arg("space"))
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
                           { return self->GetVectorPtr();})

    .def("Add", [](shared_ptr<LF> self, shared_ptr<LinearFormIntegrator> lfi)
          { 
            self->AddIntegrator (lfi);
            return self; 
          },
         py::arg("integrator"))
    
    .def("__iadd__",[](shared_ptr<LF> self, shared_ptr<LinearFormIntegrator> lfi)
         { (*self)+=lfi; return self; })

    .def_property_readonly("integrators", [](shared_ptr<LF> self)
                           { return MakePyTuple (self->Integrators()); })

    .def("Assemble", [](shared_ptr<LF> self)
         { self->Assemble(glh); }, py::call_guard<py::gil_scoped_release>())
    
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
          })

    ;

  /////////////////////////////// Preconditioner /////////////////////////////////////////////

  auto prec_class = py::class_<Preconditioner, shared_ptr<Preconditioner>, BaseMatrix>(m, "Preconditioner");
  prec_class
    .def(py::init([prec_class](shared_ptr<BilinearForm> bfa, const string & type, py::kwargs kwargs)
         {
           auto flags = CreateFlagsFromKwArgs(prec_class,kwargs);
           auto creator = GetPreconditionerClasses().GetPreconditioner(type);
           if (creator == nullptr)
             throw Exception(string("nothing known about preconditioner '") + type + "'");
           return creator->creatorbf(bfa, flags, "noname-pre");
         }),
         py::arg("bf"), py::arg("type"))

    .def_static("__flags_doc__", [] ()
                {
                  return py::dict
                    (
                     py::arg("inverse") = "Inverse type used in Preconditioner",
                     py::arg("test") = "bool = False\n"
                     "  Computes condition number for preconditioner, if testout file\n"
                     "  is set, prints eigenvalues to file."
                     );
                })
    .def ("Test", [](Preconditioner &pre) { pre.Test();}, py::call_guard<py::gil_scoped_release>())
    .def ("Update", [](Preconditioner &pre) { pre.Update();}, py::call_guard<py::gil_scoped_release>())
    .def_property_readonly("mat", [](Preconditioner &self)
                   {
                     return self.GetMatrixPtr();
                   })
    ;

  auto prec_multigrid = py::class_<MGPreconditioner, shared_ptr<MGPreconditioner>, Preconditioner>
    (m,"MultiGridPreconditioner");
  prec_multigrid
    .def(py::init([prec_multigrid](shared_ptr<BilinearForm> bfa, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(prec_multigrid, kwargs);
                    return make_shared<MGPreconditioner>(bfa,flags);
                  }), py::arg("bf"))
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
                  return mg_flags;
                })
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  py::class_<NumProc, NGS_Object, shared_ptr<NumProc>> (m, "NumProc")
    .def("Do", [](NumProc & self)
         {
           self.Do(glh);
         }, py::call_guard<py::gil_scoped_release>())
    ;

  py::class_<PyNumProc, NumProc, shared_ptr<PyNumProc>> (m, "PyNumProc")
    /*
    .def("__init__",
         [](NumProc *instance, shared_ptr<PDE> pde, Flags & flags)
                           {
                             new (instance) PyNumProc(pde, flags);
                           })
    */
    .def(py::init<> ([](shared_ptr<PDE> pde, Flags & flags)
                     { return new PyNumProc(pde, flags); }))
    .def_property_readonly("pde", [](NumProc &self) { return self.GetPDE(); })
    .def("Do", [](NumProc & self, LocalHeap & lh)
                               {
                                 self.Do(lh);
                               }, py::call_guard<py::gil_scoped_release>())
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
                             ngs_comm = MPI_COMM_WORLD;

                             //cout << "Rank = " << MyMPI_GetId(ngs_comm) << "/"
                             //     << MyMPI_GetNTasks(ngs_comm) << endl;

                             NGSOStream::SetGlobalActive (MyMPI_GetId()==0);
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
                                })

    .def("Add", [](shared_ptr<PDE> self, const string & name, double val)
                                {
                                  self->AddConstant (name, val);
                                })

    .def("Add", [](shared_ptr<PDE> self, shared_ptr<FESpace> space)
                                {
                                  self->AddFESpace (space->GetName(), space);
                                })

    .def("Add", [](shared_ptr<PDE> self, shared_ptr<GridFunction> gf)
                                {
                                  self->AddGridFunction (gf->GetName(), gf);
                                })

    .def("Add", [](shared_ptr<PDE> self, shared_ptr<BilinearForm> bf)
                                {
                                  self->AddBilinearForm (bf->GetName(), bf);
                                })

    .def("Add", [](shared_ptr<PDE> self, shared_ptr<LinearForm> lf)
                                {
                                  self->AddLinearForm (lf->GetName(), lf);
                                })

    .def("Add", [](shared_ptr<PDE> self, shared_ptr<Preconditioner> pre)
                                {
                                  self->AddPreconditioner (pre->GetName(), pre);
                                })

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
                                })

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
                                })

    .def("SetCurveIntegrator", [](shared_ptr<PDE> self, const string & filename, shared_ptr<LinearFormIntegrator> lfi)
          {
            self->SetLineIntegratorCurvePointInfo(filename, lfi.get());
          })

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
           shared_ptr<MeshAccess> ma, 
           VorB vb, int order, py::object definedon,
	   bool region_wise, bool element_wise)
        {
          static Timer t("Integrate CF"); RegionTimer reg(t);
          // static mutex addcomplex_mutex;
          BitArray mask;
          {
            py::gil_scoped_acquire aquire;
            py::extract<Region> defon_region(definedon);
            if (defon_region.check())
              {
                vb = VorB(defon_region());
                mask = BitArray(defon_region().Mask());
              }
          }
          if(!mask.Size()){
            mask = BitArray(ma->GetNRegions(vb));
            mask.Set();
          }
          int dim = cf->Dimension();
          if((region_wise || element_wise) && dim != 1)
            throw Exception("region_wise and element_wise only implemented for 1 dimensional coefficientfunctions");
          
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
                     AsAtomic(sum(i)) += hsum(i);
                   if(region_wise)
                     AsAtomic(region_sum(el.GetIndex())) += hsum(0);
                   if (element_wise)
                     element_sum(el.Nr()) = hsum(0);
                 });
              py::object result;
              if (region_wise) {
#ifdef PARALLEL
                Vector<> rs2(ma->GetNRegions(vb));
                MPI_Allreduce(&region_sum(0), &rs2(0), ma->GetNRegions(vb), MPI_DOUBLE, MPI_SUM, ngs_comm);
                region_sum = rs2;
#endif
                result = py::list(py::cast(region_sum));
              }
              else if (element_wise)
                result = py::cast(element_sum);
              else if(dim==1) {
                sum(0) = MyMPI_AllReduce(sum(0));
                result = py::cast(sum(0));
              }
              else {
#ifdef PARALLEL
                Vector<> gsum(dim);
                MPI_Allreduce(&sum(0), &gsum(0), dim, MPI_DOUBLE, MPI_SUM, ngs_comm);
                sum = gsum;
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
                     MyAtomicAdd (sum(i), hsum(i));
                   if(region_wise)
                     MyAtomicAdd (region_sum(el.GetIndex()), hsum(0));
                   if (element_wise)
                     element_sum(el.Nr()) = hsum(0);
                 });
              
              py::object result;
              if (region_wise) {
#ifdef PARALLEL
                Vector<Complex> rs2(ma->GetNRegions(vb));
                MPI_Allreduce(&region_sum(0), &rs2(0), ma->GetNRegions(vb), MPI_Traits<Complex>::MPIType(), MPI_SUM, ngs_comm);
                region_sum = rs2;
#endif
                result = py::list(py::cast(region_sum));
              }
              else if (element_wise)
                result = py::cast(element_sum);
              else if(dim==1) {
                sum(0) = MyMPI_AllReduce(sum(0));
                result = py::cast(sum(0));
              }
              else {
#ifdef PARALLEL
                Vector<Complex> gsum(dim);
                MPI_Allreduce(&sum(0), &gsum(0), dim, MPI_Traits<Complex>::MPIType(), MPI_SUM, ngs_comm);
                sum = gsum;
#endif
                result = py::cast(sum);
              }
              return result;
            }
        },
	py::arg("cf"), py::arg("mesh"), py::arg("VOL_or_BND")=VOL, 
	py::arg("order")=5,
	py::arg("definedon")=DummyArgument(),
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
  
  m.def("SymbolicLFI",
          [](spCF cf, VorB vb, bool element_boundary,
             bool skeleton, py::object definedon,
             IntegrationRule ir, int bonus_intorder, py::object definedonelem,
             bool simd_evaluate, VorB element_vb) 
           {
             py::extract<Region> defon_region(definedon);
             if (defon_region.check())
               vb = VorB(defon_region());

             if (element_boundary) element_vb = BND;
             
             shared_ptr<LinearFormIntegrator> lfi;
             if (!skeleton)
               lfi = make_shared<SymbolicLinearFormIntegrator> (cf, vb, element_vb);
             else
               lfi = make_shared<SymbolicFacetLinearFormIntegrator> (cf, vb /* , element_boundary */);
             
             if (py::extract<py::list> (definedon).check())
               {
                 Array<int> defon = makeCArray<int> (definedon);
                 for (int & d : defon) d--;
                 lfi -> SetDefinedOn (defon); 
               }

             lfi->SetSimdEvaluate (simd_evaluate);
             // lfi -> SetDefinedOn (makeCArray<int> (definedon));

             if (defon_region.check())
               lfi->SetDefinedOn(defon_region().Mask());
             lfi -> SetBonusIntegrationOrder(bonus_intorder);
	     if (ir.Size())
               {
                 cout << IM(1) << "WARNING: Setting the integration rule for all element types is deprecated, use LFI.SetIntegrationRule(ELEMENT_TYPE, IntegrationRule) instead!" << endl;
                 dynamic_pointer_cast<SymbolicLinearFormIntegrator>
		   (lfi)->SetIntegrationRule(ir);                   
               }

             if (! py::extract<DummyArgument> (definedonelem).check())
               lfi -> SetDefinedOnElements (py::extract<shared_ptr<BitArray>>(definedonelem)());

             return shared_ptr<LinearFormIntegrator>(lfi);
           },
           py::arg("form"),
           py::arg("VOL_or_BND")=VOL,
           py::arg("element_boundary")=false,
           py::arg("skeleton")=false,           
           py::arg("definedon")=DummyArgument(),
	   py::arg("intrule")=IntegrationRule(),
           py::arg("bonus_intorder")=0,
           py::arg("definedonelements")=DummyArgument(),
           py::arg("simd_evaluate")=true,
           py::arg("element_vb")=VOL        
          );

  m.def("SymbolicBFI",
          [](spCF cf, VorB vb, bool element_boundary,
             bool skeleton, py::object definedon,
             IntegrationRule ir, int bonus_intorder, py::object definedonelem,
             bool simd_evaluate, VorB element_vb)
           {
             py::extract<Region> defon_region(definedon);
             if (defon_region.check())
               vb = VorB(defon_region());

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
             
             if (py::extract<py::list> (definedon).check())
               {
                 Array<int> defon = makeCArray<int> (definedon);
                 for (int & d : defon) d--;
                 bfi -> SetDefinedOn (defon); 
               }
             // bfi -> SetDefinedOn (makeCArray<int> (definedon));
             bfi -> SetBonusIntegrationOrder(bonus_intorder);
             if (defon_region.check())
               bfi->SetDefinedOn(defon_region().Mask());

             if (ir.Size())
               {
                 cout << IM(1) << "WARNING: Setting the integration rule for all element types is deprecated, use BFI.SetIntegrationRule(ELEMENT_TYPE, IntegrationRule) instead!" << endl;
                 dynamic_pointer_cast<SymbolicBilinearFormIntegrator> (bfi)
                   ->SetIntegrationRule(ir);
               }

             bfi->SetSimdEvaluate (simd_evaluate);
             
             if (! py::extract<DummyArgument> (definedonelem).check())
               bfi -> SetDefinedOnElements (py::extract<shared_ptr<BitArray>>(definedonelem)());
             return shared_ptr<BilinearFormIntegrator>(bfi);
           },
        py::arg("form"), py::arg("VOL_or_BND")=VOL,
        py::arg("element_boundary")=false,
        py::arg("skeleton")=false,
        py::arg("definedon")=DummyArgument(),
        py::arg("intrule")=IntegrationRule(),
        py::arg("bonus_intorder")=0,
        py::arg("definedonelements")=DummyArgument(),
        py::arg("simd_evaluate")=true,
        py::arg("element_vb")=VOL
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
           VorB element_vb)
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
             return bfi;
           },
        py::arg("form"), py::arg("VOL_or_BND")=VOL,
        py::arg("definedon")=DummyArgument(), py::arg("element_boundary")=false,
        py::arg("bonus_intorder")=0,
        py::arg("definedonelements")=DummyArgument(),
        py::arg("simd_evaluate")=true,
        py::arg("element_vb")=VOL        
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
              auto spaces = makeCArraySharedPtr<shared_ptr<FESpace>> (spaces_list);
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
                      py::list names_list, string filename, int subdivision, int only_element)
         -> shared_ptr<BaseVTKOutput>
         {
           Array<shared_ptr<CoefficientFunction> > coefs
             = makeCArraySharedPtr<shared_ptr<CoefficientFunction>> (coefs_list);
           Array<string > names
             = makeCArray<string> (names_list);
           shared_ptr<BaseVTKOutput> ret;
           if (ma->GetDimension() == 2)
             ret = make_shared<VTKOutput<2>> (ma, coefs, names, filename, subdivision, only_element);
           else
             ret = make_shared<VTKOutput<3>> (ma, coefs, names, filename, subdivision, only_element);
           return ret;
         }),
         py::arg("ma"),
         py::arg("coefs")= py::list(),
         py::arg("names") = py::list(),
         py::arg("filename") = "vtkout",
         py::arg("subdivision") = 0,
         py::arg("only_element") = -1
         )
     .def("Do", [](shared_ptr<BaseVTKOutput> self)
          { 
            self->Do(glh);
          },
          py::call_guard<py::gil_scoped_release>())
     .def("Do", [](shared_ptr<BaseVTKOutput> self, const BitArray * drawelems)
          { 
            self->Do(glh, drawelems);
          },
          py::arg("drawelems"),
          py::call_guard<py::gil_scoped_release>())
     ;

  /////////////////////////////////////////////////////////////////////////////////////
}




PYBIND11_MODULE(libngcomp, m) {
  m.attr("__name__") = "comp";
  ExportNgcomp(m);
}

#endif // NGS_PYTHON
