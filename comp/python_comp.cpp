#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <boost/python/slice.hpp>
#include <boost/python/iterator.hpp>
#include <comp.hpp>

#ifdef PARALLEL
#include <mpi4py/mpi4py.h>
#endif

#include <regex>

using namespace ngcomp;

using ngfem::ELEMENT_TYPE;

typedef GridFunction GF;
typedef PyWrapperDerived<GF, CoefficientFunction> PyGF;
typedef PyWrapper<FESpace> PyFES;

template <typename T>
struct PythonTupleFromFlatArray {
  static PyObject* convert(FlatArray<T> ar)
    {
      bp::list res;
      for(int i = 0; i < ar.Size(); i++) 
        res.append (ar[i]);
      bp::tuple tup(res);
      return bp::incref(tup.ptr());
    }
};

template <typename T>
struct PythonTupleFromArray {
  static PyObject* convert(const Array<T> & ar)
    {
      bp::list res;
      for(int i = 0; i < ar.Size(); i++) 
        res.append (ar[i]);
      bp::tuple tup(res);
      return bp::incref(tup.ptr());
    }
};


template <typename T> void PyExportArray ()
{
  boost::python::to_python_converter< FlatArray<T>, PythonTupleFromFlatArray<T> >();
  boost::python::to_python_converter< Array<T>, PythonTupleFromArray<T> >();
}




template <> class cl_NonElement<ElementId>
{
public:
  static ElementId Val() { return ElementId(VOL,-1); }
};


class PyNumProc : public NumProc
{

public:
  PyNumProc (shared_ptr<PDE> pde, const Flags & flags) : NumProc (pde, flags) { ; }
  shared_ptr<PDE> GetPDE() const { return shared_ptr<PDE> (pde); }
  // virtual void Do (LocalHeap & lh) { cout << "should not be called" << endl; }
};

class NumProcWrap : public PyNumProc, public bp::wrapper<PyNumProc> {
public:
  NumProcWrap (shared_ptr<PDE> pde, const Flags & flags) : PyNumProc(pde, flags) { ; }
  virtual void Do(LocalHeap & lh)  {
    // cout << "numproc wrap - do" << endl;
    AcquireGIL gil_lock;
    try
      {
        this->get_override("Do")(boost::ref(lh));
      }
    catch (bp::error_already_set const &) {
      cout << "caught a python error:" << endl;
      PyErr_Print();
    }
  }
};

typedef PyWrapperDerived<ProxyFunction, CoefficientFunction> PyProxyFunction;

MSVC2015_UPDATE3_GET_PTR_FIX(NumProcWrap)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::BaseVTKOutput)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::BilinearForm)
MSVC2015_UPDATE3_GET_PTR_FIX(PyWrapper<ngcomp::BilinearForm>)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::CompoundFESpace)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::FESpace)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::GridFunction)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::HCurlHighOrderFESpace)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::LinearForm)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::MeshAccess)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::NGS_Object)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::NumProc)
MSVC2015_UPDATE3_GET_PTR_FIX(ngcomp::Preconditioner)

bp::object MakeProxyFunction2 (const FESpace & fes,
                              bool testfunction,
                              const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock)
{
  auto compspace = dynamic_cast<const CompoundFESpace*> (&fes);
  if (compspace)
    {
      bp::list l;
      int nspace = compspace->GetNSpaces();
      for (int i = 0; i < nspace; i++)
        {
          l.append (MakeProxyFunction2 ( *(*compspace)[i], testfunction,
                                         [&] (shared_ptr<ProxyFunction> proxy)
                                         {
                                           auto block_eval = make_shared<CompoundDifferentialOperator> (proxy->Evaluator(), i);
                                           shared_ptr <CompoundDifferentialOperator> block_deriv_eval = nullptr;
                                           if (proxy->DerivEvaluator() != nullptr)
                                             block_deriv_eval = make_shared<CompoundDifferentialOperator> (proxy->DerivEvaluator(), i);
                                           shared_ptr <CompoundDifferentialOperator> block_trace_eval = nullptr;
                                           if (proxy->TraceEvaluator() != nullptr)
                                             block_trace_eval = make_shared<CompoundDifferentialOperator> (proxy->TraceEvaluator(), i);
                                           shared_ptr <CompoundDifferentialOperator> block_trace_deriv_eval = nullptr;
                                           if (proxy->TraceDerivEvaluator() != nullptr)
                                             block_trace_deriv_eval = make_shared<CompoundDifferentialOperator> (proxy->TraceDerivEvaluator(), i);
                                           auto block_proxy = make_shared<ProxyFunction> (/* &fes, */ testfunction, fes.IsComplex(),                                                                                          block_eval, block_deriv_eval, block_trace_eval, block_trace_deriv_eval);

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

  /*
  shared_ptr<CoefficientFunction> proxy =
    addblock(make_shared<ProxyFunction> (testfunction, fes.IsComplex(),
                                         fes.GetEvaluator(),
                                         fes.GetFluxEvaluator(),
                                         fes.GetEvaluator(BND),
                                         fes.GetFluxEvaluator(BND)
                                         ));
  */
  auto proxy = make_shared<ProxyFunction>  (testfunction, fes.IsComplex(),
                                            fes.GetEvaluator(),
                                            fes.GetFluxEvaluator(),
                                            fes.GetEvaluator(BND),
                                            fes.GetFluxEvaluator(BND));
  auto add_diffops = fes.GetAdditionalEvaluators();
  for (int i = 0; i < add_diffops.Size(); i++)
    proxy->SetAdditionalEvaluator (add_diffops.GetName(i), add_diffops[i]);

  proxy = addblock(proxy);
  return bp::object(PyProxyFunction(proxy));
}

bp::object MakeProxyFunction (const FESpace & fes,
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


void NGS_DLL_HEADER ExportNgcomp()
{
  bp::docstring_options local_docstring_options(true, true, false);
  
  std::string nested_name = "comp";
  if( bp::scope() )
    nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".comp");
  bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
  
  cout << IM(1) << "exporting comp as " << nested_name << endl;
  bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
  parent.attr("comp") = module;
  
  bp::scope local_scope(module);

  //////////////////////////////////////////////////////////////////////////////////////////

  bp::enum_<VorB>("VorB")
    .value("VOL", VOL)
    .value("BND", BND)
    .export_values()
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  bp::enum_<COUPLING_TYPE> ("COUPLING_TYPE")
    .value("UNUSED_DOF", UNUSED_DOF)
    .value("LOCAL_DOF", LOCAL_DOF)
    .value("INTERFACE_DOF", INTERFACE_DOF)
    .value("NONWIREBASKET_DOF", NONWIREBASKET_DOF)
    .value("WIREBASKET_DOF", WIREBASKET_DOF)
    .value("EXTERNAL_DOF", EXTERNAL_DOF)
    .value("ANY_DOF", ANY_DOF)
    // .export_values()
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  bp::class_<ElementId> ("ElementId", 
                         "an element identifier containing element number and Volume/Boundary flag",
                         bp::no_init)
    .def(bp::init<VorB,int>())
    .def(bp::init<int>())
    .def("__str__", &ToString<ElementId>)
    .add_property("nr", &ElementId::Nr, "the element number")    
    .def("IsVolume", &ElementId::IsVolume, "is it a boundary element ?")
    .def("IsBoundary", &ElementId::IsBoundary, "is it a volume element ?")
    .def(bp::self!=bp::self)
    .def("__eq__" , FunctionPointer( [](ElementId &self, ElementId &other)
                                    { return !(self!=other); }) )
    .def("__hash__" , &ElementId::Nr)
    ;
  
  bp::def("BndElementId", FunctionPointer([] (int nr) { return ElementId(BND,nr); }),
          (bp::arg("nr")),
          "creates an element-id for a boundary element")
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  bp::class_<ElementRange,bp::bases<IntRange>> ("ElementRange",bp::init<const MeshAccess&,VorB,IntRange>())
    .def(PyDefIterable2<ElementRange>())
    ;

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<FESpace::ElementRange>);
  bp::class_<FESpace::ElementRange,shared_ptr<FESpace::ElementRange>, bp::bases<IntRange>, boost::noncopyable> ("FESpaceElementRange",bp::no_init)
    // .def(bp::init<const FESpace::ElementRange&>())
    // .def(bp::init<FESpace::ElementRange&&>())
    .def(PyDefIterable3<FESpace::ElementRange>())
    ;


  //////////////////////////////////////////////////////////////////////////////////////////

  bp::class_<Ngs_Element,bp::bases<ElementId>>("Ngs_Element", bp::no_init)
    .add_property("vertices", FunctionPointer([](Ngs_Element &el) {return bp::tuple(Array<int>(el.Vertices()));} ),
                  "list of global vertex numbers")
    .add_property("edges", FunctionPointer([](Ngs_Element &el) { return bp::tuple(Array<int>(el.Edges()));} ),
                  "list of global edge numbers")
    .add_property("faces", FunctionPointer([](Ngs_Element &el) { return bp::tuple(Array<int>(el.Faces()));} ),
                  "list of global face numbers")
    .add_property("type", &Ngs_Element::GetType, "geometric shape of element")
    .add_property("index", &Ngs_Element::GetIndex, "material or boundary condition index")
    .add_property("mat", FunctionPointer([](Ngs_Element & el)
                                         { return el.GetMaterial() ? *el.GetMaterial() : ""; }),
                  "material or boundary condition label")
    ;

  bp::class_<FESpace::Element,bp::bases<Ngs_Element>>("FESpaceElement", bp::no_init)
    .add_property("dofs",
                  FunctionPointer
                  ([](FESpace::Element & el) 
                   {
                     Array<int> tmp (el.GetDofs());
                     return bp::tuple(tmp);
                     // return bp::tuple(Array<int>(el.GetDofs()));} ))
                   }),
                  "degrees of freedom of element"
                  )

    .def("GetLH", FunctionPointer([](FESpace::Element & el) -> LocalHeap & 
                                  {
                                    return el.GetLH();
                                  }),
         bp::return_value_policy<bp::reference_existing_object>()
         )
    
    .def("GetFE", FunctionPointer([](FESpace::Element & el) -> const FiniteElement & 
                                  {
                                    return el.GetFE();
                                  }),
         bp::return_value_policy<bp::reference_existing_object>(),
         "the finite element containing shape functions"
         )

    .def("GetTrafo", FunctionPointer([](FESpace::Element & el) -> const ElementTransformation & 
                                     {
                                       return el.GetTrafo();
                                     }),
         bp::return_value_policy<bp::reference_existing_object>(),
         "the transformation from reference element to physical element"
         )

    ;
  //////////////////////////////////////////////////////////////////////////////////////////


  bp::class_<GlobalDummyVariables> ("GlobalVariables", bp::no_init)
    .add_property("msg_level", 
                 &GlobalDummyVariables::GetMsgLevel,
                 &GlobalDummyVariables::SetMsgLevel)
    .add_property("testout", 
                 &GlobalDummyVariables::GetTestoutFile,
                 &GlobalDummyVariables::SetTestoutFile)
    /*
    .add_property("pajetrace",
		  &GlobalDummyVariables::GetTestoutFile,
		  FunctionPointer([] (GlobalDummyVariables&, bool use)
				  { TaskManager::SetPajeTrace(use); }));
    */
    .add_property("pajetrace",
		  &GlobalDummyVariables::GetTestoutFile,
		  FunctionPointer([] (GlobalDummyVariables&, size_t size)
				  {
                                    TaskManager::SetPajeTrace(size > 0);
                                    PajeTrace::SetMaxTracefileSize(size);
                                  }))
    .add_property("numthreads",
		  FunctionPointer([] (GlobalDummyVariables&)
				  {
                                    return TaskManager::GetMaxThreads ();
                                  }),
		  FunctionPointer([] (GlobalDummyVariables&, int numthreads)
				  {
                                    TaskManager::SetNumThreads (numthreads);
                                  }))
    // &GlobalDummyVariables::SetTestoutFile)
    ;

  bp::scope().attr("ngsglobals") = bp::object(bp::ptr(&globvar));

  //////////////////////////////////////////////////////////////////////////////////

  PyExportArray<string>();

  struct MeshAccess_pickle_suite : bp::pickle_suite
  {
    static
    bp::tuple getinitargs(const MeshAccess & ma)
    {
      return bp::make_tuple(); 
    }

    static
    bp::tuple getstate(bp::object o)
    {
      auto & ma = bp::extract<MeshAccess const&>(o)();
      stringstream str;
      ma.SaveMesh(str);
      return bp::make_tuple (o.attr("__dict__"), str.str());
    }
    
    static
    void setstate(bp::object o, bp::tuple state)
    {
      auto & ma = bp::extract<MeshAccess&>(o)();

      /*
      if (len(state) != 2)
        {
          PyErr_SetObject(PyExc_ValueError,
                          ("expected 2-item tuple in call to __setstate__; got %s"
                           % state).ptr()
                          );
          throw_error_already_set();
        }
      */

      bp::dict d = bp::extract<bp::dict>(o.attr("__dict__"))();
      d.update(state[0]);
      string s = bp::extract<string>(state[1]);
      stringstream str(s);
      ma.LoadMesh (str);
    }

    static bool getstate_manages_dict() { return true; }
  };


  struct FESpace_pickle_suite : bp::pickle_suite
  {
    static
    bp::tuple getinitargs(bp::object obj)
    {
      auto fes = bp::extract<PyFES>(obj)().Get();
      bp::object m (fes->GetMeshAccess());
      bp::object flags = obj.attr("__dict__")["flags"];
      flags["dim"] = fes->GetDimension();
      return bp::make_tuple(fes->type, m, flags, fes->GetOrder(), fes->IsComplex());
    }

    static
    bp::tuple getstate(bp::object o)
    {
      // auto & fes = bp::extract<FESpace const&>(o)();
      return bp::make_tuple (o.attr("__dict__")); // , str.str());
    }
    
    static
    void setstate(bp::object o, bp::tuple state)
    {
      // auto & fes = bp::extract<FESpace&>(o)();
      bp::dict d = bp::extract<bp::dict>(o.attr("__dict__"))();
      d.update(state[0]);
    }

    static bool getstate_manages_dict() { return true; }
  };
  


  //////////////////////////////////////////////////////////////////////////////////////////

  bp::class_<Region> ("Region", "a subset of volume or boundary elements", bp::no_init)
    .def(bp::init<shared_ptr<MeshAccess>,VorB,string>())
    .def("Mask", FunctionPointer([](Region & reg)->BitArray { return reg.Mask(); }))
    .def(bp::self + bp::self)
    .def(bp::self + string())
    .def(bp::self - bp::self)
    .def(bp::self - string())
    .def(~bp::self)
    ;

  bp::implicitly_convertible <Region, BitArray> ();


  //////////////////////////////////////////////////////////////////////////////////////////
  
  
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<MeshAccess>);
  bp::class_<MeshAccess, shared_ptr<MeshAccess>>("Mesh", 
                                                 "the mesh")
    .def(bp::init<shared_ptr<netgen::Mesh>>())
    .def_pickle(MeshAccess_pickle_suite())
    .def("__ngsid__", FunctionPointer( [] ( MeshAccess & self)
        { return reinterpret_cast<std::uintptr_t>(&self); } ) )
    
#ifndef PARALLEL
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](const string & filename)
                           { 
                             return make_shared<MeshAccess> (filename);
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("filename"))
          ))

#else

    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](const string & filename,
                              bp::object py_mpicomm)
                           { 
                             PyObject * py_mpicomm_ptr = py_mpicomm.ptr();
                             if (py_mpicomm_ptr != Py_None)
                               {
                                 MPI_Comm * comm = PyMPIComm_Get (py_mpicomm_ptr);
                                 ngs_comm = *comm;
                               }
                             else
                               ngs_comm = MPI_COMM_WORLD;

                             NGSOStream::SetGlobalActive (MyMPI_GetId()==0);
                             return make_shared<MeshAccess> (filename, ngs_comm);
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("filename"), bp::arg("mpicomm")=bp::object())
          ))
#endif

    
    .def("__eq__", FunctionPointer
         ( [] (shared_ptr<MeshAccess> self, shared_ptr<MeshAccess> other)
           {
             return self == other;
           }))

    .def("LoadMesh", static_cast<void(MeshAccess::*)(const string &)>(&MeshAccess::LoadMesh),
         "Load mesh from file")
    
    .def("Elements", static_cast<ElementRange(MeshAccess::*)(VorB)const> (&MeshAccess::Elements),
         (bp::arg("VOL_or_BND")=VOL))

    .def("__getitem__", static_cast<Ngs_Element(MeshAccess::*)(ElementId)const> (&MeshAccess::operator[]))

    .def ("GetNE", static_cast<int(MeshAccess::*)(VorB)const> (&MeshAccess::GetNE))
    .add_property ("nv", &MeshAccess::GetNV, "number of vertices")
    .add_property ("ne",  static_cast<int(MeshAccess::*)()const> (&MeshAccess::GetNE), "number of volume elements")
    .add_property ("dim", &MeshAccess::GetDimension, "mesh dimension")
    .add_property ("ngmesh", &MeshAccess::GetNetgenMesh, "netgen mesh")
    .def ("GetTrafo", 
          static_cast<ElementTransformation&(MeshAccess::*)(ElementId,Allocator&)const>
          (&MeshAccess::GetTrafo), 
          bp::return_value_policy<bp::reference_existing_object>())

    .def ("GetTrafo", FunctionPointer([](MeshAccess & ma, ElementId id)
                                      {
                                        return &ma.GetTrafo(id, global_alloc);
                                      }),
          bp::return_value_policy<bp::manage_new_object>())

    // .def("SetDeformation", &MeshAccess::SetDeformation)
    .def("SetDeformation", FunctionPointer
	 ([](MeshAccess & ma, PyGF gf)
          { ma.SetDeformation(gf.Get()); }))
    .def("SetRadialPML", &MeshAccess::SetRadialPML)
    .def("UnsetDeformation", FunctionPointer
	 ([](MeshAccess & ma){ ma.SetDeformation(nullptr);}))
    
    .def("GetMaterials", FunctionPointer
	 ([](const MeshAccess & ma)
	  {
	    Array<string> materials(ma.GetNDomains());
	    for (int i : materials.Range())
	      materials[i] = ma.GetDomainMaterial(i);
	    return bp::list(materials);
	  }),
         (bp::arg("self")),
         "returns list of materials"
         )

    .def("Materials", FunctionPointer
	 ([](shared_ptr<MeshAccess> ma, string pattern) 
	  {
            return new Region (ma, VOL, pattern);
	  }),
         (bp::arg("self"), bp::arg("pattern")),
         "returns mesh-region matching the given regex pattern",
         bp::return_value_policy<bp::manage_new_object>()
         )
    
    .def("GetBoundaries", FunctionPointer
	 ([](const MeshAccess & ma)
	  {
	    Array<string> materials(ma.GetNBoundaries());
	    for (int i : materials.Range())
	      materials[i] = ma.GetBCNumBCName(i);
	    return bp::list(materials);
	  }),
         (bp::arg("self")),
         "returns list of boundary conditions"
         )

    .def("Boundaries", FunctionPointer
	 ([](shared_ptr<MeshAccess> ma, string pattern)
	  {
            return new Region (ma, BND, pattern);
	  }),
         (bp::arg("self"), bp::arg("pattern")),
         "returns boundary mesh-region matching the given regex pattern",
         bp::return_value_policy<bp::manage_new_object>()
         )

    .def("Refine", FunctionPointer
         ([](MeshAccess & ma)
          {
            ma.Refine();
          }),
         "local mesh refinement based on marked elements, uses element-bisection algorithm")

    .def("RefineHP", FunctionPointer
         ([](MeshAccess & ma, int levels, double factor)
          {
            Ng_HPRefinement(levels, factor);
            ma.UpdateBuffers();
          }),
         (bp::arg("self"), bp::arg("levels"), bp::arg("factor")=0.125),
         "geometric mesh refinement towards marked vertices and edges, uses factor for placement of new points"
         )

    .def("SetRefinementFlag", &MeshAccess::SetRefinementFlag,
         "set refinementflag for mesh-refinement")

    .def("GetParentElement", &MeshAccess::GetParentElement)
    .def("GetParentVertices", FunctionPointer
         ([](MeshAccess & ma, int vnum)
          {
            Array<int> parents(2);
            ma.GetParentNodes (vnum, &parents[0]);
            return bp::tuple(parents);
          }))
    
    .def("Curve", FunctionPointer
         ([](MeshAccess & ma, int order)
          {
            Ng_HighOrder(order);
          }),
         (bp::arg("self"),bp::arg("order")))



    .def("__call__", FunctionPointer
         ([](MeshAccess & ma, double x, double y, double z, VorB vb) 
          {
            IntegrationPoint ip;
            int elnr;
            if (vb == VOL)
              elnr = ma.FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            else
              elnr = ma.FindSurfaceElementOfPoint(Vec<3>(x, y, z), ip, true);
              
            if (elnr < 0) throw Exception ("point out of domain");

            ElementTransformation & trafo = ma.GetTrafo(ElementId(vb, elnr), global_alloc);
            BaseMappedIntegrationPoint & mip = trafo(ip, global_alloc);
            mip.SetOwnsTrafo(true);
            return &mip;
          } 
          ), 
         (bp::arg("self"), bp::arg("x") = 0.0, bp::arg("y") = 0.0, bp::arg("z") = 0.0,
          bp::arg("VOL_or_BND") = VOL),
         bp::return_value_policy<bp::manage_new_object>()
         )

    .def("Contains", FunctionPointer
         ([](MeshAccess & ma, double x, double y, double z) 
          {
            IntegrationPoint ip;
            int elnr = ma.FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            return (elnr >= 0);
          }), 
         (bp::arg("self"), bp::arg("x") = 0.0, bp::arg("y") = 0.0, bp::arg("z") = 0.0)
         )

    ;

  //////////////////////////////////////////////////////////////////////////////////////////
  
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<NGS_Object>);
  bp::class_<NGS_Object, shared_ptr<NGS_Object>, boost::noncopyable>("NGS_Object", bp::no_init)
    .add_property("name", FunctionPointer
                  ([](const NGS_Object & self)->string { return self.GetName();}))
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  typedef PyWrapper<CoefficientFunction> PyCF;

  bp::class_<PyProxyFunction, bp::bases<PyCF> > ("ProxyFunction",
                         // bp::init<FESpace*,bool, shared_ptr<DifferentialOperator>, shared_ptr<DifferentialOperator>>()
                         bp::no_init)
    .def("Deriv", FunctionPointer
         ([](const PyProxyFunction self)
          { return PyProxyFunction(self->Deriv()); }),
         "take canonical derivative (grad, curl, div)")
    .def("Trace", FunctionPointer
         ([](const PyProxyFunction self)
          { return PyProxyFunction(self->Trace()); }),
         "take canonical boundary trace")
    .def("Other", FunctionPointer
         ([](const PyProxyFunction self, bp::object bnd) 
          {
            if (bp::extract<double> (bnd).check())
              return PyProxyFunction(self->Other(make_shared<ConstantCoefficientFunction>(bp::extract<double> (bnd)())));
            if (bp::extract<PyCF> (bnd).check())
              return PyProxyFunction(self->Other(bp::extract<PyCF> (bnd)().Get()));
            else
              return PyProxyFunction(self->Other(nullptr));
          }),
         "take value from neighbour element (DG)",
          (bp::arg("self"), bp::arg("bnd") = bp::object())
         )
    .add_property("derivname", FunctionPointer
                  ([](const PyProxyFunction self) -> string
                   {
                     if (!self->Deriv()) return "";
                     return self->DerivEvaluator()->Name();
                   }))
    .def("Operator", FunctionPointer
         ([] (const PyProxyFunction self, string name) -> bp::object 
          {
            auto op = self->GetAdditionalProxy(name);
            if (op)
              return bp::object(PyProxyFunction(op));
            return bp::object();
          }))
    ;




  struct OrderProxy 
  {
    FESpace & fes;
    OrderProxy (FESpace & afes) : fes(afes) { ; }
  };


  bp::class_<OrderProxy> ("OrderProxy", bp::no_init)
    .def("__setitem__", FunctionPointer
         ([] (OrderProxy & self, ElementId ei, int o) 
          {
            cout << "set order of el " << ei << " to order " << o << endl;
            cout << "(not implemented)" << endl;
          }))

    .def("__setitem__", FunctionPointer
         ([] (OrderProxy & self, ELEMENT_TYPE et, int o) 
          {
            cout << "set order of eltype " << et << " to order " << o << endl;
            self.fes.SetBonusOrder (et, o - self.fes.GetOrder());

            LocalHeap lh (100000, "FESpace::Update-heap", true);
            self.fes.Update(lh);
            self.fes.FinalizeUpdate(lh);
          }))
    
    .def("__setitem__", FunctionPointer
         ([] (OrderProxy & self, NODE_TYPE nt, int o) 
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
          }))

    .def("__setitem__", FunctionPointer
         ([] (OrderProxy & self, NODE_TYPE nt, int nr, int o) 
          {
            cout << "set order of " << nt << " " << nr << " to " << o << endl;
            cout << "(not implemented)" << endl;
          }))

    .def("__setitem__", FunctionPointer
         ([] (OrderProxy & self, bp::tuple tup, int o) 
          {
            NODE_TYPE nt = bp::extract<NODE_TYPE>(tup[0])();
            int nr = bp::extract<int>(tup[1])();
            cout << "set order of " << nt << " " << nr << " to " << o << endl;
            cout << "(not implemented)" << endl;
          }))
    

    

    /*
    .def("__setitem__", FunctionPointer([] (OrderProxy & self, bp::slice inds, int o) 
                                        {
                                          cout << "set order to slice, o = " <<o << endl;
                                          auto ndof = self.fes.GetNDof();
                                          bp::object indices = inds.attr("indices")(ndof);
                                          int start = bp::extract<int> (indices[0]);
                                          int stop = bp::extract<int> (indices[1]);
                                          int step = bp::extract<int> (indices[2]);
                                          cout << "start = " << start << ", stop = " << stop << ", step = " << step << endl;
                                        }))

    .def("__setitem__", FunctionPointer([] (OrderProxy & self, bp::list inds, int o) 
                                        {
                                          cout << "set order list" << endl;

                                          for (int i = 0; i < len(inds); i++)
                                            cout << bp::extract<int> (inds[i]) << endl;
                                        }))
    */

    /*
    .def("__setitem__", FunctionPointer([] (OrderProxy & self, bp::object generator, int o) 
                                        {
                                          cout << "general setitem called" << endl;

                                          if (bp::extract<int> (generator).check())
                                            {
                                              cout << " set order, int" << endl;
                                              return;
                                            }

                                          if (bp::extract<ElementId> (generator).check())
                                            {
                                              cout << " set order, elid" << endl;
                                              return;
                                            }
                                          if (bp::extract<bp::slice> (generator).check())
                                            {
                                              cout << " set order, slice" << endl;
                                              return;
                                            }
                                          
                                          cout << "set order from generator" << endl;
                                          try
                                            {
                                              auto iter = generator.attr("__iter__")();
                                              while (1)
                                                {
                                                  auto el = iter.attr("__next__")();
                                                  cout << bp::extract<int> (el) << " ";
                                                }
                                            }
                                          catch (bp::error_already_set&) 
                                            { 
                                              if (PyErr_ExceptionMatches (PyExc_StopIteration))
                                                {
                                                  cout << endl;
                                                  PyErr_Clear();
                                                }
                                              else
                                                {
                                                  cout << "some other error" << endl;
                                                }
                                            };
                                        }))
    */
    ;


  static size_t global_heapsize = 1000000;
  static LocalHeap glh(global_heapsize, "python-comp lh", true);
  bp::def("SetHeapSize", FunctionPointer([](size_t heapsize)
                                         {
                                           if (heapsize > global_heapsize)
                                             {
                                               global_heapsize = heapsize;
                                               glh = LocalHeap (heapsize, "python-comp lh", true);
                                             }
                                         }));

  bp::def("SetTestoutFile", FunctionPointer([](string filename)
                                            {
                                              testout = new ofstream (filename);
                                            }));
  
  //////////////////////////////////////////////////////////////////////////////////////////
  bp::class_<PyFES>("FESpace",  "a finite element space", bp::no_init)


    .def("__dummy_init__", bp::make_constructor 
         (FunctionPointer ([](shared_ptr<MeshAccess> ma, const string & type, 
                              bp::dict bpflags, int order, bool is_complex,
                              bp::object dirichlet, bp::object definedon, int dim)
                           {
                             Flags flags = bp::extract<Flags> (bpflags)();

                             if (order > -1) flags.SetFlag ("order", order);
                             if (dim > -1) flags.SetFlag ("dim", dim);
                             if (is_complex) flags.SetFlag ("complex");

                             bp::extract<bp::list> dirlist(dirichlet);
                             if (dirlist.check())
                               flags.SetFlag("dirichlet", makeCArray<double>(dirlist()));

                             bp::extract<string> dirstring(dirichlet);
                             if (dirstring.check())
                               {
                                 std::regex pattern(dirstring());
                                 Array<double> dirlist;
                                 for (int i = 0; i < ma->GetNBoundaries(); i++)
                                   if (std::regex_match (ma->GetBCNumBCName(i), pattern))
                                     dirlist.Append (i+1);
                                 flags.SetFlag("dirichlet", dirlist);
                               }

                             bp::extract<string> definedon_string(definedon);
                             if (definedon_string.check())
                               {
                                 regex definedon_pattern(definedon_string());
                                 Array<double> defonlist;
                                 for (int i = 0; i < ma->GetNDomains(); i++)
                                   if (regex_match(ma->GetDomainMaterial(i), definedon_pattern))
                                     defonlist.Append(i+1);
                                 flags.SetFlag ("definedon", defonlist);
                               }
                             bp::extract<bp::list> definedon_list(definedon);
                             if (definedon_list.check())
                               flags.SetFlag ("definedon", makeCArray<double> (definedon));
                             bp::extract<Region> definedon_reg(definedon);
                             if (definedon_reg.check() && definedon_reg().IsVolume())
                               {
                                 Array<double> defonlist;
                                 for (int i = 0; i < definedon_reg().Mask().Size(); i++)
                                   if (definedon_reg().Mask().Test(i))
                                     defonlist.Append(i+1);
                                 flags.SetFlag ("definedon", defonlist);
                               }
                             
                             
                             auto fes = CreateFESpace (type, ma, flags); 
                             LocalHeap lh (1000000, "FESpace::Update-heap");
                             fes->Update(lh);
                             fes->FinalizeUpdate(lh);
                             return new PyFES(fes);
                             })
          ))

    // the raw - constructor
    .def("__init__", 
         FunctionPointer ([](bp::object self, const string & type, bp::object bp_ma, 
                             bp::dict bp_flags, int order, bool is_complex,
                             bp::object dirichlet, bp::object definedon, int dim)
                          {
                            shared_ptr<MeshAccess> ma = bp::extract<shared_ptr<MeshAccess>>(bp_ma)();
                            auto ret = self.attr("__dummy_init__")(ma, type, bp_flags, order, is_complex, dirichlet, definedon, dim);
                            self.attr("__dict__")["flags"] = bp_flags;
                            return ret;   
                           }),
         bp::default_call_policies(),        // need it to use arguments
         (bp::arg("type"), bp::arg("mesh"), bp::arg("flags") = bp::dict(), 
           bp::arg("order")=-1, 
           bp::arg("complex")=false, 
           bp::arg("dirichlet")= bp::object(),
           bp::arg("definedon")=bp::object(),
          bp::arg("dim")=-1 ),
         "allowed types are: 'h1ho', 'l2ho', 'hcurlho', 'hdivho' etc."
         )
    

    
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](bp::list lspaces, bp::dict bpflags)
                           {
                             Flags flags = bp::extract<Flags> (bpflags)();

                             auto spaces = makeCArrayUnpackWrapper<PyWrapper<FESpace>> (lspaces);
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
                             
                             shared_ptr<FESpace> fes = make_shared<CompoundFESpace> (spaces[0]->GetMeshAccess(), spaces, flags);
                             LocalHeap lh (1000000, "FESpace::Update-heap");
                             fes->Update(lh);
                             fes->FinalizeUpdate(lh);
                             return new PyFES(fes);
                           }),
          bp::default_call_policies(),       
          (bp::arg("spaces"), bp::arg("flags") = bp::dict())),
         "construct compound-FESpace from list of component spaces"
         )
    .def_pickle(FESpace_pickle_suite())
    .def("__ngsid__", FunctionPointer( [] ( PyFES & self)
        { return reinterpret_cast<std::uintptr_t>(self.Get().get()); } ) )
    .def("Update", FunctionPointer([](PyFES & self, int heapsize)
                                   { 
                                     LocalHeap lh (heapsize, "FESpace::Update-heap");
                                     self->Update(lh);
                                     self->FinalizeUpdate(lh);
                                   }),
         (bp::arg("self"),bp::arg("heapsize")=1000000),
         "update space after mesh-refinement")

    .add_property ("ndof", FunctionPointer([](PyFES & self) { return self->GetNDof(); }), 
                   "number of degrees of freedom")

    .add_property ("ndofglobal", FunctionPointer([](PyFES & self) { return self->GetNDofGlobal(); }), 
                   "global number of dofs on MPI-distributed mesh")
    .def("__str__", &ToString<FESpace>)

    // .add_property("mesh", FunctionPointer ([](FESpace & self) -> shared_ptr<MeshAccess>
    // { return self.GetMeshAccess(); }))

    .add_property("order", FunctionPointer([] (PyFES & self) { return OrderProxy(*self.Get()); }),
                  "proxy to set order for individual nodes")
    .add_property("globalorder", FunctionPointer([] (PyFES & self) { return self->GetOrder(); }),
                  "query global order of space")    
    .add_property("type", FunctionPointer([] (PyFES & self) { return self->type; }),
                  "type of finite element space")    

    .def("Elements", 
         FunctionPointer([](PyFES & self, VorB vb, int heapsize) 
                         {
                           return make_shared<FESpace::ElementRange> (self->Elements(vb, heapsize));
                         }),
         (bp::arg("self"),bp::arg("VOL_or_BND")=VOL,bp::arg("heapsize")=10000))

    .def("Elements", 
         FunctionPointer([](PyFES & self, VorB vb, LocalHeap & lh) 
                         {
                           return make_shared<FESpace::ElementRange> (self->Elements(vb, lh));
                         }),
         (bp::arg("self"), bp::arg("VOL_or_BND")=VOL, bp::arg("heap")))

    /*
    .def("Elements", 
         FunctionPointer([](FESpace & self, VorB vb, LocalHeap & lh, int heapsize) 
                         {
                           cout << "lh.avail = " << lh.Available() << endl;
                           return make_shared<FESpace::ElementRange> (self.Elements(vb, heapsize));
                         }),
         (bp::arg("self"),bp::arg("VOL_or_BND")=VOL, 
          bp::arg("heap")=LocalHeap(0), bp::arg("heapsize")=10000))
    */

    .def("GetDofNrs", FunctionPointer([](PyFES & self, ElementId ei) 
                                   {
                                     Array<int> tmp; self->GetDofNrs(ei,tmp); 
                                     return bp::tuple (tmp); 
                                   }))

    .def("CouplingType", FunctionPointer ([](PyFES & self, int dofnr) -> COUPLING_TYPE
                                          { return self.Get()->GetDofCouplingType(dofnr); }),
         // (bp::arg("self"),bp::arg("dofnr")),
         "get coupling type of a degree of freedom"
         )
    .def("SetCouplingType", FunctionPointer ([](PyFES & self, int dofnr, COUPLING_TYPE ct) 
                                             { return self.Get()->SetDofCouplingType(dofnr,ct); }),
         (bp::arg("self"),bp::arg("dofnr")),
         "set coupling type of a degree of freedom"
         )

    /*
    .def ("GetFE", 
          static_cast<FiniteElement&(FESpace::*)(ElementId,Allocator&)const>
          (&FESpace::GetFE), 
          bp::return_value_policy<bp::reference_existing_object>())
    */
    .def ("GetFE", FunctionPointer([](PyFES & self, ElementId ei) -> bp::object
                                   {
                                     Allocator alloc;

                                     auto fe = shared_ptr<FiniteElement> (&self->GetFE(ei, alloc));

                                     auto scalfe = dynamic_pointer_cast<BaseScalarFiniteElement> (fe);
                                     if (scalfe) return bp::object(scalfe);

                                     return bp::object(fe);

                                   }))
    
    .def ("GetFE", FunctionPointer([](PyFES & self, ElementId ei, LocalHeap & lh)
                                   {
                                     return &self->GetFE(ei, lh);
                                   }),
          bp::return_value_policy<bp::reference_existing_object>())


    .def("FreeDofs", FunctionPointer
         ( [] (const PyFES &self, bool coupling) -> const BitArray &{ return *self->GetFreeDofs(coupling); } ),
         bp::return_value_policy<bp::reference_existing_object>(),
         (bp::arg("self"), 
          bp::arg("coupling")=false))

    .def("Range", FunctionPointer
         ( [] (const PyFES & self, int comp) -> bp::slice
           {
             auto compspace = dynamic_pointer_cast<CompoundFESpace> (self.Get());
             if (!compspace)
               bp::exec("'Range' is available only for product spaces");
             IntRange r = compspace->GetRange(comp);
             return bp::slice(r.First(), r.Next());
           }))

    .add_property("components", FunctionPointer
                  ([](PyFES & self)-> bp::tuple
                   { 
                     auto compspace = dynamic_pointer_cast<CompoundFESpace> (self.Get());
                     if (!compspace)
                       bp::exec("'components' is available only for product spaces");
                     bp::list vecs;
                     for (int i = 0; i < compspace -> GetNSpaces(); i++) 
                       vecs.append( PyFES((*compspace)[i]) );
                     return bp::tuple(vecs);
                   }),
                  "list of gridfunctions for compound gridfunction")

    .def("TrialFunction", FunctionPointer
         ( [] (const PyFES & self) 
           {
             return MakeProxyFunction (*self.Get(), false);
           }),
         (bp::args("self")))
    .def("TestFunction", FunctionPointer
         ( [] (const PyFES & self) 
           {
             return MakeProxyFunction (*self.Get(), true);
           }),
         (bp::args("self")))

    .def("SolveM", FunctionPointer
        ( [] (const PyFES & self,
              PyCF rho, shared_ptr<BaseVector> vec, int heapsize)
          {
            // LocalHeap lh(heapsize, "solveM - lh", true);
            if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
              }
            self->SolveM(*rho.Get(), *vec, glh);
          }),
        (bp::args("self"), 
         bp::args("rho"), bp::args("vec"), bp::args("heapsize")=1000000))
        
    .def("__eq__", FunctionPointer
         ( [] (PyFES self, PyFES other)
           {
             return self.Get() == other.Get();
           }))
    ;
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<HCurlHighOrderFESpace>);
  bp::class_<HCurlHighOrderFESpace, shared_ptr<HCurlHighOrderFESpace>, bp::bases<PyFES>,boost::noncopyable>
    ("HCurlFunctionsWrap",bp::no_init)
    .def("CreateGradient", FunctionPointer([](PyFES & self) {
	  auto hcurl = dynamic_pointer_cast<HCurlHighOrderFESpace>(self.Get());
	  auto fesh1 = hcurl->CreateGradientSpace();
	  shared_ptr<BaseMatrix> grad = hcurl->CreateGradient(*fesh1);
	  auto fes = new PyFES(fesh1);
	  return bp::make_tuple(grad, fes);
	}))
    ;
  
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<CompoundFESpace>);
  bp::class_<CompoundFESpace, shared_ptr<CompoundFESpace>, bp::bases<PyFES>, boost::noncopyable>
    ("CompoundFESpace", bp::no_init)
    .def("Range", &CompoundFESpace::GetRange)
    ;

  //////////////////////////////////////////////////////////////////////////////////////////
  

  struct GF_pickle_suite : bp::pickle_suite
  {
    static
    bp::tuple getinitargs(bp::object obj)
    {
      auto gf = bp::extract<PyGF>(obj)().Get();
      bp::object space = obj.attr("__dict__")["space"];
      return bp::make_tuple(space, gf->GetName());
    }

    static
    bp::tuple getstate(bp::object obj)
    {
      auto gf = bp::extract<PyGF>(obj)().Get();
      bp::object bp_vec(gf->GetVectorPtr());
      return bp::make_tuple (obj.attr("__dict__"), bp_vec);
    }
    
    static
    void setstate(bp::object obj, bp::tuple state)
    {
      auto gf = bp::extract<PyGF>(obj)().Get();
      bp::dict d = bp::extract<bp::dict>(obj.attr("__dict__"))();
      d.update(state[0]);
      gf->GetVector() = *bp::extract<shared_ptr<BaseVector>> (state[1])();
    }

    static bool getstate_manages_dict() { return true; }
  };
  


  
  bp::class_<PyGF, bp::bases<PyCF>>
    ("GridFunction",  "a field approximated in some finite element space", bp::no_init)


    
    // raw - constructor
    .def("__init__",
         FunctionPointer ([](bp::object self, bp::object bp_fespace, string name, bp::object multidim)
                          {
                            auto fespace = bp::extract<PyFES>(bp_fespace)();

                            auto ret = 
                              bp::make_constructor
                              (FunctionPointer ([](PyFES fespace, string name, bp::object multidim)
                              {
                                Flags flags;
                                flags.SetFlag ("novisual");
                                if (bp::extract<int>(multidim).check())
                                  flags.SetFlag ("multidim", bp::extract<int>(multidim)());
                                auto gf = CreateGridFunction (fespace.Get(), name, flags);
                                gf->Update();
                                return new PyGF(gf);
                              }))(self, fespace, name, multidim);
                            
                            self.attr("__dict__")["space"] = bp_fespace;
                            return ret;   
                          }),
         (bp::arg("space"), bp::arg("name")="gfu", bp::arg("multidim")=bp::object()),
         "creates a gridfunction in finite element space"
         )
    .def_pickle(GF_pickle_suite())    
    .def("__ngsid__", FunctionPointer( [] ( PyGF self)
        { return reinterpret_cast<std::uintptr_t>(self.Get().get()); } ) )
    .def("__str__", &ToString<GF>)
    .add_property("space", FunctionPointer([](bp::object self) -> bp::object
                                           {
                                             bp::dict d = bp::extract<bp::dict>(self.attr("__dict__"))();
                                             // if gridfunction is created from python, it has the space attribute
                                             if (d.has_key("space"))
                                               return d.get("space");

                                             // if not, make a new python space object from C++ space
                                             return bp::object(PyFES(bp::extract<PyGF>(self)()->GetFESpace()));
                                           }),
                  "the finite element space")
    // .add_property ("space", &GF::GetFESpace, "the finite element spaces")
    .def("Update", FunctionPointer ([](PyGF self) { self->Update(); }),
         "update vector size to finite element space dimension after mesh refinement")
    
    .def("Save", FunctionPointer([](PyGF self, string filename)
                                 {
                                   ofstream out(filename, ios::binary);
                                   self->Save(out);
                                 }))
    .def("Load", FunctionPointer([](PyGF self, string filename)
                                 {
                                   ifstream in(filename, ios::binary);
                                   self->Load(in);
                                 }))
    
    .def("Set", FunctionPointer
         ([](PyGF self, PyCF cf,
             bool boundary, bp::object definedon, int heapsize, bp::object heap)
          {
            Region * reg = nullptr;
            if (bp::extract<Region&> (definedon).check())
              reg = &bp::extract<Region&>(definedon)();
            
            if (bp::extract<LocalHeap&> (heap).check())
              {
                LocalHeap & lh = bp::extract<LocalHeap&> (heap)();
                if (reg)
                  SetValues (cf.Get(), *self.Get(), *reg, NULL, lh);
                else
                  SetValues (cf.Get(), *self.Get(), boundary, NULL, lh);
                return;
              }

            if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
              }
            // LocalHeap lh(heapsize, "GridFunction::Set-lh", true);
            if (reg)
              SetValues (cf.Get(), *self.Get(), *reg, NULL, glh);
            else
              SetValues (cf.Get(), *self.Get(), boundary, NULL, glh);
          }),
          bp::default_call_policies(),        // need it to use arguments
         (bp::arg("self"),bp::arg("coefficient"),
          bp::arg("boundary")=false,
          bp::arg("definedon")=bp::object(),
          bp::arg("heapsize")=1000000, bp::arg("heap")=bp::object()),
         "Set values"
      )


    .add_property("components", FunctionPointer
                  ([](PyGF self)-> bp::tuple
                   { 
                     bp::list vecs;
                     for (int i = 0; i < self->GetNComponents(); i++)
                       vecs.append(PyGF(self->GetComponent(i)));
                     return bp::tuple(vecs);
                   }),
                  "list of gridfunctions for compound gridfunction")

    .add_property("vec",
                  FunctionPointer([](PyGF self) { return self->GetVectorPtr(); }),
                  "coefficient vector")

    .add_property("vecs", FunctionPointer
                  ([](PyGF self)-> bp::list
                   { 
                     bp::list vecs;
                     for (int i = 0; i < self->GetMultiDim(); i++)
                       vecs.append(self->GetVectorPtr(i));
                     return vecs;
                   }),
                  "list of coefficient vectors for multi-dim gridfunction")

    /*
    .def("CF", FunctionPointer
         ([](shared_ptr<GF> self) -> shared_ptr<CoefficientFunction>
          {
            return make_shared<GridFunctionCoefficientFunction> (self);
          }))

    .def("CF", FunctionPointer
         ([](shared_ptr<GF> self, shared_ptr<DifferentialOperator> diffop)
          -> shared_ptr<CoefficientFunction>
          {
            return make_shared<GridFunctionCoefficientFunction> (self, diffop);
          }))
    */
    .def("Deriv", FunctionPointer
         ([](PyGF self) -> PyCF
          {
            auto sp = make_shared<GridFunctionCoefficientFunction> (self.Get(),
                                                                    self->GetFESpace()->GetFluxEvaluator(),
                                                                    self->GetFESpace()->GetFluxEvaluator(BND));
            // sp->SetComplex(self->GetFESpace()->IsComplex()); 
            sp->SetDimensions(sp->Dimensions());
            return PyCF (sp);
          }))

    .def("Operator", FunctionPointer
         ([](PyGF self, string name) -> bp::object // shared_ptr<CoefficientFunction>
          {
            if (self->GetFESpace()->GetAdditionalEvaluators().Used(name))
              {
                auto diffop = self->GetFESpace()->GetAdditionalEvaluators()[name];
                cout << "diffop is " << typeid(*diffop).name() << endl;
                PyCF coef(make_shared<GridFunctionCoefficientFunction> (self.Get(), diffop));
                return bp::object(coef);
              }
            return bp::object(); //  shared_ptr<CoefficientFunction>();
          }))

    
    .add_property("derivname", FunctionPointer
                  ([](PyGF self) -> string
                   {
                     auto deriv = self->GetFESpace()->GetFluxEvaluator();
                     if (!deriv) return "";
                     return deriv->Name();
                   }))

    .def("__call__", FunctionPointer
         ([](PyGF self, double x, double y, double z)
          {
            auto space = self->GetFESpace();
            auto evaluator = space->GetEvaluator();
            LocalHeap lh(10000, "ngcomp::GridFunction::Eval");

            IntegrationPoint ip;
            int elnr = space->GetMeshAccess()->FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            if (elnr < 0) throw Exception ("point out of domain");

            const FiniteElement & fel = space->GetFE(elnr, lh);

            Array<int> dnums(fel.GetNDof(), lh);
            space->GetDofNrs(elnr, dnums);
            auto & trafo = space->GetMeshAccess()->GetTrafo(elnr, false, lh);

            if (space->IsComplex())
              {
                Vector<Complex> elvec(fel.GetNDof()*space->GetDimension());
                Vector<Complex> values(evaluator->Dim());
                self->GetElementVector(dnums, elvec);

                evaluator->Apply(fel, trafo(ip, lh), elvec, values, lh);
                return (values.Size() > 1) ? bp::object(values) : bp::object(values(0));
              }
            else
              {
                Vector<> elvec(fel.GetNDof()*space->GetDimension());
                Vector<> values(evaluator->Dim());
                self->GetElementVector(dnums, elvec);

                evaluator->Apply(fel, trafo(ip, lh), elvec, values, lh);
                return (values.Size() > 1) ? bp::object(values) : bp::object(values(0));
              }
          }
          ), (bp::arg("self"), bp::arg("x") = 0.0, bp::arg("y") = 0.0, bp::arg("z") = 0.0))


   .def("__call__", FunctionPointer
        ([](PyGF self, const BaseMappedIntegrationPoint & mip)
          {
            auto space = self->GetFESpace();

            ElementId ei = mip.GetTransformation().GetElementId();
            // auto evaluator = space->GetEvaluator(ei.IsBoundary());
            auto evaluator = space->GetEvaluator(VorB(ei));
            LocalHeap lh(10000, "ngcomp::GridFunction::Eval");

            // int elnr = mip.GetTransformation().GetElementNr();
            const FiniteElement & fel = space->GetFE(ei, lh);

            Array<int> dnums(fel.GetNDof());
            space->GetDofNrs(ei, dnums);

            if (space->IsComplex())
              {
                Vector<Complex> elvec(fel.GetNDof()*space->GetDimension());
                Vector<Complex> values(evaluator->Dim());
                self->GetElementVector(dnums, elvec);

                evaluator->Apply(fel, mip, elvec, values, lh);
                return (values.Size() > 1) ? bp::object(values) : bp::object(values(0));
              }
            else
              {
                Vector<> elvec(fel.GetNDof()*space->GetDimension());
                Vector<> values(evaluator->Dim());
                self->GetElementVector(dnums, elvec);
                evaluator->Apply(fel, mip, elvec, values, lh);
                return (values.Size() > 1) ? bp::object(values) : bp::object(values(0));
              }
          }), 
        (bp::arg("self"), bp::arg("mip")))
    

    .def("D", FunctionPointer
         ([](PyGF self, const double &x, const double &y, const double &z)
          {
            const FESpace & space = *self->GetFESpace();
            IntegrationPoint ip;
            int dim_mesh = space.GetMeshAccess()->GetDimension();
            auto evaluator = space.GetFluxEvaluator();
            cout << evaluator->Name() << endl;
            int dim = evaluator->Dim();
            LocalHeap lh(10000, "ngcomp::GridFunction::Eval");
            int elnr = space.GetMeshAccess()->FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            Array<int> dnums;
            space.GetDofNrs(elnr, dnums);
            const FiniteElement & fel = space.GetFE(elnr, lh);
            if (space.IsComplex())
              {
                Vector<Complex> elvec;
                Vector<Complex> values(dim);
                elvec.SetSize(fel.GetNDof());
                self->GetElementVector(dnums, elvec);
                if (dim_mesh == 2)
                  {
                    MappedIntegrationPoint<2, 2> mip(ip, space.GetMeshAccess()->GetTrafo(elnr, false, lh));
                    evaluator->Apply(fel, mip, elvec, values, lh);
                  }
                else if (dim_mesh == 3)
                  {
                    MappedIntegrationPoint<3, 3> mip(ip, space.GetMeshAccess()->GetTrafo(elnr, false, lh));
                    evaluator->Apply(fel, mip, elvec, values, lh);
                  }
                if (dim > 1)
                  return bp::object(values);
                else
                  return bp::object(values(0));
              }
            else
              {
                Vector<> elvec;
                Vector<> values(dim);
                elvec.SetSize(fel.GetNDof());
                self->GetElementVector(dnums, elvec);
                if (dim_mesh == 2)
                  {
                    MappedIntegrationPoint<2, 2> mip(ip, space.GetMeshAccess()->GetTrafo(elnr, false, lh));
                    evaluator->Apply(fel, mip, elvec, values, lh);
                  }
                else if (dim_mesh == 3)
                  {
                    MappedIntegrationPoint<3, 3> mip(ip, space.GetMeshAccess()->GetTrafo(elnr, false, lh));
                    evaluator->Apply(fel, mip, elvec, values, lh);
                  }
                if (dim > 1)
                  return bp::object(values);
                else
                  return bp::object(values(0));
              }
          }
          ), (bp::arg("self"), bp::arg("x") = 0.0, bp::arg("y") = 0.0, bp::arg("z") = 0.0))


    .def("CF", FunctionPointer
         ([](PyGF self, shared_ptr<DifferentialOperator> diffop) -> PyCF
          {
            if (!diffop->Boundary())
              return PyCF(make_shared<GridFunctionCoefficientFunction> (self.Get(), diffop));
            else
              return PyCF(make_shared<GridFunctionCoefficientFunction> (self.Get(), nullptr, diffop));
          }))
    
    ;



  //////////////////////////////////////////////////////////////////////////////////////////

  PyExportArray<shared_ptr<BilinearFormIntegrator>> ();

  // typedef BilinearForm BF;
  typedef PyWrapper<BilinearForm> PyBF;
  bp::class_<PyBF>("BilinearForm", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](PyFES fespace, string name, 
                              bool symmetric, bp::dict bpflags)
                           {
                             Flags flags = bp::extract<Flags> (bpflags)();
                             if (symmetric) flags.SetFlag("symmetric");
                             return new PyBF(CreateBilinearForm (fespace.Get(), name, flags));
                           }),
          bp::default_call_policies(),        // need it to use arguments
          (bp::arg("space"),
           bp::arg("name")="bfa", 
           bp::arg("symmetric") = false,
           bp::arg("flags") = bp::dict())))
    
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](PyFES fespace, PyFES fespace2,
                              string name, bool symmetric, bp::dict bpflags)
                           {
                             Flags flags = bp::extract<Flags> (bpflags)();
                             if (symmetric) flags.SetFlag("symmetric");
                             return new PyBF(CreateBilinearForm (fespace.Get(), fespace2.Get(), name, flags));
                           }),
          bp::default_call_policies(),        // need it to use arguments
          (bp::arg("space"),
           bp::arg("space2"),
           bp::arg("name")="bfa", 
           bp::arg("symmetric") = false,
           bp::arg("flags") = bp::dict())))

    .def("__str__", FunctionPointer( []( PyBF & self ) { return ToString<BilinearForm>(*self.Get()); } ))

    .def("Add", FunctionPointer ([](PyBF & self, PyWrapper<BilinearFormIntegrator> bfi) -> PyBF&
                                 { self->AddIntegrator (bfi.Get()); return self; }),
         bp::return_value_policy<bp::reference_existing_object>(),
         "add integrator to bilinear-form")
    
    .def("__iadd__",FunctionPointer
                  ([](PyBF self, PyWrapper<BilinearFormIntegrator> & other) { *self.Get()+=other.Get(); return self; } ))

    .add_property("integrators", FunctionPointer
                  ([](PyBF & self)
                   {
                     bp::list igts;
                     for (auto igt : self->Integrators())
                       igts.append (PyWrapper<BilinearFormIntegrator> (igt));
                     return igts;
                     // return bp::object (self->Integrators());
                   } ))
    
    .def("Assemble", FunctionPointer([](PyBF & self, int heapsize, bool reallocate)
                                     {
                                       // LocalHeap lh (heapsize, "BilinearForm::Assemble-heap", true);
                                       if (heapsize > global_heapsize)
                                         {
                                           global_heapsize = heapsize;
                                           glh = LocalHeap(heapsize, "python-comp lh", true);
                                         }
                                       self->ReAssemble(glh,reallocate);
                                     }),
         (bp::arg("self")=NULL,bp::arg("heapsize")=1000000,bp::arg("reallocate")=false))

    // .add_property("mat", static_cast<shared_ptr<BaseMatrix>(BilinearForm::*)()const> (&BilinearForm::GetMatrixPtr))
    .add_property("mat", FunctionPointer([](PyBF & self)
                                         {
                                           auto mat = self->GetMatrixPtr();
                                           if (!mat)
                                             bp::exec("raise RuntimeError('matrix not ready - assemble bilinearform first')\n");
                                           return mat;
                                         }))

    .def("__getitem__", FunctionPointer( [](PyBF & self, bp::tuple t)
                                         {
                                           int ind1 = bp::extract<int>(t[0])();
                                           int ind2 = bp::extract<int>(t[1])();
                                           cout << "get bf, ind = " << ind1 << "," << ind2 << endl;
                                         }))
    

    .add_property("components", FunctionPointer
                  ([](PyBF & self)-> bp::list
                   { 
                     bp::list bfs;
                     auto fes = dynamic_pointer_cast<CompoundFESpace> (self->GetFESpace());
                     if (!fes)
                       bp::exec("raise RuntimeError('not a compound-fespace')\n");
                       
                     int ncomp = fes->GetNSpaces();
                     for (int i = 0; i < ncomp; i++)
                       // bfs.append(shared_ptr<BilinearForm> (new ComponentBilinearForm(self.Get().get(), i, ncomp)));
                       bfs.append(PyWrapper<BilinearForm> (make_shared<ComponentBilinearForm>(self.Get(), i, ncomp)));
                     return bfs;
                   }),
                  "list of components for bilinearforms on compound-space")

    .def("__call__", FunctionPointer
         ([](PyBF & self, const GridFunction & u, const GridFunction & v)
          {
            auto au = self->GetMatrix().CreateVector();
            au = self->GetMatrix() * u.GetVector();
            return InnerProduct (au, v.GetVector());
          }))

    // .def("Energy", &BilinearForm::Energy)
    .def("Energy",FunctionPointer
         ([](PyBF & self, BaseVector & x)
          {
            return self->Energy(x);
          }))
    .def("Apply", FunctionPointer
	 ([](PyBF & self, BaseVector & x, BaseVector & y, int heapsize)
	  {
	    // static LocalHeap lh (heapsize, "BilinearForm::Apply", true);
            if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
              }
	    self->ApplyMatrix (x, y, glh);
	  }),
         (bp::arg("self")=NULL,bp::arg("x"),bp::arg("y"),bp::arg("heapsize")=1000000))
    .def("ComputeInternal", FunctionPointer
	 ([](PyBF & self, BaseVector & u, BaseVector & f, int heapsize)
	  {
	    LocalHeap lh (heapsize, "BilinearForm::ComputeInternal");
	    self->ComputeInternal (u ,f ,lh );
	  }),
         (bp::arg("self")=NULL,bp::arg("u"),bp::arg("f"),bp::arg("heapsize")=1000000))

    .def("AssembleLinearization", FunctionPointer
	 ([](PyBF & self, BaseVector & ulin, int heapsize)
	  {
	    // LocalHeap lh (heapsize, "BilinearForm::Assemble-heap", true);
            if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
              }
	    self->AssembleLinearization (ulin, glh);
	  }),
         (bp::arg("self")=NULL,bp::arg("ulin"),bp::arg("heapsize")=1000000))

    .def("Flux", FunctionPointer
         ([](PyBF & self, shared_ptr<GridFunction> gf) -> PyCF
          {
            return PyCF(make_shared<GridFunctionCoefficientFunction> (gf, self->GetIntegrator(0)));
          }))
    .add_property("harmonic_extension", FunctionPointer
                  ([](PyBF & self)
                   {
                     return shared_ptr<BaseMatrix> (&self->GetHarmonicExtension(),
                                                    &NOOP_Deleter);
                   })
                  )
    .add_property("harmonic_extension_trans", FunctionPointer
                  ([](PyBF & self)
                   {
                     return shared_ptr<BaseMatrix> (&self->GetHarmonicExtensionTrans(),
                                                    &NOOP_Deleter);
                   })
                  )
    .add_property("inner_solve", FunctionPointer
                  ([](PyBF & self)
                   {
                     return shared_ptr<BaseMatrix> (&self->GetInnerSolve(),
                                                    &NOOP_Deleter);
                   })
                  )
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  PyExportArray<shared_ptr<LinearFormIntegrator>> ();

  // typedef LinearForm LF;
  typedef PyWrapper<LinearForm> PyLF;
  bp::class_<PyLF>("LinearForm", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](PyFES fespace, string name, Flags flags) // -> shared_ptr<LinearForm>
                           {
                             auto f = CreateLinearForm (fespace.Get(), name, flags);
                             f->AllocateVector();
                             return new PyLF(f);
                           }),
          bp::default_call_policies(),        // need it to use arguments
          (bp::arg("space"), bp::arg("name")="lff", bp::arg("flags") = bp::dict()))
         )
    .def("__str__", FunctionPointer( []( PyLF & self ) { return ToString<LinearForm>(*self.Get()); } ))

    .add_property("vec", FunctionPointer([] (PyLF self) { return self->GetVectorPtr();}))

    .def("Add", FunctionPointer
         ([](PyLF self, PyWrapper<LinearFormIntegrator> lfi)
          { 
            self->AddIntegrator (lfi.Get());
            return self; 
          }),
         bp::default_call_policies(),        // need it to use arguments
         (bp::arg("self"), bp::arg("integrator")))

    .def("__iadd__",FunctionPointer
                  ([](PyLF self, PyWrapper<LinearFormIntegrator> & other) { *self.Get()+=other.Get(); return self; } ))


    .add_property("integrators", FunctionPointer
                  ([](PyLF self)
                   {
                     bp::list igts;
                     for (auto igt : self->Integrators())
                       igts.append (PyWrapper<LinearFormIntegrator> (igt));
                     return igts;
                     // return bp::object (self->Integrators());
                   } ))

    .def("Assemble", FunctionPointer
         ([](PyLF self, int heapsize)
          { 
            // LocalHeap lh(heapsize, "LinearForm::Assemble-heap", true);
            if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
              }
            self->Assemble(glh);
          }),
         (bp::arg("self")=NULL,bp::arg("heapsize")=1000000))

    .add_property("components", FunctionPointer
                  ([](PyLF self)-> bp::list
                   { 
                     bp::list lfs;
                     auto fes = dynamic_pointer_cast<CompoundFESpace> (self->GetFESpace());
                     if (!fes)
                       bp::exec("raise RuntimeError('not a compound-fespace')\n");
                       
                     int ncomp = fes->GetNSpaces();
                     for (int i = 0; i < ncomp; i++)
                       // lfs.append(shared_ptr<LinearForm> (new ComponentLinearForm(self.Get().get(), i, ncomp)));
                       lfs.append(PyWrapper<LinearForm> (make_shared<ComponentLinearForm>(self.Get(), i, ncomp)));
                     return lfs;
                   }),
                  "list of components for linearforms on compound-space")
    
    .def("__call__", FunctionPointer
         ([](PyLF & self, const GridFunction & v)
          {
            return InnerProduct (self->GetVector(), v.GetVector());
          }))

    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  typedef Preconditioner PRE;
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<PRE>);
  bp::class_<PRE, shared_ptr<PRE>, boost::noncopyable>("Preconditioner", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](PyWrapper<BilinearForm> bfa, const string & type,
                              Flags flags)
                           { 
                             auto creator = GetPreconditionerClasses().GetPreconditioner(type);
                             if (creator == nullptr)
                               throw Exception(string("nothing known about preconditioner '") + type + "'");
                             return creator->creatorbf(bfa.Get(), flags, "noname-pre");
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("bf"), bp::arg("type"), bp::arg("flags")=bp::dict())
          ))

    .def ("Update", &Preconditioner::Update)
    .add_property("mat", FunctionPointer
                  ([](shared_ptr<Preconditioner> self) 
                   {
                     return shared_ptr<BaseMatrix> (const_cast<BaseMatrix*> (&self->GetMatrix()),
                                                    NOOP_Deleter);
                   }))
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<NumProc>);
  bp::class_<NumProc, shared_ptr<NumProc>,bp::bases<NGS_Object>,boost::noncopyable> ("NumProc", bp::no_init)
    .def("Do", FunctionPointer([](NumProc & self, int heapsize)
                               {
                                 LocalHeap lh (heapsize, "NumProc::Do-heap");
                                 self.Do(lh);
                               }),
         (bp::arg("self")=NULL,bp::arg("heapsize")=1000000))
    ;

  // die geht
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<NumProcWrap>);
  bp::class_<NumProcWrap,shared_ptr<NumProcWrap>, bp::bases<NumProc>,boost::noncopyable>("PyNumProc", bp::init<shared_ptr<PDE>, const Flags&>())
    .def("Do", bp::pure_virtual(&PyNumProc::Do)) 
    .add_property("pde", &PyNumProc::GetPDE)
    ;
  
  bp::implicitly_convertible 
    <shared_ptr<NumProcWrap>, shared_ptr<NumProc> >(); 


  //////////////////////////////////////////////////////////////////////////////////////////

  PyExportSymbolTable<PyFES> ();
  PyExportSymbolTable<PyCF> ();
  PyExportSymbolTable<PyGF> ();
  PyExportSymbolTable<PyBF> ();
  PyExportSymbolTable<PyLF> ();
  PyExportSymbolTable<shared_ptr<Preconditioner>> ();
  PyExportSymbolTable<shared_ptr<NumProc>> ();
  PyExportSymbolTable<double> ();
  PyExportSymbolTable<shared_ptr<double>> ();

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<PDE>);  
  bp::class_<PDE,shared_ptr<PDE>> ("PDE", bp::init<>())

    // .def(bp::init<const string&>())

#ifndef PARALLEL
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](const string & filename)
                           { 
                             return LoadPDE (filename);
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("filename"))
          ))

#else

    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](const string & filename,
                              bp::object py_mpicomm)
                           { 
                             PyObject * py_mpicomm_ptr = py_mpicomm.ptr();
                             if (py_mpicomm_ptr != Py_None)
                               {
                                 MPI_Comm * comm = PyMPIComm_Get (py_mpicomm_ptr);
                                 ngs_comm = *comm;
                               }
                             else
                               ngs_comm = MPI_COMM_WORLD;

                             cout << "Rank = " << MyMPI_GetId(ngs_comm) << "/"
                                  << MyMPI_GetNTasks(ngs_comm) << endl;

                             NGSOStream::SetGlobalActive (MyMPI_GetId()==0);
                             return LoadPDE (filename);
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("filename"), bp::arg("mpicomm")=bp::object())
          ))
#endif



    .def("LoadSolution", &PDE::LoadSolution,
         (bp::arg("filename"), bp::arg("ascii")=false)
         )

    
    /*
    .def("Load", 
         // static_cast<void(PDE::*)(const string &, const bool, const bool)> 
         // (&PDE::LoadPDE),
         FunctionPointer ([](shared_ptr<PDE> pde, const string & filename)
                          { 
                            LoadPDE (pde, filename);
                          }))
    */

    .def("__str__", &ToString<PDE>)

    .def("Mesh",  &PDE::GetMeshAccess,
         (bp::arg("meshnr")=0))

    .def("Solve", &PDE::Solve)


    .def("Add", FunctionPointer([](PDE & self, shared_ptr<MeshAccess> mesh)
                                {
                                  self.AddMeshAccess (mesh);
                                }))

    .def("Add", FunctionPointer([](PDE & self, const string & name, double val)
                                {
                                  self.AddConstant (name, val);
                                }))

    .def("Add", FunctionPointer([](PDE & self, PyWrapper<FESpace> space)
                                {
                                  self.AddFESpace (space->GetName(), space.Get());
                                }))

    .def("Add", FunctionPointer([](PDE & self, PyWrapperDerived<GridFunction, CoefficientFunction> gf)
                                {
                                  self.AddGridFunction (gf->GetName(), gf.Get());
                                }))

    .def("Add", FunctionPointer([](PDE & self, PyWrapper<BilinearForm> bf)
                                {
                                  self.AddBilinearForm (bf->GetName(), bf.Get());
                                }))

    .def("Add", FunctionPointer([](PDE & self, PyWrapper<LinearForm> lf)
                                {
                                  self.AddLinearForm (lf->GetName(), lf.Get());
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<Preconditioner> pre)
                                {
                                  self.AddPreconditioner (pre->GetName(), pre);
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<NumProcWrap> np)
                                {
                                  cout << "add pynumproc" << endl;
                                  self.AddNumProc ("pynumproc", np);
                                }))
    
    .def("Add", FunctionPointer([](PDE & self, shared_ptr<NumProc> np)
                                {
				  static int cnt = 0;
				  cnt++;
				  string name = "np_from_py" + ToString(cnt);
                                  self.AddNumProc (name, np);
                                }))

    .def("Add", FunctionPointer([](PDE & self, const bp::list &l)
                                {
                                  for (int i=0; i<bp::len(l); i++)
                                    {
                                      bp::extract<shared_ptr<PyNumProc>> np(l[i]);
                                      if(np.check())
                                        {
                                          self.AddNumProc (np()->GetName(), np());
                                          continue;
                                        }
                                      
                                      bp::extract<shared_ptr<NumProc>> pnp(l[i]);
                                      if(np.check())
                                        {
                                          self.AddNumProc (pnp()->GetName(), pnp());
                                          continue;
                                        }
                                      
                                      bp::extract<shared_ptr<GridFunction>> gf(l[i]);
                                      if(gf.check())
                                        {
                                          self.AddGridFunction (gf()->GetName(), gf());
                                          continue;
                                        }
                                      
                                      bp::extract<shared_ptr<BilinearForm>> bf(l[i]);
                                      if(gf.check())
                                        {
                                          self.AddBilinearForm (bf()->GetName(), bf());
                                          continue;
                                        }
                                      
                                      bp::extract<shared_ptr<LinearForm>> lf(l[i]);
                                      if(gf.check())
                                        {
                                          self.AddLinearForm (lf()->GetName(), lf());
                                          continue;
                                        }
                                      
                                      bp::extract<shared_ptr<Preconditioner>> pre(l[i]);
                                      if(gf.check())
                                        {
                                          self.AddPreconditioner (pre()->GetName(), pre());
                                          continue;
                                        }
                                      
                                      cout << "warning: unknown object at position " << i << endl;
                                    }
                                }))

    .def("SetCurveIntegrator", FunctionPointer
         ([](PDE & self, const string & filename, PyWrapper<LinearFormIntegrator> lfi)
          {
            self.SetLineIntegratorCurvePointInfo(filename, lfi.Get().get());
          }))

    .add_property ("constants", FunctionPointer([](PDE & self) { return bp::object(self.GetConstantTable()); }))
    .add_property ("variables", FunctionPointer([](PDE & self) { return bp::object(self.GetVariableTable()); }))
    .add_property ("coefficients", FunctionPointer([](PDE & self) { return bp::object(self.GetCoefficientTable()); }))
    .add_property ("spaces", FunctionPointer([](PDE & self) {
          auto table = self.GetSpaceTable();
          SymbolTable<PyFES> pytable;
          for ( auto i : Range(table.Size() ))
                pytable.Set(table.GetName(i), PyFES(table[i]));
          return bp::object(pytable);
          }))
    .add_property ("gridfunctions", FunctionPointer([](PDE & self) {
          auto table = self.GetGridFunctionTable();
          SymbolTable<PyGF> pytable;
          for ( auto i : Range(table.Size() ))
                pytable.Set(table.GetName(i), PyGF(table[i]));
          return bp::object(pytable);
          }))
    .add_property ("bilinearforms", FunctionPointer([](PDE & self) {
          auto table = self.GetBilinearFormTable();
          SymbolTable<PyBF> pytable;
          for ( auto i : Range(table.Size() ))
                pytable.Set(table.GetName(i), PyBF(table[i]));
          return bp::object(pytable);
          }))
    .add_property ("linearforms", FunctionPointer([](PDE & self) {
          auto table = self.GetLinearFormTable();
          SymbolTable<PyLF> pytable;
          for ( auto i : Range(table.Size() ))
                pytable.Set(table.GetName(i), PyLF(table[i]));
          return bp::object(pytable);
          }))
    .add_property ("preconditioners", FunctionPointer([](PDE & self) { return bp::object(self.GetPreconditionerTable()); }))
    .add_property ("numprocs", FunctionPointer([](PDE & self) { return bp::object(self.GetNumProcTable()); }))
    ;
  
  bp::def("Integrate", 
          FunctionPointer([](PyCF cf,
                             shared_ptr<MeshAccess> ma, 
                             VorB vb, int order, 
                             bool region_wise, bool element_wise)
                          {
                            static mutex addcomplex_mutex;
                            LocalHeap lh(1000000, "lh-Integrate");
                            
                            if (!cf->IsComplex())
                              {
                                atomic<double> sum(0.0);
                                Vector<> region_sum(ma->GetNRegions(vb));
                                Vector<> element_sum(element_wise ? ma->GetNE(vb) : 0);
                                region_sum = 0;
                                element_sum = 0;

                                bool use_simd = true;
                                
                                ma->IterateElements
                                  (vb, lh, [&] (Ngs_Element el, LocalHeap & lh)
                                   {
                                     auto & trafo = ma->GetTrafo (el, lh);
                                     double hsum = 0.0;
                                     bool this_simd = use_simd;
                                     
                                     if (this_simd)
                                       {
                                         try
                                           {
                                             SIMD_IntegrationRule ir(trafo.GetElementType(), order);
                                             auto & mir = trafo(ir, lh);
                                             AFlatMatrix<> values(1, ir.GetNIP(), lh);
                                             cf.Get() -> Evaluate (mir, values);
                                             SIMD<double> vsum = 0;
                                             for (int i = 0; i < values.VWidth(); i++)
                                               vsum += mir[i].GetWeight() * values.Get(0,i);
                                             hsum = HSum(vsum);
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
                                         FlatMatrix<> values(ir.Size(), 1, lh);
                                         cf.Get() -> Evaluate (mir, values);
                                         for (int i = 0; i < values.Height(); i++)
                                           hsum += mir[i].GetWeight() * values(i,0);
                                       }

                                     sum += hsum;
                                     double & rsum = region_sum(el.GetIndex());
				     AsAtomic(rsum) += hsum;
                                     if (element_wise)
                                       element_sum(el.Nr()) = hsum;
                                   });
                                bp::object result;
                                if (region_wise)
                                  result = bp::list(bp::object(region_sum));
                                else if (element_wise)
                                  result = bp::object(element_sum);
                                else
                                  result = bp::object(sum.load());
                                return result;
                              }
                            else
                              {
                                Complex sum = 0;
                                Vector<Complex> region_sum(ma->GetNRegions(vb));
                                Vector<Complex> element_sum(element_wise ? ma->GetNE(vb) : 0);
                                region_sum = 0;
                                element_sum = 0;
                                
                                ma->IterateElements
                                  (vb, lh, [&] (Ngs_Element el, LocalHeap & lh)
                                   {
                                     auto & trafo = ma->GetTrafo (el, lh);
                                     IntegrationRule ir(trafo.GetElementType(), order);
                                     BaseMappedIntegrationRule & mir = trafo(ir, lh);
                                     FlatMatrix<Complex> values(ir.Size(), 1, lh);
                                     cf.Get() -> Evaluate (mir, values);
                                     Complex hsum = 0;
                                     for (int i = 0; i < values.Height(); i++)
                                       hsum += mir[i].GetWeight() * values(i,0);
                                     {
                                       lock_guard<mutex> guard(addcomplex_mutex);
                                       sum += hsum;
                                       region_sum(el.GetIndex()) += hsum;
                                     }
                                     if (element_wise)
                                       element_sum(el.Nr()) = hsum;
                                   });
                                bp::object result;
                                if (region_wise)
                                  result = bp::list(bp::object(region_sum));
                                else if (element_wise)
                                  result = bp::object(element_sum);
                                else
                                  result = bp::object(sum);
                                return result;
                              }
                          }),
          (bp::arg("cf"), bp::arg("mesh"), bp::arg("VOL_or_BND")=VOL, 
           bp::arg("order")=5, 
           bp::arg("region_wise")=false,
           bp::arg("element_wise")=false))
    ;
  






  bp::def("SymbolicLFI", FunctionPointer
          ([](PyCF cf, VorB vb, bool element_boundary,
              bool skeleton, bp::object definedon) 
           {
             bp::extract<Region> defon_region(definedon);
             if (defon_region.check())
               vb = VorB(defon_region());

             shared_ptr<LinearFormIntegrator> lfi;
             if (!skeleton)
               lfi = make_shared<SymbolicLinearFormIntegrator> (cf.Get(), vb, element_boundary);
             else
               lfi = make_shared<SymbolicFacetLinearFormIntegrator> (cf.Get(), vb /* , element_boundary */);
             
             if (bp::extract<bp::list> (definedon).check())
               lfi -> SetDefinedOn (makeCArray<int> (definedon));
             if (defon_region.check())
               lfi->SetDefinedOn(defon_region().Mask());

             return PyWrapper<LinearFormIntegrator>(lfi);
           }),
          (bp::args("form"),
           bp::args("VOL_or_BND")=VOL,
           bp::args("element_boundary")=false,
           bp::args("skeleton")=false,           
           bp::arg("definedon")=bp::object())
          );

  bp::def("SymbolicBFI", FunctionPointer
          ([](PyCF cf, VorB vb, bool element_boundary,
              bool skeleton, bp::object definedon)
           {
             bp::extract<Region> defon_region(definedon);
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
               bfi = make_shared<SymbolicBilinearFormIntegrator> (cf.Get(), vb, element_boundary);
             else
               bfi = make_shared<SymbolicFacetBilinearFormIntegrator> (cf.Get(), vb, element_boundary);
             
             if (bp::extract<bp::list> (definedon).check())
               bfi -> SetDefinedOn (makeCArray<int> (definedon));

             if (defon_region.check())
               {
                 cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
                 bfi->SetDefinedOn(defon_region().Mask());
               }
             
             return PyWrapper<BilinearFormIntegrator>(bfi);
           }),
          (bp::args("form"), bp::args("VOL_or_BND")=VOL,
           bp::args("element_boundary")=false,
           bp::args("skeleton")=false,
           bp::arg("definedon")=bp::object())
          );

  bp::def("SymbolicEnergy", FunctionPointer
          ([](PyCF cf, VorB vb, bp::object definedon) -> PyWrapper<BilinearFormIntegrator>
           {
             bp::extract<Region> defon_region(definedon);
             if (defon_region.check())
               vb = VorB(defon_region());

             auto bfi = make_shared<SymbolicEnergy> (cf.Get(), vb);
             
             if (defon_region.check())
               {
                 cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
                 bfi->SetDefinedOn(defon_region().Mask());
               }
             
             /*
             bp::extract<bp::list> defon_list(definedon);
             if (defon_list.check())
               {
                 BitArray bits(bp::len (defon_list));
                 bits.Clear();
                 bool all_booleans = true;
                 for (int i : Range(bits))
                   {
                     cout << "class = " << defon_list().attr("__class__") << endl;
                     bp::extract<bool> extbool(defon_list()[i]);
                     if (extbool.check())
                       {
                         if (extbool()) bits.Set(i);
                       }
                     else
                       all_booleans = false;
                   }
                 cout << "bits: " << bits << endl;
                 cout << "allbool = " << all_booleans << endl;
               }
             */
             return PyWrapper<BilinearFormIntegrator>(bfi);
           }),
          (bp::args("self"), bp::args("VOL_or_BND")=VOL, bp::args("definedon")=bp::object())
          );


  /*
  bp::def("IntegrateLF", 
          FunctionPointer
          ([](shared_ptr<LinearForm> lf, 
              shared_ptr<CoefficientFunction> cf)
           {
             lf->AllocateVector();
             lf->GetVector() = 0.0;

             Array<ProxyFunction*> proxies;
             cf->TraverseTree( [&] (CoefficientFunction & nodecf)
                               {
                                 auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
                                 if (proxy && !proxies.Contains(proxy))
                                   proxies.Append (proxy);
                               });
             
             LocalHeap lh1(1000000, "lh-Integrate");

             // for (auto el : lf->GetFESpace()->Elements(VOL, lh))
             IterateElements 
               (*lf->GetFESpace(), VOL, lh1,
                [&] (FESpace::Element el, LocalHeap & lh)
               {
                 const FiniteElement & fel = el.GetFE();
                 auto & trafo = lf->GetMeshAccess()->GetTrafo (el, lh);
                 IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
                 BaseMappedIntegrationRule & mir = trafo(ir, lh);
                 FlatVector<> elvec(fel.GetNDof(), lh);
                 FlatVector<> elvec1(fel.GetNDof(), lh);

                 FlatMatrix<> values(ir.Size(), cf->Dimension(), lh);
                 ProxyUserData ud;
                 trafo.userdata = &ud;

                 elvec = 0;
                 for (auto proxy : proxies)
                   {
                     FlatMatrix<> proxyvalues(ir.Size(), proxy->Dimension(), lh);
                     for (int k = 0; k < proxy->Dimension(); k++)
                       {
                         ud.testfunction = proxy;
                         ud.test_comp = k;
                         
                         cf -> Evaluate (mir, values);
                         for (int i = 0; i < mir.Size(); i++)
                           values.Row(i) *= mir[i].GetWeight();
                         proxyvalues.Col(k) = values.Col(0);
                       }

                     proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
                     elvec += elvec1;
                   }
                 lf->AddElementVector (el.GetDofs(), elvec);
               });
           }));
           


  bp::def("IntegrateBF", 
          FunctionPointer
          ([](shared_ptr<BilinearForm> bf1, 
              shared_ptr<CoefficientFunction> cf)
           {
             auto bf = dynamic_pointer_cast<S_BilinearForm<double>> (bf1);
             bf->GetMatrix().SetZero();

             Array<ProxyFunction*> trial_proxies, test_proxies;
             cf->TraverseTree( [&] (CoefficientFunction & nodecf)
                               {
                                 auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
                                 if (proxy) 
                                   {
                                     if (proxy->IsTestFunction())
                                       {
                                         if (!test_proxies.Contains(proxy))
                                           test_proxies.Append (proxy);
                                       }
                                     else
                                       {                                         
                                         if (!trial_proxies.Contains(proxy))
                                           trial_proxies.Append (proxy);
                                       }
                                   }
                               });

             ProxyUserData ud;
             LocalHeap lh(1000000, "lh-Integrate");

             // IterateElements (*lf->GetFESpace(), VOL, lh,
             for (auto el : bf->GetFESpace()->Elements(VOL, lh))
               {
                 const FiniteElement & fel = el.GetFE();
                 auto & trafo = bf->GetMeshAccess()->GetTrafo (el, lh);
                 trafo.userdata = &ud;
                 IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
                 BaseMappedIntegrationRule & mir = trafo(ir, lh);
                 FlatMatrix<> elmat(fel.GetNDof(), lh);

                 FlatMatrix<> values(ir.Size(), 1, lh);

                 elmat = 0;

                 for (int i = 0; i < mir.Size(); i++)
                   {
                     auto & mip = mir[i];
                     
                     for (auto proxy1 : trial_proxies)
                       for (auto proxy2 : test_proxies)
                         {
                           HeapReset hr(lh);

                           FlatMatrix<> proxyvalues(proxy2->Dimension(), 
                                                    proxy1->Dimension(), 
                                                    lh);
                           for (int k = 0; k < proxy1->Dimension(); k++)
                             for (int l = 0; l < proxy2->Dimension(); l++)
                               {
                                 ud.trialfunction = proxy1;
                                 ud.trial_comp = k;
                                 ud.testfunction = proxy2;
                                 ud.test_comp = l;
                                 proxyvalues(l,k) = 
                                   mip.GetWeight() * cf -> Evaluate (mip);
                               }
                           
                           FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), fel.GetNDof(), lh);
                           FlatMatrix<double,ColMajor> dbmat1(proxy1->Dimension(), fel.GetNDof(), lh);
                           FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), fel.GetNDof(), lh);

                           proxy1->Evaluator()->CalcMatrix(fel, mip, bmat1, lh);
                           proxy2->Evaluator()->CalcMatrix(fel, mip, bmat2, lh);
                           dbmat1 = proxyvalues * bmat1;
                           elmat += Trans (bmat2) * dbmat1;
                         }
                   }
                 bf->AddElementMatrix (el.GetDofs(), el.GetDofs(), elmat, el, lh);
               }
           }));
  */

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<BaseVTKOutput>);
  bp::class_<BaseVTKOutput, shared_ptr<BaseVTKOutput>,  boost::noncopyable>("VTKOutput", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](shared_ptr<MeshAccess> ma, bp::list coefs_list,
                              bp::list names_list, string filename, int subdivision, int only_element)
                           {
                             Array<shared_ptr<CoefficientFunction> > coefs
                               = makeCArrayUnpackWrapper<PyCF> (coefs_list);
                             Array<string > names
                               = makeCArray<string> (names_list);
                             shared_ptr<BaseVTKOutput> ret;
                             if (ma->GetDimension() == 2)
                               ret = make_shared<VTKOutput<2>> (ma, coefs, names, filename, subdivision, only_element);
                             else
                               ret = make_shared<VTKOutput<3>> (ma, coefs, names, filename, subdivision, only_element);
                             return ret;
                           }),

          bp::default_call_policies(),     // need it to use named arguments
          (
            bp::arg("ma"),
            bp::arg("coefs")= bp::list(),
            bp::arg("names") = bp::list(),
            bp::arg("filename") = "vtkout",
            bp::arg("subdivision") = 0,
            bp::arg("only_element") = -1
            )
           )
      )

    .def("Do", FunctionPointer([](BaseVTKOutput & self, int heapsize)
                               { 
                                 LocalHeap lh (heapsize, "VTKOutput-heap");
                                 self.Do(lh);
                               }),
         (bp::arg("self"),bp::arg("heapsize")=1000000))
    .def("Do", FunctionPointer([](BaseVTKOutput & self, const BitArray * drawelems, int heapsize)
                               { 
                                 LocalHeap lh (heapsize, "VTKOutput-heap");
                                 self.Do(lh, drawelems);
                               }),
         (bp::arg("self"),bp::arg("drawelems"),bp::arg("heapsize")=1000000))
    
    ;



#ifdef PARALLEL
  import_mpi4py();
#endif
}





BOOST_PYTHON_MODULE(libngcomp) 
{
  ExportNgcomp();
}



#endif // NGS_PYTHON
