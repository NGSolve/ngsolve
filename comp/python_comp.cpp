#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <comp.hpp>
using namespace ngcomp;



template <> class cl_NonElement<ElementId>
{
public:
  static ElementId Val() { return ElementId(VOL,-1); }
};




void ExportNgcomp() 
{
  std::string nested_name = "ngcomp";
  if( bp::scope() )
    nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".ngcomp");
  
  bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
  
  cout << "exporting ngstd as " << nested_name << endl;
  bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
  parent.attr("ngcomp") = module;
  
  bp::scope local_scope(module);
    



  bp::enum_<VorB>("VorB")
    .value("VOL", VOL)
    .value("BND", BND)
    .export_values()
    ;

  bp::enum_<COUPLING_TYPE> ("COUPLING_TYPE")
    .value("UNUSED_DOF", UNUSED_DOF)
    .value("LOCAL_DOF", LOCAL_DOF)
    .value("INTERFACE_DOF", INTERFACE_DOF)
    .value("NONWIREBASKET_DOF", NONWIREBASKET_DOF)
    .value("WIREBASKET_DOF", WIREBASKET_DOF)
    .value("EXTERNAL_DOF", EXTERNAL_DOF)
    .value("ANY_DOF", ANY_DOF)
    .export_values()
    ;

  bp::class_<ElementId> ("ElementId", bp::init<VorB,int>())
    .def("__str__", &ToString<ElementId>)
    .add_property("nr", &ElementId::Nr)    
    .def("IsVolume", &ElementId::IsVolume)
    .def("IsBoundary", &ElementId::IsBoundary)
    ;
  //
  bp::class_<ElementRange,bp::bases<IntRange>> ("ElementRange",bp::init<VorB,IntRange>())
    .def(PyDefIterable2<ElementRange>())
    // .def("__iter__", bp::iterator<ElementRange>())
    ;

  bp::class_<Ngs_Element>("Ngs_Element", bp::no_init)
    .add_property("vertices", FunctionPointer([](Ngs_Element &el) {return bp::tuple(Array<int>(el.Vertices()));} ))
    .add_property("edges", FunctionPointer([](Ngs_Element &el) { return bp::tuple(Array<int>(el.Edges()));} ))
    .add_property("faces", FunctionPointer([](Ngs_Element &el) { return bp::tuple(Array<int>(el.Faces()));} ))
    .add_property("type", &Ngs_Element::GetType)
    ;

  bp::class_<MeshAccess>("Mesh", "meshclass doc", bp::init<string>())
    .def("LoadMesh", static_cast<void(MeshAccess::*)(const string &)>(&MeshAccess::LoadMesh))
    .def("Elements", static_cast<ElementRange(MeshAccess::*)(VorB)const> (&MeshAccess::Elements))
    .def("GetElement", static_cast<Ngs_Element(MeshAccess::*)(ElementId)const> (&MeshAccess::GetElement))
    .def("__getitem__", static_cast<Ngs_Element(MeshAccess::*)(ElementId)const> (&MeshAccess::operator[]))

    .def ("GetNE", static_cast<int(MeshAccess::*)(VorB)const> (&MeshAccess::GetNE))
    .add_property ("nv", &MeshAccess::GetNV, "number of vertices")
    .add_property ("ne", static_cast<int(MeshAccess::*)()const> (&MeshAccess::GetNE))

    .def ("GetTrafo", 
          static_cast<ElementTransformation&(MeshAccess::*)(ElementId,LocalHeap&)const>
          (&MeshAccess::GetTrafo), 
          bp::return_value_policy<bp::reference_existing_object>())

    // first attempts, but keep for a moment ...
    .def("GetElement", static_cast< Ngs_Element (MeshAccess::*)(int, bool)const> (&MeshAccess::GetElement),
         (bp::arg("arg1")=NULL,bp::arg("arg2")=0))
    ;

  
  bp::class_<NGS_Object, shared_ptr<NGS_Object>,  boost::noncopyable>("NGS_Object", bp::no_init)
    .add_property("name", FunctionPointer([](const NGS_Object & self)->string
                                          {
                                            return self.GetName();
                                          }))
    ;


  bp::class_<FESpace, shared_ptr<FESpace>,  boost::noncopyable>("FESpace", bp::no_init)
    .def("Update", FunctionPointer([](FESpace & self, int heapsize)
                                   {
                                     LocalHeap lh (heapsize, "FESpace::Update-heap");
                                     self.Update(lh);
                                     self.FinalizeUpdate(lh);
                                   }),
         (bp::arg("self")=NULL,bp::arg("heapsize")=1000000))
    
    .def("GetDofNrs", FunctionPointer([](FESpace & self, ElementId i) 
                                   {
                                     Array<int> tmp; self.GetDofNrs(i,tmp); 
                                     return bp::tuple (tmp); 
                                   }))

    .def("CouplingType", &FESpace::GetDofCouplingType)
    .add_property ("ndof", FunctionPointer([](FESpace & self) { return self.GetNDof(); }))
    .def(PyDefToString<FESpace>())
    .def ("GetFE", 
          static_cast<const FiniteElement&(FESpace::*)(ElementId,LocalHeap&)const>
          (&FESpace::GetFE), 
          bp::return_value_policy<bp::reference_existing_object>())

    .def("FreeDofs",
         FunctionPointer( [] (const FESpace &self, bool coupling) -> const BitArray &{ return *self.GetFreeDofs(coupling); } ),
         bp::return_value_policy<bp::reference_existing_object>(),
         (bp::arg("self"), bp::arg("coupling")=0))
    ;
  
  typedef GridFunction GF;
  bp::class_<GF, shared_ptr<GF>, boost::noncopyable>
    ("GridFunction",  "a field approximated in some finite element space", bp::no_init)

    .def("__init__", bp::make_constructor
         (FunctionPointer ([](string name, shared_ptr<FESpace> fespace)
                           {
                             Flags flags;
                             GridFunction * gf = CreateGridFunction (fespace, name, flags);
                             return shared_ptr<GridFunction> (gf);
                           })),
         // (bp::arg("self"), bp::arg("name"), bp::arg("fespace")), // how ?
         "creates a gridfunction in finite element space"
         )

    .def("__str__", &ToString<GF>)
    .add_property("space", &GF::FESpacePtr, "the finite element spaces")
    .def("Update", 
         FunctionPointer([](GF & self) {self.Update();}),
         "update vector size to finite element space dimension after mesh refinement"
         )

    .def("Set",
         FunctionPointer([](GF & self, const CoefficientFunction & cf)
                         {
                          LocalHeap lh(1000000, "tmplh");
                          SetValues (cf, self, false, NULL, lh);
                        }))

    .def("Vector", 
         FunctionPointer([](GF & self)->BaseVector& { return self.GetVector(); }),
         bp::return_value_policy<bp::reference_existing_object>(),
         "vector of coefficients"
         )

    .add_property("vec",
         FunctionPointer([](GF & self)-> shared_ptr<BaseVector> { return shared_ptr<BaseVector>(&self.GetVector(), &NOOP_Deleter); }),
         "vector of coefficients of first multidim dimension"
         )

    .add_property("vecs",
         FunctionPointer([](GF & self)-> bp::list { 
             bp::list vecs;
             for (int i=0; i<self.GetMultiDim(); i++) {
                 vecs.append(shared_ptr<BaseVector>(&self.GetVector(i), &NOOP_Deleter));
             }
             return vecs;
             }),
         "list of vectors of coefficients"
         )
    ;
  
  /*
  bp::def("CreateGF",
          FunctionPointer
          ([](string name, shared_ptr<FESpace> fespace)
           {
             Flags flags;
             GridFunction * gf = CreateGridFunction (fespace, name, flags);
             return shared_ptr<GridFunction> (gf);
           }),
          (bp::arg("name"), bp::arg("fespace")),
          "creates a gridfunction in finite element space");
  */

  bp::class_<BilinearForm, shared_ptr<BilinearForm>, boost::noncopyable>("BilinearForm", bp::no_init)
    .def(PyDefToString<BilinearForm>())
    .def("Matrix", FunctionPointer([](shared_ptr<BilinearForm> self)->BaseMatrix& { return self->GetMatrix(); }),
         bp::return_value_policy<bp::reference_existing_object>())
    ;

  bp::class_<LinearForm, shared_ptr<LinearForm>, boost::noncopyable>("LinearForm", bp::no_init)
    .def(PyDefToString<LinearForm>())
    .def("Vector", FunctionPointer([](shared_ptr<LinearForm> self)->BaseVector& { return self->GetVector(); }),
         bp::return_value_policy<bp::reference_existing_object>())
    ;


}





void ExportNgstd();
void ExportNgbla();
void ExportNgfem();

BOOST_PYTHON_MODULE(libngcomp) {
  ExportNgstd();
  ExportNgbla();
  ExportNgfem();
  ExportNgcomp();
}






#endif // NGS_PYTHON
