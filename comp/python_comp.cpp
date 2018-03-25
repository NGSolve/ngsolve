#ifdef NGS_PYTHON

#include <regex>

#include "../ngstd/python_ngstd.hpp"
#include <comp.hpp>

#include "hdivdivfespace.hpp"
#include "hdivdivsurfacespace.hpp"
#include "numberfespace.hpp"
using namespace ngcomp;

using ngfem::ELEMENT_TYPE;

typedef GridFunction GF;

template <> class cl_NonElement<ElementId>
{
public:
  static ElementId Val() { return ElementId(VOL,-1); }
};


template <typename T>
py::tuple MakePyTuple (const BaseArrayObject<T> & ao)
{
  size_t s = ao.Size();
  py::tuple tup(s);
  for (size_t i = 0; i < s; i++)
    tup[i] = ao[i];
  return tup;
}

inline auto Nr2Vert(size_t nr) {  return NodeId(NT_VERTEX,nr); };
inline auto Nr2Edge(size_t nr) {  return NodeId(NT_EDGE,nr); };
inline auto Nr2Face(size_t nr) {  return NodeId(NT_FACE,nr); };


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

template<typename FESPACE>
shared_ptr<FESPACE> fesUnpickle(py::tuple state)
{
  auto fes = CreateFESpace(state[0].cast<string>(),
                           state[1].cast<shared_ptr<MeshAccess>>(),
                           state[2].cast<Flags>());
  LocalHeap lh (1000000, "FESpace::Update-heap");
  fes->Update(lh);
  fes->FinalizeUpdate(lh);
  return dynamic_pointer_cast<FESPACE>(fes);
};

void ExportPml(py::module &m)
{
  typedef CoefficientFunction CF;
  typedef PML_Transformation PML;
  py::class_<PML, shared_ptr<PML>>(m, "PML", R"raw_string(Base PML object

can only be created by generator functions. Use PML(x, [y, z]) to evaluate the scaling.)raw_string")
    .def("__call__",  [](py::args varargs) {
                      auto self = py::extract<shared_ptr<PML>>(varargs[0])();
                      int dim = self->GetDimension();
                      Vector<double> hpoint(dim);
                      hpoint = 0.;
                      for (int i : Range(min(int(py::len(varargs)-1),dim)))
                        hpoint[i] = py::extract<double>(varargs[i+1])();
                      Vector<Complex> point(dim);
                      Matrix<Complex> jac(dim,dim);
                      self->MapPointV(hpoint,point,jac);
                      return point;
                    },"map a point")
    .def("__str__", [] (shared_ptr<PML> self) { return ToString(*self); } )
    .def("call_jacobian",  [](py::args varargs) {
                      auto self = py::extract<shared_ptr<PML>>(varargs[0])();
                      int dim = self->GetDimension();
                      Vector<double> hpoint(dim);
                      hpoint = 0.;
                      for (int i : Range(min(int(py::len(varargs)-1),dim)))
                        hpoint[i] = py::extract<double>(varargs[i+1])();
                      Vector<Complex> point(dim);
                      Matrix<Complex> jac(dim,dim);
                      self->MapPointV(hpoint,point,jac);
                      return jac;
                    },"evaluate PML jacobian at point x, [y, z]")
    .def_property_readonly("dim", [] (shared_ptr<PML> self) {return self->GetDimension(); },
        "dimension")
    .def_property_readonly("PML_CF", [](shared_ptr<PML> self)->shared_ptr<CF> {
        return make_shared<PML_CF> (self);
      },
      "the scaling as coefficient function")
    .def_property_readonly("Jac_CF", [](shared_ptr<PML>self)->shared_ptr<CF> {
        return make_shared<PML_Jac> (self);
      },
      "the jacobian of the PML as coefficient function")
    .def_property_readonly("Det_CF", [](shared_ptr<PML> self)->shared_ptr<CF> {
        return make_shared<PML_Det> (self);
      },
      "the determinant of the jacobian as coefficient function")
    .def_property_readonly("JacInv_CF", [](shared_ptr<PML> self)->shared_ptr<CF> {
        return make_shared<PML_JacInv> (self);
      },
      "the inverse of the jacobian as coefficient function")
    .def("__add__", [](shared_ptr<PML> pml1, shared_ptr<PML> pml2)
         -> shared_ptr<PML>
         {
        int dim = pml1->GetDimension();
        if (pml2->GetDimension() != dim)
          throw Exception("Dimensions do not match");
        switch (dim)
          {
          case 1:
            return make_shared<SumPML<1>> (pml1,pml2);
          case 2:
            return make_shared<SumPML<2>> (pml1,pml2);
          case 3:
            return make_shared<SumPML<3>> (pml1,pml2);
          }
        throw Exception("No valid dimension");
      })
    ;

  m.def("Radial", [](py::object _origin, double rad, Complex alpha) -> shared_ptr<PML>{
      Vector<double> origin;
      int dim = 0;
      if (py::extract<double>(_origin).check())
        {
          dim = 1;
          origin.SetSize(1);
          origin(0)=py::extract<double>(_origin)();
        }
      else if (py::extract<py::tuple>(_origin).check())
        {
          py::tuple torigin(_origin);
          dim = py::len(torigin);
          origin.SetSize(dim);
          for (int j : Range(dim))
            origin(j)=py::extract<double>(torigin[j])();
        }
      switch (dim)
        {
        case 1:
          return make_shared<RadialPML_Transformation<1>> (rad,alpha,origin);
        case 2:
          return make_shared<RadialPML_Transformation<2>> (rad,alpha,origin);
        case 3:
          return make_shared<RadialPML_Transformation<3>> (rad,alpha,origin);
        }
      throw Exception("No valid dimension");
    },
    py::arg("origin"),py::arg("rad")=1,py::arg("alpha")=Complex(0,1),
    R"raw_string(radial pml transformation

origin is a list/tuple determining the dimenson)raw_string");

  m.def("Custom", [](shared_ptr<CF> trafo, shared_ptr<CF> jac) -> shared_ptr<PML>{
      switch (trafo->Dimension())
        {
        case 1:
          return make_shared<CustomPML_Transformation<1>> (trafo,jac);
        case 2:
          return make_shared<CustomPML_Transformation<2>> (trafo,jac);
        case 3:
          return make_shared<CustomPML_Transformation<3>> (trafo,jac);
        }
      throw Exception("No valid dimension");
    },
    py::arg("trafo"),py::arg("jac"),
    R"raw_string(custom pml transformation

trafo and jac are coefficient functions of the scaling and the jacobian)raw_string")
    ;
  m.def("Cartesian", [](py::object mins,py::object maxs, Complex alpha) -> shared_ptr<PML>{
      int dim = 0;
      Matrix<double> bounds;
      if (py::extract<double>(mins).check())
        {
          dim = 1;
          bounds.SetSize(dim,2);
          bounds = 0.;
          bounds(0,0)=py::extract<double>(mins)();
        }
      else if (py::extract<py::tuple>(mins).check())
        {
          py::tuple tmins(mins);
          dim = py::len(tmins);
          bounds.SetSize(dim,2);
          bounds = 0.;
          for (int j : Range(dim))
            bounds(j,0)=py::extract<double>(tmins[j])();
        }

      if (py::extract<double>(maxs).check())
        bounds(0,1)=py::extract<double>(maxs)();

          else if (py::extract<py::tuple>(maxs).check())
          {
            py::tuple tmax(maxs);
            for (int j : Range(min(int(py::len(tmax)),dim)))
              bounds(j,1)=py::extract<double>(tmax[j])();
          }
          switch (dim)
          {
            case 1:
              return make_shared<CartesianPML_Transformation<1>> (bounds,alpha);
            case 2:
              return make_shared<CartesianPML_Transformation<2>> (bounds,alpha);
            case 3:
              return make_shared<CartesianPML_Transformation<3>> (bounds,alpha);
           }
          throw Exception("No valid dimension");
        },
        py::arg("mins"),py::arg("maxs"), py::arg("alpha")=Complex(0,1),
        R"raw_string(cartesian pml transformation

mins and maxs are tuples/lists determining the dimension)raw_string")
    ;
  m.def("HalfSpace", [](py::object point,py::object normal, Complex alpha) -> shared_ptr<PML>{
          int dim = 0;
          Vector<double> vpoint;
          Vector<double> vnormal;
          if (py::extract<double>(point).check())
          {
            dim = 1;
            vpoint.SetSize(dim);
            vpoint = 0.;
            vnormal.SetSize(dim);
            vnormal = 0.;
            vpoint(0)=py::extract<double>(point)();
          }
          else if (py::extract<py::tuple>(point).check())
          {
            py::tuple tpoint(point);
            dim = py::len(tpoint);
            vpoint.SetSize(dim);
            vnormal.SetSize(dim);
            vpoint = 0.;
            vnormal = 0.;
            for (int j : Range(dim))
              vpoint(j)=py::extract<double>(tpoint[j])();
          }

          if(py::extract<double>(normal).check())
          {
            dim = 1;
            vnormal(0)=py::extract<double>(normal)();
          }
          else if (py::extract<py::tuple>(normal).check())
          {
            py::tuple tnormal(normal);
            dim = py::len(tnormal);
            for (int j : Range(min(int(py::len(tnormal)),dim)))
              vnormal(j)=py::extract<double>(tnormal[j])();
          }
          switch (dim)
          {
            case 1:
              return make_shared<HalfSpacePML_Transformation<1>> (vpoint,vnormal,alpha);
            case 2:
              return make_shared<HalfSpacePML_Transformation<2>> (vpoint,vnormal,alpha);
            case 3:
              return make_shared<HalfSpacePML_Transformation<3>> (vpoint,vnormal,alpha);
          }
          throw Exception("No valid dimension");
        },
        py::arg("point"),py::arg("normal"), py::arg("alpha")=Complex(0,1),
        R"raw_string(half space pml

scales orthogonal to specified plane in direction of normal point and normal are given as tuples/lists determining the dimension)raw_string")
    ;
    m.def("BrickRadial", [](py::object mins,py::object maxs,py::object _origin, Complex alpha) {
          int dim = 0;
          Matrix<double> bounds;
          if (py::extract<double>(mins).check())
          {
            dim = 1;
            bounds.SetSize(dim,2);
            bounds = 0.;
            bounds(0,0)=py::extract<double>(mins)();
          }
          else if (py::extract<py::tuple>(mins).check())
          {
            py::tuple tmins(mins);
            dim = py::len(tmins);
            bounds.SetSize(dim,2);
            bounds = 0.;
            for (int j : Range(dim))
              bounds(j,0)=py::extract<double>(tmins[j])();
          }

          if (py::extract<double>(maxs).check())
              bounds(0,1)=py::extract<double>(maxs)();

          else if (py::extract<py::tuple>(maxs).check())
          {
            py::tuple tmax(maxs);
            for (int j : Range(min(int(py::len(tmax)),dim)))
              bounds(j,1)=py::extract<double>(tmax[j])();
          }
          Vector<double> vorigin(dim);
          vorigin = 0.;
          if (py::extract<double>(_origin).check())
          {
            vorigin(0)=py::extract<double>(_origin)();
          }
          else if (py::extract<py::tuple>(_origin).check())
          {
            py::tuple torigin(_origin);
            for (int j : Range(min(int(py::len(torigin)),dim)))
              vorigin(j)=py::extract<double>(torigin[j])();
          }
          switch (dim)
          {
            case 1:
              return shared_ptr<PML>(make_shared<BrickRadialPML_Transformation<1>> (bounds,alpha,vorigin));
            case 2:
              return shared_ptr<PML>(make_shared<BrickRadialPML_Transformation<2>> (bounds,alpha,vorigin));
            case 3:
              return shared_ptr<PML>(make_shared<BrickRadialPML_Transformation<3>> (bounds,alpha,vorigin));
          }
          throw Exception("No valid dimension");
        },
        py::arg("mins"),py::arg("maxs"), py::arg("origin")=py::make_tuple(0.,0.,0.),py::arg("alpha")=Complex(0,1),
        R"raw_string(radial pml on a brick

mins, maxs and origin are given as tuples/lists)raw_string")
      ;
    m.def("Compound", [](shared_ptr<PML> pml1,shared_ptr<PML> pml2,py::object dims1,py::object dims2)
          ->shared_ptr<PML>
          {
          int dim1 = pml1->GetDimension();
          int dim2 = pml2->GetDimension();
          int dim = dim1 + dim2;
          Vector<int> vdims1;
          Vector<int> vdims2;
          
          if (py::extract<double>(dims1).check())
          {
            vdims1.SetSize(1);
            vdims1=py::extract<double>(dims1)();
          }
          else if (py::extract<py::tuple>(dims1).check())
          {
            py::tuple tdims1(dims1);
            vdims1.SetSize(py::len(tdims1));
            for (int j : Range(py::len(tdims1)))
              vdims1(j)=py::extract<double>(tdims1[j])();
          }
          else 
          {
            vdims1.SetSize(dim1);
            for (int j : Range(dim1))
              vdims1(j)=j+1;
          }
          if (py::extract<double>(dims2).check())
          {
            vdims2.SetSize(1);
            vdims2=py::extract<double>(dims2)();
          }
          else if (py::extract<py::tuple>(dims2).check())
          {
            py::tuple tdims2(dims2);
            vdims2.SetSize(py::len(tdims2));
            for (int j : Range(py::len(tdims2)))
              vdims2(j)=py::extract<double>(tdims2[j])();
          }
          else
          {
            vdims2.SetSize(dim2);
            for (int j : Range(dim2))
              vdims2(j)=j+dim1+1;
          }
          if (vdims1.Size()!=dim1 || vdims2.Size()!=dim2)
          {
            throw Exception("Dimensions do not match");
          }
          switch (dim)
          {
            case 1:
              if (dim1==1)
                return pml1;
              else
                return pml2;
            case 2:
              switch(dim1)
              {
                case 0:
                  return shared_ptr<PML>(make_shared<CompoundPML<2,0,2>> (pml1,pml2,vdims1,vdims2));
                case 1:
                  return shared_ptr<PML>(make_shared<CompoundPML<2,1,1>> (pml1,pml2,vdims1,vdims2));
                case 2:
                  return shared_ptr<PML>(make_shared<CompoundPML<2,2,0>> (pml1,pml2,vdims1,vdims2));
              }
            case 3:
              switch(dim1)
              {
                case 0:
                  return shared_ptr<PML>(make_shared<CompoundPML<3,0,3>> (pml1,pml2,vdims1,vdims2));
                case 1:
                  return shared_ptr<PML>(make_shared<CompoundPML<3,1,2>> (pml1,pml2,vdims1,vdims2));
                case 2:
                  return shared_ptr<PML>(make_shared<CompoundPML<3,2,1>> (pml1,pml2,vdims1,vdims2));
                case 3:
                  return shared_ptr<PML>(make_shared<CompoundPML<3,3,0>> (pml1,pml2,vdims1,vdims2));
              }
          }
          throw Exception("No valid dimension");
        },
        py::arg("pml1"),py::arg("pml2"), 
        py::arg("dims1")=DummyArgument(),py::arg("dims2")=DummyArgument(),
        R"raw_string(tensor product of two pml transformations

        dimensions are optional, given as tuples/lists and start with 1)raw_string")
      ;
}



void NGS_DLL_HEADER ExportNgcomp(py::module &m)
{

  py::module pml = m.def_submodule("pml", "module for perfectly matched layers");
  ExportPml(pml);
  //////////////////////////////////////////////////////////////////////////////////////////

  py::enum_<VorB>(m, "VorB", "Enum specifying the codimension. VOL is volume, BND is boundary and BBND is codimension 2 (edges in 3D, points in 2D)")
    .value("VOL", VOL)
    .value("BND", BND)
    .value("BBND", BBND)
    .value("BBBND", BBBND)
    .export_values()
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  py::enum_<COUPLING_TYPE> (m, "COUPLING_TYPE", docu_string(R"raw_string(
Enum specifying the coupling type of a degree of freedom, each dof is
either UNUSED_DOF, LOCAL_DOF, INTERFACE_DOF or WIREBASKET_DOF, other values
are provided as combinations of these:

UNUSED_DOF: Dof is not used, i.e the slave dofs in a :any:`Periodic` finite
    element space.

LOCAL_DOF: Inner degree of freedom, will be eliminated by static
    condensation.

INTERFACE_DOF: Degree of freedom between two elements, these will not be
    eliminated by static condensation, but not be put into the wirebasket
    system for i.e. a bddc :any:`Preconditioner`.

NONWIREBASKET_DOF: Either a LOCAL_DOF or an INTERFACE_DOF

WIREBASKET_DOF: Degree of freedom coupling with many elements (more than
    one). These will be put into the system for a bddc preconditioner.
    The HCurl space also treats degrees of freedom of badly shaped
    elements as WIREBASKET_DOFs.

EXTERNAL_DOF: Either INTERFACE_DOF or WIREBASKET_DOF

ANY_DOF: Any used dof (LOCAL_DOF or INTERFACE_DOF or WIREBASKET_DOF)

)raw_string"))
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

  py::class_<ElementId> (m, "ElementId", 
                         "an element identifier containing element number and Volume/Boundary flag")
    .def(py::init<VorB,size_t>())
    .def(py::init<size_t>())
    .def(py::init<Ngs_Element>())
    .def("__str__", &ToString<ElementId>)
    .def_property_readonly("nr", &ElementId::Nr, "the element number")    
    .def("VB", &ElementId::VB, "VorB of element")
    .def(py::self!=py::self)
    .def(py::self==py::self)
    .def("__hash__" , &ElementId::Nr)
    ;
  
  m.def("BndElementId",[] (int nr) { return ElementId(BND,nr); },
          py::arg("nr"),
          "creates an element-id for a boundary element")
    ;

  py::class_<NodeId> (m, "NodeId",
                      "an node identifier containing node type and node nr")
    .def(py::init<NODE_TYPE,size_t>())
    .def("__str__", &ToString<NodeId>)
    .def("__repr__", &ToString<NodeId>)
    /*
    .def("__repr__", [](NodeId & self)
         { return string("NodeId(")+ToString(self.GetType())+","+ToString(self.GetNr())+")"; })
    */
    .def(py::self!=py::self)
    .def(py::self==py::self)
    .def("__hash__" , &NodeId::GetNr)    
    .def_property_readonly("type", &NodeId::GetType, "the node type")        
    .def_property_readonly("nr", &NodeId::GetNr, "the node number")    
    ;

  class MeshNode : public NodeId
  {
    const MeshAccess & ma;
  public:
    MeshNode (NodeId _ni, const MeshAccess & _ma)
      : NodeId(_ni), ma(_ma) { ; }
    auto & Mesh() { return ma; }
    MeshNode operator++ (int) { return MeshNode(NodeId::operator++(0),ma); }
    MeshNode operator++ () { return MeshNode(NodeId::operator++(), ma); }
  };

  py::class_<MeshNode, NodeId> (m, "MeshNode", "a node within a mesh")
    .def_property_readonly("vertices", [](MeshNode & node) -> py::tuple
                           {
                             const MeshAccess & mesh = node.Mesh();
                             if (node.GetType() == NT_EDGE)
                               {
                                 /*
                                 auto verts = node.Mesh().GetEdgePNums(node.GetNr());
                                 py::tuple tup(2);
                                 for (size_t i = 0; i < 2; i++)
                                   tup[i] = py::cast(NodeId(NT_VERTEX, verts[i]));
                                 return tup;
                                 */
                                 
                                 auto Nr2MeshVert = [&mesh] (size_t nr)
                                   { return MeshNode(NodeId(NT_VERTEX,nr), mesh); };
                                 auto verts = node.Mesh().GetEdgePNums(node.GetNr());
                                 return MakePyTuple (Substitute(ArrayObject(verts), Nr2MeshVert));
                               }
                             if (node.GetType() == NT_FACE)
                               {
                                 /*
                                 auto verts = node.Mesh().GetFacePNums(node.GetNr());
                                 py::tuple tup(verts.Size());
                                 for (size_t i = 0; i < verts.Size(); i++)
                                   tup[i] = py::cast(NodeId(NT_VERTEX, verts[i]));
                                 return tup;
                                 */
                                 auto verts = node.Mesh().GetFacePNums(node.GetNr());
                                 return MakePyTuple (Substitute(verts, Nr2Vert));
                               }
                             throw py::type_error("vertices only available for edge and face nodes\n");
                           },
                           "tuple of global vertex numbers")
    
    ;

  /*
  py::class_<ngstd::T_Range<NodeId>> py_intrange (m, "NodeRange");
  py_intrange.def("__iter__", [] (ngstd::T_Range<NodeId> & r)
      { return py::make_iterator(r.begin(), r.end()); },
      py::keep_alive<0,1>()
    );
  */
  py::class_<ngstd::T_Range<NodeId>> (m, "NodeRange")
    .def("__len__", &T_Range<NodeId>::Size)
    .def("__iter__", [] (ngstd::T_Range<NodeId> & r)
         { return py::make_iterator(r.begin(), r.end()); },
         py::keep_alive<0,1>())
    ;

  py::class_<ngstd::T_Range<MeshNode>> (m, "MeshNodeRange")
    .def("__len__", &T_Range<MeshNode>::Size)
    .def("__iter__", [] (ngstd::T_Range<MeshNode> & r)
         { return py::make_iterator(r.begin(), r.end()); },
         py::keep_alive<0,1>())
    ;
  
  py::enum_<ORDER_POLICY>(m, "ORDER_POLICY")
    .value("CONSTANT", CONSTANT_ORDER)
    .value("NODETYPE", NODE_TYPE_ORDER)
    .value("VARIABLE", VARIABLE_ORDER)
    .value("OLDSTYLE", OLDSTYLE_ORDER)
    ;


  //////////////////////////////////////////////////////////////////////////////////////////


  py::class_<FlatArray<NodeId> > class_flatarrayNI (m, "FlatArrayNI");
  PyDefVector<FlatArray<NodeId>, NodeId>(m, class_flatarrayNI);
  PyDefToString<FlatArray<NodeId> >(m, class_flatarrayNI);
  class_flatarrayNI.def(py::init<size_t, NodeId *>());

  py::class_<Array<NodeId>, FlatArray<NodeId> >(m, "ArrayNI")
    .def(py::init<size_t>())
    ;

  
  // TODO: make tuple not doing the right thing
  py::class_<Ngs_Element>(m, "Ngs_Element")
    .def_property_readonly("nr", &Ngs_Element::Nr, "the element number")    
    .def("VB", &Ngs_Element::VB, "VorB of element")   
    .def_property_readonly("vertices", [](Ngs_Element &el)
                           {
                             return MakePyTuple(Substitute(el.Vertices(), Nr2Vert));
                           },
                           "tuple of global vertex numbers")
    .def_property_readonly("edges", [](Ngs_Element &el)
                           {
                             return MakePyTuple(Substitute(el.Edges(), Nr2Edge));
                           },
                           "tuple of global edge numbers")
    .def_property_readonly("faces", [](Ngs_Element &el)
                           {
                             return MakePyTuple(Substitute(el.Faces(), Nr2Face));
                           },
                           "tuple of global face numbers")
    .def_property_readonly("type", [](Ngs_Element &self)
        { return self.GetType(); },
        "geometric shape of element")
    .def_property_readonly("index", [](Ngs_Element &self)
        { return self.GetIndex(); },
        "material or boundary condition index")
    .def_property_readonly("mat", [](Ngs_Element & el)
                           { return el.GetMaterial(); },
                           "material or boundary condition label")
    ;

  py::implicitly_convertible <Ngs_Element, ElementId> ();
  // py::implicitly_convertible <ElementId, Ngs_Element> ();

  py::class_<FESpace::Element,Ngs_Element>(m, "FESpaceElement")
    .def_property_readonly("dofs",
                  [](FESpace::Element & el) 
                   {
                     py::list res;
                     Array<int> tmp (el.GetDofs());
                     for (int i : tmp)
                       res.append(py::cast(i));
                     return res;
                   },
                  "degrees of freedom of element"
                  )

    .def("GetLH",[](FESpace::Element & el) -> LocalHeap & 
                                  {
                                    return el.GetLH();
                                  },
         py::return_value_policy::reference
         )
    
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

  py::class_<Region> (m, "Region", "a subset of volume or boundary elements")
    .def(py::init<shared_ptr<MeshAccess>,VorB,string>())
    .def(py::init<shared_ptr<MeshAccess>,VorB,BitArray>())
    .def("Mask",[](Region & reg)->BitArray { return reg.Mask(); })
    .def("VB", [](Region & reg) { return VorB(reg); })
    .def(py::self + py::self)
    .def(py::self + string())
    .def(py::self - py::self)
    .def(py::self - string())
    .def(~py::self)
    ;

  py::implicitly_convertible <Region, BitArray> ();


  //////////////////////////////////////////////////////////////////////////////////////////
  
  
  typedef PML_Transformation PML;
  py::class_<MeshAccess, shared_ptr<MeshAccess>>(m, "Mesh", docu_string(R"raw_string(
NGSolve interface to the Netgen mesh. Provides access and functionality
to use the mesh for finite element calculations.

Parameters

mesh (netgen.Mesh): a mesh generated from Netgen


)raw_string") , py::dynamic_attr())
    .def(py::init<shared_ptr<netgen::Mesh>>())
    
#ifndef PARALLEL
    .def(py::init([](const string & filename)
                  { return make_shared<MeshAccess>(filename); }),
          py::arg("filename"))
#else
    .def(py::init([](const string & filename, PyMPI_Comm c)
                  {
                    ngs_comm = c.comm;
                    NGSOStream::SetGlobalActive (MyMPI_GetId(c.comm)==0);                      
                    return make_shared<MeshAccess>(filename, c.comm);
                  }),
         py::arg("filename"), py::arg("comm")=PyMPI_Comm(MPI_COMM_WORLD))
#endif

    
    .def("__eq__",
         [] (shared_ptr<MeshAccess> self, shared_ptr<MeshAccess> other)
           {
             return self == other;
           })
    .def(py::pickle([](const MeshAccess& ma)
                    {
                      stringstream ss;
                      ma.SaveMesh(ss);
                      return py::make_tuple(ss.str());
                    },
                    [] (py::tuple state)
                    {
                      auto ma = make_shared<MeshAccess>();
                      stringstream ss(state[0].cast<string>());
                      ma->LoadMesh(ss);
                      return ma;
                    }))
    .def("LoadMesh", static_cast<void(MeshAccess::*)(const string &)>(&MeshAccess::LoadMesh),
         "Load mesh from file")
    
    .def("Elements", static_cast<ElementRange(MeshAccess::*)(VorB)const> (&MeshAccess::Elements),
	 (py::arg("VOL_or_BND")=VOL), docu_string("Returns an iterator over ElementIds on VorB"))

    .def("__getitem__", static_cast<Ngs_Element(MeshAccess::*)(ElementId)const> (&MeshAccess::operator[]))
    .def("__getitem__", [](MeshAccess & self, NodeId ni) { return MeshNode(ni, self); })

    .def ("GetNE", static_cast<size_t(MeshAccess::*)(VorB)const> (&MeshAccess::GetNE), docu_string("Number of elements of codimension VorB."))
    .def_property_readonly ("nv", &MeshAccess::GetNV, "Number of vertices")
    .def_property_readonly ("ne",  static_cast<size_t(MeshAccess::*)()const> (&MeshAccess::GetNE), "Number of volume elements")
    .def_property_readonly ("nedge", &MeshAccess::GetNEdges, "Number of edges")
    .def_property_readonly ("nface", &MeshAccess::GetNFaces, "Number of faces")    
    .def_property_readonly ("dim", &MeshAccess::GetDimension, "Mesh dimension")
    .def_property_readonly ("ngmesh", &MeshAccess::GetNetgenMesh, "Get the Netgen mesh")

    .def_property_readonly ("vertices", [] (shared_ptr<MeshAccess> mesh)
          {
            // return mesh->Nodes(NT_VERTEX);
            return T_Range<MeshNode> (MeshNode(NodeId(NT_VERTEX, 0), *mesh),
                                      MeshNode(NodeId(NT_VERTEX, mesh->GetNNodes(NT_VERTEX)), *mesh));
          }, "iterable of mesh vertices")
    .def_property_readonly ("edges", [] (shared_ptr<MeshAccess> mesh)
          {
            return T_Range<MeshNode> (MeshNode(NodeId(NT_EDGE, 0), *mesh),
                                      MeshNode(NodeId(NT_EDGE, mesh->GetNNodes(NT_EDGE)), *mesh));
            /*
            return T_Range<NodeId> (NodeId(NT_EDGE, 0),
                                    NodeId(NT_EDGE, mesh->GetNNodes(NT_EDGE)));
            */
            // return mesh->Nodes(NT_EDGE);
          }, "iterable of mesh edges")
    .def_property_readonly ("faces", [] (shared_ptr<MeshAccess> mesh)
          {
            // return mesh->Nodes(NT_FACE);
            return T_Range<MeshNode> (MeshNode(NodeId(NT_FACE, 0), *mesh),
                                      MeshNode(NodeId(NT_FACE, mesh->GetNNodes(NT_FACE)), *mesh));
          }, "iterable of mesh faces")
    
    .def ("GetTrafo", 
          static_cast<ElementTransformation&(MeshAccess::*)(ElementId,Allocator&)const>
          (&MeshAccess::GetTrafo), 
          py::return_value_policy::reference)

    .def ("GetTrafo",
          [](MeshAccess & ma, ElementId id)
          { return &ma.GetTrafo(id, global_alloc); },
          py::return_value_policy::take_ownership)

    .def("SetDeformation", 
	 [](MeshAccess & ma, shared_ptr<GF> gf)
         { ma.SetDeformation(gf); },
         docu_string("Deform the mesh with the given GridFunction"))

    .def("SetPML", 
	 [](MeshAccess & ma,  shared_ptr<PML> apml, py::object definedon)
          {
            if (py::extract<int>(definedon).check())
              {
                ma.SetPML(apml, py::extract<int>(definedon)()-1);
              }

            if (py::isinstance<py::str>(definedon))
              {
                std::regex pattern(definedon.cast<string>());
                for (int i = 0; i < ma.GetNDomains(); i++)
                  if (std::regex_match (ma.GetMaterial(VOL,i), pattern))
                    ma.SetPML(apml, i);
              }
          },
         py::arg("pmltrafo"),py::arg("definedon"),
         "set PML transformation on domain"
         )
    .def("UnSetPML", [](MeshAccess & ma, py::object definedon)
          {
            if (py::extract<int>(definedon).check())
                ma.UnSetPML(py::extract<int>(definedon)()-1);

            if (py::isinstance<py::str>(definedon))
              {
                std::regex pattern(definedon.cast<string>());
                for (int i = 0; i < ma.GetNDomains(); i++)
                  if (std::regex_match (ma.GetMaterial(VOL,i), pattern))
                    ma.UnSetPML(i);
              }
          })
    
    .def("GetPMLTrafos", [](MeshAccess & ma) 
      {
        py::list pml_trafos(ma.GetNDomains());
        for (int i : Range(ma.GetNDomains()))
        {
          if (ma.GetPMLTrafos()[i])
            pml_trafos[i] = shared_ptr<PML>(ma.GetPMLTrafos()[i]);
          else
            pml_trafos[i] = py::none();
        }
        return pml_trafos;
      },
        "returns list of pml transformations"
    )
    .def("GetPMLTrafo", [](MeshAccess & ma, int domnr) {
        if (ma.GetPMLTrafos()[domnr])
     	  return ma.GetPMLTrafos()[domnr-1];
        else
          throw Exception("No PML Trafo set"); 
        },
        py::arg("dom")=1,
        "returns pml transformation on domain dom"
        )

    .def("UnsetDeformation", [](MeshAccess & ma){ ma.SetDeformation(nullptr);})
    
    .def("GetMaterials",
	 [](const MeshAccess & ma)
	  {
            return MakePyTuple(ma.GetMaterials(VOL));
            // return py::cast(ma.GetMaterials(VOL));
	  },
	 "Returns list of materials"
         )

    .def("Materials",
	 [](shared_ptr<MeshAccess> ma, string pattern) 
	  {
            return Region (ma, VOL, pattern);
	  },
         py::arg("pattern"),
	 "Returns mesh-region matching the given regex pattern"
         )
    .def("Materials",
	 [](shared_ptr<MeshAccess> ma, vector<int> domains)
	  {
            BitArray mask(ma->GetNDomains());
            mask.Clear();
            for (auto i : domains)
              if (i >= 0 && i < mask.Size())
                mask.Set(i);
              else
                throw Exception ("index "+ToString(i)+" out of range [0,"+ToString(mask.Size())+")");
            return Region (ma, VOL, mask);
	  },
         py::arg("domains"),
	 "Generate mesh-region by domain numbers"
         )
    
    
    .def("GetBoundaries",
	 [](const MeshAccess & ma)
	  {
            return MakePyTuple(ma.GetMaterials(BND));
	  },
	 "Returns list of boundary conditions"
         )

    .def("Boundaries",
	 [](shared_ptr<MeshAccess> ma, string pattern)
	  {
            return Region (ma, BND, pattern);
	  },
         py::arg("pattern"),
	 "Returns boundary mesh-region matching the given regex pattern"
         )
    .def("Boundaries",
	 [](shared_ptr<MeshAccess> ma, vector<int> bnds)
	  {
            BitArray mask(ma->GetNBoundaries());
            mask.Clear();
            for (auto i : bnds)
              if (i >= 0 && i < mask.Size())
                mask.Set(i);
              else
                throw Exception ("boundary index "+ToString(i)+" out of range [0,"+ToString(mask.Size())+")");
            return Region (ma, BND, mask);
	  },
         py::arg("bnds"),
	 "Generate boundary mesh-region by boundary condition numbers"
         )

    .def("GetBBoundaries",
	 [](const MeshAccess & ma)
	  {
            return MakePyTuple(ma.GetMaterials(BBND));
	  },
	 "Returns list of boundary conditions for co dimension 2"
	 )
    .def("BBoundaries", [](shared_ptr<MeshAccess> ma, string pattern)
	  {
	    return Region (ma, BBND, pattern);
	  },
	 (py::arg("self"), py::arg("pattern")),
	 "Returns co dim 2 boundary mesh-region matching the given regex pattern"
	 )

    // TODO: explain how to mark elements
    .def("Refine",
         [](MeshAccess & ma)
          {
            ma.Refine();
          },
	 "Local mesh refinement based on marked elements, uses element-bisection algorithm")

    // TODO: explain how to mark vertices and edges, explain how factor is used
    .def("RefineHP",
         [](MeshAccess & ma, int levels, double factor)
          {
            Ng_HPRefinement(levels, factor);
            ma.UpdateBuffers();
          },
         py::arg("levels"), py::arg("factor")=0.125,
	 "Geometric mesh refinement towards marked vertices and edges, uses factor for placement of new points"
         )

    // TODO: Docu string says nothing... what does refinement flag do?
    .def("SetRefinementFlag", &MeshAccess::SetRefinementFlag,
	 "Set refinementflag for mesh-refinement")

    // TODO: Docu
    .def("GetParentElement", &MeshAccess::GetParentElement)
    // TODO: Docu
    .def("GetParentVertices", [](MeshAccess & ma, int vnum)
          {
            Array<int> parents(2);
            ma.GetParentNodes (vnum, &parents[0]);
            return py::make_tuple(parents[0], parents[1]);
          })

    .def("SetElementOrder",
         [](MeshAccess & ma, ElementId id, int order)
         {
           ma.SetElOrder(id.Nr(), order);
         })
    
    // TODO: Docu
    .def("Curve",
         [](MeshAccess & ma, int order)
          {
            Ng_HighOrder(order);
          },
         py::arg("order"))

    .def("__call__",
         [](MeshAccess & ma, double x, double y, double z, VorB vb) 
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
          , 
         py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0,
         py::arg("VOL_or_BND") = VOL,
	 docu_string("Get a MappedIntegrationPoint in the point (x,y,z) on the matching volume (VorB=VOL, default) or surface (VorB=BND) element. BBND elements aren't supported"))

    .def("Contains",
         [](MeshAccess & ma, double x, double y, double z) 
          {
            IntegrationPoint ip;
            int elnr = ma.FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            return (elnr >= 0);
          }, 
         py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0
	 ,"Checks if the point (x,y,z) is in the meshed domain (is inside a volume element)")

    ;

  //////////////////////////////////////////////////////////////////////////////////////////
  
  py::class_<NGS_Object, shared_ptr<NGS_Object>>(m, "NGS_Object")
    .def_property_readonly("name", [](const NGS_Object & self)->string { return self.GetName();})
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


  static size_t global_heapsize = 1000000;
  static LocalHeap glh(global_heapsize, "python-comp lh", true);
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

  auto fesPickle = [](const FESpace& fes)
    {
      auto flags = fes.GetFlags();
      auto mesh = fes.GetMeshAccess();
      auto type = fes.type;
      // TODO: pickle order policies
      return py::make_tuple(type,mesh,flags);
    };

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
                    shared_ptr<FESpace> fes = make_shared<CompoundFESpace> (spaces[0]->GetMeshAccess(), spaces, flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                    //                              py::cast(*instance).attr("flags") = bpflags;
                  }),
                  py::arg("spaces"),
                  "construct compound-FESpace from list of component spaces"
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
    // .def("__str__", &ToString<FESpace>)
    .def("__str__", [] (shared_ptr<FESpace> self) { return ToString(*self); } )
    .def("__timing__", [] (shared_ptr<FESpace> self) { return py::cast(self->Timing()); })

    .def_property_readonly("mesh",
                           [](shared_ptr<FESpace> self) -> shared_ptr<MeshAccess>
                           { return self->GetMeshAccess(); })

    .def_property_readonly("order", [] (shared_ptr<FESpace> self) { return OrderProxy(*self); },
                  "proxy to set order for individual nodes")
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
         [](shared_ptr<FESpace> self, ELEMENT_TYPE et, py::object order, py::object order_left, py::object order_right)
         {
           if (py::isinstance<py::int_> (order))
             {
               self->SetOrderLeft (et, order.cast<py::int_>());
               self->SetOrderRight (et, order.cast<py::int_>());
             }
           if (py::isinstance<py::int_> (order_left))
             self->SetOrderLeft (et, order_left.cast<py::int_>());
           if (py::isinstance<py::int_> (order_right))
             self->SetOrderRight (et, order_right.cast<int>());
         },
         py::arg("element_type"),
         py::arg("order")=DummyArgument(),
         py::arg("order_left")=DummyArgument(),
         py::arg("order_right")=DummyArgument()
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
           Array<int> tmp; self->GetDofNrs(ei,tmp);
           return MakePyTuple(tmp);           
         })

    .def("GetDofNrs", [](shared_ptr<FESpace> self, NodeId ni)
         {
           Array<int> tmp; self->GetDofNrs(ni,tmp);
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
            
            return py::cast(fe);
          })
          
    .def ("GetFE", [](shared_ptr<FESpace> self, ElementId ei, LocalHeap & lh)
          {
            return shared_ptr<FiniteElement>(&self->GetFE(ei, lh), NOOP_Deleter);
          },
          py::return_value_policy::reference)
    
    .def("FreeDofs",
         [] (const shared_ptr<FESpace>self, bool coupling)
         { return self->GetFreeDofs(coupling); },
         py::arg("coupling")=false)

    .def("ParallelDofs",
         [] (const shared_ptr<FESpace>self)
         { return self->GetParallelDofs(); })

    .def("Range",
         [] (const shared_ptr<FESpace> self, int comp) -> py::slice
         {
           auto compspace = dynamic_pointer_cast<CompoundFESpace> (self);
           if (!compspace)
             throw py::type_error("'Range' is available only for product spaces");
           IntRange r = compspace->GetRange(comp);
           return py::slice(py::int_(r.First()), py::int_(r.Next()),1);
         })

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
                  "list of gridfunctions for compound gridfunction")

    .def("TrialFunction",
         [] (const shared_ptr<FESpace> self)
         {
           return MakeProxyFunction (self, false);
         },
         docu_string("Gives a proxy to be used as a trialfunction in :any:`Symbolic Integrators<symbolic-integrators>`"))
    
    .def("TestFunction",
         [] (const shared_ptr<FESpace> self)
           {
             return MakeProxyFunction (self, true);
           },
         docu_string("Gives a proxy to be used as a testfunction for :any:`Symbolic Integrators<symbolic-integrators>`"))

    .def("TnT",
         [] (const shared_ptr<FESpace> self)
         {
           return std::make_tuple(MakeProxyFunction (self, false), MakeProxyFunction (self, true));
         },
         docu_string("Gives a tuple of trial and testfunction"))

    .def("SolveM",
         [] (const shared_ptr<FESpace> self,
             spCF rho, BaseVector& vec, int heapsize)
          {
            if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
                bool first_time = true;
                if (first_time)
                  { first_time = false; cerr << "warning: use SetHeapSize(size) instead of heapsize=size" << endl; }
              }
            self->SolveM(*rho, vec, glh);
          },
         py::arg("rho"), py::arg("vec"), py::arg("heapsize")=1000000)
        
    .def("__eq__",
         [] (shared_ptr<FESpace> self, shared_ptr<FESpace> other)
         {
           return self == other;
         })
    ;

  auto hcurl = py::class_<HCurlHighOrderFESpace, shared_ptr<HCurlHighOrderFESpace>,FESpace>
    (m, "HCurl");
  hcurl
    .def(py::init([hcurl](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(hcurl, kwargs, info);
                    auto fes = make_shared<HCurlHighOrderFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle, (shared_ptr<HCurlHighOrderFESpace>(*)(py::tuple))
                    fesUnpickle<HCurlHighOrderFESpace>))
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

  py::class_<CompoundFESpace, shared_ptr<CompoundFESpace>, FESpace>
    (m,"CompoundFESpace")
    .def(py::pickle([] (py::object pyfes)
                    {
                      auto fes = py::cast<shared_ptr<CompoundFESpace>>(pyfes);
                      auto flags = fes->GetFlags();
                      py::list lst;
                      for(auto i : Range(fes->GetNSpaces()))
                        {
                          lst.append((*fes)[i]);
                        }
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
  
  auto hdiv = py::class_<HDivHighOrderFESpace, shared_ptr<HDivHighOrderFESpace>,FESpace>
    (m, "HDiv");
  hdiv
    .def(py::init([hdiv](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(hdiv, kwargs, info);
                    auto fes = make_shared<HDivHighOrderFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<HDivHighOrderFESpace>(*)(py::tuple))
                    fesUnpickle<HDivHighOrderFESpace>))
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
    .def("Average",
         [] (shared_ptr<HDivHighOrderFESpace> hdivfes, BaseVector & bv)
         {
           auto & pairs = hdivfes->GetDCPairs();
           auto fu = bv.FV<double>();
           for (auto pair : pairs)
             {
               auto f1 = pair[0];
               auto f2 = pair[1];
               if (f2 != -1)
                 {
                   double mean = 0.5 * (fu(f1) + fu(f2));
                   fu(f1) = fu(f2) = mean;
                 }
               else if (f1 != -1)
                 fu(f1) = 0.0;
             }
         },
         py::arg("vector"))
    ;


  auto h1 = py::class_<H1HighOrderFESpace, shared_ptr<H1HighOrderFESpace>,FESpace>
    (m, "H1");
  h1
    .def(py::init([h1](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(h1, kwargs, info);
                    auto fes = make_shared<H1HighOrderFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<H1HighOrderFESpace>(*)(py::tuple))
                    fesUnpickle<H1HighOrderFESpace>))
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  return flags_doc;
                })
    ;

  auto vectorh1 = py::class_<VectorH1FESpace, shared_ptr<VectorH1FESpace>,FESpace>
    (m, "VectorH1");
  vectorh1
    .def(py::init([vectorh1](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(vectorh1, kwargs, info);
                    auto fes = make_shared<VectorH1FESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<VectorH1FESpace>(*)(py::tuple))
                    fesUnpickle<VectorH1FESpace>))
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  return flags_doc;
                })
    ;

  auto l2 = py::class_<L2HighOrderFESpace, shared_ptr<L2HighOrderFESpace>,FESpace>
    (m, "L2");
  l2
    .def(py::init([l2](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(l2, kwargs, info);
                    auto fes = make_shared<L2HighOrderFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<L2HighOrderFESpace>(*)(py::tuple))
                    fesUnpickle<L2HighOrderFESpace>))
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  return flags_doc;
                })
    ;

  auto vectorl2 = py::class_<VectorL2FESpace, shared_ptr<VectorL2FESpace>,FESpace>
    (m, "VectorL2");
  vectorl2
    .def(py::init([vectorl2](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(vectorl2, kwargs, info);
                    auto fes = make_shared<VectorL2FESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<VectorL2FESpace>(*)(py::tuple))
                    fesUnpickle<VectorL2FESpace>))
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  return flags_doc;
                })
    ;

  auto l2surface = py::class_<L2SurfaceHighOrderFESpace, shared_ptr<L2SurfaceHighOrderFESpace>,FESpace>
    (m, "SurfaceL2");
  l2surface
    .def(py::init([l2surface](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(l2surface, kwargs, info);
                    auto fes = make_shared<L2SurfaceHighOrderFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<L2SurfaceHighOrderFESpace>(*)(py::tuple))
                    fesUnpickle<L2SurfaceHighOrderFESpace>))
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  return flags_doc;
                })
    ;

  auto numberfes = py::class_<NumberFESpace, shared_ptr<NumberFESpace>,FESpace>
    (m, "NumberSpace");
  numberfes
    .def(py::init([numberfes](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(numberfes, kwargs, info);
                    auto fes = make_shared<NumberFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<NumberFESpace>(*)(py::tuple))
                    fesUnpickle<NumberFESpace>))
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  return flags_doc;
                })
    ;

  auto hdivdiv = py::class_<HDivDivFESpace, shared_ptr<HDivDivFESpace>,FESpace>
    (m, "HDivDiv");
  hdivdiv
    .def(py::init([hdivdiv](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(hdivdiv, kwargs, info);
                    auto fes = make_shared<HDivDivFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<HDivDivFESpace>(*)(py::tuple))
                    fesUnpickle<HDivDivFESpace>))
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


  auto hdivdivsurf = py::class_<HDivDivSurfaceSpace, shared_ptr<HDivDivSurfaceSpace>,FESpace>
    (m, "HDivDivSurface");
  hdivdivsurf
    .def(py::init([hdivdivsurf](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(hdivdivsurf, kwargs, info);
                    auto fes = make_shared<HDivDivSurfaceSpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<HDivDivSurfaceSpace>(*)(py::tuple))
                    fesUnpickle<HDivDivSurfaceSpace>))
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

  auto vectorfacet = py::class_<VectorFacetFESpace, shared_ptr<VectorFacetFESpace>,FESpace>
    (m, "VectorFacet");
  vectorfacet
    .def(py::init([vectorfacet](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(vectorfacet, kwargs, info);
                    auto fes = make_shared<VectorFacetFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<VectorFacetFESpace>(*)(py::tuple))
                    fesUnpickle<VectorFacetFESpace>))
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  return flags_doc;
                })
    ;

  auto facetfes = py::class_<FacetFESpace, shared_ptr<FacetFESpace>,FESpace>
    (m, "FacetFESpace");
  facetfes
    .def(py::init([facetfes](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(facetfes, kwargs, info);
                    auto fes = make_shared<FacetFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<FacetFESpace>(*)(py::tuple))
                    fesUnpickle<FacetFESpace>))
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  return flags_doc;
                })
    ;

  auto facetsurface = py::class_<FacetSurfaceFESpace, shared_ptr<FacetSurfaceFESpace>,FESpace>
    (m, "FacetSurface");
  facetsurface
    .def(py::init([facetsurface](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(facetsurface, kwargs, info);
                    auto fes = make_shared<FacetSurfaceFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<FacetSurfaceFESpace>(*)(py::tuple))
                    fesUnpickle<FacetSurfaceFESpace>))
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  return flags_doc;
                })
    ;
   
  auto hdivsurf = py::class_<HDivHighOrderSurfaceFESpace, shared_ptr<HDivHighOrderSurfaceFESpace>,FESpace>
    (m, "HDivSurface");
  hdivsurf
    .def(py::init([hdivsurf](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
                  {
                    py::list info;
                    info.append(ma);
                    auto flags = CreateFlagsFromKwArgs(hdivsurf, kwargs, info);
                    auto fes = make_shared<HDivHighOrderSurfaceFESpace>(ma,flags);
                    fes->Update(glh);
                    fes->FinalizeUpdate(glh);
                    return fes;
                  }),py::arg("mesh"))
    .def(py::pickle(fesPickle,(shared_ptr<HDivHighOrderSurfaceFESpace>(*)(py::tuple))
                    fesUnpickle<HDivHighOrderSurfaceFESpace>))
    .def_static("__flags_doc__", [] ()
                {
                  auto flags_doc = py::cast<py::dict>(py::module::import("ngsolve").
                                                  attr("FESpace").
                                                  attr("__flags_doc__")());
                  flags_doc["discontinuous"] = "bool = False\n"
                    "  Create discontinuous HDivSurface space";
                  return flags_doc;
                })
    .def("Average",
         [] (shared_ptr<HDivHighOrderSurfaceFESpace> hdivsurffes, BaseVector & bv)
         {
           auto & pairs = hdivsurffes->GetDCPairs();
           auto fu = bv.FV<double>();
           for (auto pair : pairs)
             {
               auto f1 = pair[0];
               auto f2 = pair[1];
               if (f2 != -1)
                 {
                   double mean = 0.5 * (fu(f1) + fu(f2));
                   fu(f1) = fu(f2) = mean;
                 }
               else if (f1 != -1)
                 fu(f1) = 0.0;
             }
         },
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
                      auto fes = gf.GetFESpace();
                      auto flags = gf.GetFlags();
                      auto name = gf.GetName();
                      py::list lst;
                      for(auto val : gf.GetVector().FVDouble())
                        lst.append(val);
                      return py::make_tuple(fes,name,flags, lst);
                    },
                    [] (py::tuple state)
                    {
                      auto gf = CreateGridFunction(state[0].cast<shared_ptr<FESpace>>(),
                                                   state[1].cast<string>(),
                                                   state[2].cast<Flags>());
                      gf->Update();
                      auto lst = state[3].cast<py::list>();
                      for(auto i : Range(py::len(lst)))
                        gf->GetVector().FVDouble()[i] = lst[i].cast<double>();
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
            VorB boundary, py::object definedon, int heapsize, py::object heap)
         {
           shared_ptr<TPHighOrderFESpace> tpspace = dynamic_pointer_cast<TPHighOrderFESpace>(self->GetFESpace());          
            Region * reg = nullptr;
            if (py::extract<Region&> (definedon).check())
              reg = &py::extract<Region&>(definedon)();
            
            if (py::extract<LocalHeap&> (heap).check())
              {
                LocalHeap & lh = py::extract<LocalHeap&> (heap)();
                if(tpspace)
                {
                  Transfer2TPMesh(cf.get(),self.get(),lh);
                  return;
                }
                if (reg)
                  SetValues (cf, *self, *reg, NULL, lh);
                else
                  SetValues (cf, *self, boundary, NULL, lh);
                return;
              }

            if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
                bool first_time = true;
                if (first_time)
                  { first_time = false; cerr << "warning: use SetHeapSize(size) instead of heapsize=size" << endl; }
              }
            // LocalHeap lh(heapsize, "GridFunction::Set-lh", true);
            if(tpspace)
            {
              Transfer2TPMesh(cf.get(),self.get(),glh);
              return;
            }            
            if (reg)
              SetValues (cf, *self, *reg, NULL, glh);
            else
              SetValues (cf, *self, boundary, NULL, glh);
         },
          py::arg("coefficient"),
          py::arg("VOL_or_BND")=VOL,
          py::arg("definedon")=DummyArgument(),
          py::arg("heapsize")=1000000, py::arg("heap")=DummyArgument(),
         "Set values"
      )


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

    /*
    .def("CF", [](shared_ptr<GF> self) -> shared_ptr<CoefficientFunction>
          {
            return make_shared<GridFunctionCoefficientFunction> (self);
          })

    .def("CF", [](shared_ptr<GF> self, shared_ptr<DifferentialOperator> diffop)
          -> shared_ptr<CoefficientFunction>
          {
            return make_shared<GridFunctionCoefficientFunction> (self, diffop);
          })
    */
    .def("Deriv",
         [](shared_ptr<GF> self) -> spCF
          {
            auto sp = make_shared<GridFunctionCoefficientFunction> (self,
                                                                    self->GetFESpace()->GetFluxEvaluator(),
                                                                    self->GetFESpace()->GetFluxEvaluator(BND));
            // sp->SetDimensions(sp->Dimensions());
            return sp;
          })

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
                return py::cast(shared_ptr<CoefficientFunction>(coef));
              }
            return py::none(); //  shared_ptr<CoefficientFunction>();
          }, py::arg("name"), py::arg("VB")=VOL)

    
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
            // auto evaluator = space->GetEvaluator(ei.IsBoundary());
            auto evaluator = space->GetEvaluator(VorB(ei));

            // int elnr = mip.GetTransformation().GetElementNr();
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


    .def("CF", 
         [](shared_ptr<GF> self, shared_ptr<DifferentialOperator> diffop) -> spCF
          {
            if (!diffop->Boundary())
              return make_shared<GridFunctionCoefficientFunction> (self, diffop);
            else
              return make_shared<GridFunctionCoefficientFunction> (self, nullptr, diffop);
          })
    ;



  //////////////////////////////////////////////////////////////////////////////////////////

//   PyExportArray<shared_ptr<BilinearFormIntegrator>> ();

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

name : string
  Give a name to the BilinearForm - not really a use case in Python

check_unused : bool
  If set prints warnings if not UNUSED_DOFS are not used

)raw_string"));
  bf_class
    .def(py::init([bf_class] (shared_ptr<FESpace> fespace, bool check_unused,
                      string name, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(bf_class,kwargs);
                    auto biform = CreateBilinearForm (fespace, name, flags);
                    biform -> SetCheckUnused (check_unused);
                    return biform;
                  }),
         py::arg("space"),
         py::arg("check_unused")=true,
         py::arg("name")="bfa")
    .def(py::init([bf_class](shared_ptr<FESpace> trial_space,
                     shared_ptr<FESpace> test_space, string name,
                     bool check_unused, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(bf_class,kwargs);
                    auto biform = CreateBilinearForm (trial_space, test_space, name, flags);
                    biform -> SetCheckUnused (check_unused);
                    return biform;
                  }),
         py::arg("trialspace"),
         py::arg("testspace"),
         py::arg("name")="bfa", py::arg("check_unused")=true)

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
		     " The full matrix is stored, even if the symmetric flag is set."
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
                   {
                     py::list igts;
                     for (auto igt : self.Integrators())
                       igts.append(igt);
                     return igts;
                   } )
    
    .def("Assemble", [](BF & self, int heapsize, bool reallocate)
                                     {
                                       if (heapsize > global_heapsize)
                                         {
                                           global_heapsize = heapsize;
                                           glh = LocalHeap(heapsize, "python-comp lh", true);
                                           bool first_time = true;
                                           if (first_time)
                                             { first_time = false; cerr << "warning: use SetHeapSize(size) instead of heapsize=size" << endl; }                                           
                                         }
                                       self.ReAssemble(glh,reallocate);
                                     }, py::call_guard<py::gil_scoped_release>(),
         py::arg("heapsize")=1000000,py::arg("reallocate")=false)

    .def_property_readonly("mat", [](BF & self)
                                         {
                                           auto mat = self.GetMatrixPtr();
                                           if (!mat)
                                             throw py::type_error("matrix not ready - assemble bilinearform first");
                                           return mat;
                                         })

    .def("__getitem__",  [](BF & self, py::tuple t)
                                         {
                                           int ind1 = py::extract<int>(t[0])();
                                           int ind2 = py::extract<int>(t[1])();
                                           cout << "get bf, ind = " << ind1 << "," << ind2 << endl;
                                         })
    

    .def_property_readonly("components", [](shared_ptr<BilinearForm> self)-> py::list
                   { 
                     py::list bfs;
                     auto fes = dynamic_pointer_cast<CompoundFESpace> (self->GetFESpace());
                     if (!fes)
                       throw py::type_error("not a compound-fespace\n");
                       
                     int ncomp = fes->GetNSpaces();
                     for (int i = 0; i < ncomp; i++)
                       // bfs.append(shared_ptr<BilinearForm> (new ComponentBilinearForm(self.get(), i, ncomp)));
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
            return self.Energy(*x);
          }, py::call_guard<py::gil_scoped_release>())
    
    .def("Apply", [](BF & self, BaseVector& x, BaseVector & y, int heapsize)
	  {
            if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
                bool first_time = true;
                if (first_time)
                  { first_time = false; cerr << "warning: use SetHeapSize(size) instead of heapsize=size" << endl; }                
              }
	    self.ApplyMatrix (x, y, glh);
	  }, py::call_guard<py::gil_scoped_release>(),
         py::arg("x"),py::arg("y"),py::arg("heapsize")=1000000,docu_string(R"raw_string(
Applies a (non-)linear variational formulation to x and stores the result in y.

Parameters

x : ngsolve.BaseVector
  input vector

y : ngsolve.BaseVector
  output vector

heapsize : int
  Size of the LocalHeap for that operation. If you get an error about not
  enough heapsize increase this value.

)raw_string"))

    .def("ComputeInternal", [](BF & self, BaseVector & u, BaseVector & f, int heapsize)
	  {
            if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
                bool first_time = true;
                if (first_time)
                  { first_time = false; cerr << "warning: use SetHeapSize(size) instead of heapsize=size" << endl; }
                
              }
	    self.ComputeInternal (u, f, glh );
	  }, py::call_guard<py::gil_scoped_release>(),
         py::arg("u"),py::arg("f"),py::arg("heapsize")=1000000)

    .def("AssembleLinearization", [](BF & self, BaseVector & ulin, int heapsize)
	  {
            if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
                bool first_time = true;
                if (first_time)
                  { first_time = false; cerr << "warning: use SetHeapSize(size) instead of heapsize=size" << endl; }
              }
	    self.AssembleLinearization (ulin, glh);
	  }, py::call_guard<py::gil_scoped_release>(),
         py::arg("ulin"),py::arg("heapsize")=1000000)

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

  //////////////////////////////////////////////////////////////////////////////////////////

//   PyExportArray<shared_ptr<LinearFormIntegrator>> ();
//

  typedef LinearForm LF;
  auto lf_class = py::class_<LF, shared_ptr<LF>, NGS_Object>(m, "LinearForm", docu_string(R"raw_string(
Used to store the left hand side of a PDE. Add integrators
(ngsolve.LFI) to it to implement your PDE.

Parameters

space : ngsolve.FESpace
  The space the linearform is defined on. Can be a compound
  FESpace for a mixed formulation.

name : string
  The name of the linearform (in python not really in use...)

flags : dict
  Additional options for the linearform, for example:

    print : bool
      Write additional debug information to testout file. This
      file must be set by ngsolve.SetTestoutFile. Use
      ngsolve.SetNumThreads(1) for serial output.

)raw_string"));
  lf_class
    .def(py::init([lf_class] (shared_ptr<FESpace> fespace, string name, py::kwargs kwargs)
                  {
                    auto flags = CreateFlagsFromKwArgs(lf_class,kwargs);
                    auto f = CreateLinearForm (fespace, name, flags);
                    f->AllocateVector();
                    return f;
                  }),
         py::arg("space"), py::arg("name")="lff")
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

    .def("__iadd__",[](shared_ptr<LF> self, shared_ptr<LinearFormIntegrator> other) { (*self)+=other; return self; })

    .def_property_readonly("integrators",  [](shared_ptr<LF> self)
                   {
                     py::list igts;
                     for (auto igt : self->Integrators())
                       igts.append (igt);
                     return igts;
                   })

    .def("Assemble",  [](shared_ptr<LF> self, int heapsize)
         {
           if (heapsize > global_heapsize)
             {
               global_heapsize = heapsize;
               glh = LocalHeap(heapsize, "python-comp lh", true);
               bool first_time = true;
               if (first_time)
                 { first_time = false; cerr << "warning: use SetHeapSize(size) instead of heapsize=size" << endl; }                
             }
           self->Assemble(glh);
         }, py::call_guard<py::gil_scoped_release>(), py::arg("heapsize")=1000000)

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

  //////////////////////////////////////////////////////////////////////////////////////////

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
    .def("Do", [](NumProc & self, int heapsize)
                               {
                                 LocalHeap lh (heapsize, "NumProc::Do-heap");
                                 self.Do(lh);
                               }, py::call_guard<py::gil_scoped_release>(),
         py::arg("heapsize")=1000000)
    ;

  py::class_<PyNumProc, NumProc, shared_ptr<PyNumProc>> (m, "PyNumProc")
    .def("__init__",
         [](NumProc *instance, shared_ptr<PDE> pde, Flags & flags)
                           {
                             new (instance) PyNumProc(pde, flags);
                           })
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
	   bool region_wise, bool element_wise, int heapsize)
                          {
                            static Timer t("Integrate CF"); RegionTimer reg(t);
                            // static mutex addcomplex_mutex;
			    if (heapsize > global_heapsize)
			      {
				global_heapsize = heapsize;
				glh = LocalHeap(heapsize, "python-comp lh", true);
                                bool first_time = true;
                                if (first_time)
                                  { first_time = false; cerr << "warning: use SetHeapSize(size) instead of heapsize=size" << endl; }                                                
			      }
                           BitArray mask(ma->GetNRegions(vb));
                           mask.Set();
                            {
                              py::gil_scoped_acquire acquire;
                              py::extract<Region> defon_region(definedon);
                              if (defon_region.check())
                                vb = VorB(defon_region());
                              if(defon_region.check())
                                for(auto i : Range(ma->GetNRegions(vb)))
                                  if(!defon_region().Mask().Test(i))
                                    mask.Clear(i);
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
#ifdef PARALLEL
				  sum(0) = MyMPI_AllReduce(sum(0));
#endif
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
#ifdef PARALLEL
				  sum(0) = MyMPI_AllReduce(sum(0));
#endif
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
	py::arg("heapsize") = 1000000,
        py::call_guard<py::gil_scoped_release>())
    ;

  m.def("SymbolicLFI",
          [](spCF cf, VorB vb, bool element_boundary,
             bool skeleton, py::object definedon,
             IntegrationRule ir, int bonus_intorder, py::object definedonelem,
             bool simd_evaluate) 
           {
             py::extract<Region> defon_region(definedon);
             if (defon_region.check())
               vb = VorB(defon_region());

             shared_ptr<LinearFormIntegrator> lfi;
             if (!skeleton)
               lfi = make_shared<SymbolicLinearFormIntegrator> (cf, vb, element_boundary);
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
           py::arg("simd_evaluate")=true
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
           int bonus_intorder, py::object definedonelem, bool simd_evaluate)
        -> shared_ptr<BilinearFormIntegrator>
           {
             py::extract<Region> defon_region(definedon);
             if (defon_region.check())
               vb = VorB(defon_region());

             auto bfi = make_shared<SymbolicEnergy> (cf, vb, element_boundary);
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
        py::arg("simd_evaluate")=true
          );


  
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
             FlatVector<> x0(spaces[0]->GetSpacialDimension(),&x0_help[0]);
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
                                shared_ptr<GF> gfustd, int heapsize )
            {
              if (heapsize > global_heapsize)
              {
                global_heapsize = heapsize;
                glh = LocalHeap(heapsize, "python-comp lh", true);
              }
              static Timer tall("comp.Transfer2StdMesh"); RegionTimer rall(tall);
              Transfer2StdMesh(gfutp.get(),gfustd.get(),glh);
              return;
             },
             py::arg("gftp"),
             py::arg("gfstd"),
             py::arg() = 1000000,
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
         py::arg("only_element") = -1,
        py::call_guard<py::gil_scoped_release>()
         )
    .def("Do", [](shared_ptr<BaseVTKOutput> self, int heapsize)
                               { 
                                 LocalHeap lh (heapsize, "VTKOutput-heap");
                                 self->Do(lh);
                               },
         py::arg("heapsize")=1000000,
        py::call_guard<py::gil_scoped_release>())
    .def("Do", [](shared_ptr<BaseVTKOutput> self, const BitArray * drawelems, int heapsize)
                               { 
                                 LocalHeap lh (heapsize, "VTKOutput-heap");
                                 self->Do(lh, drawelems);
                               },
         py::arg("drawelems"),py::arg("heapsize")=1000000,
        py::call_guard<py::gil_scoped_release>())
    
    ;

  /////////////////////////////////////////////////////////////////////////////////////
}




PYBIND11_MODULE(libngcomp, m) {
  m.attr("__name__") = "comp";
  ExportNgcomp(m);
}

#endif // NGS_PYTHON
