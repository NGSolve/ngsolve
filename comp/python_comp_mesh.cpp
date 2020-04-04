#ifdef NGS_PYTHON

#include <regex>

#include <core/python_ngcore.hpp>

#include "../ngstd/python_ngstd.hpp"
#include <comp.hpp>
using namespace ngcomp;

inline auto Nr2Vert(size_t nr) {  return NodeId(NT_VERTEX,nr); };
inline auto Nr2Edge(size_t nr) {  return NodeId(NT_EDGE,nr); };
inline auto Nr2Face(size_t nr) {  return NodeId(NT_FACE,nr); };
inline auto Nr2VolElement(size_t nr) {  return ElementId(VOL,nr); };

void ExportPml(py::module &m);

void ExportNgcompMesh (py::module &m)
{
  py::module pml = m.def_submodule("pml", "module for perfectly matched layers");
  ExportPml(pml);


  py::enum_<VorB>(m, "VorB", "Enum specifying the codimension. VOL is volume, BND is boundary and BBND is codimension 2 (edges in 3D, points in 2D)")
    .value("VOL", VOL)
    .value("BND", BND)
    .value("BBND", BBND)
    .value("BBBND", BBBND)
    .export_values()
    ;
  
  py::class_<ElementId> (m, "ElementId", 
                         docu_string(R"raw_string(
An element identifier containing element number and Volume/Boundary flag

3 __init__ overloads:

1)

Parameters:

vb : ngsolve.comp.VorB
  input Volume or Boundary (VOL, BND, BBND, BBBND)

nr : int
  input element number


2)

Parameters:

nr : int
  input element number


3)

Parameters:

el : ngcomp::Ngs_Element
  input Ngs element

)raw_string"))
    .def(py::init<VorB,size_t>(), py::arg("vb"), py::arg("nr"))
    .def(py::init<size_t>(), py::arg("nr"))
    .def(py::init<Ngs_Element>(), py::arg("el"))
    .def("__str__", &ToString<ElementId>)
    .def_property_readonly("nr", &ElementId::Nr, "the element number")    
    .def("VB", &ElementId::VB, "VorB of element")
    .def_property_readonly("valid", [] (ElementId ei) { return ei.Nr() != -1; }, "is element valid")
    .def(py::self!=py::self)
    .def(py::self==py::self)
    .def("__hash__" , &ElementId::Nr)
    ;
  
  m.def("BndElementId",[] (int nr) { return ElementId(BND,nr); },
          py::arg("nr"), docu_string(R"raw_string(
Creates an element-id for a boundary element

Parameters:

nr : int
  input Bnd element number

)raw_string"))
    ;


  
  py::class_<NodeId> (m, "NodeId",
                      "an node identifier containing node type and node nr")
    .def(py::init<NODE_TYPE,size_t>(), py::arg("type"), py::arg("nr"))
    .def("__str__", &ToString<NodeId>)
    .def("__repr__", &ToString<NodeId>)
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
                             auto& mesh = node.Mesh();
                             switch (StdNodeType(node.GetType(), mesh.GetDimension()))
                               {
                               case NT_EDGE:
                                 return MakePyTuple(Substitute(ArrayObject(mesh.GetEdgePNums(node.GetNr())), Nr2Vert));
                               case NT_FACE:
                                 return MakePyTuple(Substitute(mesh.GetFacePNums(node.GetNr()), Nr2Vert));
                               case NT_CELL:
                                 return MakePyTuple(Substitute(mesh.GetElPNums(ElementId(VOL, node.GetNr())),
                                                               Nr2Vert));
                               default:
                                 throw py::type_error("vertices only available for edge, face and cell nodes\n");
                               }
                           }, "tuple of global vertex numbers")
    .def_property_readonly("edges",[](MeshNode & node) -> py::tuple
                           {
                             auto& mesh = node.Mesh();
                             switch(StdNodeType(node.GetType(), mesh.GetDimension()))
                               {
                               case NT_VERTEX:
                                 {
                                   Array<int> enums;
                                   for (auto el : mesh.GetVertexElements(node.GetNr()))
                                     for (auto edge : mesh.GetElement(ElementId(VOL,el)).Edges())
                                       if (!enums.Contains(edge))
                                         enums.Append(edge);
                                   QuickSort (enums);
                                   return MakePyTuple(Substitute(enums, Nr2Edge));
                                 }
                               case NT_FACE:
                                 return MakePyTuple(Substitute(mesh.GetFaceEdges(node.GetNr()), Nr2Edge));
                               case NT_CELL:
                                 return MakePyTuple(Substitute(mesh.GetElEdges(ElementId(VOL,node.GetNr())),
                                                               Nr2Edge));
                               default:
                                 throw py::type_error("edges only available for face and cell nodes\n");
                               }
                           }, "tuple of global edge numbers")
    .def_property_readonly("faces", [](MeshNode & node) -> py::tuple
                           {
                             auto & mesh = node.Mesh();
                             if (StdNodeType(node.GetType(), mesh.GetDimension()) == NT_VERTEX)  
                               {
                                 Array<int> fnums;
                                 for (auto el : mesh.GetVertexElements(node.GetNr()))
                                   for (auto face : mesh.GetElement(ElementId(VOL,el)).Faces())
                                     if (!fnums.Contains(face))
                                       fnums.Append(face);
                                 QuickSort (fnums);
                                 return MakePyTuple(Substitute(fnums, Nr2Face));
                               }
                             if (node.GetType() == NT_EDGE)
                               {
                                 Array<int> faces;
                                 mesh.GetEdgeFaces(node.GetNr(), faces);
                                 return MakePyTuple(Substitute(faces, Nr2Face));
                               }
                             if (node.GetType() == NT_CELL)
                               return MakePyTuple(Substitute(mesh.GetElFacets(ElementId(VOL, node.GetNr())),
                                                             Nr2Face));
                             throw py::type_error("faces only available for cell nodes\n");
                           }, "tuple of global face numbers")
    
    .def_property_readonly("point", [](MeshNode & node) -> py::tuple
                           {
                             auto & mesh = node.Mesh();
                             if (node.GetType() == NT_VERTEX)
                               switch (mesh.GetDimension())
                                 {
                                 case 1:
                                   {
                                     auto p = mesh.GetPoint<1>(node.GetNr());
                                     return py::make_tuple(p(0));
                                   }
                                 case 2:
                                   {
                                     auto p = mesh.GetPoint<2>(node.GetNr());
                                     return py::make_tuple(p(0), p(1));
                                   }
                                 case 3:
                                   {
                                     auto p = mesh.GetPoint<3>(node.GetNr());
                                     return py::make_tuple(p(0), p(1), p(2));
                                   }
                                 }
                             throw py::type_error("point only available for vertex nodes\n");
                           }, "vertex coordinates")


    .def_property_readonly("elements",[](MeshNode & node) -> py::tuple
                           {
                             auto& mesh = node.Mesh();
                             switch(node.GetType())
                               {
                               case NT_VERTEX:
                                 return MakePyTuple(Substitute(mesh.GetVertexElements(node.GetNr()), Nr2VolElement));
                               case NT_EDGE:
                                 {
                                   Array<int> elnums;
                                   mesh.GetEdgeElements(node.GetNr(), elnums);
                                   return MakePyTuple(Substitute(elnums, Nr2VolElement));
                                 }
                               case NT_FACE:
                                 {
                                   Array<int> elnums;
                                   mesh.GetFaceElements(node.GetNr(), elnums);
                                   return MakePyTuple(Substitute(elnums, Nr2VolElement));
                                 }
                               default:
                                 throw py::type_error("elements only available for vertex nodes\n");
                               }
                           }, "tuple of global element-ids")
    
    ;

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




  py::class_<Ngs_Element>(m, "Ngs_Element")
    .def_property_readonly("nr", &Ngs_Element::Nr, "the element number")    
    .def("VB", &Ngs_Element::VB, "VorB of element")
    .def_property_readonly("valid", [] (ElementId ei) { return ei.Nr() != -1; }, "is element valid")    
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
    .def_property_readonly("facets", [](Ngs_Element &el)
                           {
                             switch (ElementTopology::GetSpaceDim(el.GetType()))
                               {
                               case 1:
                                 return MakePyTuple(Substitute(el.Vertices(), Nr2Vert));
                               case 2:
                                 return MakePyTuple(Substitute(el.Edges(), Nr2Edge));
                               case 3:
                                 return MakePyTuple(Substitute(el.Faces(), Nr2Face));
                               default:
                                 throw Exception ("Illegal dimension in Ngs_Element.faces");
                               }
                           },
                           "tuple of global face, edge or vertex numbers")
    .def_property_readonly("elementnode", [](Ngs_Element &el)
                           {
                             switch (ElementTopology::GetSpaceDim(el.GetType()))
                               {
                               case 1:
                                 return Nr2Edge (el.Edges()[0]);
                               case 2:
                                 return Nr2Face (el.Faces()[0]);
                               case 3:
                                 return NodeId(NT_CELL, el.Nr());
                               default:
                                 throw Exception ("Illegal dimension in Ngs_Element.elementnode");
                               }
                           },
                           "inner node, i.e. cell, face or edge node for 3D/2D/1D")
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
  
  //////////////////////////////////////////////////////////////////////////////////////////

  py::class_<Region> (m, "Region", "a subset of volume or boundary elements")
    .def(py::init<shared_ptr<MeshAccess>,VorB,string>(), py::arg("mesh"), py::arg("vb"), py::arg("name"))
    .def(py::init<shared_ptr<MeshAccess>,VorB,BitArray>(), py::arg("mesh"), py::arg("vb"), py::arg("mask"))
    .def("Mask",[](Region & reg)->BitArray { return reg.Mask(); }, "BitArray mask of the region")
    .def("VB", [](Region & reg) { return VorB(reg); }, "VorB of the region")
    .def(py::self + py::self)
    .def(py::self + string())
    .def(py::self - py::self)
    .def(py::self - string())
    .def(~py::self)
    ;

  py::implicitly_convertible <Region, BitArray> ();


  //////////////////////////////////////////////////////////////////////////////////////////
  
  
  typedef PML_Transformation PML;
  
  py::class_<MeshAccess, shared_ptr<MeshAccess>> mesh_access(m, "Mesh", docu_string(R"raw_string(
NGSolve interface to the Netgen mesh. Provides access and functionality
to use the mesh for finite element calculations.

Parameters:

mesh (netgen.Mesh): a mesh generated from Netgen


)raw_string") , py::dynamic_attr());
    
  mesh_access
    .def(py::init<shared_ptr<netgen::Mesh>>(),
         py::arg("ngmesh"),
         "Make an NGSolve-mesh from a Netgen-mesh")

    .def(py::init([](const string & filename, NgMPI_Comm comm)
                  {
		    // MPI_Comm comm = c ? c->comm : ngs_comm;
                    NGSOStream::SetGlobalActive (comm.Rank()==0);
                    return make_shared<MeshAccess>(filename, comm);
                  }),
         py::arg("filename"), py::arg("comm")=NgMPI_Comm{},
         "Load a mesh from file.\n"
         "In MPI-parallel mode the mesh is distributed over the MPI-group given by the communicator (WIP!)")
    
    .def("__eq__",
         [] (shared_ptr<MeshAccess> self, shared_ptr<MeshAccess> other)
         { return self == other; }, py::arg("mesh"))
     .def_property_readonly("comm", [](const MeshAccess& ma) -> NgMPI_Comm
                            { return ma.GetCommunicator(); },
                            "MPI-communicator the Mesh lives in")
   
    .def(NGSPickle<MeshAccess>())
    /*
    .def("LoadMesh", static_cast<void(MeshAccess::*)(const string &)>(&MeshAccess::LoadMesh),
         "Load mesh from file")
    */
    
    .def("Elements", static_cast<ElementRange(MeshAccess::*)(VorB)const> (&MeshAccess::Elements),
	 (py::arg("VOL_or_BND")=VOL),
         docu_string("Return an iterator over elements on VOL/BND"))

    .def("__getitem__", static_cast<Ngs_Element(MeshAccess::*)(ElementId)const> (&MeshAccess::operator[]),
         "Return Ngs_Element from given ElementId")
    
    .def("__getitem__", [](MeshAccess & self, NodeId ni)
         {
           if (ni.GetNr() >= self.GetNNodes(ni.GetType()))
             throw py::index_error("illegal node number");
           return MeshNode(ni, self);
         },
         "Return MeshNode from given NodeId")

    .def ("GetNE", static_cast<size_t(MeshAccess::*)(VorB)const> (&MeshAccess::GetNE),
          docu_string("Return number of elements of codimension VorB."))
    
    .def_property_readonly ("nv", &MeshAccess::GetNV, "number of vertices")
    .def_property_readonly ("ne",  static_cast<size_t(MeshAccess::*)()const> (&MeshAccess::GetNE), "number of volume elements")
    .def_property_readonly ("nedge", &MeshAccess::GetNEdges, "number of edges")
    .def_property_readonly ("nface", &MeshAccess::GetNFaces, "number of faces")
    .def_property_readonly ("nfacet", &MeshAccess::GetNFacets, "number of facets")
    .def ("nnodes", &MeshAccess::GetNNodes, "number of nodes given type")
    .def_property_readonly ("dim", &MeshAccess::GetDimension, "mesh dimension")
    .def_property_readonly ("ngmesh", &MeshAccess::GetNetgenMesh, "the Netgen mesh")

    
    .def_property_readonly ("vertices", [] (shared_ptr<MeshAccess> mesh)
          {
            return T_Range<MeshNode> (MeshNode(NodeId(NT_VERTEX, 0), *mesh),
                                      MeshNode(NodeId(NT_VERTEX, mesh->GetNNodes(NT_VERTEX)), *mesh));
          }, "iterable of mesh vertices")
    
    .def_property_readonly ("edges", [] (shared_ptr<MeshAccess> mesh)
          {
            return T_Range<MeshNode> (MeshNode(NodeId(NT_EDGE, 0), *mesh),
                                      MeshNode(NodeId(NT_EDGE, mesh->GetNNodes(NT_EDGE)), *mesh));
          }, "iterable of mesh edges")
    
    .def_property_readonly ("faces", [] (shared_ptr<MeshAccess> mesh)
          {
            return T_Range<MeshNode> (MeshNode(NodeId(NT_FACE, 0), *mesh),
                                      MeshNode(NodeId(NT_FACE, mesh->GetNNodes(NT_FACE)), *mesh));
          }, "iterable of mesh faces")

    .def_property_readonly ("facets", [] (shared_ptr<MeshAccess> mesh)
          {
            auto nt = StdNodeType(NT_FACET, mesh->GetDimension());
            return T_Range<MeshNode> (MeshNode(NodeId(nt, 0), *mesh),
                                      MeshNode(NodeId(nt, mesh->GetNNodes(nt)), *mesh));
          }, "iterable of mesh facets")

    .def("nodes", [] (shared_ptr<MeshAccess> mesh, NODE_TYPE type)
         {
           return T_Range<MeshNode> (MeshNode(NodeId(type, 0), *mesh),
                                     MeshNode(NodeId(type, mesh->GetNNodes(type)),*mesh));
         }, py::arg("node_type"), "iterable of mesh nodes of type node_type")

    /*
    .def ("GetTrafo", 
          static_cast<ElementTransformation&(MeshAccess::*)(ElementId,Allocator&)const>
          (&MeshAccess::GetTrafo), 
          py::return_value_policy::reference)
    */

    .def("GetPeriodicNodePairs", [](MeshAccess& self, NODE_TYPE type)
         {
           py::list pairs;
           for(auto idnr : Range(self.GetNPeriodicIdentifications()))
             {
               const auto& periodic_nodes = self.GetPeriodicNodes(type, idnr);
               for(auto pair : periodic_nodes)
                 pairs.append(py::make_tuple(py::make_tuple(pair[0], pair[1]),idnr));
             }
           return pairs;
         }, "returns list of periodic nodes with their identification number as [((master_nr, slave_nr),idnr),...]")
    
    .def ("GetTrafo",
          [](MeshAccess & ma, ElementId id)
          { return &ma.GetTrafo(id, global_alloc); }, py::arg("eid"),
          py::return_value_policy::take_ownership, "returns element transformation of given element id")

    .def("SetDeformation", 
	 [](MeshAccess & ma, shared_ptr<GridFunction> gf)
         { ma.SetDeformation(gf); }, py::arg("gf"),
         docu_string("Deform the mesh with the given GridFunction"))

    .def("UnsetDeformation", [](MeshAccess & ma){ ma.SetDeformation(nullptr);}, "Unset the deformation")

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
         "Set PML transformation on domain"
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
          }, py::arg("definedon"), "Unset PML transformation on domain")
    
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
        "Return list of pml transformations"
    )
    .def("GetPMLTrafo", [](MeshAccess & ma, int domnr) {
        if (ma.GetPMLTrafos()[domnr])
     	  return ma.GetPMLTrafos()[domnr-1];
        else
          throw Exception("No PML Trafo set"); 
        },
        py::arg("dom")=1,
        "Return pml transformation on domain dom"
        )

    .def("GetMaterials",
	 [](const MeshAccess & ma)
         {
           return MakePyTuple(ma.GetMaterials(VOL));
         },
         "Return list of material names")

    .def("Materials",
	 [](const shared_ptr<MeshAccess> & ma, string pattern) 
	  {
            return Region (ma, VOL, pattern);
	  },
         py::arg("pattern"),
	 docu_string("Return mesh-region matching the given regex pattern"))
    
    .def("Materials",
	 [](const shared_ptr<MeshAccess> & ma, vector<int> domains)
	  {
            cout << "warning: Materials( [int list] ) is deprecated, pls generate Region" << endl;
            BitArray mask(ma->GetNDomains());
            mask.Clear();
            for (auto i : domains)
              if (i >= 0 && i < mask.Size())
                mask.SetBit(i);
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
	 "Return list of boundary condition names")

    .def("Boundaries",
	 [](const shared_ptr<MeshAccess> & ma, string pattern)
	  {
            return Region (ma, BND, pattern);
	  },
         py::arg("pattern"),
	 "Return boundary mesh-region matching the given regex pattern")

    .def("Boundaries",
	 [](const shared_ptr<MeshAccess> & ma, vector<int> bnds)
	  {
            cout << "warning: Boundaries( [int list] ) is deprecated, pls generate Region" << endl;            
            BitArray mask(ma->GetNBoundaries());
            mask.Clear();
            for (auto i : bnds)
              if (i >= 0 && i < mask.Size())
                mask.SetBit(i);
              else
                throw Exception ("boundary index "+ToString(i)+" out of range [0,"+ToString(mask.Size())+")");
            return Region (ma, BND, mask);
	  },
         py::arg("bnds"),
	 "Generate boundary mesh-region by boundary condition numbers")

    .def("GetBBoundaries",
	 [](const MeshAccess & ma)
	  {
            return MakePyTuple(ma.GetMaterials(BBND));
	  },
	 "Return list of boundary conditions for co dimension 2")
    
    .def("BBoundaries", [](const shared_ptr<MeshAccess> & ma, string pattern)
	  {
	    return Region (ma, BBND, pattern);
	  },
	 (py::arg("self"), py::arg("pattern")),
	 "Return co dim 2 boundary mesh-region matching the given regex pattern")

    .def("GetBBBoundaries",
	 [](const MeshAccess & ma)
	  {
            return MakePyTuple(ma.GetMaterials(BBBND));
	  },
	 "Return list of boundary conditions for co dimension 3")

    .def("BBBoundaries", [](const shared_ptr<MeshAccess> & ma, string pattern)
	  {
	    return Region (ma, BBBND, pattern);
	  },
	 (py::arg("self"), py::arg("pattern")),
	 "Return co dim 3 boundary mesh-region matching the given regex pattern")

    // TODO: explain how to mark elements
    .def("Refine",
         [](MeshAccess & ma, bool mark_surface_elements)
          {
            if (!mark_surface_elements)
              {
                for (ElementId ei : ma.Elements(BND))
                  ma.SetRefinementFlag(ei, false);
              }
            ma.Refine();
          },py::call_guard<py::gil_scoped_release>(),
         py::arg("mark_surface_elements")=false,
	 "Local mesh refinement based on marked elements, uses element-bisection algorithm")

    .def("RefineHP",
         [](MeshAccess & ma, int levels, double factor)
          {
            ma.HPRefinement(levels, factor);
            // Ng_HPRefinement(levels, factor);
            // ma.UpdateBuffers();
          },
         py::arg("levels"), py::arg("factor")=0.125,
	 "Geometric mesh refinement towards marked vertices and edges, uses factor for placement of new points")

    .def("_updateBuffers", &MeshAccess::UpdateBuffers, "Update NGSolve mesh information, needs to be called if Netgen mesh changes")
    .def("SetRefinementFlag", &MeshAccess::SetRefinementFlag,
         py::arg("ei"), py::arg("refine"),
	 "Set refinementflag for mesh-refinement")

    .def("SetRefinementFlags", [&](MeshAccess & ma, std::vector<bool> flags)
         {
           for (ElementId ei : ma.Elements(VOL))
             ma.SetRefinementFlag (ei, flags[ei.Nr()]);
         }, py::arg("refine"),
	 "Set refinementflags for mesh-refinement")
    
    .def("GetParentElement", static_cast<ElementId(MeshAccess::*)(ElementId)const> (&MeshAccess::GetParentElement),
         py::arg("ei"),
         "Return parent element id on refined mesh")

    .def("GetParentVertices", [](MeshAccess & ma, int vnum)
          {
            auto parents = ma.GetParentNodes (vnum);
            return py::make_tuple(parents[0], parents[1]);
          },
         py::arg("vnum"),
         "Return parent vertex numbers on refined mesh")
    .def("GetHPElementLevel", &MeshAccess::GetHPElementLevel,
         py::arg("ei"),
         "THIS FUNCTION IS WIP!\n Return HP-refinement level of element")
    .def("SetElementOrder",
         [](MeshAccess & ma, ElementId id, int order)
         {
           ma.SetElOrder(id.Nr(), order);
         }, py::arg("eid"), py::arg("order"), "For backward compatibility, not recommended to use")
    
    .def("Curve", &MeshAccess::Curve,
         py::arg("order"),
         "Curve the mesh elements for geometry approximation of given order")

    .def("GetCurveOrder", &MeshAccess::GetCurveOrder,
	 "")

    .def("Contains",
         [](MeshAccess & ma, double x, double y, double z) 
          {
            IntegrationPoint ip;
            int elnr = ma.FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            return (elnr >= 0);
          }, 
         py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0
	 ,"Check if the point (x,y,z) is in the meshed domain (is inside a volume element)")

    ;
    PyDefVectorized(mesh_access, "__call__",
         [](MeshAccess* ma, double x, double y, double z, VorB vb)
          {
            IntegrationPoint ip;
            int elnr;
            if (vb == VOL)
              elnr = ma->FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            else
              elnr = ma->FindSurfaceElementOfPoint(Vec<3>(x, y, z), ip, true);
            return MeshPoint { ip(0), ip(1), ip(2), ma, vb, elnr };
          },
         py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0,
         py::arg("VOL_or_BND") = VOL,
	 docu_string("Get a MappedIntegrationPoint in the point (x,y,z) on the matching volume (VorB=VOL, default) or surface (VorB=BND) element. BBND elements aren't supported"));

  
    m.def("BoundaryFromVolumeCF", 
          [] (shared_ptr<CoefficientFunction> vol_cf)
          { return MakeBoundaryFromVolumeCoefficientFunction(vol_cf); }, py::arg("vol_cf"),
          docu_string(R"raw_string(Allows the evaluation of volumetric functions on the boundary.

When evaluated on a boundary element, this function searches for the associated
volume element, transforms the local coordinates, and evaluates the function in the
volume. A typical use case is to visualize L2-functions, or mechanical stresses at
the boundary.

It is different from the boundary Trace()-operator. The trace provides a function
which is defined by boundary degrees of freedom only. E.g. the trace of an H(div)
function is only the normal component, while the BoundaryFromVolumeCF gives the
whole function. Obviously, the Trace() function is cheaper to evaluate.

If called on an interface, it evaluates from one side (which one is not specified).
If the function is only defined on one side, this side will be taken. One can use
a domain-wise CF to define a function only locally:
uloc = CoefficientFunction( [None, None, u, None] )
)raw_string")
          );

}







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
         }, py::arg("pml"))
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

origin is a list/tuple with as many entries as dimenson)raw_string");

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


#endif
