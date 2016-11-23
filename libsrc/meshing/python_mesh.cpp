#ifdef NG_PYTHON

#include <../general/ngpython.hpp>

#include <mystdlib.h>
#include "meshing.hpp"
#include <csg.hpp>
#include <geometry2d.hpp>
#include <../interface/writeuser.hpp>


using namespace netgen;

namespace netgen
{
  extern shared_ptr<NetgenGeometry> ng_geometry;
}


template <typename T, int BASE = 0, typename TIND = int>
void ExportArray (py::module &m)
{
  using TA = Array<T,BASE,TIND>;
  string name = string("Array_") + typeid(T).name();
  py::class_<Array<T,BASE,TIND>>(m, name.c_str())
    .def ("__len__", [] ( Array<T,BASE,TIND> &self ) { return self.Size(); } )
    .def ("__getitem__", 
          FunctionPointer ([](Array<T,BASE,TIND> & self, TIND i) -> T&
                           {
                             if (i < BASE || i >= BASE+self.Size())
                               throw py::index_error();
                             return self[i];
                           }),
          py::return_value_policy::reference)
    .def("__iter__", [] ( TA & self) {
	return py::make_iterator (self.begin(),self.end());
      }, py::keep_alive<0,1>()) // keep array alive while iterator is used

    ;
}

void TranslateException (const NgException & ex)
{
  string err = string("Netgen exception: ")+ex.What();
  PyErr_SetString(PyExc_RuntimeError, err.c_str());
}


DLL_HEADER void ExportNetgenMeshing(py::module &m) 
{
  
  
  py::class_<PointIndex>(m, "PointId")
    .def(py::init<int>())
    .def("__repr__", &ToString<PointIndex>)
    .def("__str__", &ToString<PointIndex>)
    .def_property_readonly("nr", &PointIndex::operator int)
    .def("__eq__" , FunctionPointer( [](PointIndex &self, PointIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](PointIndex &self ) { return static_cast<int>(self); }) )
    ;

  py::class_<ElementIndex>(m, "ElementId3D")
    .def(py::init<int>())
    .def("__repr__", &ToString<ElementIndex>)
    .def("__str__", &ToString<ElementIndex>)
    .def_property_readonly("nr", &ElementIndex::operator int)
    .def("__eq__" , FunctionPointer( [](ElementIndex &self, ElementIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](ElementIndex &self ) { return static_cast<int>(self); }) )
    ;


  py::class_<SurfaceElementIndex>(m, "ElementId2D")
    .def(py::init<int>())
    .def("__repr__", &ToString<SurfaceElementIndex>)
    .def("__str__", &ToString<SurfaceElementIndex>)
    .def_property_readonly("nr", &SurfaceElementIndex::operator int)
    .def("__eq__" , FunctionPointer( [](SurfaceElementIndex &self, SurfaceElementIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](SurfaceElementIndex &self ) { return static_cast<int>(self); }) )
    ;

  py::class_<SegmentIndex>(m, "ElementId1D")
    .def(py::init<int>())
    .def("__repr__", &ToString<SegmentIndex>)
    .def("__str__", &ToString<SegmentIndex>)
    .def_property_readonly("nr", &SegmentIndex::operator int)
    .def("__eq__" , FunctionPointer( [](SegmentIndex &self, SegmentIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](SegmentIndex &self ) { return static_cast<int>(self); }) )
    ;



  /*  
  py::class_<Point<3>> ("Point")
    .def(py::init<double,double,double>())
    ;
  */

  py::class_<MeshPoint /* ,py::bases<Point<3>> */ >(m, "MeshPoint")
    .def(py::init<Point<3>>())
    .def("__str__", &ToString<MeshPoint>)
    .def("__repr__", &ToString<MeshPoint>)
    .def_property_readonly("p", FunctionPointer([](const MeshPoint & self)
                                       {
                                         py::list l;
                                         l.append ( py::cast(self[0]) );
                                         l.append ( py::cast(self[1]) );
                                         l.append ( py::cast(self[2]) );
                                         return py::tuple(l);
                                       }))
    .def("__getitem__", FunctionPointer([](const MeshPoint & self, int index) {
	  if(index<0 || index>2)
              throw py::index_error();
	  return self[index];
	}))
    ;
  
  py::class_<Element>(m, "Element3D")
    .def("__init__", [](Element *instance, int index, py::list vertices)
                           {
                             if (py::len(vertices) == 4)
                               {
                                 new (instance) Element(TET);
                                 for (int i = 0; i < 4; i++)
                                   (*instance)[i] = py::extract<PointIndex>(vertices[i])();
                                 instance->SetIndex(index);
                               }
                             else if (py::len(vertices) == 6)
                               {
                                 new (instance) Element(PRISM);
                                 for (int i = 0; i < 6; i++)
                                   (*instance)[i] = py::extract<PointIndex>(vertices[i])();
                                 instance->SetIndex(index);
                               }
                             else if (py::len(vertices) == 8)
                               {
                                 new (instance) Element(HEX);
                                 for (int i = 0; i < 8; i++)
                                   (*instance)[i] = py::extract<PointIndex>(vertices[i])();
                                 instance->SetIndex(index);
                               }
                             else
                               throw NgException ("cannot create element");                             
                           },
          py::arg("index")=1,py::arg("vertices"),
         "create volume element"
         )
    .def("__repr__", &ToString<Element>)
    .def_property("index", &Element::GetIndex, &Element::SetIndex)
    .def_property_readonly("vertices", 
                  FunctionPointer ([](const Element & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < self.GetNV(); i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    ;

  py::class_<Element2d>(m, "Element2D")
    .def("__init__", 
         [](Element2d *instance, int index, py::list vertices)
                           {
                             if (py::len(vertices) == 3)
                               {
                                 new (instance) Element2d(TRIG);
                                 for (int i = 0; i < 3; i++)
                                   (*instance)[i] = py::extract<PointIndex>(vertices[i])();
                                 instance->SetIndex(index);
                               }
                             else
                               {
                                 new (instance) Element2d(QUAD);
                                 for (int i = 0; i < 4; i++)
                                   (*instance)[i] = py::extract<PointIndex>(vertices[i])();
                                 instance->SetIndex(index);
                               }
                               
                           },
          py::arg("index")=1,py::arg("vertices"),
         "create surface element"
         )
    .def_property("index", &Element2d::GetIndex, &Element2d::SetIndex)
    .def_property_readonly("vertices",
                  FunctionPointer([](const Element2d & self) -> py::list
                                  {
                                    py::list li;
                                    for (int i = 0; i < self.GetNV(); i++)
                                      li.append(py::cast(self[i]));
                                    return li;
                                  }))
    ;

  py::class_<Segment>(m, "Element1D")
    .def("__init__",
         [](Segment *instance, py::list vertices, py::list surfaces, int index)
                           {
                             new (instance) Segment();
                             for (int i = 0; i < 2; i++)
                               (*instance)[i] = py::extract<PointIndex>(vertices[i])();
                             instance -> si = index;
			     // needed for codim2 in 3d
			     instance -> edgenr = index;
                             if (len(surfaces))
                               {
                                 instance->surfnr1 = py::extract<int>(surfaces[0])();
                                 instance->surfnr2 = py::extract<int>(surfaces[1])();
                               }
                           },
          py::arg("vertices"),
           py::arg("surfaces")=py::list(),
           py::arg("index")=1,
         "create segment element"
         )
    .def("__repr__", &ToString<Segment>)
    .def_property_readonly("vertices", 
                  FunctionPointer ([](const Segment & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < 2; i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    .def_property_readonly("surfaces", 
                  FunctionPointer ([](const Segment & self) -> py::list
                                   {
                                     py::list li;
                                     li.append (py::cast(self.surfnr1));
                                     li.append (py::cast(self.surfnr2));
                                     return li;
                                   }))
    .def_property_readonly("index", FunctionPointer([](const Segment &self) -> size_t
		  {
		    return self.si;
		  }))
    .def_property_readonly("edgenr", FunctionPointer([](const Segment & self) -> size_t
						     {
						       return self.edgenr;
						     }))
    ;


  py::class_<Element0d>(m, "Element0D")
    .def("__init__",
         [](Element0d *instance, PointIndex vertex, int index)
                           {
                             new (instance) Element0d;
                             instance->pnum = vertex;
                             instance->index = index;
                           },
         py::arg("vertex"),
         py::arg("index")=1,
         "create point element"
         )
    .def("__repr__", &ToString<Element0d>)
    .def_property_readonly("vertices", 
                  FunctionPointer ([](const Element0d & self) -> py::list
                                   {
                                     py::list li;
                                     li.append (py::cast(self.pnum));
                                     return li;
                                   }))
    ;
  
  
  


  py::class_<FaceDescriptor>(m, "FaceDescriptor")
    .def(py::init<const FaceDescriptor&>())
    .def("__init__", 
         [](FaceDescriptor *instance, int surfnr, int domin, int domout, int bc)
                           {
                             new (instance) FaceDescriptor();
                             instance->SetSurfNr(surfnr);
                             instance->SetDomainIn(domin);
                             instance->SetDomainOut(domout);
                             instance->SetBCProperty(bc);
                           },
           py::arg("surfnr")=1, 
           py::arg("domin")=1,
           py::arg("domout")=py::int_(0),
           py::arg("bc")=py::int_(0),
         "create facedescriptor")
    .def("__str__", &ToString<FaceDescriptor>)
    .def("__repr__", &ToString<FaceDescriptor>)
    .def_property("surfnr", &FaceDescriptor::SurfNr, &FaceDescriptor::SetSurfNr)
    .def_property("domin", &FaceDescriptor::DomainIn, &FaceDescriptor::SetDomainIn)
    .def_property("domout", &FaceDescriptor::DomainOut, &FaceDescriptor::SetDomainOut)
    .def_property("bc", &FaceDescriptor::BCProperty, &FaceDescriptor::SetBCProperty)
    .def_property_readonly("bcname", FunctionPointer ([](FaceDescriptor & self) -> string { return self.GetBCName(); }))
    .def("SetSurfaceColor", [](FaceDescriptor & self, py::list color )
          {
            Vec3d c;
            c.X() = py::extract<double>(color[0])();
            c.Y() = py::extract<double>(color[1])();
            c.Z() = py::extract<double>(color[2])();
            self.SetSurfColour(c);
          })
    ;

  

  ExportArray<Element>(m);
  ExportArray<Element2d>(m);
  ExportArray<Segment>(m);
  ExportArray<Element0d>(m);
  ExportArray<MeshPoint,PointIndex::BASE,PointIndex>(m);
  ExportArray<FaceDescriptor>(m);

  py::implicitly_convertible< int, PointIndex>();
  
  py::class_<Mesh,shared_ptr<Mesh>>(m, "Mesh")
    // .def(py::init<>("create empty mesh"))

    .def("__init__",
         [](Mesh *instance, int dim)
                           {
                             new (instance) Mesh();
                             instance->SetDimension(dim);
                           },
           py::arg("dim")=3
          )

    
    .def("__str__", &ToString<Mesh>)
    .def("Load",  FunctionPointer 
	 ([](Mesh & self, const string & filename)
	  {
            istream * infile;
            if (filename.find(".vol.gz") != string::npos)
              infile = new igzstream (filename.c_str());
            else
              infile = new ifstream (filename.c_str());
	    // ifstream input(filename);
#ifdef PARALLEL
	    // int id;
	    MPI_Comm_rank(MPI_COMM_WORLD, &id);
	    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

	    if (id == 0)
	      {
		self.Load(*infile);	      
		self.Distribute();
	      }
	    else
	      {
		self.SendRecvMesh();
	      }
#else
	    self.Load(*infile);
#endif
	    for (int i = 0; i < geometryregister.Size(); i++)
	      {
		NetgenGeometry * hgeom = geometryregister[i]->LoadFromMeshFile (*infile);
		if (hgeom)
		  {
		    ng_geometry.reset (hgeom);
		    break;
		  }
	      }
	  }))
    // static_cast<void(Mesh::*)(const string & name)>(&Mesh::Load))
    .def("Save", static_cast<void(Mesh::*)(const string & name)const>(&Mesh::Save))
    .def("Export",
         [] (Mesh & self, string filename, string format)
          {
            if (WriteUserFormat (format, self, *self.GetGeometry(), filename))
              {
                string err = string ("nothing known about format")+format;
                Array<const char*> names, extensions;
                RegisterUserFormats (names, extensions);
                err += "\navailable formats are:\n";
                for (auto name : names)
                  err += string("'") + name + "'\n";
                throw NgException (err);
              }
          },
         py::arg("filename"), py::arg("format"))
    
    .def_property("dim", &Mesh::GetDimension, &Mesh::SetDimension)

    .def("Elements3D", 
         static_cast<Array<Element>&(Mesh::*)()> (&Mesh::VolumeElements),
         py::return_value_policy::reference)

    .def("Elements2D", 
         static_cast<Array<Element2d>&(Mesh::*)()> (&Mesh::SurfaceElements),
         py::return_value_policy::reference)

    .def("Elements1D", 
         static_cast<Array<Segment>&(Mesh::*)()> (&Mesh::LineSegments),
         py::return_value_policy::reference)

    .def("Elements0D", FunctionPointer([] (Mesh & self) -> Array<Element0d>&
                                       {
                                         return self.pointelements;
                                       } ),
         py::return_value_policy::reference)

    .def("Points", 
         static_cast<Mesh::T_POINTS&(Mesh::*)()> (&Mesh::Points),
         py::return_value_policy::reference)

    .def("FaceDescriptor", static_cast<FaceDescriptor&(Mesh::*)(int)> (&Mesh::GetFaceDescriptor),
         py::return_value_policy::reference)
    .def("GetNFaceDescriptors", &Mesh::GetNFD)

    .def("GetNCD2Names", &Mesh::GetNCD2Names)
    

    .def("__getitem__", FunctionPointer ([](const Mesh & self, PointIndex pi)
                                         {
                                           return self[pi];
                                         }))

    .def ("Add", FunctionPointer ([](Mesh & self, MeshPoint p)
                                  {
                                    return self.AddPoint (Point3d(p));
                                  }))

    .def ("Add", FunctionPointer ([](Mesh & self, const Element & el)
                                  {
                                    return self.AddVolumeElement (el);
                                  }))

    .def ("Add", FunctionPointer ([](Mesh & self, const Element2d & el)
                                  {
                                    return self.AddSurfaceElement (el);
                                  }))

    .def ("Add", FunctionPointer ([](Mesh & self, const Segment & el)
                                  {
                                    return self.AddSegment (el);
                                  }))
    
    .def ("Add", FunctionPointer ([](Mesh & self, const Element0d & el)
                                  {
                                    return self.pointelements.Append (el);
                                  }))

    .def ("Add", FunctionPointer ([](Mesh & self, const FaceDescriptor & fd)
                                  {
                                    return self.AddFaceDescriptor (fd);
                                  }))

    .def ("SetBCName", &Mesh::SetBCName)
    .def ("GetBCName", FunctionPointer([](Mesh & self, int bc)->string 
                                       { return self.GetBCName(bc); }))
    .def ("SetMaterial", &Mesh::SetMaterial)
    .def ("GetMaterial", FunctionPointer([](Mesh & self, int domnr)
                                         { return string(self.GetMaterial(domnr)); }))

    .def ("GetCD2Name", &Mesh::GetCD2Name)
    .def("SetCD2Name", &Mesh::SetCD2Name)

    .def ("AddPointIdentification", [](Mesh & self, py::object pindex1, py::object pindex2, int identnr, int type)
                           {
			     if(py::extract<PointIndex>(pindex1).check() && py::extract<PointIndex>(pindex2).check())
			       {
				 self.GetIdentifications().Add (py::extract<PointIndex>(pindex1)(), py::extract<PointIndex>(pindex2)(), identnr);
				 self.GetIdentifications().SetType(identnr, Identifications::ID_TYPE(type)); // type = 2 ... periodic
			       }
                           },
          //py::default_call_policies(),
          py::arg("pid1"),
           py::arg("pid2"),
           py::arg("identnr"),
           py::arg("type"))
    .def ("GenerateVolumeMesh", 
          [](Mesh & self, py::object pymp)
           {
             cout << "generate vol mesh" << endl;

             MeshingParameters mp;
             if (py::extract<MeshingParameters>(pymp).check())
               mp = py::extract<MeshingParameters>(pymp)();
             else
               {
                 mp.optsteps3d = 5;
               }
             MeshVolume (mp, self);
             OptimizeVolume (mp, self);
           },
          py::arg("mp")=NGDummyArgument())

   .def ("OptimizeVolumeMesh", FunctionPointer
         ([](Mesh & self)
          {
            MeshingParameters mp;
            mp.optsteps3d = 5;
            OptimizeVolume (mp, self);
          }))

    .def ("Refine", FunctionPointer
          ([](Mesh & self)
           {
             if (self.GetGeometry())
               self.GetGeometry()->GetRefinement().Refine(self);
             else
               Refinement().Refine(self);
           }))

    .def ("SetGeometry", FunctionPointer
          ([](Mesh & self, shared_ptr<CSGeometry> geo)
           {
             self.SetGeometry(geo);
           }))

    // TODO: fix this dependency on libgeom2d.so
//     .def ("SetGeometry", FunctionPointer
//           ([](Mesh & self, shared_ptr<SplineGeometry2d> geo)
//            {
//              self.SetGeometry(geo);
//            }))

    .def ("BuildSearchTree", &Mesh::BuildElementSearchTree)

    .def ("BoundaryLayer", FunctionPointer 
          ([](Mesh & self, int bc, py::list thicknesses, int volnr, py::list materials)
           {
             int n = py::len(thicknesses);
             BoundaryLayerParameters blp;

             for (int i = 1; i <= self.GetNFD(); i++)
               if (self.GetFaceDescriptor(i).BCProperty() == bc)
                   blp.surfid.Append (i);

             cout << "add layer at surfaces: " << blp.surfid << endl;

             blp.prismlayers = n;
             blp.growthfactor = 1.0;

             // find max domain nr
             int maxind = 0;
             for (ElementIndex ei = 0; ei < self.GetNE(); ei++)
               maxind = max (maxind, self[ei].GetIndex());
             cout << "maxind = " << maxind << endl;
             for ( int i=0; i<n; i++ )
               {
                 blp.heights.Append( py::extract<double>(thicknesses[i])()) ;
                 blp.new_matnrs.Append( maxind+1+i );
                 self.SetMaterial (maxind+1+i, py::extract<string>(materials[i])().c_str());
               }
             blp.bulk_matnr = volnr;
             GenerateBoundaryLayer (self, blp);
           }
           ))

    .def ("BoundaryLayer", FunctionPointer
          ([](Mesh & self, int bc, double thickness, int volnr, string material)
           {
             BoundaryLayerParameters blp;

             for (int i = 1; i <= self.GetNFD(); i++)
               if (self.GetFaceDescriptor(i).BCProperty() == bc)
                   blp.surfid.Append (i);

             cout << "add layer at surfaces: " << blp.surfid << endl;

             blp.prismlayers = 1;
             blp.hfirst = thickness;
             blp.growthfactor = 1.0;

             // find max domain nr
             int maxind = 0;
             for (ElementIndex ei = 0; ei < self.GetNE(); ei++)
               maxind = max (maxind, self[ei].GetIndex());
             cout << "maxind = " << maxind << endl;
             self.SetMaterial (maxind+1, material.c_str());
             blp.new_matnr = maxind+1;
             blp.bulk_matnr = volnr;
             GenerateBoundaryLayer (self, blp);
           }
           ))

    .def ("Scale", FunctionPointer([](Mesh & self, double factor)
				   {
				     for(auto i = 0; i<self.GetNP();i++)
				       self.Point(i).Scale(factor);
				   }))
                                            
    ;
  

  typedef MeshingParameters MP;
  py::class_<MP> (m, "MeshingParameters")
    .def(py::init<>())
    .def("__init__",
         [](MP *instance, double maxh, bool quad_dominated, int optsteps2d, int optsteps3d)
                           {
                             new (instance) MeshingParameters;
                             instance->maxh = maxh;
                             instance->quad = int(quad_dominated);
                             instance->optsteps2d = optsteps2d;
                             instance->optsteps3d = optsteps3d;
                           },
           py::arg("maxh")=1000,
           py::arg("quad_dominated")=false,
           py::arg("optsteps2d") = 3,
           py::arg("optsteps3d") = 3
           ,
         "create meshing parameters"
          )
    .def("__str__", &ToString<MP>)
    .def_property("maxh", 
                  FunctionPointer ([](const MP & mp ) { return mp.maxh; }),
                  FunctionPointer ([](MP & mp, double maxh) { return mp.maxh = maxh; }))
    .def("RestrictH", FunctionPointer
         ([](MP & mp, double x, double y, double z, double h)
          {
            mp.meshsize_points.Append ( MeshingParameters::MeshSizePoint (Point<3> (x,y,z), h));
          }),
         py::arg("x"), py::arg("y"), py::arg("z"), py::arg("h")
         )
    ;

  m.def("SetTestoutFile", FunctionPointer ([] (const string & filename)
                                             {
                                               delete testout;
                                               testout = new ofstream (filename);
                                             }));

  m.def("SetMessageImportance", FunctionPointer ([] (int importance)
                                                   {
                                                     int old = printmessage_importance;
                                                     printmessage_importance = importance;
                                                     return old;
                                                   }));
}

PYBIND11_PLUGIN(libmesh) {
  py::module m("mesh", "pybind mesh");
  ExportNetgenMeshing(m);
  return m.ptr();
}
#endif




