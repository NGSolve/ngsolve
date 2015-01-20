#ifdef NG_PYTHON

#include <boost/python.hpp>
#include <boost/python/slice.hpp>
#include <../general/ngpython.hpp>

#include <mystdlib.h>
#include "meshing.hpp"


using namespace netgen;
namespace bp = boost::python;

namespace netgen
{
  extern shared_ptr<NetgenGeometry> ng_geometry;
}


template <typename T, int BASE = 0, typename TIND = int>
void ExportArray ()
{
  string name = string("Array_") + typeid(T).name();
  bp::class_<Array<T,BASE,TIND>,boost::noncopyable>(name.c_str())
    .def ("__len__", &Array<T,BASE,TIND>::Size)
    .def ("__getitem__", 
          FunctionPointer ([](Array<T,BASE,TIND> & self, TIND i) -> T&
                           {
                             if (i < BASE || i >= BASE+self.Size())
                               bp::exec("raise IndexError()\n");
                             return self[i];
                           }),
          bp::return_value_policy<bp::reference_existing_object>())

    .def ("__iter__", 
          bp::range (FunctionPointer([](Array<T,BASE,TIND> & self) { return &self[BASE]; }),
                     FunctionPointer([](Array<T,BASE,TIND> & self) { return &self[BASE+self.Size()]; })))

    ;
}



DLL_HEADER void ExportNetgenMeshing() 
{
  
  ModuleScope module("meshing");

  bp::class_<PointIndex>("PointId", bp::init<int>())
    .def("__repr__", &ToString<PointIndex>)
    .def("__str__", &ToString<PointIndex>)
    .add_property("nr", &PointIndex::operator int)
    ;

  /*  
  bp::class_<Point<3>> ("Point")
    .def(bp::init<double,double,double>())
    ;
  */

  bp::class_<MeshPoint /* ,bp::bases<Point<3>> */ >("MeshPoint")
    // .def(bp::init<Point<3>>())
    .add_property("p", FunctionPointer([](const MeshPoint & self)
                                       {
                                         bp::list l;
                                         l.append ( self[0] );
                                         l.append ( self[1] );
                                         l.append ( self[2] );
                                         return bp::tuple(l);
                                       }))
    ;
  
  bp::class_<Element>("Element3D")
    .add_property("index", &Element::GetIndex, &Element::SetIndex)
    .add_property("vertices", 
                  FunctionPointer ([](const Element & self) -> bp::list
                                   {
                                     bp::list li;
                                     for (int i = 0; i < self.GetNV(); i++)
                                       li.append (self[i]);
                                     return li;
                                   }))
    ;

  bp::class_<Element2d>("Element2D")
    .add_property("index", &Element2d::GetIndex, &Element2d::SetIndex)
    .add_property("vertices",
                  FunctionPointer([](const Element2d & self) -> bp::list
                                  {
                                    bp::list li;
                                    for (int i = 0; i < self.GetNV(); i++)
                                      li.append(self[i]);
                                    return li;
                                  }))
    ;
  ExportArray<Element>();
  ExportArray<Element2d>();
  ExportArray<MeshPoint,PointIndex::BASE,PointIndex>();
  ;
  
  
  bp::class_<Mesh,shared_ptr<Mesh>,boost::noncopyable>("Mesh", bp::no_init)
    .def(bp::init<>("create empty mesh"))
    .def("__str__", &ToString<Mesh>)
    .def("Load",  FunctionPointer 
	 ([](Mesh & self, const string & filename)
	  {
	    ifstream input(filename);
	    self.Load(input);
	    for (int i = 0; i < geometryregister.Size(); i++)
	      {
		NetgenGeometry * hgeom = geometryregister[i]->LoadFromMeshFile (input);
		if (hgeom)
		  {
		    ng_geometry.reset (hgeom);
		    break;
		  }
	      }
	  }))
    // static_cast<void(Mesh::*)(const string & name)>(&Mesh::Load))
    .def("Save", static_cast<void(Mesh::*)(const string & name)const>(&Mesh::Save))

    .def("Elements3D", 
         static_cast<Array<Element>&(Mesh::*)()> (&Mesh::VolumeElements),
         bp::return_value_policy<bp::reference_existing_object>())

    .def("Elements2D", 
         static_cast<Array<Element2d>&(Mesh::*)()> (&Mesh::SurfaceElements),
         bp::return_value_policy<bp::reference_existing_object>())

    .def("Points", 
         static_cast<Mesh::T_POINTS&(Mesh::*)()> (&Mesh::Points),
         bp::return_value_policy<bp::reference_existing_object>())


    .def("__getitem__", FunctionPointer ([](const Mesh & self, PointIndex pi)
                                         {
                                           return self[pi];
                                         }))

    .def ("Add", FunctionPointer ([](Mesh & self, MeshPoint p)
                                  {
                                    return self.AddPoint (Point3d(p));
                                  }))
    /*
    .def("__init__", bp::make_constructor
         (FunctionPointer ([]()
                           {
                             cout << "create new mesh" << endl;
                             auto tmp = make_shared<Mesh>();
                             return tmp;
                           })),
         "create empty mesh"
      )
    */
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
                                            
    ;
  

  typedef MeshingParameters MP;
  bp::class_<MP> ("MeshingParameters", bp::init<>())
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](double maxh)
                           {
                             auto tmp = new MeshingParameters;
                             tmp->maxh = maxh;
                             return tmp;
                           }),
          bp::default_call_policies(),        // need it to use arguments
          (bp::arg("maxh")=1000)),
         "create meshing parameters"
         )
    .def("__str__", &ToString<MP>)
    .add_property("maxh", 
                  FunctionPointer ([](const MP & mp ) { return mp.maxh; }),
                  FunctionPointer ([](MP & mp, double maxh) { return mp.maxh = maxh; }))
                  
    ;

  bp::def("SetTestoutFile", FunctionPointer ([] (const string & filename)
                                             {
                                               delete testout;
                                               testout = new ofstream (filename);
                                             }));

  bp::def("SetMessageImportance", FunctionPointer ([] (int importance)
                                                   {
                                                     int old = printmessage_importance;
                                                     printmessage_importance = importance;
                                                     return old;
                                                   }));
}



BOOST_PYTHON_MODULE(libmesh) {
  ExportNetgenMeshing();
}



#endif




