//  GEOM PARTITION : partition algorithm
//
//  Copyright (C) 2003  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS 
// 
//  This library is free software; you can redistribute it and/or 
//  modify it under the terms of the GNU Lesser General Public 
//  License as published by the Free Software Foundation; either 
//  version 2.1 of the License. 
// 
//  This library is distributed in the hope that it will be useful, 
//  but WITHOUT ANY WARRANTY; without even the implied warranty of 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
//  Lesser General Public License for more details. 
// 
//  You should have received a copy of the GNU Lesser General Public 
//  License along with this library; if not, write to the Free Software 
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA 
// 
//  See http://www.opencascade.org/SALOME/ or email : webmaster.salome@opencascade.org 
//
//
//
//  File   : Partition_Inter2d.hxx
//  Module : GEOM

#ifndef _Partition_Inter2d_HeaderFile
#define _Partition_Inter2d_HeaderFile

#ifndef _Handle_BRepAlgo_AsDes_HeaderFile
#include <Handle_BRepAlgo_AsDes.hxx>
#endif
#ifndef _Standard_Real_HeaderFile
#include <Standard_Real.hxx>
#endif
#ifndef _Standard_Boolean_HeaderFile
#include <Standard_Boolean.hxx>
#endif
class BRepAlgo_AsDes;
class TopoDS_Face;
class TopTools_MapOfShape;
class TopoDS_Vertex;
class TopTools_ListOfShape;
class TopoDS_Edge;


#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif

class Partition_Inter2d  {

public:

   void* operator new(size_t,void* anAddress) 
   {
      return anAddress;
   }
   void* operator new(size_t size) 
   { 
      return Standard::Allocate(size); 
   }
   void  operator delete(void *anAddress) 
   { 
      if (anAddress) Standard::Free((Standard_Address&)anAddress); 
   }
   // Methods PUBLIC
   // 
   static  void CompletPart2d(const Handle(BRepAlgo_AsDes)& AsDes,const TopoDS_Face& F,const TopTools_MapOfShape& NewEdges) ;
   static  TopoDS_Vertex FindEndVertex(const TopTools_ListOfShape& VertList,const Standard_Real f,const Standard_Real l,const TopoDS_Edge& E,Standard_Boolean& First,Standard_Real& DU) ;
   static  TopoDS_Vertex AddVonE(const TopoDS_Vertex& V,const TopoDS_Edge& E1,const TopoDS_Edge& E2,const Handle(BRepAlgo_AsDes)& AsDes,const TopoDS_Face& F) ;
   static  Standard_Real GetTolerance(const TopoDS_Vertex& theV,const Standard_Real theU,const TopoDS_Edge& theE,const Handle(BRepAlgo_AsDes)& theAsDes) ;




protected:

   // Methods PROTECTED
   // 


   // Fields PROTECTED
   //


private: 

   // Methods PRIVATE
   // 


   // Fields PRIVATE
   //


};





// other Inline functions and methods (like "C++: function call" methods)
//


#endif
