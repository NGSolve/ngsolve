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
//  File   : Partition_Loop.hxx
//  Module : GEOM

#ifndef _Partition_Loop_HeaderFile
#define _Partition_Loop_HeaderFile

#ifndef _TopoDS_Face_HeaderFile
#include <TopoDS_Face.hxx>
#endif
#ifndef _TopTools_ListOfShape_HeaderFile
#include <TopTools_ListOfShape.hxx>
#endif
#ifndef _TopTools_DataMapOfShapeListOfShape_HeaderFile
#include <TopTools_DataMapOfShapeListOfShape.hxx>
#endif
class TopoDS_Face;
class TopoDS_Edge;
class TopTools_ListOfShape;


#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif

class Partition_Loop  {

public:

   inline void* operator new(size_t,void* anAddress) 
   {
      return anAddress;
   }
   inline void* operator new(size_t size) 
   { 
      return Standard::Allocate(size); 
   }
   inline void  operator delete(void *anAddress) 
   { 
      if (anAddress) Standard::Free((Standard_Address&)anAddress); 
   }
   //    inline void  operator delete(void *anAddress, size_t size) 
   //      { 
   //        if (anAddress) Standard::Free((Standard_Address&)anAddress,size); 
   //      }
   // Methods PUBLIC
   // 
   Partition_Loop();
   void Init(const TopoDS_Face& F) ;
   void AddConstEdge(const TopoDS_Edge& E) ;
   void Perform() ;
   const TopTools_ListOfShape& NewWires() const;
   void WiresToFaces() ;
   const TopTools_ListOfShape& NewFaces() const;




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
   TopoDS_Face myFace;
   TopTools_ListOfShape myConstEdges;
   TopTools_ListOfShape myNewWires;
   TopTools_ListOfShape myNewFaces;


};





// other inline functions and methods (like "C++: function call" methods)
//


#endif
