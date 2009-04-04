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
//  File   : Partition_Inter3d.hxx
//  Module : GEOM

#ifndef _Partition_Inter3d_HeaderFile
#define _Partition_Inter3d_HeaderFile

#ifndef _Handle_BRepAlgo_AsDes_HeaderFile
#include <Handle_BRepAlgo_AsDes.hxx>
#endif
#ifndef _TopTools_DataMapOfShapeListOfShape_HeaderFile
#include <TopTools_DataMapOfShapeListOfShape.hxx>
#endif
#ifndef _TopTools_MapOfShape_HeaderFile
#include <TopTools_MapOfShape.hxx>
#endif
#ifndef _TopTools_DataMapOfShapeShape_HeaderFile
#include <TopTools_DataMapOfShapeShape.hxx>
#endif
#ifndef _Standard_Boolean_HeaderFile
#include <Standard_Boolean.hxx>
#endif
class BRepAlgo_AsDes;
class TopTools_ListOfShape;
class TopTools_DataMapOfShapeShape;
class TopoDS_Face;
class TopTools_MapOfShape;
class TopoDS_Shape;
class TopoDS_Vertex;
class TopoDS_Edge;


#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif

class Partition_Inter3d  {

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
   Partition_Inter3d();
   Partition_Inter3d(const Handle(BRepAlgo_AsDes)& AsDes);
   void CompletPart3d(const TopTools_ListOfShape& SetOfFaces1,const TopTools_DataMapOfShapeShape& FaceShapeMap) ;
   void FacesPartition(const TopoDS_Face& F1,const TopoDS_Face& F2) ;
   Standard_Boolean IsDone(const TopoDS_Face& F1,const TopoDS_Face& F2) const;
   TopTools_MapOfShape& TouchedFaces() ;
   Handle_BRepAlgo_AsDes AsDes() const;
   TopTools_MapOfShape& NewEdges() ;
   Standard_Boolean HasSameDomainF(const TopoDS_Shape& F) const;
   Standard_Boolean IsSameDomainF(const TopoDS_Shape& F1,const TopoDS_Shape& F2) const;
   const TopTools_ListOfShape& SameDomain(const TopoDS_Face& F) const;
   TopoDS_Vertex ReplaceSameDomainV(const TopoDS_Vertex& V,const TopoDS_Edge& E) const;
   Handle_BRepAlgo_AsDes SectionEdgesAD() const;
   Standard_Boolean IsSectionEdge(const TopoDS_Edge& E) const;
   Standard_Boolean HasSectionEdge(const TopoDS_Face& F) const;
   Standard_Boolean IsSplitOn(const TopoDS_Edge& NewE,const TopoDS_Edge& OldE,const TopoDS_Face& F) const;
   const TopTools_ListOfShape& SectionEdgeFaces(const TopoDS_Edge& SecE) const;





protected:

   // Methods PROTECTED
   // 


   // Fields PROTECTED
   //


private: 

   // Methods PRIVATE
   // 
   void Inter3D(const TopoDS_Face& F1,const TopoDS_Face& F2,TopTools_ListOfShape& LInt) ;
   void StorePart3d(const TopoDS_Face& F1,const TopoDS_Face& F2,const TopTools_ListOfShape& LInt1) ;
   void SetDone(const TopoDS_Face& F1,const TopoDS_Face& F2) ;
   void Affiche(const TopTools_ListOfShape& SetOfFaces) const;


   // Fields PRIVATE
   //
   Handle_BRepAlgo_AsDes myAsDes;
   TopTools_DataMapOfShapeListOfShape myDone;
   TopTools_MapOfShape myTouched;
   TopTools_MapOfShape myNewEdges;
   Handle_BRepAlgo_AsDes mySectionEdgesAD;
   TopTools_DataMapOfShapeListOfShape mySameDomainFM;
   TopTools_DataMapOfShapeShape mySameDomainVM;


};





// other Inline functions and methods (like "C++: function call" methods)
//


#endif
