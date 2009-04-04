//  GEOM PARTITION : partition algorithm
//
//  Copyright (C) 2003  CEA/DEN, EDF R&D
//
//
//
//  File   : Partition_Spliter.hxx
//  Module : GEOM

#ifndef _Partition_Spliter_HeaderFile
#define _Partition_Spliter_HeaderFile

#ifndef _TopAbs_ShapeEnum_HeaderFile
#include <TopAbs_ShapeEnum.hxx>
#endif
#ifndef _TopoDS_Compound_HeaderFile
#include <TopoDS_Compound.hxx>
#endif
#ifndef _BRep_Builder_HeaderFile
#include <BRep_Builder.hxx>
#endif
#ifndef _TopTools_ListOfShape_HeaderFile
#include <TopTools_ListOfShape.hxx>
#endif
#ifndef _TopTools_MapOfShape_HeaderFile
#include <TopTools_MapOfShape.hxx>
#endif
#ifndef _TopTools_DataMapOfShapeShape_HeaderFile
#include <TopTools_DataMapOfShapeShape.hxx>
#endif
#ifndef _Handle_BRepAlgo_AsDes_HeaderFile
#include <Handle_BRepAlgo_AsDes.hxx>
#endif
#ifndef _BRepAlgo_Image_HeaderFile
#include <BRepAlgo_Image.hxx>
#endif
#ifndef _Partition_Inter3d_HeaderFile
#include "Partition_Inter3d.hxx"
#endif
#ifndef _TopTools_MapOfOrientedShape_HeaderFile
#include <TopTools_MapOfOrientedShape.hxx>
#endif
#ifndef _Standard_Boolean_HeaderFile
#include <Standard_Boolean.hxx>
#endif
class BRepAlgo_AsDes;
class TopoDS_Shape;
class TopTools_ListOfShape;
class TopoDS_Edge;


#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif

class Partition_Spliter  {

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
   Partition_Spliter();
   void AddShape(const TopoDS_Shape& S) ;
   void AddTool(const TopoDS_Shape& S) ;
   void Compute(const TopAbs_ShapeEnum Limit = TopAbs_SHAPE) ;
   void KeepShapesInside(const TopoDS_Shape& S) ;
   void RemoveShapesInside(const TopoDS_Shape& S) ;
   TopoDS_Shape Shape() const;
   void Clear() ;





protected:

   // Methods PROTECTED
   // 


   // Fields PROTECTED
   //


private: 

   // Methods PRIVATE
   // 
   void MakeSolids(const TopoDS_Shape& Solid,TopTools_ListOfShape& Shells) ;
   void MakeShells(const TopoDS_Shape& S,TopTools_ListOfShape& NS) ;
   TopoDS_Shape MakeFaces(const TopoDS_Shape& S) ;
   void MakeEdges(const TopoDS_Edge& E,const TopTools_ListOfShape& VOnE,TopTools_ListOfShape& NE) const;
   TopoDS_Shape FindFacesInside(const TopoDS_Shape& S,const Standard_Boolean CheckClosed = Standard_False,const Standard_Boolean All = Standard_False) ;
   Standard_Boolean CheckTool(const TopoDS_Shape& S) ;
   void MergeEqualEdges(const TopTools_ListOfShape& LE) ;
   static  Standard_Boolean IsInside(const TopoDS_Shape& S1,const TopoDS_Shape& S2) ;
   TopoDS_Shape GetOriginalShape(const TopoDS_Shape& aShape) const;
   void FindToolsToReconstruct() ;


   // Fields PRIVATE
   //
   TopAbs_ShapeEnum myDoneStep;
   TopoDS_Compound myShape;
   BRep_Builder myBuilder;
   TopTools_ListOfShape myListShapes;
   TopTools_MapOfShape myMapFaces;
   TopTools_MapOfShape myMapTools;
   TopTools_MapOfShape myEqualEdges;
   TopTools_MapOfShape myNewSection;
   TopTools_MapOfShape myClosedShapes;
   TopTools_MapOfShape mySharedFaces;
   TopTools_MapOfShape myWrappingSolid;
   TopTools_DataMapOfShapeShape myFaceShapeMap;
   TopTools_DataMapOfShapeShape myInternalFaces;
   TopTools_DataMapOfShapeShape myIntNotClFaces;
   Handle_BRepAlgo_AsDes myAsDes;
   BRepAlgo_Image myImagesFaces;
   BRepAlgo_Image myImagesEdges;
   BRepAlgo_Image myImageShape;
   Partition_Inter3d myInter3d;
   TopTools_MapOfOrientedShape myAddedFacesMap;


};





// other Inline functions and methods (like "C++: function call" methods)
//


#endif
