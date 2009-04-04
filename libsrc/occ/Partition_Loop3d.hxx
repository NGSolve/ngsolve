//  GEOM PARTITION : partition algorithm
//
//  Copyright (C) 2003  CEA/DEN, EDF R&D
//
//
//
//  File   : Partition_Loop3d.hxx
//  Module : GEOM

#ifndef _Partition_Loop3d_HeaderFile
#define _Partition_Loop3d_HeaderFile

#ifndef _TopTools_ListOfShape_HeaderFile
#include <TopTools_ListOfShape.hxx>
#endif
#ifndef _TopTools_IndexedDataMapOfShapeListOfShape_HeaderFile
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#endif
#ifndef _Standard_Boolean_HeaderFile
#include <Standard_Boolean.hxx>
#endif
#ifndef _Standard_Real_HeaderFile
#include <Standard_Real.hxx>
#endif
class TopoDS_Shape;
class TopTools_ListOfShape;
class TopTools_MapOfOrientedShape;
class TopoDS_Edge;
class TopoDS_Face;
class gp_Vec;


#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif

class Partition_Loop3d  {

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
   Partition_Loop3d();
   void AddConstFaces(const TopoDS_Shape& S) ;
   void AddSectionFaces(const TopoDS_Shape& S) ;
   const TopTools_ListOfShape& MakeShells(const TopTools_MapOfOrientedShape& AvoidFacesMap) ;
   static  Standard_Boolean IsInside(const TopoDS_Edge& E,const TopoDS_Face& F1,const TopoDS_Face& F2,const Standard_Boolean CountDot,Standard_Real& Dot,Standard_Boolean& GoodOri) ;
   static  gp_Vec Normal(const TopoDS_Edge& E,const TopoDS_Face& F) ;





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
   TopTools_ListOfShape myNewShells;
   TopTools_ListOfShape myFaces;
   TopTools_IndexedDataMapOfShapeListOfShape myEFMap;


};





// other Inline functions and methods (like "C++: function call" methods)
//


#endif
