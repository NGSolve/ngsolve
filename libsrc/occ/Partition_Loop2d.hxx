//  GEOM PARTITION : partition algorithm
//
//  Copyright (C) 2003  CEA/DEN, EDF R&D
//
//
//
//  File   : Partition_Loop2d.hxx
//  Module : GEOM

#ifndef _Partition_Loop2d_HeaderFile
#define _Partition_Loop2d_HeaderFile

#ifndef _TopoDS_Face_HeaderFile
#include <TopoDS_Face.hxx>
#endif
#ifndef _TopAbs_Orientation_HeaderFile
#include <TopAbs_Orientation.hxx>
#endif
#ifndef _TopTools_ListOfShape_HeaderFile
#include <TopTools_ListOfShape.hxx>
#endif
#ifndef _TopTools_MapOfShape_HeaderFile
#include <TopTools_MapOfShape.hxx>
#endif
class TopoDS_Face;
class TopoDS_Edge;
class TopTools_ListOfShape;
class BRepAlgo_Image;


#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif

class Partition_Loop2d  {

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
   Partition_Loop2d();
   void Init(const TopoDS_Face& F) ;
   void AddConstEdge(const TopoDS_Edge& E) ;
   void AddSectionEdge(const TopoDS_Edge& E) ;
   void Perform() ;
   const TopTools_ListOfShape& NewWires() const;
   void WiresToFaces(const BRepAlgo_Image& EdgeImage) ;
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
   TopAbs_Orientation myFaceOri;
   TopTools_ListOfShape myConstEdges;
   TopTools_ListOfShape myNewWires;
   TopTools_ListOfShape myNewFaces;
   TopTools_ListOfShape myInternalWL;
   TopTools_MapOfShape mySectionEdges;


};





// other Inline functions and methods (like "C++: function call" methods)
//


#endif
