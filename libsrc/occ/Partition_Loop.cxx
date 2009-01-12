#ifdef OCCGEOMETRY

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
//  File   : Partition_Loop.cxx
//  Author : Benedicte MARTIN
//  Module : GEOM
//  $Header: /cvs/netgen/netgen/libsrc/occ/Partition_Loop.cxx,v 1.6 2008/03/31 14:20:28 wabro Exp $

//using namespace std;
#include <cstdio>
#include <climits>

#include "Partition_Loop.ixx"

#include "utilities.h"

#include <BRep_Builder.hxx>
#include <BRepAlgo_FaceRestrictor.hxx>
#include <BRep_Tool.hxx>

#include <Geom2d_Curve.hxx>
#include <Geom_Surface.hxx>

#include <TopTools_SequenceOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopTools_MapIteratorOfMapOfShape.hxx>
#include <TopTools_MapOfOrientedShape.hxx>
#include <TopTools_DataMapOfShapeShape.hxx>
#include <TopTools_DataMapIteratorOfDataMapOfShapeListOfShape.hxx>

#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Iterator.hxx>

#include <Precision.hxx>
#include <BRep_TVertex.hxx>
#include <BRep_TEdge.hxx>

#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>

static char* name = new char[100];
static int nbe = 0;

//=======================================================================
//function : Partition_Loop
//purpose  : 
//=======================================================================
Partition_Loop::Partition_Loop()
{
}

//=======================================================================
//function : Init
//purpose  : 
//=======================================================================
void Partition_Loop::Init(const TopoDS_Face& F)
{
  myConstEdges.Clear(); 
  myNewWires  .Clear();
  myNewFaces  .Clear();
  myFace = F;
}

//=======================================================================
//function : AddConstEdge
//purpose  : 
//=======================================================================
void Partition_Loop::AddConstEdge (const TopoDS_Edge& E)
{
  myConstEdges.Append(E);
}


//=======================================================================
//function : FindDelta
//purpose  : 
//=======================================================================
static Standard_Real FindDelta(TopTools_ListOfShape& LE,
			       const TopoDS_Face& F)
{
  Standard_Real dist, f, l;
  Standard_Real d = Precision::Infinite();
  TopTools_ListIteratorOfListOfShape itl;

  for ( itl.Initialize(LE); itl.More(); itl.Next()) {
    const TopoDS_Edge& E = TopoDS::Edge(itl.Value());
    Handle(Geom2d_Curve) C = BRep_Tool::CurveOnSurface(E,F,f,l);
    gp_Pnt2d p = C->Value(f);
    gp_Pnt2d pp = C->Value(l);
    Standard_Real d1 = p.Distance(pp);
    if (d1<d) { d=d1;}
  }
  dist = d ;
  return dist;
}

//=======================================================================
//function : SelectEdge
//purpose  : Find the edge <NE> connected <CE> by the vertex <CV> in the list <LE>.
//           <NE> Is erased  of the list. If <CE> is too in the list <LE> 
//			 with the same orientation, it's erased of the list 
//=======================================================================
static Standard_Boolean  SelectEdge(const TopoDS_Face&    F,
				    const TopoDS_Edge&    CE,
				    const TopoDS_Vertex&  CV,
				    TopoDS_Edge&          NE,
				    TopTools_ListOfShape& LE)
{
  TopTools_ListIteratorOfListOfShape itl;
  NE.Nullify();
  for ( itl.Initialize(LE); itl.More(); itl.Next()) {
    if (itl.Value().IsEqual(CE)) {
      LE.Remove(itl);
      break;
    }
  }

  if (LE.Extent() > 1) {
    //--------------------------------------------------------------
    // Several possible edges.   
    // - Test the edges differents of CE 
    //--------------------------------------------------------------
    Standard_Real   cf, cl, f, l;
    TopoDS_Face FForward = F;
    Handle(Geom2d_Curve) Cc, C;
    FForward.Orientation(TopAbs_FORWARD);
			
    Cc = BRep_Tool::CurveOnSurface(CE,FForward,cf,cl);
    Standard_Real dist,distmin  = 100*BRep_Tool::Tolerance(CV);
    Standard_Real uc,u;
    if (CE.Orientation () == TopAbs_FORWARD) uc = cl;
    else                                     uc = cf;

    gp_Pnt2d P2,PV = Cc->Value(uc); 

    Standard_Real delta = FindDelta(LE,FForward);

    for ( itl.Initialize(LE); itl.More(); itl.Next()) {
      const TopoDS_Edge& E = TopoDS::Edge(itl.Value());
      if (!E.IsSame(CE)) {
	C = BRep_Tool::CurveOnSurface(E,FForward,f,l);
	if (E.Orientation () == TopAbs_FORWARD) u = f;
	else                                    u = l;
	P2 = C->Value(u);
	dist = PV.Distance(P2);
	if (dist <= distmin){
	  distmin = dist;
	}
				
      }
    }

    Standard_Real anglemax = - PI;
    TopoDS_Edge   SelectedEdge;	
    for ( itl.Initialize(LE); itl.More(); itl.Next()) {
      const TopoDS_Edge& E = TopoDS::Edge(itl.Value());
      if (!E.IsSame(CE)) {
	C = BRep_Tool::CurveOnSurface(E,FForward,f,l);
	if (E.Orientation () == TopAbs_FORWARD) u = f;
	else                                    u = l;
	P2 = C->Value(u);
	dist = PV.Distance(P2);
	if (dist <= distmin + (1./3)*delta){ 
	  gp_Pnt2d PC, P;
	  gp_Vec2d CTg1, CTg2, Tg1, Tg2;
	  Cc->D2(uc, PC, CTg1, CTg2);
	  C->D2(u, P, Tg1, Tg2);

	  Standard_Real angle;

	  if (CE.Orientation () == TopAbs_REVERSED && E.Orientation () == TopAbs_FORWARD) {
	    angle = CTg1.Angle(Tg1.Reversed());
	  }
	  else if (CE.Orientation () == TopAbs_FORWARD && E.Orientation () == TopAbs_REVERSED) {
	    angle = (CTg1.Reversed()).Angle(Tg1);
	  }
	  else if (CE.Orientation () == TopAbs_REVERSED && E.Orientation () == TopAbs_REVERSED) {
	    angle = CTg1.Angle(Tg1);
	  }
	  else if (CE.Orientation () == TopAbs_FORWARD && E.Orientation () == TopAbs_FORWARD) {
	    angle = (CTg1.Reversed()).Angle(Tg1.Reversed());
	  }
	  if (angle >= anglemax) {
	    anglemax = angle ;
	    SelectedEdge = E;	
	  }
	}
      }
    }
    for ( itl.Initialize(LE); itl.More(); itl.Next()) {
      const TopoDS_Edge& E = TopoDS::Edge(itl.Value());
      if (E.IsEqual(SelectedEdge)) {
	NE = TopoDS::Edge(E);
	LE.Remove(itl);
	break;
      }
    }					
  }
  else if (LE.Extent() == 1) {
    NE = TopoDS::Edge(LE.First());
    LE.RemoveFirst();
  }
  else {
    return Standard_False;
  }
  return Standard_True;
}

//=======================================================================
//function : SamePnt2d
//purpose  : 
//=======================================================================
static Standard_Boolean  SamePnt2d(TopoDS_Vertex  V,
				   TopoDS_Edge&   E1,
				   TopoDS_Edge&   E2,
				   TopoDS_Face&   F)
{
  Standard_Real   f1,f2,l1,l2;
  gp_Pnt2d        P1,P2;
  TopoDS_Shape aLocalF = F.Oriented(TopAbs_FORWARD);
  TopoDS_Face FF = TopoDS::Face(aLocalF);
  Handle(Geom2d_Curve) C1 = BRep_Tool::CurveOnSurface(E1,FF,f1,l1);  
  Handle(Geom2d_Curve) C2 = BRep_Tool::CurveOnSurface(E2,FF,f2,l2);  
  if (E1.Orientation () == TopAbs_FORWARD) P1 = C1->Value(f1);
  else                                     P1 = C1->Value(l1);
  
  if (E2.Orientation () == TopAbs_FORWARD) P2 = C2->Value(l2);
  else                                     P2 = C2->Value(f2);
  Standard_Real Tol  = 100*BRep_Tool::Tolerance(V);
  Standard_Real Dist = P1.Distance(P2);
  return Dist < Tol; 
}

//=======================================================================
//function : PurgeNewEdges
//purpose  : 
//=======================================================================
static void  PurgeNewEdges(TopTools_ListOfShape& ConstEdges,
			   const TopTools_MapOfOrientedShape&          UsedEdges)
{
  TopTools_ListIteratorOfListOfShape it(ConstEdges);
  while ( it.More()) {
    const TopoDS_Shape& NE = it.Value();
    if (!UsedEdges.Contains(NE)) {
      ConstEdges.Remove(it);
    }
    else {
      it.Next();
    }
  }  
}

//=======================================================================
//function : StoreInMVE
//purpose  : 
//=======================================================================
static void StoreInMVE (const TopoDS_Face& F,
			TopoDS_Edge& E,
			TopTools_DataMapOfShapeListOfShape& MVE )

{ 
  TopoDS_Vertex V1, V2;
  TopTools_ListOfShape Empty;

  TopExp::Vertices(E,V1,V2);
  if (!MVE.IsBound(V1)) {
    MVE.Bind(V1,Empty);
  }
  MVE(V1).Append(E);
	
  if (!MVE.IsBound(V2)) {
    MVE.Bind(V2,Empty);
  }
  MVE(V2).Append(E);
}

//=======================================================================
//function : Perform
//purpose  : 
//=======================================================================
void Partition_Loop::Perform()
{

  TopTools_DataMapOfShapeListOfShape MVE;
  TopTools_DataMapIteratorOfDataMapOfShapeListOfShape Mapit, Mapit1;  
  TopTools_ListIteratorOfListOfShape                  itl;
  TopoDS_Vertex                                       V1,V2;

  //-----------------------------------
  // Construction map vertex => edges
  //-----------------------------------
  for (itl.Initialize(myConstEdges); itl.More(); itl.Next()) {
    TopoDS_Edge& E = TopoDS::Edge(itl.Value());
    StoreInMVE(myFace,E,MVE);
  }

  //----------------------------------------------
  // Construction of all the wires and of all the new faces. 
  //----------------------------------------------
  TopTools_MapOfOrientedShape UsedEdges;

  while (!MVE.IsEmpty()) {
    TopoDS_Vertex    VF,CV;
    TopoDS_Edge      CE,NE,EF;
    TopoDS_Wire      NW;
    BRep_Builder     B;
    Standard_Boolean End= Standard_False;

    B.MakeWire(NW);
    //--------------------------------
    // EF first edge.
    //--------------------------------
    Mapit.Initialize(MVE);
    EF = CE = TopoDS::Edge(Mapit.Value().First());

    TopExp::Vertices(CE,V1,V2);
    //--------------------------------
    // VF first vertex 
    //--------------------------------
    if (CE.Orientation() == TopAbs_FORWARD) { 
      CV = VF = V1;
    }
    else  { 
      CV = VF = V2;
    }
    if (!MVE.IsBound(CV)) continue;
    for ( itl.Initialize(MVE(CV)); itl.More(); itl.Next()) {
      if (itl.Value().IsEqual(CE)) {
	MVE(CV).Remove(itl);
	break;
      }
    }

    int i = 0;
    while (!End) { 
      //-------------------------------
      // Construction of a wire.
      //-------------------------------
      TopExp::Vertices(CE,V1,V2);
      if (!CV.IsSame(V1)) CV = V1; else CV = V2; 
      B.Add (NW,CE);
      UsedEdges.Add(CE);

      //--------------
      // stop test
      //--------------			
      if (!MVE.IsBound(CV) || MVE(CV).IsEmpty() || CV.IsSame(VF) ) {
	if (CV.IsSame(VF)) {
	  if (MVE(CV).Extent() == 1 ) MVE.UnBind(CV);
	  else {
	    for ( itl.Initialize(MVE(CV)); itl.More(); itl.Next()) {
	      if (itl.Value().IsEqual(CE)) {
		MVE(CV).Remove(itl);
		break;
	      }
	    }
	  }
	}
	End=Standard_True;
      } 

      //--------------
      // select edge
      //--------------
      else {
	Standard_Boolean find = SelectEdge(myFace,CE,CV,NE,MVE(CV));
	if (find) {
	  CE=NE;
	  if (MVE(CV).IsEmpty()) MVE.UnBind(CV);
	  if (CE.IsNull() ) {
	    MESSAGE ( " CE is  NULL !!! " )
	    End=Standard_True;
	  }
	}
	else {
	  MESSAGE ( " edge doesn't exist " )
	  End=Standard_True;
	}
      }
    }

    //-----------------------------
    // Test if the wire is closed  
    //-----------------------------
    if (VF.IsSame(CV) && SamePnt2d(VF,EF,CE,myFace)) {
    }
    else{
      MESSAGE ( "wire not closed" )
    }
    myNewWires.Append (NW);			
  }

  PurgeNewEdges(myConstEdges,UsedEdges);

}


//=======================================================================
//function : NewWires
//purpose  : 
//=======================================================================
const TopTools_ListOfShape&  Partition_Loop::NewWires() const 
{  
  return myNewWires;
}

//=======================================================================
//function : NewFaces
//purpose  : 
//=======================================================================
const TopTools_ListOfShape&  Partition_Loop::NewFaces() const 
{  
  return myNewFaces;
}
 
//=======================================================================
//function : WiresToFaces
//purpose  : 
//=======================================================================
void  Partition_Loop::WiresToFaces() 
{  
  if (!myNewWires.IsEmpty()) {
    BRepAlgo_FaceRestrictor FR;

    TopAbs_Orientation OriF = myFace.Orientation();
    TopoDS_Shape aLocalS = myFace.Oriented(TopAbs_FORWARD);

    FR.Init (TopoDS::Face(aLocalS),Standard_False);
    TopTools_ListIteratorOfListOfShape it(myNewWires);
    for (; it.More(); it.Next()) {
      FR.Add(TopoDS::Wire(it.Value()));
    }

    FR.Perform();
    
    if (FR.IsDone()) {
      for (; FR.More(); FR.Next()) {
	myNewFaces.Append(FR.Current().Oriented(OriF));
      }
    }
  }
}


#endif
