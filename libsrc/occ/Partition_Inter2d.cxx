#ifdef OCCGEOMETRY

//  GEOM PARTITION : partition algorithm
//
//  Copyright (C) 2003  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R& D, LEG, PRINCIPIA R& D, BUREAU VERITAS
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
//  File   : Partition_Inter2d.cxx
//  Author : Benedicte MARTIN
//  Module : GEOM
//  $Header: /cvs/netgen/netgen/libsrc/occ/Partition_Inter2d.cxx,v 1.5 2008/03/31 14:20:28 wabro Exp $

//using namespace std;

#include "Partition_Inter2d.ixx"

#include "utilities.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgo_AsDes.hxx>
#include <BRepLib_MakeVertex.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopOpeBRepDS_Transition.hxx>
#include <TopOpeBRep_EdgesIntersector.hxx>
#include <TopOpeBRep_Point2d.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapIteratorOfMapOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Pnt.hxx>

#ifdef DEB
static Standard_Boolean TestEdges = 0;
static Standard_Integer NbF2d = 0;
static Standard_Integer NbE2d = 0;
#endif

//=======================================================================
//function : getOtherShape
//purpose  :
//=======================================================================

static TopoDS_Shape getOtherShape(const TopoDS_Shape&         theS,
                                  const TopTools_ListOfShape& theSList)
{
  TopTools_ListIteratorOfListOfShape anIt( theSList );
  for ( ; anIt.More(); anIt.Next() )
    if (!theS.IsSame( anIt.Value() ))
      return anIt.Value();

  return TopoDS_Shape();
}

//=======================================================================
//function : findVOnE
//purpose  : on theE, find a vertex close to theV, such that an edge
//           passing through it is an itersection of theF1 and theF2.
//           theE intersects theE2 at theV
//=======================================================================

static Standard_Boolean findVOnE(const TopoDS_Vertex &         theV,
                                 const TopoDS_Edge&            theE,
                                 const TopoDS_Edge&            theE2,
                                 const TopoDS_Shape&           theF1,
                                 const TopoDS_Shape&           theF2,
                                 const Handle(BRepAlgo_AsDes)& theAsDes,
                                 TopoDS_Vertex &               theFoundV)
{
  Standard_Real MinDist2 = ::RealLast();
  gp_Pnt P;

  // check all vertices on theE
  const TopTools_ListOfShape& aVList = theAsDes->Descendant( theE );
  TopTools_ListIteratorOfListOfShape anIt( aVList );
  if (anIt.More())
    P = BRep_Tool::Pnt( theV );
  for ( ; anIt.More(); anIt.Next() )
  {
    // check by distance
    TopoDS_Vertex & V = TopoDS::Vertex( anIt.Value() );
    Standard_Real dist2 = P.SquareDistance( BRep_Tool::Pnt( V ));
    if (dist2 < MinDist2)
      MinDist2 = dist2;
    else
      continue;

    // V is a candidate if among edges passing through V there is one
    // which is an intersection of theF1 and theF2
    TopTools_ListIteratorOfListOfShape anEIt( theAsDes->Ascendant( V ));
    Standard_Boolean isOk = Standard_False;
    for (  ; !isOk && anEIt.More(); anEIt.Next() )
    {
      const TopoDS_Shape & E2 = anEIt.Value();
      if ( theE2.IsSame( E2 ))
        continue;
      const TopTools_ListOfShape & aFList = theAsDes->Ascendant( E2 );
      if (aFList.IsEmpty())
        continue;
      if ( theF1.IsSame( aFList.First() ))
        isOk = theF2.IsSame( aFList.Last() );
      else
        isOk = theF2.IsSame( aFList.First() ) && theF1.IsSame( aFList.Last() );
    }
    if (isOk)
      theFoundV = V;
  }

  if (theFoundV.IsNull())
    return Standard_False;

  // check that MinDist2 is not too large
  Standard_Real f, l;
  TopLoc_Location L;
  Handle(Geom_Curve) aCurve = BRep_Tool::Curve( theE, L, f, l );
  gp_Pnt P1 = aCurve->Value( f );
  gp_Pnt P2 = aCurve->Value( 0.3 * f + 0.7 * l );
  //gp_Pnt P2 = aCurve->Value( 0.5 * ( f + l ));
  if (MinDist2 > P1.SquareDistance( P2 ))
    return Standard_False;

#ifdef DEB
  MESSAGE("findVOnE: found MinDist = " << sqrt (MinDist2));
#endif

  return Standard_True;
}

//=======================================================================
//function : AddVonE
//purpose  : Put V in AsDes as intersection of E1 and E2.
//           Check that vertex equal to V already exists on one
//           of edges, in  such  a  case,  V  is  not added but
//           existing vertex is updated to  be on E1 and E2 and
//           is returned insead of V.
//=======================================================================

TopoDS_Vertex Partition_Inter2d::AddVonE(const TopoDS_Vertex& theV,
                                         const TopoDS_Edge&   E1,
                                         const TopoDS_Edge&   E2,
                                         const Handle(BRepAlgo_AsDes)& AsDes,
                                         const TopoDS_Face&   theF)

{
  //-------------------------------------------------------------
  // test if the points of intersection already exist. If not,
  // add as descendants of the edges.
  // nb: theses points are only vertices of intersection.
  //-------------------------------------------------------------
  const TopTools_ListOfShape& VOnE1 = AsDes->Descendant(E1);
  const TopTools_ListOfShape& VOnE2 = AsDes->Descendant(E2);
  gp_Pnt                      P1,P2;
  TopoDS_Vertex               V1,V2;
  TopTools_ListIteratorOfListOfShape it;
  BRep_Builder                       B;
  TopAbs_Orientation                 O1,O2;
  Standard_Real                      U1,U2;
  Standard_Real                      Tol,Tol1,Tol2;
  Standard_Boolean                   OnE1,OnE2;

  TopoDS_Vertex V    = theV;

  U1 = BRep_Tool::Parameter(V,E1);
  U2 = BRep_Tool::Parameter(V,E2);
  O1 = V.Orientation();
  O2 = O1;
  P1  = BRep_Tool::Pnt(V);
  Tol = BRep_Tool::Tolerance( V );
  OnE1 = OnE2 = Standard_False;

  //-----------------------------------------------------------------
  // Search if the point of intersection is a vertex of E1.
  //-----------------------------------------------------------------
  for (it.Initialize(VOnE1); it.More(); it.Next()) {
    const TopoDS_Vertex& CV = TopoDS::Vertex( it.Value() );
    if (V.IsSame( CV )) {
      V1   = V;
      OnE1 = Standard_True;
      break;
    }
    P2 = BRep_Tool::Pnt( CV );
    Tol1 = 1.1*(Tol + BRep_Tool::Tolerance( CV ));
    if (P1.SquareDistance(P2) <= Tol1*Tol1) {
      V    = CV;
      V1   = V;
      OnE1 = Standard_True;
      break;
    }
  }
  if (OnE1) {
    //-----------------------------------------------------------------
    // Search if the vertex found is still on E2.
    //-----------------------------------------------------------------
    for (it.Initialize(VOnE2); it.More(); it.Next()) {
      if (V.IsSame( it.Value() )) {
        OnE2 = Standard_True;
        V2   = V;
        break;
      }
    }
  }
  if (!OnE2) {
    for (it.Initialize(VOnE2); it.More(); it.Next()) {
      //-----------------------------------------------------------------
      // Search if the point of intersection is a vertex of E2.
      //-----------------------------------------------------------------
      const TopoDS_Vertex& CV = TopoDS::Vertex( it.Value() );
      P2 = BRep_Tool::Pnt( CV );
      Tol2 = 1.1*(Tol + BRep_Tool::Tolerance( CV ));
      if (P1.SquareDistance(P2) <= Tol2*Tol2) {
        V  = CV;
        V2 = V;
        OnE2 = Standard_True;
        break;
      }
    }
  }


  if (!OnE1 && !OnE2 && !theF.IsNull())
  {
    // if 3 faces intersects each others, 3 new edges on them must pass
    // through one vertex but real intersection points of each
    // pair of edges are sometimes more far than a tolerance.
    // Try to analitically find vertices that E1 and E2 must pass trough

    TopoDS_Shape F1 = getOtherShape( theF, AsDes->Ascendant( E1 ));
    TopoDS_Shape F2 = getOtherShape( theF, AsDes->Ascendant( E2 ));
    if (!F1.IsNull() && !F2.IsNull() && !F1.IsSame( F2 ))
    {
      OnE1 = findVOnE ( theV, E1, E2, F1, F2, AsDes, V1 );
      OnE2 = findVOnE ( theV, E2, E1, F1, F2, AsDes, V2 );
      if (OnE2) V = V2;
      if (OnE1) V = V1;
    }
  }

  if (OnE1 && OnE2) {
    if (!V1.IsSame(V2)) {
      // replace V1 with V2 on all edges V1 is on
      Standard_Real UV1;
      TopoDS_Edge   EWE1;
      TopoDS_Vertex VI;
      const TopTools_ListOfShape& EdgeWithV1 = AsDes->Ascendant(V1);

      for (it.Initialize(EdgeWithV1); it.More(); it.Next()) {
        EWE1  = TopoDS::Edge(it.Value());
        VI = V1;
        VI.Orientation(TopAbs_INTERNAL);
        UV1 = BRep_Tool::Parameter(VI,EWE1);
        VI = V2;
        VI.Orientation(TopAbs_INTERNAL);
        B.UpdateVertex( VI, UV1, EWE1, GetTolerance( VI, UV1, EWE1, AsDes));
      }
      AsDes->Replace(V1,V2);
      V = V2;
    }
  }

  // add existing vertices instead of new ones
  if (!OnE1) {
    if (OnE2) {
      V.Orientation(TopAbs_INTERNAL);
      B.UpdateVertex (V, U1, E1, GetTolerance( V, U1, E1, AsDes));
    }
    V.Orientation(O1);
    AsDes->Add(E1,V);
  }
  if (!OnE2) {
    if (OnE1) {
      V.Orientation(TopAbs_INTERNAL);
      B.UpdateVertex (V, U2, E2, GetTolerance( V, U2, E2, AsDes ));
    }
    V.Orientation(O2);
    AsDes->Add(E2,V);
  }

  return V;
}

//=======================================================================
//function : FindEndVertex
//purpose  : Returns a vertex  from  <VertList> having parameter on
//           <E>  closest  to  <f>  or  <l>.  <isFirst>  is True if
//           found vertex is closer  to <f>. <DU> returns parameter
//           difference.
//=======================================================================

TopoDS_Vertex Partition_Inter2d::FindEndVertex(const TopTools_ListOfShape& LV,
                                               const Standard_Real f,
                                               const Standard_Real l,
                                               const TopoDS_Edge&  E,
                                               Standard_Boolean&   isFirst,
                                               Standard_Real&      minDU)
{
  TopoDS_Vertex endV;
  Standard_Real U, endU, min;
  minDU = 1.e10;

  TopTools_ListIteratorOfListOfShape it;
  it.Initialize(LV);
  for (; it.More(); it.Next()) {
    const TopoDS_Vertex& v = TopoDS::Vertex(it.Value());
    U = BRep_Tool::Parameter(v, E);
    min = Min( Abs(U-f), Abs(U-l) );
    if (min < minDU) {
      endV = v;
      endU = U;
      minDU = min;
    }
  }
  if (Abs(endU-f) < Abs(endU-l))
    isFirst = Standard_True;
  else
    isFirst = Standard_False;

  return endV;
}

//=======================================================================
//function : treatClosed
//purpose  : add second vertex to closed edge. Vertex is one of <LV1>
//=======================================================================

static void treatClosed (const TopoDS_Edge& E1,
                          const Standard_Real f,
                          const Standard_Real l,
                          TopTools_ListOfShape& LV1,
                          TopTools_ListOfShape& /*LV2*/)
{
  Standard_Boolean isFirst=0;
  Standard_Real    minDU = 1.e10;
  TopoDS_Vertex endV;
  endV = Partition_Inter2d::FindEndVertex(LV1, f,l, E1, isFirst,minDU);

  if (minDU > Precision::PConfusion())
    return; // not end point

  Standard_Real newU;
  if (isFirst)
    newU = f + (l - f);
  else
    newU = l - (l - f);

  // update end parameter
  BRep_Builder B;
  endV.Orientation(TopAbs_INTERNAL);
  B.UpdateVertex(endV,newU,E1,BRep_Tool::Tolerance(endV));
}

//=======================================================================
//function : EdgesPartition
//purpose  :
//=======================================================================

static void EdgesPartition(const TopoDS_Face&            F,
                           const TopoDS_Edge&            E1,
                           const TopoDS_Edge&            E2,
                           const Handle(BRepAlgo_AsDes)& AsDes,
                           const TopTools_MapOfShape&    NewEdges,
                           const Standard_Boolean        WithOri)
{

  Standard_Real f[3],l[3];
  Standard_Real MilTol2;
  Standard_Real Tol = Max (BRep_Tool::Tolerance(E1),
                           BRep_Tool::Tolerance(E2));
  MilTol2 = Tol * Tol * 10;

  BRep_Tool::Range(E1, f[1], l[1]);
  BRep_Tool::Range(E2, f[2], l[2]);

  BRepAdaptor_Curve CE1(E1,F);
  BRepAdaptor_Curve CE2(E2,F);

  TopoDS_Edge                 EI[3]; EI[1] = E1; EI[2] = E2;
  TopTools_ListOfShape        LV1; // new vertices at intersections on E1
  TopTools_ListOfShape        LV2; // ... on E2
  BRep_Builder                B;

  // if E1 and E2 are results of intersection of F and two connex faces then
  // no need to intersect edges, they can contact by vertices only
  // (encounted an exception in TopOpeBRep_EdgesIntersector in such a case)
  Standard_Boolean intersect = Standard_True;
  TopTools_IndexedMapOfShape ME;
  TopExp::MapShapes(F, TopAbs_EDGE, ME);
  if (!ME.Contains(E1) && ! ME.Contains(E2)) { // if E1 and E2 are new on F
    TopoDS_Shape F1, F2;
    const TopTools_ListOfShape& LF1 = AsDes->Ascendant( E1 );
    F1 = F.IsSame( LF1.First() ) ? LF1.Last() : LF1.First();
    const TopTools_ListOfShape& LF2 = AsDes->Ascendant( E2 );
    F2 = F.IsSame( LF2.First() ) ? LF2.Last() : LF2.First();
    if (!F.IsSame(F2) && !F.IsSame(F1) ) {
      TopExp_Explorer exp(F2, TopAbs_EDGE);
      TopExp::MapShapes(F1, TopAbs_EDGE, ME);
      for (; exp.More(); exp.Next()) {
        if (ME.Contains( exp.Current())) {
          intersect = Standard_False;
          break;
        }
      }
    }
  }

  if (intersect) {
    //------------------------------------------------------
    // compute the points of Intersection in 2D
    //-----------------------------------------------------
    // i.e. fill LV1 and LV2
    TopOpeBRep_EdgesIntersector EInter;
    EInter.SetFaces(F,F);
    Standard_Real TolDub = 1.e-7;
    EInter.ForceTolerances(TolDub,TolDub);
    Standard_Boolean reducesegments = Standard_False;
    EInter.Perform (E1,E2,reducesegments);

    Standard_Boolean rejectreducedsegmentpoints = Standard_False;
    EInter.InitPoint(rejectreducedsegmentpoints);
    for ( ; EInter.MorePoint(); EInter.NextPoint() )
    {
      const TopOpeBRep_Point2d& P2D = EInter.Point();
      const gp_Pnt&    P    = P2D.Value();
      TopoDS_Vertex    V    = BRepLib_MakeVertex(P);

      //-------------------------
      // control the point found.
      //-------------------------
      gp_Pnt P1 = CE1.Value(P2D.Parameter(1));
      gp_Pnt P2 = CE2.Value(P2D.Parameter(2));
      Standard_Real sqd1 = P1.SquareDistance(P);
      Standard_Real sqd2 = P2.SquareDistance(P);
      if (sqd1 > MilTol2 || sqd2 > MilTol2  )
        continue;

      // add a new vertex to the both edges
      Standard_Real toler = Max( Tol, sqrt( Max( sqd1, sqd2 )));
      Standard_Integer i;
      for (i = 1; i <= 2; i++) {
        Standard_Real U = P2D.Parameter(i);
        V.Orientation(TopAbs_INTERNAL);
        B.UpdateVertex( V,U,EI[i], toler);
        TopAbs_Orientation OO = TopAbs_REVERSED;
        if (WithOri) {
          if (P2D.IsVertex(i))
            OO = P2D.Vertex(i).Orientation();
          else if (P2D.Transition(i).Before() == TopAbs_OUT) {
            OO = TopAbs_FORWARD;
          }
          V.Orientation(OO);
          if (i == 1) LV1.Append(V);
          else        LV2.Append(V);
        }
      }
    }
  } // if (intersect)

  //----------------------------------
  // Test the extremities of the edges.
  //----------------------------------
  // add to LV* vertices for vertex-vertex closeness
  Standard_Real U1,U2;
  Standard_Real TolConf2, TolConf;
  TopoDS_Vertex V1[2],V2[2];
  TopExp::Vertices(E1,V1[0],V1[1]);
  TopExp::Vertices(E2,V2[0],V2[1]);

  Standard_Integer i,j,k;
  for (j = 0; j < 2; j++) {
    if (V1[j].IsNull()) continue;
    for ( k = 0; k < 2; k++) {
      if (V2[k].IsNull()) continue;
      gp_Pnt P1 = BRep_Tool::Pnt(V1[j]);
      gp_Pnt P2 = BRep_Tool::Pnt(V2[k]);
      TolConf = BRep_Tool::Tolerance(V1[j]) + BRep_Tool::Tolerance(V2[k]);
      TolConf = Max (Tol, TolConf);
      TolConf2 = TolConf * TolConf;
      if (!intersect)
        TolConf2 *= 100;
      Standard_Real SqDist = P1.SquareDistance(P2);

      if (SqDist <= TolConf2) {
        TopoDS_Vertex V = BRepLib_MakeVertex(P1);
        V.Orientation(TopAbs_INTERNAL);
        U1 = (j == 0) ? f[1] : l[1];
        U2 = (k == 0) ? f[2] : l[2];
        B.UpdateVertex(V,U1,E1,TolConf);
        B.UpdateVertex(V,U2,E2,TolConf);
        LV1.Prepend(V.Oriented(V1[j].Orientation()));
        LV2.Prepend(V.Oriented(V2[k].Orientation()));
      }
    }
  }

  Standard_Boolean AffichPurge = Standard_False;

  if ( LV1.IsEmpty()) return;

  //----------------------------------
  // Purge of all the vertices.
  //----------------------------------
  // remove one of close vertices
  TopTools_ListIteratorOfListOfShape it1LV1,it1LV2,it2LV1;
  gp_Pnt P1,P2;
  Standard_Boolean Purge = Standard_True;

  while (Purge) {
    i = 1;
    Purge = Standard_False;
    for (it1LV1.Initialize(LV1),it1LV2.Initialize(LV2);
         it1LV1.More();
         it1LV1.Next(),it1LV2.Next()) {
      j = 1;
      it2LV1.Initialize(LV1);
      while (j < i) {
        const TopoDS_Vertex& VE1 = TopoDS::Vertex(it1LV1.Value());
        const TopoDS_Vertex& VE2 = TopoDS::Vertex(it2LV1.Value());
        Standard_Real Tol1 = BRep_Tool::Tolerance( VE1 );
        Standard_Real Tol2 = BRep_Tool::Tolerance( VE2 );
        P1 = BRep_Tool::Pnt( VE1 );
        P2 = BRep_Tool::Pnt( VE2 );
        if (P1.IsEqual(P2, Tol1 + Tol2)) {
          LV1.Remove(it1LV1);
          LV2.Remove(it1LV2);
          Purge = Standard_True;
          break;
        }
        j++;
        it2LV1.Next();
      }
      if (Purge) break;
      i++;
    }
  }

  // care of new closed edges, they always intersect with seam at end
  if (V1[0].IsSame( V1[1] ) && NewEdges.Contains(E1) )
    treatClosed (E1, f[1], l[1], LV1, LV2);
  if (V2[0].IsSame( V2[1] ) && NewEdges.Contains(E2) )
    treatClosed (E2, f[2], l[2], LV2, LV1);

  //----------------
  // Stocking vertex
  //----------------

  for ( it1LV1.Initialize( LV1 ); it1LV1.More(); it1LV1.Next())
    Partition_Inter2d::AddVonE (TopoDS::Vertex( it1LV1.Value()),
                                E1, E2, AsDes, F);
}

//=======================================================================
//function : CompletPart2d
//purpose  : Computes the intersections between the edges stored
//           is AsDes as descendants of <F> . Intersections is computed
//           between two edges if one of them is bound in NewEdges.
//=======================================================================

void Partition_Inter2d::CompletPart2d (const Handle(BRepAlgo_AsDes)&   AsDes,
                                       const TopoDS_Face&              F,
                                       const TopTools_MapOfShape&      NewEdges)
{

#ifdef DEB
  NbF2d++;
  NbE2d = 0;
#endif

  //Do not intersect the edges of a face
  TopTools_IndexedMapOfShape EdgesOfFace;
  TopExp::MapShapes( F, TopAbs_EDGE , EdgesOfFace);

  //-------------------------------------------------------------------
  // compute the intersection2D on the faces touched by the intersection3D
  //-------------------------------------------------------------------
  TopTools_ListIteratorOfListOfShape it1LE ;
  TopTools_ListIteratorOfListOfShape it2LE ;

  //-----------------------------------------------
  // Intersection edge-edge.
  //-----------------------------------------------
  const TopTools_ListOfShape&        LE = AsDes->Descendant(F);
  TopoDS_Vertex                      V1,V2;
  Standard_Integer                   j, i = 1;

  TopoDS_Face FF = F;
  FF.Orientation(TopAbs_FORWARD);

  for ( it1LE.Initialize(LE) ; it1LE.More(); it1LE.Next()) {
    const TopoDS_Edge& E1 = TopoDS::Edge(it1LE.Value());
    j = 1;
    it2LE.Initialize(LE);

    while (j < i && it2LE.More()) {
      const TopoDS_Edge& E2 = TopoDS::Edge(it2LE.Value());
      //----------------------------------------------------------
      // Intersections of the new edges obtained by intersection
      // between them and with the restrictions edges
      //----------------------------------------------------------
      if ( (!EdgesOfFace.Contains(E1) || !EdgesOfFace.Contains(E2)) &&
           (NewEdges.Contains(E1) || NewEdges.Contains(E2)) ) {
        EdgesPartition(FF,E1,E2,AsDes,NewEdges,Standard_True);
      }
      it2LE.Next();
      j++;
    }
    i++;
  }
}

//=======================================================================
//function : GetTolerance
//purpose  : Returns  tolerance  theV   must   have  atfer  its
//           addition to theE with  theU parameter. theAsDes is
//           used to find pcurves of theE
//=======================================================================

Standard_Real Partition_Inter2d::GetTolerance
                         (const TopoDS_Vertex &         theV,
                          const Standard_Real           theU,
                          const TopoDS_Edge &           theE,
                          const Handle(BRepAlgo_AsDes)& theAsDes)
{
  Standard_Real aTol = BRep_Tool::Tolerance( theV );
  gp_Pnt aPnt = BRep_Tool::Pnt( theV );

  // check point on 3D curve
  Standard_Real f,l;
  Handle(Geom_Curve) C = BRep_Tool::Curve( theE, f, l );
  if (!C.IsNull())
    aTol = Max ( aTol, aPnt.Distance( C->Value( theU )));

  // check points on pcurves
  const TopTools_ListOfShape& aFList = theAsDes->Ascendant( theE );
  TopTools_ListIteratorOfListOfShape aFIt( aFList );
  for (  ; aFIt.More(); aFIt.Next() )
  {
    const TopoDS_Face& F = TopoDS::Face( aFIt.Value() );
    Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface( theE, F, f, l );
    if (!pcurve.IsNull())
    {
      gp_Pnt2d aPnt2d = pcurve->Value( theU );
      TopLoc_Location L;
      Handle(Geom_Surface) S = BRep_Tool::Surface( F, L );
      gp_Pnt aPntOnS = S->Value( aPnt2d.X(), aPnt2d.Y() );
      if (!L.IsIdentity())
        aPntOnS.Transform( L.Transformation() );
      aTol = Max ( aTol, aPnt.Distance( aPntOnS ));
    }
  }

  return aTol;
}

#endif
