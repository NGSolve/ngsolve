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
//  File   : Partition_Inter3d.cxx
//  Author : Benedicte MARTIN
//  Module : GEOM
//  $Header: /cvs/netgen/netgen/libsrc/occ/Partition_Inter3d.cxx,v 1.6 2008/03/31 14:20:28 wabro Exp $

//using namespace std;
#include <climits>

#include "Partition_Inter2d.hxx"
#include "Partition_Inter3d.ixx"
#include "utilities.h"

#include <BRepAlgo_AsDes.hxx>
#include <BRepAlgo_Image.hxx>
#include <BRepLib.hxx>
#include <BRepOffset_Tool.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>

#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>

#include <TopOpeBRepTool_BoxSort.hxx>
#include <TopTools_DataMapIteratorOfDataMapOfShapeListOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>

#ifdef DEB
#include <DBRep.hxx>
#endif

#include <BRepLib_MakeVertex.hxx>
#include <BRepTools.hxx>
#include <Extrema_ExtPS.hxx>
#include <Extrema_POnSurf.hxx>
#include <Geom2dAPI_ProjectPointOnCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Precision.hxx>
#include <TColStd_MapOfInteger.hxx>
#include <TopOpeBRepBuild_Builder.hxx>
#include <TopOpeBRepDS_BuildTool.hxx>
#include <TopOpeBRepDS_CurveExplorer.hxx>
#include <TopOpeBRepDS_HDataStructure.hxx>
#include <TopOpeBRepDS_Interference.hxx>
#include <TopOpeBRepDS_PointIterator.hxx>
#include <TopOpeBRepDS_Transition.hxx>
#include <TopOpeBRepTool_CurveTool.hxx>
#include <TopOpeBRepTool_GeomTool.hxx>
#include <TopOpeBRepTool_OutCurveType.hxx>
#include <TopOpeBRep_DSFiller.hxx>
#include <TopTools_DataMapIteratorOfDataMapOfShapeShape.hxx>
#include <stdio.h>

//=======================================================================
//function : Partition_Inter3d
//purpose  : 
//=======================================================================

Partition_Inter3d::Partition_Inter3d()
{
}
//=======================================================================
//function : Partition_Inter3d
//purpose  : 
//=======================================================================

Partition_Inter3d::Partition_Inter3d(const Handle(BRepAlgo_AsDes)& AsDes)
  :myAsDes(AsDes)
{
  mySectionEdgesAD = new BRepAlgo_AsDes;
}

//=======================================================================
//function : CompletPart3d
//purpose  : FaceShapeMap is just to know the shape a face belongs to
//=======================================================================

void Partition_Inter3d::CompletPart3d(const TopTools_ListOfShape& SetOfFaces1,
				      const TopTools_DataMapOfShapeShape& FaceShapeMap)
{
  if (myAsDes.IsNull())
    myAsDes = new BRepAlgo_AsDes;
  
  TopTools_ListIteratorOfListOfShape it;

  //---------------------------------------------------------------
  // Construction of bounding boxes.
  //---------------------------------------------------------------

  BRep_Builder B;
  TopoDS_Compound CompOS;
  B.MakeCompound(CompOS);
  for (it.Initialize(SetOfFaces1); it.More(); it.Next())
    B.Add(CompOS, it.Value());
    
  TopOpeBRepTool_BoxSort BOS;
  BOS.AddBoxesMakeCOB(CompOS,TopAbs_FACE);

  for (it.Initialize(SetOfFaces1); it.More(); it.Next()) {
    TopoDS_Face F1 = TopoDS::Face(it.Value());
    
    // avoid intersecting faces of one shape
    TopoDS_Shape S1;
    if (FaceShapeMap.IsBound(F1)) S1 = FaceShapeMap.Find(F1);

    // to filter faces sharing an edge
    TopTools_IndexedMapOfShape EM;
    TopExp::MapShapes( F1, TopAbs_EDGE, EM);
    
    TColStd_ListIteratorOfListOfInteger itLI = BOS.Compare(F1);
    for (; itLI.More(); itLI.Next()) {
      TopoDS_Face F2 = TopoDS::Face(BOS.TouchedShape(itLI));
      if (F1.IsSame(F2) || IsDone(F1,F2))
	continue;

      TopoDS_Shape S2;
      if (FaceShapeMap.IsBound(F2)) S2 = FaceShapeMap.Find(F2);
      if (!S1.IsNull() && S1.IsSame(S2))
	continue; // descendants of one shape

      TopExp_Explorer expE (F2, TopAbs_EDGE);
      for ( ; expE.More(); expE.Next())
	if (EM.Contains( expE.Current() ))
	  break;
      if (expE.More())
      {
        // faces have a common edge, check if they are a tool and a face
        // generated by the tool in another shape; in that case they are
        // to be intersected
        TopLoc_Location L1, L2;
        Handle(Geom_Surface) S1 = BRep_Tool::Surface( F1, L1 );
        Handle(Geom_Surface) S2 = BRep_Tool::Surface( F2, L2 );
        if ( S1 != S2 || L1 != L2 )
          continue;
      }

      F1.Orientation(TopAbs_FORWARD);
      F2.Orientation(TopAbs_FORWARD);
      FacesPartition(F1,F2);	  
    }

    // mark as modified a face which has at least one new edge
    if (!myAsDes->HasDescendant( F1 ))
      continue;
    TopTools_ListIteratorOfListOfShape itE (myAsDes->Descendant( F1 ));
    for ( ; itE.More(); itE.Next()) {
      if (myNewEdges.Contains( itE.Value())) {
	myTouched.Add( F1 );
	break;
      }
    }
  }
}

//=======================================================================
//function : PutInBounds
//purpose  : 
//=======================================================================

static void PutInBounds (const TopoDS_Face&          F,
			 const TopoDS_Edge&          E,
			 Handle(Geom2d_Curve)&       C2d)
{
  Standard_Real   umin,umax,vmin,vmax;
  Standard_Real   f,l;
  BRep_Tool::Range(E,f,l);

  TopLoc_Location L; // Recup S avec la location pour eviter la copie.
  Handle (Geom_Surface) S   = BRep_Tool::Surface(F,L);

  if (S->IsKind(STANDARD_TYPE(Geom_RectangularTrimmedSurface))) {
    S = (*(Handle_Geom_RectangularTrimmedSurface*)&S)->BasisSurface();
  }
  if (!S->IsUPeriodic() && !S->IsVPeriodic())
    return;

  BRepTools::UVBounds(F,umin,umax,vmin,vmax);

  gp_Pnt2d Pf = C2d->Value(f);
  gp_Pnt2d Pl = C2d->Value(l);
  const Standard_Real Um = 0.34*f + 0.66*l;
  gp_Pnt2d Pm = C2d->Value( Um );

  // sometimes on shpere, pcurve is out of domain by V though S is
  // UPeriodic, sometimes it is in domain but nontheless it has
  // wrong position.
  // Check pcurve position by 3D point
  if (S->IsKind(STANDARD_TYPE( Geom_SphericalSurface )))
  {
    // get point on the surface
    gp_Pnt Ps = S->Value( Pm.X(), Pm.Y() );
    // get point on the edge
    Handle(Geom_Curve) C = BRep_Tool::Curve( E, f, l );
    gp_Pnt Pc = C->Value( Um );
    // compare points
    Standard_Real TolE = BRep_Tool::Tolerance( E );
    if ( Pc.SquareDistance( Ps ) * 0.95 < TolE * TolE )
      return; // OK

    // find good UV for Pc: project Pc on S
    GeomAdaptor_Surface  SA (S);
    Extrema_ExtPS anExtPS (Pc, SA,
                           SA.UResolution( TolE ), SA.VResolution( TolE ));
    if (anExtPS.IsDone())
    {
      Standard_Integer i, nbExt = anExtPS.NbExt();
      Extrema_POnSurf aPOnSurf;
      for (i = 1; i <= nbExt; ++i )
	if (anExtPS.Value( i ) <= TolE)               // V6.3
	  // if (anExtPS.SquareDistance( i ) <= TolE)   // V6.5
	  {
          aPOnSurf = anExtPS.Point( i );
          break;
        }
      if (i <= nbExt) {
        // a point found
        Standard_Real u, v;
        aPOnSurf.Parameter( u, v );
        gp_Pnt2d aGoodPm ( u, v );
        C2d->Translate( Pm , aGoodPm );
      }
    }
  }

  //---------------
  // Recadre en U.
  //---------------
  if (S->IsUPeriodic()) {
    Standard_Real period  = S->UPeriod();
    Standard_Real eps     = period*1.e-6;
    Standard_Real minC    = Min(Pf.X(),Pl.X()); minC = Min(minC,Pm.X());
    Standard_Real maxC    = Max(Pf.X(),Pl.X()); maxC = Max(maxC,Pm.X());
    Standard_Real du = 0.;
    if (minC< umin - eps) {
      du = (int((umin - minC)/period) + 1)*period;
    }
    if (minC > umax + eps) {
      du = -(int((minC - umax)/period) + 1)*period;
    }
    if (du != 0) {
      gp_Vec2d T1(du,0.);
      C2d->Translate(T1);
      minC += du; maxC += du;
    }
    // Ajuste au mieux la courbe dans le domaine.
    if (maxC > umax +100*eps) {
      Standard_Real d1 = maxC - umax;
      Standard_Real d2 = umin - minC + period;
      if (d2 < d1) du =-period;
      if ( du != 0.) {
	gp_Vec2d T2(du,0.);
	C2d->Translate(T2);
      }
    }
  }
  //------------------
  // Recadre en V.
  //------------------
  if (S->IsVPeriodic()) {
    Standard_Real period  = S->VPeriod();
    Standard_Real eps     = period*1.e-6;
    Standard_Real minC    = Min(Pf.Y(),Pl.Y()); minC = Min(minC,Pm.Y());
    Standard_Real maxC    = Max(Pf.Y(),Pl.Y()); maxC = Max(maxC,Pm.Y());
    Standard_Real dv = 0.;
    if (minC< vmin - eps) {
      dv = (int((vmin - minC)/period) + 1)*period;
    }
    if (minC > vmax + eps) {
      dv = -(int((minC - vmax)/period) + 1)*period;
    }
    if (dv != 0) {
      gp_Vec2d T1(0.,dv);
      C2d->Translate(T1);
      minC += dv; maxC += dv;
    }
    // Ajuste au mieux la courbe dans le domaine.
    if (maxC > vmax +100*eps) {
      Standard_Real d1 = maxC - vmax;
      Standard_Real d2 = vmin - minC + period;
      if (d2 < d1) dv =-period;
      if ( dv != 0.) {
	gp_Vec2d T2(0.,dv);
	C2d->Translate(T2);
      }
    }
  }
}

//=======================================================================
//function : Inter3D
//purpose  : 
//=======================================================================

void Partition_Inter3d::Inter3D(const TopoDS_Face& F1,
				const TopoDS_Face& F2,
				TopTools_ListOfShape& L)
{
  BRep_Builder B;
  
  // fill the data Structure
  Handle(TopOpeBRepDS_HDataStructure) DatStr = new TopOpeBRepDS_HDataStructure();
  TopOpeBRep_DSFiller DSFiller;
  DSFiller.Insert(F1,F2,DatStr);

  // define the GeomTool used by the DSFiller :
  // compute BSpline of degree 1 on intersection curves.
  Standard_Real tol3dAPPROX = 1e-7;
  Standard_Real tol2dAPPROX = 1e-7;
  TopOpeBRepTool_GeomTool GT2 (TopOpeBRepTool_APPROX);  
  GT2.SetTolerances(tol3dAPPROX,tol2dAPPROX);
  TopOpeBRepDS_BuildTool  BT(GT2);

  // Perform Section
  TopOpeBRepBuild_Builder TopB(BT);
  TopB.Perform(DatStr);

  // ===============
  // Store new edges
  // ===============
  
  L.Clear();
  TopOpeBRepDS_CurveExplorer cex(DatStr->DS());
  for (; cex.More(); cex.Next()) {
    const TopOpeBRepDS_Curve& CDS = cex.Curve();
    Standard_Integer ic = cex.Index();
    Handle(Geom2d_Curve) pc1 = CDS.Curve1();
    Handle(Geom2d_Curve) pc2 = CDS.Curve2();
    
    TopTools_ListIteratorOfListOfShape itLE = TopB.NewEdges(ic);
    while (itLE.More()) {
      TopoDS_Edge E = TopoDS::Edge(itLE.Value());
      
      PutInBounds (F1,E,pc1);
      PutInBounds (F2,E,pc2);
      
      B.UpdateEdge (E,pc1,F1,0.);
      B.UpdateEdge (E,pc2,F2,0.);
      
      L.Append (E);
      
      itLE.Next();
      if (itLE.More()) {
	pc1 = Handle(Geom2d_Curve)::DownCast(pc1->Copy());
	pc2 = Handle(Geom2d_Curve)::DownCast(pc2->Copy());
      }
    }
  }

  // ========================
  // store same domain faces 
  // ========================


  if ( DatStr->HasSameDomain( F1 ))
  {
    TopTools_ListOfShape emptyList;
    if (!mySameDomainFM.IsBound(F1))
      mySameDomainFM.Bind(F1,emptyList);
    if (!mySameDomainFM.IsBound(F2))
      mySameDomainFM.Bind(F2,emptyList);
    mySameDomainFM(F1).Append(F2);
    mySameDomainFM(F2).Append(F1);
  }

  // ====================
  // Store section edges
  // ====================

  const TopOpeBRepDS_DataStructure& DS = DatStr->DS();
  Standard_Integer j,i,nse = DS.NbSectionEdges();
  if (nse == 0) return;

    
  TopoDS_Vertex V, sdeV1, sdeV2;
  TopTools_MapOfShape MV;
  TopTools_ListOfShape LSE; // list of section edges
  TopoDS_Face dummyF;
  
  for (i = 1; i <= nse; i++)
  {
    const TopoDS_Edge & se = DS.SectionEdge(i);
    if (! TopB.IsSplit(se,TopAbs_ON))
      continue;
    LSE.Append( se );

    // add vertices where section edges interferes with other
    // edges as its descendant in myAsDes
    
    TopoDS_Edge sde, oe; // same domain, other edge
    if (DatStr->HasSameDomain(se)) {
      sde = TopoDS::Edge( DatStr->SameDomain(se).Value() );
      TopExp::Vertices( sde, sdeV1, sdeV2);
    }
    TColStd_MapOfInteger MIV; // indices of added edges
    TopOpeBRepDS_PointIterator itP (DS.ShapeInterferences( se ));
    itP.SupportKind( TopOpeBRepDS_EDGE );
    // loop on intersections of se
    for (; itP.More(); itP.Next()) {
      oe = TopoDS::Edge( DS.Shape( itP.Support()));
      if (itP.IsVertex()) {
        // there is a vertex at intersection
	if ( !MIV.Add( itP.Current() ))
	  continue;
	V = TopoDS::Vertex( DS.Shape( itP.Current()));
	if ( !sde.IsNull() && (V.IsSame(sdeV1) || V.IsSame(sdeV2)) )
	  oe = sde;
	V = ReplaceSameDomainV( V , oe );
	V.Orientation( TopAbs_INTERNAL);
	B.UpdateVertex( V, itP.Parameter(), se, 0.); // AddVonE() sets real U
      }
      else {
        // create a new vertex at the intersection point
	const TopOpeBRepDS_Point& DSP = DS.Point( itP.Current());
	V = BRepLib_MakeVertex( DSP.Point() );
	V.Orientation( TopAbs_INTERNAL);
	B.UpdateVertex( V, itP.Parameter(), se, DSP.Tolerance());
	// make V be on the other edge
	TopOpeBRepDS_PointIterator itOP (DS.ShapeInterferences( oe ));
	for (; itOP.More(); itOP.Next()) {
	  const TopOpeBRepDS_Point& ODSP = DS.Point( itOP.Current());
	  if ( DSP.IsEqual (ODSP)) {
	    B.UpdateVertex( V, itOP.Parameter(), TopoDS::Edge(oe), ODSP.Tolerance());
	    break;
	  }
	}
      }
      // add V on the both intersecting edges
      TopoDS_Vertex addedV = Partition_Inter2d::AddVonE( V,se,oe,myAsDes,dummyF);
      if (!addedV.IsSame( V ))
	mySameDomainVM.Bind (V, addedV); // equal vertex is already there

      MV.Add( addedV ); // to ease storage of vertices of ON splits
    }
  }

  // add section edge to the face it intersects and find
  // splits ON that do not have same domain pair
  
  TopB.SplitSectionEdges(); // let TopB find ON splits

  TopTools_MapOfShape SPM; // map of ON splits
  TopTools_IndexedMapOfShape ME[2];
  TopExp::MapShapes( F1, TopAbs_EDGE, ME[1]);
  TopExp::MapShapes( F2, TopAbs_EDGE, ME[0]);

  TopTools_ListIteratorOfListOfShape itSP, itLSE (LSE);
  while ( itLSE.More() ) {

    TopoDS_Edge se = TopoDS::Edge( itLSE.Value() );

    // move itLSE to the next se
    Standard_Integer ancRank = DS.AncestorRank(se);
    if (ME[ancRank-1].Contains( se ))
    {
      LSE.Remove( itLSE ); // se is an edge of face it intersects
      continue;
    }
    else
    {
      itLSE.Next();
    }

    const TopoDS_Face& F = (ancRank == 1) ? F2 : F1;

    // add se to face but dont add twice
    TopTools_ListIteratorOfListOfShape itE( myAsDes->Descendant( F ));
    if (myAsDes->HasDescendant( F )) {
      for ( ; itE.More(); itE.Next())
	if (se.IsSame( itE.Value() ))
	  break;
    }
    if (!itE.More())
    {
      myAsDes->Add( F, se );

      // check se pcurve on F
      Standard_Real tol, f,l, umin=1e100, umax=-1e100;
      Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface( se, F, f,l);
      if (pc.IsNull()) {
	itSP.Initialize( TopB.Splits(se,TopAbs_ON) );
	for ( ; itSP.More(); itSP.Next()) {
	  const TopoDS_Edge& E = TopoDS::Edge ( itSP.Value());
	  BRep_Tool::Range(E, f, l);
	  umin = Min( umin, f);
	  umax = Max( umax, l);
	}
	Handle(Geom_Curve) C3d = BRep_Tool::Curve( se, f, l);
	if (umin < umax) // sometimes umin == umax for closed edge
	  C3d = new Geom_TrimmedCurve( C3d, umin, umax);
	pc = TopOpeBRepTool_CurveTool::MakePCurveOnFace (F,C3d,tol);
	if (pc.IsNull()) {
	  MESSAGE (" CANT BUILD PCURVE ");
	}
	B.UpdateEdge( se, pc, F, tol);
      }
    }

    // to detect splits that do not have same domain pair
    // ie which split a face into parts and not pass by its boundary
    itSP.Initialize( TopB.Splits(se,TopAbs_ON) );
    for ( ; itSP.More(); itSP.Next()) {
      const TopoDS_Shape& SP = itSP.Value();
      if (!SPM.Add( SP ))
	SPM.Remove( SP );
    }
  }

  // store vertices of ON splits and bind section edges to faces
  
  for (itLSE.Initialize (LSE); itLSE.More(); itLSE.Next())
  {
    const TopoDS_Shape& se = itLSE.Value();

    Standard_Integer ancRank = DS.AncestorRank(se);
    TopoDS_Face F = (ancRank == 1) ? F2 : F1;

    // add vertices of ON splits which have no same domain pair
    Standard_Boolean added = Standard_False;
    itSP.Initialize( TopB.Splits(se,TopAbs_ON) );
    for ( ; itSP.More(); itSP.Next())
    {
      if (!SPM.Contains( itSP.Value() ))
	continue;
      
      const TopoDS_Edge& S = TopoDS::Edge ( itSP.Value());

      added = Standard_True;
      mySectionEdgesAD->Add( F, se );
      
      TopoDS_Vertex VS[2];
      TopExp::Vertices (S, VS[0], VS[1]);
      for (j=0; j<2; ++j)
      {
	if (mySameDomainVM.IsBound( VS[j] ))
	  VS[j] = TopoDS::Vertex( mySameDomainVM( VS[j] ));
	if ( !MV.Contains( VS[j] )) {
	  // find equal vertex on se - point interference
	  gp_Pnt P1 = BRep_Tool::Pnt( VS[j] );
	  TopTools_ListIteratorOfListOfShape itV( myAsDes->Descendant(se) );
	  for (; itV.More(); itV.Next()) {
	    V = TopoDS::Vertex( itV.Value() );
            if ( V.IsSame( VS[j] ))
              break;
	    gp_Pnt P2 = BRep_Tool::Pnt( V );
	    if (P1.IsEqual( P2, Precision::Confusion())) {
	      mySameDomainVM.Bind (VS[j], V);
	      VS[j] = V;
	      break;
	    }
	  }
	  if (!itV.More())  // no interferences with edges
	    myAsDes->Add( se, VS[j]);
	}

        // add ends of ON splits to F in order to detect later
        // if a split is on face in IsSplitOn()
	mySectionEdgesAD->Add( F, VS[j]);
      }
      // in the descendants of F, first go ends of an ON split and
      // then a split itself
      mySectionEdgesAD->Add( F, S );
    }
    if (!added)
      mySectionEdgesAD->Add( F, se );
    
    myNewEdges.Add( se );
  }
}

//=======================================================================
//function : FacesPartition
//purpose  : 
//=======================================================================

void Partition_Inter3d::FacesPartition(const TopoDS_Face& F1,
				       const TopoDS_Face& F2)
     //(const TopTools_DataMapOfShapeListOfShape& /*SetOfFaces2*/)
{
  TopTools_ListOfShape LInt;

  Inter3D (F1,F2,LInt);
  
  StorePart3d (F1,F2,LInt);
}

//=======================================================================
//function : SetDone
//purpose  : 
//=======================================================================

void Partition_Inter3d::SetDone(const TopoDS_Face& F1, 
				const TopoDS_Face& F2)
{
  if (!myDone.IsBound(F1)) {
    TopTools_ListOfShape emptyList;
    myDone.Bind(F1,emptyList);
  }
  myDone(F1).Append(F2);
  if (!myDone.IsBound(F2)) {
    TopTools_ListOfShape emptyList;
    myDone.Bind(F2,emptyList);
  }
  myDone(F2).Append(F1);
}

//=======================================================================
//function : IsDone
//purpose  : 
//=======================================================================

Standard_Boolean Partition_Inter3d::IsDone(const TopoDS_Face& F1, 
					   const TopoDS_Face& F2) 

  const 
{
  if (myDone.IsBound(F1)) {
    TopTools_ListIteratorOfListOfShape it (myDone(F1));
    for (; it.More(); it.Next()) {
      if (it.Value().IsSame(F2)) return Standard_True;
    }
  }
  return Standard_False;
}

//=======================================================================
//function : StorePart3d
//purpose  : 
//=======================================================================

void Partition_Inter3d::StorePart3d(const TopoDS_Face& F1, 
				    const TopoDS_Face& F2, 
				    const TopTools_ListOfShape& LInt)
{
  if (!LInt.IsEmpty()) {
    myAsDes->Add( F1,LInt);
    myAsDes->Add( F2,LInt);

    TopTools_ListIteratorOfListOfShape it(LInt);
    for (; it.More(); it.Next()) {

      TopoDS_Edge E = TopoDS::Edge(it.Value());

      BRep_Builder B;
      B.SameParameter(E,Standard_False);
      BRepLib::SameParameter(E,1.0e-7);
      
      myNewEdges.Add(E);
    }
  }
  SetDone(F1,F2);
}

//=======================================================================
//function : TouchedFaces
//purpose  : 
//=======================================================================

TopTools_MapOfShape& Partition_Inter3d::TouchedFaces()
{
  return myTouched;
}

//=======================================================================
//function : AsDes
//purpose  : 
//=======================================================================

Handle(BRepAlgo_AsDes) Partition_Inter3d::AsDes() const 
{
  return myAsDes;
}

//=======================================================================
//function : NewEdges
//purpose  : 
//=======================================================================

TopTools_MapOfShape& Partition_Inter3d::NewEdges() 
{
  return myNewEdges;
}

//=======================================================================
//function : Affiche
//purpose  : 
//=======================================================================

void Partition_Inter3d::Affiche(const TopTools_ListOfShape& SetOfFaces) const
{
#ifdef DEB
  char PSection[1024];
  char *section=PSection;
  Standard_Integer i = 0;
  Standard_Real j=1;
  TopTools_ListOfShape aList;
  TopTools_ListIteratorOfListOfShape it;
  for (it.Initialize(SetOfFaces); it.More(); it.Next()) {
    const TopoDS_Shape& OS = it.Value();
    aList=myAsDes->Descendant(OS);
    MESSAGE ( " the number of items stored in the list " << j << " :  " << aList.Extent() )
    j++;
    TopTools_ListIteratorOfListOfShape itaList;
    for (itaList.Initialize(aList); itaList.More(); itaList.Next()) {
      const TopoDS_Shape& SS = itaList.Value();
      i++;
      sprintf(PSection,"section_%d",i);
      DBRep::Set(section,SS);  
    }
  }
#endif
}

//=======================================================================
//function : SameDomain
//purpose  : 
//=======================================================================

const TopTools_ListOfShape& Partition_Inter3d::SameDomain(const TopoDS_Face& F) const
{
  if (mySameDomainFM.IsBound( F ))
    return mySameDomainFM (F);

  static TopTools_ListOfShape emptyList;
  return emptyList;
}

//=======================================================================
//function : HasSameDomainF
//purpose  : Return true if F has same domain faces
//=======================================================================

Standard_Boolean Partition_Inter3d::HasSameDomainF(const TopoDS_Shape& F) const
{
  return mySameDomainFM.IsBound( F );
}

//=======================================================================
//function : IsSameDomain
//purpose  : Return true if F1 and F2 are same domain faces
//=======================================================================

Standard_Boolean Partition_Inter3d::IsSameDomainF(const TopoDS_Shape& F1,
						 const TopoDS_Shape& F2) const
{
  if (mySameDomainFM.IsBound( F1 )) {
    TopTools_ListIteratorOfListOfShape it (mySameDomainFM( F1 ));
    for (; it.More(); it.Next()) 
      if (F2.IsSame( it.Value()))
	return Standard_True;
  }
  return F1.IsSame( F2 );
}

//=======================================================================
//function : ReplaceSameDomainV
//purpose  : return same domain vertex of  V if it was replaced
//           and make this vertex to be on E too, else return V
//=======================================================================

TopoDS_Vertex Partition_Inter3d::ReplaceSameDomainV(const TopoDS_Vertex& V,
						    const TopoDS_Edge&   E) const
{
  TopoDS_Vertex SDV = V;
  if (mySameDomainVM.IsBound( V )) {

    TopoDS_Vertex V1,V2;
    TopExp::Vertices(E,V1,V2);
    Standard_Boolean isClosed = V1.IsSame( V2 ) && V.IsSame(V1);

    SDV = TopoDS::Vertex( mySameDomainVM(V) );
    Standard_Real tol = BRep_Tool::Tolerance( V );
    BRep_Builder B;
    SDV.Orientation( V.Orientation());

    if (isClosed) {
      Standard_Real f, l;
      BRep_Tool::Range (E, f, l);
      Standard_Boolean isFirst = IsEqual( BRep_Tool::Parameter(V,E), f );
      B.UpdateVertex(SDV, (isFirst ? f : l), E, tol);
      SDV.Reverse();
      B.UpdateVertex(SDV, (isFirst ? l : f), E, tol);
    }
    else
      B.UpdateVertex (SDV, BRep_Tool::Parameter(V,E), E, tol);
      
  }
  return SDV;
}

//=======================================================================
//function : SectionEdgesAD
//purpose  : 
//=======================================================================

Handle(BRepAlgo_AsDes) Partition_Inter3d::SectionEdgesAD() const
{
  return mySectionEdgesAD;
}

//=======================================================================
//function : IsSectionEdge
//purpose  : return True if  E  is  an  edge  of  a face and it
//           intersects an other face
//=======================================================================

Standard_Boolean
  Partition_Inter3d::IsSectionEdge(const TopoDS_Edge& E) const
{
  return mySectionEdgesAD->HasAscendant(E);
}

//=======================================================================
//function : HasSectionEdge
//purpose  : return True if an  edge  of  F intersects an other
//           face or F is intersected by edge of an other face
//=======================================================================

Standard_Boolean
  Partition_Inter3d::HasSectionEdge(const TopoDS_Face& F) const
{
  return mySectionEdgesAD->HasDescendant(F);
}

//=======================================================================
//function : IsSplitOn
//purpose  : return True if NewE is split of OldE on F
//=======================================================================

Standard_Boolean
  Partition_Inter3d::IsSplitOn(const TopoDS_Edge& NewE,
			       const TopoDS_Edge& OldE,
			       const TopoDS_Face& F) const
{
  if (! mySectionEdgesAD->HasDescendant(F))
    return Standard_False;

  TopTools_ListIteratorOfListOfShape itE ( mySectionEdgesAD->Descendant(F) );
  for ( ; itE.More(); itE.Next()) {
    if ( itE.Value().ShapeType() != TopAbs_EDGE ||
	! OldE.IsSame ( itE.Value() ))
      continue;
    // an edge encountered, its vertices and a split come next
    itE.Next();
    if (!itE.More()) break;
    const TopoDS_Shape& V3 = itE.Value();
    if (V3.ShapeType() != TopAbs_VERTEX) continue;
    itE.Next();
    if (!itE.More()) break;
    const TopoDS_Shape& V4 = itE.Value();
    if (V4.ShapeType() != TopAbs_VERTEX) continue;

    TopoDS_Vertex V1, V2;
    TopExp::Vertices( OldE, V1, V2);
    
    if ( V1.IsSame(V2) &&
	(V1.IsSame(V3) || V1.IsSame(V4)) ) {
      // closed old edge; use the split for the test 
      itE.Next();
      if (!itE.More()) break;
      const TopoDS_Edge& split = TopoDS::Edge( itE.Value() );
      // check distance at middle point of NewE
      Standard_Real f1,l1, f2,l2;
      Handle(Geom2d_Curve) PC1 = BRep_Tool::CurveOnSurface( split, F ,f1,l1);
      if (!PC1.IsNull()) {
	Handle(Geom2d_Curve) PC2 = BRep_Tool::CurveOnSurface(NewE, F ,f2,l2);
	gp_Pnt2d P = PC2->Value( 0.5*(f2+l2) );
	Geom2dAPI_ProjectPointOnCurve proj (P, PC1, f1, l1);
	if (proj.NbPoints() &&
	    proj.LowerDistance() <= Precision::Confusion())
	  return Standard_True;
      }
      else {
        Handle(Geom_Curve) C1 = BRep_Tool::Curve( split ,f1,l1);
	Handle(Geom_Curve) C2 = BRep_Tool::Curve( NewE  ,f2,l2);
	gp_Pnt P = C2->Value( 0.5*(f2+l2) );
	GeomAPI_ProjectPointOnCurve proj (P, C1, f1, l1);
	if (proj.NbPoints() &&
	    proj.LowerDistance() <= Precision::Confusion())
	  return Standard_True;
      }
    }
    else {
      Standard_Real u3 = BRep_Tool::Parameter( TopoDS::Vertex(V3), OldE);
      Standard_Real u4 = BRep_Tool::Parameter( TopoDS::Vertex(V4), OldE);

      Standard_Real f,l, u;
      BRep_Tool::Range( NewE, f,l);
      u = 0.5*(f+l);
      f = Min(u3,u4);
      l = Max(u3,u4);

      if (u <= l && u >= f)
        return Standard_True;
    }
  }
  return Standard_False;
}

//=======================================================================
//function : SectionEdgeFaces
//purpose  : return faces cut by section edge
//=======================================================================

const TopTools_ListOfShape&
  Partition_Inter3d::SectionEdgeFaces(const TopoDS_Edge& SecE) const
{
  return mySectionEdgesAD->Ascendant( SecE );
}

#endif
