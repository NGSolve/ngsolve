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
//  File   : Partition_Spliter.cxx
//  Author : Benedicte MARTIN
//  Module : GEOM
//  $Header: /cvs/netgen/netgen/libsrc/occ/Partition_Spliter.cxx,v 1.7 2008/03/31 14:20:28 wabro Exp $

//using namespace std;
#include <climits>
#include "Partition_Inter2d.hxx"
#include "Partition_Inter3d.hxx"
#include "Partition_Loop2d.hxx"
#include "Partition_Loop3d.hxx"
#include "Partition_Spliter.ixx"

#include "utilities.h"

#include <Precision.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>

#include <TopTools_DataMapIteratorOfDataMapOfShapeListOfShape.hxx>
#include <TopTools_DataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapIteratorOfMapOfShape.hxx>
#include <TopTools_SequenceOfShape.hxx>

#include <Geom2d_Curve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>

#include <BRepBndLib.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepLib.hxx>
#include <BRep_Tool.hxx>

#include <Extrema_ExtPC.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <TopOpeBRepTool_CurveTool.hxx>

#ifdef DEB
//# define PART_PERF
#endif

#ifdef PART_PERF
# include <OSD_Chronometer.hxx>
#endif

//=======================================================================
//function : isClosed
//purpose  : check id a shape is closed, ie is a solid or a closed shell
//=======================================================================

static Standard_Boolean isClosed(const TopoDS_Shape& theShape)
{
  Standard_Boolean isClosed = (theShape.ShapeType() == TopAbs_SOLID);

  if (!isClosed && theShape.ShapeType() == TopAbs_SHELL) {
    TopTools_IndexedDataMapOfShapeListOfShape MEF;
    TopExp::MapShapesAndAncestors(theShape, TopAbs_EDGE, TopAbs_FACE, MEF);
    for (Standard_Integer i=1;  isClosed && i<=MEF.Extent();  ++i)
      isClosed = ( MEF(i).Extent() != 1 );
  }
  
  return isClosed;
}

//=======================================================================
//function : Partition_Spliter
//purpose  : constructor
//=======================================================================

Partition_Spliter::Partition_Spliter()
{
  myAsDes = new BRepAlgo_AsDes;
  Clear();
}

//=======================================================================
//function : AddTool
//purpose  : add cutting tool that will _NOT_ be in result
//=======================================================================

void Partition_Spliter::AddTool(const TopoDS_Shape& S)
{
  if (S.ShapeType() < TopAbs_SOLID) { // compound or compsolid
    TopoDS_Iterator it (S);
    for (; it.More(); it.Next())
    {
      AddTool( it.Value());
      myFaceShapeMap.Bind( it.Value(), S ); // to know compound by shape
    }
    return;
  }

  for (TopExp_Explorer exp(S,TopAbs_FACE); exp.More(); exp.Next())
  {
    myMapTools.Add(exp.Current());
    myFaceShapeMap.Bind( exp.Current(), S );
  }
  if (isClosed( S ))
    myClosedShapes.Add( S );
}

//=======================================================================
//function : AddShape
//purpose  : add object Shape to be splited
//=======================================================================

void Partition_Spliter::AddShape(const TopoDS_Shape& S)
{
  if (S.ShapeType() < TopAbs_SOLID) { // compound or compsolid
    TopoDS_Iterator it (S);
    for (; it.More(); it.Next())
    {
      AddShape( it.Value());
      myFaceShapeMap.Bind( it.Value(), S ); // to know compound by shape
    }
    return;
  }

  TopExp_Explorer exp(S,TopAbs_FACE);
  if (!exp.More()) { // do not split edges and vertices
    //myBuilder.Add( myShape, S );
    return;
  }

  Standard_Integer nbFacesBefore = myMapFaces.Extent(); // not to add twice the same S
  for (; exp.More(); exp.Next()) {
    const TopoDS_Shape & aFace = exp.Current();
    if ( ! myFaceShapeMap.IsBound( aFace )) // keep shape of tool face added as object
      myFaceShapeMap.Bind( aFace, S );
    if (myMapFaces.Add( aFace ))
      myImagesFaces.SetRoot( aFace );
  }

  if (nbFacesBefore == myMapFaces.Extent())
    return;

  // solids must be processed before all
  if (S.ShapeType() == TopAbs_SOLID)
    myListShapes.Prepend(S);
  else
    myListShapes.Append(S);

  if (isClosed( S ))
    myClosedShapes.Add( S );

}

//=======================================================================
//function : Shape
//purpose  : return resulting compound
//=======================================================================

TopoDS_Shape Partition_Spliter::Shape() const
{
  return myShape;
}

//=======================================================================
//function : Clear
//purpose  : clear fields
//=======================================================================

void Partition_Spliter::Clear()
{
  myDoneStep = TopAbs_SHAPE;
  
  myListShapes.Clear();
  myMapFaces.Clear();
  myMapTools.Clear();
  myEqualEdges.Clear();
  myNewSection.Clear();
  myClosedShapes.Clear();
  mySharedFaces.Clear();
  myWrappingSolid.Clear();
  myFaceShapeMap.Clear();
  
  myInternalFaces.Clear();
  myIntNotClFaces.Clear();
  
  myAsDes->Clear();
  myImagesFaces.Clear();
  myImagesEdges.Clear();
  myImageShape.Clear();
  
  //  myInter3d = Partition_Inter3d(myAsDes);
  Partition_Inter3d hinter3d (myAsDes);
  myInter3d = hinter3d;
  
  myAddedFacesMap.Clear();

}

//=======================================================================
//function : Compute
//purpose  : produce a result
//=======================================================================

void Partition_Spliter::Compute(const TopAbs_ShapeEnum Limit)
{
  if ((Limit != TopAbs_SHAPE && myDoneStep == Limit) ||
      (Limit == TopAbs_SHAPE && myDoneStep == TopAbs_SOLID))
    return;
  
  myBuilder.MakeCompound( myShape );
  
  TopTools_MapIteratorOfMapOfShape it;
  TopTools_ListIteratorOfListOfShape itl;
  TopExp_Explorer exp;

#ifdef PART_PERF
  OSD_Chronometer aCron;
#endif

  if (myDoneStep > TopAbs_VERTEX) {

    TopTools_ListOfShape aListFaces;
    aListFaces = myImagesFaces.Roots();
    for (it.Initialize(myMapTools); it.More(); it.Next())
      aListFaces.Append(it.Key());

#ifdef PART_PERF
    aCron.Start();
#endif

    //-----------------------------------------------
    // Intersection between faces
    //-----------------------------------------------
    // result is in myAsDes as a map Face - list of new edges;
    // special care is done for section edges, same domain faces and vertices:
    // data about them is inside myInter3d
    myInter3d.CompletPart3d(aListFaces, myFaceShapeMap);

#ifdef PART_PERF
    MESSAGE("+++ CompletPart3d()");
    aCron.Show( cout );
    aCron.Reset();
    aCron.Start();
#endif
    //-----------------------------------------------
    // Intersection of edges
    //-----------------------------------------------

    // add tool faces which must be reconstructed to myMapFaces too
    FindToolsToReconstruct();

#ifdef PART_PERF
    MESSAGE("+++ FindToolsToReconstruct()");
    aCron.Show( cout );
    aCron.Reset();
    aCron.Start();
#endif

    // add existing vertices to edges of object faces in myAsDes
    TopTools_MapOfShape DoneEM;
    for ( it.Initialize(myMapFaces); it.More(); it.Next()) {
      const TopoDS_Shape& F  = it.Key();
      TopoDS_Face FForward = TopoDS::Face(F.Oriented(TopAbs_FORWARD));
      for (exp.Init(FForward,TopAbs_EDGE); exp.More(); exp.Next()) {
        const TopoDS_Edge& E = TopoDS::Edge( exp.Current() );
        myAsDes->Add(FForward,E);
        if (DoneEM.Add(E)) {
          TopoDS_Iterator itV(E);
          for (; itV.More(); itV.Next()) {
            const TopoDS_Vertex& V = TopoDS::Vertex( itV.Value());
            myAsDes->Add(E, myInter3d.ReplaceSameDomainV( V, E ));
          }
        }
      }
    }

    // intersect edges that are descendants of a face in myAsDes
    TopTools_MapOfShape& Modif = myInter3d.TouchedFaces();
    for ( it.Initialize(Modif); it.More(); it.Next()) {
      const TopoDS_Face& F  = TopoDS::Face(it.Key());
      Partition_Inter2d::CompletPart2d (myAsDes, F, myInter3d.NewEdges());
    }
    // now myAsDes contains also new vertices made at edge intersection as
    // descendant of edges both new and old

    myDoneStep = TopAbs_VERTEX;
    
#ifdef PART_PERF
    MESSAGE("+++ CompletPart2d()");
    aCron.Show( cout );
    aCron.Reset();
    aCron.Start();
#endif
  } //   if (myDoneStep > TopAbs_VERTEX)
  
  if (Limit == TopAbs_VERTEX) {
    // add new vertices to myShape
    for ( it.Initialize( myInter3d.NewEdges() ); it.More(); it.Next()) {
      if (! myAsDes->HasDescendant( it.Key() ))
        continue;
      itl.Initialize( myAsDes->Descendant( it.Key() ));
      for (; itl.More(); itl.Next()) 
        myBuilder.Add ( myShape, itl.Value() );
    }
    return;
  }
  

  if (myDoneStep > TopAbs_EDGE) {

    //-----------------------------------------------
    //-----------------------------------------------
    // ------- Reconstruction of all the edges.------
    //-----------------------------------------------
    //-----------------------------------------------

    // ==============
    // cut new edges
    // ==============
    TopTools_ListOfShape LSE; // all edge splits
    for ( it.Initialize(myInter3d.NewEdges()); it.More(); it.Next()) {

      TopoDS_Vertex V1,V2;
      TopoDS_Edge EE = TopoDS::Edge(it.Key());

      TopTools_ListOfShape aListV, aListF;
      aListV = myAsDes->Descendant(EE); // intersection vertices
      aListF = myAsDes->Ascendant(EE);  // intersected faces

      if (aListV.IsEmpty())
        continue;  // new edge does not intersect any other edge

      // Add end vertices to new edges only if 
      // one face is Tool and the other is Shape
      Standard_Boolean isTool1 = ! myMapFaces.Contains( aListF.First() );
      Standard_Boolean isTool2 = ! myMapFaces.Contains( aListF.Last() );
      if (isTool1 || isTool2)
      {
        TopExp::Vertices(EE,V1,V2);
	Standard_Real Tol = Max (BRep_Tool::Tolerance( V1 ),
				 BRep_Tool::Tolerance( V2 ));

        gp_Pnt P1 = BRep_Tool::Pnt(V1);
        gp_Pnt P2 = BRep_Tool::Pnt(V2);
        Standard_Boolean AddV1 = Standard_True;
        Standard_Boolean AddV2 = Standard_True;

        // add only if there is no intersection at end vertex
        for (itl.Initialize(aListV); itl.More(); itl.Next()) {
          const TopoDS_Vertex& Ve = TopoDS::Vertex(itl.Value()) ;
          Standard_Real Tol2 = Max ( Tol, BRep_Tool::Tolerance( Ve ));
          Tol2 *= Tol2;
          gp_Pnt P = BRep_Tool::Pnt(Ve);
          if (AddV1 && P.SquareDistance(P1) <= Tol2)
            AddV1 = Standard_False;

          if (AddV2 && P.SquareDistance(P2) <= Tol2) 
            AddV2 = Standard_False;
        }

        if (AddV1) {
          aListV.Append(V1);
          myAsDes->Add(EE,V1);
        }

        if (AddV2) {
          aListV.Append(V2);
          myAsDes->Add(EE,V2);
        }
      }

      // cut new edges
      Standard_Integer NbV=aListV.Extent() ;
      if (NbV>1 || (NbV==1 && V1.IsSame(V2)) ) {
        TopTools_ListOfShape LNE;
        MakeEdges (EE,aListV, LNE);
        myImagesEdges.Bind(EE,LNE);
	LSE.Append( LNE );
      }
    }

    // ==============
    // cut old edges
    // ==============
    for ( it.Initialize(myMapFaces); it.More(); it.Next()) {
      for (exp.Init( it.Key(), TopAbs_EDGE); exp.More(); exp.Next()) {
        const TopoDS_Edge& EE = TopoDS::Edge( exp.Current() );
        if ( myImagesEdges.HasImage( EE ))
          continue;
        TopTools_ListOfShape  LNE;
        const TopTools_ListOfShape& aListVV = myAsDes->Descendant(EE);
        MakeEdges (EE, aListVV, LNE);
        myImagesEdges.Bind(EE,LNE);
	LSE.Append( LNE );
      }
    }
#ifdef PART_PERF
    MESSAGE("+++ Cut Edges");
    aCron.Show( cout );
    aCron.Reset();
    aCron.Start();
#endif

    // process same domain section edges
    MergeEqualEdges( LSE );
    
    myDoneStep = TopAbs_EDGE;
    
#ifdef PART_PERF
    MESSAGE("+++ MergeEqualEdges()");
    aCron.Show( cout );
    aCron.Reset();
    aCron.Start();
#endif
  }  //   if (myDoneStep > TopAbs_EDGE) 

  if (Limit == TopAbs_EDGE) {
    // add splits of old edges
    TopTools_ListIteratorOfListOfShape itNE;
    for (itl.Initialize( myListShapes );itl.More();itl.Next()) {
      if (myMapTools.Contains( itl.Value() ))
        continue; // skip tool faces
      for ( exp.Init( itl.Value(), TopAbs_EDGE ); exp.More(); exp.Next()) {
	itNE.Initialize( myImagesEdges.Image( exp.Current() ));
	for ( ; itNE.More(); itNE.Next())
	  myBuilder.Add ( myShape, itNE.Value() );
      }
    }
    // add splits of new edges
    for ( it.Initialize( myInter3d.NewEdges() ); it.More(); it.Next()) {
      itNE.Initialize( myImagesEdges.Image( it.Key() ));
      for (; itNE.More(); itNE.Next())
        myBuilder.Add ( myShape, itNE.Value() );
    }
    return;
  }
  
  
  //-----------------------------------------------
  // split faces
  //-----------------------------------------------

  if (myDoneStep > TopAbs_FACE) {
    
    for (itl.Initialize(myListShapes);itl.More();itl.Next()) {
      TopoDS_Shape FacesComp = MakeFaces ( itl.Value());
      // there is a cunning here: myImagesFaces keeps faces made by Loop2d
      // but some of them may be replaced with splits of same domain face
      // and myImageShape keeps ultimate result
      myImageShape.Bind( itl.Value(), FacesComp );
    }
    
    myDoneStep = TopAbs_FACE;
#ifdef PART_PERF
    MESSAGE("+++ MakeFaces()");
    aCron.Show( cout );
    aCron.Reset();
    aCron.Start();
#endif
  }
  
  if (Limit == TopAbs_WIRE ||
      Limit == TopAbs_FACE)   {
    for (itl.Initialize(myListShapes);itl.More();itl.Next()) {
      if ( myMapTools.Contains( itl.Value() ))
	continue; // no result needed for a tool face
      const TopoDS_Shape& FacesComp = myImageShape.Image( itl.Value() ).First();
      for ( exp.Init( FacesComp, Limit); exp.More(); exp.Next())
	myBuilder.Add ( myShape, exp.Current());
    }
    return;
  }

  
  //-----------------------------------------------
  // split and add solids and shells
  //-----------------------------------------------

  Standard_Boolean makeSolids = (Limit == TopAbs_SHAPE ||
				 Limit < TopAbs_SHELL);
  for (itl.Initialize(myListShapes);itl.More();itl.Next())
  {
    const TopoDS_Shape & S = itl.Value();
    if (S.ShapeType() > TopAbs_SHELL)
      continue;

    TopTools_ListOfShape NSL; // new shape list
    MakeShells (S , NSL);
    if (makeSolids && S.ShapeType() == TopAbs_SOLID )
      MakeSolids( S, NSL );

    // store new shells or solids
    TopTools_ListIteratorOfListOfShape itNSL (NSL);
    for ( ; itNSL.More(); itNSL.Next()) 
      myBuilder.Add (myShape, itNSL.Value());
  }
#ifdef PART_PERF
    MESSAGE("+++ MakeShells()");
    aCron.Show( cout );
#endif

  //-----------------------------------------------
  // add split faces
  //-----------------------------------------------

  for (itl.Initialize(myListShapes);itl.More();itl.Next())
  {
    const TopoDS_Shape & S = itl.Value();
    if (S.ShapeType() != TopAbs_FACE ||
        myMapTools.Contains( S ))
      continue; 
    TopoDS_Iterator itS( myImageShape.Image(S).First() );
    for (; itS.More(); itS.Next())
      if (! myAddedFacesMap.Contains( itS.Value() ))
        myBuilder.Add (myShape, itS.Value());
  }

  myDoneStep = makeSolids ? TopAbs_SOLID : TopAbs_SHELL;
  
}

//=======================================================================
//function : MakeSolids
//purpose  : make solids out of Shells
//=======================================================================

void Partition_Spliter::MakeSolids(const TopoDS_Shape &   theSolid,
                                   TopTools_ListOfShape & theShellList)
{
  // for a solid wrapping other shells or solids without intersection,
  // it is necessary to find shells making holes in it

  TopTools_ListOfShape aNewSolids; // result
  TopTools_ListOfShape aHoleShells;
  TopoDS_Shape anInfinitePointShape;

  Standard_Boolean isWrapping = myWrappingSolid.Contains( theSolid );
  if (!isWrapping && !theShellList.IsEmpty())
  {
    // check if theSolid initially has internal shells
    TopoDS_Iterator aShellExp (theSolid);
    aShellExp.Next();
    isWrapping = aShellExp.More();
  }
  
  TopTools_ListIteratorOfListOfShape aShellIt(theShellList);
  for ( ; aShellIt.More(); aShellIt.Next())
  {
    const TopoDS_Shape & aShell = aShellIt.Value();

    // check if a shell is a hole
    if (isWrapping && IsInside (anInfinitePointShape, aShell))
      aHoleShells.Append( aShell );
    else
    {
      // make a solid from a shell
      TopoDS_Solid Solid;
      myBuilder.MakeSolid( Solid );
      myBuilder.Add (Solid, aShell);

      aNewSolids.Append (Solid);
    }
  }

  // find an outer shell most close to each hole shell
  TopTools_DataMapOfShapeShape aInOutMap;
  for (aShellIt.Initialize( aHoleShells ); aShellIt.More(); aShellIt.Next())
  {
    const TopoDS_Shape & aHole = aShellIt.Value();
    TopTools_ListIteratorOfListOfShape aSolisIt (aNewSolids);
    for ( ; aSolisIt.More(); aSolisIt.Next())
    {
      const TopoDS_Shape & aSolid = aSolisIt.Value();
      if (! IsInside( aHole, aSolid ))
        continue;

      if ( aInOutMap.IsBound (aHole))
      {
        const TopoDS_Shape & aSolid2 = aInOutMap( aHole );
        if ( IsInside( aSolid, aSolid2 ))
        {
          aInOutMap.UnBind( aHole );
          aInOutMap.Bind ( aHole, aSolid );
        }
      }
      else
        aInOutMap.Bind ( aHole, aSolid );
    }

    // add aHole to a solid
    if (aInOutMap.IsBound( aHole ))
      myBuilder.Add ( aInOutMap( aHole ), aHole );
  }

  theShellList.Clear();
  theShellList.Append( aNewSolids );
}
 
//=======================================================================
//function : FindFacesInside
//purpose  : return compound of faces  of other shapes that are
//           inside <theShape>. 
//           <theShape> is an object shape.
//           <CheckClosed> makes avoid faces that do not form a
//           closed shell
//           <All> makes return already added faces
//=======================================================================

TopoDS_Shape Partition_Spliter::FindFacesInside(const TopoDS_Shape& theShape,
						const Standard_Boolean CheckClosed,
						const Standard_Boolean All)
{
  // ================================================
  // check if internal faces have been already found
  // ================================================
  TopExp_Explorer expl;
  if (myInternalFaces.IsBound( theShape ))
  {
    TopoDS_Shape aIntFComp = myInternalFaces.Find ( theShape );
    TopoDS_Shape aIntRemFComp = myIntNotClFaces.Find ( theShape );

    expl.Init( aIntRemFComp, TopAbs_FACE);
    if (CheckClosed || !expl.More())
      return aIntFComp;

    TopoDS_Compound C;
    myBuilder.MakeCompound( C );
    // add removed faces
    for (; expl.More(); expl.Next())
      myBuilder.Add( C, expl.Current() );
    // add good internal faces
    for (expl.Init( aIntFComp, TopAbs_FACE); expl.More(); expl.Next())
      myBuilder.Add( C, expl.Current() );
    return C;
  }

  // ===================================
  // get data for internal faces search
  // ===================================

  // compound of split faces of theShape 
  const TopoDS_Shape& CSF = myImageShape.Image(theShape).First();

  TopTools_MapOfShape MSE, MFP;
  TopTools_DataMapOfShapeListOfShape DMSEFP;
  TopTools_MapIteratorOfMapOfShape itm;
  TopTools_ListOfShape EmptyL;

  // MSE filling: map of new section edges of CSF
  for (expl.Init(CSF,TopAbs_EDGE); expl.More(); expl.Next()) {
    const TopoDS_Shape & resE = expl.Current() ;
    if (myNewSection.Contains( resE )) // only new edges
      MSE.Add(resE);
  }

  // DMEF: map edge of CSF - faces of CSF
  TopTools_IndexedDataMapOfShapeListOfShape DMEF;
  TopExp::MapShapesAndAncestors(CSF, TopAbs_EDGE, TopAbs_FACE, DMEF);

  // Fill
  // 1.  MFP - a map of faces to process: map of resulting faces except
  // those of theShape; we`ll add to C those of them which are inside CSF
  // 2.  DMSEFP - edge of MSE => faces of MFP
  TopTools_ListIteratorOfListOfShape itl;
  for (itl.Initialize(myListShapes);itl.More();itl.Next()) {
    const TopoDS_Shape& aShape = itl.Value();
    if ( theShape.IsSame( aShape )) continue;
    // fill maps
    // iterate on split faces of aShape
    TopoDS_Iterator itF ( myImageShape.Image(aShape).First() );
    for ( ; itF.More(); itF.Next()) {
      const TopoDS_Shape& sf = itF.Value();
      MFP.Add(sf);
      // iterate on edges of split faces of aShape,
      // add to DMSEFP edges that are new
      for (expl.Init( sf, TopAbs_EDGE ); expl.More(); expl.Next()) {
	TopoDS_Shape se = expl.Current();
	if ( MSE.Contains(se)) {// section edge
	  if (!DMSEFP.IsBound(se)) 
	    DMSEFP.Bind(se,EmptyL);
	  DMSEFP(se).Append(sf);
	}
      }
    }
  }

  // add tool faces having section edges on faces of theShape to MFP and DMSEFP;
  // such tool faces need not to be reconstructed and so they are not in myListShapes
  for (itm.Initialize(myMapTools); itm.More(); itm.Next())
  {
    const TopoDS_Shape & aToolFace = itm.Key();
    if (myMapFaces.Contains( aToolFace ))
      continue;
    MFP.Add(aToolFace);
    for (expl.Init( aToolFace, TopAbs_EDGE ); expl.More(); expl.Next()) {
      TopoDS_Shape se = expl.Current();
      if ( MSE.Contains( se )) {// section edge
        if (!DMSEFP.IsBound( se )) 
          DMSEFP.Bind( se, EmptyL );
        DMSEFP( se ).Append( aToolFace );
      }
    }
  }
  

  // ===========================
  // find faces inside theShape
  // ===========================

  Standard_Boolean skipAlreadyAdded = Standard_False;
  Standard_Boolean GoodOri, inside;
  Standard_Real dot;
  TopTools_ListOfShape KeepFaces;
  TopTools_DataMapIteratorOfDataMapOfShapeListOfShape Mapit;

  // iterate on section edges, check faces of other shapes
  // sharing section edges and put internal faces to KeepFaces
  for (Mapit.Initialize(DMSEFP); Mapit.More() ; Mapit.Next() ) {
    // a new edge of theShape
    const TopoDS_Edge& E = TopoDS::Edge (Mapit.Key());
    // an original edge of which E is a split
    const TopoDS_Edge& OrigE = TopoDS::Edge ( myImagesEdges.Root( E ));
    // does OrigE itself splits a face
    Standard_Boolean isSectionE = myInter3d.IsSectionEdge ( OrigE );

    // split faces of other shapes sharing E
    TopTools_ListOfShape& LSF = DMSEFP.ChangeFind(E);
    itl.Initialize( LSF );
    while (itl.More()) {
      // a split faces of other shape
      TopoDS_Face aFace1 = TopoDS::Face(itl.Value());
      // remove aFace1 form DMSEFP and MFP
      LSF.Remove( itl ); // == itl.Next();
      if (!MFP.Remove( aFace1 ))
	continue; // was not is MFP ( i.e already checked)
      // check if aFace1 was already added to 2 shells
      if (!All &&
	  myAddedFacesMap.Contains( aFace1 ) &&
	  myAddedFacesMap.Contains( aFace1.Reversed() )) {
	skipAlreadyAdded = Standard_True;
	continue;
      }

      // find another face which originates from the same face as aFace1:
      // usually aFace2 is internal if aFace1 is not and vice versa

      TopoDS_Shape anOrigFace = aFace1;
      if (myImagesFaces.IsImage(aFace1))
        anOrigFace = myImagesFaces.Root(aFace1);
      TopoDS_Shape aFace2;
      if ( !isSectionE ) {
        while (itl.More()) {
          aFace2 = itl.Value();
          if (!MFP.Contains( aFace2 )) {
            LSF.Remove( itl );
            continue;
          }
          if (anOrigFace.IsSame( myImagesFaces.Root( aFace2 )))
            break;
          itl.Next();
        }
        if (itl.More()) { // aFace2 found, remove it from maps
          LSF.Remove( itl );
          MFP.Remove(aFace2);
        }
        else
          aFace2.Nullify();
        itl.Initialize( LSF );
      }

      // check that anOrigFace is not same domain with CSF faces it intersects

      const TopTools_ListOfShape& FL = DMEF.FindFromKey(E); //faces of CSF sharing E
      const TopoDS_Shape& origF1 = myImagesFaces.Root(FL.First());
      const TopoDS_Shape& origF2 = myImagesFaces.Root(FL.Last());
      Standard_Boolean sameDom1 = anOrigFace.IsSame( origF1 );
      Standard_Boolean sameDom2 = anOrigFace.IsSame( origF2 );
      if (!(sameDom1 || sameDom2) && myInter3d.HasSameDomainF( anOrigFace )) {
	sameDom1 = myInter3d.IsSameDomainF( anOrigFace, origF1);
        if (origF1 == origF2)
          sameDom2 = sameDom1;
        else
          myInter3d.IsSameDomainF( anOrigFace, origF2);
      }
      if (sameDom1 && sameDom2)
	continue;
      if ((sameDom1 || sameDom2)) {
	inside = Partition_Loop3d::IsInside (E,
					     TopoDS::Face(FL.First()),
					     TopoDS::Face(FL.Last()),
					     1, dot, GoodOri);
	if (inside || (dot + Precision::Angular() >= 1.0))
	  continue; // E is convex between origF1 and origF2 or they are tangent
      }


      // keep one of found faces

      //face of CSF sharing E
      const TopoDS_Shape& aShapeFace = sameDom1 ? FL.Last() : FL.First();
      // analyse aFace1 state
      inside = Partition_Loop3d::IsInside (E, TopoDS::Face(aShapeFace), aFace1,
					   1, dot, GoodOri);
      if (inside && isSectionE)
      {
        // aFace1 must be tested with both adjacent faces of CSF
        const TopoDS_Shape& aShapeFace2 = sameDom1 ? FL.First() : FL.Last();
        if (aShapeFace2 != aShapeFace)
          inside = Partition_Loop3d::IsInside (E, TopoDS::Face(aShapeFace2), aFace1,
                                               1, dot, GoodOri);
      }

      // store internal face
      if (inside)
        KeepFaces.Append(aFace1);

      else if (!aFace2.IsNull())
      {
        if (dot + Precision::Angular() >= 1.0)
        {
          // aFace2 state is not clear, it will be analysed alone,
          // put it back to the maps
          MFP.Add( aFace2 );
          LSF.Append( aFace2 );
        }
        else
          KeepFaces.Append(aFace2);
      }
    }
  }

  // ===================================================
  // add not distributed faces connected with KeepFaces
  // ===================================================

  // ultimate list of internal faces
  TopTools_ListOfShape KeptFaces;

  // add to MFP not split tool faces as well, they may be connected with
  // tool faces interfering with theShape
  for ( itm.Initialize(myMapTools); itm.More(); itm.Next() ) {
    const TopoDS_Shape& aToolFace = itm.Key();
    if (!myImageShape.HasImage(aToolFace))
      MFP.Add (aToolFace);
  }

  if (MFP.IsEmpty())
    KeptFaces.Append (KeepFaces);

  while (!KeepFaces.IsEmpty())
  {
    // KeepEdges : map of edges of faces kept last time
    TopTools_IndexedMapOfShape KeepEdges;
    for ( itl.Initialize(KeepFaces); itl.More(); itl.Next() ) {
      TopExp::MapShapes( itl.Value(), TopAbs_EDGE, KeepEdges);
      KeptFaces.Append( itl.Value() );
    }

    KeepFaces.Clear();

    // keep faces connected with already kept faces by KeepEdges
    for ( itm.Initialize(MFP); itm.More(); itm.Next() ) {
      const TopoDS_Shape& FP = itm.Key();
      for (expl.Init(FP,TopAbs_EDGE); expl.More(); expl.Next()) {
        const TopoDS_Shape& se = expl.Current();
        if (!MSE.Contains(se) && KeepEdges.Contains(se) ) {
          KeepFaces.Append(FP);
          MFP.Remove(FP);
          break;
        }
      }
    }
  }

  // ===============================================================
  // here MFP contains faces outer of theShape and those of shapes
  // which do not interfere with theShape at all and between which
  // there may be those wrapped by theShape and whose faces may be
  // needed to be returned as well
  // ===============================================================

  Standard_Boolean isSolid = (theShape.ShapeType() == TopAbs_SOLID);
  if (All || isSolid)  // All is for sub-result removal
  {
    // loop on not used faces; checked faces will be removed from MFP
    // during the loop
    for ( itm.Initialize( MFP ); itm.More(); itm.Next() ) {
      const TopoDS_Shape & aFace = itm.Key();

      // a shape which aFace originates from
      TopoDS_Shape anOrigShape = GetOriginalShape( aFace );

      // find out if all split faces of anOrigShape are not in MFP
      // and by the way remove them from MFP
      Standard_Boolean isAllOut = Standard_True;
      TopoDS_Shape aSplitFaces = anOrigShape;
      if (myImageShape.HasImage(anOrigShape))
        aSplitFaces = myImageShape.Image(anOrigShape).First();

      TopTools_ListOfShape aSplitFaceL; // faces candidate to be kept
      for (expl.Init( aSplitFaces, TopAbs_FACE ); expl.More(); expl.Next())
      {
        const TopoDS_Shape & aSpFace = expl.Current();
        // a tool face which became object has image but the whole tool shape has not
        if (myImageShape.HasImage( aSpFace ))
        {
          TopExp_Explorer exF (myImageShape.Image( aSpFace ).First(), TopAbs_FACE );
          for ( ; exF.More(); exF.Next() )
          {
            aSplitFaceL.Append( exF.Current() );
            if ( ! MFP.Remove( exF.Current() ) && isAllOut )
              // a shared face might be removed from MFP during a prev loop
              isAllOut = mySharedFaces.Contains( exF.Current() );
          }
        }
        else
        {
          aSplitFaceL.Append( aSpFace );
          if ( ! MFP.Remove( aSpFace ) && isAllOut)
            // a shared face might be removed from MFP during a prev loop
            isAllOut = mySharedFaces.Contains( aSpFace );
        }
      }
      itm.Initialize( MFP ); // iterate remaining faces
      if ( !isAllOut )
        continue;

      // classify anOrigShape against theShape
      if (IsInside (anOrigShape, theShape))
      {
        if (isSolid && myClosedShapes.Contains( anOrigShape ))
          // to make a special care at solid reconstruction
          myWrappingSolid.Add ( theShape );

        // keep faces of an internal shape anOrigShape
        KeptFaces.Append( aSplitFaceL );
      }
    }
  }

  // ====================================================
  // check if kept faces form a shell without free edges
  // ====================================================

  DMEF.Clear();  // edge - kept faces
  MFP.Clear(); // reuse it for wrong faces
  if (CheckClosed) {
    for (itl.Initialize(KeptFaces); itl.More(); itl.Next() ) 
      TopExp::MapShapesAndAncestors(itl.Value(), TopAbs_EDGE, TopAbs_FACE, DMEF);

    Standard_Integer i, nb = DMEF.Extent();
    Standard_Boolean isClosed = Standard_False;
    while (!isClosed) {
      isClosed = Standard_True;
      for (i=1;  isClosed && i<=nb;  ++i) {
        const TopoDS_Shape& E = DMEF.FindKey( i );
        if (! BRep_Tool::Degenerated( TopoDS::Edge( E )) &&
            ! MSE.Contains( E ))
          isClosed = ( DMEF(i).Extent() != 1 );
      }
      if (!isClosed) {
        const TopoDS_Shape& F = DMEF.FindFromIndex( i-1 ).First(); // bad face
        MFP.Add( F ); 
        // remove bad face from DMEF
        for (expl.Init( F, TopAbs_EDGE); expl.More(); expl.Next()) {
	  const TopoDS_Shape& E = expl.Current();
          TopTools_ListOfShape& FL = DMEF.ChangeFromKey( E );
          for (itl.Initialize( FL ); itl.More(); itl.Next() ) {
            if ( F.IsSame( itl.Value() )) {
              FL.Remove( itl );
              break;
            }
          }
        }
      }
    }
  }

  // ==============
  // make a result
  // ==============

  TopoDS_Compound C;
  // compound of removed internal faces
  TopoDS_Compound CNotCl;

  myBuilder.MakeCompound(C);
  myBuilder.MakeCompound(CNotCl);

  // add to compounds
  for (itl.Initialize(KeptFaces); itl.More(); itl.Next() )
  {
    TopoDS_Shape & aIntFace = itl.Value();
    if (! MFP.Contains( aIntFace ))
      myBuilder.Add( C, aIntFace);
    else
      myBuilder.Add( CNotCl, aIntFace);
  }

  if (!skipAlreadyAdded && CheckClosed)
  {
    myInternalFaces.Bind( theShape, C );
    myIntNotClFaces.Bind( theShape, CNotCl );
  }

  return C;
}

//=======================================================================
//function : MakeShell
//purpose  : split S into compound of shells
//=======================================================================

void Partition_Spliter::MakeShells(const TopoDS_Shape& S,
                                   TopTools_ListOfShape& NS)
{
  Partition_Loop3d ShellMaker;
  // get compound of split faces of S
  const TopoDS_Shape& FacesComp = myImageShape.Image(S).First();
  ShellMaker.AddConstFaces( FacesComp );
  // add split faces inside S
  if (myClosedShapes.Contains( S )) {
    TopoDS_Shape InternalFacesComp = FindFacesInside(S, Standard_True);
    ShellMaker.AddSectionFaces( InternalFacesComp );
  }
  
  NS = ShellMaker.MakeShells( myAddedFacesMap );

  // Add faces added to new shell to myAddedFacesMap:
  // avoid rebuilding twice commont part of 2 solids.
  TopTools_ListIteratorOfListOfShape itS(NS);
  while ( itS.More()) {
    TopExp_Explorer expF (itS.Value(), TopAbs_FACE);
    for (; expF.More(); expF.Next())
      myAddedFacesMap.Add (expF.Current());
    
    itS.Next();
  }
}

//=======================================================================
//function : findEqual
//purpose  : compare edges of EL1 against edges of EL2,
//           Result is in EMM binding EL1 edges to list of equal edges.
//           Edges are considered equall only if they have same vertices.
//           <addSame>==True makes consider same edges as equal
//           Put in <AllEqMap> all equal edges
//=======================================================================

static void findEqual (const TopTools_ListOfShape& EL1,
		       const TopTools_ListOfShape& EL2,
		       const Standard_Boolean addSame,
		       TopTools_DataMapOfShapeListOfShape& EEM,
		       TopTools_MapOfShape& AllEqMap)
{
  // map vertices to edges for EL2
  TopTools_DataMapOfShapeListOfShape VEM;
  TopTools_ListIteratorOfListOfShape itE1, itE2(EL2);
  TopoDS_Iterator itV;
  TopTools_ListOfShape emptyL;
  for (; itE2.More(); itE2.Next()) {
    for (itV.Initialize( itE2.Value() ); itV.More(); itV.Next()) {
      const TopoDS_Shape& V = itV.Value(); 
      if (! VEM.IsBound( V ) )
	VEM.Bind( V, emptyL);
      VEM( V ).Append( itE2.Value());
    }
  }

  gp_Vec D1, D2;
  gp_Pnt P;
  Standard_Real f,l,u,tol;
  Handle(Geom_Curve) C1, C2;
  Extrema_ExtPC Extrema;
  TopoDS_Vertex V1, V2, V3, V4;

  AllEqMap.Clear();
  
  for (itE1.Initialize(EL1); itE1.More(); itE1.Next()) {
    const TopoDS_Edge& E1 = TopoDS::Edge( itE1.Value());
    if (BRep_Tool::Degenerated( E1 ) || AllEqMap.Contains (E1))
      continue;
    TopExp::Vertices( E1, V1, V2 );

    if (VEM.IsBound(V1))
      itE2.Initialize( VEM(V1) );
    for (; itE2.More(); itE2.Next()) {
      const TopoDS_Edge& E2 = TopoDS::Edge( itE2.Value());
      if (BRep_Tool::Degenerated( E2 ) || AllEqMap.Contains (E2))
        continue;

      if (E1.IsSame(E2)) {
	if (!addSame)
	  continue;
      }
      else {
	TopExp::Vertices( E2, V3, V4);
	if (!V2.IsSame(V4) && !V2.IsSame(V3))
	  continue;
	// E1 and E2 have same vertices
	// check D1 at end points.
        C2 = BRep_Tool::Curve( E2, f,l);
        C1 = BRep_Tool::Curve( E1, f,l);
	u = BRep_Tool::Parameter(V1,E1);
        C1->D1(u, P, D1);
	u = BRep_Tool::Parameter(V1.IsSame(V3) ? V3 : V4, E2);
	C2->D1(u, P, D2);
        D1.Normalize(); D2.Normalize();
        if (Abs(D1*D2) + Precision::Angular() < 1.0)
          continue;
	if (! V1.IsSame(V2)) {
	  u = BRep_Tool::Parameter(V2,E1);
	  C1->D1(u, P, D1);
	  u = BRep_Tool::Parameter(V2.IsSame(V3) ? V3 : V4, E2);
	  C2->D1(u, P, D2);
	  D1.Normalize(); D2.Normalize();
	  if (Abs(D1*D2) + Precision::Angular() < 1.0)
	    continue;
	}
        // check distance at a couple of internal points
        tol = Max(BRep_Tool::Tolerance(E1),
                  BRep_Tool::Tolerance(E2));
        GeomAdaptor_Curve AC1(C1);
        Extrema.Initialize(AC1,f,l);
	Standard_Boolean ok = Standard_True, hasMin = Standard_False;
	BRep_Tool::Range( E2, f, l);
        Standard_Integer i=1, nbi=3;
        for (; i<nbi && ok; ++i) {
          Extrema.Perform( C2->Value( f+(l-f)*i/nbi ));
          Standard_Integer j=1, nbj=Extrema.NbExt();
          for (; j<=nbj && ok; ++j) {
            if (Extrema.IsMin(j)) {
	      hasMin = Standard_True;
	      ok = Extrema.Value(j) <= tol;  // V6.3
	      // ok = Extrema.SquareDistance(j) <= tol;  // V6.5
	    }
          }
        }
        if ( !hasMin || !ok)
          continue;
      }
      // bind E2 to E1 in EEM
      if (!EEM.IsBound(E1)) {
        EEM.Bind (E1, emptyL);
	AllEqMap.Add (E1);
      }
      EEM(E1).Append(E2);
      AllEqMap.Add (E2);
    }
  }
}

//=======================================================================
//function : MakeFaces
//purpose  : split faces of S, return compound of new faces
//=======================================================================

TopoDS_Shape Partition_Spliter::MakeFaces (const TopoDS_Shape& S)
{
  TopoDS_Compound C;
  myBuilder.MakeCompound(C);
  
  TopTools_ListIteratorOfListOfShape itl, itNE;
  
  TopExp_Explorer exp(S,TopAbs_FACE);
  for (; exp.More(); exp.Next()) {

    const TopoDS_Face& F = TopoDS::Face(exp.Current());

    TopTools_ListOfShape LNF;
    
    if (myImagesFaces.HasImage( F )) {
      myImagesFaces.LastImage( F, LNF );
      TopAbs_Orientation oriF = F.Orientation();
      for ( itl.Initialize( LNF ); itl.More(); itl.Next())
	itl.Value().Orientation( oriF );
    }
    else {

      Partition_Loop2d loops;
      loops.Init(F);

      TopTools_IndexedMapOfShape EM;
      TopExp::MapShapes( F, TopAbs_EDGE, EM);

      TopTools_MapOfShape AddedEqualM, EqualSeamM;
      Standard_Boolean needRebuild = Standard_False;

      // add splits to loops

      // LE: old edges + new not splitted edges
      const TopTools_ListOfShape& LE = myAsDes->Descendant(F);
      for (itl.Initialize(LE); itl.More(); itl.Next()) {
	const TopoDS_Edge& E = TopoDS::Edge( itl.Value() );

	Standard_Boolean isSectionE = myInter3d.IsSectionEdge(E);
	Standard_Boolean isNewE = !EM.Contains( E );

	// LSE: list of split edges
	TopTools_ListOfShape LSE;
	myImagesEdges.LastImage(E,LSE); // splits of E or E itself

	for (itNE.Initialize(LSE); itNE.More(); itNE.Next()) {

	  TopoDS_Edge NE = TopoDS::Edge( itNE.Value() );
	  Standard_Boolean isSameE = NE.IsSame ( E );
	  
	  if ( isNewE || isSectionE || !isSameE) {
	    if (AddedEqualM.Contains( NE )) {
              // a seam must be twice in a loop
              if (!BRep_Tool::IsClosed( E, F ) || !EqualSeamM.Add( NE ))
                continue;
            }

	    if (isNewE) {
	      if (isSectionE) {
		if ( ! myInter3d.IsSplitOn( NE, E, F) )
		  continue;
	      }
	      else {
		TopoDS_Vertex V1,V2;
		TopExp::Vertices(NE,V1,V2);
		const TopTools_ListOfShape& EL1 = myAsDes->Ascendant(V1);
		const TopTools_ListOfShape& EL2 = myAsDes->Ascendant(V2);
		if ( EL1.Extent() < 2 && EL2.Extent() < 2 )
		  continue;
	      }
	    }
	    else {
	      NE.Orientation( E.Orientation());
	      if (!isSameE) {
		// orient NE because it may be a split of other edge
		Standard_Real f,l,u;
		Handle(Geom_Curve) C3d  = BRep_Tool::Curve( E,f,l );
		Handle(Geom_Curve) NC3d = BRep_Tool::Curve( NE,f,l);
		if ( C3d != NC3d) {
		  gp_Vec D1, ND1;  gp_Pnt P;
		  TopoDS_Vertex V = TopExp::FirstVertex(NE);
		  u = BRep_Tool::Parameter(V,NE);
		  NC3d->D1 (u, P, ND1);
		  u = BRep_Tool::Parameter(V,E);
		  C3d ->D1 (u, P, D1);
		  if (ND1.Dot(D1) < 0)
		    NE.Reverse();
		}
	      }
	    }
	    if (myEqualEdges.Contains( NE ))
              AddedEqualM.Add( NE );

	    needRebuild = Standard_True;
	  }

	  if (isNewE || isSectionE)
	    myNewSection.Add( NE );

	  if (isNewE) 
	    loops.AddSectionEdge(NE);
	  else
	    loops.AddConstEdge(NE);
	}
      }

      //-------------------
      // Build the faces.
      //-------------------
      
      if (needRebuild) {
	
        loops.Perform();
        loops.WiresToFaces(myImagesEdges);

        LNF = loops.NewFaces();

        myImagesFaces.Bind(F,LNF);

        // replace the result faces that have already been built
        // during same domain faces reconstruction done earlier
        if (myInter3d.HasSameDomainF( F ))
        {
          // build map edge to same domain faces: EFM
          TopTools_IndexedDataMapOfShapeListOfShape EFM;
          TopTools_MapOfShape SDFM; // avoid doubling
          itl.Initialize( myInter3d.SameDomain( F ));
          for (; itl.More(); itl.Next()) {
            if ( !myImagesFaces.HasImage( itl.Value() ))
              continue;
            // loop on splits of a SD face
            TopTools_ListIteratorOfListOfShape itNF;
            itNF.Initialize (myImagesFaces.Image( itl.Value() ));
            for ( ; itNF.More(); itNF.Next()) {
              TopoDS_Shape SDF = itNF.Value();
              if (myImagesFaces.HasImage( SDF )) // already replaced
                SDF = myImagesFaces.Image( SDF ).First();
              if (SDFM.Add (SDF))
                TopExp::MapShapesAndAncestors(SDF, TopAbs_EDGE, TopAbs_FACE, EFM);
            }
          }
          // do replace faces in the LNF
          TopTools_ListOfShape LOF;
          if ( !EFM.IsEmpty() )
            itl.Initialize( LNF );
          while (itl.More()) {
            const TopoDS_Shape& NF = itl.Value();
            TopExp_Explorer expE ( NF, TopAbs_EDGE );
            const TopoDS_Edge& E  = TopoDS::Edge (expE.Current());
            if (EFM.Contains( E )) {
              const TopTools_ListOfShape& SDFL = EFM.FindFromKey( E );
              TopoDS_Shape SDF = SDFL.First();
              Standard_Boolean GoodOri;
              Standard_Real dot;
              Partition_Loop3d::IsInside (E, TopoDS::Face(NF), TopoDS::Face(SDF),
                                          1, dot, GoodOri);
              if (dot < 0)
              {
                // NF and SDF are on different side of E
                if (SDFL.Extent() == 1) {
                  itl.Next();
                  continue;
                }
                else
                  SDF = SDFL.Last(); // next face must be on the same side
              }
              gp_Vec V1 = Partition_Loop3d::Normal( E, TopoDS::Face( NF ));
              gp_Vec V2 = Partition_Loop3d::Normal( E, TopoDS::Face( SDF ));
              if (V1*V2 < 0)
                SDF.Reverse();

              if (!myImagesFaces.HasImage( NF ))
                myImagesFaces.Bind( NF, SDF );

              // mySharedFaces is used in FindFacesInside()
              mySharedFaces.Add( SDF );

              LOF.Prepend ( SDF );
              LNF.Remove (itl);
            }
            else
              itl.Next();
          }

          LNF.Append (LOF);
        }
      } // if (needRebuild)
      
      else {
	LNF.Append( F );
	myImagesFaces.Bind(F,LNF);
      }
    } // if (myImagesFaces.HasImage( F ))

    // fill the resulting compound
    for (itl.Initialize(LNF); itl.More(); itl.Next())
      myBuilder.Add ( C, itl.Value());
    
  } // loop on faces of S

  return C;
}


//=======================================================================
//function : Tri
//purpose  : 
//=======================================================================

static void Tri(const TopoDS_Edge&        E,
		TopTools_SequenceOfShape& Seq,
                const Partition_Inter3d & theInter3d)
{
  Standard_Boolean Invert   = Standard_True;
  Standard_Real    U1,U2;
  TopoDS_Vertex    V1,V2;

  while (Invert) {
    Invert = Standard_False;
    for ( Standard_Integer i = 1; i < Seq.Length(); i++) {
      
      V1 = TopoDS::Vertex(Seq.Value(i));
      V2 = TopoDS::Vertex(Seq.Value(i+1));
      
      V1.Orientation(TopAbs_INTERNAL);
      V2.Orientation(TopAbs_INTERNAL);
      
      U1 = BRep_Tool::Parameter(V1,E);
      U2 = BRep_Tool::Parameter(V2,E);
      
      if (IsEqual(U1,U2)) {
        if (theInter3d.ReplaceSameDomainV( V1, E ).IsSame( V1 ))
          Seq.Remove(i+1); // remove V2
        else
          Seq.Remove(i);
	i--;
	continue;
      }
      if (U2 < U1) {
	Seq.Exchange(i,i+1);
	Invert = Standard_True;
      }
    }
  }
}

//=======================================================================
//function : MakeEdges
//purpose  : cut E by vertices VOnE, return list of new edges NE
//=======================================================================

void Partition_Spliter::MakeEdges (const TopoDS_Edge& E,
                                   const TopTools_ListOfShape& VOnE,
                                   TopTools_ListOfShape& NE   ) const
{
  TopoDS_Edge WE = E;
  WE.Orientation(TopAbs_FORWARD);

  Standard_Real    U1,U2, f, l;
  TopoDS_Vertex    V1,V2,VF,VL;

  BRep_Tool::Range(WE,f,l);
  TopExp::Vertices(WE,VF,VL);

  if (VOnE.Extent() < 3) { // do not rebuild not cut edge
    if (( VF.IsSame( VOnE.First() ) && VL.IsSame( VOnE.Last() )) ||
	VL.IsSame( VOnE.First() ) && VF.IsSame( VOnE.Last() )  ) {
      NE.Append( E );
      return;
    }
  }

  TopTools_SequenceOfShape SV;
  TopTools_ListIteratorOfListOfShape itv(VOnE);
  TopTools_MapOfOrientedShape VM( VOnE.Extent() );
  for (; itv.More(); itv.Next())
    if ( VM.Add( itv.Value() ))
      SV.Append(itv.Value());

  Tri( WE, SV, myInter3d );

  if (SV.Length() < 3) { // do not rebuild not cut edge
    if (( VF.IsSame( SV.First() ) && VL.IsSame( SV.Last() )) ||
	VL.IsSame( SV.First() ) && VF.IsSame( SV.Last() )  ) {
      NE.Append( E );
      return;
    }
  }

  Standard_Integer iVer, NbVer = SV.Length();


  //----------------------------------------------------------------
  // Construction of the new edges .
  //----------------------------------------------------------------

  if (VF.IsSame(VL)) { // closed edge
    if (NbVer==1) 
      SV.Append( SV.First() );
    else if (!SV.First().IsSame(SV.Last())) {
      Standard_Boolean isFirst=0;
      Standard_Real    minDU = 1.e10;
      TopoDS_Vertex endV = Partition_Inter2d::FindEndVertex(VOnE, f,l, E, isFirst,minDU);
      if (endV.IsSame(SV.First()))
	SV.Append(endV);
      else if (endV.IsSame(SV.Last()))
	SV.Prepend(endV);
      else
	MESSAGE ("END VERTEX IS IN SEQUNCE MIDDLE");
    }
    NbVer = SV.Length();
  }

  for (iVer=1; iVer < NbVer; iVer++) {
    V1  = TopoDS::Vertex(SV(iVer));
    V2  = TopoDS::Vertex(SV(iVer+1));
    
    TopoDS_Shape NewEdge = WE.EmptyCopied();
    V1.Orientation(TopAbs_FORWARD);
    myBuilder.Add  (NewEdge,V1);
    V2.Orientation(TopAbs_REVERSED);
    myBuilder.Add  (NewEdge,V2);
    
    if (iVer==1)
      U1 = f;
    else 	{
      V1.Orientation(TopAbs_INTERNAL);
      U1=BRep_Tool::Parameter(V1,WE);
    }
    if (iVer+1 == NbVer)
      U2 = l;
    else	{
      V2.Orientation(TopAbs_INTERNAL);
      U2=BRep_Tool::Parameter(V2,WE);
    }
    if (Abs(U1-U2) <= Precision::PConfusion()) {
      MESSAGE( "MakeEdges(), EQUAL PARAMETERS OF DIFFERENT VERTICES");
      continue;
    }
    TopoDS_Edge EE=TopoDS::Edge(NewEdge);
    myBuilder.Range (EE,U1,U2);

    TopoDS_Edge NEdge = TopoDS::Edge(NewEdge);
    myBuilder.SameParameter(NEdge,Standard_False);

    Standard_Real tol = 1.0e-2;
    Standard_Boolean flag = BRep_Tool::SameParameter(NEdge);
    if (!flag) {
      BRepLib::SameParameter(NEdge,tol);
    }
    NE.Append(NEdge.Oriented(E.Orientation()));
  }
}

//=======================================================================
//function : MergeEqualEdges
//purpose  : find equal edges,  choose  ones  to  keep and make
//           them have pcurves on all faces they are shared by
//=======================================================================

void Partition_Spliter::MergeEqualEdges (const TopTools_ListOfShape& LSE)
{
  // find equal edges
  // map: edge - equal edges
  TopTools_DataMapOfShapeListOfShape EEM( LSE.Extent() );
  findEqual (LSE, LSE, 0, EEM, myEqualEdges);

  TopTools_ListOfShape EEL; // list of equal edges
  TopTools_DataMapIteratorOfDataMapOfShapeListOfShape itM (EEM);
  for ( ; itM.More(); itM.Next()) {
    EEL = itM.Value();
    EEL.Append( itM.Key() );

    // choose an edge to keep, section edges have priority
    TopoDS_Edge EKeep;
    TopTools_ListIteratorOfListOfShape itEE (EEL);
    for (; itEE.More(); itEE.Next()) {
      EKeep = TopoDS::Edge( itEE.Value() );
      const TopoDS_Edge& EKeepOrig = TopoDS::Edge( myImagesEdges.Root( EKeep ));
      if (myInter3d.IsSectionEdge( EKeepOrig ))
        break;
    }

    // update edge images and build pcurves
    Standard_Real f,l, tol;
    for (itEE.Initialize (EEL); itEE.More(); itEE.Next()) {
      const TopoDS_Edge& E = TopoDS::Edge( itEE.Value() );
      if ( E.IsSame( EKeep )) 
        continue;

      // 1. build pcurves of the kept edge on faces where replaced edges exist
      const TopoDS_Edge& EReplOrig = TopoDS::Edge( myImagesEdges.Root( E ));
      TopTools_ListOfShape FL;
      FL = myAsDes->Ascendant( EReplOrig );
      Standard_Integer iFace, iFirstSectionFace = FL.Extent() + 1;
      // add faces where the replaced edge is a section edge
      if (myInter3d.IsSectionEdge( EReplOrig )) {
        TopTools_ListIteratorOfListOfShape seIt;
        seIt.Initialize( myInter3d.SectionEdgeFaces ( EReplOrig ));
        for ( ; seIt.More(); seIt.Next())
          FL.Append( seIt.Value() );
      }
      // loop on faces
      TopTools_ListIteratorOfListOfShape itF (FL);
      for ( iFace = 1 ; itF.More(); itF.Next(), ++iFace ) {
        const TopoDS_Face& F = TopoDS::Face( itF.Value());

        Handle(Geom2d_Curve) pc = BRep_Tool::CurveOnSurface( EKeep, F, f,l);
        if (pc.IsNull()) {
          Handle(Geom_Curve) C3d = BRep_Tool::Curve( EKeep, f, l);
          C3d = new Geom_TrimmedCurve( C3d, f,l);
          pc = TopOpeBRepTool_CurveTool::MakePCurveOnFace (F,C3d,tol);
          if (pc.IsNull()) {
            MESSAGE (" CANT BUILD PCURVE ");
          }
          myBuilder.UpdateEdge( EKeep, pc, F, tol);
        }

        if (iFace >= iFirstSectionFace ||
            !BRep_Tool::IsClosed( EReplOrig, F ))
          continue;

        // build the second pcurve for a seam
        TopoDS_Vertex V = TopExp::FirstVertex( EKeep );
        Standard_Real Ukeep = BRep_Tool::Parameter( V, EKeep );
        Standard_Real Urepl = BRep_Tool::Parameter( V, E );

        TopoDS_Edge EReplRev = E;
        EReplRev.Reverse();
        Handle(Geom2d_Curve) pcRepl1 = BRep_Tool::CurveOnSurface( E, F, f,l);
        Handle(Geom2d_Curve) pcRepl2 = BRep_Tool::CurveOnSurface( EReplRev, F, f,l);

        gp_Pnt2d p1r, p2r, pk;
        p1r = pcRepl1->Value( Urepl );
        p2r = pcRepl2->Value( Urepl );
        pk  = pc->Value( Ukeep );

        // suppose that pk is equal to either p1r or p2r
        Standard_Boolean isUPeriod =
          ( Abs( p1r.X() - p2r.X() ) > Abs( p1r.Y() - p2r.Y() ));
        Standard_Boolean is1Equal;
        if (isUPeriod)
          is1Equal = ( Abs( p1r.X() - pk.X() ) < Abs( p2r.X() - pk.X() ));
        else
          is1Equal = ( Abs( p1r.Y() - pk.Y() ) < Abs( p2r.Y() - pk.Y() ));

        Handle(Geom2d_Curve) pc2 = Handle(Geom2d_Curve)::DownCast
          ( pc->Translated( pk, is1Equal ? p2r : p1r ) );

        if (E.Orientation() == TopAbs_REVERSED)
          is1Equal = !is1Equal;

        if (is1Equal)
          myBuilder.UpdateEdge( EKeep, pc, pc2, F, tol);
        else
          myBuilder.UpdateEdge( EKeep, pc2, pc, F, tol);

      } // loop on a Faces where a replaced edge exists


      // 2. update edge images according to replacement
      if (myImagesEdges.HasImage( E ))
        myImagesEdges.Remove( E );
      myImagesEdges.Bind( E, EKeep );

    } // loop on a list of equal edges EEL
  } // loop on a map of equal edges EEM
}

//=======================================================================
//function : KeepShapesInside
//purpose  : remove shapes that are outside of S from resul
//=======================================================================

void Partition_Spliter::KeepShapesInside (const TopoDS_Shape& S)
{
  TopoDS_Iterator it;
  if (S.ShapeType() < TopAbs_SOLID) { // compound or compsolid
    for (it.Initialize( S ); it.More(); it.Next())
      KeepShapesInside( it.Value());
    return;
  }

  Standard_Boolean isTool = Standard_False;
  if (!myImageShape.HasImage( S )) {
    isTool = CheckTool( S );
    if (!isTool) return;
  }

  // build map of internal faces
  TopTools_IndexedMapOfShape MIF;
  TopoDS_Shape IntFacesComp = FindFacesInside( S, Standard_False, Standard_True);
  TopExp::MapShapes( IntFacesComp, TopAbs_FACE, MIF );

  TopoDS_Compound C;
  myBuilder.MakeCompound(C);

  TopAbs_ShapeEnum anInternalShapeType = TopAbs_SHAPE;
  if (!MIF.IsEmpty())
  {
    // leave in the result only those shapes having a face in MIF
    for (it.Initialize( myShape ); it.More(); it.Next()) {
      const TopoDS_Shape & aResShape = it.Value();
      TopExp_Explorer expResF( aResShape, TopAbs_FACE );
      for (; expResF.More(); expResF.Next()) {
        if ( MIF.Contains( expResF.Current())) {
          myBuilder.Add( C, aResShape );
          if (aResShape.ShapeType() < anInternalShapeType)
            anInternalShapeType = aResShape.ShapeType();
          break;
        }
      }
    }
  }

  // may be S was not split by internal faces then it is missing
  // in myShape, add it
  if (!isTool &&
      (anInternalShapeType > TopAbs_SOLID || S.ShapeType() > TopAbs_SOLID))
  {
    TopTools_IndexedMapOfShape MSF; // map of split faces of S
    TopExp::MapShapes( myImageShape.Image(S).First(), TopAbs_FACE, MSF);

    // find a shape having all faces in MSF
    for (it.Initialize( myShape ); it.More(); it.Next()) {
      TopExp_Explorer expResF( it.Value(), TopAbs_FACE );
      for (; expResF.More(); expResF.Next()) {
        if (! MSF.Contains( expResF.Current())) 
          break;
      }
      if (! expResF.More()) {
        myBuilder.Add( C, it.Value() );
        break;
      }
    }
  }

  myShape = C;
}

//=======================================================================
//function : RemoveShapesInside
//purpose  : remove shapes that are inside S from resul
//=======================================================================

void Partition_Spliter::RemoveShapesInside (const TopoDS_Shape& S)
{
  TopoDS_Iterator it;
  if (S.ShapeType() < TopAbs_SOLID) { // compound or compsolid
    for (it.Initialize( S ); it.More(); it.Next())
      RemoveShapesInside( it.Value());
    return;
  }
  Standard_Boolean isTool = Standard_False;
  if (!myImageShape.HasImage( S )) {
    isTool = CheckTool( S );
    if (!isTool) return;
  }

  TopoDS_Shape IntFacesComp = FindFacesInside( S, Standard_False, Standard_True);
  TopTools_IndexedMapOfShape MIF; // map of internal faces
  TopExp::MapShapes( IntFacesComp, TopAbs_FACE, MIF);

  if (MIF.IsEmpty()) return;

  // add to MIF split faces of S
  if (myImageShape.HasImage(S))
    TopExp::MapShapes( myImageShape.Image(S).First(), TopAbs_FACE, MIF);

  // leave in the result only those shapes not having all face in MIF
  
  TopoDS_Compound C;
  myBuilder.MakeCompound(C);

  // RMF : faces of removed shapes that encounter once
  TopTools_MapOfShape RFM;
  
  for (it.Initialize( myShape ); it.More(); it.Next()) {
    
    TopExp_Explorer expResF( it.Value(), TopAbs_FACE );
    for (; expResF.More(); expResF.Next())
      if (!MIF.Contains( expResF.Current()))
	break;

    if (expResF.More())
      // add shape to result
      myBuilder.Add( C, it.Value() );
    else 
      // add faces of a removed shape to RFM
      for (expResF.ReInit(); expResF.More(); expResF.Next()) {
	const TopoDS_Shape& F = expResF.Current();
	if ( ! RFM.Remove ( F ))
	  RFM.Add( F );
      }
  }

  if (!isTool) {

    // rebuild S, it must remain in the result

    Standard_Boolean isClosed = Standard_False;
    switch (S.ShapeType()) {
    case TopAbs_SOLID :
      isClosed = Standard_True; break;
    case TopAbs_SHELL: {
      TopTools_IndexedDataMapOfShapeListOfShape MEF;
      TopExp::MapShapesAndAncestors(S, TopAbs_EDGE, TopAbs_FACE, MEF);
      Standard_Integer i;
      for (i=1;  isClosed && i<=MEF.Extent();  ++i) 
        isClosed = ( MEF(i).Extent() != 1 );
      break;
    }
    default:
      isClosed = Standard_False;
    }
    if (isClosed) {

      // add to a new shape external faces of removed shapes, ie those in RFM

      TopoDS_Shell Shell;
      myBuilder.MakeShell( Shell );

      // exclude redundant internal face with edges encounterd only once
      TopTools_IndexedDataMapOfShapeListOfShape MEF;
      TopTools_MapIteratorOfMapOfShape itF (RFM);
      for ( ; itF.More(); itF.Next()) 
        TopExp::MapShapesAndAncestors(itF.Key(), TopAbs_EDGE, TopAbs_FACE, MEF);

      // add only faces forming a closed shell
      for (itF.Reset() ; itF.More(); itF.Next())
      {
        TopExp_Explorer expE (itF.Key(), TopAbs_EDGE);
        for (; expE.More(); expE.Next())
          if (MEF.FindFromKey(expE.Current()).Extent() == 1)
            break;
        if (!expE.More())
          myBuilder.Add( Shell, itF.Key());
      }

      if (S.ShapeType() == TopAbs_SOLID) {
        TopoDS_Solid Solid;
        myBuilder.MakeSolid( Solid );
        myBuilder.Add (Solid, Shell);
        myBuilder.Add (C, Solid);
      }
      else
        myBuilder.Add (C, Shell);
    }
    else {
      if (myImageShape.HasImage( S )) {
        for (it.Initialize( myImageShape.Image(S).First()); it.More(); it.Next())
          myBuilder.Add (C, it.Value());
      }
    }
  }
  
  myShape = C;
}

//=======================================================================
//function : CheckTool
//purpose  : Return True if <S>  is  a tool shape. Prepare tool
//           faces of <S> for the search of internal faces.
//=======================================================================

Standard_Boolean Partition_Spliter::CheckTool(const TopoDS_Shape& S)
{
  // suppose S has not an image
  
  Standard_Boolean isTool = Standard_False;
  TopoDS_Compound C;
  myBuilder.MakeCompound( C );

  TopExp_Explorer expF( S, TopAbs_FACE);
  for (; expF.More(); expF.Next()) {

    const TopoDS_Face& F = TopoDS::Face( expF.Current() );
    if (myMapTools.Contains( F ))
      isTool = Standard_True;
    else
      continue;

    if (myImagesFaces.HasImage( F )) {
      // F has been reconstructed
      TopAbs_Orientation Fori = F.Orientation();
      TopTools_ListOfShape LNF;
      myImagesFaces.LastImage( F, LNF);
      TopTools_ListIteratorOfListOfShape itF (LNF);
      for ( ; itF.More(); itF.Next())
	myBuilder.Add( C, itF.Value().Oriented(Fori) );
      continue;
    }
    
    Standard_Boolean hasSectionE = myInter3d.HasSectionEdge( F );
    Standard_Boolean hasNewE     = myAsDes->HasDescendant( F );
    if (!hasSectionE && !hasNewE)
    {
      // F intersects nothing
      myBuilder.Add( C, F );
      continue;
    }
    
    // make an image for F
    
    TopoDS_Face NF = F;
    NF.Orientation(TopAbs_FORWARD);
    NF = TopoDS::Face( NF.EmptyCopied() ); // make a copy
    TopoDS_Wire NW;
    myBuilder.MakeWire( NW );

    // add edges, as less as possible
    TopTools_ListOfShape NEL;
    TopTools_ListIteratorOfListOfShape itNE;
    if (hasSectionE) {
      // add section edges
      TopExp_Explorer expE;
      for ( ; expE.More(); expE.Next()) {
	if (! myImagesEdges.HasImage( expE.Current() ))
	  continue;
	myImagesEdges.LastImage( expE.Current(), NEL );
	for ( itNE.Initialize( NEL ); itNE.More(); itNE.Next())
	  myBuilder.Add ( NW, itNE.Value());
      }
    }
    if (hasNewE) {
      // add new adges
      NEL = myAsDes->Descendant( F );
      for ( itNE.Initialize( NEL ); itNE.More(); itNE.Next()) {
	TopTools_ListOfShape SEL; // splits
	myImagesEdges.LastImage( itNE.Value(), SEL );
	TopTools_ListIteratorOfListOfShape itSE (SEL);
	for ( ; itSE.More(); itSE.Next()) 
	  myBuilder.Add ( NW, itSE.Value());
      }
    }
    myBuilder.Add( NF, NW );
    myBuilder.Add (C, NF);
    
    NF.Orientation( F.Orientation() ); // NF is most probably invalid
    myImagesFaces.Bind (F, NF);
  }
  if (isTool)
    myImageShape.Bind (S, C);

  return isTool;
}

//=======================================================================
//function : IsInside
//purpose  : Return True if the first vertex of S1 inside S2.
//           If S1.IsNull(), check infinite point against S2.
//=======================================================================

Standard_Boolean Partition_Spliter::IsInside (const TopoDS_Shape& theS1,
                                              const TopoDS_Shape& theS2)
{
  BRepClass3d_SolidClassifier aClassifier( theS2 );

  TopExp_Explorer expl( theS1, TopAbs_VERTEX );
  if (!expl.More())
    aClassifier.PerformInfinitePoint( ::RealSmall());
  else
  {
    const TopoDS_Vertex & aVertex = TopoDS::Vertex( expl.Current() );
    aClassifier.Perform (BRep_Tool::Pnt( aVertex ),
                         BRep_Tool::Tolerance( aVertex ));
  }

  return ( aClassifier.State() == TopAbs_IN );
}

//=======================================================================
//function : GetOriginalShape
//purpose  : Return the  shape  aShape  originates from. aShape
//           should be a face or more complex result shape
//=======================================================================

TopoDS_Shape Partition_Spliter::GetOriginalShape(const TopoDS_Shape& theShape) const
{
  TopoDS_Shape anOrigShape;

  TopExp_Explorer expl( theShape, TopAbs_FACE);
  if (expl.More())
  {

    TopoDS_Shape aFace = expl.Current();
    if (myImagesFaces.IsImage( aFace ))
      aFace = myImagesFaces.Root( aFace );
    anOrigShape = myFaceShapeMap.Find( aFace );
  }
  return anOrigShape;
}

//=======================================================================
//function : FindToolsToReconstruct
//purpose  : find and store  as  objects  tools which interfere
//           with  solids   or   are   inside   solids  without
//           an interference
//=======================================================================

void Partition_Spliter::FindToolsToReconstruct()
{
  if (myMapTools.IsEmpty())
    return;

  Standard_Integer nbFoundTools = 0;

  // build edge - face map in order to detect interference with section edges
  TopTools_IndexedDataMapOfShapeListOfShape EFM;
  TopTools_MapIteratorOfMapOfShape aMapIt;
  for (aMapIt.Initialize(myMapTools); aMapIt.More(); aMapIt.Next())
    TopExp::MapShapesAndAncestors( aMapIt.Key(), TopAbs_EDGE, TopAbs_FACE, EFM);
  for (aMapIt.Initialize(myMapFaces); aMapIt.More(); aMapIt.Next())
    TopExp::MapShapesAndAncestors( aMapIt.Key(), TopAbs_EDGE, TopAbs_FACE, EFM);

  TopTools_MapOfShape aCurrentSolids, aCheckedShapes;

  // faces cut by new edges
  TopTools_MapOfShape & aSectionFaces = myInter3d.TouchedFaces();

  // keep solids interfering with each other in aCurrentSolids map
  // and add tool faces intersecting solids as object shapes

  TopTools_ListIteratorOfListOfShape itS, itF, itCF, itE;
  for (itS.Initialize( myListShapes ); itS.More(); itS.Next()) {
    TopExp_Explorer expSo (itS.Value(), TopAbs_SOLID);
    for (; expSo.More(); expSo.Next()) {

      // check if a solid has been already processed
      const TopoDS_Shape & aSo = expSo.Current();
      if (!aCheckedShapes.Add( aSo ))
        continue;
      aCurrentSolids.Add( aSo );

      // faces to check
      TopTools_ListOfShape aFacesToCheck;
      TopExp_Explorer exp( aSo, TopAbs_FACE );
      for ( ; exp.More(); exp.Next())
        aFacesToCheck.Append ( exp.Current());

      // add other shapes interefering with a solid.
      // iterate faces to check while appending new ones
      for (itCF.Initialize (aFacesToCheck) ; itCF.More(); itCF.Next())
      {
        const TopoDS_Shape& aCheckFace = itCF.Value();
//         if (!aCheckedShapes.Add( aCheckFace ))
//           continue;

        // find faces interfering with aCheckFace
        TopTools_ListOfShape anIntFaces;

        // ** 1. faces intersecting aCheckFace with creation of new edges on it
        if ( myAsDes->HasDescendant( aCheckFace ))
        {
          // new edges on aCheckFace
          const TopTools_ListOfShape& NEL = myAsDes->Descendant( aCheckFace );
          for (itE.Initialize( NEL); itE.More(); itE.Next())
          {
            const TopoDS_Shape & aNewEdge = itE.Value();
            if (!aCheckedShapes.Add( aNewEdge ))
              continue;

            // faces interfering by aNewEdge
            itF.Initialize (myAsDes->Ascendant( aNewEdge ));
            for (; itF.More(); itF.Next())
              if (aCheckFace != itF.Value())
                anIntFaces.Append( itF.Value() );

            // ** 2. faces having section edge aNewEdge on aFacesToCheck
            if (EFM.Contains( aNewEdge))
            {
              itF.Initialize ( EFM.FindFromKey (itE.Value()));
              for (; itF.More(); itF.Next())
                if (aCheckFace != itF.Value())
                  anIntFaces.Append( itF.Value() );
            }
          }
        }

        // ** 3. faces cut by edges of aCheckFace
        TopExp_Explorer expE (aCheckFace, TopAbs_EDGE);
        for ( ; expE.More(); expE.Next())
        {
          const TopoDS_Shape & aCheckEdge = expE.Current();
          if (aCheckedShapes.Add( aCheckEdge ) &&
              myInter3d.IsSectionEdge( TopoDS::Edge( aCheckEdge )))
          {
            itF.Initialize( myInter3d.SectionEdgeFaces( TopoDS::Edge( aCheckEdge )));
            for (; itF.More(); itF.Next()) 
              if (aCheckFace != itF.Value())
                anIntFaces.Append( itF.Value() );
          }
        }

        // process faces interfering with aCheckFace and shapes they
        // belong to
        for (itF.Initialize (anIntFaces); itF.More(); itF.Next())
        {
          const TopoDS_Shape & F = itF.Value();
          if (! aCheckedShapes.Add( F ))
            continue;

          Standard_Boolean isTool = myMapTools.Contains( F );
          if (isTool && 
              myFaceShapeMap( aCheckFace ).ShapeType() == TopAbs_SOLID )
          {
            // a tool interfering with a solid
            if (aSectionFaces.Contains( F ))
              AddShape( F );
            ++ nbFoundTools;
            if (nbFoundTools == myMapTools.Extent())
              return;
          }

          const TopoDS_Shape & S = myFaceShapeMap( F );
          if (aCheckedShapes.Add( S ))
          {
            // a new shape interefering with aCurrentSolids is found
            if (!isTool && S.ShapeType() == TopAbs_SOLID)
              aCurrentSolids.Add ( S );
            // add faces to aFacesToCheck list
            for ( exp.Init( S, TopAbs_FACE ); exp.More(); exp.Next())
              aFacesToCheck.Append ( exp.Current() );
          }
        }
      } // loop on aFacesToCheck

      // Here aCurrentSolids contains all solids interfering with each other.
      // aCheckedShapes contains all faces belonging to shapes included
      // in or interfering with aCurrentSolids or previously checked solids.
      // Test if tool faces that do not interefere with other shapes are
      // wrapped by any of aCurrentSolids

      TopTools_MapIteratorOfMapOfShape aSolidIt (aCurrentSolids);
      for ( ; aSolidIt.More(); aSolidIt.Next())
      {
        const TopoDS_Shape & aSolid = aSolidIt.Key();
        TopTools_MapOfShape aCheckedTools( myMapTools.Extent() );

        TopTools_MapIteratorOfMapOfShape aToolIt (myMapTools);
        for ( ; aToolIt.More(); aToolIt.Next())
        {
          const TopoDS_Shape & aToolFace = aToolIt.Key();
          if (aCheckedShapes.Contains( aToolFace ) || // already found
              aCheckedTools.Contains( aToolFace ))    // checked against aSolid
            continue;

          const TopoDS_Shape & aToolShape = myFaceShapeMap( aToolFace );
          TopExp_Explorer aToolFaceIt( aToolShape, TopAbs_FACE );
          
          Standard_Boolean isInside = IsInside( aToolShape, aSolid );
          for ( ; aToolFaceIt.More(); aToolFaceIt.Next() )
          {
            const TopoDS_Shape & aTool = aToolFaceIt.Current();
            aCheckedTools.Add( aTool );
            if (isInside)
            {
              if (aSectionFaces.Contains( aTool ))
                AddShape( aTool );
              ++ nbFoundTools;
              if (nbFoundTools == myMapTools.Extent())
                return;
              aCheckedShapes.Add( aTool );
            }
          }
        }
      }
      
    } // loop on solid shapes
  }
}

#endif
