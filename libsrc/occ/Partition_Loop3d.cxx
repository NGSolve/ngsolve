#ifdef OCCGEOMETRY

//  GEOM PARTITION : partition algorithm
//
//  Copyright (C) 2003  CEA/DEN, EDF R&D
//
//
//
//  File   : Partition_Loop3d.cxx
//  Module : GEOM

//using namespace std;
#include <climits>
#include "Partition_Loop3d.ixx"

#include <TopExp_Explorer.hxx>
#include <TopExp.hxx>
#include <BRep_Builder.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS.hxx>
#include <TopTools_MapIteratorOfMapOfShape.hxx>
#include <gp_Vec.hxx>
#include <gp_Pnt.hxx>
#include <Geom2d_Curve.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec2d.hxx>
#include <gp_Dir2d.hxx>
#include <Geom_Curve.hxx>

//=======================================================================
//function : Partition_Loop3d
//purpose  : 
//=======================================================================

Partition_Loop3d::Partition_Loop3d()
{
}

//=======================================================================
//function : AddConstFaces
//purpose  : Add faces of <S> as unique faces in the result.
//=======================================================================

void Partition_Loop3d::AddConstFaces(const TopoDS_Shape& S) 
{
  TopExp_Explorer FaceExp(S, TopAbs_FACE);
  for (; FaceExp.More(); FaceExp.Next())
    myFaces.Append( FaceExp.Current() );

  TopExp::MapShapesAndAncestors(S, TopAbs_EDGE, TopAbs_FACE, myEFMap);
}

//=======================================================================
//function : AddSectionFaces
//purpose  : Add faces of <S> as double faces in the result.
//=======================================================================

void Partition_Loop3d::AddSectionFaces(const TopoDS_Shape& S) 
{
  AddConstFaces( S );
  AddConstFaces( S.Reversed() );
}

//=======================================================================
//function : MakeShells
//purpose  : Make and return shells. 
//           <AvoidFacesMap> can contain faces that must not be
//           added to result shells.
//=======================================================================

const TopTools_ListOfShape&
  Partition_Loop3d::MakeShells (const TopTools_MapOfOrientedShape& AvoidFacesMap)
{
  myNewShells.Clear();
  
  BRep_Builder Builder;
  TopTools_MapOfShape CheckedEdgesMap;
  TopTools_MapOfOrientedShape AddedFacesMap;
  
  TopTools_ListIteratorOfListOfShape itF (myFaces);
  for (; itF.More(); itF.Next())
  {
    const TopoDS_Shape& FF = itF.Value();
    if (AvoidFacesMap.Contains( FF ) ||
	! AddedFacesMap.Add( FF ) )
      continue;

    // make a new shell
    TopoDS_Shell Shell;
    Builder.MakeShell(Shell);
    Builder.Add(Shell,FF);

    // clear the maps from shapes added to previous Shell
    TopTools_MapIteratorOfMapOfShape itEM (CheckedEdgesMap);
    for (; itEM.More(); itEM.Next()) {
      TopTools_ListOfShape& FL = myEFMap.ChangeFromKey( itEM.Key());
      TopTools_ListIteratorOfListOfShape it (FL);
      while ( it.More()) {
        if (AddedFacesMap.Contains( it.Value()))
          FL.Remove( it );
        else
          it.Next();
      }
    }
    CheckedEdgesMap.Clear();

    
    // loop on faces added to Shell; add their neighbor faces to Shell and so on
    TopoDS_Iterator itAddedF (Shell);
    for (; itAddedF.More(); itAddedF.Next())
    {
      const TopoDS_Face& F = TopoDS::Face (itAddedF.Value());

      // loop on edges of F; find a good neighbor face of F by E
      TopExp_Explorer EdgeExp(F, TopAbs_EDGE);
      for (; EdgeExp.More(); EdgeExp.Next())
      {
        const TopoDS_Edge& E = TopoDS::Edge( EdgeExp.Current());
	if (! CheckedEdgesMap.Add( E ))
	  continue;

	// candidate faces list
        const TopTools_ListOfShape& FL = myEFMap.ChangeFromKey(E);
        if (FL.IsEmpty())
          continue;
	// select one of neighbors
        TopoDS_Face SelF;
        if (FL.Extent() == 2) {
          if (! F.IsSame( FL.First() ))
            SelF = TopoDS::Face( FL.First() );
          else if (!F.IsSame( FL.Last() ))
            SelF = TopoDS::Face( FL.Last() );
        }
        else {
          // check if a face already added to Shell shares E
	  TopTools_ListIteratorOfListOfShape it (FL);
          Standard_Boolean found = Standard_False;
          for (; !found && it.More(); it.Next())
            if (F != it.Value())
              found = AddedFacesMap.Contains( it.Value() );
          if (found)
            continue;
          // select basing on geometrical check
          Standard_Boolean GoodOri, inside;
          Standard_Real dot, MaxDot = -100;
          TopTools_ListOfShape TangFL; // tangent faces
          for ( it.Initialize( FL ) ; it.More(); it.Next()) {
            const TopoDS_Face& NeighborF = TopoDS::Face( it.Value());
            if (NeighborF.IsSame( F ))
              continue;
            inside = Partition_Loop3d::IsInside( E, F, NeighborF, 1, dot, GoodOri);
            if (!GoodOri)
              continue;
            if (!inside)
              dot = -dot - 3;
            if (dot < MaxDot)
              continue;
            if ( IsEqual( dot, MaxDot))
              TangFL.Append(SelF);
            else
              TangFL.Clear();
            MaxDot = dot;
            SelF = NeighborF;
          }
          if (!TangFL.IsEmpty()) {
            for (it.Initialize( TangFL ); it.More(); it.Next()) {
              const TopoDS_Face& NeighborF = TopoDS::Face( it.Value());
              if (Partition_Loop3d:: IsInside( E, SelF , NeighborF, 0, dot, GoodOri))
                SelF = NeighborF;
            }
          }
        }
        if (!SelF.IsNull() &&
	    AddedFacesMap.Add( SelF ) &&
	    !AvoidFacesMap.Contains( SelF )) 
          Builder.Add( Shell, SelF);

      } // loop on edges of F
      
    } // loop on the faces added to Shell

    // Shell is complete
    myNewShells.Append( Shell );

  } // loop on myFaces


  // prepare to the next call
  myFaces.Clear();
  myEFMap.Clear();

  return myNewShells;
}



//=======================================================================
//function : Normal
//purpose  : 
//=======================================================================

gp_Vec Partition_Loop3d::Normal(const TopoDS_Edge& E,
				const TopoDS_Face& F)
{
  gp_Vec Norm, V1, V2;
  Standard_Real First, Last;
  gp_Pnt Ps;

  Handle(Geom2d_Curve) C2d = BRep_Tool::CurveOnSurface (E, F, First, Last);
  Handle(Geom_Surface) Sf = BRep_Tool::Surface(F);

  gp_Pnt2d p = C2d->Value( 0.5*(First+Last) );
  Sf->D1(p.X(), p.Y(), Ps, V1, V2);
  Norm = V1.Crossed(V2);

  if (F.Orientation() == TopAbs_REVERSED ) 
    Norm.Reverse();

  return Norm;
}

//=======================================================================
//function : NextNormal
//purpose  : find normal to F at point a little inside F near the middle of E
//warning  : E must be properly oriented in F.
//=======================================================================

static gp_Vec NextNormal(const TopoDS_Edge& E,
			 const TopoDS_Face& F)
{
  Standard_Real First, Last;

  Handle(Geom2d_Curve) C2d = BRep_Tool::CurveOnSurface (E, F, First, Last);
  Handle(Geom_Surface) Sf = BRep_Tool::Surface(F);

  gp_Pnt2d p;
  gp_Vec2d v;
  C2d->D1( 0.5*(First+Last), p, v);
  if (E.Orientation() != F.Orientation())
    v.Reverse();
  gp_Dir2d dir( -v.Y(), v.X() ); // dir inside F
  
  Standard_Real duv = 1e-6; // this is not Ok and may give incorrect result if
  // resolutionUV of compared faces is very different. To have a good result,
  //it is necessary to get normal to faces at points equidistant from E in 3D
  
  p.SetX( p.X() + dir.X()*duv );
  p.SetY( p.Y() + dir.Y()*duv );
  
  gp_Pnt Ps;
  gp_Vec Norm, V1, V2, VV1, VV2;
  Sf->D1( p.X(), p.Y(), Ps, V1, V2);
  Norm = V1.Crossed(V2);

  if (F.Orientation() == TopAbs_REVERSED ) 
    Norm.Reverse();

  return Norm;
}


//=======================================================================
//function : FindEinF
//purpose  : find E in F
//=======================================================================

static TopoDS_Edge FindEinF(const TopoDS_Edge& E,
			    const TopoDS_Face& F)
{
  TopExp_Explorer expl (F, TopAbs_EDGE);
  for (; expl.More(); expl.Next()) 
    if( E.IsSame( expl.Current() ))
      return TopoDS::Edge(expl.Current());
  TopoDS_Edge nullE;
  return nullE;
}

//=======================================================================
//function : IsInside
//purpose  : check if <F2> is inside <F1> by edge <E>.
//           if <CountDot>, compute <Dot>: scalar production of
//           normalized  vectors  pointing  inside  faces,  and
//           check if faces are oriented well for sewing
//=======================================================================

Standard_Boolean Partition_Loop3d::IsInside(const TopoDS_Edge& E,
					    const TopoDS_Face& F1,
					    const TopoDS_Face& F2,
					    const Standard_Boolean CountDot,
					    Standard_Real& Dot,
					    Standard_Boolean& GoodOri) 
{
  Standard_Real f, l;
  gp_Pnt P;
  gp_Vec Vc1, Vc2, Vin1, Vin2, Nf1, Nf2;
  Handle(Geom_Curve) Curve = BRep_Tool::Curve(E,f,l);
  Curve->D1( 0.5*(f + l), P, Vc2);
  TopoDS_Edge E1, E2 = FindEinF (E, F2);
  if (E2.Orientation() == TopAbs_REVERSED ) Vc2.Reverse();

  Nf1 = Normal(E,F1);
  Nf2 = Normal(E,F2);

  Standard_Real sin =
    Nf1.CrossSquareMagnitude(Nf2) / Nf1.SquareMagnitude() / Nf2.SquareMagnitude();
  Standard_Boolean tangent = sin < 0.001;

  Standard_Boolean inside = 0;
  if (tangent) {
    E1 = FindEinF (E, F1);
    gp_Vec NNf1 = NextNormal(E1,F1);
    gp_Vec NNf2 = NextNormal(E2,F2);
    Vin2 = NNf2.Crossed(Vc2);
    inside = Vin2 * NNf1 < 0;
  }
  else {
    Vin2 = Nf2.Crossed(Vc2);
    inside = Vin2 * Nf1 < 0;
  }
  
  if (!CountDot) return inside;

  if (tangent)
    Vin2 = Nf2.Crossed(Vc2);
  else
    E1 = FindEinF (E, F1);
    
  Vc1 = Vc2;
  if (E1.Orientation() != E2.Orientation()) 
    Vc1.Reverse();
  Vin1 = Nf1.Crossed(Vc1);

  if (tangent) {
    Standard_Real N1N2 = Nf1 * Nf2;
    GoodOri = (Vin2 * Vin1 < 0) ? N1N2 > 0 : N1N2 < 0;
  }
  else {
    Standard_Real V1N2 = Vin1 * Nf2;
    GoodOri = ( inside ? V1N2 <= 0 : V1N2 >= 0);
  }

  Vin1.Normalize();
  Vin2.Normalize();
  
  Dot = Vin2 * Vin1;
  
  return inside;
}


#endif
