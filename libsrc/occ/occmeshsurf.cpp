#ifdef OCCGEOMETRY

#include <mystdlib.h>

#include <occgeom.hpp>
#include <meshing.hpp>
#include <GeomLProp_SLProps.hxx>
#include <ShapeAnalysis_Surface.hxx>


namespace netgen
{
#include "occmeshsurf.hpp"


  bool glob_testout(false);

  void OCCSurface :: GetNormalVector (const Point<3> & p, 
				      const PointGeomInfo & geominfo,
				      Vec<3> & n) const
  {
    gp_Pnt pnt;
    gp_Vec du, dv;

    /*
      double gu = geominfo.u;
      double gv = geominfo.v;

      if (fabs (gu) < 1e-3) gu = 0;
      if (fabs (gv) < 1e-3) gv = 0;

      occface->D1(gu,gv,pnt,du,dv);
    */

    /*
      occface->D1(geominfo.u,geominfo.v,pnt,du,dv);

      n = Cross (Vec<3>(du.X(), du.Y(), du.Z()),
      Vec<3>(dv.X(), dv.Y(), dv.Z()));
      n.Normalize();
    */



    GeomLProp_SLProps lprop(occface,geominfo.u,geominfo.v,1,1e-5);
    double setu=geominfo.u,setv=geominfo.v;

    if(lprop.D1U().Magnitude() < 1e-5 || lprop.D1V().Magnitude() < 1e-5)
      {
	double ustep = 0.01*(umax-umin);
	double vstep = 0.01*(vmax-vmin);

	n=0;

	while(setu < umax && (lprop.D1U().Magnitude() < 1e-5 || lprop.D1V().Magnitude() < 1e-5))
	  setu += ustep;
	if(setu < umax)
	  {
	    lprop.SetParameters(setu,setv);
	    n(0)+=lprop.Normal().X();
	    n(1)+=lprop.Normal().Y();
	    n(2)+=lprop.Normal().Z();
	  }
	setu = geominfo.u;
	while(setu > umin && (lprop.D1U().Magnitude() < 1e-5 || lprop.D1V().Magnitude() < 1e-5))
	  setu -= ustep;
	if(setu > umin)
	  {
	    lprop.SetParameters(setu,setv);
	    n(0)+=lprop.Normal().X();
	    n(1)+=lprop.Normal().Y();
	    n(2)+=lprop.Normal().Z();
	  }
	setu = geominfo.u;

	while(setv < vmax && (lprop.D1U().Magnitude() < 1e-5 || lprop.D1V().Magnitude() < 1e-5))
	  setv += ustep;
	if(setv < vmax)
	  {
	    lprop.SetParameters(setu,setv);
	    n(0)+=lprop.Normal().X();
	    n(1)+=lprop.Normal().Y();
	    n(2)+=lprop.Normal().Z();
	  }
	setv = geominfo.v;
	while(setv > vmin && (lprop.D1U().Magnitude() < 1e-5 || lprop.D1V().Magnitude() < 1e-5))
	  setv -= ustep;
	if(setv > vmin)
	  {
	    lprop.SetParameters(setu,setv);
	    n(0)+=lprop.Normal().X();
	    n(1)+=lprop.Normal().Y();
	    n(2)+=lprop.Normal().Z();
	  }
	setv = geominfo.v;

	n.Normalize();
      }
    else
      {
	n(0)=lprop.Normal().X();
	n(1)=lprop.Normal().Y();
	n(2)=lprop.Normal().Z();
      }

    if(glob_testout)
      {
	(*testout) << "u " << geominfo.u << " v " << geominfo.v 
		   << " du " << lprop.D1U().X() << " "<< lprop.D1U().Y() << " "<< lprop.D1U().Z()
		   << " dv " << lprop.D1V().X() << " "<< lprop.D1V().Y() << " "<< lprop.D1V().Z() << endl;
      }



    if (orient == TopAbs_REVERSED) n = -1*n;
    //  (*testout) << "GetNormalVector" << endl;
  }


  void OCCSurface :: DefineTangentialPlane (const Point<3> & ap1,
					    const PointGeomInfo & geominfo1,
					    const Point<3> & ap2,
					    const PointGeomInfo & geominfo2)
  {
    if (projecttype == PLANESPACE)
      {
	p1 = ap1; p2 = ap2;

	//cout << "p1 = " << p1 << endl;
	//cout << "p2 = " << p2 << endl;
      
	GetNormalVector (p1, geominfo1, ez);
      
	ex = p2 - p1;
	ex -= (ex * ez) * ez;
	ex.Normalize();
	ey = Cross (ez, ex); 

	GetNormalVector (p2, geominfo2, n2);
  
	nmid = 0.5*(n2+ez);
      
	ez = nmid;
	ez.Normalize(); 
      
	ex = (p2 - p1).Normalize();
	ez -= (ez * ex) * ex;
	ez.Normalize();
	ey = Cross (ez, ex);
	nmid = ez;
	//cout << "ex " << ex << " ey " << ey << " ez " << ez << endl;
      }
    else
      {
	if ( (geominfo1.u < umin) ||
	     (geominfo1.u > umax) ||
	     (geominfo2.u < umin) ||
	     (geominfo2.u > umax) ||
	     (geominfo1.v < vmin) ||
	     (geominfo1.v > vmax) ||
	     (geominfo2.v < vmin) ||
	     (geominfo2.v > vmax) ) throw UVBoundsException();
	  

	p1 = ap1; p2 = ap2;
	psp1 = Point<2>(geominfo1.u, geominfo1.v);
	psp2 = Point<2>(geominfo2.u, geominfo2.v);
      
	Vec<3> n;
	GetNormalVector (p1, geominfo1, n);

	gp_Pnt pnt;
	gp_Vec du, dv;
	occface->D1 (geominfo1.u, geominfo1.v, pnt, du, dv);

	DenseMatrix D1(3,2), D1T(2,3), DDTinv(2,2);
	D1(0,0) = du.X(); D1(1,0) = du.Y(); D1(2,0) = du.Z();
	D1(0,1) = dv.X(); D1(1,1) = dv.Y(); D1(2,1) = dv.Z();

	/*
	  (*testout) << "DefineTangentialPlane" << endl
	  << "---------------------" << endl;
	  (*testout) << "D1 = " << endl << D1 << endl;
	*/

	Transpose (D1, D1T);
	DenseMatrix D1TD1(3,3);

	D1TD1 = D1T*D1;
	if (D1TD1.Det() == 0) throw SingularMatrixException();
      
	CalcInverse (D1TD1, DDTinv);
	DenseMatrix Y(3,2);
	Vec<3> y1 = (ap2-ap1).Normalize();
	Vec<3> y2 = Cross(n, y1).Normalize();
	for (int i = 0; i < 3; i++)
	  {
	    Y(i,0) = y1(i);
	    Y(i,1) = y2(i);
	  }

	DenseMatrix A(2,2);
	A = DDTinv * D1T * Y;
	DenseMatrix Ainv(2,2);

	if (A.Det() == 0) throw SingularMatrixException();

	CalcInverse (A, Ainv);

	for (int i = 0; i < 2; i++)
	  for (int j = 0; j < 2; j++)
	    {
	      Amat(i,j) = A(i,j);
	      Amatinv(i,j) = Ainv(i,j);
	    }

	Vec<2> temp = Amatinv * (psp2-psp1);
      

	double r = temp.Length();
	//      double alpha = -acos (temp(0)/r);
	double alpha = -atan2 (temp(1),temp(0));
	DenseMatrix R(2,2);
	R(0,0) = cos (alpha);
	R(1,0) = -sin (alpha);
	R(0,1) = sin (alpha);
	R(1,1) = cos (alpha);


	A = A*R;

	if (A.Det() == 0) throw SingularMatrixException();

	CalcInverse (A, Ainv);
    

	for (int i = 0; i < 2; i++)
	  for (int j = 0; j < 2; j++)
	    {
	      Amat(i,j) = A(i,j);
	      Amatinv(i,j) = Ainv(i,j);
	    }

	temp = Amatinv * (psp2-psp1);
      
      };
 
  }


  void OCCSurface :: ToPlane (const Point<3> & p3d,
			      const PointGeomInfo & geominfo,
			      Point<2> & pplane, 
			      double h, int & zone) const
  {
    if (projecttype == PLANESPACE)
      {
	Vec<3> p1p, n;
	GetNormalVector (p3d, geominfo, n);
      
	p1p = p3d - p1;
	pplane(0) = (p1p * ex) / h;
	pplane(1) = (p1p * ey) / h;
      
	if (n * nmid < 0)
	  zone = -1;
	else
	  zone = 0;

	/*
	  if(zone == -1)
	  {
	  (*testout) << "zone = -1 for " << p3d << " 2D: " << pplane << " n " << n << " nmid " << nmid << endl;
	  glob_testout = true;
	  GetNormalVector (p3d, geominfo, n);
	  glob_testout = false;
	  }
	*/
      }
    else
      {
	pplane = Point<2>(geominfo.u, geominfo.v);
	//      (*testout) << "(u,v) = " << geominfo.u << ", " << geominfo.v << endl;
	pplane = Point<2> (1/h * (Amatinv * (pplane-psp1)));
	//      pplane = Point<2> (h * (Amatinv * (pplane-psp1)));
	//      pplane = Point<2> (1/h * ((pplane-psp1)));

	zone = 0;
      };
  }	


  void OCCSurface :: FromPlane (const Point<2> & pplane, 
				Point<3> & p3d,
				PointGeomInfo & gi,
				double h) 
  { 
    if (projecttype == PLANESPACE)
      {
	//      cout << "2d   : " << pplane << endl;
	p3d = p1 + (h * pplane(0)) * ex + (h * pplane(1)) * ey;
	//      cout << "3d   : " << p3d << endl;
	Project (p3d, gi);  
	//      cout << "proj : " << p3d << endl;
      }
    else
      {
	//      Point<2> pspnew = Point<2>(1/h * (Amat * Vec<2>(pplane)) + Vec<2>(psp1));
	Point<2> pspnew = Point<2>(h * (Amat * Vec<2>(pplane)) + Vec<2>(psp1));
	//      Point<2> pspnew = Point<2>(h * (Vec<2>(pplane)) + Vec<2>(psp1));
	gi.u = pspnew(0);
	gi.v = pspnew(1);
	gi.trignum = 1;
	gp_Pnt val = occface->Value (gi.u, gi.v);
	p3d = Point<3> (val.X(), val.Y(), val.Z());
      };
  }



  void OCCSurface :: Project (Point<3> & p, PointGeomInfo & gi)
  {
    //   static int cnt = 0;
    //  if (cnt++ % 1000 == 0) cout << "********************************************** OCCSurfce :: Project, cnt = " << cnt << endl;
  
    gp_Pnt pnt(p(0), p(1), p(2));

    //(*testout) << "pnt = " << pnt.X() << ", " << pnt.Y() << ", " << pnt.Z() << endl;


    /*
    GeomAPI_ProjectPointOnSurf proj(pnt, occface, umin, umax, vmin, vmax);

    if (!proj.NbPoints())
      {
	cout << "Project Point on Surface FAIL" << endl;
	throw UVBoundsException();
      }
    */

    



    /*
      cout << "NP = " << proj.NbPoints() << endl;

      for (int i = 1; i <= proj.NbPoints(); i++)
      {
      gp_Pnt pnt2 = proj.Point(i);
      Point<3> p2 = Point<3> (pnt2.X(), pnt2.Y(), pnt2.Z());
      cout << i << ". p = " << p2 << ", dist = " << (p2-p).Length() << endl;
      }
    */

    /*
    pnt = proj.NearestPoint();
    proj.LowerDistanceParameters (gi.u, gi.v);
    */

    double u,v;
    Handle( ShapeAnalysis_Surface ) su = new ShapeAnalysis_Surface( occface );
    gp_Pnt2d suval = su->ValueOfUV ( pnt, BRep_Tool::Tolerance( topods_face ) );
    suval.Coord( u, v);
    pnt = occface->Value( u, v );
    
    //(*testout) << "pnt(proj) = " << pnt.X() << ", " << pnt.Y() << ", " << pnt.Z() << endl;
    gi.u = u;
    gi.v = v;
    

    gi.trignum = 1;

    p = Point<3> (pnt.X(), pnt.Y(), pnt.Z());
  } 


  Meshing2OCCSurfaces :: Meshing2OCCSurfaces (const TopoDS_Shape & asurf,
					      const Box<3> & abb, int aprojecttype)
    : Meshing2(mparam, Box<3>(abb.PMin(), abb.PMax())), surface(TopoDS::Face(asurf), aprojecttype)
  {
    ;
  }


  void Meshing2OCCSurfaces :: DefineTransformation (const Point3d & p1, const Point3d & p2,
						    const PointGeomInfo * geominfo1,
						    const PointGeomInfo * geominfo2)
  {
    ((OCCSurface&)surface).DefineTangentialPlane (p1, *geominfo1, p2, *geominfo2);
  }
 
  void Meshing2OCCSurfaces :: TransformToPlain (const Point3d & locpoint, 
						const MultiPointGeomInfo & geominfo,
						Point2d & planepoint, 
						double h, int & zone)
  {
    Point<2> hp;
    surface.ToPlane (locpoint, geominfo.GetPGI(1), hp, h, zone);
    planepoint.X() = hp(0);
    planepoint.Y() = hp(1);
  }

  int Meshing2OCCSurfaces :: TransformFromPlain (Point2d & planepoint,
						 Point3d & locpoint,
						 PointGeomInfo & gi,
						 double h)
  {
    Point<3> hp;
    Point<2> hp2 (planepoint.X(), planepoint.Y());
    surface.FromPlane (hp2, hp, gi, h);
    locpoint = hp;
    return 0;
  }



  double Meshing2OCCSurfaces :: CalcLocalH (const Point3d & p, double gh) const
  {
    return gh;
  }






  MeshOptimize2dOCCSurfaces :: MeshOptimize2dOCCSurfaces (const OCCGeometry & ageometry)
    : MeshOptimize2d(), geometry(ageometry)
  {
    ;
  }


  void MeshOptimize2dOCCSurfaces :: ProjectPoint (INDEX surfind, Point<3> & p) const
  {
    geometry.Project (surfind, p);
  }


  int MeshOptimize2dOCCSurfaces :: ProjectPointGI (INDEX surfind, Point<3> & p, PointGeomInfo & gi) const
  {
    double u = gi.u;
    double v = gi.v;

    Point<3> hp = p;
    if (geometry.FastProject (surfind, hp, u, v))
      {
	p = hp;
	return 1;
      }
    ProjectPoint (surfind, p); 
    return CalcPointGeomInfo (surfind, gi, p); 
  }


  void MeshOptimize2dOCCSurfaces :: ProjectPoint2 (INDEX surfind, INDEX surfind2, 
						   Point<3> & p) const
  {
    TopExp_Explorer exp0, exp1;
    bool done = false;
    Handle(Geom_Curve) c;

    for (exp0.Init(geometry.fmap(surfind), TopAbs_EDGE); !done && exp0.More(); exp0.Next())
      for (exp1.Init(geometry.fmap(surfind2), TopAbs_EDGE); !done && exp1.More(); exp1.Next())
	{
	  if (TopoDS::Edge(exp0.Current()).IsSame(TopoDS::Edge(exp1.Current())))
	    {
	      done = true;
	      double s0, s1;
	      c = BRep_Tool::Curve(TopoDS::Edge(exp0.Current()), s0, s1);
	    }
	}
  
    gp_Pnt pnt(p(0), p(1), p(2));
    GeomAPI_ProjectPointOnCurve proj(pnt, c);
    pnt = proj.NearestPoint();  
    p(0) = pnt.X();
    p(1) = pnt.Y();
    p(2) = pnt.Z();
	
  }

  void MeshOptimize2dOCCSurfaces :: 
  GetNormalVector(INDEX surfind, const Point<3> & p, PointGeomInfo & geominfo, Vec<3> & n) const
  {
    gp_Pnt pnt;
    gp_Vec du, dv;

    Handle(Geom_Surface) occface;
    occface = BRep_Tool::Surface(TopoDS::Face(geometry.fmap(surfind)));

    occface->D1(geominfo.u,geominfo.v,pnt,du,dv);

    n = Cross (Vec<3>(du.X(), du.Y(), du.Z()),
	       Vec<3>(dv.X(), dv.Y(), dv.Z()));
    n.Normalize();

    if (geometry.fmap(surfind).Orientation() == TopAbs_REVERSED) n = -1*n;  
  
    //  GetNormalVector (surfind, p, n);
  }


  void MeshOptimize2dOCCSurfaces :: 
  GetNormalVector(INDEX surfind, const Point<3> & p, Vec<3> & n) const
  {
    //  static int cnt = 0;
    //  if (cnt++ % 1000 == 0) cout << "GetNV cnt = " << cnt << endl;
    Standard_Real u,v;

    gp_Pnt pnt(p(0), p(1), p(2));

    Handle(Geom_Surface) occface;
    occface = BRep_Tool::Surface(TopoDS::Face(geometry.fmap(surfind)));

    /*
    GeomAPI_ProjectPointOnSurf proj(pnt, occface);

    if (proj.NbPoints() < 1)
      {
	cout << "ERROR: OCCSurface :: GetNormalVector: GeomAPI_ProjectPointOnSurf failed!"
	     << endl;
	cout << p << endl;
	return;
      }
 
    proj.LowerDistanceParameters (u, v);
    */
    
    Handle( ShapeAnalysis_Surface ) su = new ShapeAnalysis_Surface( occface );
    gp_Pnt2d suval = su->ValueOfUV ( pnt, BRep_Tool::Tolerance( TopoDS::Face(geometry.fmap(surfind)) ) );
    suval.Coord( u, v);
    pnt = occface->Value( u, v );
    
    

    gp_Vec du, dv;
    occface->D1(u,v,pnt,du,dv);

    /*
      if (!occface->IsCNu (1) || !occface->IsCNv (1))
      (*testout) << "SurfOpt: Differentiation FAIL" << endl;
    */

    n = Cross (Vec3d(du.X(), du.Y(), du.Z()),
	       Vec3d(dv.X(), dv.Y(), dv.Z()));
    n.Normalize();

    if (geometry.fmap(surfind).Orientation() == TopAbs_REVERSED) n = -1*n;  
  }


  int MeshOptimize2dOCCSurfaces :: 
  CalcPointGeomInfo(int surfind, PointGeomInfo& gi, const Point<3> & p) const
  {
    Standard_Real u,v;

    gp_Pnt pnt(p(0), p(1), p(2));

    Handle(Geom_Surface) occface;
    occface = BRep_Tool::Surface(TopoDS::Face(geometry.fmap(surfind)));

    /*
    GeomAPI_ProjectPointOnSurf proj(pnt, occface);

    if (proj.NbPoints() < 1)
      {
	cout << "ERROR: OCCSurface :: GetNormalVector: GeomAPI_ProjectPointOnSurf failed!"
	     << endl;
	cout << p << endl;
	return 0;
      }
 
    proj.LowerDistanceParameters (u, v);  
    */

    Handle( ShapeAnalysis_Surface ) su = new ShapeAnalysis_Surface( occface );
    gp_Pnt2d suval = su->ValueOfUV ( pnt, BRep_Tool::Tolerance( TopoDS::Face(geometry.fmap(surfind)) ) );
    suval.Coord( u, v);
    //pnt = occface->Value( u, v );
    

    gi.u = u;
    gi.v = v;
    return 1;
  }






  OCCRefinementSurfaces :: OCCRefinementSurfaces (const OCCGeometry & ageometry)
    : Refinement(), geometry(ageometry)
  {
    ;
  }

  OCCRefinementSurfaces :: ~OCCRefinementSurfaces ()
  {
    ;
  }

  /*
    inline double Det3 (double a00, double a01, double a02,
    double a10, double a11, double a12,
    double a20, double a21, double a22)
    {
    return a00*a11*a22 + a01*a12*a20 + a10*a21*a02 - a20*a11*a02 - a10*a01*a22 - a21*a12*a00;
    }

    bool ProjectToSurface (gp_Pnt & p, Handle(Geom_Surface) surface, double& u, double& v)
    {
    gp_Pnt x = surface->Value (u,v);

    if (p.SquareDistance(x) <= sqr(PROJECTION_TOLERANCE)) return true;

    gp_Vec du, dv;

    surface->D1(u,v,x,du,dv);

    int count = 0;

    gp_Pnt xold;
    gp_Vec n;
    double det, lambda, mu;

    do {
    count++;

    n = du^dv;

    det = Det3 (n.X(), du.X(), dv.X(),
    n.Y(), du.Y(), dv.Y(),
    n.Z(), du.Z(), dv.Z());

    if (det < 1e-15) return false; 

    lambda = Det3 (n.X(), p.X()-x.X(), dv.X(),
    n.Y(), p.Y()-x.Y(), dv.Y(),
    n.Z(), p.Z()-x.Z(), dv.Z())/det;

    mu     = Det3 (n.X(), du.X(), p.X()-x.X(),
    n.Y(), du.Y(), p.Y()-x.Y(),
    n.Z(), du.Z(), p.Z()-x.Z())/det;
  
    u += lambda;
    v += mu;

    xold = x;
    surface->D1(u,v,x,du,dv);

    } while (xold.SquareDistance(x) > sqr(PROJECTION_TOLERANCE) || count > 50);

    if (count > 50) return false;

    p = x;

    return true;
    }
  */

  void OCCRefinementSurfaces :: 
  PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
		int surfi, 
		const PointGeomInfo & gi1, 
		const PointGeomInfo & gi2,
		Point<3> & newp, PointGeomInfo & newgi) const
  {
    Point<3> hnewp;
    hnewp = p1+secpoint*(p2-p1);

    if (surfi > 0)
      {
      
	double u = gi1.u+secpoint*(gi2.u-gi1.u);
	double v = gi1.v+secpoint*(gi2.v-gi1.v);
 
	if (!geometry.FastProject (surfi, hnewp, u, v))
	  {
	  //  cout << "Fast projection to surface fails! Using OCC projection" << endl;
	    geometry.Project (surfi, hnewp);
	  }

	newgi.trignum = 1;
        newgi.u = u;
        newgi.v = v;
      }
  
    newp = hnewp;
  }


  void OCCRefinementSurfaces :: 
  PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
		int surfi1, int surfi2, 
		const EdgePointGeomInfo & ap1, 
		const EdgePointGeomInfo & ap2,
		Point<3> & newp, EdgePointGeomInfo & newgi) const
  {
    double s0, s1;

    Point<3> hnewp = p1+secpoint*(p2-p1);
    gp_Pnt pnt(hnewp(0), hnewp(1), hnewp(2));
    GeomAPI_ProjectPointOnCurve proj(pnt, BRep_Tool::Curve(TopoDS::Edge(geometry.emap(ap1.edgenr)), s0, s1));
    pnt = proj.NearestPoint();
    hnewp = Point<3> (pnt.X(), pnt.Y(), pnt.Z());
    newp = hnewp;
    newgi = ap1;
  };


  void OCCRefinementSurfaces :: ProjectToSurface (Point<3> & p, int surfi) const
  {
    if (surfi > 0)
      geometry.Project (surfi, p);
  };

  void OCCRefinementSurfaces :: ProjectToSurface (Point<3> & p, int surfi, PointGeomInfo & gi) const
  {
    if (surfi > 0)
      if (!geometry.FastProject (surfi, p, gi.u, gi.v))
	{
	  cout << "Fast projection to surface fails! Using OCC projection" << endl;
	  geometry.Project (surfi, p);
	}
  };



}


#endif
