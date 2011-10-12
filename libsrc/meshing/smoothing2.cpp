#include <mystdlib.h>

#include "meshing.hpp"
#include <opti.hpp>

namespace netgen
{

  static const double c_trig = 0.14433756;      // sqrt(3.0) / 12
  static const double c_trig4 = 0.57735026;     // sqrt(3.0) / 3

  inline double CalcTriangleBadness (double x2, double x3, double y3, 
				     double metricweight, double h)
  {
    // badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1 
    // p1 = (0, 0), p2 = (x2, 0), p3 = (x3, y3);

    double cir_2 = (x2*x2 + x3*x3 + y3*y3 - x2*x3);
    double area = x2 * y3;
    
    if (area <= 1e-24 * cir_2)
      return 1e10;
    
    double badness = c_trig4 * cir_2 / area - 1;
    
    if (metricweight > 0)
      {
	// add:  metricweight * (area / h^2 + h^2 / area - 2)

	double areahh = area / (h * h);
	badness += metricweight * (areahh + 1 / areahh - 2);
      }
    return badness;
  }

  
  inline void CalcTriangleBadness (double x2, double x3, double y3, double metricweight,
				   double h, double & badness, double & g1x, double & g1y)
  {
    // old: badness = sqrt(3.0) /36 * circumference^2 / area - 1 
    // badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1 
    // p1 = (0, 0), p2 = (x2, 0), p3 = (x3, y3);


    double cir_2 = 2* (x2*x2 + x3*x3 + y3*y3 - x2*x3);
    double area = 0.5 * x2 * y3;

    if (area <= 1e-24 * cir_2)
      {
	g1x = 0;
	g1y = 0;
	badness = 1e10;
	return;
      }

    badness = c_trig * cir_2 / area - 1;

    double c1 = -2 * c_trig / area;
    double c2 = 0.5 * c_trig * cir_2 / (area * area);
    g1x = c1 * (x2 + x3) + c2 * y3;
    g1y = c1 * (y3)      + c2 * (x2-x3); 

    if (metricweight > 0)
      {
	// area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
	// add:  metricweight * (area / h^2 + h^2 / area - 2)
      
	area = x2 * y3;
	double dareax1 = -y3; 
	double dareay1 = x3 - x2; 

	double areahh = area / (h * h);
	double fac = metricweight * (areahh - 1 / areahh) / area;

	badness += metricweight * (areahh + 1 / areahh - 2);
	g1x += fac * dareax1;
	g1y += fac * dareay1; 
      }
  }





  double CalcTriangleBadness (const Point3d & p1, 
			      const Point3d & p2, 
			      const Point3d & p3,
			      double metricweight,
			      double h)
  {
    // badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1 
    // p1 = (0, 0), p2 = (x2, 0), p3 = (x3, y3);

    Vec3d e12(p1,p2);
    Vec3d e13(p1,p3);
    Vec3d e23(p2,p3);
  
    double l12_2 = e12.Length2();
    double l13_2 = e13.Length2();
    double l23_2 = e23.Length2();

    double cir_2 = l12_2 + l13_2 + l23_2;
    Vec3d area_v = Cross (e12, e13);
    double area = 0.5 * area_v.Length();

    if (area <= 1e-24 * cir_2)
      return 1e10;

    double badness = c_trig * cir_2 / area - 1;

    if (metricweight > 0)
      {
	// area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
	// add:  metricweight * (area / h^2 + h^2 / area - 2)

	const double areahh = area / (h * h);
	badness += metricweight * (areahh + 1 / areahh - 2);
      }

    return badness;
  }


  double CalcTriangleBadness (const Point3d & p1, 
			      const Point3d & p2, 
			      const Point3d & p3,
			      const Vec3d & n,
			      double metricweight,
			      double h)
  {
    Vec3d v1 (p1, p2);
    Vec3d v2 (p1, p3);

    Vec3d e1 = v1;
    Vec3d e2 = v2;

    e1 -= (e1 * n) * n;
    e1 /= (e1.Length() + 1e-24);
    e2 = Cross (n, e1);

    return CalcTriangleBadness ( (e1 * v1), (e1 * v2), (e2 * v2), 
				 metricweight, h);
  }



  class Opti2dLocalData
  {
  public:
    const MeshOptimize2d * meshthis;
    MeshPoint sp1; 
    PointGeomInfo gi1;
    Vec<3> normal, t1, t2;
    Array<SurfaceElementIndex> locelements;
    Array<int> locrots;
    Array<double> lochs;
  // static int locerr2;
    double locmetricweight;
    double loch;
    int surfi, surfi2;
    int uselocalh;
  public:
    Opti2dLocalData ()
    {
      locmetricweight = 0;
    }
  };


  class Opti2SurfaceMinFunction : public MinFunction
  {
    const Mesh & mesh;
    Opti2dLocalData & ld;
  public:
    Opti2SurfaceMinFunction (const Mesh & amesh,
			     Opti2dLocalData & ald)
      : mesh(amesh), ld(ald)
    { } ;
    virtual double FuncGrad (const Vector & x, Vector & g) const;
    virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
    virtual double Func (const Vector & x) const;
  };
  
  double Opti2SurfaceMinFunction :: 
  Func (const Vector & x) const
  {
    Vector g(x.Size());
    return FuncGrad (x, g);
  }


  double Opti2SurfaceMinFunction :: 
  FuncGrad (const Vector & x, Vector & grad) const
  {
    Vec<3> n, vgrad;
    Point<3> pp1;
    double g1x, g1y;
    double badness, hbadness;

    vgrad = 0;
    badness = 0;

    ld.meshthis -> GetNormalVector (ld.surfi, ld.sp1, ld.gi1, n);
    pp1 = ld.sp1 + x(0) * ld.t1 + x(1) * ld.t2;

    //  meshthis -> ProjectPoint (surfi, pp1);
    // meshthis -> GetNormalVector (surfi, pp1, n);

    for (int j = 0; j < ld.locelements.Size(); j++)
      {
	int roti = ld.locrots[j];
	const Element2d & bel = mesh[ld.locelements[j]];

	Vec<3> e1 = mesh[bel.PNumMod(roti + 1)] - pp1;
	Vec<3> e2 = mesh[bel.PNumMod(roti + 2)] - pp1;

	if (ld.uselocalh) ld.loch = ld.lochs[j];

	double e1l = e1.Length();
	if (Determinant(e1, e2, n) > 1e-8 * e1l * e2.Length())
	  {
	    e1 /= e1l;
	    double e1e2 = e1 * e2;
            e2 -= e1e2 * e1;
	    double e2l = e2.Length();

	    CalcTriangleBadness ( e1l, e1e2, e2l, ld.locmetricweight, ld.loch,
				  hbadness, g1x, g1y);

	    badness += hbadness;
            vgrad += g1x * e1 + (g1y/e2l) * e2;
	  }
	else
	  {
	    // (*testout) << "very very bad badness" << endl;
	    badness += 1e8;
	  }
      }

    vgrad -=  (vgrad * n) * n;

    grad(0) = vgrad * ld.t1;
    grad(1) = vgrad * ld.t2;
    return badness;
  }




  double Opti2SurfaceMinFunction :: 
  FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
  {
    Vec<3> n, vgrad;
    Point<3> pp1;
    double g1x, g1y;
    double badness, hbadness;

    vgrad = 0;
    badness = 0;

    ld.meshthis -> GetNormalVector (ld.surfi, ld.sp1, ld.gi1, n);

    pp1 = ld.sp1 + x(0) * ld.t1 + x(1) * ld.t2;

    for (int j = 0; j < ld.locelements.Size(); j++)
      {
	int roti = ld.locrots[j];

	const Element2d & bel = mesh[ld.locelements[j]];

	Vec<3> e1 = mesh[bel.PNumMod(roti + 1)] - pp1;
	Vec<3> e2 = mesh[bel.PNumMod(roti + 2)] - pp1;

	if (ld.uselocalh) ld.loch = ld.lochs[j];

	double e1l = e1.Length();
	if (Determinant(e1, e2, n) > 1e-8 * e1l * e2.Length())
	  {
	    e1 /= e1l;
	    double e1e2 = e1 * e2;
	    e2 -= e1e2 * e1;
	    double e2l = e2.Length();
	    CalcTriangleBadness ( e1l, e1e2, e2l, ld.locmetricweight, ld.loch,
				  hbadness, g1x, g1y);

	    badness += hbadness;
            vgrad += g1x * e1 + (g1y / e2l) * e2;
	  }
	else
	  {
	    // (*testout) << "very very bad badness" << endl;
	    badness += 1e8;
	  }
      }

    vgrad -= (vgrad * n) * n;
    deriv = dir(0) * (vgrad*ld.t1) + dir(1) * (vgrad*ld.t2);

    return badness;
  }












  class Opti2EdgeMinFunction : public MinFunction
  {
    const Mesh & mesh;
    Opti2dLocalData & ld;

  public:
    Opti2EdgeMinFunction (const Mesh & amesh,
			  Opti2dLocalData & ald)
      : mesh(amesh), ld(ald) { } ;

    virtual double FuncGrad (const Vector & x, Vector & g) const;
    virtual double Func (const Vector & x) const;
  };

  double Opti2EdgeMinFunction :: Func (const Vector & x) const
  {
    Vector g(x.Size());
    return FuncGrad (x, g);
  }

  double Opti2EdgeMinFunction :: FuncGrad (const Vector & x, Vector & grad) const
  {
    int j, rot;
    Vec<3> n1, n2, v1, v2, e1, e2, vgrad;
    Point<3> pp1;
    Vec<2> g1;
    double badness, hbadness;

    vgrad = 0.0;
    badness = 0;

    pp1 = ld.sp1 + x(0) * ld.t1;
    ld.meshthis -> ProjectPoint2 (ld.surfi, ld.surfi2, pp1);

    for (j = 0; j < ld.locelements.Size(); j++)
      {
	rot = ld.locrots[j];
	const Element2d & bel = mesh[ld.locelements[j]];

	v1 = mesh[bel.PNumMod(rot + 1)] - pp1;
	v2 = mesh[bel.PNumMod(rot + 2)] - pp1;

	e1 = v1;
	e2 = v2;
	e1 /= e1.Length();
	e2 -= (e1 * e2) * e1;
	e2 /= e2.Length();

	if (ld.uselocalh) ld.loch = ld.lochs[j];
	CalcTriangleBadness ( (e1 * v1), (e1 * v2), (e2 * v2), ld.locmetricweight, ld.loch,
			      hbadness, g1(0), g1(1));

	badness += hbadness;
        vgrad += g1(0) * e1 + g1(1) * e2;
      }

    ld.meshthis -> GetNormalVector (ld.surfi, pp1, n1);
    ld.meshthis -> GetNormalVector (ld.surfi2, pp1, n2);

    v1 = Cross (n1, n2);
    v1.Normalize();

    grad(0) = (vgrad * v1) * (ld.t1 * v1);

    return badness;
  }




  class Opti2SurfaceMinFunctionJacobian : public MinFunction
  {
    const Mesh & mesh;
    Opti2dLocalData & ld;

  public:
    Opti2SurfaceMinFunctionJacobian (const Mesh & amesh,
				     Opti2dLocalData & ald)
      : mesh(amesh), ld(ald)
    { } ;
    virtual double FuncGrad (const Vector & x, Vector & g) const;
    virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
    virtual double Func (const Vector & x) const;
  };
  
  double Opti2SurfaceMinFunctionJacobian :: 
  Func (const Vector & x) const
  {
    Vector g(x.Size());
    return FuncGrad (x, g);
  }


  double Opti2SurfaceMinFunctionJacobian :: 
  FuncGrad (const Vector & x, Vector & grad) const
  {
    // from 2d:

    int lpi, gpi;
    Vec<3> n, vgrad;
    Point<3> pp1;
    Vec2d g1, vdir;
    double badness, hbad, hderiv;

    vgrad = 0;
    badness = 0;

    ld.meshthis -> GetNormalVector (ld.surfi, ld.sp1, ld.gi1, n);

    pp1 = ld.sp1 + x(0) * ld.t1 + x(1) * ld.t2;

    //  meshthis -> ProjectPoint (surfi, pp1);
    //  meshthis -> GetNormalVector (surfi, pp1, n);

    static Array<Point2d> pts2d;
    pts2d.SetSize(mesh.GetNP());

    grad = 0;

    for (int j = 1; j <= ld.locelements.Size(); j++)
      {
	lpi = ld.locrots.Get(j);
	const Element2d & bel = 
	  mesh[ld.locelements.Get(j)];
      
	gpi = bel.PNum(lpi);

	for (int k = 1; k <= bel.GetNP(); k++)
	  {
	    PointIndex pi = bel.PNum(k);
	    pts2d.Elem(pi) = Point2d (ld.t1 * (mesh.Point(pi) - ld.sp1), 
				      ld.t2 * (mesh.Point(pi) - ld.sp1)); 
	  }				    
	pts2d.Elem(gpi) = Point2d (x(0), x(1));
      

	for (int k = 1; k <= 2; k++)
	  {
	    if (k == 1)
	      vdir = Vec2d (1, 0);
	    else
	      vdir = Vec2d (0, 1);
	  
	    hbad = bel.
	      CalcJacobianBadnessDirDeriv (pts2d, lpi, vdir, hderiv);
            
	    grad(k-1) += hderiv;
	    if (k == 1)
	      badness += hbad;
	  }
      }


    /*
      vgrad.Add (-(vgrad * n), n);

      grad.Elem(1) = vgrad * t1;
      grad.Elem(2) = vgrad * t2;
    */
    return badness;
  }




  double Opti2SurfaceMinFunctionJacobian :: 
  FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
  {
    // from 2d:

    int j, k, lpi, gpi;
    Vec<3> n, vgrad;
    Point<3> pp1;
    Vec2d g1, vdir;
    double badness, hbad, hderiv;

    vgrad = 0;
    badness = 0;

    ld.meshthis -> GetNormalVector (ld.surfi, ld.sp1, ld.gi1, n);

    // pp1 = sp1;
    //    pp1.Add2 (x.Get(1), t1, x.Get(2), t2);
    pp1 = ld.sp1 + x(0) * ld.t1 + x(1) * ld.t2;

    static Array<Point2d> pts2d;
    pts2d.SetSize(mesh.GetNP());

    deriv = 0;

    for (j = 1; j <= ld.locelements.Size(); j++)
      {
	lpi = ld.locrots.Get(j);
	const Element2d & bel = 
	  mesh[ld.locelements.Get(j)];
      
	gpi = bel.PNum(lpi);

	for (k = 1; k <= bel.GetNP(); k++)
	  {
	    PointIndex pi = bel.PNum(k);
	    pts2d.Elem(pi) = Point2d (ld.t1 * (mesh.Point(pi) - ld.sp1), 
				      ld.t2 * (mesh.Point(pi) - ld.sp1)); 
	  }				    
	pts2d.Elem(gpi) = Point2d (x(0), x(1));
      

	vdir = Vec2d (dir(0), dir(1));
	  
	hbad = bel.
	  CalcJacobianBadnessDirDeriv (pts2d, lpi, vdir, hderiv);
      
	deriv += hderiv;
	badness += hbad;
      }


    return badness;
  }







  MeshOptimize2d dummy;

  MeshOptimize2d :: MeshOptimize2d ()
  {
    SetFaceIndex (0);
    SetImproveEdges (0);
    SetMetricWeight (0);
    SetWriteStatus (1);
  }


  void MeshOptimize2d :: SelectSurfaceOfPoint (const Point<3> & p,
					       const PointGeomInfo & gi)
  {
    ;
  }

  void MeshOptimize2d :: ImproveMesh (Mesh & mesh, const MeshingParameters & mp)
  {
    if (!faceindex)
      {
	PrintMessage (3, "Smoothing");

	for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
	  {
	    ImproveMesh (mesh, mp);
	    if (multithread.terminate)
	      throw NgException ("Meshing stopped");
	  }
	faceindex = 0;
	return;
      }
 

    static int timer = NgProfiler::CreateTimer ("MeshSmoothing 2D");
    static int timer1 = NgProfiler::CreateTimer ("MeshSmoothing 2D start");

    NgProfiler::RegionTimer reg (timer);
    NgProfiler::StartTimer (timer1);

    CheckMeshApproximation (mesh);

    Opti2dLocalData ld;

    // SurfaceElementIndex sei;

    Array<SurfaceElementIndex> seia;
    mesh.GetSurfaceElementsOfFace (faceindex, seia);

    bool mixed = 0;
    for (int i = 0; i < seia.Size(); i++)
      if (mesh[seia[i]].GetNP() != 3)
	{
	  mixed = 1;
	  break;
	}

    int loci;
    double fact;
    int moveisok;

    PointGeomInfo ngi;
    Point3d origp;



    Vec3d n1, n2;
    Vector x(2), xedge(1);

    Array<MeshPoint, PointIndex::BASE> savepoints(mesh.GetNP());

    ld.uselocalh = mp.uselocalh;


    Array<int, PointIndex::BASE> nelementsonpoint(mesh.GetNP());
    nelementsonpoint = 0;
    for (int i = 0; i < seia.Size(); i++)
      {
	const Element2d & el = mesh[seia[i]];
	for (int j = 0; j < el.GetNP(); j++)
	  nelementsonpoint[el[j]]++;
      }

    TABLE<SurfaceElementIndex,PointIndex::BASE> elementsonpoint(nelementsonpoint);

    for (int i = 0; i < seia.Size(); i++)
      {
	const Element2d & el = mesh[seia[i]];
	for (int j = 0; j < el.GetNP(); j++)
	  elementsonpoint.Add (el[j], seia[i]);
      }


    ld.loch = mp.maxh;
    ld.locmetricweight = metricweight;
    ld.meshthis = this;



    Opti2SurfaceMinFunction surfminf(mesh, ld);
    Opti2EdgeMinFunction edgeminf(mesh, ld);
    Opti2SurfaceMinFunctionJacobian surfminfj(mesh, ld);

    OptiParameters par;
    par.maxit_linsearch = 8;
    par.maxit_bfgs = 5;

    /*
    int i, j, k;

      if (improveedges)
      for (i = 1; i <= mesh.GetNP(); i++)
      if (mesh.PointType(i) == EDGEPOINT)
      {
      continue;
      PrintDot ();
      sp1 = mesh.Point(i);
	  
      locelements.SetSize(0);
      locrots.SetSize (0);
      lochs.SetSize (0);
      surfi = surfi2 = surfi3 = 0;
	  
      for (j = 0; j < elementsonpoint[i].Size(); j++)
      {
      sei = elementsonpoint[i][j];
      const Element2d * bel = &mesh[sei];
	      
      if (!surfi)
      surfi = mesh.GetFaceDescriptor(bel->GetIndex()).SurfNr();
      else if (surfi != mesh.GetFaceDescriptor(bel->GetIndex()).SurfNr())
      {
      if (surfi2 != 0 && surfi2 != 
      mesh.GetFaceDescriptor(bel->GetIndex()).SurfNr())
      surfi3 = mesh.GetFaceDescriptor(bel->GetIndex()).SurfNr();
      else
      surfi2 = mesh.GetFaceDescriptor(bel->GetIndex()).SurfNr();
      }
	      
      locelements.Append (sei);
	      
      if (bel->PNum(1) == i)
      locrots.Append (1);
      else if (bel->PNum(2) == i)
      locrots.Append (2);
      else
      locrots.Append (3);

      if (uselocalh)
      {
      Point3d pmid = Center (mesh.Point(bel->PNum(1)),
      mesh.Point(bel->PNum(2)),
      mesh.Point(bel->PNum(3)));
      lochs.Append (mesh.GetH(pmid));
      }
      }
	  
      if (surfi2 && !surfi3)
      {
      GetNormalVector (surfi, sp1, n1);
      GetNormalVector (surfi2, sp1, n2);
      t1 = Cross (n1, n2);
	      
      xedge = 0;
      BFGS (xedge, edgeminf, par, 1e-6);
	      
      mesh.Point(i).X() += xedge.Get(1) * t1.X();
      mesh.Point(i).Y() += xedge.Get(1) * t1.Y();
      mesh.Point(i).Z() += xedge.Get(1) * t1.Z();
      ProjectPoint2 (surfi, surfi2, mesh.Point(i));
      }
      }
    */


    bool printeddot = 0;
    char plotchar = '.';
    int modplot = 1;
    if (mesh.GetNP() > 1000)
      {
	plotchar = '+';
	modplot = 10;
      }
    if (mesh.GetNP() > 10000)
      {
	plotchar = 'o';
	modplot = 100;
      }
    int cnt = 0;


    NgProfiler::StopTimer (timer1);


    for (PointIndex pi = PointIndex::BASE; pi < mesh.GetNP()+PointIndex::BASE; pi++)
      if (mesh[pi].Type() == SURFACEPOINT)
	{
	  if (multithread.terminate)
	    throw NgException ("Meshing stopped");
	
	  cnt++;
	  if (cnt % modplot == 0 && writestatus)
	    {
	      printeddot = 1;
	      PrintDot (plotchar);
	    }

	  if (elementsonpoint[pi].Size() == 0)
	    continue;

	
	  ld.sp1 = mesh[pi];

	  Element2d & hel = mesh[elementsonpoint[pi][0]];

	  int hpi = 0;
	  for (int j = 1; j <= hel.GetNP(); j++)
	    if (hel.PNum(j) == pi)
	      {
		hpi = j;
		break;
	      }

	  ld.gi1 = hel.GeomInfoPi(hpi);
	  SelectSurfaceOfPoint (ld.sp1, ld.gi1);
	  
	  ld.locelements.SetSize(0);
	  ld.locrots.SetSize (0);
	  ld.lochs.SetSize (0);
	
	  for (int j = 0; j < elementsonpoint[pi].Size(); j++)
	    {
	      SurfaceElementIndex sei = elementsonpoint[pi][j];
	      const Element2d & bel = mesh[sei];
	      ld.surfi = mesh.GetFaceDescriptor(bel.GetIndex()).SurfNr();
	    
	      ld.locelements.Append (sei);
	    
	      for (int k = 1; k <= bel.GetNP(); k++)
		if (bel.PNum(k) == pi)
		  {
		    ld.locrots.Append (k);
		    break;
		  }
	      
	      if (ld.uselocalh)
		{
		  Point3d pmid = Center (mesh[bel[0]], mesh[bel[1]], mesh[bel[2]]);
		  ld.lochs.Append (mesh.GetH(pmid));
		}
	    }
	  
	  GetNormalVector (ld.surfi, ld.sp1, ld.gi1, ld.normal);
	  ld.t1 = ld.normal.GetNormal ();
	  ld.t2 = Cross (ld.normal, ld.t1);
	  
	  // save points, and project to tangential plane
	  for (int j = 0; j < ld.locelements.Size(); j++)
	    {
	      const Element2d & el = mesh[ld.locelements[j]];
	      for (int k = 0; k < el.GetNP(); k++)
		savepoints[el[k]] = mesh[el[k]];
	    }

	  for (int j = 0; j < ld.locelements.Size(); j++)
	    {
	      const Element2d & el = mesh[ld.locelements[j]];
	      for (int k = 0; k < el.GetNP(); k++)
		{
		  PointIndex hhpi = el[k];
		  double lam = ld.normal * (mesh[hhpi] - ld.sp1);
		  mesh[hhpi] -= lam * ld.normal;
		}
	    }
	  
	  x = 0;
	  par.typx = ld.lochs[0];

	  if (mixed)
	    {
	      //	      (*testout) << "vorher : " << surfminfj.Func (x) << endl;
	      BFGS (x, surfminfj, par, 1e-6);
	      //	      (*testout) << "nachher: " << surfminfj.Func (x) << endl;
	      //	      (*testout) << "x = " << x << endl;
	    }
	  else
	    {
	      //	      (*testout) << "vorher : " << surfminf.Func (x) << endl;
	      BFGS (x, surfminf, par, 1e-6);
	      //	      (*testout) << "nachher: " << surfminf.Func (x) << endl;
	      //	      (*testout) << "x = " << x << endl;
	    }

	
	  origp = mesh[pi];
	  loci = 1;
	  fact = 1;
	  moveisok = 0;

	  // restore other points
	  for (int j = 0; j < ld.locelements.Size(); j++)
	    {
	      const Element2d & el = mesh[ld.locelements[j]];
	      for (int k = 0; k < el.GetNP(); k++)
		{
		  PointIndex hhpi = el[k];
		  if (hhpi != pi) mesh[hhpi] = savepoints[hhpi];
		}
	    }

	  
	  //optimizer loop (if whole distance is not possible, move only a bit!!!!)
	  while (loci <= 5 && !moveisok)
	    {
	      loci ++;
              /*
	      mesh[pi].X() = origp.X() + (x.Get(1) * t1.X() + x.Get(2) * t2.X())*fact;
	      mesh[pi].Y() = origp.Y() + (x.Get(1) * t1.Y() + x.Get(2) * t2.Y())*fact;
	      mesh[pi].Z() = origp.Z() + (x.Get(1) * t1.Z() + x.Get(2) * t2.Z())*fact;
              */
              Vec<3> hv = x(0) * ld.t1 + x(1) * ld.t2;
              Point3d hnp = origp + Vec3d (hv);
              mesh[pi](0) = hnp.X();
              mesh[pi](1) = hnp.Y();
              mesh[pi](2) = hnp.Z();

	      fact = fact/2.;

	      // ProjectPoint (surfi, mesh[pi]);
	      // moveisok = CalcPointGeomInfo(surfi, ngi, mesh[pi]); 

	      ngi = ld.gi1;
	      moveisok = ProjectPointGI (ld.surfi, mesh[pi], ngi);
	      // point lies on same chart in stlsurface
	    
	      if (moveisok)
		{
		  for (int j = 0; j < ld.locelements.Size(); j++)
		    mesh[ld.locelements[j]].GeomInfoPi(ld.locrots[j]) = ngi;
		}
	      else
		{
		  mesh[pi] = Point<3> (origp);
		}
	    
	    }
	}

    if (printeddot)
      PrintDot ('\n');
  
    CheckMeshApproximation (mesh);
    mesh.SetNextTimeStamp();
  }

  void MeshOptimize2d :: GetNormalVector(INDEX /* surfind */, const Point<3> & p, Vec<3> & nv) const
  {
    nv = Vec<3> (0, 0, 1);
  }

  void MeshOptimize2d :: GetNormalVector(INDEX surfind, const Point<3> & p, PointGeomInfo & gi, Vec<3> & n) const
  {
    GetNormalVector (surfind, p, n);
  }
}
