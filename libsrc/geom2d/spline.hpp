#ifndef FILE_SPLINE_HPP
#define FILE_SPLINE_HPP

/**************************************************************************/
/* File:   spline.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   24. Jul. 96                                                    */
/**************************************************************************/

namespace netgen
{


void CalcPartition (double l, double h, double h1, double h2,
		    double hcurve, double elto0, Array<double> & points);

/*
  Spline curves for 2D mesh generation
*/


/// Geometry point
template < int D >
class GeomPoint : public Point<D>
{
public:
  /// refinement factor at point
  double refatpoint;
  /// max mesh-size at point
  double hmax;
  /// hp-refinement
  bool hpref;

  ///
  GeomPoint () { ; }

  ///
  GeomPoint (const Point<D> & ap, double aref = 1, bool ahpref=false)
  : Point<D>(ap), refatpoint(aref), hpref(ahpref) { ; }
};




/// base class for 2d - segment
template < int D >
class SplineSeg
{
public:
  /// left domain
  int leftdom;
  /// right domain
  int rightdom;
  /// refinement at line
  double reffak;
  /// maximal h;
  double hmax;
  /// boundary condition number
  int bc;
  /// copy spline mesh from other spline (-1.. do not copy)
  int copyfrom;
  /// perfrom anisotropic refinement (hp-refinement) to edge
  bool hpref_left;
  /// perfrom anisotropic refinement (hp-refinement) to edge
  bool hpref_right;
  ///
  int layer;


  SplineSeg ()
  {
    layer = 1;
  }

  /// calculates length of curve
  virtual double Length () const;
  /// returns point at curve, 0 <= t <= 1
  virtual Point<D> GetPoint (double t) const = 0;
  /// returns a (not necessarily unit-length) tangent vector for 0 <= t <= 1
  virtual Vec<D> GetTangent (const double t) const
  { cerr << "GetTangent not implemented for spline base-class"  << endl; Vec<D> dummy; return dummy;}
  virtual void GetDerivatives (const double t, 
			       Point<D> & point,
			       Vec<D> & first,
			       Vec<D> & second) const {;}
  /// partitionizes curve
  void Partition (double h, double elto0,
		  Mesh & mesh, Point3dTree & searchtree, int segnr) const;
  /// returns initial point on curve
  virtual const GeomPoint<D> & StartPI () const = 0;
  /// returns terminal point on curve
  virtual const GeomPoint<D> & EndPI () const = 0;
  /** writes curve description for fepp:
      for implicitly given quadratic curves, the 6 coefficients of
      the polynomial
      $$ a x^2 + b y^2 + c x y + d x + e y + f = 0 $$
      are written to ost */
  void PrintCoeff (ostream & ost) const;

  virtual void GetCoeff (Vector & coeffs) const = 0;

  virtual void GetPoints (int n, Array<Point<D> > & points);

  /** calculates (2D) lineintersections:
      for lines $$ a x + b y + c = 0 $$ the interecting points are calculated
      and stored in points */
  virtual void LineIntersections (const double a, const double b, const double c,
				  Array < Point<D> > & points, const double eps) const
  {points.SetSize(0);}

  virtual double MaxCurvature(void) const = 0;

  virtual string GetType(void) const {return "splinebase";}

  virtual void Project (const Point<D> point, Point<D> & point_on_curve, double & t) const
  { cerr << "Project not implemented for spline base-class" << endl;}

  virtual void GetRawData (Array<double> & data) const
  { cerr << "GetRawData not implemented for spline base-class" << endl;}

};


/// Straight line form p1 to p2
template< int D >
class LineSeg : public SplineSeg<D>
{
  ///
  GeomPoint<D> p1, p2;
public:
  ///
  LineSeg (const GeomPoint<D> & ap1, const GeomPoint<D> & ap2);
  ///
  virtual double Length () const;
  ///
  inline virtual Point<D> GetPoint (double t) const;
  ///
  virtual Vec<D> GetTangent (const double t) const;

  
  virtual void GetDerivatives (const double t, 
			       Point<D> & point,
			       Vec<D> & first,
			       Vec<D> & second) const;
  ///
  virtual const GeomPoint<D> & StartPI () const { return p1; };
  ///
  virtual const GeomPoint<D> & EndPI () const { return p2; }
  ///
  virtual void GetCoeff (Vector & coeffs) const;

  virtual string GetType(void) const {return "line";}

  virtual void LineIntersections (const double a, const double b, const double c,
				  Array < Point<D> > & points, const double eps) const;

  virtual double MaxCurvature(void) const {return 0;}

  virtual void Project (const Point<D> point, Point<D> & point_on_curve, double & t) const;

  virtual void GetRawData (Array<double> & data) const;
};


/// curve given by a rational, quadratic spline (including ellipses)
template< int D >
class SplineSeg3 : public SplineSeg<D>
{
  ///
  GeomPoint<D> p1, p2, p3;

  mutable double proj_latest_t;
public:
  ///
  SplineSeg3 (const GeomPoint<D> & ap1, 
	      const GeomPoint<D> & ap2, 
	      const GeomPoint<D> & ap3);
  ///
  inline virtual Point<D> GetPoint (double t) const;
  ///
  virtual Vec<D> GetTangent (const double t) const;

  
  virtual void GetDerivatives (const double t, 
			       Point<D> & point,
			       Vec<D> & first,
			       Vec<D> & second) const;
  ///
  virtual const GeomPoint<D> & StartPI () const { return p1; };
  ///
  virtual const GeomPoint<D> & EndPI () const { return p3; }
  ///
  virtual void GetCoeff (Vector & coeffs) const;

  virtual string GetType(void) const {return "spline3";}

  const GeomPoint<D> & TangentPoint (void) const { return p2; }

  virtual void LineIntersections (const double a, const double b, const double c,
				  Array < Point<D> > & points, const double eps) const;

  virtual double MaxCurvature(void) const;

  virtual void Project (const Point<D> point, Point<D> & point_on_curve, double & t) const;

  virtual void GetRawData (Array<double> & data) const;
};


// Gundolf Haase  8/26/97
/// A circle
template < int D >
class CircleSeg : public SplineSeg<D>
{
  ///
private:
  GeomPoint<D>	p1, p2, p3;
  //const GeomPoint<D>	&p1, &p2, &p3;
  Point<D>		pm;
  double		radius, w1,w3;
public:
  ///
  CircleSeg (const GeomPoint<D> & ap1, 
	     const GeomPoint<D> & ap2, 
	     const GeomPoint<D> & ap3);
  ///
  virtual Point<D> GetPoint (double t) const;
  ///
  virtual const GeomPoint<D> & StartPI () const { return p1; }
  ///
  virtual const GeomPoint<D> & EndPI () const { return p3; }
  ///
  virtual void GetCoeff (Vector & coeffs) const;
  ///
  double Radius() const { return radius; }
  ///
  double StartAngle() const { return w1; }
  ///
  double EndAngle() const { return w3; }
  ///
  const Point<D> & MidPoint(void) const {return pm; }

  virtual string GetType(void) const {return "circle";}

  virtual void LineIntersections (const double a, const double b, const double c,
				  Array < Point<D> > & points, const double eps) const;

  virtual double MaxCurvature(void) const {return 1./radius;}
};






/// 
template<int D>
class DiscretePointsSeg : public SplineSeg<D>
{
  Array<Point<D> > pts;
  GeomPoint<D> p1, p2;
public:
  ///
  DiscretePointsSeg (const Array<Point<D> > & apts);
  ///
  virtual ~DiscretePointsSeg ();
  ///
  virtual Point<D> GetPoint (double t) const;
  ///
  virtual const GeomPoint<D> & StartPI () const { return p1; };
  ///
  virtual const GeomPoint<D> & EndPI () const { return p2; }
  ///
  virtual void GetCoeff (Vector & coeffs) const {;}

  virtual double MaxCurvature(void) const {return 1;}
};







// calculates length of spline-curve
template<int D>
double SplineSeg<D> :: Length () const
{
  Point<D> p, pold;

  int i, n = 100;
  double dt = 1.0 / n;

  pold = GetPoint (0);

  double l = 0;
  for (i = 1; i <= n; i++)
    {
      p = GetPoint (i * dt);
      l += Dist (p, pold);
      pold = p;
    }
  return l;
}



// partitionizes spline curve
template<int D>
void SplineSeg<D> :: Partition (double h, double elto0,
				Mesh & mesh, Point3dTree & searchtree, int segnr) const
{
  int i, j;
  double l; // , r1, r2, ra;
  double lold, dt, frac;
  int n = 100;
  Point<D> p, pold, mark, oldmark;
  Array<double> curvepoints;
  double edgelength, edgelengthold;
  l = Length();

  double h1 = min (StartPI().hmax, h/StartPI().refatpoint);
  double h2 = min (EndPI().hmax, h/EndPI().refatpoint);
  double hcurve = min (hmax, h/reffak);

  //  cout << "Partition, l = " << l << ", h = " << h << endl;
  CalcPartition (l, h, h1, h2, hcurve, elto0, curvepoints);
  //  cout << "curvepoints = " << curvepoints << endl;

  dt = 1.0 / n;

  l = 0;
  j = 1;

  pold = GetPoint (0);
  lold = 0;
  oldmark = pold;
  edgelengthold = 0;
  Array<int> locsearch;

  for (i = 1; i <= n; i++)
    {
      p = GetPoint (i*dt);
      l = lold + Dist (p, pold);
      while (j < curvepoints.Size() && (l >= curvepoints[j] || i == n))
	{
	  frac = (curvepoints[j]-lold) / (l-lold);
	  edgelength = i*dt + (frac-1)*dt;
	  // mark = pold + frac * (p-pold);
	  mark = GetPoint (edgelength);
	  
	  // cout << "mark = " << mark << " =?= " << GetPoint (edgelength) << endl;

	  {
	    PointIndex pi1 = -1, pi2 = -1;
	  
	    Point3d mark3(mark(0), mark(1), 0);
	    Point3d oldmark3(oldmark(0), oldmark(1), 0);

	    Vec<3> v (1e-4*h, 1e-4*h, 1e-4*h);
	    searchtree.GetIntersecting (oldmark3 - v, oldmark3 + v, locsearch);

	    for (int k = 0; k < locsearch.Size(); k++)
	      if ( mesh[PointIndex(locsearch[k])].GetLayer() == layer)
		pi1 = locsearch[k];
	    // if (locsearch.Size()) pi1 = locsearch[0];
	      
	    searchtree.GetIntersecting (mark3 - v, mark3 + v, locsearch);
	    for (int k = 0; k < locsearch.Size(); k++)
	      if ( mesh[PointIndex(locsearch[k])].GetLayer() == layer)
		pi2 = locsearch[k];
	    // if (locsearch.Size()) pi2 = locsearch[0];

	    /*	    
	      for (PointIndex pk = PointIndex::BASE; 
	      pk < mesh.GetNP()+PointIndex::BASE; pk++)
	      {
	      if (Dist (mesh[pk], oldmark3) < 1e-4 * h) pi1 = pk;
	      if (Dist (mesh[pk], mark3) < 1e-4 * h) pi2 = pk;
	      }
	    */
	    

	    //	    cout << "pi1 = " << pi1 << endl;
	    //	    cout << "pi2 = " << pi2 << endl;
	    
	    if (pi1 == -1)
	      {
		pi1 = mesh.AddPoint(oldmark3, layer);
		searchtree.Insert (oldmark3, pi1);
	      }
	    if (pi2 == -1)
	      {
		pi2 = mesh.AddPoint(mark3, layer);
		searchtree.Insert (mark3, pi2);
	      }

	    Segment seg;
	    seg.edgenr = segnr;
	    seg.si = bc; // segnr;
	    seg[0] = pi1;
	    seg[1] = pi2;
	    seg.domin = leftdom;
	    seg.domout = rightdom;
	    seg.epgeominfo[0].edgenr = segnr;
	    seg.epgeominfo[0].dist = edgelengthold;
	    seg.epgeominfo[1].edgenr = segnr;
	    seg.epgeominfo[1].dist = edgelength;
	    seg.singedge_left = hpref_left;
	    seg.singedge_right = hpref_right;
	    mesh.AddSegment (seg);
	  }
	
	  oldmark = mark;
	  edgelengthold = edgelength;
	  j++;
	}
    
      pold = p;
      lold = l;
    }
}


template<int D>
void SplineSeg<D> :: GetPoints (int n, Array<Point<D> > & points)
{
  points.SetSize (n);
  if (n >= 2)
    for (int i = 0; i < n; i++)
      points[i] = GetPoint(double(i) / (n-1));
}


template<int D>
void SplineSeg<D> :: PrintCoeff (ostream & ost) const
{
  Vector u(6);

  GetCoeff(u);

  for ( int i=0; i<6; i++)
    ost << u[i] << "  ";
  ost << endl;
}



/* 
   Implementation of line-segment from p1 to p2
*/


template<int D>
LineSeg<D> :: LineSeg (const GeomPoint<D> & ap1, 
		       const GeomPoint<D> & ap2)
  : p1(ap1), p2(ap2)
{
  ;
}


template<int D>
inline Point<D> LineSeg<D> :: GetPoint (double t) const
{
  return p1 + t * (p2 - p1);
}

template<int D>
Vec<D> LineSeg<D> :: GetTangent (const double t) const
{
  return p2-p1;
}

template<int D>
void LineSeg<D> :: GetDerivatives (const double t, 
				   Point<D> & point,
				   Vec<D> & first,
				   Vec<D> & second) const
{
  first = p2 - p1;
  point = p1 + t * first;
  second = 0;
}


template<int D>
double LineSeg<D> :: Length () const
{
  return Dist (p1, p2);
}


template<int D>
void LineSeg<D> :: GetCoeff (Vector & coeffs) const
{
  coeffs.SetSize(6);

  double dx = p2(0) - p1(0);
  double dy = p2(1) - p1(1);

  coeffs[0] = coeffs[1] = coeffs[2] = 0;
  coeffs[3] = -dy;
  coeffs[4] = dx;
  coeffs[5] = -dx * p1(1) + dy * p1(0);
}



template<int D>
void LineSeg<D> :: LineIntersections (const double a, const double b, const double c,
				      Array < Point<D> > & points, const double eps) const
{
  points.SetSize(0);

  double denom = -a*p2(0)+a*p1(0)-b*p2(1)+b*p1(1);
  if(fabs(denom) < 1e-20)
    return;

  double t = (a*p1(0)+b*p1(1)+c)/denom;
  if((t > -eps) && (t <  1.+eps))
    points.Append(GetPoint(t));
}



template<int D>
void LineSeg<D> :: Project (const Point<D> point, Point<D> & point_on_curve, double & t) const
{
  Vec<D> v = p2-p1;
  double l = v.Length();
  v *= 1./l;
  t = (point-p1)*v;

  if(t<0) t = 0;
  if(t>l) t = l;

  point_on_curve = p1+t*v;

  t *= 1./l;
}


template<int D>
void LineSeg<D> :: GetRawData (Array<double> & data) const
{
  data.Append(2);
  for(int i=0; i<D; i++)
    data.Append(p1[i]);
  for(int i=0; i<D; i++)
    data.Append(p2[i]);
}





template<int D>
SplineSeg3<D> :: SplineSeg3 (const GeomPoint<D> & ap1, 
			     const GeomPoint<D> & ap2,
			     const GeomPoint<D> & ap3)
  : p1(ap1), p2(ap2), p3(ap3)
{
  proj_latest_t = 0.5;
}

template<int D>
inline Point<D> SplineSeg3<D> :: GetPoint (double t) const
{
  double x, y, w;
  double b1, b2, b3;

  b1 = (1-t)*(1-t);
  b2 = sqrt(2.0) * t * (1-t);
  b3 = t * t;

  x = p1(0) * b1 + p2(0) * b2 + p3(0) * b3;
  y = p1(1) * b1 + p2(1) * b2 + p3(1) * b3;
  w = b1 + b2 + b3;

  if(D==3)
    {
      double z = p1(2) * b1 + p2(2) * b2 + p3(2) * b3;
      return Point<D> (x/w, y/w, z/w);
    }
  else
    return Point<D> (x/w, y/w);
}




template<int D>
Vec<D> SplineSeg3<D> :: GetTangent (const double t) const
{
  const double b1 = (1.-t)*((sqrt(2.)-2.)*t-sqrt(2.));
  const double b2 = sqrt(2.)*(1.-2.*t);
  const double b3 = t*((sqrt(2.)-2)*t+2.);


  Vec<D> retval;
  for(int i=0; i<D; i++) 
    retval(i) = b1*p1(i) + b2*p2(i) + b3*p3(i);

  return retval;

}


template<int D>
void SplineSeg3<D> :: GetCoeff (Vector & u) const
{
  DenseMatrix a(6, 6);
  DenseMatrix ata(6, 6);
  Vector f(6);

  u.SetSize(6);

  //  ata.SetSymmetric(1);

  double t = 0;
  for (int i = 0; i < 5; i++, t += 0.25)
    {
      Point<D> p = GetPoint (t);
      a(i, 0) = p(0) * p(0);
      a(i, 1) = p(1) * p(1);
      a(i, 2) = p(0) * p(1);
      a(i, 3) = p(0);
      a(i, 4) = p(1);
      a(i, 5) = 1;
    }
  a(5, 0) = 1;

  CalcAtA (a, ata);

  u = 0;
  u(5) = 1;
  a.MultTrans (u, f);
  ata.Solve (f, u);

  // the sign
  Point<D> p0 = GetPoint(0);
  Vec<D> ht = GetTangent(0);
  Vec<2> tang(ht(0), ht(1));

  double gradx = 2.*u(0)*p0(0) + u(2)*p0(1) + u(3);
  double grady = 2.*u(1)*p0(1) + u(2)*p0(0) + u(4);
  Vec<2> gradn (grady, -gradx);
  
  if (tang * gradn < 0) u *= -1;
}

/*
template<int D>
double SplineSeg3<D> :: MaxCurvature(void) const
{
  Vec<D> v1 = p1-p2;
  Vec<D> v2 = p3-p2;
  double l1 = v1.Length();
  double l2 = v2.Length();
  (*testout) << "v1 " << v1 << " v2 " << v2 << endl;

  double cosalpha = v1*v2/(l1*l2);

  (*testout) << "cosalpha " << cosalpha << endl;

  return sqrt(cosalpha + 1.)/(min2(l1,l2)*(1.-cosalpha));
}
*/
  

template<int D>
void SplineSeg3<D> :: LineIntersections (const double a, const double b, const double c,
					 Array < Point<D> > & points, const double eps) const
{
  points.SetSize(0);

  double t;

  const double c1 = a*p1(0) - sqrt(2.)*a*p2(0) + a*p3(0) 
    + b*p1(1) - sqrt(2.)*b*p2(1) + b*p3(1) 
    + (2.-sqrt(2.))*c;
  const double c2 = -2.*a*p1(0) + sqrt(2.)*a*p2(0) -2.*b*p1(1) + sqrt(2.)*b*p2(1) + (sqrt(2.)-2.)*c;
  const double c3 = a*p1(0) + b*p1(1) + c;

  if(fabs(c1) < 1e-20)
    {
      if(fabs(c2) < 1e-20)
	return;

      t = -c3/c2;
      if((t > -eps) && (t < 1.+eps))
	points.Append(GetPoint(t));
      return;
    }

  const double discr = c2*c2-4.*c1*c3;

  if(discr < 0)
    return;

  if(fabs(discr/(c1*c1)) < 1e-14)
    {
      t = -0.5*c2/c1;
      if((t > -eps) && (t < 1.+eps))
	points.Append(GetPoint(t));
      return;
    }

  t = (-c2 + sqrt(discr))/(2.*c1);
  if((t > -eps) && (t < 1.+eps))
    points.Append(GetPoint(t));

  t = (-c2 - sqrt(discr))/(2.*c1);
  if((t > -eps) && (t < 1.+eps))
    points.Append(GetPoint(t));
}


template < int D >
void SplineSeg3<D> :: GetRawData (Array<double> & data) const
{
  data.Append(3);
  for(int i=0; i<D; i++)
    data.Append(p1[i]);
  for(int i=0; i<D; i++)
    data.Append(p2[i]);
  for(int i=0; i<D; i++)
    data.Append(p3[i]);
}


//########################################################################
//		circlesegment

template<int D>
CircleSeg<D> :: CircleSeg (const GeomPoint<D> & ap1, 
			   const GeomPoint<D> & ap2,
			   const GeomPoint<D> & ap3)
  : p1(ap1), p2(ap2), p3(ap3)
{
  Vec<D> v1,v2;
  
  v1 = p1 - p2;
  v2 = p3 - p2;
  
  Point<D> p1t(p1+v1);
  Point<D> p2t(p3+v2);

  // works only in 2D!!!!!!!!!
    
  Line2d g1t,g2t;

  g1t.P1() = Point<2>(p1(0),p1(1));
  g1t.P2() = Point<2>(p1t(0),p1t(1));
  g2t.P1() = Point<2>(p3(0),p3(1));
  g2t.P2() = Point<2>(p2t(0),p2t(1));

  Point<2> mp = CrossPoint (g1t,g2t);

  pm(0) = mp(0); pm(1) = mp(1);
  radius  = Dist(pm,StartPI());
  Vec2d auxv;
  auxv.X() = p1(0)-pm(0); auxv.Y() = p1(1)-pm(1);
  w1      = Angle(auxv);
  auxv.X() = p3(0)-pm(0); auxv.Y() = p3(1)-pm(1);
  w3      = Angle(auxv);
  if ( fabs(w3-w1) > M_PI )
    {  
      if ( w3>M_PI )   w3 -= 2*M_PI;
      if ( w1>M_PI )   w1 -= 2*M_PI;
    }
}
 

template<int D>
Point<D> CircleSeg<D> :: GetPoint (double t) const
{
  if (t >= 1.0)  { return p3; }
     
  double phi = StartAngle() + t*(EndAngle()-StartAngle());
  Vec<D>  tmp(cos(phi),sin(phi));
     
  return pm + Radius()*tmp;
}
  
template<int D>
void CircleSeg<D> :: GetCoeff (Vector & coeff) const
{ 
  coeff[0] = coeff[1] = 1.0;
  coeff[2] = 0.0;
  coeff[3] = -2.0 * pm[0];
  coeff[4] = -2.0 * pm[1];
  coeff[5] = sqr(pm[0]) + sqr(pm[1]) - sqr(Radius());
}

  
template<int D>
void CircleSeg<D> :: LineIntersections (const double a, const double b, const double c,
					Array < Point<D> > & points, const double eps) const
{
  points.SetSize(0);

  double px=0,py=0;

  if(fabs(b) > 1e-20)
    py = -c/b;
  else
    px = -c/a;

  const double c1 = a*a + b*b;
  const double c2 = 2. * ( a*(py-pm(1)) - b*(px-pm(0)));
  const double c3 = pow(px-pm(0),2) + pow(py-pm(1),2) - pow(Radius(),2);
    
  const double discr = c2*c2 - 4*c1*c3;

  if(discr < 0)
    return;

  Array<double> t;

  if(fabs(discr) < 1e-20)
    t.Append(-0.5*c2/c1);
  else
    {
      t.Append((-c2+sqrt(discr))/(2.*c1));
      t.Append((-c2-sqrt(discr))/(2.*c1));
    }

  for(int i=0; i<t.Size(); i++)
    {
      Point<D> p (px-t[i]*b,py+t[i]*a);

      double angle = atan2(p(1),p(0))+M_PI;

      if(angle > StartAngle()-eps && angle < EndAngle()+eps)
	points.Append(p);
    }
}




template<int D>
DiscretePointsSeg<D> ::   DiscretePointsSeg (const Array<Point<D> > & apts)
  : pts (apts)
{ 
  for(int i=0; i<D; i++)
    {
      p1(i) = apts[0](i);
      p2(i) = apts.Last()(i);
    }
  p1.refatpoint = true;
  p2.refatpoint = true;
}


template<int D>
DiscretePointsSeg<D> :: ~DiscretePointsSeg ()
{ ; }

template<int D>
Point<D> DiscretePointsSeg<D> :: GetPoint (double t) const
{
  double t1 = t * (pts.Size()-1);
  int segnr = int(t1);
  if (segnr < 0) segnr = 0;
  if (segnr >= pts.Size()) segnr = pts.Size()-1;

  double rest = t1 - segnr;
    
  return pts[segnr] + rest*Vec<D>(pts[segnr+1]-pts[segnr]);
}

  

typedef GeomPoint<2> GeomPoint2d;
typedef SplineSeg<2> SplineSegment;
typedef LineSeg<2> LineSegment;
typedef SplineSeg3<2> SplineSegment3;
typedef CircleSeg<2> CircleSegment;
typedef DiscretePointsSeg<2> DiscretePointsSegment;


}




#endif
