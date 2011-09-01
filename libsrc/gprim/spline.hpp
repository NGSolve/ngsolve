#ifndef FILE_SPLINE_HPP
#define FILE_SPLINE_HPP

/**************************************************************************/
/* File:   spline.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   24. Jul. 96                                                    */
/**************************************************************************/

namespace netgen
{



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
      : Point<D>(ap), refatpoint(aref), hmax(1e99), hpref(ahpref) { ; }
  };




  /// base class for 2d - segment
  template < int D >
  class SplineSeg
  {
  public:
    SplineSeg () { ; }
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

    virtual void GetPoints (int n, Array<Point<D> > & points) const;

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
    double weight;
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

  
    DLL_HEADER virtual void GetDerivatives (const double t, 
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

    DLL_HEADER virtual void LineIntersections (const double a, const double b, const double c,
				    Array < Point<D> > & points, const double eps) const;

    DLL_HEADER virtual double MaxCurvature(void) const;

    DLL_HEADER virtual void Project (const Point<D> point, Point<D> & point_on_curve, double & t) const;

    DLL_HEADER virtual void GetRawData (Array<double> & data) const;
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
    GeomPoint<D> p1n, p2n;
  public:
    ///
    DiscretePointsSeg (const Array<Point<D> > & apts);
    ///
    virtual ~DiscretePointsSeg ();
    ///
    virtual Point<D> GetPoint (double t) const;
    ///
    virtual const GeomPoint<D> & StartPI () const { return p1n; };
    ///
    virtual const GeomPoint<D> & EndPI () const { return p2n; }
    ///
    virtual void GetCoeff (Vector & coeffs) const {;}

    virtual double MaxCurvature(void) const {return 1;}
  };






  // calculates length of spline-curve
  template<int D>
  double SplineSeg<D> :: Length () const
  {
    int n = 100;
    double dt = 1.0 / n;

    Point<D> pold = GetPoint (0);

    double l = 0;
    for (int i = 1; i <= n; i++)
      {
	Point<D> p = GetPoint (i * dt);
	l += Dist (p, pold);
	pold = p;
      }

    return l;
  }


  template<int D>
  void SplineSeg<D> :: GetPoints (int n, Array<Point<D> > & points) const
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
  DiscretePointsSeg<D> ::   DiscretePointsSeg (const Array<Point<D> > & apts)
    : pts (apts)
  { 
    for(int i=0; i<D; i++)
      {
	p1n(i) = apts[0](i);
	p2n(i) = apts.Last()(i);
      }
    p1n.refatpoint = 1;
    p2n.refatpoint = 1;
    p1n.hmax = 1e99;
    p2n.hmax = 1e99;
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







  

  // *************************************
  // Template for B-Splines of order ORDER
  // thx to Gerhard Kitzler
  // *************************************
  
  template<int D, int ORDER>
  class BSplineSeg : public SplineSeg<D>
  {
    Array<Point<D> > pts;
    GeomPoint<D> p1n, p2n;    
    Array<int> ti; 

  public:
    ///
    BSplineSeg (const Array<Point<D> > & apts);
    ///
    virtual ~BSplineSeg();
    ///
    virtual Point<D> GetPoint (double t) const;
    ///
    virtual const GeomPoint<D> & StartPI () const { return p1n; };
    ///
    virtual const GeomPoint<D> & EndPI () const { return p2n; }
    ///
    virtual void GetCoeff (Vector & coeffs) const {;}

    virtual double MaxCurvature(void) const {return 1;}
  };

  // Constructor
  template<int D,int ORDER>
  BSplineSeg<D,ORDER> :: BSplineSeg (const Array<Point<D> > & apts)
    : pts (apts)
  { 
    /*
    for(int i=0; i<D; i++)
      {
	p1n(i) = apts[0](i);
	p2n(i) = apts.Last()(i);
      }
    */
    p1n = apts[0];
    p2n = apts.Last();

    /*
    p1n.refatpoint = 1;
    p2n.refatpoint = 1;
    p1n.hmax = 1e99;
    p2n.hmax = 1e99;
    */

    int m=pts.Size()+ORDER;
    ti.SetSize(m);
    // b.SetSize(m-1);
    ti=0;    
    //    b=0.0;
    for(int i=ORDER;i<m-ORDER+1;i++)
      ti[i]=i-ORDER+1;   
    for(int i=m-ORDER+1;i<m;i++)
      ti[i]=m-2*ORDER+1;
  }
  // Destructor
  template<int D,int ORDER>
  BSplineSeg<D, ORDER> :: ~BSplineSeg ()
  { ; }


  // GetPoint Method...(evaluation of BSpline Curve)
  template<int D,int ORDER>
  Point<D> BSplineSeg<D,ORDER> :: GetPoint (double t_in) const
  {    
    int m=pts.Size()+ORDER;           

    double t = t_in * (m-2*ORDER+1);    

    double b[ORDER];
    
    int interval_nr = int(t)+ORDER-1;    
    if (interval_nr < ORDER-1) interval_nr = ORDER-1;
    if (interval_nr > m-ORDER-1) interval_nr = m-ORDER-1;

    b[ORDER-1] = 1.0;
    
    for(int degree=1;degree<ORDER;degree++)
      for (int k = 0; k <= degree; k++)
	{
	  int j = interval_nr-degree+k;
	  double bnew = 0;

	  if (k != 0) 
	    bnew += (t-ti[j]) / ( ti[j+degree]-ti[j] ) * b[k-degree+ORDER-1];
	  if (k != degree)
	    bnew += (ti[j+degree+1]-t) / ( ti[j+degree+1]-ti[j+1] ) * b[k-degree+ORDER];
	  b[k-degree+ORDER-1] = bnew;
	}

    Point<D> p = 0.0;
    for(int i=0; i < ORDER; i++) 
      p += b[i] * Vec<D> (pts[i+interval_nr-ORDER+1]);
    return p;
  }



}


#endif
