/*

Spline curve for Mesh generator

*/

#include <mystdlib.h>
#include <linalg.hpp>
#include <gprim.hpp>
#include "spline.hpp"

namespace netgen
{

  // just for testing (JS)
  template <int D>
  void ProjectTrivial (const SplineSeg3<D> & seg, 
                       const Point<D> point, Point<D> & point_on_curve, double & t)
  {
    double mindist = -1;
    for (int i = 0; i <= 1000; i++)
      {
        double ht = double(i)/1000;
        Point<D> p = seg.GetPoint(ht);
        double dist = Dist2 (p, point);
        if (i == 0 || dist < mindist)
          {
            mindist = dist;
            t = ht;
          }
      }
    point_on_curve = seg.GetPoint(t);
  }


  template <> 
  void CircleSeg<3> :: LineIntersections (const double a, const double b, const double c,
					  Array < Point<3> > & points, const double eps) const
  {
    cerr << "CircleSeg<3>::LineIntersections not implemented" << endl;
  }
  
  template <> 
  void CircleSeg<2> :: LineIntersections (const double a, const double b, const double c,
					  Array < Point<2> > & points, const double eps) const
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
	Point<2> p (px-t[i]*b,py+t[i]*a);

	double angle = atan2(p(1),p(0))+M_PI;

	if(angle > StartAngle()-eps && angle < EndAngle()+eps)
	  points.Append(p);
      }
  }




  template<int D>
  SplineSeg3<D> :: SplineSeg3 (const GeomPoint<D> & ap1, 
			       const GeomPoint<D> & ap2,
			       const GeomPoint<D> & ap3)
    : p1(ap1), p2(ap2), p3(ap3)
  {
    weight = Dist (p1, p3) / sqrt (0.5 * (Dist2 (p1, p2) + Dist2 (p2, p3)));
    // weight = sqrt(2);
    // cout << "weight = " << weight << endl;
    proj_latest_t = 0.5;
  }

  template<int D>
  inline Point<D> SplineSeg3<D> :: GetPoint (double t) const
  {
    double x, y, w;
    double b1, b2, b3;

    b1 = (1-t)*(1-t);
    b2 = weight * t * (1-t);
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
    const double b1 = (1.-t)*((weight-2.)*t-weight);
    const double b2 = weight*(1.-2.*t);
    const double b3 = t*((weight-2)*t+2.);


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

  


  template<int D>
  void SplineSeg3<D> :: Project (const Point<D> point, Point<D> & point_on_curve, double & t) const
  {
    double t_old = -1;

    if(proj_latest_t > 0. && proj_latest_t < 1.)
      t = proj_latest_t;
    else
      t = 0.5;
	
    Point<D> phi;
    Vec<D> phip,phipp,phimp;
    
    int i=0;

    while(t > -0.5 && t < 1.5 && i<20 && fabs(t-t_old) > 1e-15 )
      {
        GetDerivatives(t,phi,phip,phipp);
	
        t_old = t;

        phimp = phi-point;

        //t = min2(max2(t-(phip*phimp)/(phipp*phimp + phip*phip),0.),1.);
        t -= (phip*phimp)/(phipp*phimp + phip*phip);

        i++;
      }
    
    //if(i<10 && t > 0. && t < 1.)
    if(i<20 && t > -0.4 && t < 1.4)
      {
        if(t < 0)
          {
            t = 0.;
          }
        if(t > 1)
          {
            t = 1.;
          }

        point_on_curve = SplineSeg3<D>::GetPoint(t);
	
        double dist = Dist(point,point_on_curve);
	
        phi =  SplineSeg3<D> ::GetPoint(0);
        double auxdist = Dist(phi,point);
        if(auxdist < dist)
          {
            t = 0.;
            point_on_curve = phi;
            dist = auxdist;
          }
        phi =  SplineSeg3<D> ::GetPoint(1);
        auxdist = Dist(phi,point);
        if(auxdist < dist)
          {
            t = 1.;
            point_on_curve = phi;
            dist = auxdist;
          }
      }
    else
      {
        double t0 = 0;
        double t1 = 0.5;
        double t2 = 1.;

        double d0,d1,d2;

	
        //(*testout) << "newtonersatz" << endl;
        while(t2-t0 > 1e-8)
          {
	    
            phi =  SplineSeg3<D> ::GetPoint(t0); d0 = Dist(phi,point);
            phi =  SplineSeg3<D> ::GetPoint(t1); d1 = Dist(phi,point);
            phi =  SplineSeg3<D> ::GetPoint(t2); d2 = Dist(phi,point);

            double a = (2.*d0 - 4.*d1 +2.*d2)/pow(t2-t0,2);

            if(a <= 0)
              {
                if(d0 < d2)
                  t2 -= 0.3*(t2-t0);
                else
                  t0 += 0.3*(t2-t0);

                t1 = 0.5*(t2+t0);
              }
            else
              {
                double b = (d1-d0-a*(t1*t1-t0*t0))/(t1-t0);

                double auxt1 = -0.5*b/a;

                if(auxt1 < t0)
                  {
                    t2 -= 0.4*(t2-t0);
                    t0 = max2(0.,t0-0.1*(t2-t0));
                  }
                else if (auxt1 > t2)
                  {
                    t0 += 0.4*(t2-t0);
                    t2 = min2(1.,t2+0.1*(t2-t0));
                  }
                else
                  {
                    t1 = auxt1;
                    auxt1 = 0.25*(t2-t0);
                    t0 = max2(0.,t1-auxt1);
                    t2 = min2(1.,t1+auxt1);
                  }
		
                t1 = 0.5*(t2+t0);
              }  

          }

	
        phi =  SplineSeg3<D> ::GetPoint(t0); d0 = Dist(phi,point);
        phi =  SplineSeg3<D> ::GetPoint(t1); d1 = Dist(phi,point);
        phi =  SplineSeg3<D> ::GetPoint(t2); d2 = Dist(phi,point);

        double mind = d0;
        t = t0;
        if(d1 < mind)
          {
            t = t1;
            mind = d1;
          }
        if(d2 < mind)
          {
            t = t2;
            mind = d2;
          }

        point_on_curve =  SplineSeg3<D> ::GetPoint(t);
      }
    //(*testout) << " latest_t " << proj_latest_t << " t " << t << endl;

    proj_latest_t = t;

    /*
    // test it by trivial sampling
    double ht;
    Point<D> hp;
    ProjectTrivial (*this, point, hp, ht);
    if (fabs (t-ht) > 1e-3)
    {
    // if (Dist2 (point, hp) < Dist2 (point, point_on_curve))
    cout << "project is wrong" << endl;
    cout << "t = " << t << ", ht = " << ht << endl;
    cout << "dist org = " << Dist(point, point_on_curve) << endl;
    cout << "dist trivial = " << Dist(point, hp) << endl;
    }
    */
  }






  template<int D>
  void SplineSeg3<D> :: GetDerivatives (const double t, 
                                        Point<D> & point,
                                        Vec<D> & first,
                                        Vec<D> & second) const
  {
    Vec<D> v1(p1), v2(p2), v3(p3);

    double b1 = (1.-t)*(1.-t);
    double b2 = weight*t*(1.-t);
    double b3 = t*t;
    double w = b1+b2+b3;
    b1 *= 1./w; b2 *= 1./w; b3 *= 1./w;

    double b1p = 2.*(t-1.);
    double b2p = weight*(1.-2.*t);
    double b3p = 2.*t;
    const double wp = b1p+b2p+b3p;
    const double fac1 = wp/w;
    b1p *= 1./w; b2p *= 1./w; b3p *= 1./w;

    const double b1pp = 2.;
    const double b2pp = -2.*weight;
    const double b3pp = 2.;
    const double wpp = b1pp+b2pp+b3pp;
    const double fac2 = (wpp*w-2.*wp*wp)/(w*w);

    for(int i=0; i<D; i++)
      point(i) = b1*p1(i) + b2*p2(i) + b3*p3(i);
    
 
    first = (b1p - b1*fac1) * v1 +
      (b2p - b2*fac1) * v2 +
      (b3p - b3*fac1) * v3;

    second = (b1pp/w - 2*b1p*fac1 - b1*fac2) * v1 +
      (b2pp/w - 2*b2p*fac1 - b2*fac2) * v2 +
      (b3pp/w - 2*b3p*fac1 - b3*fac2) * v3;
  }



  template<>
  double SplineSeg3<2> :: MaxCurvature(void) const
  {
    Vec<2> v1 = p1-p2;
    Vec<2> v2 = p3-p2;
    double l1 = v1.Length();
    double l2 = v2.Length();
        
    double cosalpha = (v1*v2)/(l1*l2);
    
            
    return sqrt(cosalpha + 1.)/(min2(l1,l2)*(1.-cosalpha));
  }

  template<>
  double SplineSeg3<3> :: MaxCurvature(void) const
  {
    Vec<3> v1 = p1-p2;
    Vec<3> v2 = p3-p2;
    double l1 = v1.Length();
    double l2 = v2.Length();
        
    double cosalpha = v1*v2/(l1*l2);
    
        
    return sqrt(cosalpha + 1.)/(min2(l1,l2)*(1.-cosalpha));
  }


  template<int D>
  void SplineSeg3<D> :: LineIntersections (const double a, const double b, const double c,
					   Array < Point<D> > & points, const double eps) const
  {
    points.SetSize(0);

    double t;

    const double c1 = a*p1(0) - weight*a*p2(0) + a*p3(0) 
      + b*p1(1) - weight*b*p2(1) + b*p3(1) 
      + (2.-weight)*c;
    const double c2 = -2.*a*p1(0) + weight*a*p2(0) -2.*b*p1(1) + weight*b*p2(1) + (weight-2.)*c;
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



  template class  SplineSeg3<2>;
  template class  SplineSeg3<3>;








}
