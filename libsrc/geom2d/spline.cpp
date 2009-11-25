/*

Spline curve for Mesh generator

*/

#include <mystdlib.h>
#include <csg.hpp>
#include <linalg.hpp>
#include <meshing.hpp>

namespace netgen
{
#include "spline.hpp"


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
    double b2 = sqrt(2.)*t*(1.-t);
    double b3 = t*t;
    double w = b1+b2+b3;
    b1 *= 1./w; b2 *= 1./w; b3 *= 1./w;

    double b1p = 2.*(t-1.);
    double b2p = sqrt(2.)*(1.-2.*t);
    double b3p = 2.*t;
    const double wp = b1p+b2p+b3p;
    const double fac1 = wp/w;
    b1p *= 1./w; b2p *= 1./w; b3p *= 1./w;

    const double b1pp = 2.;
    const double b2pp = -2.*sqrt(2.);
    const double b3pp = 2.;
    const double wpp = b1pp+b2pp+b3pp;
    const double fac2 = (wpp*w-2.*wp*wp)/(w*w);

    for(int i=0; i<D; i++)
      point(i) = b1*p1(i) + b2*p2(i) + b3*p3(i);
    
 
    first = (b1p - b1*fac1) * v1 +
      (b2p - b2*fac1) * v2 +
      (b3p - b3*fac1) * v3;

    /*
      second = (b1pp/w - b1p*fac1 - b1*fac2) * v1 +
      (b2pp/w - b2p*fac1 - b2*fac2) * v2 +
      (b3pp/w - b3p*fac1 - b3*fac2) * v3;
    */
    // JS: 2 was missing
    second = (b1pp/w - 2*b1p*fac1 - b1*fac2) * v1 +
      (b2pp/w - 2*b2p*fac1 - b2*fac2) * v2 +
      (b3pp/w - 2*b3p*fac1 - b3*fac2) * v3;

  }











  void CalcPartition (double l, double h, double h1, double h2,
		      double hcurve, double elto0, Array<double> & points)
  {
    // cout << "calcpart, h = " << h << ", h1 = " << h1 << ", h2 = " << h2 << ", hcurve = " << hcurve << endl;
    int i, j, n, nel;
    double sum, t, dt, fun, fperel, oldf, f;

    n = 1000;

    points.SetSize (0);

    sum = 0;
    dt = l / n;
    t = 0.5 * dt;
    for (i = 1; i <= n; i++)
      {
	fun = min3 (hcurve, t/elto0 + h1, (l-t)/elto0 + h2);
	sum += dt / fun;
	t += dt;
      }

    nel = int (sum+1);
    fperel = sum / nel;

    points.Append (0);

    i = 1;
    oldf = 0;
    t = 0.5 * dt;
    for (j = 1; j <= n && i < nel; j++)
      {
	fun = min3 (hcurve, t/elto0 + h1, (l-t)/elto0 + h2);

	f = oldf + dt / fun;

	while (f > i * fperel && i < nel)
	  {
	    points.Append ( (l/n) * (j-1 +  (i * fperel - oldf) / (f - oldf)) );
	    i++;
	  }
	oldf = f;
	t += dt;
      }
    points.Append (l);
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




  template class  SplineSeg3<2>;
  template class  SplineSeg3<3>;

}
