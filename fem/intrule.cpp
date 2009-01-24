/*********************************************************************/
/* File:   intrule.cpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Mar. 2000                                             */
/*********************************************************************/

/* 
   Finite Element Integration Rules
*/

#include <fem.hpp>
   
namespace ngfem
{
  using namespace ngfem;
  
  /*
  IntegrationPoint :: 
  IntegrationPoint (const double api[3], double aw)
  {
    glob_nr = -1;
    nr = -1;
    pi[0] = api[0];
    pi[1] = api[1];
    pi[2] = api[2];
    weight = aw;
    precomputed_geometry = 0;
  }
  
  IntegrationPoint :: 
  IntegrationPoint (double p1, double p2, double p3, double aw)
  {
    glob_nr = -1;
    nr = -1;
    pi[0] = p1;
    pi[1] = p2;
    pi[2] = p3;
    weight = aw;
    precomputed_geometry = 0;
  }
  
  IntegrationPoint :: 
  IntegrationPoint (const FlatVector<double> & ap, double aw)
  {
    glob_nr = -1;
    nr = -1;
    pi[0] = (ap.Size() >= 1) ? ap(0) : 0;
    pi[1] = (ap.Size() >= 2) ? ap(1) : 0;
    pi[2] = (ap.Size() >= 3) ? ap(2) : 0;
    weight = aw;
    precomputed_geometry = 0;
  }

  IntegrationPoint & IntegrationPoint :: operator=(const IntegrationPoint & aip)
  {
    glob_nr = aip.IPNr();
    nr = aip.Nr();
    pi[0] = aip(0);
    pi[1] = aip(1);
    pi[2] = aip(2);
    weight = aip.Weight();
    precomputed_geometry = aip.precomputed_geometry;
    return *this;
  }
*/
  
  ostream & operator<< (ostream & ost, const IntegrationPoint & ip)
  {
    ost << "IP globnr " << ip.glob_nr << " locnr = " << ip.nr << ": (" 
	<< ip.pi[0] << ", " << ip.pi[1] << ", " << ip.pi[2] 
	<< "), weight = " << ip.weight << endl;
    return ost;
  }
  


  template <int S, int R, typename SCAL>
  SpecificIntegrationPoint<S,R,SCAL> :: 
  SpecificIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    LocalHeap & lh)
    : DimSpecificIntegrationPoint<R,SCAL> (aip, aeltrans)
  {
    this->eltrans.CalcPointJacobian (this->ip, this->point, dxdxi, lh);

    if (S == R)
      {
	det = Det (dxdxi); 
	if(det == 0) 
	  { 
	    cout << " dxdxi " << dxdxi << endl; 
	    cout << " GetJacobieDet is ZERO !!! " << endl; 
	    *testout << " GetJacobieDet is ZERO !!! " << endl; 
            *testout << "ip = " << this -> ip << endl;
            *testout << "point = " << this->point << endl;
            *testout << "dxdxi = " << dxdxi << endl;
	    exit(0); 
	  } 
	if (det < 0 && 0)
	  {
	    cout << "A,det<0" << endl;
	    (*testout) << "A,det<0" << endl;
	  }
	dxidx = Inv (dxdxi);
      }
    else
      {
	if (R == 3)
	  {
            normalvec = Cross (Vec<3,SCAL> (dxdxi.Col(0)),
                               Vec<3,SCAL> (dxdxi.Col(1)));
            det = L2Norm (normalvec);
            normalvec /= det;
            
            /*
            cout << "normalvec 1 = " << normalvec << endl;
            cout << "det 1 = " << det << endl;

	    det = sqrt ( sqr (dxdxi(1,0) * dxdxi(2,1) -
			      dxdxi(2,0) * dxdxi(1,1)) +
			 sqr (dxdxi(0,0) * dxdxi(2,1) -
			      dxdxi(2,0) * dxdxi(0,1)) +
			 sqr (dxdxi(1,0) * dxdxi(0,1) -
			      dxdxi(0,0) * dxdxi(1,1)) );

	    normalvec(0) = dxdxi(1,0) * dxdxi(2,1) - dxdxi(2,0) * dxdxi(1,1);
	    normalvec(1) = dxdxi(2,0) * dxdxi(0,1) - dxdxi(0,0) * dxdxi(2,1);
	    normalvec(2) = dxdxi(0,0) * dxdxi(1,1) - dxdxi(1,0) * dxdxi(0,1);
	    normalvec /= L2Norm (normalvec);

            cout << "normalvec 2 = " << normalvec << endl;
            cout << "det 2 = " << det << endl;
            */
	  }
	else
	  {
	    det = sqrt ( sqr (dxdxi(0,0)) + sqr (dxdxi(1,0)));

	    normalvec(0) = -dxdxi(1,0) / det;
	    normalvec(1) = dxdxi(0,0) / det;
	  }
	
	Mat<S,S> ata, iata;
	
	ata = Trans (dxdxi) * dxdxi;
	iata = Inv (ata);
	dxidx = iata * Trans (dxdxi);
      }
  }


  /*
  template <int S, int R, typename SCAL>
  SpecificIntegrationPoint<S,R,SCAL> :: 
  SpecificIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,			    
			    const Vec<R, SCAL> & ax,
			    const Mat<R, S, SCAL> & adxdxi, 
			    LocalHeap & lh)
    : DimSpecificIntegrationPoint<R,SCAL> (aip, aeltrans)
  {
    this->point = ax;
    dxdxi = adxdxi;
      
    if (S == R)
      {
	det = Det (dxdxi);
	dxidx = Inv (dxdxi);
      }
    else
      {
	if (R == 3)
	  {
            normalvec = Cross (Vec<3> (dxdxi.Col(0)),
                               Vec<3> (dxdxi.Col(1)));
            det = L2Norm (normalvec);
            normalvec /= det;
	  }
	else
	  {
	    det = sqrt ( sqr (dxdxi(0,0)) + sqr (dxdxi(1,0)));

	    normalvec(0) = -dxdxi(1,0) / det;
	    normalvec(1) = dxdxi(0,0) / det;
	  }
	
	Mat<S,S> ata, iata;
	
	ata = Trans (dxdxi) * dxdxi;
	iata = Inv (ata);
	dxidx = iata * Trans (dxdxi);
      }
  }
  */
  

  template <int S, int R, typename SCAL>
  void SpecificIntegrationPoint<S,R,SCAL> :: 
  CalcHesse (Mat<2> & ddx1, Mat<2> & ddx2) const
  {
    double eps = 1e-6;

    Mat<2> jacr, jacl;
    for (int dir = 0; dir < 2; dir++)
      {
	IntegrationPoint ipr = this->ip;
	IntegrationPoint ipl = this->ip;
	ipr(dir) += eps;
	ipl(dir) -= eps;
	this->eltrans.CalcJacobian (ipr, jacr);    
	this->eltrans.CalcJacobian (ipl, jacl);    

	for (int j = 0; j < 2; j++)
	  {
	    ddx1(dir,j) = (jacr(0,j) - jacl(0,j) ) / (2*eps);
	    ddx2(dir,j) = (jacr(1,j) - jacl(1,j) ) / (2*eps);
	  }
      }
  }



  template <int S, int R, typename SCAL>
  void SpecificIntegrationPoint<S,R,SCAL> :: 
  CalcHesse (Mat<3> & ddx1, Mat<3> & ddx2, Mat<3> & ddx3) const
  {
    double eps = 1e-6;

    Mat<3> jacr, jacl;
    for (int dir = 0; dir < 3; dir++)
      {
	IntegrationPoint ipr = this->ip;
	IntegrationPoint ipl = this->ip;
	ipr(dir) += eps;
	ipl(dir) -= eps;
	this->eltrans.CalcJacobian (ipr, jacr);    
	this->eltrans.CalcJacobian (ipl, jacl);    

	for (int j = 0; j < 3; j++)
	  {
	    ddx1(dir,j) = (jacr(0,j) - jacl(0,j) ) / (2*eps);
	    ddx2(dir,j) = (jacr(1,j) - jacl(1,j) ) / (2*eps);
	    ddx3(dir,j) = (jacr(2,j) - jacl(2,j) ) / (2*eps);
	  }
      }
  }

  template <int S, int R, typename SCAL>
  void SpecificIntegrationPoint<S,R,SCAL> :: 
  CalcHesse (Mat<2> & ddx1, Mat<2> & ddx2, Mat<2> & ddx3) const
  {
    this -> CalcHesse(ddx1, ddx2); 
  }

  
  //   template <int DIMS, int DIMR, typename SCAL>
  //   void SpecificIntegrationPoint<DIMS,DIMR,SCAL> :: SetTV ( const Vec<DIMR,SCAL> & vec )
  //   {
  //     for(int i=0; i<DIMR; i++)
  //       dxdxi(i,0) = vec(i);
  //   }




 



  /*
    template <int DIMS, int DIMR, typename SCAL>
    void SpecificIntegrationPoint<DIMS,DIMR,SCAL> :: SetNV ( const Vec<DIMR,SCAL> & vec)
    {
    normalvec = vec;
    }
  */

  template class SpecificIntegrationPoint<1,1>;
  template class SpecificIntegrationPoint<2,2>;
  template class SpecificIntegrationPoint<3,3>;
  template class SpecificIntegrationPoint<1,2>;
  template class SpecificIntegrationPoint<2,3>;
  template class SpecificIntegrationPoint<1,3>;
  
  IntegrationRule :: IntegrationRule ()
  {
    ;
  }
  
  IntegrationRule :: IntegrationRule (int nips)
  {
    ipts.SetAllocSize (nips);
  }
  
  IntegrationRule :: ~IntegrationRule ()
  {
    for (int i = 0; i < ipts.Size(); i++)
      delete ipts[i];
  }
  
  
  void IntegrationRule :: 
  AddIntegrationPoint (IntegrationPoint * ip)
  {
    ipts.Append (ip);
  }

  
  // computes Gaussean integration formula on (0,1) with n points
  // in: Numerical algs in C (or, was it the Fortran book ?)
  void ComputeGaussRule (int n, 
			 Array<double> & xi, Array<double> & wi)
  {
    //  cout << "compute gauss rule, n = " << n << endl;
    xi.SetSize (n);
    wi.SetSize (n);
    
    int m = (n+1)/2;
    double p1, p2, p3;
    double pp, z, z1;
    for (int i = 1; i <= m; i++)
      {
	z = cos ( M_PI * (i - 0.25) / (n + 0.5));
	
	while(1)
	  {
	    p1 = 1;
	    p2 = 0;
	    for (int j = 1; j <= n; j++)
	      {
		p3 = p2;
		p2 = p1;
		p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
	      }
	    // p1 is legendre polynomial
	    
	    pp = n * (z*p1-p2) / (z*z - 1);
	    z1 = z;
	    z = z1-p1/pp;
	    
	    if (fabs (z - z1) < 1e-14) break;
	  }
	
	xi[i-1] = 0.5 * (1 - z);
	xi[n-i] = 0.5 * (1 + z);
	wi[i-1] = wi[n-i] = 1.0 / ( (1  - z * z) * pp * pp);
      }
    //  (*testout) << "Gauss points with n = " << n << ":" << endl;
    //  for (i = 1; i <= n; i++)
    //    (*testout) << xi.Elem(i) << ",  w= " << wi.Elem(i) << endl;
  }





  double  gammln(double xx)
  {
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
			  -1.231739516,0.120858003e-2,-0.536382e-5};
    
    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (int j=0;j<=5;j++) {
      x += 1.0;
      ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
  }



  const double  EPS = 3.0e-14;
  const int  MAXIT = 10;

  void ComputeGaussJacobiRule (int n, 
			       Array<double> & x, 
			       Array<double> & w,
			       double alf, double bet)
  {
    x.SetSize (n);
    w.SetSize (n);

    int i,its,j;
    double alfbet,an,bn,r1,r2,r3;
    double a,b,c,p1,p2,p3,pp,temp,z,z1;
    
    for (i=1;i<=n;i++) {
      if (i == 1) {
	an=alf/n;
	bn=bet/n;
	r1=(1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n);
	r2=1.0+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
	z=1.0-r1/r2;
      } else if (i == 2) {
	r1=(4.1+alf)/((1.0+alf)*(1.0+0.156*alf));
	r2=1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n;
	r3=1.0+0.012*bet*(1.0+0.25*fabs(alf))/n;
	z -= (1.0-z)*r1*r2*r3;
      } else if (i == 3) {
	r1=(1.67+0.28*alf)/(1.0+0.37*alf);
	r2=1.0+0.22*(n-8.0)/n;
	r3=1.0+8.0*bet/((6.28+bet)*n*n);
	z -= (x[0]-z)*r1*r2*r3;
      } else if (i == n-1) {
	r1=(1.0+0.235*bet)/(0.766+0.119*bet);
	r2=1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0)));
	r3=1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
	z += (z-x[n-4])*r1*r2*r3;
      } else if (i == n) {
	r1=(1.0+0.37*bet)/(1.67+0.28*bet);
	r2=1.0/(1.0+0.22*(n-8.0)/n);
	r3=1.0/(1.0+8.0*alf/((6.28+alf)*n*n));
	z += (z-x[n-3])*r1*r2*r3;
      } else {
	z=3.0*x[i-2]-3.0*x[i-3]+x[i-4];
      }
      alfbet=alf+bet;
      for (its=1;its<=MAXIT;its++) {
	temp=2.0+alfbet;
	p1=(alf-bet+temp*z)/2.0;
	p2=1.0;
	for (j=2;j<=n;j++) {
	  p3=p2;
	  p2=p1;
	  temp=2*j+alfbet;
	  a=2*j*(j+alfbet)*(temp-2.0);
	  b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z);
	  c=2.0*(j-1+alf)*(j-1+bet)*temp;
	  p1=(b*p2-c*p3)/a;
	}
	pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z));
	z1=z;
	z=z1-p1/pp;
	if (fabs(z-z1) <= EPS) break;
      }
      if (its > MAXIT) 
	cout << "too many iterations in gaujac";
      x[i-1]=z;


      /*
      // is not double precision accurate !!!
      w[i-1]=exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-
      gammln(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2);    
      */

      /*
	w[i-1]=
	(alf+n-1)! * (beta+n-1) ! /  (  (n)!  (n+alpha+beta)! )
	*temp*pow(2.0,alfbet)/(pp*p2);    
	*/

      if (bet == 0.0)
	{
	  w[i-1] = 1.0 / ( (n+alf) * n )
	    *temp*pow(2.0,alfbet)/(pp*p2);    
	}
      else
	w[i-1]=
	  exp(  gammln(alf+n) + gammln(bet+n) - gammln(n+1.0) - gammln(n+alfbet+1.0))
	  *temp*pow(2.0,alfbet)/(pp*p2);    
    }
    
    for (int i = 0; i < n; i++)
      {
	w[i] *= 0.5 * pow(1-x[i],-alf) * pow(1+x[i], -bet);
	x[i] = 0.5 * (x[i]+1);
      }
  }
    


  
  void ComputeHermiteRule (int n, 
                           Array<double> & x,
                           Array<double> & w)
  {
    const double EPS = 3e-14;
    const double PIM4 = 1.0 / pow (M_PI, 0.25);

    x.SetSize (n);
    w.SetSize (n);

    int its, j, m;
    double p1, p2, p3, pp, z, z1;

    m = (n+1)/2;
    for (int i = 1; i <= m; i++)
      {
        if (i == 1) 
          {
            z = sqrt ( (double) (2*n+1)) - 1.85575*pow( (double)(2*n+1), -0.16667);
          }
        else if (i == 2)
          {
            z -= 1.14*pow( (double)n, 0.426)/z;
          } 
        else if (i == 3)
          {
            z = 1.86*z-0.86*x[0];
          }
        else if (i == 4)
          {
            z = 1.91*z-0.91*x[1];
          } 
        else
          z = 2.0*z-x[i-3];


        for (its = 1; its <= 1000; its++)
          {
            p1 = PIM4;
            p2 = 0;
            for (int j = 1; j <= n; j++)
              {
                p3 = p2;
                p2 = p1;
                p1 = z * sqrt(2.0/j) * p2 - sqrt (((double) (j-1))/j)*p3;
              }
            
            pp = sqrt ( (double)2*n) * p2;
            z1 = z;
            z = z1-p1/pp;
            if (fabs (z-z1) <= EPS) break;

            if (its > 20) cout << "too many steps" << endl;
          }

        x[i-1] = z;
        x[n-i] = -z;
        w[i-1] = 2.0 / (pp*pp);
        w[n-i] = w[i-1];
      }
  }






#ifdef OLD

  TensorProductIntegrationRule ::   
  TensorProductIntegrationRule (const FiniteElement & el, const IntegrationRule & aseg_rule)
    : type(el.ElementType()), seg_rule(aseg_rule) 
  {
    useho = 0;
    nx = seg_rule.GetNIP();


    if (el.ElementType() == ET_TET)
      {
        const H1HighOrderFiniteElement<3> * h1ho = 
          dynamic_cast<const H1HighOrderFiniteElement<3>*> (&el);

	n = nx*nx*nx;

	if (h1ho)
	  {
	    useho = 1;
	    
	    for (int i = 0; i < 4; i++)
	      vnums[i] = h1ho->GetVNums()[i];
	    
	    for (int i = 0; i < 4; i++)
	      sort[i] = i;
	    
	    for (int i = 0; i < 4; i++)
	      for (int j = 0; j < 3; j++)
		if (vnums[sort[j]] < vnums[sort[j+1]])
		  swap (sort[j], sort[j+1]);
	    
	    for (int i = 0; i < 4; i++)
	      isort[sort[i]] = i;
	    
	    const double refpts[][3] =
	      { { 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 },
		{ 0, 0, 0 } };
	    
	    for (int i = 0; i < 4; i++)
	      for (int j = 0; j < 3; j++)
		vertices[i](j) = refpts[sort[i]][j];
	  }
      }

    if (el.ElementType() == ET_TRIG)
      {
        const H1HighOrderFiniteElement<2> * h1ho = 
          dynamic_cast<const H1HighOrderFiniteElement<2>*> (&el);

	n = nx * nx;

	if (h1ho)
	  {
	    useho = 1;
	    
	    for (int i = 0; i < 3; i++)
	      vnums[i] = h1ho->GetVNums()[i];
	    
	    for (int i = 0; i < 3; i++)
	      sort[i] = i;
	
	    for (int i = 0; i < 3; i++)
	      for (int j = 0; j < 2; j++)
		if (vnums[sort[j]] < vnums[sort[j+1]])
		  swap (sort[j], sort[j+1]);
	    
	    for (int i = 0; i < 3; i++)
	      isort[sort[i]] = i;

	    const double refpts[][3] =
	      { { 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 0 } };

	    for (int i = 0; i < 3; i++)
	      for (int j = 0; j < 3; j++)
		vertices[i](j) = refpts[sort[i]][j];
	  }
      }

  }








  TensorProductIntegrationRule ::   
  TensorProductIntegrationRule (ELEMENT_TYPE atype, const int * avnums, const IntegrationRule & aseg_rule)
    : type(atype), seg_rule(aseg_rule) 
  {
    useho = 1;
    nx = seg_rule.GetNIP();

    if (type == ET_TET)
      {
	n = nx * nx * nx;
	useho = 1;

	for (int i = 0; i < 4; i++)
	  vnums[i] = avnums[i];

	for (int i = 0; i < 4; i++)
	  sort[i] = i;
	
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j < 3; j++)
	    if (vnums[sort[j]] < vnums[sort[j+1]])
	      swap (sort[j], sort[j+1]);

	for (int i = 0; i < 4; i++)
	  isort[sort[i]] = i;

	const double refpts[][3] =
	  { { 1, 0, 0 },
	    { 0, 1, 0 },
	    { 0, 0, 1 },
	    { 0, 0, 0 } };

	for (int i = 0; i < 4; i++)
	  for (int j = 0; j < 3; j++)
	    vertices[i](j) = refpts[sort[i]][j];
      }

    if (type == ET_TRIG)
      {
	n = nx * nx;
	useho = 1;

	for (int i = 0; i < 3; i++)
	  vnums[i] = avnums[i];

	for (int i = 0; i < 3; i++)
	  sort[i] = i;
	
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 2; j++)
	    if (vnums[sort[j]] < vnums[sort[j+1]])
	      swap (sort[j], sort[j+1]);

	for (int i = 0; i < 3; i++)
	  isort[sort[i]] = i;

	const double refpts[][3] =
	  { { 1, 0, 0 },
	    { 0, 1, 0 },
	    { 0, 0, 0 } };

	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 3; j++)
	    vertices[i](j) = refpts[sort[i]][j];
      }


  }


  IntegrationPoint 
  TensorProductIntegrationRule :: GetIP(int i) const
  {
    int ix = i % nx;
    i /= nx;
    int iy = i % nx;
    i /= nx;
    int iz = i;

    return GetIP (ix, iy, iz);
  }


  IntegrationPoint 
  TensorProductIntegrationRule :: GetIP(int ix, int iy, int iz) const
  {
    IntegrationPoint ip;
    switch (type)
      {
      case ET_TET:
	{
	  double xi = seg_rule[ix](0);
	  double eta = seg_rule[iy](0);
	  double zeta = seg_rule[iz](0);

	  if (useho)
	    {
	      double lami[4];
	      
	      lami[3] = zeta;
	      lami[2] = eta * (1-zeta);
	      lami[1] = xi * (1-eta)*(1-zeta);
	      lami[0] = 1-lami[1]-lami[2]-lami[3];

	      for (int j = 0; j < 3; j++)
		ip(j) = lami[isort[j]];
	    }
	  else
	    {
	      ip(0) = xi * (1-eta)*(1-zeta);
	      ip(1) = eta * (1-zeta);
	      ip(2) = zeta;
	    }

	  ip.SetWeight(seg_rule[ix].Weight() * (1-eta)*(1-zeta) *
		       seg_rule[iy].Weight() * (1-zeta) *
		       seg_rule[iz].Weight());
	  
	  ip.SetNr (ix + nx * iy + nx*nx * iz);
	  break;
	}


      case ET_TRIG:
	{
	  double xi = seg_rule[ix](0);
	  double eta = seg_rule[iy](0);

	  if (useho)
	    {
	      double lami[3];
	      
	      lami[2] = eta;
	      lami[1] = xi * (1-eta);
	      lami[0] = 1-lami[1]-lami[2];

	      for (int j = 0; j < 2; j++)
		ip(j) = lami[isort[j]];
	    }
	  else
	    {
	      ip(0) = xi * (1-eta);
	      ip(1) = eta;
	    }

	  ip.SetWeight(seg_rule[ix].Weight() * (1-eta) *
		       seg_rule[iy].Weight());

	  break;
	}
      };

    ip.precomputed_geometry = 1;
    return ip;
  }



  void
  TensorProductIntegrationRule :: GetIP_DTet_DHex(int i, Mat<3> & dtet_dhex) const
  {
    int ix = i % nx;
    i /= nx;
    int iy = i % nx;
    i /= nx;
    int iz = i;

    switch (type)
      {
      case ET_TET:
	{
	  double xi = seg_rule[ix](0);
	  double eta = seg_rule[iy](0);
	  double zeta = seg_rule[iz](0);

	  if (useho)
	    {
	      //	      double lami[4];
	      double dlami[4][4];
	      
	      /*
		lami[3] = zeta;
		lami[2] = eta * (1-zeta);
		lami[1] = xi * (1-eta)*(1-zeta);
		lami[0] = 1-lami[1]-lami[2]-lami[3];
	      */

	      dlami[3][0] = 0;
	      dlami[3][1] = 0;
	      dlami[3][2] = 1;
	      
	      dlami[2][0] = 0;
	      dlami[2][1] = 1-zeta;
	      dlami[2][2] = -eta;
	      
	      dlami[1][0] = (1-eta)*(1-zeta);
	      dlami[1][1] = -xi * (1-zeta);
	      dlami[1][2] = -xi * (1-eta);
	      
	      dlami[0][0] = -(1-eta)*(1-zeta);
	      dlami[0][1] = -(1-xi) * (1-zeta);
	      dlami[0][2] = -(1-xi) * (1-eta);
	      
	      
	      for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
		  dtet_dhex(j,k) = dlami[isort[j]][k];
	    }

	  break;
	}
      };
    // return ip;    
  }


  void
  TensorProductIntegrationRule :: GetIP_DTrig_DQuad(int i, Mat<2> & dtrig_dquad) const
  {
    int ix = i % nx;
    i /= nx;
    int iy = i;

    switch (type)
      {
      case ET_TRIG:
	{
	  double xi = seg_rule[ix](0);
	  double eta = seg_rule[iy](0);

	  if (useho)
	    {
	      // 	      double lami[3];
	      // 	      lami[2] = eta;
	      // 	      lami[1] = xi * (1-eta);
	      // 	      lami[0] = 1-lami[1]-lami[2];

	      double dlami[3][2];
	      
	      dlami[2][0] = 0;
	      dlami[2][1] = 1;
	      
	      dlami[1][0] = (1-eta);
	      dlami[1][1] = -xi;
	      
	      dlami[0][0] = -(1-eta);
	      dlami[0][1] = -(1-xi);
	      
	      for (int j = 0; j < 2; j++)
		for (int k = 0; k < 2; k++)
		  dtrig_dquad(j,k) = dlami[isort[j]][k];
	    }

	  break;
	}
      };
  }


  void TensorProductIntegrationRule :: GetWeights (FlatArray<double> weights)
  {
    switch (type)
      {
      case ET_TET:
	{
	  for (int iz = 0, i = 0; iz < nx; iz++)
	    for (int iy = 0; iy < nx; iy++)
	      {
		double eta = seg_rule[iy](0);
		double zeta = seg_rule[iz](0);
		
		double fac_yz = 
		  (1-eta)*(1-zeta) *
		  seg_rule[iy].Weight() * (1-zeta) *
		  seg_rule[iz].Weight();
		
		for (int ix = 0; ix < nx; ix++, i++)
		  weights[i] = seg_rule[ix].Weight() * fac_yz;
	      }
	  break;
	}
      case ET_TRIG:
	{
	  for (int iy = 0, i = 0; iy < nx; iy++)
	    {
	      double eta = seg_rule[iy](0);
	      double fac_y = (1-eta)* seg_rule[iy].Weight();
	      for (int ix = 0; ix < nx; ix++, i++)
		weights[i] = seg_rule[ix].Weight() * fac_y;
	    }
	  break;
	}
      }
  }
#endif

  











  template <int D>
  IntegrationRuleTP<D> :: IntegrationRuleTP (const ElementTransformation & eltrans,
                                             int order, bool compute_mapping, LocalHeap & lh)
  {
    int nip;


    switch (eltrans.GetElementType())
      {
      case ET_TRIG:
        {
          irx = &SelectIntegrationRuleJacobi10 (order);
          iry = &SelectIntegrationRule (ET_SEGM, order);
          
          int sort[3];
          eltrans.GetSort (FlatArray<int> (3, sort) );
          int isort[3];
          for (int i = 0; i < 3; i++) isort[sort[i]] = i;
          
          nip = irx->GetNIP() * iry->GetNIP();
          xi.SetSize(nip);
          weight.SetSize(nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++, ii++)
              {
                double 
                  x = (*irx)[i1](0), 
                  y = (*iry)[i2](0);
                
                double lami[] = { x, 
                                  y * (1-x), 
                                  (1-x) * (1-y) };
                
                xi[ii](0) = lami[isort[0]];
                xi[ii](1) = lami[isort[1]];
                weight[ii] = (*irx)[i1].Weight()*(*iry)[i2].Weight()*(1-x);
              }

          // trig permutation transformation
          double dlamdx[3][2] = 
            { { 1, 0 },
              { 0, 1 },
              { -1, -1 }};
        
          Mat<2> trans2;
          for (int k = 0; k < 2; k++)
            for (int l = 0; l < 2; l++)
              trans2(l,k) = dlamdx[sort[l]][k];
        
          // trig permutation plus quad -> trig mapping
          dxdxi_duffy.SetSize(nip);
          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            {
              double x = (*irx)[i1](0);
              double invx = 1.0 / (1-x);
            
              for (int i2 = 0; i2 < iry->GetNIP(); i2++, ii++)
                {
                  double y = (*iry)[i2](0);
                  
                  Mat<2> trans3, trans, hmat;
                  
                  // quad -> trig transform
                  trans3(0,0) = 1;
                  trans3(0,1) = 0;
                  
                  trans3(1,0) = y*invx;
                  trans3(1,1) = 1*invx;
                  
                  dxdxi_duffy[ii] = trans3 * trans2;
                }
            }
          break;
        }


      case ET_QUAD:
        {
          irx = &SelectIntegrationRule (ET_SEGM, order);
          iry = &SelectIntegrationRule (ET_SEGM, order);

          nip = irx->GetNIP() * iry->GetNIP();

          xi.SetSize(nip);
          weight.SetSize(nip);
          dxdxi_duffy.SetSize(nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++, ii++)
              {
                xi[ii](0) = (*irx)[i1](0);
                xi[ii](1) = (*iry)[i2](0);
                weight[ii] = (*irx)[i1].Weight()*(*iry)[i2].Weight();
              }
        
          Mat<2> id;
          id = 0;
          id(0,0) = id(1,1) = 1;
        
          for (int i = 0; i < nip; i++)
            dxdxi_duffy[i] = id;

          break;
        }


      case ET_TET:
        {
          irx = &SelectIntegrationRuleJacobi20 (order);
          iry = &SelectIntegrationRuleJacobi10 (order);
          irz = &SelectIntegrationRule (ET_SEGM, order);
          
          int sort[4];
          eltrans.GetSort (FlatArray<int> (4, sort) );
          int isort[4];
          for (int i = 0; i < 4; i++) isort[sort[i]] = i;
          
          nip = irx->GetNIP() * iry->GetNIP() * irz->GetNIP();

          xi.SetSize(nip);
          weight.SetSize(nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
                  double 
                    x = (*irx)[i1](0), 
                    y = (*iry)[i2](0), 
                    z = (*irz)[i3](0);

                  double lami[] = { x, 
                                    y * (1-x), 
                                    z * (1-x) * (1-y), 
                                    (1-x)*(1-y)*(1-z) };
                
                  xi[ii](0) = lami[isort[0]];
                  xi[ii](1) = lami[isort[1]];
                  xi[ii](2) = lami[isort[2]];
                  weight[ii] = (*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight() * 
                    sqr(1-x) * (1-y);
                }


          // tet permutation transformation
          double dlamdx[4][3] = 
            { { 1, 0, 0 },
              { 0, 1, 0 },
              { 0, 0, 1 },
              { -1, -1, -1 }};
        
          Mat<3> trans2;
          for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
              trans2(l,k) = dlamdx[sort[l]][k];
        
        
          // tet permutation plus hex -> tet mapping
          dxdxi_duffy.SetSize(nip);
          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            {
              double x = (*irx)[i1](0);
              double invx = 1.0 / (1-x);
            
              for (int i2 = 0; i2 < iry->GetNIP(); i2++)
                {
                  double y = (*iry)[i2](0);
                  double invxy = 1.0 / ( (1-x) * (1-y) );
                
                  for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                    {
                      double z = (*irz)[i3](0);
                    
                      Mat<3> trans3, trans, hmat;
                    
                      // hex -> tet transform
                      trans3(0,0) = 1;
                      trans3(0,1) = 0;
                      trans3(0,2) = 0;
                    
                      trans3(1,0) = y*invx;
                      trans3(1,1) = 1*invx;
                      trans3(1,2) = 0;
                    
                      trans3(2,0) = z*invxy;
                      trans3(2,1) = z*invxy;
                      trans3(2,2) = 1*invxy;
                    
                      dxdxi_duffy[ii] = trans3 * trans2;
                    }
                }
            }
          break;
        }


      case ET_PRISM:
        {
          irx = &SelectIntegrationRule (ET_SEGM, order);
          iry = &SelectIntegrationRule (ET_SEGM, order+1);
          irz = &SelectIntegrationRule (ET_SEGM, order);

          int sort[6], isort[6];

          eltrans.GetSort (FlatArray<int> (6, sort) );
          for (int i = 0; i < 6; i++) isort[sort[i]] = i;


          nip = irx->GetNIP() * iry->GetNIP() * irz->GetNIP();

          xi.SetSize(nip);
          weight.SetSize(nip);
          dxdxi_duffy.SetSize(nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
                  double 
                    z = (*irx)[i1](0), 
                    x = (*iry)[i2](0), 
                    y = (*irz)[i3](0);
 
                  double lami[] = { x, y *(1-x), (1-x)*(1-y) };
                
                  xi[ii](0) = lami[isort[0]];
                  xi[ii](1) = lami[isort[1]];
                  xi[ii](2) = z;
                  weight[ii] = (*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight() * (1-x);
                }


          // trig permutation transformation
          double dlamdx[3][2] = 
            { { 1, 0 },
              { 0, 1 },
              { -1, -1 }};
        
          Mat<3> trans2;
          trans2 = 0.0;
          for (int k = 0; k < 2; k++)
            for (int l = 0; l < 2; l++)
              trans2(l,k) = dlamdx[sort[l]][k];
          trans2(2,2) = 1.0;
        
          // prism permutation plus hex -> tet mapping
          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
                  double x = (*iry)[i2](0);
                  double invx = 1.0 / (1-x);
                  double y = (*irz)[i3](0);
                
                  Mat<3> trans3;
                
                  // hex -> prism transform
                  trans3 = 0.0;
                
                  trans3(1,0) = 1;
                  trans3(2,0) = y*invx;
                  trans3(2,1) = 1*invx;
                
                  trans3(0,2) = 1;

                  dxdxi_duffy[ii] = trans3 * trans2;
                }

          break;
        }

      case ET_HEX:
        {
          irx = &SelectIntegrationRule (ET_SEGM, order);
          iry = &SelectIntegrationRule (ET_SEGM, order);
          irz = &SelectIntegrationRule (ET_SEGM, order);

          nip = irx->GetNIP() * iry->GetNIP() * irz->GetNIP();

          xi.SetSize(nip);
          weight.SetSize(nip);
          dxdxi_duffy.SetSize(nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
                  xi[ii](0) = (*irx)[i1](0);
                  xi[ii](1) = (*iry)[i2](0);
                  xi[ii](2) = (*irz)[i3](0);
                  weight[ii] = (*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight();
                }
        
          Mat<3> id;
          id = 0;
          id(0,0) = id(1,1) = id(2,2) = 1;
        
          for (int i = 0; i < nip; i++)
            dxdxi_duffy[i] = id;

          break;
        }
      }

    if (compute_mapping && !eltrans.Boundary())
      {
        x.SetSize(nip);
        dxdxi.SetSize(nip);
        eltrans.CalcMultiPointJacobian (xi, x, dxdxi, lh);
      }
  }




  template <int D>
  IntegrationRuleTP<D> :: IntegrationRuleTP (ELEMENT_TYPE eltype, FlatArray<int> sort, 
                                             NODE_TYPE nt, int nodenr, int order, LocalHeap & lh)
  {
    int nip;

    static IntegrationRule ir0, ir1;
    if (ir0.GetNIP() == 0)
      {
        ir0.AddIntegrationPoint (new IntegrationPoint (0.0, 0, 0, 1.0));
        ir1.AddIntegrationPoint (new IntegrationPoint (1.0, 0, 0, 1.0));
      }

    switch (eltype)
      {
      case ET_TRIG:
        {
          if (nt == NT_EDGE)
            {
              // int sort[3];
              // eltrans.GetSort (FlatArray<int> (3, sort) );
              int isort[3];
              for (int i = 0; i < 3; i++) isort[sort[i]] = i;


              const EDGE & edge = ElementTopology::GetEdges (ET_TRIG)[nodenr];
              EDGE sedge;
              sedge[0] = isort[edge[0]];
              sedge[1] = isort[edge[1]];
              if (sedge[0] > sedge[1]) swap (sedge[0], sedge[1]);

              irx = 0; iry = 0;
              if (sedge[0] == 1 && sedge[1] == 2)
                {
                  irx = &ir0;
                  iry = &SelectIntegrationRule (ET_SEGM, order);
                }
                  
              if (sedge[0] == 0 && sedge[1] == 1)
                {
                  irx = &SelectIntegrationRule (ET_SEGM, order);
                  iry = &ir1;
                }

              if (sedge[0] == 0 && sedge[1] == 2)
                {
                  irx = &SelectIntegrationRule (ET_SEGM, order);
                  iry = &ir0;
                }

              nip = irx->GetNIP() * iry->GetNIP();
              xi.SetSize(nip);
              weight.SetSize(nip);
              
              for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
                for (int i2 = 0; i2 < iry->GetNIP(); i2++, ii++)
                  {
                    double 
                      x = (*irx)[i1](0), 
                      y = (*iry)[i2](0);
                    
                    double lami[] = { x, 
                                      y * (1-x), 
                                      (1-x) * (1-y) };
                    
                    xi[ii](0) = lami[isort[0]];
                    xi[ii](1) = lami[isort[1]];
                    weight[ii] = (*irx)[i1].Weight()*(*iry)[i2].Weight(); // *(1-x);
                  }
              
              // trig permutation transformation
              double dlamdx[3][2] = 
                { { 1, 0 },
                  { 0, 1 },
                  { -1, -1 }};
              
              Mat<2> trans2;
              for (int k = 0; k < 2; k++)
                for (int l = 0; l < 2; l++)
                  trans2(l,k) = dlamdx[sort[l]][k];
              
              // trig permutation plus quad -> trig mapping
              dxdxi_duffy.SetSize(nip);
              for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
                {
                  double x = (*irx)[i1](0);
                  double invx = 1.0 / (1-x);
                  
                  for (int i2 = 0; i2 < iry->GetNIP(); i2++, ii++)
                    {
                      double y = (*iry)[i2](0);
                      
                      Mat<2> trans3, trans, hmat;
                      
                      // quad -> trig transform
                      trans3(0,0) = 1;
                      trans3(0,1) = 0;
                      
                      trans3(1,0) = y*invx;
                      trans3(1,1) = 1*invx;
                      
                      dxdxi_duffy[ii] = trans3 * trans2;
                    }
                }
              break;
            }
          
        }


      case ET_TET:
        {
          if (nt != NT_FACE) break;
          // int sort[4];
          // eltrans.GetSort (FlatArray<int> (4, sort) );
          int isort[4];
          for (int i = 0; i < 4; i++) isort[sort[i]] = i;
          
          const FACE & face = ElementTopology::GetFaces (ET_TET)[nodenr];
          FACE sface;
          sface[0] = isort[face[0]];
          sface[1] = isort[face[1]];
          sface[2] = isort[face[2]];
          if (sface[0] > sface[1]) swap (sface[0], sface[1]);
          if (sface[1] > sface[2]) swap (sface[1], sface[2]);
          if (sface[0] > sface[1]) swap (sface[0], sface[1]);

          irx = iry = irz = 0;
          int powx = 0, powy = 0;
          if (sface[0] == 1 && sface[1] == 2 && sface[2] == 3)
            {
              powy = 1;
              irx = &ir0;
              iry = &SelectIntegrationRule (ET_SEGM, order+1);
              irz = &SelectIntegrationRule (ET_SEGM, order);
            }
          if (sface[0] == 0 && sface[1] == 2 && sface[2] == 3)
            {
              powx = 1;
              irx = &SelectIntegrationRule (ET_SEGM, order+1);
              iry = &ir0;
              irz = &SelectIntegrationRule (ET_SEGM, order);
            }
          if (sface[0] == 0 && sface[1] == 1 && sface[2] == 3)
            {
              powx = 1;
              irx = &SelectIntegrationRule (ET_SEGM, order+1);
              iry = &SelectIntegrationRule (ET_SEGM, order);
              irz = &ir0;
            }
          if (sface[0] == 0 && sface[1] == 1 && sface[2] == 2)
            {
              powx = 1;
              irx = &SelectIntegrationRule (ET_SEGM, order+1);
              iry = &SelectIntegrationRule (ET_SEGM, order);
              irz = &ir1;
            }

          nip = irx->GetNIP() * iry->GetNIP() * irz->GetNIP();

          xi.SetSize(nip);
          weight.SetSize(nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
                  double 
                    x = (*irx)[i1](0), 
                    y = (*iry)[i2](0), 
                    z = (*irz)[i3](0);

                  double lami[] = { x, 
                                    y * (1-x), 
                                    z * (1-x) * (1-y), 
                                    (1-x)*(1-y)*(1-z) };
                
                  xi[ii](0) = lami[isort[0]];
                  xi[ii](1) = lami[isort[1]];
                  xi[ii](2) = lami[isort[2]];
                  weight[ii] = (*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight() * 
                    pow (1-x, powx) * pow(1-y, powy);
                }


          // tet permutation transformation
          double dlamdx[4][3] = 
            { { 1, 0, 0 },
              { 0, 1, 0 },
              { 0, 0, 1 },
              { -1, -1, -1 }};
        
          Mat<3> trans2;
          for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
              trans2(l,k) = dlamdx[sort[l]][k];
        
        
          // tet permutation plus hex -> tet mapping
          dxdxi_duffy.SetSize(nip);
          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            {
              double x = (*irx)[i1](0);
              double invx = 1.0 / (1-x);
            
              for (int i2 = 0; i2 < iry->GetNIP(); i2++)
                {
                  double y = (*iry)[i2](0);
                  double invxy = 1.0 / ( (1-x) * (1-y) );
                
                  for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                    {
                      double z = (*irz)[i3](0);
                    
                      Mat<3> trans3, trans, hmat;
                    
                      // hex -> tet transform
                      trans3(0,0) = 1;
                      trans3(0,1) = 0;
                      trans3(0,2) = 0;
                    
                      trans3(1,0) = y*invx;
                      trans3(1,1) = 1*invx;
                      trans3(1,2) = 0;
                    
                      trans3(2,0) = z*invxy;
                      trans3(2,1) = z*invxy;
                      trans3(2,2) = 1*invxy;
                    
                      dxdxi_duffy[ii] = trans3 * trans2;
                    }
                }
            }
          break;
        }
      }
  }




  template class IntegrationRuleTP<1>;
  template class IntegrationRuleTP<2>;
  template class IntegrationRuleTP<3>;








  IntegrationRules :: IntegrationRules ()
    {
    int i, p;
    IntegrationRule * rule;    
    IntegrationPoint * ip;


    Array<double> xi;
    Array<double> wi;

    // ************************************
    // ** Segment integration rules
    // ************************************


    for (p = 0; p < 20; p++)
      GenerateIntegrationRule (ET_SEGM, p);

    static double qf_segm_lumping_points[][3] = 
      { 
	{ 0 },
	{ 1 },
      };
    static double qf_segm_lumping_weights[] = 
      { 0.5, 0.5 } ;
    for (i = 0; i < 2; i++)
      {
	ip = new IntegrationPoint (qf_segm_lumping_points[i],
				   qf_segm_lumping_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (segmentpoints.Append (ip)-1);
	segmentlumping.AddIntegrationPoint (ip);
      }


    // ************************************
    // ** Triangle integration rules
    // ************************************

    trigrules.SetSize (7);

    static double qf_trig_order1_points[][3] = 
      {
	{ 1.0/3.0, 1.0/3.0 },
      };
    
    static double qf_trig_order1_weights[] = 
      {  0.5} ;


    
    rule = new IntegrationRule (1);
    ip = new IntegrationPoint (qf_trig_order1_points[0],
			       qf_trig_order1_weights[0]);
    ip->SetNr (0);
    ip->SetGlobNr (trigpoints.Append (ip)-1);
    rule->AddIntegrationPoint (ip);
    trigrules[0] = rule;


    rule = new IntegrationRule (1);
    ip = new IntegrationPoint (qf_trig_order1_points[0],
			       qf_trig_order1_weights[0]);
    ip->SetNr (0);
    ip->SetGlobNr (trigpoints.Append (ip)-1);
    rule->AddIntegrationPoint (ip);
    trigrules[1] = rule;


    static double qf_tria_order2_points[][3] = 
      {
	{ 0,   0.5 },
	{ 0.5, 0,  },
	{ 0.5, 0.5 }
      };
   
    static double qf_tria_order2_weights[] = 
      {
	1.0/6.0, 1.0/6.0 , 1.0/6.0
      };
    
    rule = new IntegrationRule (3);
    for (i = 0; i < 3; i++)
      {
	ip = new IntegrationPoint (qf_tria_order2_points[i],
				   qf_tria_order2_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (trigpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    trigrules[2] = rule;


    static double qf_trig_order4_points[][3] = 
      {
	{ 0.816847572980459, 0.091576213509771, },
	{ 0.091576213509771, 0.816847572980459, },
	{ 0.091576213509771, 0.091576213509771, },
	{ 0.108103018168070, 0.445948490915965, },
	{ 0.445948490915965, 0.108103018168070, },
	{ 0.445948490915965, 0.445948490915965 }
      };
    
    
    static double qf_trig_order4_weights[] = 
      {
	0.054975871827661, 0.054975871827661, 0.054975871827661,
	0.111690794839005, 0.111690794839005, 0.111690794839005
      };

    rule = new IntegrationRule (6);
    for (i = 0; i < 6; i++)
      {
	ip = new IntegrationPoint (qf_trig_order4_points[i],
				   qf_trig_order4_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (trigpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    trigrules[3] = rule;


    rule = new IntegrationRule (6);
    for (i = 0; i < 6; i++)
      {
	ip = new IntegrationPoint (qf_trig_order4_points[i],
				   qf_trig_order4_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (trigpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    trigrules[4] = rule;



    static double qf_trig_order6_points[][3] = 
      {
	{ 0.873821971016996, 0.063089014491502, },
	{ 0.063089014491502, 0.873821971016996, },
	{ 0.063089014491502, 0.063089014491502, },
	{ 0.501426509658179, 0.249286745170910, },
	{ 0.249286745170910, 0.501426509658179, },
	{ 0.249286745170910, 0.249286745170910, },

	{ 0.636502499121399, 0.310352451033785, },
	{ 0.310352451033785, 0.053145049844816, },
	{ 0.053145049844816, 0.636502499121399, },
	{ 0.636502499121399, 0.053145049844816, },
	{ 0.310352451033785, 0.636502499121399, },
	{ 0.053145049844816, 0.310352451033785, }
      };

    static double qf_trig_order6_weights[] = 
      {
	0.025422453185103, 0.025422453185103, 0.025422453185103,
	0.058393137863189, 0.058393137863189, 0.058393137863189,

	0.041425537809187, 0.041425537809187, 0.041425537809187,
	0.041425537809187, 0.041425537809187, 0.041425537809187 
      };
  
    rule = new IntegrationRule (12);
    for (i = 0; i < 12; i++)
      {
	ip = new IntegrationPoint (qf_trig_order6_points[i],
				   qf_trig_order6_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (trigpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    trigrules[5] = rule;

    rule = new IntegrationRule (12);
    for (i = 0; i < 12; i++)
      {
	ip = new IntegrationPoint (qf_trig_order6_points[i],
				   qf_trig_order6_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (trigpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    trigrules[6] = rule;

    for (p = 7; p <= 10; p++)
      GenerateIntegrationRule (ET_TRIG, p);



    static double qf_trig_lumping_points[][3] = 
      {
	{ 1, 0 },
	{ 0, 1, },
	{ 0, 0, }
      };

    static double qf_trig_lumping_weights[] = 
      {
	1.0/6.0, 1.0/6.0 , 1.0/6.0
      };
    
    for (i = 0; i < 3; i++)
      {
	ip = new IntegrationPoint (qf_trig_lumping_points[i],
				   qf_trig_lumping_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (trigpoints.Append (ip)-1);
	triglumping.AddIntegrationPoint (ip);
      }



    static double qf_trig_lumping2_points[][3] = 
      {
	{ 1, 0 },
	{ 0, 1, },
	{ 0, 0, },
	{ 0, 0.5 },
	{ 0.5, 0 },
	{ 0.5, 0.5 }
      };
       
    static double qf_trig_lumping2_weights[] = 
      {
	1.0/12.0, 1.0/12.0 , 1.0/12.0, 
	1.0/12.0, 1.0/12.0 , 1.0/12.0 
      };
    
    for (i = 0; i < 6; i++)
      {
	ip = new IntegrationPoint (qf_trig_lumping2_points[i],
				   qf_trig_lumping2_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (trigpoints.Append (ip)-1);
	triglumping2.AddIntegrationPoint (ip);
      }





    trignodalrules.SetSize(12);
    for (p = 1; p <= trignodalrules.Size(); p++)
      {
	rule = new IntegrationRule ((p+1)*(p+2)/2);
	int nelp = (p*p+3*p+2) / 2;
      
	int lami[3];
	double xi[3];
	int i = 0;
	for (lami[0] = 0; lami[0] <= p; lami[0]++)
	  for (lami[1] = 0; lami[1] <= p-lami[0]; lami[1]++)
	    {
	      lami[2] = p - lami[0] - lami[1];
	      
	      for (int n = 0; n < 3; n++)
		xi[n] = double(lami[n]) / p;
	    
	      ip = new IntegrationPoint (xi, 1.0 / (2.0 * nelp));
	      ip->SetNr (i); i++;
	      ip->SetGlobNr (trigpoints.Append (ip)-1);
	      rule->AddIntegrationPoint (ip);
	    }
	trignodalrules[p-1] = rule;
      }
  



    // ************************************
    // ** Quadrilateral integration rules
    // ************************************
    
    for (p = 0; p <= 10; p++)
      GenerateIntegrationRule (ET_QUAD, p);


    {
      static double qf_quad_lumping_points[][4] = 
	{
	  { 0, 0 },
	  { 1, 0 },
	  { 1, 1 },
	  { 0, 1 }
	};
      
      static double qf_quad_lumping_weights[] = 
	{
	  0.25, 0.25, 0.25, 0.25, 
	};
    
      for (i = 0; i < 4; i++)
	{
	  ip = new IntegrationPoint (qf_quad_lumping_points[i],
				     qf_quad_lumping_weights[i]);
	  ip->SetNr (i);
	  ip->SetGlobNr (quadpoints.Append (ip)-1);
	  quadlumping.AddIntegrationPoint (ip);
	}

    }










    // ************************************
    // ** Tetrahedral integration rules
    // ************************************


    tetrules.SetSize(5);

    static double qf_tetra_order1_points[][3] = 
      { 
	{ 0.25, 0.25, 0.25 },
      };
    
    static double qf_tetra_order1_weights[] = 
      {
	1.0/6.0
      };
    
    rule = new IntegrationRule (1);

    ip = new IntegrationPoint (qf_tetra_order1_points[0],
			       qf_tetra_order1_weights[0]);
    ip->SetNr (0);
    ip->SetGlobNr (tetpoints.Append (ip)-1);
    rule->AddIntegrationPoint (ip);

    tetrules[0] = rule;
    




    static double qf_tetra_order2_points[][3] = 
      {
	{ 0.585410196624969, 0.138196601125011, 0.138196601125011 },
	{ 0.138196601125011, 0.585410196624969, 0.138196601125011 },
	{ 0.138196601125011, 0.138196601125011, 0.585410196624969 },
	{ 0.138196601125011, 0.138196601125011, 0.138196601125011 }
      };
    
    static double qf_tetra_order2_weights[] = 
      { 1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0 };
    
    /*
      rule = new IntegrationRule (4);
      for (i = 0; i < 4; i++)
      {
      ip = new IntegrationPoint (qf_tetra_order2_points[i],
      qf_tetra_order2_weights[i]);
      ip->SetNr (i);
      ip->SetGlobNr (tetpoints.Append (ip)-1);
      rule->AddIntegrationPoint (ip);
      }

      tetrules[1] = rule;
    */
    rule = new IntegrationRule (1);
    for (i = 0; i < 1; i++)
      {
	ip = new IntegrationPoint (qf_tetra_order1_points[i],
				   qf_tetra_order1_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (tetpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    tetrules[1] = rule;




    rule = new IntegrationRule (4);
    for (i = 0; i < 4; i++)
      {
	ip = new IntegrationPoint (qf_tetra_order2_points[i],
				   qf_tetra_order2_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (tetpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    tetrules[2] = rule;




    static double qf_tetra_order5_points[][3] = 
      {
	{ 0.454496295874350,   0.454496295874350,   0.045503704125650 },
	{ 0.454496295874350,   0.045503704125650,   0.454496295874350 },
	{ 0.454496295874350,   0.045503704125650,   0.045503704125650 },
	{ 0.045503704125650,   0.454496295874350,   0.454496295874350 },
	{ 0.045503704125650,   0.454496295874350,   0.045503704125650 },
	{ 0.045503704125650,   0.045503704125650,   0.454496295874350 },
	{ 0.310885919263301,   0.310885919263301,   0.310885919263301 },   
	{ 0.067342242210098,   0.310885919263301,   0.310885919263301 }, 
	{ 0.310885919263301,   0.067342242210098,   0.310885919263301 },  
	{ 0.310885919263301,   0.310885919263301,   0.067342242210098 },
      
	{ 0.092735250310891,   0.092735250310891,   0.092735250310891 },
	{ 0.721794249067326,   0.092735250310891,   0.092735250310891 },
	{ 0.092735250310891,   0.721794249067326,   0.092735250310891 },
	{ 0.092735250310891,   0.092735250310891,   0.721794249067326 }
      };

    static double qf_tetra_order5_weights[] = 
      {
	0.007091003462847, 0.007091003462847, 0.007091003462847,
	0.007091003462847, 0.007091003462847, 0.007091003462847,
	0.018781320953003, 0.018781320953003, 0.018781320953003, 0.018781320953003,
	0.012248840519394, 0.012248840519394, 0.012248840519394, 0.012248840519394
      };

    rule = new IntegrationRule (14);
    for (i = 0; i < 14; i++)
      {
	ip = new IntegrationPoint (qf_tetra_order5_points[i],
				   qf_tetra_order5_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (tetpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    tetrules[3] = rule;

    rule = new IntegrationRule (14);
    for (i = 0; i < 14; i++)
      {
	ip = new IntegrationPoint (qf_tetra_order5_points[i],
				   qf_tetra_order5_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (tetpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    tetrules[4] = rule;

    for (p = 5; p <= 8; p++)
      GenerateIntegrationRule (ET_TET, p);



    tetnodalrules.SetSize(8);
    for (p = 3; p <= tetnodalrules.Size(); p++)
      {
	int nelp = (p*p*p + 6 * p * p + 11 * p + 6) / 6;
	rule = new IntegrationRule (nelp);
      
	int lami[4];
	double xi[4];
	int i = 0;
	for (lami[0] = 0; lami[0] <= p; lami[0]++)
	  for (lami[1] = 0; lami[1] <= p-lami[0]; lami[1]++)
	    for (lami[2] = 0; lami[2] <= p-lami[0]-lami[1]; lami[2]++)
	      {
		lami[3] = p - lami[0] - lami[1] - lami[2];
	      
		int n;
		for (n = 0; n < 4; n++)
		  xi[n] = double(lami[n]) / p;

		ip = new IntegrationPoint (xi, 1.0 / (6.0 * nelp));
		ip->SetNr (i); i++;
		ip->SetGlobNr (tetpoints.Append (ip)-1);
		rule->AddIntegrationPoint (ip);
	      }
	tetnodalrules[p-1] = rule;
      }

    double tet1pts[][3] = 
      { { 1, 0, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1 },
	{ 0, 0, 0 } };
      
    rule = new IntegrationRule (4);
    for (i = 0; i < 4; i++)
      {
	ip = new IntegrationPoint (tet1pts[i], 1.0 / (6.0 * 4));
	ip->SetNr (i);
	ip->SetGlobNr (tetpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    tetnodalrules[0] = rule;

    rule = new IntegrationRule (4);
    for (i = 0; i < 4; i++)
      {
	ip = new IntegrationPoint (tet1pts[i], 1.0 / (6.0 * 4));
	ip->SetNr (i);
	ip->SetGlobNr (tetpoints.Append (ip)-1);
	rule->AddIntegrationPoint (ip);
      }

    tetnodalrules[1] = rule;
  




    static double qf_tet_lumping_points[][3] = 
      {
	{ 0, 0, 0,},
	{ 1, 0, 0,},
	{ 0, 1, 0,},
	{ 0, 0, 1,},
      };
       
    static double qf_tet_lumping_weights[] = 
      {
	1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0, 
      };
    
    for (i = 0; i < 4; i++)
      {
	ip = new IntegrationPoint (qf_tet_lumping_points[i],
				   qf_tet_lumping_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (tetpoints.Append (ip)-1);
	tetlumping.AddIntegrationPoint (ip);
      }



    // ************************************
    // ** Prismatic integration rules
    // ************************************

    for (p = 0; p <= 6; p++)
      GenerateIntegrationRule (ET_PRISM, p);



    static double qf_prismfacemidpoint[][3] = 
      {
	{ 0.5, 0, 0.5 },
	{ 0.5, 0.5, 0.5 },
	{ 0, 0.5, 0.5 },
      };

    for (i = 0; i < 3; i++)
      {
	ip = new IntegrationPoint (qf_prismfacemidpoint[i], 0.5);
	ip->SetNr (i);
	ip->SetGlobNr (prismpoints.Append (ip)-1);
	prismfacemidpoint.AddIntegrationPoint (ip);
      }

    static double qf_prism_lumping_points[][3] = 
      {
	{ 1, 0, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 0 },
	{ 1, 0, 1 },
	{ 0, 1, 1 },
	{ 0, 0, 1 }
      };
       
    static double qf_prism_lumping_weights[] = 
      {
	1.0/12.0, 1.0/12.0, 1.0/12.0, 1.0/12.0, 
	1.0/12.0, 1.0/12.0, 
      };
    
    for (i = 0; i < 6; i++)
      {
	ip = new IntegrationPoint (qf_prism_lumping_points[i],
				   qf_prism_lumping_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (prismpoints.Append (ip)-1);
	prismlumping.AddIntegrationPoint (ip);
      }







    // ************************************
    // ** Pyramid integration rules
    // ************************************

    /*

    static double qf_pyramidz_order1_points_weight[][2] = 
    {
    { 0.75,  1.0/3.0 }
    };

    static double qf_pyramidz_order3_points_weight[][2] = 
    {
    { 0.455848155988775, 0.100785882079825 },
    { 0.877485177344559, 0.232547451253508 }
    };
    static double qf_pyramidz_order5_points_weight[][2] = 
    {
    { 0.294997790111502, 0.029950703008581 },
    { 0.652996233961648, 0.146246269259866 },
    { 0.927005975926850, 0.157136361064887 }
    };

    pyramidrules.SetSize(order3d);
    for (k = 0; k < pyramidrules.Size(); k++)
    {
    const IntegrationRule & quadrule = *quadrules[k];
    double (*zrule) [2];
    int nz;
    switch (k+1)
    {
    case 1: 
    zrule = qf_pyramidz_order1_points_weight; 
    nz = 1;
    break;
    case 2: case 3:
    zrule = qf_pyramidz_order3_points_weight; 
    nz = 2;
    break;
    // case 4: case 5:
    default:
    zrule = qf_pyramidz_order5_points_weight; 
    nz = 3;
    break;
    }

    IntegrationRule * pyramidrule = new IntegrationRule();
	
    double point[3], weight;
	
    for (i = 0; i < quadrule.GetNIP(); i++)
    for (j = 0; j < nz; j++)
    {
    const IntegrationPoint & ipquad = quadrule.GetIP(i);
    point[0] = zrule[j][0] * ipquad.Point()[0];
    point[1] = zrule[j][0] * ipquad.Point()[1];
    point[2] = 1-zrule[j][0];
    weight = zrule[j][1] * ipquad.Weight();
	      
    ip = new IntegrationPoint (point, weight);
    ip->SetGlobNr (pyramidpoints.Append (ip)-1);
    pyramidrule->AddIntegrationPoint (ip);
    }
    pyramidrules[k] = pyramidrule;
    }
    */    


    for (p = 0; p <= 6; p++)
      GenerateIntegrationRule (ET_PYRAMID, p);

    static double qf_pyramid_lumping_points[][3] = 
      {
	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{ 1, 1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1-1e-14 },
      };
       
    static double qf_pyramid_lumping_weights[] = 
      {
	1.0/15.0,
	1.0/15.0,
	1.0/15.0,
	1.0/15.0,
	1.0/15.0,
      };    // not optimal, JS !!!
    
    for (i = 0; i < 5; i++)
      {
	ip = new IntegrationPoint (qf_pyramid_lumping_points[i],
				   qf_pyramid_lumping_weights[i]);
	ip->SetNr (i);
	ip->SetGlobNr (pyramidpoints.Append (ip)-1);
	pyramidlumping.AddIntegrationPoint (ip);
      }


    // ************************************
    // ** Hexaeder integration rules
    // ************************************

    for (p = 0; p <= 6; p++)
      GenerateIntegrationRule (ET_HEX, p);

    /*
      cout << "Check trig intrule:";
      for (int order = 0; order < 10; order++)
      {
      cout << "order = " << order << endl;
      const IntegrationRule & rule = *trigrules[order];
      for (int ox = 0; ox <= 4; ox++)
      for (int oy = 0; ox+oy <= 4; oy++)
      {
      double sum = 0;
      for (int j = 0; j < rule.GetNIP(); j++)
      {
      const IntegrationPoint & ip = rule[j];
      sum += ip.Weight() * pow (ip(0), ox) * pow (ip(1), oy);
      }
      cout << "\\int x^" << ox << " y^" << oy << " = " << sum << endl;
      }
      }
    */      

  }




  IntegrationRules :: ~IntegrationRules ()
  {
    for (int i = 0; i < segmentrules.Size(); i++)
      delete segmentrules[i];

    for (int i = 0; i < trigrules.Size(); i++)
      delete trigrules[i];

    for (int i = 0; i < trignodalrules.Size(); i++)
      delete trignodalrules[i];

    for (int i = 0; i < quadrules.Size(); i++)
      delete quadrules[i];

    for (int i = 0; i < tetrules.Size(); i++)
      delete tetrules[i];

    for (int i = 0; i < tetnodalrules.Size(); i++)
      delete tetnodalrules[i];

    for (int i = 0; i < prismrules.Size(); i++)
      delete prismrules[i];

    for (int i = 0; i < pyramidrules.Size(); i++)
      delete pyramidrules[i];

    for (int i = 0; i < hexrules.Size(); i++)
      delete hexrules[i];
  }


  const Array<IntegrationPoint*> & IntegrationRules ::
  GetIntegrationPoints (ELEMENT_TYPE eltyp) const
  {
    switch (eltyp)
      {
      case ET_SEGM:
	return segmentpoints;
      case ET_TRIG:
	return trigpoints;
      case ET_QUAD:
	return quadpoints;
      case ET_TET:
	return tetpoints;
      case ET_PRISM:
	return prismpoints;
      case ET_PYRAMID:
	return pyramidpoints;
      case ET_HEX:
	return hexpoints;
      }

    stringstream str;
    str<< "no intpoints available for element type " 
       << ElementTopology::GetElementName(eltyp) << endl;
    throw Exception (str.str());
  }



  const IntegrationRule & IntegrationRules :: 
  SelectIntegrationRule (ELEMENT_TYPE eltyp, int order) const
  {
    const Array<IntegrationRule*> * ira;

    switch (eltyp)
      {
      case ET_SEGM:
	ira = &segmentrules; break;
      case ET_TRIG:
	ira = &trigrules; break;
      case ET_QUAD:
	ira = &quadrules; break;
      case ET_TET:
	ira = &tetrules; break;
      case ET_PYRAMID:
	ira = &pyramidrules; break;
      case ET_PRISM:
	ira = &prismrules; break;
      case ET_HEX:
	ira = &hexrules; break;
      default:
	{
	  stringstream str;
	  str << "no integration rules for element " << int(eltyp) << endl;
	  throw Exception (str.str());
	}
      }

    if (order < 0) 
      { order = 0; }

    if (order >= ira->Size() || (*ira)[order] == 0)
      {
        return const_cast<IntegrationRules&> (*this).
	  GenerateIntegrationRule (eltyp, order);
      }

    return *((*ira)[order]);
  }
 

  const IntegrationRule & IntegrationRules :: SelectIntegrationRuleJacobi10 (int order) const
  {
    const Array<IntegrationRule*> * ira;
  
    ira = &jacobirules10; 

    if (order < 0) { order = 0; }

    if (order >= ira->Size() || (*ira)[order] == 0)
      {
	return const_cast<IntegrationRules&> (*this).
	  GenerateIntegrationRuleJacobi10 (order);
      }

    return *((*ira)[order]);
  }
 

  const IntegrationRule & IntegrationRules :: SelectIntegrationRuleJacobi20 (int order) const
  {
    const Array<IntegrationRule*> * ira;
  
    ira = &jacobirules20; 

    if (order < 0) { order = 0; }

    if (order >= ira->Size() || (*ira)[order] == 0)
      {
	return const_cast<IntegrationRules&> (*this).
	  GenerateIntegrationRuleJacobi20 (order);
      }

    return *((*ira)[order]);
  }
 






  const IntegrationRule & IntegrationRules :: 
  GenerateIntegrationRule (ELEMENT_TYPE eltyp, int order)
  {
    Array<IntegrationRule*> * ira;

    if (eltyp == ET_QUAD || eltyp == ET_TRIG)
      {
        GenerateIntegrationRule (ET_SEGM, order);
      }

    if (eltyp == ET_TET || eltyp == ET_PRISM || eltyp == ET_HEX || eltyp == ET_PYRAMID) 
      {        
        GenerateIntegrationRule (ET_SEGM, order);
        GenerateIntegrationRule (ET_SEGM, order+2);
        GenerateIntegrationRule (ET_TRIG, order);
        GenerateIntegrationRule (ET_QUAD, order);
      }

#pragma omp critical(genintrule)
    {
  
      switch (eltyp)
	{
	case ET_SEGM:
	  ira = &segmentrules; break;
	case ET_TRIG:
	  ira = &trigrules; break;
	case ET_QUAD:
	  ira = &quadrules; break;
	case ET_TET:
	  ira = &tetrules; break;
	case ET_PYRAMID:
	  ira = &pyramidrules; break;
	case ET_PRISM:
	  ira = &prismrules; break;
	case ET_HEX:
	  ira = &hexrules; break;
	default:
	  {
	    stringstream str;
	    str << "no integration rules for element " << int(eltyp) << endl;
	    throw Exception (str.str());
	  }
	}

      if (ira -> Size() < order+1)
	{
	  int oldsize = ira -> Size();
	  ira -> SetSize (order+1);
	  for (int i = oldsize; i < order+1; i++)
	    (*ira)[i] = 0;
	}

      if ( (*ira)[order] == 0)
	{
	  switch (eltyp)
	    {
	    case ET_SEGM:
	      {
		Array<double> xi, wi;
		ComputeGaussRule (order/2+1, xi, wi);
		IntegrationRule * rule = new IntegrationRule (xi.Size());
		double xip[3] = { 0, 0, 0 };
		for (int j = 0; j < xi.Size(); j++)
		  {
		    xip[0] = xi[j];
		    IntegrationPoint * ip = new IntegrationPoint (xip, wi[j]);
		    ip->SetNr (j);
		    if (order < 20)
		      ip->SetGlobNr (segmentpoints.Append (ip)-1);
		    rule->AddIntegrationPoint (ip);
		  }
		segmentrules[order] = rule;	      
		break;
	      }

	    case ET_TRIG:
	      {
		Array<double> xx, wx, xy, wy;
		ComputeGaussJacobiRule (order/2+1, xx, wx, 1, 0);
		// ComputeGaussRule (order/2+2, xx, wx);
		ComputeGaussRule (order/2+1, xy, wy);

		IntegrationRule * trigrule = new IntegrationRule(xx.Size()*xy.Size());
	      
		int ii = 0;
		for (int i = 0; i < xx.Size(); i++)
		  for (int j = 0; j < xy.Size(); j++)
		    {
		      IntegrationPoint * ip = 
			new IntegrationPoint (xx[i], xy[j]*(1-xx[i]), 0,
					      wx[i]*wy[j]*(1-xx[i]));
		      ip->SetNr (ii); ii++;
		      if (order <= 10)
			ip->SetGlobNr (trigpoints.Append (ip)-1);
		      trigrule->AddIntegrationPoint (ip);
		    }
		trigrules[order] = trigrule;
		break;


		/*
		  const IntegrationRule & segmrule = SelectIntegrationRule (ET_SEGM, order+1);
		  IntegrationRule * trigrule = new IntegrationRule(segmrule.GetNIP()*segmrule.GetNIP());
	
		  double point[3], weight;
	      
		  int ii = 0;
		  for (int i = 0; i < segmrule.GetNIP(); i++)
		  for (int j = 0; j < segmrule.GetNIP(); j++)
		  {
		  const IntegrationPoint & ipsegmi = segmrule.GetIP(i);
		  const IntegrationPoint & ipsegmj = segmrule.GetIP(j);
		    
		  point[0] = ipsegmi.Point()[0];
		  point[1] = ipsegmj.Point()[0]*(1-point[0]);
		  point[2] = 0;
		    
		  weight = ipsegmi.Weight() *
		  ipsegmj.Weight() * (1-point[0]);
		    
		  IntegrationPoint * ip = new IntegrationPoint (point, weight);
		  ip->SetNr (ii); ii++;
		  if (order <= 10)
		  ip->SetGlobNr (trigpoints.Append (ip)-1);
		  trigrule->AddIntegrationPoint (ip);
		  }
		  trigrules[order] = trigrule;
		  break;
		*/
	      }


	    case ET_QUAD:
	      {
		const IntegrationRule & segmrule = SelectIntegrationRule (ET_SEGM, order);
		IntegrationRule * quadrule = new IntegrationRule(segmrule.GetNIP()*segmrule.GetNIP());
	      
		double point[3], weight;

		int ii = 0;
		for (int i = 0; i < segmrule.GetNIP(); i++)
		  for (int j = 0; j < segmrule.GetNIP(); j++)
		    {
		      const IntegrationPoint & ipsegm1 = segmrule.GetIP(i);
		      const IntegrationPoint & ipsegm2 = segmrule.GetIP(j);
		      point[0] = ipsegm1.Point()[0];
		      point[1] = ipsegm2.Point()[0];
		      point[2] = 0;
		      weight = ipsegm1.Weight() * ipsegm2.Weight();
		    
		      IntegrationPoint * ip = new IntegrationPoint (point, weight);
		      ip->SetNr (ii); ii++;
		      if (order <= 10)
			ip->SetGlobNr (quadpoints.Append (ip)-1);
		      quadrule->AddIntegrationPoint (ip);
		    }
		quadrules[order] = quadrule;
		break;
	      }
  

	    case ET_TET:
	      {
		// tet-rules by degenerated tensor product rules:

		Array<double> xx, wx, xy, wy, xz, wz;
		ComputeGaussRule (order/2+1, xz, wz);
		ComputeGaussJacobiRule (order/2+1, xy, wy, 1, 0);
		ComputeGaussJacobiRule (order/2+1, xx, wx, 2, 0);

		IntegrationRule * tetrule = new IntegrationRule(xx.Size()*xy.Size()*xz.Size());
	      
		int ii = 0;
		for (int i = 0; i < xx.Size(); i++)
		  for (int j = 0; j < xy.Size(); j++)
		    for (int k = 0; k < xz.Size(); k++)
		      {
			IntegrationPoint * ip = 
			  new IntegrationPoint (xx[i], 
						xy[j]*(1-xx[i]),
						xz[k]*(1-xx[i])*(1-xy[j]),
						wx[i]*wy[j]*wz[k]*sqr(1-xx[i])*(1-xy[j]));
			ip->SetNr (ii); ii++;
			if (order <= 6)
			  ip->SetGlobNr (tetpoints.Append (ip)-1);
			tetrule->AddIntegrationPoint (ip);
		      }
		tetrules[order] = tetrule;
		break;
	      }

	    case ET_HEX:
	      {
		const IntegrationRule & segmrule = SelectIntegrationRule (ET_SEGM, order);

		IntegrationRule * hexrule = 
		  new IntegrationRule(segmrule.GetNIP()*segmrule.GetNIP()*segmrule.GetNIP());
	
		double point[3], weight;
		int ii = 0;
		for (int i = 0; i < segmrule.GetNIP(); i++)
		  for (int j = 0; j < segmrule.GetNIP(); j++)
		    for (int l = 0; l < segmrule.GetNIP(); l++)
		      {
			const IntegrationPoint & ipsegm1 = segmrule.GetIP(i);
			const IntegrationPoint & ipsegm2 = segmrule.GetIP(j);
			const IntegrationPoint & ipsegm3 = segmrule.GetIP(l);
		      
			point[0] = ipsegm1.Point()[0];
			point[1] = ipsegm2.Point()[0];
			point[2] = ipsegm3.Point()[0];
			weight = ipsegm1.Weight() * ipsegm2.Weight() * ipsegm3.Weight();
		      
			IntegrationPoint * ip = 
			  new IntegrationPoint (point, weight);
		      
			ip->SetNr (ii); ii++;
			if (order <= 6)
			  ip->SetGlobNr (hexpoints.Append (ip)-1);
			hexrule->AddIntegrationPoint (ip);
		      }
		hexrules[order] = hexrule;
		break;
	      }


	    case ET_PRISM:
	      {
		const IntegrationRule & segmrule = SelectIntegrationRule (ET_SEGM, order);
		const IntegrationRule & trigrule = SelectIntegrationRule (ET_TRIG, order);

		IntegrationRule * prismrule = 
		  new IntegrationRule(segmrule.GetNIP()*trigrule.GetNIP());
      
		double point[3], weight;
	      
		int ii = 0;
		for (int i = 0; i < segmrule.GetNIP(); i++)
		  for (int j = 0; j < trigrule.GetNIP(); j++)
		    {
		      const IntegrationPoint & ipsegm = segmrule.GetIP(i);
		      const IntegrationPoint & iptrig = trigrule.GetIP(j);
		      point[0] = iptrig.Point()[0];
		      point[1] = iptrig.Point()[1];
		      point[2] = ipsegm.Point()[0];
		      weight = iptrig.Weight() * ipsegm.Weight();
	      
		      IntegrationPoint * ip = 
			new IntegrationPoint (point, weight);
		      ip->SetNr (ii); ii++;
		      if (order <= 6)
			ip->SetGlobNr (prismpoints.Append (ip)-1);
		      prismrule->AddIntegrationPoint (ip);
		    }
		prismrules[order] = prismrule;

                /*
		Array<double> xxz, wxz, xy, wy;
                ComputeGaussRule (order/2+1, xxz, wxz);
                ComputeGaussJacobiRule (order/2+1, xy, wy, 1, 0);

		IntegrationRule * prismrule = new IntegrationRule(xxz.Size()*xy.Size()*xxz.Size());
	      
		int ii = 0;
		for (int i = 0; i < xxz.Size(); i++)
		  for (int j = 0; j < xy.Size(); j++)
		    for (int k = 0; k < xxz.Size(); k++)
		      {
			IntegrationPoint * ip = 
			  new IntegrationPoint (xxz[i], 
						xy[j]*(1-xxz[i]),
						xxz[k],
						wxz[i]*wy[j]*wxz[k]*(1-xxz[i]));
			ip->SetNr (ii); ii++;
			if (order <= 6)
			  ip->SetGlobNr (prismpoints.Append (ip)-1);
			prismrule->AddIntegrationPoint (ip);
		      }
		prismrules[order] = prismrule;
                */
		break;
	      }


	    case ET_PYRAMID:
	      {
		const IntegrationRule & quadrule = SelectIntegrationRule (ET_QUAD, order);
		const IntegrationRule & segrule = SelectIntegrationRule (ET_SEGM, order+2);

		IntegrationRule * pyramidrule = 
		  new IntegrationRule(quadrule.GetNIP()*segrule.GetNIP());
	
		double point[3], weight;
	      
		int ii = 0;
		for (int i = 0; i < quadrule.GetNIP(); i++)
		  for (int j = 0; j < segrule.GetNIP(); j++)
		    {
		      const IntegrationPoint & ipquad = quadrule.GetIP(i);
		      const IntegrationPoint & ipseg = segrule.GetIP(j);
		      point[0] = (1-ipseg(0)) * ipquad(0);
		      point[1] = (1-ipseg(0)) * ipquad(1);
		      point[2] = ipseg(0);
		      weight = ipseg.Weight() * sqr (1-ipseg(0)) * ipquad.Weight();
		    
		      IntegrationPoint * ip = new IntegrationPoint (point, weight);
		      ip->SetNr (ii); ii++;
		      if (order <= 6)
			ip->SetGlobNr (pyramidpoints.Append (ip)-1);
		      pyramidrule->AddIntegrationPoint (ip);
		    }
		pyramidrules[order] = pyramidrule;
		break;
	      }
	    }
	}

      if ( (*ira)[order] == 0)
	{
	  stringstream str;
	  str << "could not generate Integration rule of order " << order 
	      << " for element type " 
	      << ElementTopology::GetElementName(eltyp) << endl;
	  throw Exception (str.str());
	}
    }

    return *(*ira)[order];
  }






  const IntegrationRule & IntegrationRules :: GenerateIntegrationRuleJacobi10 (int order)
  {
    Array<IntegrationRule*> * ira;
    ira = &jacobirules10; 

#pragma omp critical(genintrule)
    {
    if (ira -> Size() < order+1)
      {
	int oldsize = ira -> Size();
	ira -> SetSize (order+1);
	for (int i = oldsize; i < order+1; i++)
	  (*ira)[i] = 0;
      }

    if ( (*ira)[order] == 0)
      {
        Array<double> xi, wi;
        // ComputeGaussRule (order/2+1, xi, wi);
        ComputeGaussJacobiRule (order/2+1, xi, wi, 1, 0);
        IntegrationRule * rule = new IntegrationRule (xi.Size());
        double xip[3] = { 0, 0, 0 };
        for (int j = 0; j < xi.Size(); j++)
          {
            xip[0] = xi[j];
            IntegrationPoint * ip = new IntegrationPoint (xip, wi[j]);
            ip->SetNr (j);
            if (order < 20)
              ip->SetGlobNr (segmentpoints.Append (ip)-1);
            rule->AddIntegrationPoint (ip);
          }
        jacobirules10[order] = rule;	      
      }

    if ( (*ira)[order] == 0)
      {
	stringstream str;
	str << "could not generate Jacobi-10 integration rule of order " << order 
	    << " for element type " << endl;
	throw Exception (str.str());
      }
    }
    return *(*ira)[order];
  }






  const IntegrationRule & IntegrationRules :: GenerateIntegrationRuleJacobi20 (int order)
  {
    Array<IntegrationRule*> * ira;
    ira = &jacobirules20; 

#pragma omp critical(genintrule)
    {
    if (ira -> Size() < order+1)
      {
	int oldsize = ira -> Size();
	ira -> SetSize (order+1);
	for (int i = oldsize; i < order+1; i++)
	  (*ira)[i] = 0;
      }

    if ( (*ira)[order] == 0)
      {
        Array<double> xi, wi;
        // ComputeGaussRule (order/2+1, xi, wi);
        ComputeGaussJacobiRule (order/2+1, xi, wi, 2, 0);
        IntegrationRule * rule = new IntegrationRule (xi.Size());
        double xip[3] = { 0, 0, 0 };
        for (int j = 0; j < xi.Size(); j++)
          {
            xip[0] = xi[j];
            IntegrationPoint * ip = new IntegrationPoint (xip, wi[j]);
            ip->SetNr (j);
            if (order < 20)
              ip->SetGlobNr (segmentpoints.Append (ip)-1);
            rule->AddIntegrationPoint (ip);
          }
        jacobirules20[order] = rule;	      
      }

    if ( (*ira)[order] == 0)
      {
	stringstream str;
	str << "could not generate Jacobi-20 integration rule of order " << order 
	    << " for element type " << endl;
	throw Exception (str.str());
      }
    }
    return *(*ira)[order];
  }



















  const IntegrationRule & IntegrationRules :: 
  SelectLumpingIntegrationRule (ELEMENT_TYPE eltyp) const
  {
    const IntegrationRule * ir = &triglumping;

    switch (eltyp)
      {
      case ET_SEGM:
	ir = &segmentlumping; break;
      case ET_TRIG:
	ir = &triglumping; break;
      case ET_QUAD:
	ir = &quadlumping; break;
      case ET_TET:
	ir = &tetlumping; break;
      case ET_PRISM:
	ir = &prismlumping; break;
      case ET_PYRAMID:      
	ir = &pyramidlumping; break;
      default:
	{
	  cout << "no lumping integration rules for element " << int(eltyp) << endl;
	  ir = &SelectIntegrationRule (eltyp, 1);
	  //	ir = &triglumping;
	}
      }

    return *ir;
  }
  



  const IntegrationRule & IntegrationRules :: 
  SelectNodalIntegrationRule (ELEMENT_TYPE eltyp, int order) const
  {
    const IntegrationRule * ir;
    ir = NULL;

    switch (eltyp)
      {
      case ET_TET:
	ir = tetnodalrules[order-1]; break;
      case ET_TRIG:
	if (order == 1)
	  ir = &triglumping;
	else if (order == 2)
	  ir = &triglumping2;
	else ir = trignodalrules[order-1]; break;
	break;
      case ET_PYRAMID:
	ir = &pyramidlumping; break;
      case ET_PRISM:
	ir = &prismlumping; break;
      case ET_QUAD:
	ir = &quadlumping; break;
      default:
	{
	  cout << "no nodal integration rules for element " << int(eltyp) << endl;
	  ir = &triglumping;
	}
      }
    if (!ir)
      {
	cout << "no nodal integration rules for element " << int(eltyp)
	     << ", order " << order << endl;
	ir = &triglumping;
      }
    return *ir;
  }
  


  int Integrator :: common_integration_order = -1;

  const IntegrationRules & GetIntegrationRules ()
  {
    static IntegrationRules intrules;
    return intrules;
  }


  const IntegrationRule & SelectIntegrationRule (ELEMENT_TYPE eltype, int order)
  {
    return GetIntegrationRules ().SelectIntegrationRule (eltype, order);
  }

  const IntegrationRule & SelectIntegrationRuleJacobi20 (int order)
  {
    return GetIntegrationRules ().SelectIntegrationRuleJacobi20 (order);
  }

  const IntegrationRule & SelectIntegrationRuleJacobi10 (int order)
  {
    return GetIntegrationRules ().SelectIntegrationRuleJacobi10 (order);
  }

 
 
  /*
    using namespace std;
    class Initabc
    {
    public:
    Initabc ()
    {
    cout << "intrule" << endl;
    cout << "gammln(1) = " << gammln(1) << ", (2) = " << gammln(2) << endl;
    ofstream out ("gamma.out");
    out.precision(15);
    for (double x = 0; x < 10; x += 0.1)
    {
    out << x << " " << gammln(x) << " " << exp(gammln(x)) << endl;
    }
    }
    };

    Initabc xy;
  */
}
