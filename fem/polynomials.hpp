#ifndef FILE_POLYNOMIALS_HPP
#define FILE_POLYNOMIALS_HPP
// log of gamma function:

#ifndef WIN32
      double min(double,double);
#endif
inline double my_exp(double x)
{
if(abs(x)<1.0e-12)
	return 1.0;
else
	return exp(x);
}

inline double my_gammln(double xx)
{  
    double tmp, ser;
    tmp  = xx + 4.5;
    tmp -= (xx - 0.5) * log(tmp);

    ser = 1.000000000190015
    + (76.18009172947146     / xx)
    - (86.50532032941677     / (xx + 1.0))
    + (24.01409824083091     / (xx + 2.0))
    - (1.231739572450155     / (xx + 3.0))
    + (0.1208650973866179e-2 / (xx + 4.0))
    - (0.5395239384953e-5    / (xx + 5.0));

    return (log(2.5066282746310005 * ser) - tmp);  
}

// exp( log( gamma(n) ) )	-> gamma function
inline double Fact (int n)
{
	return my_exp(my_gammln(n+1.0));
}
/*
	Laguerre polynomials L_n^alpha
    orthogonal w.r.t \int_0^\infty  e^{-x} x^alpha  u(x) v(x) dx

	Recurrence: L_n^alpha = ( 2n + alpha - 1 - x )/n * L_{n-1}^alpha - (n + alpha - 1)/n L_{n-2}^alpha
    
   H_0 = 1
   H_1 = -x + alpha +1
   H_2 = x^2/2 - ( alpha + 2 ) x + ( alpha + 2 ) ( alpha + 1 ) / 2
   H_3 = -x^3/6 + ( alpha + 3) / 2 x^3 - ( alpha + 2 ) ( alpha + 3 ) / 2 x + ( alpha + 3 ) ( alpha + 2 ) ( alpha + 1 ) / 6

	|| L_n^alpha ||^2 = Gamma(alpha+n+1) / n!    ( = (alpha+n)! / n! for alpha natural)
*/
template <class S, class T>
inline void LaguerrePolynomial (int n, S x, double alpha, T & values)
{
    S p1 = 1.0, p2 = 0.0, p3;
    if (n >= 0) 
        p2 = values[0] = 1.0;
    if (n >= 1)
        p1 = values[1] = -x+alpha+1;
    for (int i  = 1; i < n; i++)
    {
        p3 = p2; p2=p1;
        p1 = 1.0 / ( i+1 ) *( (2*i+alpha+1-x) * p2 - (i+alpha) * p3 );
        values[i+1] = p1;
    }
}

/*

*/
template <class S,class T>
inline void LegendreFunction(int n, S x, T & values)
{
	values = 0.0;
	for(int mm=0;mm<=n;mm++)
	{
		// P_mm^mm to start recurrence
		values(mm,mm) = exp( my_gammln(2*mm+1) - (mm)*log(2) - my_gammln(mm+1)  )*pow((1.0-x*x),mm/2.0);
		if(mm % 2 == 1)
			values(mm,mm)*=-1;
		for(int nn = mm;nn<n;nn++)
		{
			values(mm,nn+1) = x*(2.0*nn+1.0)/(nn-mm+1.0)*values(mm,nn) - (nn+mm)/(nn-mm+1.0)*values(mm,nn-1);			
		}
	}		
}

inline void TestLegendreFuncs()
{
	Matrix<> LegendreFuncs(4,4);
	ofstream out11("legendrefuncs11.txt");
	ofstream out21("legendrefuncs21.txt");
	ofstream out22("legendrefuncs22.txt");
	ofstream out31("legendrefuncs31.txt");
	ofstream out32("legendrefuncs32.txt");
	ofstream out33("legendrefuncs33.txt");
	for(double x=-1;x<=1;x+=0.05)
	{
		LegendreFunction(LegendreFuncs.Height()-1,x,LegendreFuncs);		
		out11 <<x<<" "<< Trans(LegendreFuncs(1,1))<<endl;
		out21 <<x<<" "<< Trans(LegendreFuncs(1,2))<<endl;
		out22 <<x<<" "<< Trans(LegendreFuncs(2,2))<<endl;
		out31 <<x<<" "<< Trans(LegendreFuncs(1,3))<<endl;
		out32 <<x<<" "<< Trans(LegendreFuncs(2,3))<<endl;
		out33 <<x<<" "<< Trans(LegendreFuncs(3,3))<<endl;		
	}
}

template <class S, class T>
inline void SphericalHarmonicReal(int n, S theta, S phi,T & values)
{
	Matrix<> legendrefuncs(n+1,n+1);
	LegendreFunction(n,cos(theta),legendrefuncs);
	values = 0.0;	
	for(int nn=0;nn<=n;nn++)
	{
		double scal = sqrt((2*nn+1)/(4*M_PI));
		values(0,nn) = scal*cos(0*phi)*legendrefuncs(0,nn);
		for(int mm=1;mm<=nn;mm++)
		{
			double scal = sqrt((2*nn+1)/(4*M_PI))*exp(0.5*my_gammln(nn-fabs(mm)+1) - 0.5*my_gammln(nn+fabs(mm)+1));			
			if(mm%2==1)
				scal*=-1;
			values(2*mm-1 ,nn) = scal*cos(mm*phi)*legendrefuncs(mm,nn);
			values(2*mm   ,nn) = scal*sin(mm*phi)*legendrefuncs(mm,nn);
		}
	}	
}

inline void TestSphericalHarmonics()
{
	Matrix<> sphericalharmonics(7,4);
	ofstream out00("SphericalHarmonics00.txt");
	ofstream out10("SphericalHarmonics10.txt");
	ofstream out11("SphericalHarmonics11.txt");
	ofstream out20("SphericalHarmonics20.txt");
	ofstream out21("SphericalHarmonics21.txt");
	ofstream out22("SphericalHarmonics22.txt");
	ofstream out30("SphericalHarmonics30.txt");
	ofstream out31("SphericalHarmonics31.txt");
	ofstream out32("SphericalHarmonics32.txt");
	ofstream out33("SphericalHarmonics33.txt");
	for(double phi=0;phi<=2.0*M_PI;phi+=0.025)
		for(double theta=0;theta<=M_PI;theta+=0.025)
	{
		SphericalHarmonicReal(sphericalharmonics.Width()-1,theta,phi,sphericalharmonics);
		out00 <<theta<<" "<<phi<<" "<< Trans(sphericalharmonics(0,0))<<" "<<0.0<<endl;
		out10 <<theta<<" "<<phi<<" "<< Trans(sphericalharmonics(0,1))<<" "<<0.0<<endl;
		out11 <<theta<<" "<<phi<<" "<< Trans(sphericalharmonics(2*1-1,1))<<" "<<Trans(sphericalharmonics(2*1,1))<<endl;
		out20 <<theta<<" "<<phi<<" "<< Trans(sphericalharmonics(0,2))<<" "<<0.0<<endl;
		out21 <<theta<<" "<<phi<<" "<< Trans(sphericalharmonics(2*1-1,2))<<" "<<Trans(sphericalharmonics(2*1,2))<<endl;
		out22 <<theta<<" "<<phi<<" "<< Trans(sphericalharmonics(2*2-1,2))<<" "<<Trans(sphericalharmonics(2*2,2))<<endl;
		out30 <<theta<<" "<<phi<<" "<< Trans(sphericalharmonics(0,3))<<" "<<0.0<<endl;
		out31 <<theta<<" "<<phi<<" "<< Trans(sphericalharmonics(2*1-1,3))<<" "<<Trans(sphericalharmonics(2*1,3))<<endl;
		out32 <<theta<<" "<<phi<<" "<< Trans(sphericalharmonics(2*2-1,3))<<" "<<Trans(sphericalharmonics(2*2,3))<<endl;
		out33 <<theta<<" "<<phi<<" "<< Trans(sphericalharmonics(2*3-1,3))<<" "<<Trans(sphericalharmonics(2*3,3))<<endl;
	}
}





/*
   Hermite polynomials H_n,
   orthogonal w.r.t. \int_R exp(-x^2) u(x) v(x) dx   
   
   Recurrence: H_n =  2 x H_{n-1} - 2 (n-1)  H_{n-2}

   H_0 = 1
   H_1 = 2*x
   H_2 = 4*x^2 - 2
   H_3 = 8*x^2 - 12 x

   || H_n ||^2 = n! 2^n sqrt(pi)
*/

template <class S, class T>
inline void HermitePolynomial (int n, S x, T & values)
{
  S p1, p2, p3;
  
  p2 = 0;
  if (n >= 0)
    p1 = values[0] = 1.0;
  for (int j=1; j<=n; j++)
    {
      p3 = p2; p2 = p1;
      p1 = 2*x*p2 - 2*(j-1)*p3;
      values[j] = p1;
    }
}




/*
   Scaled hermite polynomials H_n,
   orthonormal w.r.t. \int_R exp(-x^2) u(x) v(x) dx   
   
   Recurrence: H_n =  sqrt( 2 / n ) x H_{n-1} - sqrt( (n-1) / n )  H_{n-2}

   ~H_0 = pi^( - 1 / 4 )
   ~H_1 = 1 / ( sqrt(2 sqrt(pi) ) ) x
   ~H_2 = H_2 / sqrt( 2! 2^2 sqrt(pi) )
   ~H_3 = H_3 / sqrt( 3! 2^3 sqrt(pi) )

   || H_n ||^2 = 1
*/
template <class S, class T>
inline void ScaledHermitePolynomial (int n, S x, T & values)
{
  S p1, p2, p3;
  
  p2 = 0;
  if (n >= 0)
    p1 = values[0] = 1.0 / sqrt(sqrt(M_PI));
  for (int j=1; j<=n; j++)
    {
      p3 = p2; p2 = p1;
      // p1 = 2*x*p2 - 2*(j-1)*p3;
      p1 = sqrt(2.0/j)*x*p2 - sqrt( (j-1.0)/j )*p3;
      values[j] = p1;
    }
}

// ~H_i(x) * exp(-x*x/2)
template <class S, class T>
inline void HermiteFunction (int n, S x, T & values)
{    
    int nmin = int(x*x/2.4-500);    
    if (nmin < 0)
    {
		S p1, p2, p3;	
		p2 = 0;
		if (n >= 0) p1 = values[0] = exp(-x*x/2) / sqrt(sqrt(M_PI));
		for (int j=1; j<=n; j++)
		{
			p3 = p2; p2 = p1;
			p1 = sqrt(2.0/j)*x*p2 - sqrt( (j-1.0)/j )*p3;
			values[j] = p1;
		}
    }
    else
    {
		HermiteFunction (min2 (nmin+1, n), sqrt(0.5) * x, values);
		S starti[2];
		for (int i = nmin; i <= min(nmin+1, n); i++)
		{
			S sum = 0.0;
			for (int j = 0; j <= i; j++)
				sum += sqrt(sqrt(M_PI)) * exp (-i/2.0*log (2.0) + 0.5 * (my_gammln(i+1)-my_gammln(j+1)-my_gammln(i-j+1))) * values[j] * values[i-j];
			starti[i-nmin] = sum;
		}	
		for (int i = 0; i < min(nmin, n+1); i++)
			values[i] = 0.0;
		S p1=0, p2=0, p3=0;
		if (nmin <= n)
			p2 = values[nmin] = starti[0];
		if (nmin+1 <= n)
			p1 = values[nmin+1] = starti[1];

		for (int j=nmin+2; j <= n; j++)
		{
			p3 = p2; p2 = p1;
			p1 = sqrt(2.0/j)*x*p2 - sqrt( (j-1.0)/j )*p3;
			values[j] = p1;
		}
    }
}  // End Hermite Function

inline void ComputeLaguerreRule(int order,double alpha,Array<double> & nodes, Array<double> & weights)
{
      nodes.SetSize(order);
      weights.SetSize(order);
      Matrix<> CM(order),eveci(order);
      Vector<ngbla::complex<double> > lami(order);

      CM=0.0;
      for(int i=0;i<order;i++)
            CM(i,i)=2.0*(i+1.0)+alpha-1.0;
      for(int i=0;i<order-1;i++)      
            CM(i+1,i)=CM(i,i+1)=sqrt((i+1.0)*(i+1.0 + alpha));            
      LapackEigenValues(CM,lami,eveci);
      for(int i=0;i<order;i++)
      {
            weights[i] = exp(my_gammln(alpha+1.0))*sqr(eveci(i,0));
            nodes[i] = lami(i).real();
      }        
}

#endif