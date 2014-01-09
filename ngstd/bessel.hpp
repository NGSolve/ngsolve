/*							j0.c
 *
 *	Bessel function of order zero
 *
 *      thanks to Thorsten Hohage
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j0();
 *
 * y = j0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval the following rational
 * approximation is used:
 *
 *
 *        2         2
 * (w - r  ) (w - r  ) P (w) / Q (w)
 *       1         2    3       8
 *
 *            2
 * where w = x  and the two r's are zeros of the function.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30       10000       4.4e-17     6.3e-18
 *    IEEE      0, 30       60000       4.2e-16     1.1e-16
 *
 */
/*							y0.c
 *
 *	Bessel function of the second kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y0();
 *
 * y = y0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind, of order
 * zero, of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval a rational approximation
 * R(x) is employed to compute
 *   y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.
 * Thus a call to j0() is required.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *  Absolute error, when y0(x) < 1; else relative error:
 *
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        9400       7.0e-17     7.9e-18
 *    IEEE      0, 30       30000       1.3e-15     1.6e-16
 *
 */
/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Note: all coefficients satisfy the relative error criterion
 * except YP, YQ which are designed for absolute error. */

//#include "mconf.h"
#define UNK
#ifdef UNK
static double PP[7] = {
  7.96936729297347051624E-4,
  8.28352392107440799803E-2,
  1.23953371646414299388E0,
  5.44725003058768775090E0,
  8.74716500199817011941E0,
  5.30324038235394892183E0,
  9.99999999999999997821E-1,
};
static double PQ[7] = {
  9.24408810558863637013E-4,
  8.56288474354474431428E-2,
  1.25352743901058953537E0,
  5.47097740330417105182E0,
  8.76190883237069594232E0,
  5.30605288235394617618E0,
  1.00000000000000000218E0,
};
#endif
#ifdef DEC
static unsigned short PP[28] = {
0035520,0164604,0140733,0054470,
0037251,0122605,0115356,0107170,
0040236,0124412,0071500,0056303,
0040656,0047737,0045720,0045263,
0041013,0172143,0045004,0142103,
0040651,0132045,0026241,0026406,
0040200,0000000,0000000,0000000,
};
static unsigned short PQ[28] = {
0035562,0052006,0070034,0134666,
0037257,0057055,0055242,0123424,
0040240,0071626,0046630,0032371,
0040657,0011077,0032013,0012731,
0041014,0030307,0050331,0006414,
0040651,0145457,0065021,0150304,
0040200,0000000,0000000,0000000,
};
#endif
#ifdef IBMPC
static unsigned short PP[28] = {
0x6b27,0x983b,0x1d30,0x3f4a,
0xd1cf,0xb35d,0x34b0,0x3fb5,
0x0b98,0x4e68,0xd521,0x3ff3,
0x0956,0xe97a,0xc9fb,0x4015,
0x9888,0x6940,0x7e8c,0x4021,
0x25a1,0xa594,0x3684,0x4015,
0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short PQ[28] = {
0x9737,0xce03,0x4a80,0x3f4e,
0x54e3,0xab54,0xebc5,0x3fb5,
0x069f,0xc9b3,0x0e72,0x3ff4,
0x62bb,0xe681,0xe247,0x4015,
0x21a1,0xea1b,0x8618,0x4021,
0x3a19,0xed42,0x3965,0x4015,
0x0000,0x0000,0x0000,0x3ff0,
};
#endif
#ifdef MIEEE
static unsigned short PP[28] = {
0x3f4a,0x1d30,0x983b,0x6b27,
0x3fb5,0x34b0,0xb35d,0xd1cf,
0x3ff3,0xd521,0x4e68,0x0b98,
0x4015,0xc9fb,0xe97a,0x0956,
0x4021,0x7e8c,0x6940,0x9888,
0x4015,0x3684,0xa594,0x25a1,
0x3ff0,0x0000,0x0000,0x0000,
};
static unsigned short PQ[28] = {
0x3f4e,0x4a80,0xce03,0x9737,
0x3fb5,0xebc5,0xab54,0x54e3,
0x3ff4,0x0e72,0xc9b3,0x069f,
0x4015,0xe247,0xe681,0x62bb,
0x4021,0x8618,0xea1b,0x21a1,
0x4015,0x3965,0xed42,0x3a19,
0x3ff0,0x0000,0x0000,0x0000,
};
#endif

#ifdef UNK
static double QP[8] = {
-1.13663838898469149931E-2,
-1.28252718670509318512E0,
-1.95539544257735972385E1,
-9.32060152123768231369E1,
-1.77681167980488050595E2,
-1.47077505154951170175E2,
-5.14105326766599330220E1,
-6.05014350600728481186E0,
};
static double QQ[7] = {
/*  1.00000000000000000000E0,*/
  6.43178256118178023184E1,
  8.56430025976980587198E2,
  3.88240183605401609683E3,
  7.24046774195652478189E3,
  5.93072701187316984827E3,
  2.06209331660327847417E3,
  2.42005740240291393179E2,
};
#endif
#ifdef DEC
static unsigned short QP[32] = {
0136472,0035021,0142451,0141115,
0140244,0024731,0150620,0105642,
0141234,0067177,0124161,0060141,
0141672,0064572,0151557,0043036,
0142061,0127141,0003127,0043517,
0142023,0011727,0060271,0144544,
0141515,0122142,0126620,0143150,
0140701,0115306,0106715,0007344,
};
static unsigned short QQ[28] = {
/*0040200,0000000,0000000,0000000,*/
0041600,0121272,0004741,0026544,
0042526,0015605,0105654,0161771,
0043162,0123155,0165644,0062645,
0043342,0041675,0167576,0130756,
0043271,0052720,0165631,0154214,
0043000,0160576,0034614,0172024,
0042162,0000570,0030500,0051235,
};
#endif
#ifdef IBMPC
static unsigned short QP[32] = {
0x384a,0x38a5,0x4742,0xbf87,
0x1174,0x3a32,0x853b,0xbff4,
0x2c0c,0xf50e,0x8dcf,0xc033,
0xe8c4,0x5a6d,0x4d2f,0xc057,
0xe8ea,0x20ca,0x35cc,0xc066,
0x392d,0xec17,0x627a,0xc062,
0x18cd,0x55b2,0xb48c,0xc049,
0xa1dd,0xd1b9,0x3358,0xc018,
};
static unsigned short QQ[28] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x25ac,0x413c,0x1457,0x4050,
0x9c7f,0xb175,0xc370,0x408a,
0x8cb5,0xbd74,0x54cd,0x40ae,
0xd63e,0xbdef,0x4877,0x40bc,
0x3b11,0x1d73,0x2aba,0x40b7,
0x9e82,0xc731,0x1c2f,0x40a0,
0x0a54,0x0628,0x402f,0x406e,
};
#endif
#ifdef MIEEE
static unsigned short QP[32] = {
0xbf87,0x4742,0x38a5,0x384a,
0xbff4,0x853b,0x3a32,0x1174,
0xc033,0x8dcf,0xf50e,0x2c0c,
0xc057,0x4d2f,0x5a6d,0xe8c4,
0xc066,0x35cc,0x20ca,0xe8ea,
0xc062,0x627a,0xec17,0x392d,
0xc049,0xb48c,0x55b2,0x18cd,
0xc018,0x3358,0xd1b9,0xa1dd,
};
static unsigned short QQ[28] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4050,0x1457,0x413c,0x25ac,
0x408a,0xc370,0xb175,0x9c7f,
0x40ae,0x54cd,0xbd74,0x8cb5,
0x40bc,0x4877,0xbdef,0xd63e,
0x40b7,0x2aba,0x1d73,0x3b11,
0x40a0,0x1c2f,0xc731,0x9e82,
0x406e,0x402f,0x0628,0x0a54,
};
#endif


#ifdef UNK
static double YP[8] = {
 1.55924367855235737965E4,
-1.46639295903971606143E7,
 5.43526477051876500413E9,
-9.82136065717911466409E11,
 8.75906394395366999549E13,
-3.46628303384729719441E15,
 4.42733268572569800351E16,
-1.84950800436986690637E16,
};
static double YQ[7] = {
/* 1.00000000000000000000E0,*/
 1.04128353664259848412E3,
 6.26107330137134956842E5,
 2.68919633393814121987E8,
 8.64002487103935000337E10,
 2.02979612750105546709E13,
 3.17157752842975028269E15,
 2.50596256172653059228E17,
};
#endif
#ifdef DEC
static unsigned short YP[32] = {
0043563,0120677,0042264,0046166,
0146137,0140371,0113444,0042260,
0050241,0175707,0100502,0063344,
0152144,0125737,0007265,0164526,
0053637,0051621,0163035,0060546,
0155105,0004416,0107306,0060023,
0056035,0045133,0030132,0000024,
0155603,0065132,0144061,0131732,
};
static unsigned short YQ[28] = {
/*0040200,0000000,0000000,0000000,*/
0042602,0024422,0135557,0162663,
0045030,0155665,0044075,0160135,
0047200,0035432,0105446,0104005,
0051240,0167331,0056063,0022743,
0053223,0127746,0025764,0012160,
0055064,0044206,0177532,0145545,
0056536,0111375,0163715,0127201,
};
#endif
#ifdef IBMPC
static unsigned short YP[32] = {
0x898f,0xe896,0x7437,0x40ce,
0x8896,0x32e4,0xf81f,0xc16b,
0x4cdd,0xf028,0x3f78,0x41f4,
0xbd2b,0xe1d6,0x957b,0xc26c,
0xac2d,0x3cc3,0xea72,0x42d3,
0xcc02,0xd1d8,0xa121,0xc328,
0x4003,0x660b,0xa94b,0x4363,
0x367b,0x5906,0x6d4b,0xc350,
};
static unsigned short YQ[28] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xfcb6,0x576d,0x4522,0x4090,
0xbc0c,0xa907,0x1b76,0x4123,
0xd101,0x5164,0x0763,0x41b0,
0x64bc,0x2b86,0x1ddb,0x4234,
0x828e,0xc57e,0x75fc,0x42b2,
0x596d,0xdfeb,0x8910,0x4326,
0xb5d0,0xbcf9,0xd25f,0x438b,
};
#endif
#ifdef MIEEE
static unsigned short YP[32] = {
0x40ce,0x7437,0xe896,0x898f,
0xc16b,0xf81f,0x32e4,0x8896,
0x41f4,0x3f78,0xf028,0x4cdd,
0xc26c,0x957b,0xe1d6,0xbd2b,
0x42d3,0xea72,0x3cc3,0xac2d,
0xc328,0xa121,0xd1d8,0xcc02,
0x4363,0xa94b,0x660b,0x4003,
0xc350,0x6d4b,0x5906,0x367b,
};
static unsigned short YQ[28] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4090,0x4522,0x576d,0xfcb6,
0x4123,0x1b76,0xa907,0xbc0c,
0x41b0,0x0763,0x5164,0xd101,
0x4234,0x1ddb,0x2b86,0x64bc,
0x42b2,0x75fc,0xc57e,0x828e,
0x4326,0x8910,0xdfeb,0x596d,
0x438b,0xd25f,0xbcf9,0xb5d0,
};
#endif

#ifdef UNK
/*  5.783185962946784521175995758455807035071 */
static double DR1 = 5.78318596294678452118E0;
/* 30.47126234366208639907816317502275584842 */
static double DR2 = 3.04712623436620863991E1;
#endif

#ifdef DEC
static unsigned short R1[] = {0040671,0007734,0001061,0056734};
#define DR1 *(double *)R1
static unsigned short R2[] = {0041363,0142445,0030416,0165567};
#define DR2 *(double *)R2
#endif

#ifdef IBMPC
static unsigned short R1[] = {0x2bbb,0x8046,0x21fb,0x4017};
#define DR1 *(double *)R1
static unsigned short R2[] = {0xdd6f,0xa621,0x78a4,0x403e};
#define DR2 *(double *)R2
#endif

#ifdef MIEEE
static unsigned short R1[] = {0x4017,0x21fb,0x8046,0x2bbb};
#define DR1 *(double *)R1
static unsigned short R2[] = {0x403e,0x78a4,0xa621,0xdd6f};
#define DR2 *(double *)R2
#endif

#ifdef UNK
static double RP[4] = {
-4.79443220978201773821E9,
 1.95617491946556577543E12,
-2.49248344360967716204E14,
 9.70862251047306323952E15,
};
static double RQ[8] = {
/* 1.00000000000000000000E0,*/
 4.99563147152651017219E2,
 1.73785401676374683123E5,
 4.84409658339962045305E7,
 1.11855537045356834862E10,
 2.11277520115489217587E12,
 3.10518229857422583814E14,
 3.18121955943204943306E16,
 1.71086294081043136091E18,
};
#endif
#ifdef DEC
static unsigned short RP[16] = {
0150216,0161235,0064344,0014450,
0052343,0135216,0035624,0144153,
0154142,0130247,0003310,0003667,
0055411,0173703,0047772,0176635,
};
static unsigned short RQ[32] = {
/*0040200,0000000,0000000,0000000,*/
0042371,0144025,0032265,0136137,
0044451,0133131,0132420,0151466,
0046470,0144641,0072540,0030636,
0050446,0126600,0045042,0044243,
0052365,0172633,0110301,0071063,
0054215,0032424,0062272,0043513,
0055742,0005013,0171731,0072335,
0057275,0170646,0036663,0013134,
};
#endif
#ifdef IBMPC
static unsigned short RP[16] = {
0x8325,0xad1c,0xdc53,0xc1f1,
0x990d,0xc772,0x7751,0x427c,
0x00f7,0xe0d9,0x5614,0xc2ec,
0x5fb4,0x69ff,0x3ef8,0x4341,
};
static unsigned short RQ[32] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xb78c,0xa696,0x3902,0x407f,
0x1a67,0x36a2,0x36cb,0x4105,
0x0634,0x2eac,0x1934,0x4187,
0x4914,0x0944,0xd5b0,0x4204,
0x2e46,0x7218,0xbeb3,0x427e,
0x48e9,0x8c97,0xa6a2,0x42f1,
0x2e9c,0x7e7b,0x4141,0x435c,
0x62cc,0xc7b6,0xbe34,0x43b7,
};
#endif
#ifdef MIEEE
static unsigned short RP[16] = {
0xc1f1,0xdc53,0xad1c,0x8325,
0x427c,0x7751,0xc772,0x990d,
0xc2ec,0x5614,0xe0d9,0x00f7,
0x4341,0x3ef8,0x69ff,0x5fb4,
};
static unsigned short RQ[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x407f,0x3902,0xa696,0xb78c,
0x4105,0x36cb,0x36a2,0x1a67,
0x4187,0x1934,0x2eac,0x0634,
0x4204,0xd5b0,0x0944,0x4914,
0x427e,0xbeb3,0x7218,0x2e46,
0x42f1,0xa6a2,0x8c97,0x48e9,
0x435c,0x4141,0x7e7b,0x2e9c,
0x43b7,0xbe34,0xc7b6,0x62cc,
};
#endif

#define ANSIPROT
#ifndef ANSIPROT
double bessj0(), polevl(), p1evl(), log(), sin(), cos(), sqrt();
#endif

double polevl(double x,double coef[],int N)
{
double ans;
int i;
double *p;

p = coef;
ans = *p++;
i = N;

do
        ans = ans * x  +  *p++;
while( --i );

return( ans );
}

/*                                                      p1evl() */
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

double p1evl(double x, double coef[],int N)
{
double ans;
double *p;
int i;

p = coef;
ans = x + *p++;
i = N-1;

do
        ans = ans * x  + *p++;
while( --i );

return( ans );
}

double TWOOPI =  6.36619772367581343075535E-1; /* 2/pi */
double THPIO4 =  2.35619449019234492885;       /* 3*pi/4 */
double SQ2OPI =  7.9788456080286535587989E-1;  /* sqrt( 2/pi ) */
double PIO4   =  7.85398163397448309616E-1;    /* pi/4 */

double bessj0(double x)
{
double w, z, p, q, xn;

if( x < 0 )
	x = -x;

if( x <= 5.0 )
	{
	z = x * x;
	if( x < 1.0e-5 )
		return( 1.0 - z/4.0 );

	p = (z - DR1) * (z - DR2);
	p = p * polevl( z, RP, 3)/p1evl( z, RQ, 8 );
	return( p );
	}

w = 5.0/x;
q = 25.0/(x*x);
p = polevl( q, PP, 6)/polevl( q, PQ, 6 );
q = polevl( q, QP, 7)/p1evl( q, QQ, 7 );
xn = x - PIO4;
p = p * cos(xn) - w * q * sin(xn);
return( p * SQ2OPI / sqrt(x) );
}
/*							bessy0() 2	*/
/* Bessel function of second kind, order zero	*/

/* Rational approximation coefficients YP[], YQ[] are used here.
 * The function computed is  bessy0(x)  -  2 * log(x) * bessj0(x) / PI,
 * whose value at x = 0 is  2 * ( log(0.5) + EUL ) / PI
 * = 0.073804295108687225.
 */


double bessy0(double x)
{
double w, z, p, q, xn;

if( x <= 5.0 )
	{
	if( x <= 0.0 )
	  throw Exception ("arg<=0 in bessy0");
	z = x * x;
	w = polevl( z, YP, 7) / p1evl( z, YQ, 7 );
	w += TWOOPI * log(x) * bessj0(x);
	return( w );
	}

w = 5.0/x;
z = 25.0 / (x * x);
p = polevl( z, PP, 6)/polevl( z, PQ, 6 );
q = polevl( z, QP, 7)/p1evl( z, QQ, 7 );
xn = x - PIO4;
p = p * sin(xn) + w * q * cos(xn);
return( p * SQ2OPI / sqrt(x) );
}

/*							j1.c
 *
 *	Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j1();
 *
 * y = j1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 24 term Chebyshev
 * expansion is used. In the second, the asymptotic
 * trigonometric representation is employed using two
 * rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       4.0e-17     1.1e-17
 *    IEEE      0, 30       30000       2.6e-16     1.1e-16
 *
 *
 */
/*							y1.c
 *
 *	Bessel function of second kind of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y1();
 *
 * y = y1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind of order one
 * of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 25 term Chebyshev
 * expansion is used, and a call to j1() is required.
 * In the second, the asymptotic trigonometric representation
 * is employed using two rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       8.6e-17     1.3e-17
 *    IEEE      0, 30       30000       1.0e-15     1.3e-16
 *
 * (error criterion relative when |y1| > 1).
 *
 */
/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#define PIO4 .78539816339744830962
#define THPIO4 2.35619449019234492885
#define SQ2OPI .79788456080286535588

// #include "mconf.h"

#ifdef UNK
static double RP1[4] = {
-8.99971225705559398224E8,
 4.52228297998194034323E11,
-7.27494245221818276015E13,
 3.68295732863852883286E15,
};
static double RQ1[8] = {
/* 1.00000000000000000000E0,*/
 6.20836478118054335476E2,
 2.56987256757748830383E5,
 8.35146791431949253037E7,
 2.21511595479792499675E10,
 4.74914122079991414898E12,
 7.84369607876235854894E14,
 8.95222336184627338078E16,
 5.32278620332680085395E18,
};
#endif
#ifdef DEC
static unsigned short RP1[16] = {
0147526,0110742,0063322,0077052,
0051722,0112720,0065034,0061530,
0153604,0052227,0033147,0105650,
0055121,0055025,0032276,0022015,
};
static unsigned short RQ1[32] = {
/*0040200,0000000,0000000,0000000,*/
0042433,0032610,0155604,0033473,
0044572,0173320,0067270,0006616,
0046637,0045246,0162225,0006606,
0050645,0004773,0157577,0053004,
0052612,0033734,0001667,0176501,
0054462,0054121,0173147,0121367,
0056237,0002777,0121451,0176007,
0057623,0136253,0131601,0044710,
};
#endif
#ifdef IBMPC
static unsigned short RP1[16] = {
0x4fc5,0x4cda,0xd23c,0xc1ca,
0x8c6b,0x0d43,0x52ba,0x425a,
0xf175,0xe6cc,0x8a92,0xc2d0,
0xc482,0xa697,0x2b42,0x432a,
};
static unsigned short RQ1[32] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x86e7,0x1b70,0x66b1,0x4083,
0x01b2,0x0dd7,0x5eda,0x410f,
0xa1b1,0xdc92,0xe954,0x4193,
0xeac1,0x7bef,0xa13f,0x4214,
0xffa8,0x8076,0x46fb,0x4291,
0xf45f,0x3ecc,0x4b0a,0x4306,
0x3f81,0xf465,0xe0bf,0x4373,
0x2939,0x7670,0x7795,0x43d2,
};
#endif
#ifdef MIEEE
static unsigned short RP1[16] = {
0xc1ca,0xd23c,0x4cda,0x4fc5,
0x425a,0x52ba,0x0d43,0x8c6b,
0xc2d0,0x8a92,0xe6cc,0xf175,
0x432a,0x2b42,0xa697,0xc482,
};
static unsigned short RQ1[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4083,0x66b1,0x1b70,0x86e7,
0x410f,0x5eda,0x0dd7,0x01b2,
0x4193,0xe954,0xdc92,0xa1b1,
0x4214,0xa13f,0x7bef,0xeac1,
0x4291,0x46fb,0x8076,0xffa8,
0x4306,0x4b0a,0x3ecc,0xf45f,
0x4373,0xe0bf,0xf465,0x3f81,
0x43d2,0x7795,0x7670,0x2939,
};
#endif

#ifdef UNK
static double PP1[7] = {
 7.62125616208173112003E-4,
 7.31397056940917570436E-2,
 1.12719608129684925192E0,
 5.11207951146807644818E0,
 8.42404590141772420927E0,
 5.21451598682361504063E0,
 1.00000000000000000254E0,
};
static double PQ1[7] = {
 5.71323128072548699714E-4,
 6.88455908754495404082E-2,
 1.10514232634061696926E0,
 5.07386386128601488557E0,
 8.39985554327604159757E0,
 5.20982848682361821619E0,
 9.99999999999999997461E-1,
};
#endif
#ifdef DEC
static unsigned short PP1[28] = {
0035507,0144542,0061543,0024326,
0037225,0145105,0017766,0022661,
0040220,0043766,0010254,0133255,
0040643,0113047,0142611,0151521,
0041006,0144344,0055351,0074261,
0040646,0156520,0120574,0006416,
0040200,0000000,0000000,0000000,
};
static unsigned short PQ1[28] = {
0035425,0142330,0115041,0165514,
0037214,0177352,0145105,0052026,
0040215,0072515,0141207,0073255,
0040642,0056427,0137222,0106405,
0041006,0062716,0166427,0165450,
0040646,0133352,0035425,0123304,
0040200,0000000,0000000,0000000,
};
#endif
#ifdef IBMPC
static unsigned short PP1[28] = {
0x651b,0x4c6c,0xf92c,0x3f48,
0xc4b6,0xa3fe,0xb948,0x3fb2,
0x96d6,0xc215,0x08fe,0x3ff2,
0x3a6a,0xf8b1,0x72c4,0x4014,
0x2f16,0x8b5d,0xd91c,0x4020,
0x81a2,0x142f,0xdbaa,0x4014,
0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short PQ1[28] = {
0x3d69,0x1344,0xb89b,0x3f42,
0xaa83,0x5948,0x9fdd,0x3fb1,
0xeed6,0xb850,0xaea9,0x3ff1,
0x51a1,0xf7d2,0x4ba2,0x4014,
0xfd65,0xdda2,0xccb9,0x4020,
0xb4d9,0x4762,0xd6dd,0x4014,
0x0000,0x0000,0x0000,0x3ff0,
};
#endif
#ifdef MIEEE
static unsigned short PP1[28] = {
0x3f48,0xf92c,0x4c6c,0x651b,
0x3fb2,0xb948,0xa3fe,0xc4b6,
0x3ff2,0x08fe,0xc215,0x96d6,
0x4014,0x72c4,0xf8b1,0x3a6a,
0x4020,0xd91c,0x8b5d,0x2f16,
0x4014,0xdbaa,0x142f,0x81a2,
0x3ff0,0x0000,0x0000,0x0000,
};
static unsigned short PQ1[28] = {
0x3f42,0xb89b,0x1344,0x3d69,
0x3fb1,0x9fdd,0x5948,0xaa83,
0x3ff1,0xaea9,0xb850,0xeed6,
0x4014,0x4ba2,0xf7d2,0x51a1,
0x4020,0xccb9,0xdda2,0xfd65,
0x4014,0xd6dd,0x4762,0xb4d9,
0x3ff0,0x0000,0x0000,0x0000,
};
#endif

#ifdef UNK
static double QP1[8] = {
 5.10862594750176621635E-2,
 4.98213872951233449420E0,
 7.58238284132545283818E1,
 3.66779609360150777800E2,
 7.10856304998926107277E2,
 5.97489612400613639965E2,
 2.11688757100572135698E2,
 2.52070205858023719784E1,
};
static double QQ1[7] = {
/* 1.00000000000000000000E0,*/
 7.42373277035675149943E1,
 1.05644886038262816351E3,
 4.98641058337653607651E3,
 9.56231892404756170795E3,
 7.99704160447350683650E3,
 2.82619278517639096600E3,
 3.36093607810698293419E2,
};
#endif
#ifdef DEC
static unsigned short QP1[32] = {
0037121,0037723,0055605,0151004,
0040637,0066656,0031554,0077264,
0041627,0122714,0153170,0161466,
0042267,0061712,0036520,0140145,
0042461,0133315,0131573,0071176,
0042425,0057525,0147500,0013201,
0042123,0130122,0061245,0154131,
0041311,0123772,0064254,0172650,
};
static unsigned short QQ1[28] = {
/*0040200,0000000,0000000,0000000,*/
0041624,0074603,0002112,0101670,
0042604,0007135,0010162,0175565,
0043233,0151510,0157757,0172010,
0043425,0064506,0112006,0104276,
0043371,0164125,0032271,0164242,
0043060,0121425,0122750,0136013,
0042250,0005773,0053472,0146267,
};
#endif
#ifdef IBMPC
static unsigned short QP1[32] = {
0xba40,0x6b70,0x27fa,0x3faa,
0x8fd6,0xc66d,0xedb5,0x4013,
0x1c67,0x9acf,0xf4b9,0x4052,
0x180d,0x47aa,0xec79,0x4076,
0x6e50,0xb66f,0x36d9,0x4086,
0x02d0,0xb9e8,0xabea,0x4082,
0xbb0b,0x4c54,0x760a,0x406a,
0x9eb5,0x4d15,0x34ff,0x4039,
};
static unsigned short QQ1[28] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x5077,0x6089,0x8f30,0x4052,
0x5f6f,0xa20e,0x81cb,0x4090,
0xfe81,0x1bfd,0x7a69,0x40b3,
0xd118,0xd280,0xad28,0x40c2,
0x3d14,0xa697,0x3d0a,0x40bf,
0x1781,0xb4bd,0x1462,0x40a6,
0x5997,0x6ae7,0x017f,0x4075,
};
#endif
#ifdef MIEEE
static unsigned short QP1[32] = {
0x3faa,0x27fa,0x6b70,0xba40,
0x4013,0xedb5,0xc66d,0x8fd6,
0x4052,0xf4b9,0x9acf,0x1c67,
0x4076,0xec79,0x47aa,0x180d,
0x4086,0x36d9,0xb66f,0x6e50,
0x4082,0xabea,0xb9e8,0x02d0,
0x406a,0x760a,0x4c54,0xbb0b,
0x4039,0x34ff,0x4d15,0x9eb5,
};
static unsigned short QQ1[28] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4052,0x8f30,0x6089,0x5077,
0x4090,0x81cb,0xa20e,0x5f6f,
0x40b3,0x7a69,0x1bfd,0xfe81,
0x40c2,0xad28,0xd280,0xd118,
0x40bf,0x3d0a,0xa697,0x3d14,
0x40a6,0x1462,0xb4bd,0x1781,
0x4075,0x017f,0x6ae7,0x5997,
};
#endif

#ifdef UNK
static double YP1[6] = {
 1.26320474790178026440E9,
-6.47355876379160291031E11,
 1.14509511541823727583E14,
-8.12770255501325109621E15,
 2.02439475713594898196E17,
-7.78877196265950026825E17,
};
static double YQ1[8] = {
/* 1.00000000000000000000E0,*/
 5.94301592346128195359E2,
 2.35564092943068577943E5,
 7.34811944459721705660E7,
 1.87601316108706159478E10,
 3.88231277496238566008E12,
 6.20557727146953693363E14,
 6.87141087355300489866E16,
 3.97270608116560655612E18,
};
#endif
#ifdef DEC
static unsigned short YP1[24] = {
0047626,0112763,0013715,0133045,
0152026,0134552,0142033,0024411,
0053720,0045245,0102210,0077565,
0155347,0000321,0136415,0102031,
0056463,0146550,0055633,0032605,
0157054,0171012,0167361,0054265,
};
static unsigned short YQ1[32] = {
/*0040200,0000000,0000000,0000000,*/
0042424,0111515,0044773,0153014,
0044546,0005405,0171307,0075774,
0046614,0023575,0047105,0063556,
0050613,0143034,0101533,0156026,
0052541,0175367,0166514,0114257,
0054415,0014466,0134350,0171154,
0056164,0017436,0025075,0022101,
0057534,0103614,0103663,0121772,
};
#endif
#ifdef IBMPC
static unsigned short YP1[24] = {
0xb6c5,0x62f9,0xd2be,0x41d2,
0x6521,0x5883,0xd72d,0xc262,
0x0fef,0xb091,0x0954,0x42da,
0xb083,0x37a1,0xe01a,0xc33c,
0x66b1,0x0b73,0x79ad,0x4386,
0x2b17,0x5dde,0x9e41,0xc3a5,
};
static unsigned short YQ1[32] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x7ac2,0xa93f,0x9269,0x4082,
0xef7f,0xbe58,0xc160,0x410c,
0xacee,0xa9c8,0x84ef,0x4191,
0x7b83,0x906b,0x78c3,0x4211,
0x9316,0xfda9,0x3f5e,0x428c,
0x1e4e,0xd71d,0xa326,0x4301,
0xa488,0xc547,0x83e3,0x436e,
0x747f,0x90f6,0x90f1,0x43cb,
};
#endif
#ifdef MIEEE
static unsigned short YP1[24] = {
0x41d2,0xd2be,0x62f9,0xb6c5,
0xc262,0xd72d,0x5883,0x6521,
0x42da,0x0954,0xb091,0x0fef,
0xc33c,0xe01a,0x37a1,0xb083,
0x4386,0x79ad,0x0b73,0x66b1,
0xc3a5,0x9e41,0x5dde,0x2b17,
};
static unsigned short YQ1[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4082,0x9269,0xa93f,0x7ac2,
0x410c,0xc160,0xbe58,0xef7f,
0x4191,0x84ef,0xa9c8,0xacee,
0x4211,0x78c3,0x906b,0x7b83,
0x428c,0x3f5e,0xfda9,0x9316,
0x4301,0xa326,0xd71d,0x1e4e,
0x436e,0x83e3,0xc547,0xa488,
0x43cb,0x90f1,0x90f6,0x747f,
};
#endif


#ifdef UNK
static double Z1 = 1.46819706421238932572E1;
static double Z2 = 4.92184563216946036703E1;
#endif

#ifdef DEC
static unsigned short DZ1[] = {0041152,0164532,0006114,0010540};
static unsigned short DZ2[] = {0041504,0157663,0001625,0020621};
#define Z1 (*(double *)DZ1)
#define Z2 (*(double *)DZ2)
#endif

#ifdef IBMPC
static unsigned short DZ1[] = {0x822c,0x4189,0x5d2b,0x402d};
static unsigned short DZ2[] = {0xa432,0x6072,0x9bf6,0x4048};
#define Z1 (*(double *)DZ1)
#define Z2 (*(double *)DZ2)
#endif

#ifdef MIEEE
static unsigned short DZ1[] = {0x402d,0x5d2b,0x4189,0x822c};
static unsigned short DZ2[] = {0x4048,0x9bf6,0x6072,0xa432};
#define Z1 (*(double *)DZ1)
#define Z2 (*(double *)DZ2)
#endif

#ifndef ANSIPROT
double bessj1(), polevl(), p1evl(), log(), sin(), cos(), sqrt();
#endif

double bessj1(double x)
{
double w, z, p, q, xn;

w = x;
if( x < 0 )
	w = -x;

if( w <= 5.0 )
	{
	z = x * x;	
	w = polevl( z, RP1, 3 ) / p1evl( z, RQ1, 8 );
	w = w * x * (z - Z1) * (z - Z2);
	return( w );
	}

w = 5.0/x;
z = w * w;
p = polevl( z, PP1, 6)/polevl( z, PQ1, 6 );
q = polevl( z, QP1, 7)/p1evl( z, QQ1, 7 );
xn = x - THPIO4;
p = p * cos(xn) - w * q * sin(xn);
return( p * SQ2OPI / sqrt(x) );
}

double bessy1(double x)
{
double w, z, p, q, xn;

if( x <= 5.0 )
	{
	if( x <= 0.0 )
	  throw Exception("arg<=0 in bessy1");
	z = x * x;
	w = x * (polevl( z, YP1, 5 ) / p1evl( z, YQ1, 8 ));
	w += TWOOPI * ( bessj1(x) * log(x)  -  1.0/x );
	return( w );
	}

w = 5.0/x;
z = w * w;
p = polevl( z, PP1, 6)/polevl( z, PQ1, 6 );
q = polevl( z, QP1, 7)/p1evl( z, QQ1, 7 );
xn = x - THPIO4;
p = p * sin(xn) + w * q * cos(xn);
return( p * SQ2OPI / sqrt(x) );
}
