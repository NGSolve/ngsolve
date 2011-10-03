/*
  Basic Linear Algebra demos
 */


#include <iostream>
#include <iomanip>
#include <complex>

#include <ctime>
#include <cmath>



// ng-soft header files
#include <bla.hpp>


using namespace std;
using namespace ngbla;



int main ()
{
  Vector<double> u(3), v(5);
  Matrix<double> m(5,3);

  u = 1.0;

  for (int i = 0; i < Height(m); i++)
    for (int j = 0; j < Width(m); j++)
      m(i,j) = i+j;

  
  // do some math:
  v = m * u;

  cout << "u = " << u << endl;
  cout << "v = " << v << endl;
  cout << "m = " << m << endl;

  // use result directly
  cout << "Trans(m) * v = " << Trans(m) * v << endl;

  // fix-size objects:
  Vec<3,double> u3;
  Mat<3,3,double> m3;

  u3 = 2.7;
  m3 = 1.0;
  cout << "m3 * u3 = " << m3 * u3 << endl;


  // own memory management:
  double data[1000];
  
  FlatVector<double> fu(100, data);
  FlatVector<double> fv(100, data+100);

  // overlay
  FlatVector<double> fw(200, data);

  fu = 1.0;
  fv = fu;

  cout << "(u,v) = " << InnerProduct (fu, fv) << endl;
  cout << "(w,w) = " << InnerProduct (fw, fw) << endl;

  // more complicated vectors
  Vector<Vec<2,Complex> > sysu(4);
  Vector<Vec<3,Complex> > sysv(3);
  Matrix<Mat<3,2,Complex> > sysm(3,4);

  for (int i = 0; i < sysu.Size(); i++)
    {
      sysu(i)(0) = Complex (1,0);
      sysu(i)(1) = Complex (0,1);
    }

  cout << "sysu = " << sysu << endl;

  for (int i = 0; i < sysm.Height(); i++)
    for (int j = 0; j < sysm.Width(); j++)
      {
	sysm(i,j) = 0.0;
	sysm(i,j)(1,0) = 2.5;
      }

  cout << "sysm = " << sysm << endl;

  sysv = sysm * sysu;

  cout << "sysv = " << sysv << endl;

  Vector<double> a(1000);
  Vector<double> b(1000);
  a = 1.0;
  b = 1.0;

  clock_t starttime, endtime;


  starttime = clock();

  double sum = 0;
  for (int i = 0; i < 1000000; i++)
    sum += InnerProduct (a, b);

  endtime = clock();
  cout << "time = " << double(endtime - starttime) / CLOCKS_PER_SEC << endl;
  cout << "sum = " << sum << endl;


  //  Vector<Complex> ac(1000);
  //  Vector<Complex> bc(1000);
  Vector<complex<double> > ac(1000);
  Vector<complex<double> > bc(1000);
  //  Vector<Vec<2> > ac(1000);
  //  Vector<Vec<2> > bc(1000);
  ac = 1.0;
  bc = 1.0;

  starttime = clock();

  Complex sumc = 0;
  for (int i = 0; i < 100000; i++)
    sumc += InnerProduct (ac, bc);

  endtime = clock();
  cout << "time = " << double(endtime - starttime) / CLOCKS_PER_SEC << endl;
  cout << "sum = " << sumc << endl;



  LocalHeap lh(10000, "demobla - localheap");
  FlatMatrix<Mat<2,2,double> > m2(2, lh), invm2(2, lh);
  m2 = 0;
  for (int i = 0; i < m2.Height(); i++)
    {
      m2(i,i)(0,0) = m2(i,i)(1,1) = i+1;
      m2(i,i)(1,0) = m2(i,i)(0,1) = 0.1;
      for (int j = 0; j < i; j++)
	m2(i,j)(1,0) = j+5;
    }
  CalcInverse (m2, invm2);

  cout << "m = " << m2 << endl;
  cout << "invm = " << invm2 << endl;
  cout << "check = " << (m2 * invm2) << endl;


  // test cholesky factors:
  Matrix<double> sm(5), invsm(5);
  Vector<double> vx(5), vy(5);
  sm = 0;
  for (int i = 0; i < Height(sm); i++)
    for (int j = 0; j < Width(sm); j++)
      sm(i,j) = i*j + 1.0 / (1+i+j);
  
  CalcInverse (sm, invsm);
  cout << "sm = " << sm << endl;
  cout << "invsm = " << invsm << endl;

  CholeskyFactors<double> cholm(sm);
  cout << cholm; // .Print (cout);
  vx = 0;
  vx(4) = 1;
  cholm.Mult (vx, vy);
  cout << "inv * e1 = " << vy << endl;

  
  Matrix<double> evecs(5);
  Vector<double> lami(5);
  CalcEigenSystem (sm, lami, evecs);
  cout << "lami = " << lami << endl 
       << "evecs = " << evecs << endl;


  
  Array<int> rows(2);
  rows[0] = 1;
  rows[1] = 3;
  cout << "rows " << rows << " of sm are " << endl << sm.Rows(rows).Cols(rows) << endl;

  Matrix<> bmat(5,5), cmat(5,5);
  bmat = 0.0;
  cmat = 1.1;
  bmat.Rows(rows) = cmat.Rows(2,3);
  bmat.Cols(rows) = cmat.Cols(2,3);
  cout << "bmat = " << endl << bmat << endl;

  {
    Matrix<> a(1000), b(1000), c(1000);
    a = 1;
    b = 2;
    LapackMultAB (a, b, c);
  }

  return 0;
}



/*
  icc7.1 ... 2.37 sec  (3.67, 2.76, ..),   Complex .. 200
  gcc3.3 ...  Complex 15
 */
