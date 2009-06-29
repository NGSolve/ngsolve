/*
  Finite Element demos

  compile with:  'make fem'
*/

#include <iostream>
#include <sstream>
#include <iomanip>

#include <time.h>
#include <math.h>
#include <complex>



// ng-soft header files
#include <fem.hpp>
#include <../fem/h1lofe.hpp>


using namespace std;
using namespace ngstd;
using namespace ngbla;
using namespace ngfem;


namespace netgen {
  ostream * testout = &cout;
}


int main ()
{
  // ***************** test integration rule

  // integrate  x*y  over triangle [(0,0), (1,0), (0,1)]
  // with integrationrule of order 10 

  const IntegrationRule & ir = 
    GetIntegrationRules().SelectIntegrationRule (ET_TRIG, 10);

  cout << "number of ipts = " << ir.GetNIP() << endl;

  double sum = 0;
  for (int i = 0; i < ir.GetNIP(); i++)
    {
      const IntegrationPoint & ip = ir[i];
      sum += ip.Weight() * ip(0) * ip(1);
    }
  
  cout << "Integral = " << sum << " =?= " << 1.0/24.0 << endl << endl;

  // **************** Compute shape functions


  H1HighOrderFE<ET_SEGM> seg(4);
  Vector<> shapess(seg.GetNDof());
  cout << "ndof = " << seg.GetNDof() << endl;
  ofstream outf("shape");

  for(double x = 0; x <= 1.0001; x += 0.1)
    {
      IntegrationPoint ip(x);
      seg.CalcShape(ip,shapess);
      for(int i=0; i<seg.GetNDof(); i++)
	outf << shapess(i) << " ";
      outf << endl;
    }
  outf.close();


  
  FE_Quad1 quad;
  Vector<> shapes(quad.GetNDof());
  for (double x = 0; x <= 1.0001; x += 0.25)
    for (double y = 0; y <= 1.0001; y += 0.25)
      {
	IntegrationPoint ip (x, y);
	quad.CalcShape (ip, shapes);
	cout << x << ", " << y << ": " << shapes << endl;
      }
  cout << endl << endl;


  // ***************** element matrix integration

  // own memory management
  LocalHeap lh(10000);

  // reference finite element
  FE_Trig1 trig;

  // integrators for (\nabla u, \nabla v) and (u,v)
  ConstantCoefficientFunction coef(1);
  LaplaceIntegrator<2> laplace(&coef);
  MassIntegrator<2> mass(&coef);


  // element geometry:
  ElementTransformation eltrans;

  eltrans.SetElement (&trig, 0, 0);
  

  // vertex coordinates  
  double pts[3][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 } };
  eltrans.AllocPointMatrix (2,3);
  eltrans.PointMatrix() = Trans (FlatMatrix<double> (3, 2, &pts[0][0]));

  cout << "PointMatrix = " << endl << eltrans.PointMatrix() << endl;

 
  FlatMatrix<double> elmat_lap (trig.GetNDof(), lh);
  laplace.AssembleElementMatrix (trig, eltrans, elmat_lap, lh);

  cout << "elmat laplace = " << endl << elmat_lap << endl;



  FlatMatrix<double> elmat_mass(trig.GetNDof(), lh);
  mass.AssembleElementMatrix (trig, eltrans, elmat_mass, lh);
  
  cout << "elmat mass = " << endl << elmat_mass << endl;

  return 0;
}
