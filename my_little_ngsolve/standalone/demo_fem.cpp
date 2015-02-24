/*
  Finite Element demos

  compile with:  'make fem'
*/


// ng-soft header files
#include <fem.hpp>
using namespace ngfem;


using FE_Quad1 = ScalarFE<ET_QUAD,1>;
using FE_Trig1 = ScalarFE<ET_TRIG,1>;



int main ()
{
  // ***************** test integration rule

  // integrate  x*y  over triangle [(0,0), (1,0), (0,1)]
  // with integrationrule of order 10 

  IntegrationRule ir(ET_TRIG, 10);

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
      seg.CalcShape (ip,shapess);
      for(int i = 0; i < seg.GetNDof(); i++)
	outf << shapess(i) << " ";
      outf << endl;
    } 

  
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
  LocalHeap lh(1000000, "demofem - localheap");

  // reference finite element
  FE_Trig1 trig_mapping;
  H1HighOrderFE<ET_TRIG> trig(2);

  // integrators for (\nabla u, \nabla v) and (u,v)
  LaplaceIntegrator<2> laplace(make_shared<ConstantCoefficientFunction>(1));
  MassIntegrator<2> mass(make_shared<ConstantCoefficientFunction>(1));


  // element geometry:
  FE_ElementTransformation<2,2> eltrans;

  // set finite elment for geometry, element-nr and material index
  eltrans.SetElement (&trig_mapping, 0, 0);

  // vertex coordinates  
  eltrans.PointMatrix() = Trans (Matrix<> ({ { 1, 0 }, { 0, 1 }, { 0, 0 } }));
  cout << "PointMatrix = " << endl << eltrans.PointMatrix() << endl;
 
  Matrix<double> elmat_lap (trig.GetNDof());
  laplace.CalcElementMatrix (trig, eltrans, elmat_lap, lh);

  cout << "elmat laplace = " << endl << elmat_lap << endl;



  Matrix<double> elmat_mass(trig.GetNDof());
  mass.CalcElementMatrix (trig, eltrans, elmat_mass, lh);
  
  cout << "elmat mass = " << endl << elmat_mass << endl;

  return 0;
}
