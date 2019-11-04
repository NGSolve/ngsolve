#include <fem.hpp>

namespace ngfem
{
  using namespace ngfem;

  
  void SpecialElement :: 
  Apply (FlatVector<double> elx, FlatVector<double> ely,
	 LocalHeap & lh) const
  {
    FlatVector<double> hx1(elx.Size(), lh);
    FlatVector<double> hx2(elx.Size(), lh);

    double dx = 1e-12 + 1e-6 * L2Norm (elx);
    for (int i = 0; i < elx.Size(); i++)
      {
	hx1 = elx;
	hx2 = elx;
	hx1(i) += dx;
	hx2(i) -= dx;
	ely(i) = (Energy (hx1, lh) - Energy(hx2, lh)) / (2 * dx);
      }
    (*testout) << "ely = " << ely << endl;
  }
  
  void SpecialElement ::
  CalcElementMatrix(FlatMatrix<double> elmat,
                    LocalHeap & lh) const
  {
    cerr << "SpecialElement::Assemble not implementd" << endl;
  }
  void SpecialElement ::
  CalcElementVector(FlatVector<double> elvec,
                    LocalHeap & lh) const
  {
    cerr << "SpecialElement::Assemble not implementd" << endl;
  }
  
}
