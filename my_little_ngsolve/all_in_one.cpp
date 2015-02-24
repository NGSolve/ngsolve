/*********************************************************************/
/* File:   all_in_one.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   4. May. 2009                                              */
/*********************************************************************/


/*

We write one numproc which creates and maintains all objects to solve
the Poisson equation.

*/

#include <solve.hpp>

using namespace ngsolve;

namespace all_in_one
{
  class NumProcAllInOne : public NumProc
  {
  protected:
    FESpace * fes;
    GridFunction *gfu;
    BilinearForm *bfa;
    LinearForm * lff;

  public:
    
    NumProcAllInOne (PDE & apde, const Flags & flags)
      : NumProc (apde)
    { 
      cout << "Constructor of NumProcAllInOne" << endl;

      Flags flags_fes;
      flags_fes.SetFlag ("order", 4);
      fes = new H1HighOrderFESpace (ma, flags_fes);

      Flags flags_gfu;
      gfu = new T_GridFunction<double> (*fes, "u", flags_gfu);

      Flags flags_bfa;
      bfa = new T_BilinearFormSymmetric<double> (*fes, "a", flags_bfa);

      BilinearFormIntegrator * bfi;
      bfi = new LaplaceIntegrator<2> (new ConstantCoefficientFunction (1));
      bfa -> AddIntegrator (bfi);

      Array<double> penalty(ma.GetNBoundaries());
      penalty = 0.0;
      penalty[0] = 1e10;

      bfi = new RobinIntegrator<2> (new DomainConstantCoefficientFunction (penalty));
      bfa -> AddIntegrator (bfi);


      bfi = new MassIntegrator<2> (new ConstantCoefficientFunction (1));
      bfa -> AddIntegrator (bfi);

      Flags flags_lff;
      lff = new T_LinearForm<double> (*fes, "f", flags_lff);

      LinearFormIntegrator * lfi;
      lfi = new SourceIntegrator<2> (new ConstantCoefficientFunction (5));
      lff -> AddIntegrator (lfi);
    }

    virtual string GetClassName () const
    {
      return "AllInOne";
    }

    virtual void Do (LocalHeap & lh)
    {
      fes -> Update(lh);
      fes -> FinalizeUpdate(lh);

      gfu -> Update();
      bfa -> Assemble(lh);
      lff -> Assemble(lh);

      const BaseMatrix & mata = bfa -> GetMatrix();
      const BaseVector & vecf = lff -> GetVector();
      BaseVector & vecu = gfu -> GetVector();

      BaseMatrix * inverse = mata.InverseMatrix(fes->GetFreeDofs());
      
      vecu = (*inverse) * vecf;

      delete inverse;
    }
  };

  


  static RegisterNumProc<NumProcAllInOne> npinit1("allinone");
}
