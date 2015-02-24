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
    shared_ptr<FESpace> fes;
    shared_ptr<GridFunction> gfu;
    shared_ptr<BilinearForm> bfa;
    shared_ptr<LinearForm> lff;

  public:
    
    NumProcAllInOne (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    { 
      cout << "Constructor of NumProcAllInOne" << endl;

      Flags flags_fes;
      flags_fes.SetFlag ("order", 4);
      fes = make_shared<H1HighOrderFESpace> (ma, flags_fes);

      Flags flags_gfu;
      gfu = make_shared<T_GridFunction<double>> (fes, "u", flags_gfu);

      Flags flags_bfa;
      bfa = make_shared<T_BilinearFormSymmetric<double>> (fes, "a", flags_bfa);

      shared_ptr<BilinearFormIntegrator> bfi = 
        make_shared<LaplaceIntegrator<2>> (make_shared<ConstantCoefficientFunction> (1));
      bfa -> AddIntegrator (bfi);

      Array<double> penalty(ma->GetNBoundaries());
      penalty = 0.0;
      penalty[0] = 1e10;

      bfi = make_shared<RobinIntegrator<2>> (make_shared<DomainConstantCoefficientFunction> (penalty));
      bfa -> AddIntegrator (bfi);


      bfi = make_shared<MassIntegrator<2>> (make_shared<ConstantCoefficientFunction> (1));
      bfa -> AddIntegrator (bfi);

      Flags flags_lff;
      lff = make_shared<T_LinearForm<double>> (fes, "f", flags_lff);

      auto lfi = make_shared<SourceIntegrator<2>> (make_shared<ConstantCoefficientFunction> (5));
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

      auto inverse = mata.InverseMatrix(fes->GetFreeDofs());
      
      vecu = *inverse * vecf;
    }
  };

  


  static RegisterNumProc<NumProcAllInOne> npinit1("allinone");
}
