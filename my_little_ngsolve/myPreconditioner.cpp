/*********************************************************************/
/* File:   myPreconditioner.cpp                                      */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/


/*

  We define a preconditioner ...

*/


#include <solve.hpp>    // provides FESpace, ...


namespace ngcomp
{
  
  class MyPreconditioner : public Preconditioner
  {
    const BilinearForm * bfa;
    BaseJacobiPrecond * jacobi;

  public:
    MyPreconditioner (const PDE & pde, const Flags & flags);
    ~MyPreconditioner ();
    virtual void Update();
    
    virtual void Mult (const BaseVector & f, BaseVector & u) const
    {
      jacobi -> Mult (f, u);

      /*
      u = 0.0;
      jacobi -> GSSmooth (u, f);
      jacobi -> GSSmoothBack (u, f);
      */
    }

  };


  MyPreconditioner :: MyPreconditioner (const PDE & pde, const Flags & flags)
    : Preconditioner (&pde, flags), jacobi(NULL)
  {
    cout << "Constructor of MyPreconditioner" << endl;

    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
  }
    
  
  MyPreconditioner :: ~MyPreconditioner ()
  {
    delete jacobi;
  }

  
  void MyPreconditioner :: Update()
  {
    delete jacobi;

    const BaseSparseMatrix & mat = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix());
    const BitArray * freedofs = bfa->GetFESpace().GetFreeDofs();

    jacobi = mat.CreateJacobiPrecond (freedofs);
  }

  
  

  static RegisterPreconditioner<MyPreconditioner> initpre ("mypreconditioner");
}

