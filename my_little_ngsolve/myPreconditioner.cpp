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
    shared_ptr<BilinearForm> bfa;
    shared_ptr<BaseJacobiPrecond> jacobi;

  public:
    MyPreconditioner (const PDE & pde, const Flags & flags, const string & aname);
    MyPreconditioner (shared_ptr<BilinearForm> & abfa, const Flags & flags, const string & aname);
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

    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa -> GetMatrix();
    }

  };


  MyPreconditioner :: MyPreconditioner (const PDE & pde, const Flags & flags, const string & aname)
    : Preconditioner (&pde, flags, aname), jacobi(NULL)
  {
    cout << "Constructor of MyPreconditioner" << endl;

    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
  }

  MyPreconditioner :: MyPreconditioner (shared_ptr<BilinearForm> & abfa, const Flags & flags, const string & aname)
    : Preconditioner (abfa, flags, name), bfa(abfa)
  { ; }

  
  MyPreconditioner :: ~MyPreconditioner ()
  {
    ;
  }

  
  void MyPreconditioner :: Update()
  {
    const BaseSparseMatrix & mat = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix());
    const BitArray * freedofs = bfa->GetFESpace()->GetFreeDofs(bfa->UsesEliminateInternal());

    jacobi = mat.CreateJacobiPrecond (freedofs);

    if (test) Test();
  }

  
  

  static RegisterPreconditioner<MyPreconditioner> initpre ("mypreconditioner");
}

