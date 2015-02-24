// ngscxx -shared dd_precond.cpp -lngcomp -o mydd.so


/*

parallel preconditioners:

Version 1:
diagonal

Version 2:
block-diagonal (inner, coupling)

Version 3:
sub-structuring 

*/


#include <comp.hpp>
using namespace ngcomp;



namespace ngcomp
{
  
  class MyDDPreconditioner : public Preconditioner
  {
    shared_ptr<BilinearForm> bfa;
    const SparseMatrix<double> * local_matrix;
    const ParallelDofs * pardofs;
    Array<double> inv_diag;
    shared_ptr<BaseMatrix> inv_inner;
    BitArray * localdofs;

  public:
    MyDDPreconditioner (const PDE & pde, const Flags & flags, const string & aname)
      : Preconditioner (&pde, flags, aname)
    {
      bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    }

    ~MyDDPreconditioner ()
    {
      ;
    }

    virtual void Update()
    {
      const ParallelMatrix * hm = dynamic_cast<const ParallelMatrix*> (&bfa->GetMatrix());
      if (!hm) throw Exception ("MyDD Preconditioner: not a parallel matrix");

      pardofs = hm -> GetParallelDofs();
      local_matrix = dynamic_cast<const SparseMatrix<double>*> (&hm->GetMatrix());
      const BitArray * freedofs = bfa->GetFESpace()->GetFreeDofs(bfa->UsesEliminateInternal());
      
      Array<double> diag (local_matrix->Height());
      for (int i = 0; i < diag.Size(); i++)
	diag[i] = (*local_matrix)(i,i);
      
      pardofs -> AllReduceDofData (diag, MPI_SUM);  // sum up local diagonals across processors
      
      inv_diag.SetSize (diag.Size());
      for (int i = 0; i < diag.Size(); i++)
	if (freedofs->Test(i))
	  inv_diag[i] = 1.0 / diag[i];
	else
	  inv_diag[i] = 0.0;

      localdofs = new BitArray (*freedofs);
      for (int i = 0; i < localdofs->Size(); i++)
	if (pardofs->GetDistantProcs(i).Size() != 0)
	  localdofs->Clear(i);

      local_matrix->SetInverseType (SPARSECHOLESKY);
      inv_inner = local_matrix->InverseMatrix (localdofs);

      if (test) Test();
    }
    
    virtual void Mult (const BaseVector & f, BaseVector & u) const
    {
     
      // diagonal preconditioner
      f.Distribute();
      
      FlatVector<double> ff = f.FV<double>();      
      FlatVector<double> fu = u.FV<double>();
      
      for (int i = 0; i < inv_diag.Size(); i++)
	fu(i) = inv_diag[i] * ff(i);

      u.SetParallelStatus (DISTRIBUTED);
      u.Cumulate();
     
      
/*
	// block-diagonal with inner block
      f.Distribute();
      
      FlatVector<double> ff = f.FV<double>();      
      FlatVector<double> fu = u.FV<double>();
      
      u = (*inv_inner) * f;
      for (int i = 0; i < inv_diag.Size(); i++)
	fu(i) += inv_diag[i] * ff(i);

      u.SetParallelStatus (DISTRIBUTED);
      u.Cumulate();
      */

      /*
      // sub-structuring
      f.Distribute();
      auto hv1 = f.CreateVector();
      auto hv2 = f.CreateVector();
      FlatVector<double> fhv1 = hv1.FV<double>();      
      FlatVector<double> fhv2 = hv2.FV<double>();
      
      
      hv2 = f;
      hv1 = (*inv_inner) * hv2;
      hv2 -= (*local_matrix) * hv1;

      for (int i = 0; i < inv_diag.Size(); i++)
	fhv1(i) = inv_diag[i] * fhv2(i);
      hv1.SetParallelStatus(DISTRIBUTED);
      hv1.Cumulate();

      hv2 = (*local_matrix) * hv1;
      hv1 -= (*inv_inner) * hv2;
      
      hv1 += (*inv_inner) * f;
      u = hv1;
      */
    }

    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa -> GetMatrix();
    }

  };



  

  static RegisterPreconditioner<MyDDPreconditioner> initpre ("mydd");
}

