#ifdef HYPRE

#ifndef FILE_HYPRE_PRECOND
#define FILE_HYPRE_PRECOND


/*********************************************************************/
/* File:   hypre_precond.hpp                                         */
/* Author: Martin Huber, Joachim Schoeberl                           */
/* Date:   July 2012                                                 */
/*********************************************************************/


#include "HYPRE_utilities.h"
// #include "HYPRE_krylov.h"

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"


namespace ngcomp
{

  
class HyprePreconditioner : public Preconditioner
{
  shared_ptr<BilinearForm> bfa;

  HYPRE_Solver precond;
  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;

  Array<int> global_nums;
  int ilower, iupper;
  shared_ptr<BitArray> freedofs;
  shared_ptr<ParallelDofs> pardofs;

public:
    
  // HyprePreconditioner (const PDE & pde, const Flags & flags, const string & name);
  HyprePreconditioner (const BaseMatrix & matrix, const shared_ptr<BitArray> afreedofs); 
  HyprePreconditioner (shared_ptr<BilinearForm> bfa, const Flags & aflags,
		       const string aname = "precond");

  ~HyprePreconditioner ();
	
  virtual void FinalizeLevel (const ngla::BaseMatrix * mat = NULL) override;
  virtual void Update() override;
  virtual void Mult (const BaseVector & f, BaseVector & u) const override;
  virtual int VHeight() const override { return pardofs->GetNDofLocal();}
  virtual int VWidth() const override { return pardofs->GetNDofLocal();}
  virtual const BaseMatrix & GetAMatrix() const override { return bfa->GetMatrix(); }
  virtual const BaseMatrix & GetMatrix() const override { return *this; }

  virtual const char * ClassName() const override
  { return "HYPRE AMG Preconditioner"; }

private:
  void Setup (const BaseMatrix & matrix);
};

}



#endif
#endif






