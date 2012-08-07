#ifdef HYPRE

#ifndef FILE_HYPRE_PRECOND
#define FILE_HYPRE_PRECOND


/*********************************************************************/
/* File:   hypre_precond.hpp                                         */
/* Author: Martin Huber, Joachim Schoeberl                           */
/* Date:   July 2012                                                 */
/*********************************************************************/


// #include "_hypre_utilities.h"
// #include "HYPRE_krylov.h"

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"


namespace ngcomp
{

  
class HyprePreconditioner : public Preconditioner
{
  const S_BilinearForm<double> * bfa;

  HYPRE_Solver precond;
  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;

  Array<int> global_nums;
  int ilower, iupper;
  const BitArray * freedofs;
  const ParallelDofs * pardofs;

public:
    
  HyprePreconditioner (const PDE & pde, const Flags & flags, const string & name);
  HyprePreconditioner (const BaseMatrix & matrix, const BitArray * afreedofs); 

  ~HyprePreconditioner ();
	
  virtual void Update();
  virtual void Mult (const BaseVector & f, BaseVector & u) const;

  virtual const char * ClassName() const
  { return "HYPRE AMG Preconditioner"; }

private:
  void Setup (const BaseMatrix & matrix);
};

}



#endif
#endif






