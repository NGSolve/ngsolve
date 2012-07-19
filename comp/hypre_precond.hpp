#ifdef HYPRE

#ifndef FILE_HYPRE_PRECOND
#define FILE_HYPRE_PRECOND


/*********************************************************************/
/* File:   hypre_precond.hpp                                         */
/* Author: Martin Huber, Joachim Schoeberl                           */
/* Date:   July 2012                                                 */
/*********************************************************************/


#include "_hypre_utilities.h"

#include "HYPRE_krylov.h"
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
  // int ne;
  // int nse;
	
  int ilower, iupper;

public:
    
  HyprePreconditioner (const PDE & pde, const Flags & flags, const string & name);
  HyprePreconditioner (const BaseMatrix & matrix); // 

  ~HyprePreconditioner ();
	
  virtual void Update();
  virtual void Mult (const BaseVector & f, BaseVector & u) const;

private:
  void Setup (const BaseMatrix & matrix); 
};

}



#endif
#endif






