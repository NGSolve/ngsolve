/*********************************************************************/
/* File:   hypre_precond.cpp                                         */
/* Author: Martin Huber, Joachim Schoeberl                           */
/* Date:   July 2012                                                 */
/*********************************************************************/



#include <solve.hpp>

/*
#include "_hypre_utilities.h"

#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
*/

extern ngsolve::PDE * pde;
// #include "hypre_precond.hpp"

namespace ngcomp
{


HyprePreconditioner :: HyprePreconditioner (const PDE & pde, const Flags & flags, const string & name)  
  : Preconditioner (&pde, flags)
{
  cout << IM(1) << "Constructor of HyprePreconditioner" << endl;
  bfa = dynamic_cast<const S_BilinearForm<double>*>
    (pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL)));
}
    

HyprePreconditioner :: HyprePreconditioner (const BaseMatrix & matrix)
  : Preconditioner (pde, Flags())
{
  Setup (matrix);
}


HyprePreconditioner :: ~HyprePreconditioner ()
{
  ;
}



void HyprePreconditioner :: Update()
{
  Setup (bfa->GetMatrix());
}



void HyprePreconditioner :: Setup (const BaseMatrix & matrix)
{
  cout << IM(1) << "Setup Hypre preconditioner" << endl;
  // const ParallelMatrix & pmat = (dynamic_cast<const ParallelMatrix&> (bfa->GetMatrix()));
  const ParallelMatrix & pmat = (dynamic_cast<const ParallelMatrix&> (matrix));
  const SparseMatrix<double> & mat = dynamic_cast< const SparseMatrix<double> &>(pmat.GetMatrix());


  const ParallelDofs * pardofs = pmat.GetParallelDofs ();
  int ndof = pardofs->GetNDof();

  int ntasks = MyMPI_GetNTasks();
  int id = MyMPI_GetId();
    

  // find global dof enumeration 
  global_nums.SetSize(ndof);
  global_nums = -1;
  int num_master_dofs = 0;
  for (int i = 0; i < ndof; i++)
    if (pardofs -> IsMasterDof (i))
      global_nums[i] = num_master_dofs++;
    
  Array<int> first_master_dof(ntasks);
  MPI_Allgather (&num_master_dofs, 1, MPI_INT, 
		 &first_master_dof[0], 1, MPI_INT, 
		 pardofs -> GetCommunicator());
    
  int num_glob_dofs = 0;
  for (int i = 0; i < ntasks; i++)
    {
      int cur = first_master_dof[i];
      first_master_dof[i] = num_glob_dofs;
      num_glob_dofs += cur;
    }
    
  first_master_dof.Append(num_glob_dofs);
  for (int i = 0; i < ndof; i++)
    if (global_nums[i] != -1)
      global_nums[i] += first_master_dof[id];
    
  ScatterDofData (global_nums, *pardofs);
	
	
  // range of my master dofs ...
  ilower = first_master_dof[id];
  iupper = first_master_dof[id+1]-1;
   
  HYPRE_IJMatrixCreate(ngs_comm, ilower, iupper, ilower, iupper, &A);
  HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(A);
   
  /*   
  const BaseMatrix & bmat = (dynamic_cast<const ParallelMatrix&> (bfa->GetMatrix())).GetMatrix();
  const SparseMatrix<double> & mat = dynamic_cast< const SparseMatrix<double> &>(bmat);
  */

  for( int i = 0; i < mat.Height(); i++)
    {
      FlatArray<int> cols = mat.GetRowIndices(i);
      FlatVector<double> values = mat.GetRowValues(i);
      
      Array<int> cols_global(cols.Size());
      for( int j = 0; j < cols.Size(); j++)
	cols_global[j] = global_nums[cols[j]];
      int row = global_nums[i];
      int size = cols.Size();
      
      HYPRE_IJMatrixAddToValues(A, 1, &size, &row, &cols_global[0], &values[0]);
    }

  HYPRE_IJMatrixAssemble(A);
  HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
  // HYPRE_IJMatrixPrint(A, "IJ.out.A");


  HYPRE_BoomerAMGCreate(&precond);
 
  HYPRE_ParVector par_b = NULL;
  HYPRE_ParVector par_x = NULL;
     
  HYPRE_BoomerAMGSetup (precond, parcsr_A, par_b, par_x);
	
  HYPRE_BoomerAMGSetMaxIter (precond,1);
}



void HyprePreconditioner :: Mult (const BaseVector & f, BaseVector & u) const
{
  f.Distribute();
  ParallelBaseVector& pu = dynamic_cast< ParallelBaseVector& > (u);
  pu.SetStatus(DISTRIBUTED);
	
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;

  HYPRE_IJVectorCreate(ngs_comm, ilower, iupper,&b);
  HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(b);

  HYPRE_IJVectorCreate(ngs_comm, ilower, iupper,&x);
  HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(x);
  
  const FlatVector<double> fvf =f.FVDouble();
  FlatVector<double> fu =u.FVDouble();
	
  int local_size = global_nums.Size();

  HYPRE_IJVectorAddToValues(b, local_size, &global_nums[0], &fvf(0));
  HYPRE_IJVectorAssemble(b);

	
  HYPRE_IJVectorGetObject(b, (void **) &par_b);

  HYPRE_IJVectorAssemble(x);
  HYPRE_IJVectorGetObject(x, (void **) &par_x);
   
  // HYPRE_IJVectorPrint(b, "IJ.out.b");

  HYPRE_BoomerAMGSolve(precond, parcsr_A, par_b, par_x);

    
  // ParallelDofs * pardofs = &(bfa->GetFESpace()).GetParallelDofs ();
  const ParallelDofs * pardofs = pu.GetParallelDofs ();
  int ndof = pardofs->GetNDof();
    
  Vector<> hu(iupper-ilower+1);
  Array<int> locind(iupper-ilower+1);

  for (int i = 0, cnt = 0; i < ndof; i++)
    if (pardofs->IsMasterDof(i))
      locind[cnt++] = global_nums[i];
  int lsize = locind.Size();
    
  HYPRE_IJVectorGetValues (x, lsize, &locind[0], &hu(0));
    
  for (int i = 0, cnt = 0; i < ndof; i++)
    if (pardofs->IsMasterDof(i))
      fu(i) = hu(cnt++);
    else
      fu(i) = 0.0;
      
  // (*testout)<<"fu = "<<fu<<endl;
   
  u.Cumulate();
}


static RegisterPreconditioner<HyprePreconditioner> init_hyprepre ("hypre");







}
