#ifdef HYPRE
/*********************************************************************/
/* File:   hypre_precond.cpp                                         */
/* Author: Martin Huber, Joachim Schoeberl                           */
/* Date:   July 2012                                                 */
/*********************************************************************/



#include <solve.hpp>

extern ngsolve::PDE * pde;
namespace ngcomp
{


  HyprePreconditioner :: HyprePreconditioner (const PDE & pde, const Flags & flags, const string & name)  
    : Preconditioner (&pde, flags)
  {
    bfa = dynamic_cast<const S_BilinearForm<double>*>
      (pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL)));
  }
  

  HyprePreconditioner :: HyprePreconditioner (const BaseMatrix & matrix, const BitArray * afreedofs)
    : Preconditioner (pde, Flags()), freedofs(afreedofs)
  {
    Setup (matrix);
  }
  

  HyprePreconditioner :: ~HyprePreconditioner ()
  {
    ;
  }
  
  
  
  void HyprePreconditioner :: Update()
  {
    freedofs = bfa->GetFESpace().GetFreeDofs (bfa->UsesEliminateInternal());
    Setup (bfa->GetMatrix());
  }
  


  void HyprePreconditioner :: Setup (const BaseMatrix & matrix)
  {
    cout << IM(1) << "Setup Hypre preconditioner" << endl;
    static Timer t("hypre setup");
    RegionTimer reg(t);

    const ParallelMatrix & pmat = (dynamic_cast<const ParallelMatrix&> (matrix));
    const SparseMatrix<double> & mat = dynamic_cast< const SparseMatrix<double> &>(pmat.GetMatrix());


    pardofs = pmat.GetParallelDofs ();
    int ndof = pardofs->GetNDof();

    int ntasks = MyMPI_GetNTasks();
    int id = MyMPI_GetId();
    

    // find global dof enumeration 
    global_nums.SetSize(ndof);
    global_nums = -1;
    int num_master_dofs = 0;
    for (int i = 0; i < ndof; i++)
      if (pardofs -> IsMasterDof (i) && (!freedofs || freedofs -> Test(i)))
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
    cout << IM(3) << "num glob dofs = " << num_glob_dofs << endl;
	
    VT_OFF();
	
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
	int row = global_nums[i];
	if (row == -1) continue;

	FlatArray<int> cols = mat.GetRowIndices(i);
	FlatVector<double> values = mat.GetRowValues(i);
      
	Array<int> cols_global;
	Array<double> values_global;

	for( int j = 0; j < cols.Size(); j++)
	  if (global_nums[cols[j]] != -1)
	    {
	      cols_global.Append (global_nums[cols[j]]);
	      values_global.Append (values[j]);
	    }

	int size = cols_global.Size();
	HYPRE_IJMatrixAddToValues(A, 1, &size, &row, &cols_global[0], &values_global[0]);
      }

    HYPRE_IJMatrixAssemble(A);
    HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
    // HYPRE_IJMatrixPrint(A, "IJ.out.A");

    HYPRE_BoomerAMGCreate(&precond);
 
    HYPRE_ParVector par_b = NULL;
    HYPRE_ParVector par_x = NULL;
     
    cout << IM(2) << "Call BoomerAMGSetup" << endl;
    HYPRE_BoomerAMGSetup (precond, parcsr_A, par_b, par_x);
	
    HYPRE_BoomerAMGSetMaxIter (precond,3);
    VT_ON();
  }



  void HyprePreconditioner :: Mult (const BaseVector & f, BaseVector & u) const
  {
    static Timer t("hypre mult");
    RegionTimer reg(t);

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
	
    Array<int> nzglobal;
    Array<double> free_f;
    for (int i = 0; i < global_nums.Size(); i++)
      if (global_nums[i] != -1)
	{
	  nzglobal.Append (global_nums[i]);
	  free_f.Append (fvf(i));
	}

    int local_size = nzglobal.Size();


    HYPRE_IJVectorAddToValues(b, local_size, &nzglobal[0], &free_f[0]);
    HYPRE_IJVectorAssemble(b);

	
    HYPRE_IJVectorGetObject(b, (void **) &par_b);

    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);
   
    // HYPRE_IJVectorPrint(b, "IJ.out.b");

    VT_OFF();
    HYPRE_BoomerAMGSolve(precond, parcsr_A, par_b, par_x);
    VT_ON();

    // HYPRE_IJVectorPrint(x, "IJ.out.x");
    
    // ParallelDofs * pardofs = &(bfa->GetFESpace()).GetParallelDofs ();
    // const ParallelDofs * pardofs = pu.GetParallelDofs ();

    int ndof = pardofs->GetNDof();
    
    Vector<> hu(iupper-ilower+1);
    Array<int> locind(iupper-ilower+1);

    for (int i = 0, cnt = 0; i < ndof; i++)
      if (pardofs->IsMasterDof(i) && global_nums[i] != -1)
	locind[cnt++] = global_nums[i];
    int lsize = locind.Size();
    
    HYPRE_IJVectorGetValues (x, lsize, &locind[0], &hu(0));
    
    for (int i = 0, cnt = 0; i < ndof; i++)
      if (pardofs->IsMasterDof(i) && global_nums[i] != -1)
	fu(i) = hu(cnt++);
      else
	fu(i) = 0.0;

    u.Cumulate();
  }


  static RegisterPreconditioner<HyprePreconditioner> init_hyprepre ("hypre");
}
#endif
