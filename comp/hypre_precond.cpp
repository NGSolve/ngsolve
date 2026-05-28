#ifdef HYPRE
/*********************************************************************/
/* File:   hypre_precond.cpp                                         */
/* Author: Martin Huber, Joachim Schoeberl                           */
/* Date:   July 2012                                                 */
/*********************************************************************/



#include <solve.hpp>
#include "hypre_precond.hpp"

// In sequential mode, hypre defines MPI_Comm as int but may not expose MPI_COMM_WORLD
// in public headers. Define it as 0 (the conventional stub value) if not defined.
#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
#endif

namespace ngcomp
{

  // Ensure HYPRE_Init is called exactly once.
  // Since this hypre may be built with MPI, we must ensure MPI is initialized first.
  static bool hypre_initialized = false;
  static void EnsureHypreInit() {
    if (!hypre_initialized) {
#ifdef HYPRE_WITH_MPI
      // Hypre built with MPI requires MPI to be initialized before HYPRE_Init()
      int mpi_initialized = 0;
      NG_MPI_Initialized(&mpi_initialized);
      if (!mpi_initialized) {
        NG_MPI_Init(nullptr, nullptr);
      }
#endif
      HYPRE_Init();
      hypre_initialized = true;
    }
  }

#ifdef HYPRE_WITH_MPI
  // Helper to convert NgMPI_Comm to hypre's MPI_Comm type.
  // In NG_MPI_WRAPPER mode, NG_MPI_Comm stores the real MPI_Comm as uintptr_t.
  // In non-wrapper mode, NG_MPI_Comm is directly MPI_Comm.
  static inline MPI_Comm ToHypreComm(NgMPI_Comm ngcomm) {
#ifdef NG_MPI_WRAPPER
    NG_MPI_Comm ng_comm = ngcomm;
    return (MPI_Comm)(ng_comm.value);
#else
    return (NG_MPI_Comm)ngcomm;
#endif
  }
#endif // HYPRE_WITH_MPI



  /*
  HyprePreconditioner :: HyprePreconditioner (const PDE & pde, const Flags & aflags, const string & aname)  
    : Preconditioner (&pde, aflags)
  {
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL));
  }
  */
  

  HyprePreconditioner :: HyprePreconditioner (const BaseMatrix & matrix, const shared_ptr<BitArray> afreedofs)
    : Preconditioner(shared_ptr<BilinearForm>(nullptr), Flags("not_register_for_auto_update"))
  {
    freedofs = afreedofs; 
    Setup (matrix);
  }
  
  HyprePreconditioner :: HyprePreconditioner (shared_ptr<BilinearForm> abfa, const Flags & aflags,
					      const string aname)
    : Preconditioner(abfa, aflags, aname)
  {
    bfa = abfa;
  }


  HyprePreconditioner :: ~HyprePreconditioner ()
  {
    ;
  }
  
  
  
  void HyprePreconditioner :: Update()
  {
    freedofs = bfa->GetFESpace()->GetFreeDofs(bfa->UsesEliminateInternal());
    Setup (bfa->GetMatrix());
  }

  void HyprePreconditioner :: FinalizeLevel (const BaseMatrix * mat)
  {
    freedofs = bfa->GetFESpace()->GetFreeDofs(bfa->UsesEliminateInternal());
    Setup (*mat);
  }
  


  void HyprePreconditioner :: Setup (const BaseMatrix & matrix)
  {
    cout << IM(1) << "Setup Hypre preconditioner" << endl;
    static Timer t("hypre setup");
    RegionTimer reg(t);

    EnsureHypreInit();

    // Determine if we are running in parallel (MPI) or sequential mode
    const SparseMatrix<double> * sp_mat = nullptr;

#ifdef HYPRE_WITH_MPI
    const ParallelMatrix * pmat = dynamic_cast<const ParallelMatrix*>(&matrix);
    if (pmat) {
      sp_mat = &dynamic_cast<const SparseMatrix<double> &>(*pmat->GetMatrix());
      pardofs = pmat->GetParallelDofs();
      parallel = (pardofs && pardofs->GetCommunicator().Size() > 1);
    } else
#endif
    {
      sp_mat = &dynamic_cast<const SparseMatrix<double> &>(matrix);
      pardofs = nullptr;
      parallel = false;
    }

    const SparseMatrix<double> & mat = *sp_mat;
    if (dynamic_cast<const SparseMatrixSymmetric<double> *>(&mat))
      throw Exception ("Please use fully stored sparse matrix for hypre (bf -nonsymmetric)");

    ndof = mat.Height();

    // find global dof enumeration 
    global_nums.SetSize(ndof);
    if(ndof)
      global_nums = -1;
    int num_master_dofs = 0;

    MPI_Comm hypre_comm = MPI_COMM_WORLD;

#ifdef HYPRE_WITH_MPI
    if (parallel) {
      NgMPI_Comm comm = pardofs->GetCommunicator();
      hypre_comm = ToHypreComm(comm);
      int ntasks = comm.Size();
      int id = comm.Rank();

      for (int i = 0; i < ndof; i++)
        if (pardofs->IsMasterDof(i) && (!freedofs || freedofs->Test(i)))
          global_nums[i] = num_master_dofs++;

      Array<int> first_master_dof(ntasks);
      NG_MPI_Allgather (&num_master_dofs, 1, NG_MPI_INT,
                        &first_master_dof[0], 1, NG_MPI_INT,
                        pardofs->GetCommunicator());
    
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

      ScatterDofData (global_nums, pardofs);
      cout << IM(3) << "num glob dofs = " << num_glob_dofs << endl;
	
      // range of my master dofs ...
      ilower = first_master_dof[id];
      iupper = first_master_dof[id+1]-1;
    } else
#endif // HYPRE_WITH_MPI
    {
      // Sequential: global numbering = local numbering for free dofs
      for (int i = 0; i < ndof; i++)
        if (!freedofs || freedofs->Test(i))
          global_nums[i] = num_master_dofs++;

      cout << IM(3) << "num dofs = " << num_master_dofs << endl;
      ilower = 0;
      iupper = num_master_dofs - 1;
    }

    HYPRE_IJMatrixCreate(hypre_comm, ilower, iupper, ilower, iupper, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);
   

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

    HYPRE_BoomerAMGCreate(&precond);
 
    HYPRE_ParVector par_b = NULL;
    HYPRE_ParVector par_x = NULL;

    int printlevel = flags.GetNumFlag("printlevel", 0);

    HYPRE_BoomerAMGSetPrintLevel(precond, printlevel);  /* 0=none, 1=setup, 2=solve, 3=both */
    HYPRE_BoomerAMGSetCoarsenType(precond, 10); /* Falgout coarsening */
    HYPRE_BoomerAMGSetRelaxType(precond, 6);  // 3 GS, 6 .. sym GS 
    HYPRE_BoomerAMGSetStrongThreshold(precond, 0.5);
    HYPRE_BoomerAMGSetInterpType(precond,6);
    HYPRE_BoomerAMGSetPMaxElmts(precond,4);
    HYPRE_BoomerAMGSetAggNumLevels(precond,1);
    HYPRE_BoomerAMGSetNumSweeps(precond, 1);   /* Sweeeps on each level */
    HYPRE_BoomerAMGSetMaxLevels(precond, 20);  /* maximum number of levels */
    HYPRE_BoomerAMGSetTol(precond, 0.0);      /* conv. tolerance */
    HYPRE_BoomerAMGSetMaxIter(precond,1);
     
    cout << IM(2) << "Call BoomerAMGSetup" << endl;
    HYPRE_BoomerAMGSetup (precond, parcsr_A, par_b, par_x);
	
  }



  void HyprePreconditioner :: Mult (const BaseVector & f, BaseVector & u) const
  {
    static Timer t("hypre mult");
    RegionTimer reg(t);

    MPI_Comm hypre_comm = MPI_COMM_WORLD;

#ifdef HYPRE_WITH_MPI
    if (parallel) {
      hypre_comm = ToHypreComm(pardofs->GetCommunicator());
      f.Distribute();
      u.SetParallelStatus(DISTRIBUTED);
    }
#endif

    HYPRE_IJVector b;
    HYPRE_ParVector par_b;
    HYPRE_IJVector x;
    HYPRE_ParVector par_x;

    HYPRE_IJVectorCreate(hypre_comm, ilower, iupper, &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVectorCreate(hypre_comm, ilower, iupper, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);
  
    const FlatVector<double> fvf = f.FVDouble();
    FlatVector<double> fu = u.FVDouble();
	
    Array<int> nzglobal;
    Array<double> free_f;
    for (auto i : Range(global_nums.Size()))
      if (global_nums[i] != -1)
	{
	  nzglobal.Append (global_nums[i]);
	  free_f.Append (fvf(i));
	}
        
    Array<int> setzi;
    Array<double> zeros;
    if (iupper >= ilower)
      for (auto k : Range(ilower, iupper+1))
	{
	  setzi.Append(k);
	  zeros.Append(0.0);
	}

    int local_size = nzglobal.Size();

    // set vector to 0
    if (setzi.Size())
      HYPRE_IJVectorSetValues(b, setzi.Size(), &setzi[0], &zeros[0]);
    // add values
    if (local_size)
      HYPRE_IJVectorAddToValues(b, local_size, &nzglobal[0], &free_f[0]);
    HYPRE_IJVectorAssemble(b);

    HYPRE_IJVectorGetObject(b, (void **) &par_b);

    if (setzi.Size())
      HYPRE_IJVectorSetValues(x, setzi.Size(), &setzi[0], &zeros[0]);
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);

    HYPRE_BoomerAMGSolve(precond, parcsr_A, par_b, par_x);

    Vector<> hu(iupper-ilower+1);
    Array<int> locind(iupper-ilower+1);

#ifdef HYPRE_WITH_MPI
    if (parallel) {
      for (int i = 0, cnt = 0; i < ndof; i++)
        if (pardofs->IsMasterDof(i) && global_nums[i] != -1)
          locind[cnt++] = global_nums[i];
    } else
#endif
    {
      for (int i = 0, cnt = 0; i < ndof; i++)
        if (global_nums[i] != -1)
          locind[cnt++] = global_nums[i];
    }
    int lsize = locind.Size();
    
    HYPRE_IJVectorGetValues (x, lsize, &locind[0], &hu(0));

#ifdef HYPRE_WITH_MPI
    if (parallel) {
      for (int i = 0, cnt = 0; i < ndof; i++)
        if (pardofs->IsMasterDof(i) && global_nums[i] != -1)
          fu(i) = hu(cnt++);
        else
          fu(i) = 0.0;
      u.Cumulate();
    } else
#endif
    {
      for (int i = 0, cnt = 0; i < ndof; i++)
        if (global_nums[i] != -1)
          fu(i) = hu(cnt++);
        else
          fu(i) = 0.0;
    }

    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);
  }


  static RegisterPreconditioner<HyprePreconditioner> init_hyprepre ("hypre");
}
#endif
