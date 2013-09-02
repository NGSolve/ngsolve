/*********************************************************************/
/* File:   paralleldofs.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   July 2011                                                 */
/*********************************************************************/

/* 
   parallel dofs
*/


#include <la.hpp>


namespace ngla
{

  ParallelDofs :: ParallelDofs (MPI_Comm acomm, Table<int> * adist_procs, 
				int dim, bool iscomplex)
    : comm(acomm), dist_procs(adist_procs)
    
  {
      int ntasks = MyMPI_GetNTasks(comm);
      int id = MyMPI_GetId(comm);
      
      ndof = dist_procs->Size();
    
      if (ntasks == 1)
	{
	  Array<int> nexdofs(ntasks), ndistprocs(ndof);
	  nexdofs = 0;
	  ndistprocs = 0;
	  
	  exchangedofs = new Table<int> (nexdofs);
	  dist_procs = new Table<int> (ndistprocs);
	  
	  ismasterdof.SetSize(ndof);
	  ismasterdof.Set();
	  return;
	}

      // transpose table
      Array<int> nexdofs(ntasks);
      nexdofs = 0;

      for (int i = 0; i < ndof; i++)
	{
	  FlatArray<int> dist = (*dist_procs)[i];
	  for (int j = 0; j < dist.Size(); j++)
	    nexdofs[dist[j]]++;
	}
      
      exchangedofs = new Table<int> (nexdofs);
      nexdofs = 0;
      
      for (int i = 0; i < ndof; i++)
	{
	  FlatArray<int> dist = (*dist_procs)[i];
	  for (int j = 0; j < dist.Size(); j++)
	    (*exchangedofs)[dist[j]][nexdofs[dist[j]]++] = i;
	}
      
      
      ismasterdof.SetSize (ndof);
      ismasterdof.Set();
      for (int i = 0; i < id; i++)
	{
	  FlatArray<int> ex = (*exchangedofs)[i];
	  for (int j = 0; j < ex.Size(); j++)
	    ismasterdof.Clear (ex[j]);
	}
      
      mpi_t.SetSize (ntasks); 
      
      MPI_Datatype mpi_type;
      
      mpi_type = iscomplex ?
	MyGetMPIType<Complex>() : MyGetMPIType<double>(); 
      
      if (dim != 1)
	{
	  MPI_Datatype htype;
	  MPI_Type_contiguous (dim, mpi_type, &htype);
	  mpi_type = htype;
	}
      
      for (int dest = 0; dest < ntasks; dest++ )
	{
	  if ( !IsExchangeProc(dest) ) continue;
	  
	  FlatArray<int> exdofs = GetExchangeDofs(dest);
	  int len_vec = exdofs.Size();
	  if (len_vec == 0) continue;
	  
	  Array<int> blocklen(len_vec);
	  blocklen = 1;
	  
	  MPI_Type_indexed (len_vec, &blocklen[0], &exdofs[0], 
			    mpi_type, &mpi_t[dest]);
	  MPI_Type_commit (&mpi_t[dest]);
	}
    }


  
  void ParallelDofs :: EnumerateGlobally (const BitArray * freedofs, 
					  Array<int> & global_nums,
					  int & num_glob_dofs) const
  {
    int ntasks = MyMPI_GetNTasks(comm);
    int id = MyMPI_GetId(comm);

    global_nums.SetSize(ndof);
    global_nums = -1;
    int num_master_dofs = 0;
    for (int i = 0; i < ndof; i++)
      if (IsMasterDof (i) && (!freedofs || (freedofs && freedofs->Test(i))))
	global_nums[i] = num_master_dofs++;
    
    Array<int> first_master_dof(ntasks);
    MPI_Allgather (&num_master_dofs, 1, MPI_INT, 
		   &first_master_dof[0], 1, MPI_INT, 
		   GetCommunicator());

    num_glob_dofs = 0;
    for (int i = 0; i < ntasks; i++)
      {
	int cur = first_master_dof[i];
	first_master_dof[i] = num_glob_dofs;
	num_glob_dofs += cur;
      }
    
    for (int i = 0; i < ndof; i++)
      if (global_nums[i] != -1)
	global_nums[i] += first_master_dof[id];
    
    ScatterDofData (global_nums); 
    

  }

}
