/*********************************************************************/
/* File:   paralleldofs.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   July 2011                                                 */
/*********************************************************************/

/* 
   parallel dofs
*/


#include <la.hpp>

#ifdef PARALLEL

namespace ngla
{

  // extern void MyFunction();  ????
  
  ParallelDofs :: ParallelDofs (MPI_Comm acomm, Table<int> && adist_procs, 
				int dim, bool iscomplex)
    : comm(acomm), dist_procs(adist_procs), es(dim), complex(iscomplex)
    
  {
    int ntasks = comm.Size();
    int id = comm.Rank();
      
    ndof = dist_procs.Size();
    
    if (ntasks == 1)
      {
	Array<int> nexdofs(ntasks), ndistprocs(ndof);
	nexdofs = 0;
	ndistprocs = 0;
	
	exchangedofs = Table<int> (nexdofs);
	dist_procs = Table<int> (ndistprocs);
	
	ismasterdof.SetSize(ndof);
	ismasterdof.Set();

        global_ndof = ndof;
	return;
      }
    
    // transpose table
    Array<int> nexdofs(ntasks);
    nexdofs = 0;
    
    for (int i = 0; i < ndof; i++)
      for (int d : dist_procs[i])
	nexdofs[d]++;
    
    exchangedofs = Table<int>(nexdofs);
    nexdofs = 0;

    for (int i = 0; i < ndof; i++)
      for (int d : dist_procs[i])
	exchangedofs[d][nexdofs[d]++] = i;
      
    ismasterdof.SetSize (ndof);
    ismasterdof.Set();

    for (int i = 0; i < id; i++)
      for (int ex : exchangedofs[i])
	ismasterdof.Clear (ex);


    mpi_t.SetSize (ntasks); 
    
    MPI_Datatype mpi_type;
    
    mpi_type = iscomplex ?
      MyGetMPIType<Complex>() : MyGetMPIType<double>(); 
    
    if (es != 1)
      {
	MPI_Datatype htype;
	MPI_Type_contiguous (es, mpi_type, &htype);
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
    
    for (int i = 0; i < ntasks; i++)
      if (IsExchangeProc (i))
	all_dist_procs.Append (i);



    size_t nlocal = 0;
    for (int i = 0; i < ndof; i++)
      if (ismasterdof.Test(i)) nlocal++;
    global_ndof = comm.AllReduce (nlocal, MPI_SUM);
  }

  ParallelDofs :: ~ParallelDofs()  {
    for (auto dest : Range(mpi_t.Size()))
      if ( IsExchangeProc(dest) )
	MPI_Type_free(&mpi_t[dest]);
  }

  shared_ptr<ParallelDofs> ParallelDofs :: SubSet (shared_ptr<BitArray> take_dofs) const
  {
    auto ndloc = this->GetNDofLocal();
    Array<size_t> s(ndloc);
    for(auto k:Range(ndloc))
      if(take_dofs->Test(k)) s[k] = this->GetDistantProcs(k).Size();
      else s[k] = 0;
    Table<int> tab(s);
    for(auto k:Range(ndloc))
      if(take_dofs->Test(k)) tab[k] = this->GetDistantProcs(k);
    return make_shared<ParallelDofs>(this->GetCommunicator(), move(tab));
  }

  void ParallelDofs :: EnumerateGlobally (shared_ptr<BitArray> freedofs, 
					  Array<int> & global_nums,
					  int & num_glob_dofs) const
  {
    int ntasks = comm.Size();
    int id = comm.Rank();

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


#else 
namespace ngla
{
  
  shared_ptr<BitArray> ParallelDofs :: MasterDofs () const
  {
    auto ismaster = make_shared<BitArray> (ndof);
    ismaster->Set();
    return ismaster;
  }
  
  void ParallelDofs ::
  EnumerateGlobally (shared_ptr<BitArray> freedofs, Array<int> & globnum, int & num_glob_dofs) const
  {
    if (!freedofs)
      {
        for (int i = 0; i < globnum.Size(); i++)
          globnum[i] = i;
        num_glob_dofs = globnum.Size();
      }
    else
      {
        int cnt = 0;
        for (int i = 0; i < globnum.Size(); i++)
          if (freedofs->Test(i))
            globnum[i] = cnt++;
          else
            globnum[i] = -1;            
        num_glob_dofs = cnt;
      }
  }
  
}
#endif
