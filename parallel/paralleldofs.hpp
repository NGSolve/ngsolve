#ifndef FILE_PARALLELDOFS
#define FILE_PARALLELDOFS


namespace ngfem
{
  class Node;
}

namespace ngcomp
{
  class MeshAccess;
}

namespace ngparallel
{
  using ngfem::Node;
  using ngcomp::MeshAccess;


#ifdef PARALLEL



  class ParallelMeshDofs : public ParallelDofs
  {
    const MeshAccess & ma;
    Array<Node> dofnodes;
  public:
    ParallelMeshDofs (const MeshAccess & ama, const Array<Node> & adofnodes, 
		      int dim = 1, bool iscomplex = false);

    const MeshAccess & GetMeshAccess() const { return ma; }
    const Array<Node> & GetDofNodes() const { return dofnodes; }
  };
  


  template <typename T>
  void ReduceDofData (FlatArray<T> data, MPI_Op op, const ParallelDofs & pardofs)
  {
    if (&pardofs == NULL) return;
    MPI_Comm comm = pardofs.GetCommunicator();
    int ntasks = MyMPI_GetNTasks (comm);
    if (ntasks <= 1) return;

    DynamicTable<T> dist_data(ntasks);

    for (int i = 0; i < pardofs.GetNDof(); i++)
      if (!pardofs.IsMasterDof(i))
	{
	  FlatArray<int> distprocs = pardofs.GetDistantProcs (i);
	  int master = ntasks;
	  for (int j = 0; j < distprocs.Size(); j++)
	    master = min (master, distprocs[j]);
	  dist_data.Add (master, data[i]);
	}
    
    Array<int> nsend(ntasks), nrecv(ntasks);
    for (int i = 0; i < ntasks; i++)
      nsend[i] = dist_data[i].Size();
    MyMPI_AllToAll (nsend, nrecv, comm);
    Table<T> recv_data(nrecv);

    Array<MPI_Request> requests;
    for (int i = 0; i < ntasks; i++)
      {
	if (nsend[i])
	  requests.Append (MyMPI_ISend (dist_data[i], i, MPI_TAG_SOLVE, comm));
	if (nrecv[i])
	  requests.Append (MyMPI_IRecv (recv_data[i], i, MPI_TAG_SOLVE, comm));
      }

    MyMPI_WaitAll (requests);

    Array<int> cnt(ntasks);
    cnt = 0;
    
    MPI_Datatype type = MyGetMPIType<T>();
    for (int i = 0; i < pardofs.GetNDof(); i++)
      if (pardofs.IsMasterDof(i))
	{
	  FlatArray<int> distprocs = pardofs.GetDistantProcs (i);
	  for (int j = 0; j < distprocs.Size(); j++)
	    MPI_Reduce_local (&recv_data[distprocs[j]][cnt[distprocs[j]]++], 
			      &data[i], 1, type, op);
	}
  }    





  template <typename T>
  void ScatterDofData (FlatArray<T> data, const ParallelDofs & pardofs)
  {
    if (&pardofs == NULL) return;
    MPI_Comm comm = pardofs.GetCommunicator();

    int ntasks = MyMPI_GetNTasks (comm);
    // int id = MyMPI_GetId (comm);
    if (ntasks <= 1) return;

    DynamicTable<T> dist_data(ntasks);
    for (int i = 0; i < pardofs.GetNDof(); i++)
      if (pardofs.IsMasterDof(i))
	{
	  FlatArray<int> distprocs = pardofs.GetDistantProcs (i);
	  for (int j = 0; j < distprocs.Size(); j++)
	    dist_data.Add (distprocs[j], data[i]);
	}

    Array<int> nsend(ntasks), nrecv(ntasks);
    for (int i = 0; i < ntasks; i++)
      nsend[i] = dist_data[i].Size();

    MyMPI_AllToAll (nsend, nrecv, comm);

    Table<T> recv_data(nrecv);

    Array<MPI_Request> requests;
    for (int i = 0; i < ntasks; i++)
      {
	if (nsend[i])
	  requests.Append (MyMPI_ISend (dist_data[i], i, MPI_TAG_SOLVE, comm));
	if (nrecv[i])
	  requests.Append (MyMPI_IRecv (recv_data[i], i, MPI_TAG_SOLVE, comm));
      }

    MyMPI_WaitAll (requests);

    Array<int> cnt(ntasks);
    cnt = 0;
    
    for (int i = 0; i < pardofs.GetNDof(); i++)
      if (!pardofs.IsMasterDof(i))
	{
	  FlatArray<int> distprocs = pardofs.GetDistantProcs (i);
	  
	  int master = ntasks;
	  for (int j = 0; j < distprocs.Size(); j++)
	    master = min (master, distprocs[j]);
	  data[i] = recv_data[master][cnt[master]++];
	}
  }    

  template <typename T>
  void AllReduceDofData (FlatArray<T> data, MPI_Op op, 
			 const ParallelDofs & pardofs)
  {
    if (&pardofs == NULL) return;
    ReduceDofData (data, op, pardofs);
    ScatterDofData (data, pardofs);
  }




#else


  class ParallelMeshDofs : public ParallelDofs 
  {
    ParallelDofs (const MeshAccess & ama, const Array<Node> & adofnodes, 
		  int dim = 1, bool iscomplex = false)
    { ndof = adofnodes.Size(); }
  };

  
#endif //PARALLEL














}



#endif
