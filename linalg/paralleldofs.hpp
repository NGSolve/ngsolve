#ifndef FILE_PARALLELDOFS
#define FILE_PARALLELDOFS

/**************************************************************************/
/* File:   paralleldofs.hpp                                               */
/* Author: Joachim Schoeberl                                              */
/* Date:   July 2011                                                      */
/**************************************************************************/



namespace ngla
{

#ifdef PARALLEL

  /**
     Handles the distribution of degrees of freedom for vectors and matrices
   */
  class ParallelDofs
  {
  protected:
    /// the communicator 
    MPI_Comm comm;
    
    /// local ndof
    int ndof;

    /// proc 2 dofs
    Table<int> * exchangedofs;
    
    /// dof 2 procs
    Table<int> * dist_procs;
    
    /// mpi-datatype to send exchange dofs
    Array<MPI_Datatype> mpi_t;
    
    /// am I the master process ?
    BitArray ismasterdof;
    
  public:
    /**
       setup parallel dofs.
       Table adist_procs must provide the distant processes for each dof.
       table 
     */
    ParallelDofs (MPI_Comm acomm, Table<int> * adist_procs, 
		  int dim = 1, bool iscomplex = false);

    /*
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
    */
    
    virtual ~ParallelDofs()  
    { 
      delete exchangedofs;
    }

    int GetNTasks() const { return exchangedofs->Size(); }

    FlatArray<int> GetExchangeDofs (int proc) const
    { return (*exchangedofs)[proc]; }

    FlatArray<int> GetDistantProcs (int dof) const
    { return (*dist_procs)[dof]; }

    bool IsMasterDof ( int localdof ) const
    { return ismasterdof.Test(localdof); }

    int GetNDofLocal () const { return ndof; }

    int GetNDofGlobal () const
    {
      if (GetNTasks() == 1) return ndof;
      int nlocal = 0;
      for (int i = 0; i < ndof; i++)
	if (ismasterdof.Test(i)) nlocal++;
      return MyMPI_AllReduce (nlocal);
    }

    bool IsExchangeProc (int proc) const
    { return (*exchangedofs)[proc].Size() != 0; }

    MPI_Datatype MyGetMPI_Type (int dest) const
    { return mpi_t[dest]; }

    MPI_Comm GetCommunicator () const { return comm; }

    int GetMasterProc (int dof) const
    {
      FlatArray<int> procs = GetDistantProcs(dof);
      int m = MyMPI_GetId(comm);
      for (int j = 0; j < procs.Size(); j++)
	m = min2(procs[j], m);
      return m;
    }

    void EnumerateGlobally (const BitArray * freedofs, Array<int> & globnum, int & num_glob_dofs) const;


    template <typename T>
    void ReduceDofData (FlatArray<T> data, MPI_Op op) const;

    template <typename T>
    void ScatterDofData (FlatArray<T> data) const;
    
    template <typename T>
    void AllReduceDofData (FlatArray<T> data, MPI_Op op) const
    {
      if (this == NULL) return;
      ReduceDofData (data, op);
      ScatterDofData (data);
    }
    
  };

#else
  class ParallelDofs 
  {
  protected:
    int ndof;
    
  public:
    
    int GetNDofGlobal () const { return ndof; }

    template <typename T>
    void ReduceDofData (FlatArray<T> data, MPI_Op op) const { ; }

    template <typename T>
    void ScatterDofData (FlatArray<T> data) const { ; } 
    
    template <typename T>
    void AllReduceDofData (FlatArray<T> data, MPI_Op op) const { ; }
  };
  
#endif



  template <typename T>
  void ReduceDofData (FlatArray<T> data, MPI_Op op, const ParallelDofs & pardofs)
  {
    pardofs.ReduceDofData(data, op);
  }

  template <typename T>
  void ScatterDofData (FlatArray<T> data, const ParallelDofs & pardofs)
  {
    pardofs.ScatterDofData (data);
  }

  template <typename T>
  void AllReduceDofData (FlatArray<T> data, MPI_Op op, 
			 const ParallelDofs & pardofs)
  { 
    pardofs.AllReduceDofData (data, op);
  }




#ifdef PARALLEL

  template <typename T>
  void ParallelDofs::ReduceDofData (FlatArray<T> data, MPI_Op op) const
  {
    if (this == NULL) return;
    int ntasks = GetNTasks ();
    if (ntasks <= 1) return;

    DynamicTable<T> dist_data(ntasks);

    for (int i = 0; i < GetNDofLocal(); i++)
      if (!IsMasterDof(i))
	{
	  FlatArray<int> distprocs = GetDistantProcs (i);
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
    for (int i = 0; i < GetNDofLocal(); i++)
      if (IsMasterDof(i))
	{
	  FlatArray<int> distprocs = GetDistantProcs (i);
	  for (int j = 0; j < distprocs.Size(); j++)
	    MPI_Reduce_local (&recv_data[distprocs[j]][cnt[distprocs[j]]++], 
			      &data[i], 1, type, op);
	}
  }    





  template <typename T>
  void ParallelDofs :: ScatterDofData (FlatArray<T> data) const
  {
    if (this == NULL) return;
    MPI_Comm comm = GetCommunicator();

    int ntasks = MyMPI_GetNTasks (comm);
    if (ntasks <= 1) return;

    DynamicTable<T> dist_data(ntasks);
    for (int i = 0; i < GetNDofLocal(); i++)
      if (IsMasterDof(i))
	{
	  FlatArray<int> distprocs = GetDistantProcs (i);
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
    
    for (int i = 0; i < GetNDofLocal(); i++)
      if (!IsMasterDof(i))
	{
	  FlatArray<int> distprocs = GetDistantProcs (i);
	  
	  int master = ntasks;
	  for (int j = 0; j < distprocs.Size(); j++)
	    master = min (master, distprocs[j]);
	  data[i] = recv_data[master][cnt[master]++];
	}
  }    


#endif //PARALLEL


}



#endif
