#ifndef FILE_PARALLELDOFS
#define FILE_PARALLELDOFS

/**************************************************************************/
/* File:   paralleldofs.hpp                                               */
/* Author: Joachim Schoeberl                                              */
/* Date:   July 2011                                                      */
/**************************************************************************/


#include <core/mpi_wrapper.hpp>
#include <ngstd.hpp>

namespace ngla
{
  using namespace ngstd;

  
#ifdef PARALLEL

  /**
     Handles the distribution of degrees of freedom for vectors and matrices
   */
  class ParallelDofs
  {
  protected:
    /// the communicator 
    NgMPI_Comm comm;
    
    /// local ndof
    size_t ndof;

    /// global ndof
    size_t global_ndof;

    /// proc 2 dofs
    Table<int> exchangedofs;
    
    /// dof 2 procs
    Table<int> dist_procs;

    /// all procs with connected dofs
    Array<int> all_dist_procs;
    
    /// mpi-datatype to send exchange dofs
    Array<NG_MPI_Datatype> mpi_t;
    
    /// am I the master process ?
    BitArray ismasterdof;

    /// entry-size
    int es;
    bool complex;
    
  public:
    /**
       setup parallel dofs.
       Table adist_procs must provide the distant processes for each dof.
       table 
     */
    ParallelDofs (NgMPI_Comm acomm, Table<int> && adist_procs, 
		  int dim = 1, bool iscomplex = false);

    shared_ptr<ParallelDofs> SubSet (shared_ptr<BitArray> take_dofs) const;

    shared_ptr<ParallelDofs> Range (IntRange range) const;
      
    virtual ~ParallelDofs();

    int GetNTasks() const { return exchangedofs.Size(); }

    FlatArray<int> GetExchangeDofs (int proc) const
    { return exchangedofs[proc]; }

    FlatArray<int> GetDistantProcs (int dof) const
    { return dist_procs[dof]; }

    FlatArray<int> GetDistantProcs () const
    { return all_dist_procs; }

    bool IsMasterDof (size_t localdof) const
    { return ismasterdof.Test(localdof); }
    
    const auto & MasterDofs () const
    { return ismasterdof; }
    
    size_t GetNDofLocal () const { return ndof; }

    size_t GetNDofGlobal () const { return global_ndof; }

    bool IsExchangeProc (int proc) const
    { return exchangedofs[proc].Size() != 0; }

    NG_MPI_Datatype GetMPI_Type (int dest) const
    { return mpi_t[dest]; }

    const NgMPI_Comm & GetCommunicator () const { return comm; }

    int GetEntrySize () const { return es; }
    bool IsComplex () const { return complex; }
    
    int GetMasterProc (int dof) const
    {
      auto dps = GetDistantProcs(dof);
      return (dps.Size() > 0) ? min2(comm.Rank(), dps[0]) : comm.Rank();
    }

    void EnumerateGlobally (shared_ptr<BitArray> freedofs, Array<int> & globnum, int & num_glob_dofs) const;


    template <typename T>
    void ReduceDofData (FlatArray<T> data, NG_MPI_Op op) const;

    template <typename T>
    void ScatterDofData (FlatArray<T> data) const;
    
    template <typename T>
    void AllReduceDofData (FlatArray<T> data, NG_MPI_Op op) const
    {
      // if (this == NULL)  // illformed C++, shall get rid of this
        // throw Exception("AllReduceDofData for null-object");
      ReduceDofData (data, op);
      ScatterDofData (data);
    }

  };

#else

  class ParallelDofs 
  {
  protected:
    int ndof;
    int es;
    bool complex;
    BitArray masterdofs;
  public:
    ParallelDofs (NgMPI_Comm acomm, Table<int> && adist_procs, 
		  int dim = 1, bool iscomplex = false)
      : ndof(adist_procs.Size()), es(dim), complex(iscomplex)
    { ; }
    
    int GetNDofLocal () const { return ndof; }
    int GetNDofGlobal () const { return ndof; }
    
    int GetNTasks() const { return 1; }
    NgMPI_Comm GetCommunicator () const { return NgMPI_Comm(NG_MPI_COMM_WORLD); }

    shared_ptr<ParallelDofs> Range (IntRange range) const
    { return nullptr; }
    
    FlatArray<int> GetExchangeDofs (int proc) const
    { return FlatArray<int> (0, nullptr); }

    FlatArray<int> GetDistantProcs (int dof) const
    { return FlatArray<int> (0, nullptr); }

    FlatArray<int> GetDistantProcs () const
    { return FlatArray<int> (0, nullptr); }      

    bool IsMasterDof (size_t localdof) const
    { return true; }

    const BitArray & MasterDofs () const;
    void EnumerateGlobally (shared_ptr<BitArray> freedofs, Array<int> & globnum, int & num_glob_dofs) const;
    
    template <typename T>
    void ReduceDofData (FlatArray<T> data, NG_MPI_Op op) const { ; }

    template <typename T>
    void ScatterDofData (FlatArray<T> data) const { ; } 
    
    template <typename T>
    void AllReduceDofData (FlatArray<T> data, NG_MPI_Op op) const { ; }

    int GetEntrySize () const { return es; }
    bool IsComplex () const { return complex; }
  };
  
#endif
  

  
  template <typename T>
  [[deprecated("use pardofs.ReduceDofData")]]
  void ReduceDofData (FlatArray<T> data, NG_MPI_Op op, const shared_ptr<ParallelDofs> & pardofs)
  {
    if (pardofs)
      pardofs->ReduceDofData(data, op);
  }
  
  template <typename T>
  [[deprecated("use pardofs.ScatterDofData")]]   
  void ScatterDofData (FlatArray<T> data, const shared_ptr<ParallelDofs> & pardofs)
  {
    if (pardofs)
      pardofs->ScatterDofData (data);
  }
  
  template <typename T>
  [[deprecated("use pardofs.AllReduceDofData")]]
  void AllReduceDofData (FlatArray<T> data, NG_MPI_Op op, 
			 const shared_ptr<ParallelDofs> & pardofs)
  {
    if (pardofs)
      pardofs->AllReduceDofData (data, op);
  }
  


#ifdef PARALLEL

  template <typename T>
  void ParallelDofs :: ReduceDofData (FlatArray<T> data, NG_MPI_Op op) const
  {
    static Timer t0("ParallelDofs :: ReduceDofData");
    RegionTimer rt(t0);

    auto comm = GetCommunicator();
    int ntasks = comm.Size();
    int rank = comm.Rank();
    if (ntasks <= 1) return;


    Array<int> nsend(ntasks), nrecv(ntasks);
    nsend = 0;
    nrecv = 0;

    /** Count send/recv size **/
    for (int i = 0; i < GetNDofLocal(); i++)
      if (auto dps = GetDistantProcs(i); dps.Size())
        {
          if (rank < dps[0])
            for (auto p : dps)
              nrecv[p]++;
          else
            nsend[dps[0]]++;
        }
    
    Table<T> send_data(nsend);
    Table<T> recv_data(nrecv);

    /** Fill send_data **/
    nsend = 0;
    for (int i = 0; i < GetNDofLocal(); i++)
      if (auto dps = GetDistantProcs(i); dps.Size())
        if (rank > dps[0])
          send_data[dps[0]][nsend[dps[0]]++] = data[i];

    NgMPI_Requests send_requests; 
    NgMPI_Requests recv_requests; 
    for (int i = 0; i < ntasks; i++)
      {
	if (nsend[i])
	  send_requests += comm.ISend(send_data[i], i, NG_MPI_TAG_SOLVE);
	if (nrecv[i])
	  recv_requests += comm.IRecv(recv_data[i], i, NG_MPI_TAG_SOLVE);
      }

    Array<int> cnt(ntasks);
    cnt = 0;
    
    NG_MPI_Datatype type = GetMPIType<T>();

    recv_requests.WaitAll();
    for (int i = 0; i < GetNDofLocal(); i++)
      if (IsMasterDof(i))
        for (auto p : GetDistantProcs (i))
          NG_MPI_Reduce_local (&recv_data[p][cnt[p]++], &data[i], 1, type, op);

    send_requests.WaitAll();    
  }    





  template <typename T>
  void ParallelDofs :: ScatterDofData (FlatArray<T> data) const
  {
    static Timer t0("ParallelDofs :: ScatterDofData");
    RegionTimer rt(t0);

    NgMPI_Comm comm = GetCommunicator();
    int ntasks = comm.Size();
    int rank = comm.Rank();
    if (ntasks <= 1) return;

    Array<int> nsend(ntasks), nrecv(ntasks);
    nsend = 0;
    nrecv = 0;

    /** Count send/recv size **/
    for (int i = 0; i < GetNDofLocal(); i++) 
      if (auto dps = GetDistantProcs(i); dps.Size() > 0)
        {
          if (rank < dps[0])
            for (auto p : dps)
              nsend[p]++;
          else
            nrecv[dps[0]]++;
        }
    
    Table<T> send_data(nsend);
    Table<T> recv_data(nrecv);

    /** Fill send_data **/
    nsend = 0;
    for (int i = 0; i < GetNDofLocal(); i++) 
      if (auto dps = GetDistantProcs(i); dps.Size() > 0)
        if (rank < dps[0])
          for (auto p : dps)
            send_data[p][nsend[p]++] = data[i];
    
    NgMPI_Requests requests;
    for (int i = 0; i < ntasks; i++)
      {
	if (nsend[i])
	  requests += comm.ISend (send_data[i], i, NG_MPI_TAG_SOLVE);
	if (nrecv[i])
	  requests += comm.IRecv (recv_data[i], i, NG_MPI_TAG_SOLVE);
      }
    requests.WaitAll();

    Array<int> cnt(ntasks);
    cnt = 0;

    for (int i = 0; i < GetNDofLocal(); i++)
      if (!IsMasterDof(i))
	{
	  int master = GetDistantProcs (i)[0];
	  data[i] = recv_data[master][cnt[master]++];
	}
  }    

#endif //PARALLEL

}

#endif
