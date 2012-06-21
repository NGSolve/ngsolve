#ifndef FILE_PARALLELDOFS
#define FILE_PARALLELDOFS

namespace ngfem
{
  class Node;
}


using namespace ngfem;

namespace ngcomp
{
  class MeshAccess;
}


namespace ngparallel
{
  
  



#ifdef PARALLEL
  class ParallelDofs
  {
  protected:
    int ndof;
    const MeshAccess & ma;

    /// these are local exhangedofs
    Table<int> * exchangedofs;

    /// mpi-datatype to send exchange dofs
    Array<MPI_Datatype> mpi_t;

    /// is this the master process ?
    BitArray ismasterdof;

    Array<Node> dofnodes;

  public:

    ParallelDofs (const MeshAccess & ama, const Array<Node> & adofnodes, int dim = 1, bool iscomplex = false);

    const Array<Node> & GetDofNodes() const 
    { return dofnodes; }

    int GetNTasks() const { return exchangedofs->Size(); }

    const FlatArray<int>  GetExchangeDofs (int proc) const
    { return (*exchangedofs)[proc]; }

    bool IsMasterDof ( int localdof ) const
    { return ismasterdof.Test(localdof); }

    int GetNDof () const { return ndof; }

    int GetNDofGlobal () const;

    const MeshAccess & GetMeshAccess () const { return ma; }

    bool IsExchangeProc ( int proc ) const
    { return (*exchangedofs)[proc].Size() != 0; }

    MPI_Datatype MyGetMPI_Type ( int dest )
    { return mpi_t[dest]; }

    MPI_Comm GetCommunicator () const { return ngs_comm; }
  };




#else

  class ParallelDofs 
  {
    int ndof;

  public:

    ParallelDofs (const MeshAccess & ama, const Array<Node> & adofnodes, 
		  int dim = 1, bool iscomplex = false)
    { ndof = adofnodes.Size(); }

    int GetNDofGlobal () const { return ndof; }

  };



  
#endif //PARALLEL
}



#endif
