#ifndef FILE_PARALLELDOFS
#define FILE_PARALLELDOFS


namespace ngparallel
{

#ifdef PARALLEL
  class ParallelDofs
  {
  protected:
    int ndof;
    const FESpace * fes;
    /// these are local exhangedofs, computed by FESpace
    Table<int> * sorted_exchangedof;

    /// mpi-datatype to send exchange dofs
    Array<MPI_Datatype> mpi_t;

    /// is this the master process??
    BitArray ismasterdof;

  public:
    ParallelDofs (int andof, Table<int> * exdofs, const FESpace * afes);
    
    const FlatArray<int>  GetSortedExchangeDofs (int proc) const
    { return (*sorted_exchangedof)[proc]; }

    bool IsMasterDof ( int localdof ) const
    { return ismasterdof.Test(localdof); }

    int GetNDof () const { return ndof; }

    int GetNDofGlobal () const;

    // should not be used
    const FESpace & GetFESpace() const { return *fes; }

    bool IsExchangeProc ( int proc ) const
    { return (*sorted_exchangedof)[proc].Size() != 0; }

    MPI_Datatype MyGetMPI_Type ( int dest )
    { return mpi_t[dest]; }
  };

#else

  class ParallelDofs 
  {
    int ndof;
    const FESpace * fes;

  public:
    ParallelDofs (int andof, Table<int> * exdofs, const FESpace * afes)
      : ndof(andof), fes(afes) 
    { ; }
    int GetNDofGlobal () const { return ndof; }
  };

#endif //PARALLEL
}



#endif
