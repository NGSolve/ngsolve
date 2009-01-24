#ifndef FILE_PARALLELDOFS
#define FILE_PARALLELDOFS



using namespace ngcomp;

#ifndef PARALLEL
class ParallelDofs
{
  // dummy class
  ParallelDofs( const FESpace & afespace )
  { ; }
  ~ParallelDofs()
  { ; }

};

#else

class ParallelDofs
{

protected:
  ///
  const FESpace & fespace;
  /// number of degrees of freedom
  int ndof;
  /// number of degrees of freedom shared with processors
  Array<int> nexdof;
  /// is dof i local or exchange-dof??
  BitArray isexchangedof;
  /// is dof i ghost dof or not??
  BitArray * isghostdof;
  /// number of dofs on distant processes
  Array<int> distndof;
  /// mpi-datatype to send exchange dofs
  // MPI_Datatype ** mpi_t;
  Array<MPI_Datatype> mpi_t;

public:
  
  /// arrays containing local dof nrs for exchangedof i, and distant dof nrs for exchangedof i
  /// ordering is such that distantdofs are sorted (right, JS ?)
  /// can we change it, such that the local is ordered for (me < you) and distant is ordered for (me > you)  (JS) ?
  /// then both pairs have same ordering !
//   Table < int> *localexchangedof, *distantexchangedof;

  /// these are local exhangedofs, sorted if (me < you), and sorted such that distant is sorted for (me > you)
  /// computed in updatempitype.
  Table < int> *sorted_exchangedof;

  /// is this the master process??
  BitArray ismasterdof;

  ParallelDofs( const FESpace & afespace );
  ~ParallelDofs();

  int GetNDof () const
  { return ndof; }

  int GetNExDofs() const
  { return nexdof[id]; } 

  int GetNExDofs( int proc ) const
  { return nexdof[proc]; } 

  void SetNExDof ( const Array<int> & anexdof ) 
  { 
    nexdof.SetSize(ntasks);
    for ( int i = 0; i < ntasks; i++ ) 
      nexdof[i] = anexdof[i];
  }

  void SetExchangeDof ( int dof )
  {  isexchangedof.Set ( (ntasks+1)*dof + ntasks );  }

  void SetExchangeDof ( int proc, int dof )
  {  isexchangedof.Set ( (ntasks+1)*dof + proc );  }

//   void SetDistantDof ( int proc, int localdof, int distantdof );

//   void SetDistantDofs ( int proc, const Array<int> & dofs );

  void Update();

//   void UpdateLowOrder(LocalHeap & lh);

//   void UpdateLowOrder ()
//   {
//     LocalHeap lh ( 100000 );
//     UpdateLowOrder ( lh );
//   }
//   void UpdateGlobal(LocalHeap & lh);

//   void UpdateGlobal ()
//   {
//     LocalHeap lh ( 100000 );
//     UpdateGlobal ( lh );
//   }

  const inline  bool IsExchangeDof ( const int dof ) const
  {  return isexchangedof.Test( (ntasks+1)*dof + ntasks ) ;  }

  const inline bool IsExchangeDof ( const int proc, const int dof ) const
  {  return isexchangedof.Test ( (ntasks+1)*dof + proc ) ;  }

//   const inline int ExDof2LocalDof ( const int proc, const int exdof ) const
//   { return (*localexchangedof) [proc][exdof]; }
//   const inline int ExDof2DistantDof ( const int proc, const int localexdof ) const
//   { return (*distantexchangedof) [proc][localexdof]; }

//   const int GetDistantDof ( const int proc, const int localdof ) const;
  
  bool IsMasterDof ( int localdof ) const
  { return ismasterdof.Test(localdof); }

  inline void SetMasterDof ( const int localdof ) 
  { ismasterdof.Set ( localdof ); }
  
//   void GetDistantDofs ( const int localdof,  Array <int> & distantdofs ) const;

  void GetExchangeProcs ( const int localdof, Array<int> & procs ) const;

  void GetExchangeProcs ( Array<int> & procs ) const;

  void GetHOExchangeProcs ( const int localdof, Array<int> & procs ) const;

  void GetHOExchangeProcs ( Array<int> & procs ) const;

  const bool IsExchangeProc ( const int proc ) const;

  const int GetReceiveDof ( const int proc, const int recvnum ) const;

  void Print() const;

//   const FlatArray<int>  GetLocalExchangeDofs ( const int proc ) const
//   { return (*localexchangedof)[proc]; }

  const FlatArray<int>  GetSortedExchangeDofs ( const int proc ) const
  { return (*sorted_exchangedof)[proc]; }

  inline void SetGhostDof ( const int dof )
  { isghostdof -> Set(dof); }

  inline void ClearGhostDof ( const int dof )
  { isghostdof -> Clear(dof); }

  inline bool IsGhostDof ( const int dof ) const
  { return isghostdof->Test(dof); }

  bool ContainsParallelDof ( const FlatArray<int> & doftable ) const;

  void UpdateCompound (  );

  inline int GetDistNDof ( const int dest ) const { return distndof[dest]; }
  inline void SetDistNDof ( const int dest, const int dndof ) { distndof[dest] = dndof; }

  inline void ClearExchangeDof ( const int dof )
  {  isexchangedof.Clear ( (ntasks+1)*dof + ntasks );  }

  inline void ClearExchangeDof ( const int proc, const int dof )
  {  isexchangedof.Clear ( (ntasks+1)*dof + proc );  }
  inline void ClearMasterDof ( const int localdof ) 
  { ismasterdof.Clear ( localdof ); }

  MPI_Datatype MyGetMPI_Type ( int dest )
  { return mpi_t[dest]; }

  void UpdateMPIType () ;
};


#endif //PARALLEL

#endif
