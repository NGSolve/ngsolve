#ifndef FILE_PARALLELDOFS
#define FILE_PARALLELDOFS


using namespace ngcomp;

#ifdef PARALLEL
 
class ParallelDofs
{

protected:
  ///
  const FESpace & fespace;

  /// number of degrees of freedom shared with processors
  Array<int> nexdof;

  /// is dof i local or exchange-dof??
  BitArray isexchangedof;

  /// is dof i ghost dof or not??
  BitArray * isghostdof;

  /// number of dofs on distant processes
  Array<int> distndof;


  /// mpi-datatype to send exchange dofs
  Array<MPI_Datatype> mpi_t;

  friend class ngcomp::FESpace;
  /// these are local exhangedofs, computed by FESpace
  Table<int> * sorted_exchangedof;

  /// is this the master process??
  BitArray ismasterdof;

public:

  ParallelDofs( const FESpace & afespace );
  ~ParallelDofs();

  int GetNDof () const { return fespace.GetNDof(); }
  const FESpace & GetFESpace() const { return fespace; }


  int GetNExDofs() const { return nexdof[id]; } 
  int GetNExDofs( int proc ) const { return nexdof[proc]; } 
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


  void Update();

  const inline  bool IsExchangeDof ( const int dof ) const
  {  return isexchangedof.Test( (ntasks+1)*dof + ntasks ) ;  }

  const inline bool IsExchangeDof ( const int proc, const int dof ) const
  {  return isexchangedof.Test ( (ntasks+1)*dof + proc ) ;  }
  
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

  inline int GetDistNDof (int dest) const { return distndof[dest]; }
  inline void SetDistNDof (int dest, int dndof) { distndof[dest] = dndof; }

  inline void ClearExchangeDof (int dof )
  {  isexchangedof.Clear ( (ntasks+1)*dof + ntasks );  }

  inline void ClearExchangeDof (int proc, int dof )
  {  isexchangedof.Clear ( (ntasks+1)*dof + proc );  }
  inline void ClearMasterDof (int localdof) 
  { ismasterdof.Clear ( localdof ); }

  MPI_Datatype MyGetMPI_Type ( int dest )
  { return mpi_t[dest]; }

  void UpdateMPIType () ;
};


#endif //PARALLEL

#endif
