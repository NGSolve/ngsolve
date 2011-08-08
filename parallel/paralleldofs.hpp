#ifndef FILE_PARALLELDOFS
#define FILE_PARALLELDOFS

namespace ngfem
{
  class Node;
}


using namespace ngfem;

namespace ngcomp
{
  class FESpace;
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
    // const FESpace * fes;
    /// these are local exhangedofs, computed by FESpace
    Table<int> * exchangedofs;

    /// mpi-datatype to send exchange dofs
    Array<MPI_Datatype> mpi_t;

    /// is this the master process??
    BitArray ismasterdof;

    Array<Node> dofnodes;

  public:
    // ParallelDofs (int andof, Table<int> * exdofs, const FESpace * afes);

    /*
    ParallelDofs (const MeshAccess & ma, const Array<Node> & adofnodes, const FESpace * afes = NULL);
    */

    ParallelDofs (const MeshAccess & ama, const Array<Node> & adofnodes, int dim = 1, bool iscomplex = false);

    const Array<Node> & GetDofNodes() const 
    { return dofnodes; }

    const FlatArray<int>  GetExchangeDofs (int proc) const
    { return (*exchangedofs)[proc]; }

    bool IsMasterDof ( int localdof ) const
    { return ismasterdof.Test(localdof); }

    int GetNDof () const { return ndof; }

    int GetNDofGlobal () const;

    // used by masterinverse
    // const FESpace & GetFESpace() const { return *fes; }

    const MeshAccess & GetMeshAccess () const { return ma; }

    bool IsExchangeProc ( int proc ) const
    { return (*exchangedofs)[proc].Size() != 0; }

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

    // ParallelDofs (const MeshAccess & ma, const Array<Node> & dofnodes, const FESpace * afes = NULL)
    // { ndof = dofnodes.Size(); }


    ParallelDofs (const MeshAccess & ama, const Array<Node> & adofnodes, int dim = 1, bool iscomplex = false)
    { ndof = adofnodes.Size(); }

    int GetNDofGlobal () const { return ndof; }

  };

#endif //PARALLEL
}



#endif
