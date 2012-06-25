#ifdef PARALLEL

#include <parallelngs.hpp> 
#include <comp.hpp>



namespace ngparallel
{
  using namespace ngcomp;


  MPI_Comm ngs_comm;
  /*
  int id = 0;
  int ntasks = 1;
  */

  ParallelDofs :: ParallelDofs (const MeshAccess & ama, 
				const Array<Node> & adofnodes, 
				int dim, bool iscomplex)
    : ma(ama), dofnodes(adofnodes)
  {
    static Timer timer ("ParallelDofs");
    RegionTimer reg(timer);

    int ntasks = MyMPI_GetNTasks();
    int id = MyMPI_GetId();

    ndof = dofnodes.Size();
    Array<int> distprocs;


    Array<int> nexdofs(ntasks), ndistprocs(ndof);
    nexdofs = 0;
    ndistprocs = 0;
    for (int i = 0; i < ndof; i++)
      {
	if (dofnodes[i].GetNr() == -1) continue;
	ma.GetDistantProcs (dofnodes[i], distprocs);
	ndistprocs[i] = distprocs.Size();
	for (int j = 0; j < distprocs.Size(); j++)
	  nexdofs[distprocs[j]]++;
      }

    exchangedofs = new Table<int> (nexdofs);
    dist_procs = new Table<int> (ndistprocs);
    nexdofs = 0;

    for (int i = 0; i < ndof; i++)
      {
	if (dofnodes[i].GetNr() == -1) continue;
	ma.GetDistantProcs (dofnodes[i], distprocs);
	for (int j = 0; j < distprocs.Size(); j++)
	  {
	    int dest = distprocs[j];
	    (*exchangedofs)[dest][nexdofs[dest]++] = i;
	    (*dist_procs)[i][j] = distprocs[j];
	  }
      }


    // *testout << "exchangedofs = " << endl << *exchangedofs << endl;

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
    
    if ( dim != 1)
      {
	MPI_Datatype htype;
	MPI_Type_contiguous (dim, mpi_type, &htype);
	mpi_type = htype;
      }

    for (int dest = 0; dest < ntasks; dest++ )
      {
	if ( !IsExchangeProc(dest) ) continue;
	
        FlatArray<int> sortedexchangedof = GetExchangeDofs(dest);
	int len_vec = sortedexchangedof.Size();
	if ( len_vec == 0 ) continue;

	Array<int> blocklen(len_vec);
	blocklen = 1;

	MPI_Type_indexed ( len_vec, &blocklen[0], &sortedexchangedof[0], 
			   mpi_type, &mpi_t[dest]);
	MPI_Type_commit ( &mpi_t[dest] );
      }
  }

  ParallelDofs :: ~ParallelDofs()
  {
    ;
  }


  int ParallelDofs :: GetNDofGlobal () const
  {
    if (GetNTasks() == 1) return ndof;
    
    int nlocal = 0;

    for (int i = 0; i < ndof; i++)
      if (ismasterdof.Test(i)) nlocal++;

    return MyMPI_AllReduce (nlocal);

    // MPI_Allreduce (&nlocal, &nglobal, 1,  MPI_INT, MPI_SUM, ngs_comm);
    // return nglobal;
  }






}

#endif
