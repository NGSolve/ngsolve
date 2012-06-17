#ifdef PARALLEL

#include <comp.hpp>
#include <parallelngs.hpp> 



namespace ngparallel
{
  using namespace ngcomp;


  MPI_Comm ngs_comm;
  int id = 0;
  int ntasks = 1;
  bool working_proc = true;




  ParallelDofs :: ParallelDofs (const MeshAccess & ama, 
				const Array<Node> & adofnodes, int dim, bool iscomplex)
    : ma(ama), dofnodes(adofnodes) 
  {
    // fes = NULL;
    ndof = dofnodes.Size();
    Array<int[2]> distnums;

    Array<int> nexdofs(ntasks);
    nexdofs = 0;
    for (int i = 0; i < ndof; i++)
      {
	if (dofnodes[i].GetNr() == -1) continue;
	ma.GetDistantNodeNums (dofnodes[i], distnums);
	for (int j = 0; j < distnums.Size(); j++)
	  {
	    int dest = distnums[j][0];
	    if (dest == 0) continue;
	    nexdofs[dest]++;
	  }
      }

    exchangedofs = new Table<int> (nexdofs);
    nexdofs = 0;

    for (int i = 0; i < ndof; i++)
      {
	if (dofnodes[i].GetNr() == -1) continue;
	ma.GetDistantNodeNums (dofnodes[i], distnums);
	for (int j = 0; j < distnums.Size(); j++)
	  {
	    int dest = distnums[j][0];
	    if (dest == 0) continue;
	    (*exchangedofs)[dest][nexdofs[dest]++] = i;
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











  int ParallelDofs :: GetNDofGlobal () const
  {
    if (ntasks == 1) return ndof;
    
    int nlocal = 0, nglobal;

    if (id > 0)
      for (int i = 0; i < ndof; i++)
	if (ismasterdof.Test(i)) nlocal++;

    MPI_Allreduce (&nlocal, &nglobal, 1,  MPI_INT, MPI_SUM, ngs_comm);
    return nglobal;
  }
}

#endif
