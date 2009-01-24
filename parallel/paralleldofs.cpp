#ifdef PARALLEL

#include <mystdlib.h>
#include <parallelngs.hpp>
#include <comp.hpp>

namespace ngparallel
{
  using namespace ngparallel;
  using namespace ngcomp;


  ParallelDofs :: ParallelDofs( const FESpace & afespace) : fespace ( afespace )
  {
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    ndof = 0;
    isexchangedof.SetSize(0);;
//     localexchangedof = 0;
//     distantexchangedof = 0;
    ismasterdof.SetSize(0);

    isghostdof = 0;

    nexdof.SetSize(ntasks);
    nexdof = 0;

    distndof.SetSize(ntasks);
    distndof = 0;
  }


  ParallelDofs :: ~ParallelDofs()
  {
//     delete localexchangedof;
//     delete distantexchangedof;
    delete isghostdof;

    /*
    for ( int dest = 0; dest < ntasks; dest++)
      if ( IsExchangeProc (dest) )
	delete this->mpi_t[dest];
    delete [] this->mpi_t;
    */
  }


//   void ParallelDofs :: SetDistantDofs ( const int proc, const Array<int> & dofs )
//   {
//     for ( int i = 0; i < dofs.Size(); i++)
//       SetDistantDof ( proc, i, dofs[i] );
//   }


//   void ParallelDofs :: SetDistantDof ( const int proc, const int localdof, const int distantdof )
//   { 
//     // exchangedofs[dest][i] ---- distantdof zu exchangedofs[id][i] = localdof

//     // find i s.t. exchangedofs[id][i] == localdof
//     int i = 0;
//     while ( i < nexdof[proc] )
//       {
// 	if ( (*localexchangedof)[proc][i] == localdof ) break;
// 	i++;
//       }
//     if ( i == nexdof[proc] ) 
//       {
// 	cout << "!!!!!!!! error in SetDistantDof - localdof is no exchange dof!!!!!!!!" << endl;
// 	*testout << "!!!!!!!! error in SetDistantDof - localdof is no exchange dof!!!!!!!!" << endl;
//       }
//     // set exchangedofs[dest][i] = distantdof
//     (*distantexchangedof) [proc][i] = distantdof;

//   }


//   const int ParallelDofs :: GetDistantDof ( const int proc, const int localdof ) const
//   {
//     // find i s.t. exchangedof[id][i] == localdof
//     int i = 0;
//     while ( i < nexdof[proc] )
//       {
// 	if ( (*localexchangedof)[proc][i] == localdof ) break;
// 	i++;
//       }
//     if ( i == nexdof[proc] ) 
//       {
//       cout << "!!!!!!!! error in GetDistantDof (" << id << ") - localdof is no exchange dof!!!!!!!!" << endl;
//       *testout << "!!!!!!!! error in GetDistantDof (" << id << ") for " << proc << " - localdof is no exchange dof!!!!!!!!" << endl;
//       }
//     // set exchangedofs[dest][i] = distantdof
//     return (*distantexchangedof) [proc][i];

//   }


  void ParallelDofs :: GetExchangeProcs ( const int localdof, Array<int> & procs ) const
  {
    procs.SetSize(0);
    for ( int dest = 0; dest < ntasks; dest++)
      if ( IsExchangeDof( dest, localdof ) && dest != id  )
	procs.Append (dest );
    return;
  }

  void ParallelDofs :: GetExchangeProcs ( Array<int> & procs ) const
  {
    procs.SetSize(0);
    if ( ntasks == 1 ) return;
    for ( int dest = 0; dest < ntasks; dest++)
      if (  (*sorted_exchangedof)[dest].Size() > 0 && dest != id )
	procs.Append (dest );
    return;

  }


  void ParallelDofs :: GetHOExchangeProcs ( const int localdof, Array<int> & procs ) const
  {
    procs.SetSize(0);
    for ( int dest = 1; dest < ntasks; dest++)
      if ( IsExchangeDof( dest, localdof ) && dest != id )
	procs.Append (dest );
    return;
  }

  void ParallelDofs :: GetHOExchangeProcs ( Array<int> & procs ) const
  {
    procs.SetSize(0);
    if ( ntasks == 1 ) return;
    for ( int dest = 1; dest < ntasks; dest++)
      if (  (*sorted_exchangedof)[dest].Size() > 0 && dest != id )
	procs.Append (dest );
    return;

  }

  const bool  ParallelDofs :: IsExchangeProc ( const int proc ) const
  {
    if ( id == proc )
      return false;
    if ( id == 0 || proc == 0 )
      return true;

    if ( (*sorted_exchangedof)[proc].Size() > 0 )
      return true;
    else 
      return false;
  }
  

  // high-order update for the local procs
  void ParallelDofs :: Update()
  {
    const MeshAccess & ma = fespace . GetMeshAccess();
    ndof = fespace.GetNDof();
    isexchangedof.SetSize( (ntasks+1)*ndof );
    isexchangedof.Clear();


     // ****************************
     // find ghost dofs
     // ****************************

     if ( isghostdof )
       delete isghostdof;
     isghostdof = new BitArray ( ndof );
     isghostdof->Set();

     for ( int el = 0; el < ma.GetNE(); el++ )
       {
	 if ( ma.IsGhostEl ( el ) ) continue;
	 Array<int> dofs;
	 fespace.GetDofNrs ( el, dofs );
	 for ( int i = 0; i < dofs.Size(); i++ )
	   isghostdof -> Clear ( dofs[i] );
       }



    
    
     // **********************
     // decide if dof is masterdof or not
     // **********************


    ismasterdof.SetSize ( ndof );
    ismasterdof.Set();

     for ( int dof = 0; dof < ndof; dof++ )
       {
	 if ( IsGhostDof ( dof ) )
	   ismasterdof.Clear ( dof );
       }

  }

  void ParallelDofs :: UpdateMPIType ()
  {
    /*
    this->mpi_t = new MPI_Datatype * [ntasks];

    MPI_Datatype MPI_T;
    if ( fespace . IsComplex() )
      MPI_T = MyGetMPIType<Complex>();
    else
      MPI_T = MyGetMPIType<double>();

    for ( int dest = 0; dest < ntasks; dest++ )
      {
	if ( !IsExchangeProc(dest) ) continue;
	
	FlatArray<int>  localexchangedof = GetLocalExchangeDofs(dest);
	int len_vec = (localexchangedof).Size();


	int * blocklen = new int [len_vec];

	for ( int i = 0; i < len_vec; i++ )
	  {
	    blocklen[i] = fespace.GetDimension();
	  }
	this->mpi_t[dest] = new MPI_Datatype;
	
	MPI_Type_indexed ( len_vec, blocklen, &((localexchangedof)[0]), MPI_T, this->mpi_t[dest]);    
	MPI_Type_commit ( this->mpi_t[dest] );
	delete [] blocklen; 
      }
    // delete MPI_T;
    */


    
    /*
      compute sorted_exchangedof
      these are local exhangedofs, sorted if (me < you), 
      and sorted such that distant is sorted for (me > you)
    */

//     Array<int> nexdofs(ntasks);
//     for (int i = 0; i < ntasks; i++)
//       nexdofs[i] = (*localexchangedof)[i].Size();

//     sorted_exchangedof = new Table<int> (nexdofs);
//     for (int dest = 0; dest < ntasks; dest++)
//       {
//         for (int j = 0; j < nexdofs[dest]; j++)
//           (*sorted_exchangedof)[dest][j] = (*localexchangedof)[dest][j];
//         if (id < dest)
//           BubbleSort ( (*sorted_exchangedof)[dest] );
//       }

    mpi_t.SetSize (ntasks); 

    cout << "update MPItype, id = " << id << endl;

    MPI_Datatype MPI_T;
    if ( fespace . IsComplex() )
      {
	switch ( fespace.GetDimension() )
	  {
	  case 2:
	    MPI_T = MyGetMPIType<Vec<2, Complex> >(); break;
	  case 3:
	    MPI_T = MyGetMPIType<Vec<3, Complex> >(); break;
	  case 4:
	    MPI_T = MyGetMPIType<Vec<4, Complex> >(); break;
	  default:
	    MPI_T = MyGetMPIType<Complex>(); break;
	  }
      }
    else
      {
	switch ( fespace.GetDimension() )
	  {
	  case 2:
	    MPI_T = MyGetMPIType<Vec<2> >(); break;
	  case 3:
	    MPI_T = MyGetMPIType<Vec<3> >(); break;
	  case 4:
	    MPI_T = MyGetMPIType<Vec<4> >(); break;
	  default:
	    MPI_T = MyGetMPIType<double>(); break;
	  }
      }

    for ( int dest = 0; dest < ntasks; dest++ )
      {
	if ( !IsExchangeProc(dest) ) continue;
	
        FlatArray<int> sortedexchangedof = GetSortedExchangeDofs(dest);
	int len_vec = sortedexchangedof.Size();
	if ( len_vec == 0 ) continue;

         Array<int> blocklen(len_vec);
	// das geht wirklich nicht --> MPI_Type entspricht schon Vec, blocklen = 1
//         blocklen = fespace.GetDimension();  // ich glaub das funktioniert nicht so -> offset in multiples of MPI_T
	blocklen = 1;
	
	MPI_Type_indexed ( len_vec, &blocklen[0], &sortedexchangedof[0], MPI_T, &mpi_t[dest]);
	MPI_Type_commit ( &mpi_t[dest] );
      }
   }






  void ParallelDofs ::  Print() const
  {
    if ( this == 0 ) return;

    (*testout) << endl <<  "DEGREES OF FREEDOM FOR PARALLEL SPACES" << endl << endl;
    (*testout) << "NDOF = "  << ndof << endl;
//       for ( int i = 0; i < ndof; i++ )
// 	{ 
// 	if ( IsExchangeDof(i) )
// 	  {

// 	    (*testout) << "exchange dof  " << i << ": " <<  endl;
// 	    for ( int dest = 0; dest < ntasks; dest ++)
// 	      if ( IsExchangeDof( dest, i ) )
// 		(*testout) << "   p" << dest << ": " << GetDistantDof ( dest, i ) << endl; 

// 	  }
// 	}

  }


//   void ParallelDofs ::  GetDistantDofs ( const int localdof, Array<int> & distantdofs ) const
//   { 
//     distantdofs.SetSize(ntasks);
//     // no distant dof on proc 0, as this one has only low order fespace
//     distantdofs[0] = -1;
//     for ( int i=1; i<ntasks; i++)
//       distantdofs[i] = GetDistantDof ( i, localdof );

//   }


  bool ParallelDofs :: ContainsParallelDof ( const FlatArray<int> & doftable ) const
  {
    if ( this == 0 ) { return false;}
    if ( ntasks == 1 ) return false;
    for ( int i = 0; i < doftable.Size(); i++ )
      if ( IsExchangeDof ( doftable[i] ) ) return true;
    return false;
  }



  void ParallelDofs :: UpdateCompound ( )
  {
    int base = 0;

    const CompoundFESpace & compoundspace = dynamic_cast<const CompoundFESpace & > (fespace);
    int nspaces = compoundspace.GetNSpaces();

    const MeshAccess & ma = fespace . GetMeshAccess();
    ndof = fespace.GetNDof();
    isexchangedof.SetSize( (ntasks+1)*ndof );
    isexchangedof.Clear();
  
    // **************************
    // find exchange dofs
    // **************************
  
    isexchangedof.SetSize ( ndof * (ntasks+1) );
    isexchangedof.Clear();
  
    base = 0;
    for ( int i = 0; i < nspaces; i++ )
      {
	for ( int dof = 0; dof < compoundspace[i]->GetNDof(); dof++ )
	  {
	    if ( compoundspace[i]->GetParallelDofs().IsExchangeDof(dof) ) 
	      SetExchangeDof ( dof + base );
	    for ( int dest = 0; dest < ntasks; dest++ )
	      if ( compoundspace[i]->GetParallelDofs().IsExchangeDof(dest, dof) )
		SetExchangeDof ( dest, dof + base );
	  }
	base += compoundspace[i]->GetNDof();
      }



     // ****************************
     // find ghost dofs
     // ****************************

    if ( isghostdof )
      delete isghostdof;
    isghostdof = new BitArray ( ndof );
    isghostdof->Set();
    
    base = 0;
    for ( int i = 0; i < nspaces; i++ )
      {
	for ( int dof = 0; dof < compoundspace[i]->GetNDof(); dof++ )
	  {
	    if ( ! compoundspace[i]->GetParallelDofs().IsGhostDof(dof) ) 
	      isghostdof -> Clear ( dof + base );
	  }
	base += compoundspace[i]->GetNDof();
      }


    
    
     // **********************
     // decide if dof is masterdof or not
     // **********************


    ismasterdof.SetSize ( ndof );
    ismasterdof.Set();

  base = 0;
    for ( int i = 0; i < nspaces; i++ )
      {
	for ( int dof = 0; dof < compoundspace[i]->GetNDof(); dof++ )
	  {
	    if ( ! compoundspace[i]->GetParallelDofs().IsMasterDof(dof) ) 
	      ismasterdof.Clear ( dof + base );
	  }
	base += compoundspace[i]->GetNDof();
      }

    // ****************************
    // copy localexchangedof, distantexchangedof
    // ****************************


    nexdof = 0;
    for ( int dest = 0; dest < ntasks; dest ++ )
      for ( int i = 0; i < nspaces; i++ )
	nexdof[dest] += compoundspace[i]->GetParallelDofs().GetNExDofs( dest );

//     localexchangedof = new Table<int> (nexdof);
//     distantexchangedof = new Table<int> (nexdof);
    sorted_exchangedof = new Table<int> (nexdof);

    Array<int> distbase(ntasks);
    distbase = 0;
  base = 0;
  Array<int> ii_comp(ntasks);
  ii_comp = 0;
    for ( int i = 0; i < nspaces; i++ )
      {
	for ( int dest = 0; dest < ntasks; dest++ )
	  {
	    for ( int ii = 0; ii < (*compoundspace[i]->GetParallelDofs().sorted_exchangedof)[dest].Size(); ii++ )
	      {
// 		(*localexchangedof)[dest][ii_comp[dest]] =  
// 		  (*compoundspace[i]->GetParallelDofs().localexchangedof)[dest][ii] + base;
		(*sorted_exchangedof)[dest][ii_comp[dest]] =  
		  (*compoundspace[i]->GetParallelDofs().sorted_exchangedof)[dest][ii] + base;
// 		(*distantexchangedof)[dest][ii_comp[dest]] =  
// 		  (*compoundspace[i]->GetParallelDofs().distantexchangedof)[dest][ii] + distbase[dest];
		ii_comp[dest]++;
	      }
	    distbase[dest] += compoundspace[i]->GetParallelDofs() . GetDistNDof (dest);
	  }
	base += compoundspace[i]->GetNDof();
      }


  }




}

#endif
