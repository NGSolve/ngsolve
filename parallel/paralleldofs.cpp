#ifdef PARALLEL


#include <comp.hpp>
#include <parallelngs.hpp> 


namespace ngparallel
{
  MPI_Group MPI_HIGHORDER_WORLD;
  MPI_Comm MPI_HIGHORDER_COMM;


  using namespace ngcomp;

  ParallelDofs :: ParallelDofs (int andof, Table<int> * exdofs, const FESpace * afes)
    : ndof(andof), fes(afes), sorted_exchangedof(exdofs), ismasterdof(ndof)
  {
    ismasterdof.Set();
    for (int i = 0; i < id; i++)
      {
	FlatArray<int> ex = (*sorted_exchangedof)[i];
	for (int j = 0; j < ex.Size(); j++)
	  ismasterdof.Clear (ex[j]);
      }

    mpi_t.SetSize (ntasks); 

    MPI_Datatype mpi_type = 
      ( fes -> IsComplex() ) ? 
      MyGetMPIType<Complex>() : MyGetMPIType<double>(); 
    
    if ( fes -> GetDimension() != 1)
      {
	MPI_Datatype htype;
	MPI_Type_contiguous (fes -> GetDimension(), mpi_type, &htype);
	mpi_type = htype;
      }

    for (int dest = 0; dest < ntasks; dest++ )
      {
	if ( !IsExchangeProc(dest) ) continue;
	
        FlatArray<int> sortedexchangedof = GetSortedExchangeDofs(dest);
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

    MPI_Allreduce (&nlocal, &nglobal, 1,  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return nglobal;
  }

  
  /*
  void ParallelDofs :: UpdateMPIType () 
  {
    ;
  }

  void ParallelDofs :: Update()
  {
    ;
  }
  void ParallelDofs :: Print() const
  {
    ;
  }
  */


#ifdef VERSION1

  ParallelDofs :: ParallelDofs(const FESpace & afespace) : fespace ( afespace )
  {
    isexchangedof.SetSize(0);;
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
    if (id == proc) return false;
    
    /*
    if ( id == 0 || proc == 0 )
      return true;
    */

    return ((*sorted_exchangedof)[proc].Size() > 0);
  }
  

  // high-order update for the local procs
  void ParallelDofs :: Update()
  {
    *testout << "update parallel dofs" << endl;
    const MeshAccess & ma = fespace.GetMeshAccess();
    int ndof = fespace.GetNDof();

    isexchangedof.SetSize( (ntasks+1)*ndof );
    isexchangedof.Clear();


    // ****************************
    // find ghost dofs
    // ****************************

    delete isghostdof;
    isghostdof = new BitArray ( ndof );
    isghostdof->Set();

    for (int el = 0; el < ma.GetNE(); el++ )
      {
	if ( ma.IsGhostEl ( el ) ) continue;
	Array<int> dofs;
	fespace.GetDofNrs ( el, dofs );
	for (int i = 0; i < dofs.Size(); i++ )
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
    mpi_t.SetSize (ntasks); 

    MPI_Datatype mpi_type = 
      ( fespace . IsComplex() ) ? 
      MyGetMPIType<Complex>() : MyGetMPIType<double>(); 
    
    if (fespace.GetDimension() != 1)
      {
	MPI_Datatype htype;
	MPI_Type_contiguous (fespace.GetDimension(), mpi_type, &htype);
	mpi_type = htype;
      }

    for (int dest = 0; dest < ntasks; dest++ )
      {
	if ( !IsExchangeProc(dest) ) continue;
	
        FlatArray<int> sortedexchangedof = GetSortedExchangeDofs(dest);
	int len_vec = sortedexchangedof.Size();
	if ( len_vec == 0 ) continue;

	Array<int> blocklen(len_vec);
	blocklen = 1;

	MPI_Type_indexed ( len_vec, &blocklen[0], &sortedexchangedof[0], 
			   mpi_type, &mpi_t[dest]);
	MPI_Type_commit ( &mpi_t[dest] );
      }


    /*
    // non-overlapping types
    
    mpi_t_no.SetSize (0);
    comm_no.SetSize (0);

    int ndof = fespace.GetNDof();
    Array<int> cnt(ndof);
    cnt = 0;
    for (int i = 1; i < sorted_exchangedof->Size(); i++)
      if (i != id)
	{
	  FlatArray<int> exdofs = (*sorted_exchangedof)[i];
	  for (int j = 0; j < exdofs.Size(); j++)
	    cnt[exdofs[j]]++;
	}
    Table<int> dof2proc (cnt);
    cnt = 0;
    for (int i = 1; i < sorted_exchangedof->Size(); i++)
      if (i != id)
	{
	  FlatArray<int> exdofs = (*sorted_exchangedof)[i];
	  for (int j = 0; j < exdofs.Size(); j++)
	    dof2proc[exdofs[j]][cnt[exdofs[j]]++] = i;
	}
    
    *testout << "dof2proc = " << endl << dof2proc << endl;

    Array<Array<int>*> groups;
    Array<int> groupnr(ndof);
    groupnr = -1;
    for (int i = 0; i < ndof; i++)
      {
	FlatArray<int> procs = dof2proc[i];
	if (!procs.Size()) continue;

	int found = -1;
	for (int j = 0; j < groups.Size(); j++)
	  {
	    if ( procs == *groups[j] )
	      {
		found = j;
		break;
	      }
	  }
	if (found == -1)
	  {
	    found = groups.Size();
	    groups.Append (new Array<int>);
	    for (int k = 0; k < procs.Size(); k++)
	      groups.Last() -> Append (procs[k]);
	  }
	groupnr[i] = found;
      }

    for (int i = 0; i < groups.Size(); i++)
      {
	groups[i] -> Append (id);
	BubbleSort (*groups[i]);
      }

    *testout << "groups:" << endl;
    for (int i = 0; i < groups.Size(); i++)
      *testout << "group(" << i << ") has procs:" << endl << *groups[i] << endl;
    *testout << "groupnr = " << endl << groupnr << endl;

    mpi_t_no.SetSize (groups.Size());
    comm_no.SetSize  (groups.Size());
    
    for (int i = 0; i < groups.Size(); i++)
      {
	MPI_Group newgroup, group_world;
	MPI_Comm_group (MPI_COMM_WORLD, &group_world);
	MPI_Group_incl (group_world, groups[i]->Size(), &(*groups[i])[0], &newgroup);
	MPI_Comm_create (MPI_COMM_WORLD, newgroup, &comm_no[i]);
      }


    for (int i = 0; i < groups.Size(); i++)
      {
	Array<int> group_members;
	for (int j = 0; j < groupnr.Size(); j++)
	  if (groupnr[j] == i)
	    group_members.Append (j);
	
	int len_vec = group_members.Size();
	if ( len_vec == 0 ) continue;

	Array<int> blocklen(len_vec);
	blocklen = 1;

	MPI_Type_indexed ( len_vec, &blocklen[0], &group_members[0], 
			   mpi_type, &mpi_t_no[i]);
	MPI_Type_commit ( &mpi_t_no[i] );
      }
    */
  }






  void ParallelDofs ::  Print() const
  {
    if ( this == 0 ) return;

    *testout << endl <<  "DEGREES OF FREEDOM FOR PARALLEL SPACES" << endl << endl;

    *testout << "nexdof = " << endl << nexdof << endl;
    *testout << "distndof = " << endl << distndof << endl;
    *testout << "sorted_exdof = " << endl << *sorted_exchangedof << endl;
    
    // (*testout) << "NDOF = "  << ndof << endl;
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

    // const MeshAccess & ma = fespace . GetMeshAccess();
    int ndof = fespace.GetNDof();
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


#endif

}

#endif
