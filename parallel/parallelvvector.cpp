#ifdef PARALLEL


#include <parallelngs.hpp>

namespace ngla
{
  using namespace ngla;
  using namespace ngparallel;



  template <typename T>
  ParallelVFlatVector<T> :: ParallelVFlatVector () 
    : paralleldofs(0)
  { 
    ;
  }


  template <typename T>
  ParallelVFlatVector<T> :: ParallelVFlatVector (int as, T * adata)
    : VFlatVector<T>(as, adata), paralleldofs(0)
  { 
    ;
  }


  template <typename T>
  ParallelVFlatVector<T> :: ParallelVFlatVector (int as, T * adata, ParallelDofs * aparalleldofs)
    : VFlatVector<T> (as, adata), paralleldofs(aparalleldofs)
  { 
    this -> SetParallelDofs ( aparalleldofs );
  }

  template <typename T>
  ParallelVFlatVector<T> :: ParallelVFlatVector (int as, T * adata,
						 ParallelDofs * aparalleldofs,
						 PARALLEL_STATUS astatus)
    : VFlatVector<T> (as, adata),
      status(astatus),
      paralleldofs(aparalleldofs)
  {
    this -> SetParallelDofs ( aparalleldofs );
  }


  template <typename T>
  ParallelVVector<T> :: ParallelVVector (int as)
    : VVector<T> (as), paralleldofs(0)
  { 
    this -> paralleldofs  = 0;
  }

  template <typename T>
  ParallelVVector<T> :: ParallelVVector (int as, ParallelDofs * aparalleldofs )
   
    : VVector<T> (as),
      paralleldofs(aparalleldofs)
  {
    if ( aparalleldofs != 0 )
      {
	this -> SetParallelDofs ( aparalleldofs );
	this->status = CUMULATED;
      }
    else
      {
	paralleldofs = 0;
	status = NOT_PARALLEL;
      }

  }

  template <typename T>
  ParallelVVector<T> :: ParallelVVector (int as, ngparallel::ParallelDofs * aparalleldofs,
					 PARALLEL_STATUS astatus ) 
    : VVector<T> (as), status(astatus)
  {
    this -> SetParallelDofs ( aparalleldofs );
  }

  template <typename T>
  ParallelVVector<T> :: ~ParallelVVector() throw()
  { 
    if ( ntasks == 1 || this -> paralleldofs == 0 ) return;
    delete this->recvvalues;
  }

  template <typename T>
  ParallelVFlatVector<T> :: ~ParallelVFlatVector() throw()
  { 
    if ( ntasks == 1 || this->paralleldofs == 0 ) return;
    delete this->recvvalues;
  }


  template <typename T>
  void ParallelVFlatVector<T> :: SetParallelDofs ( ParallelDofs * aparalleldofs, const ARRAY<int> * procs  )
  {
    this -> paralleldofs = aparalleldofs;
    if ( this -> paralleldofs == 0 ) return;


    if ( ntasks == 1 ) return;

    if ( this->recvvector_size.Size() ) return;

    this->sendvector_size . SetSize(ntasks);
    this->recvvector_size . SetSize(ntasks);
    recvvector_size = 0;

    ARRAY<int> exprocs(0);
    if ( procs == 0 )
      {
	for ( int dest = 0; dest < ntasks; dest++ )
	  if ( this->paralleldofs->IsExchangeProc(dest) )
	    exprocs.Append(dest);
      }
    else
      {
	for ( int idest = 0; idest < procs->Size(); idest++ )
	  if ( this->paralleldofs->IsExchangeProc( (*procs)[idest] ) )
	    exprocs.Append ( (*procs)[idest] );
      }
  
    for ( int idest = 0; idest < exprocs.Size(); idest++)
      {
	int dest = exprocs[idest];

	FlatArray<int>  sortedexchangedof = this->paralleldofs->GetSortedExchangeDofs(dest);
	int len_vec = sortedexchangedof.Size();
	int len_fvec = len_vec* this->entrysize;
	this->sendvector_size[dest] = len_fvec;

	MPI_Request intrequest;
	MPI_Isend( &len_fvec, 1, MPI_INT, dest, 10, MPI_COMM_WORLD, &intrequest);
      }

    for ( int idest = 0; idest < exprocs.Size(); idest++)
      {
	int dest = exprocs[idest];
	if ( ! this -> paralleldofs -> IsExchangeProc (dest) ) continue;
      
	MPI_Request intrequest;
	MPI_Status status;
	MPI_Irecv( &recvvector_size[dest], 1, MPI_INT, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &intrequest);
	MPI_Wait ( &intrequest, &status);
      }

    this->recvvalues = new Table<T>(recvvector_size); 
 }



  template <typename T>
  void ParallelVVector<T> :: SetParallelDofs ( ParallelDofs * aparalleldofs, const ARRAY<int> * procs )
  {
    this -> paralleldofs = aparalleldofs;

    if ( this -> paralleldofs == 0 ) return;
    if ( this->recvvector_size.Size() ) return;
  
#ifdef SCALASCA
#pragma pomp inst begin (vvector_setparalleldofs)
#endif
  
    sendvector_size . SetSize(ntasks);
    recvvector_size . SetSize(ntasks);
    recvvector_size = 0;
  
    ARRAY<MPI_Request> sendintrequest(ntasks);
    ARRAY<MPI_Request> recvintrequest(ntasks);

    ARRAY<int> exprocs(0);
    if ( procs == 0 )
      {
	for ( int dest = 0; dest < ntasks; dest++ )
	  if ( this->paralleldofs->IsExchangeProc(dest) )
	    exprocs.Append(dest);
      }
    else
      {
	for ( int idest = 0; idest < procs->Size(); idest++ )
	  if ( this->paralleldofs->IsExchangeProc( (*procs)[idest] ) )
	    exprocs.Append ( (*procs)[idest] );
      }

    for ( int idest = 0; idest < exprocs.Size(); idest++)
      {
	int dest = exprocs[idest];

	sendvector_size[dest] = 
          this->entrysize * this->paralleldofs->GetSortedExchangeDofs(dest).Size();

	MPI_Isend(&sendvector_size[dest], 1, MPI_INT, dest, 100, MPI_COMM_WORLD, &sendintrequest[dest]);
      }
  
    for ( int idest = 0; idest < exprocs.Size(); idest++)
      {
	int dest = exprocs[idest];
	MPI_Irecv( &recvvector_size[dest], 1, MPI_INT, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &recvintrequest[dest]);
      }

    for (int idest = 0; idest < exprocs.Size(); idest++)
      {
	MPI_Status status;
	int dest = exprocs[idest];
	MPI_Wait ( &recvintrequest[dest], &status);
	MPI_Wait ( &sendintrequest[dest], &status );
      }
 
    this -> recvvalues = new Table<T> (recvvector_size);

#ifdef SCALASCA
#pragma pomp inst end (vvector_setparalleldofs)
#endif
  }


  /// values from reduceprocs are added up,
  /// vectors in sendtoprocs are set to the cumulated values
  /// default pointer 0 means send to proc 0
  template <typename T>
  void ParallelVVector<T> :: AllReduce ( ARRAY<int> * reduceprocs, ARRAY<int> * sendtoprocs ) const 
  {
    // in case of one process only, return
    if ( status != DISTRIBUTED ) return;
    (*testout) << "ALLREDUCE" << endl;
    
#ifdef SCALASCA
#pragma pomp inst begin (vvector_allreduce)
#endif

    MPI_Status status;

    ARRAY<int> exprocs(0);
    
    // find which processors to communicate with
    for ( int i = 0; i < reduceprocs->Size(); i++)
      if ( paralleldofs->IsExchangeProc((*reduceprocs)[i]) )
	exprocs.Append((*reduceprocs)[i]);
    
    int nexprocs = exprocs.Size();
    
    ParallelVVector<T> * constvec = const_cast<ParallelVVector<T> * > (this);
    
    ARRAY<MPI_Request> sendrequest(nexprocs), recvrequest(nexprocs);

    ARRAY<int> sendto_exprocs(0);
    if ( sendtoprocs )
      {
	for ( int i = 0; i < sendtoprocs->Size(); i++ )
	  if ( paralleldofs->IsExchangeProc((*sendtoprocs)[i]) 
	       || (*sendtoprocs)[i] == id )
	    sendto_exprocs.Append((*sendtoprocs)[i] );
      }
    else
      sendto_exprocs.Append(0);

    // if the vectors are distributed, reduce
    if ( reduceprocs->Contains(id) && this->Status() == DISTRIBUTED )
      {
	// send 
	for ( int idest = 0; idest < nexprocs; idest ++ ) 
          constvec->ISend ( exprocs[idest], sendrequest[idest] );

	// receive	
	for ( int isender=0; isender < nexprocs; isender++)
          constvec -> IRecvVec ( exprocs[isender], recvrequest[isender] );

        // warten bis alles abgeschickt, danach darf man erst aufaddieren (JS)
	for ( int isender = 0;  isender < nexprocs; isender ++)
	  MPI_Wait ( &sendrequest[isender], &status);


	// cumulate
	// MPI_Waitany --> wait for first receive, not necessarily the one with smallest id
	for ( int cntexproc=0; cntexproc < nexprocs; cntexproc++ )
	  {
	    int sender, isender;
	    MPI_Waitany ( nexprocs, &recvrequest[0], &isender, &status); 
	    sender = exprocs[isender];

	    constvec->AddRecvValues(sender);
	  } 
      }

    constvec->SetStatus(CUMULATED);

    // +++++++++++++++
    // 
    // now send vector to the sendto-procs

    if ( reduceprocs->Contains(id) )
      {
	nexprocs = sendto_exprocs.Size();

	for ( int idest = 0; idest < nexprocs; idest++ )
	  {
	    int dest = sendto_exprocs[idest];
	    if ( ! paralleldofs->IsExchangeProc(dest) ) continue;
	    constvec->Send ( dest );
	  }
      }
    else if ( sendto_exprocs.Contains(id) )
      {
	for ( int isender = 0; isender < nexprocs; isender ++)
	  {
	    int sender = exprocs[isender];
	    int ii = 0; 
	    if ( ! paralleldofs->IsExchangeProc ( sender ) ) continue; 
	    constvec -> IRecvVec ( sender, recvrequest[isender] );
	    MPI_Wait( &recvrequest[isender], &status);

	    constvec -> AddRecvValues(sender);
	  }
      }
#ifdef SCALASCA
#pragma pomp inst end (vvector_horeduce)
#endif
  }
  






  template <typename T>
  void ParallelVVector<T> :: Distribute() const
  {
    if (  status != CUMULATED ) return;

    *testout << "distribute! " << endl;
    // *testout << "vector before distributing " << *this << endl;

    ParallelVVector<T> * constvec = const_cast<ParallelVVector<T> * > (this);
    constvec->SetStatus(DISTRIBUTED);

    for ( int dof = 0; dof < paralleldofs->GetNDof(); dof ++ )
      if ( ! paralleldofs->IsMasterDof ( dof ) )
	(*constvec)(dof) = 0;

    //   *testout << "distributed vector " << *constvec << endl;
  }



  /// values from reduceprocs are added up,
  /// vectors in sendtoprocs are set to the cumulated values
  /// default pointer 0 means send to proc 0
  template <typename T>
  void ParallelVFlatVector<T> :: AllReduce ( ARRAY<int> * reduceprocs, ARRAY<int> * sendtoprocs ) const 
  {
    // in case of one process only, return
    if ( status != DISTRIBUTED ) return;

#ifdef SCALASCA
#pragma pomp inst begin (vvector_allreduce)
#endif

    ARRAY<int> exprocs(0);
    int nexprocs;
    MPI_Status status;

    // find which processors to communicate with
    for ( int i = 0; i < reduceprocs->Size(); i++)
      if ( paralleldofs->IsExchangeProc((*reduceprocs)[i]) )
	exprocs.Append((*reduceprocs)[i]);
    
    nexprocs = exprocs.Size();
    int cntexproc = 0;
    
    ParallelVFlatVector<T> * constvec = const_cast<ParallelVFlatVector<T> * > (this);
    
    MPI_Request * sendrequest, *recvrequest;
    sendrequest = new MPI_Request[nexprocs];
    recvrequest = new MPI_Request[nexprocs];
    
    ARRAY<int> sendto_exprocs(0);
    if ( sendtoprocs )
      {
	for ( int i = 0; i < sendtoprocs->Size(); i++ )
	  if ( paralleldofs->IsExchangeProc((*sendtoprocs)[i]) )
	    sendto_exprocs.Append((*sendtoprocs)[i] );
      }
    else
      sendto_exprocs.Append(0);

    // if the vectors are distributed, reduce
    if ( this->status == DISTRIBUTED && reduceprocs->Contains(id) )
      {
	(*testout) << "reduce! high order" << endl;
	constvec->SetStatus(CUMULATED);
	
	// send 
	for ( int idest = 0; idest < nexprocs; idest ++ ) 
	  {
	    int dest = exprocs[idest];

	    constvec->ISend ( dest, sendrequest[idest] );
	  }

	// receive	
	for ( int isender=0; isender < nexprocs; isender++)
	  {
	    int sender = exprocs[isender];
	    constvec -> IRecvVec ( sender, recvrequest[isender] );
	  }
    

	// cumulate
    
	// MPI_Waitany --> wait for first receive, not necessarily the one with smallest id

	for ( cntexproc=0; cntexproc < nexprocs; cntexproc++ )
	  {
	    int sender, isender;
	    MPI_Waitany ( nexprocs, &recvrequest[0], &isender, &status); 
	    sender = exprocs[isender];
	    MPI_Wait ( sendrequest+isender, &status );

	    constvec->AddRecvValues(sender);
	  } 
	
      }


    // +++++++++++++++
    // 
    // now send vector to the sendto-procs


    if ( reduceprocs->Contains(id) )
      {
	nexprocs = sendto_exprocs.Size();
	delete [] sendrequest; delete [] recvrequest;
	sendrequest = new MPI_Request[nexprocs];
	recvrequest = new MPI_Request[nexprocs];
	    
	for ( int idest = 0; idest < nexprocs; idest++ )
	  {
	    int dest = sendto_exprocs[idest];
	    if ( ! paralleldofs->IsExchangeProc(dest) ) continue;
	    constvec->Send ( dest );
	  }
      }
    else if ( sendto_exprocs. Contains(id) )
      {
	for ( int isender = 0; isender < nexprocs; isender ++)
	  {
	    int sender = exprocs[isender];
	    int ii = 0; 
	    if ( ! paralleldofs->IsExchangeProc ( sender ) ) continue; 
	    constvec -> IRecvVec ( sender, recvrequest[isender] );
	    MPI_Status status;
	    MPI_Wait(recvrequest+isender, &status);
	    
	    constvec -> AddRecvValues(sender);
	  }
      }

#ifdef SCALASCA
#pragma pomp inst end (vvector_horeduce)
#endif
    delete [] recvrequest; delete []sendrequest;
  }
  




  template <typename T>
  void ParallelVFlatVector<T> :: Distribute() const
  {
    if ( status != CUMULATED ) return;
    ParallelVFlatVector * constflatvec = const_cast<ParallelVFlatVector<T> * > (this);
    T * constvec = const_cast<T * > (this->Data());
    constflatvec->SetStatus(DISTRIBUTED);
    
    for ( int dof = 0; dof < paralleldofs->GetNDof(); dof ++ )
      if ( ! paralleldofs->IsMasterDof ( dof ) )
	constvec[dof] = 0;
   
  }





  template <typename T>
  void ParallelVVector<T> :: ISend ( const int dest, MPI_Request & request )
  {
    MPI_Datatype mpi_t = this->paralleldofs->MyGetMPI_Type(dest);
    MPI_Isend( & (this->data[0]), 1, mpi_t, dest, 2000, MPI_COMM_WORLD, &request);
  }



  template <typename T>
  void ParallelVFlatVector<T> :: ISend ( const int dest, MPI_Request & request )
  {
    MPI_Datatype mpi_t = this->paralleldofs->MyGetMPI_Type(dest);
    MPI_Isend( & (this->data[0]), 1, mpi_t, dest, 2000, MPI_COMM_WORLD, &request);
  }



  template <typename T>
  void ParallelVVector<T> :: Send ( const int dest )
  {
    MPI_Datatype mpi_t = this->paralleldofs->MyGetMPI_Type(dest);
    MPI_Send( & (this->data[0]), 1, mpi_t, dest, 2001, MPI_COMM_WORLD);
  }



  template <typename T>
  void ParallelVFlatVector<T> :: Send ( const int dest )
  {
    MPI_Datatype mpi_t = this->paralleldofs->MyGetMPI_Type(dest);
    MPI_Send( & (this->data[0]), 1, mpi_t, dest, 2001, MPI_COMM_WORLD);
  }



  template <typename T>
  void ParallelVVector<T> :: IRecvVec ( const int dest, MPI_Request & request )
  {
    MPI_Datatype MPI_T = MyGetMPIType<T> ();
    MPI_Irecv( &( (*this->recvvalues)[dest][0]), recvvector_size[dest], MPI_T, dest, 
	       MPI_ANY_TAG, MPI_COMM_WORLD, &request);
  }

  template <typename T>
  void ParallelVFlatVector<T> :: IRecvVec ( const int dest, MPI_Request & request )
  {
    MPI_Datatype MPI_T = MyGetMPIType<T> ();
    MPI_Irecv(&( (*this->recvvalues)[dest][0]), recvvector_size[dest], MPI_T, 
	      dest, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
  }

template <typename T>
void ParallelVFlatVector<T> :: PrintParallelDofs () const
{ 
  paralleldofs->Print(); 
}

template <typename T>
void ParallelVVector<T> :: PrintParallelDofs () const
{ 
  paralleldofs->Print(); 
}


template <typename T>
void ParallelVFlatVector<T> :: SetStatus ( PARALLEL_STATUS astatus )
{ 
  status = astatus; 
}

template <typename T>
void ParallelVVector<T> :: SetStatus ( PARALLEL_STATUS astatus ) 
{ status = astatus; }



template <typename T>
BaseVector * ParallelVVector<T> :: CreateVector ( const ARRAY<int> * procs ) const
{
  ParallelVVector<T> * parvec;
  
  parvec = new ParallelVVector<T> (this->size);
  parvec->SetStatus( this->Status() );
  
  // *testout << "procs... " << procs << endl;
  if( this->Status() != NOT_PARALLEL )
    {
      ParallelDofs * aparalleldofs = this-> GetParallelDofs();
      parvec->SetParallelDofs ( aparalleldofs, procs );
    }
  return parvec;
  
}



  template <typename T>
  BaseVector * ParallelVFlatVector<T> :: CreateVector ( const ARRAY<int> * procs ) const
  {
    ParallelVVector<T> * parvec;

    parvec = new ParallelVVector<T> (this->size);
    parvec->SetStatus( this->Status() );
    
    if( this->Status() != NOT_PARALLEL )
      {
	ParallelDofs * aparalleldofs = this-> GetParallelDofs();
	parvec->SetParallelDofs ( aparalleldofs, procs );
      }
    return parvec;
  }


  template < class T >
  ostream & ParallelVVector<T> :: Print (ostream & ost) const
  {
    ost << "addr = " << &this->FV()(0) << endl;
    const PARALLEL_STATUS status = this->Status();
    if ( status == NOT_PARALLEL )
      ost << "NOT PARALLEL" << endl;
    else if ( status == DISTRIBUTED )
      ost << "DISTRIBUTED" <<endl;
    else if ( status == CUMULATED )
      ost << "CUMULATED" << endl;
    return (ost << this->FV() << endl);
  }
  
  template < class T >
  ostream & ParallelVFlatVector<T> :: Print (ostream & ost) const
  {
    ost << "addr = " << &this->FV()(0) << endl;
    const PARALLEL_STATUS status = this->Status();
    if ( status == NOT_PARALLEL )
      ost << "NOT PARALLEL" << endl;
    else if ( status == DISTRIBUTED )
      ost << "DISTRIBUTED" <<endl;
    else if ( status == CUMULATED )
      ost << "CUMULATED" << endl;
    return (ost << this->FV() << endl);
  }

  template <class T>
  void ParallelVVector<T> :: AddRecvValues( int sender )
  {
    FlatArray<int> exdofs = paralleldofs->GetSortedExchangeDofs(sender);
    for (int i = 0; i < exdofs.Size(); i++)
      (*this) (exdofs[i]) +=  (*this->recvvalues)[sender][i];
  }


  template <class T>
  void ParallelVFlatVector<T> :: AddRecvValues( int sender )
  {
    FlatArray<int> exdofs = paralleldofs->GetSortedExchangeDofs(sender);
    for (int i = 0; i < exdofs.Size(); i++)
      this->Data()[exdofs[i]] += (*this->recvvalues)[sender][i];
  }


template class ParallelVFlatVector<double>;
template class ParallelVFlatVector<Complex>;
template class ParallelVFlatVector<Vec<1,double> >;
template class ParallelVFlatVector<Vec<1,Complex> >;
template class ParallelVFlatVector<Vec<2,double> >;
template class ParallelVFlatVector<Vec<2,Complex> >;
template class ParallelVFlatVector<Vec<3,double> >;
template class ParallelVFlatVector<Vec<3,Complex> >;
template class ParallelVFlatVector<Vec<4,double> >;
template class ParallelVFlatVector<Vec<4,Complex> >;

template class ParallelVFlatVector<Vec<5,double> >;
template class ParallelVFlatVector<Vec<5,Complex> >;
template class ParallelVFlatVector<Vec<6,double> >;
template class ParallelVFlatVector<Vec<6,Complex> >;
template class ParallelVFlatVector<Vec<7,double> >;
template class ParallelVFlatVector<Vec<7,Complex> >;
template class ParallelVFlatVector<Vec<8,double> >;
template class ParallelVFlatVector<Vec<8,Complex> >;
template class ParallelVFlatVector<Vec<9,double> >;
template class ParallelVFlatVector<Vec<9,Complex> >;

template class ParallelVFlatVector<Vec<12,double> >;
template class ParallelVFlatVector<Vec<12,Complex> >;
template class ParallelVFlatVector<Vec<18,double> >;
template class ParallelVFlatVector<Vec<18,Complex> >;
template class ParallelVFlatVector<Vec<24,double> >;
template class ParallelVFlatVector<Vec<24,Complex> >;


template class ParallelVVector<double>;
template class ParallelVVector<Complex>;
template class ParallelVVector<Vec<1,double> >;
template class ParallelVVector<Vec<1,Complex> >;
template class ParallelVVector<Vec<2,double> >;
template class ParallelVVector<Vec<2,Complex> >;
template class ParallelVVector<Vec<3,double> >;
template class ParallelVVector<Vec<3,Complex> >;
template class ParallelVVector<Vec<4,double> >;
template class ParallelVVector<Vec<4,Complex> >;

template class ParallelVVector<Vec<5,double> >;
template class ParallelVVector<Vec<5,Complex> >;
template class ParallelVVector<Vec<6,double> >;
template class ParallelVVector<Vec<6,Complex> >;
template class ParallelVVector<Vec<7,double> >;
template class ParallelVVector<Vec<7,Complex> >;
template class ParallelVVector<Vec<8,double> >;
template class ParallelVVector<Vec<8,Complex> >;
template class ParallelVVector<Vec<9,double> >;
template class ParallelVVector<Vec<9,Complex> >;
template class ParallelVVector<Vec<10,double> >;
template class ParallelVVector<Vec<10,Complex> >;
template class ParallelVVector<Vec<11,double> >;
template class ParallelVVector<Vec<11,Complex> >;
template class ParallelVVector<Vec<15,double> >;
template class ParallelVVector<Vec<15,Complex> >;
template class ParallelVVector<Vec<13,double> >;
template class ParallelVVector<Vec<13,Complex> >;
template class ParallelVVector<Vec<14,double> >;
template class ParallelVVector<Vec<14,Complex> >;

template class ParallelVVector<Vec<12,double> >;
template class ParallelVVector<Vec<12,Complex> >;
template class ParallelVVector<Vec<18,double> >;
template class ParallelVVector<Vec<18,Complex> >;
template class ParallelVVector<Vec<24,double> >;
template class ParallelVVector<Vec<24,Complex> >;

  
  //  template class ParallelVFlatVector<FlatVector<double> >;
  // template class ParallelVVector<FlatVector<double> >;

}
#endif

