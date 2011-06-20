#ifdef PARALLEL


#include <parallelngs.hpp>
#include <la.hpp>


namespace ngla
{
  using namespace ngparallel;


  BaseVector & ParallelBaseVector :: SetScalar (double scal)
  {
    FVDouble() = scal;
    if ( IsParallelVector() )
      this->SetStatus(CUMULATED);
    else
      this->SetStatus(NOT_PARALLEL);
    return *this;
  }

  BaseVector & ParallelBaseVector :: SetScalar (Complex scal)
  {
    FVComplex() = scal;
    if ( IsParallelVector() )
      this->SetStatus(CUMULATED);
    else
      this->SetStatus(NOT_PARALLEL);
    return *this;
  }

  BaseVector & ParallelBaseVector :: Set (double scal, const BaseVector & v)
  {
    FVDouble() = scal * v.FVDouble();
    const ParallelBaseVector * parv = dynamic_cast<const ParallelBaseVector *> (&v);

    if ( parv && parv->IsParallelVector() )
      {
	this->SetParallelDofs (parv->GetParallelDofs());
	this->SetStatus(parv->Status());
      }
    else 
      {
	this->SetParallelDofs(0);
	this->SetStatus(NOT_PARALLEL);
      }
    return *this;
  }

  BaseVector & ParallelBaseVector :: Set (Complex scal, const BaseVector & v)
  {
    FVComplex() = scal * v.FVComplex();
    const ParallelBaseVector * parv = dynamic_cast<const ParallelBaseVector *> (&v);

    if ( parv->IsParallelVector() )
      this->SetParallelDofs(parv->GetParallelDofs());
    else
      this->SetParallelDofs(0);

    this->SetStatus(parv->Status());
    return *this;
  }
    
  BaseVector & ParallelBaseVector :: Add (double scal, const BaseVector & v)
  {
    const ParallelBaseVector * parv = dynamic_cast<const ParallelBaseVector *> (&v);

    if ( (*this).Status() != parv->Status() )
      {
        if ( (*this).Status() == DISTRIBUTED )
	  Cumulate();
          // AllReduce(&hoprocs);
        else 
	  parv -> Cumulate();
          // parv->AllReduce(&hoprocs);
      }
    FVDouble() += scal * parv->FVDouble();
    return *this;
  }

  BaseVector & ParallelBaseVector :: Add (Complex scal, const BaseVector & v)
  {
    const ParallelBaseVector * parv = dynamic_cast<const ParallelBaseVector *> (&v);

    if ( (*this).Status() != parv->Status() )
      {
        if ( (*this).Status() == DISTRIBUTED ) 
	  Cumulate();
        else 
	  parv->Cumulate();
      }

    FVComplex() += scal * parv->FVComplex();
    return *this;
  }




  void ParallelBaseVector :: Cumulate () const
  {
    if ( status != DISTRIBUTED ) return;
    (*testout) << "CUMULATE" << endl;

    
    MPI_Status status;
    Array<int> exprocs;
    
    // find which processors to communicate with
    for ( int i = 1; i < ntasks; i++)
      // if ( paralleldofs->IsExchangeProc(i) )
      if (paralleldofs -> GetSortedExchangeDofs (i).Size())
	exprocs.Append(i);
    
    int nexprocs = exprocs.Size();
    
    ParallelBaseVector * constvec = const_cast<ParallelBaseVector * > (this);
    
    Array<MPI_Request> sendrequest(nexprocs), recvrequest(nexprocs);

    // if the vectors are distributed, reduce
    if ( id >= 1)
      {
	for ( int idest = 0; idest < nexprocs; idest ++ ) 
          constvec->ISend ( exprocs[idest], sendrequest[idest] );

	for ( int isender=0; isender < nexprocs; isender++)
          constvec -> IRecvVec ( exprocs[isender], recvrequest[isender] );

        // wait till everything is sent 
	for ( int isender = 0;  isender < nexprocs; isender ++)
	  MPI_Wait ( &sendrequest[isender], &status);

	// cumulate
	// MPI_Waitany --> wait for first receive, not necessarily the one with smallest id
	for ( int cntexproc=0; cntexproc < nexprocs; cntexproc++ )
	  {
	    int isender;
	    MPI_Waitany ( nexprocs, &recvrequest[0], &isender, &status); 
	    constvec->AddRecvValues(exprocs[isender]);
	  } 
      }

    SetStatus(CUMULATED);
  }









  template <class SCAL>
  SCAL S_ParallelBaseVector<SCAL> :: InnerProduct (const BaseVector & v2) const
  {
    const ParallelBaseVector * parv2 = dynamic_cast<const ParallelBaseVector *> (&v2);

    SCAL globalsum = 0;
    if ( id == 0 && ntasks > 1 )
      {
	if (this->Status() == parv2->Status() && this->Status() == DISTRIBUTED )
	  {
	    Cumulate();
	    // this->AllReduce(&hoprocs);
	  }
	// two cumulated vectors -- distribute one
	else if ( this->Status() == parv2->Status() && this->Status() == CUMULATED )
	  this->Distribute();
	MyMPI_Recv ( globalsum, 1 );
      }
    else
      {
	// not parallel
	if ( this->Status() == NOT_PARALLEL && parv2->Status() == NOT_PARALLEL )
	  return ngbla::InnerProduct (this->FVScal(), 
				      dynamic_cast<const S_BaseVector<SCAL>&>(*parv2).FVScal());
	// two distributed vectors -- cumulate one
	else if ( this->Status() == parv2->Status() && this->Status() == DISTRIBUTED )
	  {
	    Cumulate();
	    // this->AllReduce(&hoprocs);
	  }
	// two cumulated vectors -- distribute one
	else if ( this->Status() == parv2->Status() && this->Status() == CUMULATED )
	  this->Distribute();

	SCAL localsum = ngbla::InnerProduct (this->FVScal(), 
					     dynamic_cast<const S_BaseVector<SCAL>&>(*parv2).FVScal());
	MPI_Datatype MPI_SCAL = MyGetMPIType<SCAL>();
      
	MPI_Allreduce ( &localsum, &globalsum, 1,  MPI_SCAL, MPI_SUM, MPI_HIGHORDER_COMM); //MPI_COMM_WORLD);
	if ( id == 1 )
	  MyMPI_Send( globalsum, 0 );
      }
  
    return globalsum;
  }




  template <>
  Complex S_ParallelBaseVector<Complex> :: InnerProduct (const BaseVector & v2) const
  {
    const ParallelBaseVector * parv2 = dynamic_cast<const ParallelBaseVector *> (&v2);
  
    // not parallel
    if ( this->Status() == NOT_PARALLEL && parv2->Status() == NOT_PARALLEL )
      return ngbla::InnerProduct (FVScal(), 
				  dynamic_cast<const S_BaseVector<Complex>&>(*parv2).FVScal());
    // two distributed vectors -- cumulate one
    else if (this->Status() == parv2->Status() && this->Status() == DISTRIBUTED )
      {
	Cumulate();
	// this->AllReduce(&hoprocs);
      }
    // two cumulated vectors -- distribute one
    else if ( this->Status() == parv2->Status() && this->Status() == CUMULATED )
      Distribute();

    Complex localsum ;
    Complex globalsum = 0;
    if ( id == 0 )
      localsum = 0;
    else 
      localsum = ngbla::InnerProduct (FVComplex(), 
				      dynamic_cast<const S_BaseVector<Complex>&>(*parv2).FVComplex());
    MPI_Allreduce ( &localsum, &globalsum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce ( &localsum, &globalsum, 2, MPI_DOUBLE, MPI_SUM, MPI_HIGHORDER_WORLD);
 
    return globalsum;
  }




  // template
  // Complex S_BaseVector<Complex> :: InnerProduct (const BaseVector & v2) const;

  //template <>
  // Complex S_BaseVector<Complex> :: InnerProduct (const BaseVector & v2) const


  template class S_ParallelBaseVector<double>;
  template class S_ParallelBaseVector<Complex>;










  template <typename T>
  ParallelVFlatVector<T> :: ParallelVFlatVector () 
  { 
    ;
  }


  template <typename T>
  ParallelVFlatVector<T> :: ParallelVFlatVector (int as, T * adata)
    : VFlatVector<T>(as, adata)
  { 
    ;
  }


  template <typename T>
  ParallelVFlatVector<T> :: ParallelVFlatVector (int as, T * adata, ParallelDofs * aparalleldofs)
    : VFlatVector<T> (as, adata)
  { 
    this -> SetParallelDofs ( aparalleldofs );
  }

  template <typename T>
  ParallelVFlatVector<T> :: 
  ParallelVFlatVector (int as, T * adata,
		       ParallelDofs * aparalleldofs,
		       PARALLEL_STATUS astatus)
    : VFlatVector<T> (as, adata)
  {
    status = astatus;
    this -> SetParallelDofs ( aparalleldofs );
  }


  template <typename T>
  ParallelVVector<T> :: ParallelVVector (int as)
    : VVector<T> (as)
  { 
    this -> paralleldofs  = 0;
  }

  template <typename T>
  ParallelVVector<T> :: ParallelVVector (int as, ParallelDofs * aparalleldofs )
   
    : VVector<T> (as)
					// paralleldofs(aparalleldofs)
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
    : VVector<T> (as)
  {
    status = astatus;
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
  void ParallelVFlatVector<T> :: SetParallelDofs ( ParallelDofs * aparalleldofs, const Array<int> * procs  )
  {
    this -> paralleldofs = aparalleldofs;
    if ( this -> paralleldofs == 0 ) return;

    Array<int> exdofs(ntasks);
    for (int i = 0; i < ntasks; i++)
      exdofs[i] = this->paralleldofs->GetSortedExchangeDofs(i).Size();
    this -> recvvalues = new Table<T> (exdofs);
  }


  template <typename T>
  void ParallelVVector<T> :: SetParallelDofs ( ParallelDofs * aparalleldofs, const Array<int> * procs )
  {
    this -> paralleldofs = aparalleldofs;
    if ( this -> paralleldofs == 0 ) return;
    
    Array<int> exdofs(ntasks);
    for (int i = 0; i < ntasks; i++)
      exdofs[i] = this->paralleldofs->GetSortedExchangeDofs(i).Size();
    this -> recvvalues = new Table<T> (exdofs);
  }





  template <typename T>
  void ParallelVVector<T> :: Distribute() const
  {
    if (status != CUMULATED) return;

    // *testout << "distribute! " << endl;
    // *testout << "vector before distributing " << *this << endl;

    this -> SetStatus(DISTRIBUTED);

    for ( int dof = 0; dof < paralleldofs->GetNDof(); dof ++ )
      if ( ! paralleldofs->IsMasterDof ( dof ) )
	(*this)(dof) = 0;

    //   *testout << "distributed vector " << *constvec << endl;
  }



  template <typename T>
  void ParallelVFlatVector<T> :: Distribute() const
  {
    if ( status != CUMULATED ) return;
    ParallelVFlatVector * constflatvec = const_cast<ParallelVFlatVector<T> * > (this);
    // T * constvec = const_cast<T * > (this->pdata);
    T * constvec = this->FV().Addr(0);
    constflatvec->SetStatus(DISTRIBUTED);
    
    for ( int dof = 0; dof < paralleldofs->GetNDof(); dof ++ )
      if ( ! paralleldofs->IsMasterDof ( dof ) )
	constvec[dof] = 0;
  }




  void ParallelBaseVector :: ISend ( int dest, MPI_Request & request ) const
  {
    MPI_Datatype mpi_t = this->paralleldofs->MyGetMPI_Type(dest);
    MPI_Isend( Memory(), 1, mpi_t, dest, 2000, MPI_COMM_WORLD, &request);
  }

  void ParallelBaseVector :: Send ( int dest ) const
  {
    MPI_Datatype mpi_t = this->paralleldofs->MyGetMPI_Type(dest);
    MPI_Send( Memory(), 1, mpi_t, dest, 2001, MPI_COMM_WORLD);
  }

  template <typename T>
  void ParallelVVector<T> :: IRecvVec ( int dest, MPI_Request & request )
  {
    MPI_Datatype MPI_T = MyGetMPIType<T> ();
    MPI_Irecv( &( (*this->recvvalues)[dest][0]), 
	       (*this->recvvalues)[dest].Size(), 
	       MPI_T, dest, 
	       MPI_ANY_TAG, MPI_COMM_WORLD, &request);
  }

  template <typename T>
  void ParallelVFlatVector<T> :: IRecvVec ( int dest, MPI_Request & request )
  {
    MPI_Datatype MPI_T = MyGetMPIType<T> ();
    MPI_Irecv(&( (*this->recvvalues)[dest][0]), 
	      (*this->recvvalues)[dest].Size(), 
	      MPI_T, 
	      dest, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
  }

  template <typename T>
  void ParallelVVector<T> :: RecvVec ( int dest)
  {
    MPI_Status status;

    MPI_Datatype MPI_T = MyGetMPIType<T> ();
    MPI_Recv( &( (*this->recvvalues)[dest][0]), 
	      (*this->recvvalues)[dest].Size(), 
	      MPI_T, dest, 
	      MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

  template <typename T>
  void ParallelVFlatVector<T> :: RecvVec ( int dest) 
  {
    MPI_Status status;

    MPI_Datatype MPI_T = MyGetMPIType<T> ();
    MPI_Recv(&( (*this->recvvalues)[dest][0]), 
	     (*this->recvvalues)[dest].Size(), 
	     MPI_T, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }





  template <typename T>
  void ParallelVFlatVector<T> :: PrintParallelDofs () const
  { 
    ; // paralleldofs->Print(); 
  }

  template <typename T>
  void ParallelVVector<T> :: PrintParallelDofs () const
  { 
    ; // paralleldofs->Print(); 
  }




  template <typename T>
  BaseVector * ParallelVVector<T> :: CreateVector ( const Array<int> * procs ) const
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
  double ParallelVVector<T> :: L2Norm () const
  {
    this->Cumulate();
    // this->AllReduce (&hoprocs);
    
    double sum = 0;

    FlatVector<double> fv = this -> FVDouble ();

    if (id > 0)
      for (int dof = 0; dof < paralleldofs->GetNDof(); dof++)
	if (paralleldofs->IsMasterDof ( dof ) )
	  sum += L2Norm2 (fv[dof]);
    
    double globalsum = 0;

    MPI_Datatype MPI_SCAL = MyGetMPIType<double>();
    MPI_Allreduce (&sum, &globalsum, 1, MPI_SCAL, MPI_SUM, MPI_COMM_WORLD);

    return sqrt (globalsum);
  }


  template <typename T>
  BaseVector * ParallelVFlatVector<T> :: CreateVector ( const Array<int> * procs ) const
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
    // if (sender != 0)
    if (1)
      {
	for (int i = 0; i < exdofs.Size(); i++)
	  (*this) (exdofs[i]) +=  (*this->recvvalues)[sender][i];
      }
    else
      {
	for (int i = 0; i < exdofs.Size(); i++)
	  (*this)(i) += (*this->recvvalues)[0][i];
      }
  }


  template <class T>
  void ParallelVFlatVector<T> :: AddRecvValues( int sender )
  {
    FlatArray<int> exdofs = paralleldofs->GetSortedExchangeDofs(sender);
    // if (sender != 0)
    if (1)
      {
	for (int i = 0; i < exdofs.Size(); i++)
	  this->FV()(exdofs[i]) += (*this->recvvalues)[sender][i];
      }
    else
      {
	for (int i = 0; i < exdofs.Size(); i++)
	  this->FV()(i) += (*this->recvvalues)[0][i];
	// (*this).Range (0, (*this->recvvalues).Size()) += (*this->recvvalues)[sender];
      }
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
  /*
    template class ParallelVFlatVector<Vec<9,double> >;
    template class ParallelVFlatVector<Vec<9,Complex> >;

    template class ParallelVFlatVector<Vec<12,double> >;
    template class ParallelVFlatVector<Vec<12,Complex> >;
    template class ParallelVFlatVector<Vec<18,double> >;
    template class ParallelVFlatVector<Vec<18,Complex> >;
    template class ParallelVFlatVector<Vec<24,double> >;
    template class ParallelVFlatVector<Vec<24,Complex> >;
  */

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

  /*
    template class ParallelVVector<Vec<9,double> >;
    template class ParallelVVector<Vec<9,Complex> >;
    template class ParallelVVector<Vec<10,double> >;
    template class ParallelVVector<Vec<10,Complex> >;
    // template class ParallelVVector<Vec<11,double> >;
    // template class ParallelVVector<Vec<11,Complex> >;
    // template class ParallelVVector<Vec<15,double> >;
    // template class ParallelVVector<Vec<15,Complex> >;
    // template class ParallelVVector<Vec<13,double> >;
    // template class ParallelVVector<Vec<13,Complex> >;
    // template class ParallelVVector<Vec<14,double> >;
    // template class ParallelVVector<Vec<14,Complex> >;

    template class ParallelVVector<Vec<12,double> >;
    template class ParallelVVector<Vec<12,Complex> >;
    template class ParallelVVector<Vec<18,double> >;
    template class ParallelVVector<Vec<18,Complex> >;
    template class ParallelVVector<Vec<24,double> >;
    template class ParallelVVector<Vec<24,Complex> >;
  */
  
  //  template class ParallelVFlatVector<FlatVector<double> >;
  // template class ParallelVVector<FlatVector<double> >;

}
#endif

