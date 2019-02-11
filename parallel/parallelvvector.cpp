// #ifdef PARALLEL


#include <parallelngs.hpp>
#include <la.hpp>


namespace ngla
{

  BaseVector & ParallelBaseVector :: SetScalar (double scal)
  {
    if (IsComplex())
      FVComplex() = scal;
    else
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
    const ParallelBaseVector * parv = dynamic_cast_ParallelBaseVector (&v);

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
    const ParallelBaseVector * parv = dynamic_cast_ParallelBaseVector (&v);

    if ( parv->IsParallelVector() )
      this->SetParallelDofs(parv->GetParallelDofs());
    else
      this->SetParallelDofs(0);

    this->SetStatus(parv->Status());
    return *this;
  }
    
  BaseVector & ParallelBaseVector :: Add (double scal, const BaseVector & v)
  {
    const ParallelBaseVector * parv = dynamic_cast_ParallelBaseVector (&v);

    if ( (*this).Status() != parv->Status() )
      {
        if ( (*this).Status() == DISTRIBUTED )
	  Cumulate();
        else 
	  parv -> Cumulate();
      }
    FVDouble() += scal * parv->FVDouble();
    return *this;
  }

  BaseVector & ParallelBaseVector :: Add (Complex scal, const BaseVector & v)
  {
    const ParallelBaseVector * parv = dynamic_cast_ParallelBaseVector (&v);

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


  void ParallelBaseVector :: PrintStatus ( ostream & ost ) const
  {
    if ( this->status == NOT_PARALLEL )
      ost << "NOT PARALLEL" << endl;
    else if ( this->status == DISTRIBUTED )
      ost << "DISTRIBUTED" << endl ;
    else if (this->status == CUMULATED )
	ost << "CUMULATED" << endl ;
  }
  

  void ParallelBaseVector :: Cumulate () const
  {
#ifdef PARALLEL
    if (status != DISTRIBUTED) return;
    
    int ntasks = paralleldofs->GetNTasks();
    Array<int> exprocs;
    for (int i = 0; i < ntasks; i++)
      if (paralleldofs -> GetExchangeDofs (i).Size())
	exprocs.Append(i);
    
    int nexprocs = exprocs.Size();
    
    ParallelBaseVector * constvec = const_cast<ParallelBaseVector * > (this);
    
    Array<MPI_Request> sendrequest(nexprocs), recvrequest(nexprocs);

    for (int idest = 0; idest < nexprocs; idest ++ ) 
      constvec->ISend (exprocs[idest], sendrequest[idest] );
    for (int isender=0; isender < nexprocs; isender++)
      constvec -> IRecvVec (exprocs[isender], recvrequest[isender] );
    
    MyMPI_WaitAll (sendrequest);
    
    // cumulate
    for (int cntexproc=0; cntexproc < nexprocs; cntexproc++)
      {
	int isender = MyMPI_WaitAny (recvrequest);
	constvec->AddRecvValues(exprocs[isender]);
      } 

    SetStatus(CUMULATED);
#endif
  }



  void ParallelBaseVector :: ISend ( int dest, MPI_Request & request ) const
  {
#ifdef PARALLEL
    MPI_Datatype mpi_t = this->paralleldofs->MyGetMPI_Type(dest);
    MPI_Isend( Memory(), 1, mpi_t, dest, MPI_TAG_SOLVE, this->paralleldofs->GetCommunicator(), &request);
#endif
  }

  /*
  void ParallelBaseVector :: Send ( int dest ) const
  {
    MPI_Datatype mpi_t = this->paralleldofs->MyGetMPI_Type(dest);
    MPI_Send( Memory(), 1, mpi_t, dest, MPI_TAG_SOLVE, ngs_comm);
  }
  */





  template <class SCAL>
  SCAL S_ParallelBaseVector<SCAL> :: InnerProduct (const BaseVector & v2) const
  {
    const ParallelBaseVector * parv2 = dynamic_cast_ParallelBaseVector(&v2);

    if ( this->Status() == NOT_PARALLEL && parv2->Status() == NOT_PARALLEL )
      return ngbla::InnerProduct (this->FVScal(), 
				  dynamic_cast<const S_BaseVector<SCAL>&>(*parv2).FVScal());

    // two distributed vectors -- cumulate one
    else if ( this->Status() == parv2->Status() && this->Status() == DISTRIBUTED )
      Cumulate();
    
    // two cumulated vectors -- distribute one
    else if ( this->Status() == parv2->Status() && this->Status() == CUMULATED )
      this->Distribute();
    
    SCAL localsum = ngbla::InnerProduct (this->FVScal(), 
					 dynamic_cast<const S_BaseVector<SCAL>&>(*parv2).FVScal());


    return MyMPI_AllReduce (localsum, MPI_SUM, paralleldofs->GetCommunicator());
  }




  template <>
  Complex S_ParallelBaseVector<Complex> :: InnerProduct (const BaseVector & v2) const
  {
    const ParallelBaseVector * parv2 = dynamic_cast_ParallelBaseVector(&v2);
  
    // not parallel
    if ( this->Status() == NOT_PARALLEL && parv2->Status() == NOT_PARALLEL )
      return ngbla::InnerProduct (FVScal(), 
				  dynamic_cast<const S_BaseVector<Complex>&>(*parv2).FVScal());
    // two distributed vectors -- cumulate one
    else if (this->Status() == parv2->Status() && this->Status() == DISTRIBUTED )
      Cumulate();

    // two cumulated vectors -- distribute one
    else if ( this->Status() == parv2->Status() && this->Status() == CUMULATED )
      Distribute();

    Complex localsum = ngbla::InnerProduct (FVComplex(), 
                                            dynamic_cast<const S_BaseVector<Complex>&>(*parv2).FVComplex());

    /*    
    Complex globalsum = 0;
#ifdef PARALLEL
    MPI_Allreduce (&localsum, &globalsum, 2, MPI_DOUBLE, MPI_SUM, ngs_comm);
#else
    globalsum = localsum;
#endif
    return globalsum;
    */
    return MyMPI_AllReduce (localsum, MPI_SUM, paralleldofs->GetCommunicator());
  }

  template class S_ParallelBaseVector<double>;
  template class S_ParallelBaseVector<Complex>;



  /*
  template <class SCAL>
  S_ParallelBaseVectorPtr<SCAL> :: S_ParallelBaseVectorPtr (int as, int aes, void * adata) throw()
    : S_BaseVectorPtr<SCAL> (as, aes, adata)
  { 
    cout << "gibt's denn das ?" << endl;
    recvvalues = NULL;
    this -> paralleldofs = 0;
    status = NOT_PARALLEL;
  }
  */

  template <class SCAL>
  S_ParallelBaseVectorPtr<SCAL> :: 
  S_ParallelBaseVectorPtr (int as, int aes, 
			   shared_ptr<ParallelDofs> apd, PARALLEL_STATUS stat) throw()
    : S_BaseVectorPtr<SCAL> (as, aes)
  { 
    recvvalues = NULL;
    if ( apd != 0 )
      {
	this -> SetParallelDofs ( apd );
	status = stat;
      }
    else
      {
	paralleldofs = 0;
	status = NOT_PARALLEL;
      }
    local_vec = make_shared<S_BaseVectorPtr<SCAL>>(as, aes, (void*)pdata);
  }

  template <class SCAL>
  S_ParallelBaseVectorPtr<SCAL> :: ~S_ParallelBaseVectorPtr ()
  {
    delete recvvalues;
  }


  template <typename SCAL>
  void S_ParallelBaseVectorPtr<SCAL> :: 
  SetParallelDofs (shared_ptr<ParallelDofs> aparalleldofs, const Array<int> * procs )
  {
    if (this->paralleldofs == aparalleldofs) return;

    this -> paralleldofs = aparalleldofs;
    if ( this -> paralleldofs == 0 ) return;
    
    int ntasks = this->paralleldofs->GetNTasks();
    Array<int> exdofs(ntasks);
    for (int i = 0; i < ntasks; i++)
      exdofs[i] = this->es * this->paralleldofs->GetExchangeDofs(i).Size();
    delete this->recvvalues;
    this -> recvvalues = new Table<TSCAL> (exdofs);
  }


  template <typename SCAL>
  void S_ParallelBaseVectorPtr<SCAL>  :: Distribute() const
  {
    if (status != CUMULATED) return;
    this->SetStatus(DISTRIBUTED);

    for ( int dof = 0; dof < paralleldofs->GetNDofLocal(); dof ++ )
      if ( ! paralleldofs->IsMasterDof ( dof ) )
	(*this)(dof) = 0;
  }



  template < class SCAL >
  ostream & S_ParallelBaseVectorPtr<SCAL> :: Print (ostream & ost) const
  {
    this->PrintStatus (ost);
    S_BaseVectorPtr<SCAL>::Print (ost);
    return ost;
  }


  template <typename SCAL>
  void S_ParallelBaseVectorPtr<SCAL> :: IRecvVec ( int dest, MPI_Request & request )
  {
#ifdef PARALLEL
    MPI_Datatype MPI_TS = MyGetMPIType<TSCAL> ();
    MPI_Irecv( &( (*recvvalues)[dest][0]), 
	       (*recvvalues)[dest].Size(), 
	       MPI_TS, dest, 
	       MPI_TAG_SOLVE, this->paralleldofs->GetCommunicator(), &request);
#endif
  }

  /*
  template <typename SCAL>
  void S_ParallelBaseVectorPtr<SCAL> :: RecvVec ( int dest)
  {
    MPI_Status status;

    MPI_Datatype MPI_TS = MyGetMPIType<TSCAL> ();
    MPI_Recv( &( (*recvvalues)[dest][0]), 
	      (*recvvalues)[dest].Size(), 
	      MPI_TS, dest, 
	      MPI_TAG_SOLVE, ngs_comm, &status);
  }
  */


  template <typename SCAL>
  void S_ParallelBaseVectorPtr<SCAL> :: AddRecvValues( int sender )
  {
    FlatArray<int> exdofs = paralleldofs->GetExchangeDofs(sender);
    FlatMatrix<SCAL> rec (exdofs.Size(), this->es, &(*this->recvvalues)[sender][0]);
    for (int i = 0; i < exdofs.Size(); i++)
      (*this) (exdofs[i]) += rec.Row(i);
  }



  template <typename SCAL>
  AutoVector S_ParallelBaseVectorPtr<SCAL> :: 
  CreateVector () const
  {
    return make_shared<S_ParallelBaseVectorPtr<TSCAL>>
      (this->size, this->es, paralleldofs, status);
    /*
    S_ParallelBaseVectorPtr<TSCAL> * parvec = 
      new S_ParallelBaseVectorPtr<TSCAL> (this->size, this->es, paralleldofs, status);
    return parvec;
    */
  }



  template <typename SCAL>
  double S_ParallelBaseVectorPtr<SCAL> :: L2Norm () const
  {
    this->Cumulate();
    
    double sum = 0;
    
    if (this->entrysize == 1)
      {
	FlatVector<double> fv = this -> FVDouble ();
	for (int dof = 0; dof < paralleldofs->GetNDofLocal(); dof++)
	  if (paralleldofs->IsMasterDof ( dof ) )
	    sum += sqr (fv[dof]);
      }
    else
      {
	FlatMatrix<double> fv (paralleldofs->GetNDofLocal(), 
			       this->entrysize, (double*)this->Memory());
	for (int dof = 0; dof < paralleldofs->GetNDofLocal(); dof++)
	  if (paralleldofs->IsMasterDof ( dof ) )
	    sum += L2Norm2 (fv.Row(dof));
      }
      
    double globsum = MyMPI_AllReduce (sum, MPI_SUM, ngs_comm);
    return sqrt (globsum);
  }


  template class S_ParallelBaseVectorPtr<double>;
  template class S_ParallelBaseVectorPtr<Complex>;
}



// #endif

