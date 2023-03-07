// #ifdef PARALLEL


#include <parallelngs.hpp>
#include <la.hpp>


namespace ngla
{

  class ParallelRangeVector : public ParallelBaseVector
  {
  protected:
    const ParallelBaseVector * orig;
    IntRange range;
  public:
    ParallelRangeVector (const ParallelBaseVector * aorig, IntRange arange)
      : BaseVector(), ParallelBaseVector(), orig(aorig), range(arange)
    {
      BaseVector::entrysize = aorig->EntrySize();
      this->size = range.Size();
      status = orig->GetParallelStatus();
      local_vec = orig->GetLocalVector()->Range(range);
    }
    
    ~ParallelRangeVector () override { ; }

    void * Memory () const override
    { return orig->GetLocalVector()->Range(range).Memory(); }

    virtual int EntrySizeScal() const throw () override
    { return orig->EntrySizeScal(); }

    virtual bool IsComplex() const override
    { return orig->IsComplex(); } 
    
    FlatVector<double> FVDouble () const override
    { return orig->GetLocalVector()->Range(range).FVDouble(); }

    FlatVector<Complex> FVComplex () const override
    { return orig->GetLocalVector()->Range(range).FVComplex(); }

    AutoVector Range (T_Range<size_t> range2) const override
    {
      T_Range<size_t> totrange(range.First()+range2.First(), range.First()+range2.Next());
      return make_unique<ParallelRangeVector> (orig, totrange);
    }
    
    AutoVector CreateVector () const override
    {
      cout << "CreateVector" << endl;
      throw Exception("CreateVector not avail");
    }
    
    void GetIndirect (FlatArray<int> ind, 
                      FlatVector<double> v) const override
    {
      local_vec -> GetIndirect (ind, v);
    }
    void GetIndirect (FlatArray<int> ind, 
                      FlatVector<Complex> v) const override
    {
      local_vec -> GetIndirect (ind, v);
    }

    void Cumulate () const override
    {
      orig->Cumulate();
      status = CUMULATED;
    }

    void Distribute() const override
    {
      orig->Distribute();
      status = DISTRIBUTED;
    }
    
    PARALLEL_STATUS GetParallelStatus () const override
    {

      return orig->GetParallelStatus();
    }
    
    virtual void SetParallelStatus (PARALLEL_STATUS stat) const override
    {
      if (stat != status)
        {
          cout << "setparallel, should changed from " << status << " to " << stat << endl;
          throw Exception("SetParallelStatus of rangevec called, change from"
                          +ToString(status) + " to " +ToString(stat));
        }
    }

    virtual void SetParallelDofs (shared_ptr<ParallelDofs> aparalleldofs) override
    {
      if (Size() != aparalleldofs->GetNDofLocal())
        throw Exception("SetParallelDofs of rangevec called, wrong sizes");
      paralleldofs = aparalleldofs;
    }


    void IRecvVec ( int dest, MPI_Request & request ) override
    {
      cout << "irecvec, and throw" << endl;
      throw Exception("ParallelRangeVector, don't know how to IRecVec");
    }
    
    void AddRecvValues( int sender ) override
    {
      cout << "addvecrec, and throw" << endl;
      throw Exception("ParallelRangeVector, don't know how to IRecVec");
    }
    
  };


  
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

  void ParallelBaseVector :: SetZero ()
  {
    local_vec->SetZero();
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
    static Timer t("ParallelVector::Add"); RegionTimer reg(t);
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
    static Timer t("ParallelVector - Cumulate");
    RegionTimer reg(t);
    
    // #ifdef PARALLEL
    if (status != DISTRIBUTED) return;
    
    // int ntasks = paralleldofs->GetNTasks();
    auto exprocs = paralleldofs->GetDistantProcs();
    
    int nexprocs = exprocs.Size();
    
    ParallelBaseVector * constvec = const_cast<ParallelBaseVector * > (this);
    
    for (int idest = 0; idest < nexprocs; idest ++ ) 
      constvec->ISend (exprocs[idest], sreqs[idest] );
    for (int isender=0; isender < nexprocs; isender++)
      constvec -> IRecvVec (exprocs[isender], rreqs[isender] );
    
    // if (rreqs.Size()) { // apparently Startall with 0 requests fails b/c invalid request ??
    //   MPI_Startall(rreqs.Size(), &rreqs[0]);
    //   MPI_Startall(sreqs.Size(), &sreqs[0]);
    // }

    MyMPI_WaitAll (sreqs);
    
    // cumulate
    for (int cntexproc=0; cntexproc < nexprocs; cntexproc++)
      {
	int isender = MyMPI_WaitAny (rreqs);
	constvec->AddRecvValues(exprocs[isender]);
      } 

    SetStatus(CUMULATED);
    // #endif
  }
  



  void ParallelBaseVector :: ISend ( int dest, MPI_Request & request ) const
  {
#ifdef PARALLEL
    MPI_Datatype mpi_t = this->paralleldofs->GetMPI_Type(dest);
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
  SCAL S_ParallelBaseVector<SCAL> :: InnerProduct (const BaseVector & v2, bool conjugate) const
  {
    static Timer t("ParallelVector - InnerProduct");
    RegionTimer reg(t);
    
    const ParallelBaseVector * parv2 = dynamic_cast_ParallelBaseVector(&v2);

    // two distributed vectors -- cumulate one
    if ( this->Status() == parv2->Status() && this->Status() == DISTRIBUTED )
      Cumulate();
    
    // two cumulated vectors -- distribute one
    else if ( this->Status() == parv2->Status() && this->Status() == CUMULATED )
      {
        // both cumulated
        if (this->EntrySize() == 1)
          {
            static Timer t("masked ip"); RegionTimer reg(t);
            FlatVector<> me = this->FVDouble();
            FlatVector<> you = parv2->FVDouble();
            const BitArray & ba = const_cast<BitArray&>(paralleldofs->MasterDofs());
	    SCAL localsum = MatKernelMaskedScalAB(me.Size(), me.Data(), 0, you.Data(), 0, ba);
	    if ( this->Status() == NOT_PARALLEL && parv2->Status() == NOT_PARALLEL )
	      { return localsum; }
	    return paralleldofs->GetCommunicator().AllReduce (localsum, MPI_SUM);
            // double sum = 0;
            // for (size_t i = 0; i < me.Size(); i++)
            //   if (ba.Test(i))
            //     sum += me(i) * you(i);
            // return paralleldofs->GetCommunicator().AllReduce (sum, MPI_SUM);
          }
        Distribute();
      }
    
    SCAL localsum = ngbla::InnerProduct (this->FVScal(), 
					 dynamic_cast<const S_BaseVector<SCAL>&>(*parv2).FVScal());

    if ( this->Status() == NOT_PARALLEL && parv2->Status() == NOT_PARALLEL )
      return localsum;

    return paralleldofs->GetCommunicator().AllReduce (localsum, MPI_SUM);
  }




  template <>
  Complex S_ParallelBaseVector<Complex> :: InnerProduct (const BaseVector & v2, bool conjugate) const
  {
    const ParallelBaseVector * parv2 = dynamic_cast_ParallelBaseVector(&v2);
  
    // two distributed vectors -- cumulate one
    if (this->Status() == parv2->Status() && this->Status() == DISTRIBUTED )
      Cumulate();

    // two cumulated vectors -- distribute one
    else if ( this->Status() == parv2->Status() && this->Status() == CUMULATED )
      Distribute();

    Complex localsum = conjugate ?
      ngbla::InnerProduct (Conj(FVComplex()), dynamic_cast<const S_BaseVector<Complex>&>(*parv2).FVComplex()) :
      ngbla::InnerProduct (FVComplex(), dynamic_cast<const S_BaseVector<Complex>&>(*parv2).FVComplex());

    // not parallel
    if ( this->Status() == NOT_PARALLEL && parv2->Status() == NOT_PARALLEL )
      return localsum;

    return paralleldofs->GetCommunicator().AllReduce (localsum, MPI_SUM);
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
  S_ParallelBaseVectorPtr<SCAL> :: 
  S_ParallelBaseVectorPtr (int as, int aes, void * adata,
			   shared_ptr<ParallelDofs> apd, PARALLEL_STATUS stat) throw()
    : S_BaseVectorPtr<SCAL> (as, aes, adata)
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
  SetParallelDofs (shared_ptr<ParallelDofs> aparalleldofs) // , const Array<int> * procs )
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

    // Initiate persistent send/recv requests for vector cumulate operation
    auto dps = paralleldofs->GetDistantProcs();
    this->sreqs.SetSize(dps.Size());
    this->rreqs.SetSize(dps.Size());
  }


  template <typename SCAL>
  void S_ParallelBaseVectorPtr<SCAL>  :: Distribute() const
  {
    if (status != CUMULATED) return;
    this->SetStatus(DISTRIBUTED);

    auto rank = paralleldofs->GetCommunicator().Rank();

    if (paralleldofs->GetEntrySize() == 1) {
      // this is a bit faster for BS = 1
      auto fv = this->template FV<SCAL>();
      for (auto p : paralleldofs->GetDistantProcs())
	if (p < rank)
	  for (auto dof : paralleldofs->GetExchangeDofs(p))
	    fv[dof] = 0;
    }
    else {
      for (auto p : paralleldofs->GetDistantProcs())
	if (p < rank)
	  for (auto dof : paralleldofs->GetExchangeDofs(p))
	    (*this)(dof) = 0;
    }

  }


  template < class SCAL >
  ostream & S_ParallelBaseVectorPtr<SCAL> :: Print (ostream & ost) const
  {
    this->PrintStatus (ost);
    S_BaseVectorPtr<SCAL>::Print (ost);
    return ost;
  }

  template < class SCAL >  
  AutoVector S_ParallelBaseVectorPtr<SCAL> :: Range (T_Range<size_t> range) const
  {
    /*
      Version 1: 
      A ParallelFlatVector
      + possible to cumulate / distribue by its own (???)
      - requies Range of pardofs, expensive ? 
      .... will provide method   Range (range, pardofs)
     */

    shared_ptr<ParallelDofs> sub_pardofs;
    if (paralleldofs)
      sub_pardofs = paralleldofs->Range(range);
    
    AutoVector locvec = S_BaseVectorPtr<SCAL>::Range (range);
    auto vec = make_unique<S_ParallelBaseVectorPtr<SCAL>> (range.Size(),
							   this->EntrySizeScal(),
                                                           locvec.Memory(), 
                                                           sub_pardofs,
                                                           this->GetParallelStatus());

    return std::move(vec);

    /*
      Version 2:
      Has pointer to original vector + Range
      - cumulate/distribute acts on original vector -> expensive
     */
    
    // return make_unique<ParallelRangeVector> (this, range);
  }


  template < class SCAL >  
  AutoVector S_ParallelBaseVectorPtr<SCAL> :: Range (DofRange range) const
  {
    AutoVector locvec = S_BaseVectorPtr<SCAL>::Range (range);
    auto vec = make_unique<S_ParallelBaseVectorPtr<SCAL>> (range.Size(),
							   this->EntrySizeScal(),
                                                           locvec.Memory(), 
                                                           range.GetParallelDofs(),
                                                           this->GetParallelStatus());
    return std::move(vec);
  }



  
  template <typename SCAL>
  void S_ParallelBaseVectorPtr<SCAL> :: IRecvVec ( int dest, MPI_Request & request )
  {
#ifdef PARALLEL
    MPI_Datatype MPI_TS = GetMPIType<TSCAL> ();
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
    return make_unique<S_ParallelBaseVectorPtr<TSCAL>>
      (this->size, this->es, paralleldofs, status);
    /*
    S_ParallelBaseVectorPtr<TSCAL> * parvec = 
      new S_ParallelBaseVectorPtr<TSCAL> (this->size, this->es, paralleldofs, status);
    return parvec;
    */
  }

  class ParallelMultiVector : public MultiVector
  {
    shared_ptr<ParallelDofs> paralleldofs;    
    
  public:
    ParallelMultiVector (shared_ptr<BaseVector> v, size_t cnt)
      : MultiVector(v, cnt)
    {
      paralleldofs = dynamic_cast<ParallelBaseVector&>(*v).GetParallelDofs();
    }

    void MakeSameStatus () const
    {
      if (Size() == 0) return;
      if ((*this)[0] -> GetParallelStatus() == DISTRIBUTED)
        for (int i = 1; i < Size(); i++)
          (*this)[i]->Distribute();

      if ((*this)[0] -> GetParallelStatus() == CUMULATED)
        for (int i = 1; i < Size(); i++)
          (*this)[i]->Cumulate();
    }
                         
    Matrix<> InnerProductD (const MultiVector & v2) const override
    {
      Matrix<double> res(v2.Size(), Size());
      if (res.Height()*res.Width()==0) return res;
      
      auto status1 = (*this)[0] -> GetParallelStatus();
      auto status2 = (v2)[0] -> GetParallelStatus();

      if (auto pv2 = dynamic_cast<const ParallelMultiVector*>(&v2))
        if ( ((status1 == CUMULATED) && (status2 == DISTRIBUTED)) ||
             ((status2 == CUMULATED) && (status1 == DISTRIBUTED)) )
          {
            ArrayMem<double*,32> ptrs1(Size());
            ArrayMem<double*,32> ptrs2(v2.Size());
            MakeSameStatus();
            pv2->MakeSameStatus();
            
            for (int i = 0; i < Size(); i++)
              ptrs1[i] = (*this)[i]->FVDouble().Data();
            for (int i = 0; i < v2.Size(); i++)
              ptrs2[i] = v2[i]->FVDouble().Data();
            
            ngbla::PairwiseInnerProduct((*this)[0]->FVDouble().Size(), ptrs2, ptrs1, res);
            
            paralleldofs->GetCommunicator()
              .AllReduce(FlatArray(res.Height()*res.Width(), res.Data()), MPI_SUM);
            
            return res;
          }

      // fallback
      return MultiVector::InnerProductD(v2);
    }

    
  };

  
  template <typename SCAL>  
  unique_ptr<MultiVector> S_ParallelBaseVectorPtr<SCAL> ::
  CreateMultiVector (size_t cnt) const
  {
    return make_unique<ParallelMultiVector> (CreateVector(), cnt);
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
      
    // double globsum = MyMPI_AllReduce (sum, MPI_SUM, paralleldofs->GetCommunicator()); // ngs_comm);
    double globsum = paralleldofs->GetCommunicator().AllReduce (sum, MPI_SUM); 
    return sqrt (globsum);
  }


  AutoVector CreateParallelVector (shared_ptr<ParallelDofs> pardofs, PARALLEL_STATUS status)
  {
    if (!pardofs)
      throw Exception ("CreateParallelVector called, but pardofs is nulltpr");

    // taken this version from the py-interface, maybe we should always create a parallel-vector 
    // #ifdef PARALLEL
    if(pardofs->IsComplex())
      return make_unique<S_ParallelBaseVectorPtr<Complex>> (pardofs->GetNDofLocal(), pardofs->GetEntrySize(), pardofs, status);
    else
      return make_unique<S_ParallelBaseVectorPtr<double>> (pardofs->GetNDofLocal(), pardofs->GetEntrySize(), pardofs, status);
    // #else
    // return CreateBaseVector (pardofs->GetNDofLocal(), pardofs->IsComplex(), pardofs->GetEntrySize());
    // #endif
  }

  
  template class S_ParallelBaseVectorPtr<double>;
  template class S_ParallelBaseVectorPtr<Complex>;
}



// #endif

