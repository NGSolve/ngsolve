#ifndef FILE_NGS_PARALLELVECTOR
#define FILE_NGS_PARALLELVECTOR

/* ************************************************************************/
/* File:   parallelvector.hpp                                             */
/* Author: Astrid Sinwel, Joachim Schoeberl                               */
/* Date:   2007,2011                                                      */
/* ************************************************************************/



#ifdef PARALLEL

namespace ngla
{
  // using ngparallel::ParallelDofs;
  // using ngla::ParallelDofs;


  class NGS_DLL_HEADER ParallelBaseVector : virtual public BaseVector
  {
  protected:
    mutable PARALLEL_STATUS status;
    
  public:
    ParallelBaseVector ()
    { ; }

    template <typename T> 
    BaseVector & operator= (const VVecExpr<T> & v)
    {
      v.AssignTo (1.0, *this);
      return *this;
    }

    virtual PARALLEL_STATUS Status () const { return status; }

    virtual void SetStatus ( PARALLEL_STATUS astatus ) const
    {
      status = astatus;
    }

    virtual PARALLEL_STATUS GetParallelStatus () const { return Status(); }
    virtual void SetParallelStatus (PARALLEL_STATUS stat) const { SetStatus (stat); }



    virtual const ParallelDofs * GetParallelDofs () const
    {
      return paralleldofs; 
    }
    
    virtual bool IsParallelVector () const
    {
      return (this->Status() != NOT_PARALLEL);
    }
    
    virtual BaseVector & SetScalar (double scal);
    virtual BaseVector & SetScalar (Complex scal);
    
    virtual BaseVector & Set (double scal, const BaseVector & v);
    virtual BaseVector & Set (Complex scal, const BaseVector & v);

    virtual BaseVector & Add (double scal, const BaseVector & v);
    virtual BaseVector & Add (Complex scal, const BaseVector & v);

    void PrintStatus ( ostream & ost ) const;


    virtual void Cumulate () const; 
    
    virtual void Distribute() const = 0;
    // { cerr << "ERROR -- Distribute called for BaseVector, is not parallel" << endl; }
    
    virtual void ISend ( int dest, MPI_Request & request ) const;
    // virtual void Send ( int dest ) const;
    
    virtual void IRecvVec ( int dest, MPI_Request & request ) = 0;
    // { cerr << "ERROR -- IRecvVec called for BaseVector, is not parallel" << endl; }

    // virtual void RecvVec ( int dest )
    // { cerr << "ERROR -- IRecvVec called for BaseVector, is not parallel" << endl; }
    
    virtual void AddRecvValues( int sender ) = 0;
    // { cerr << "ERROR -- AddRecvValues called for BaseVector, is not parallel" << endl; }

    virtual void SetParallelDofs (const ParallelDofs * aparalleldofs, 
				  const Array<int> * procs = 0) = 0;
    /*
    { 
      if ( aparalleldofs == 0 ) return;
      cerr << "ERROR -- SetParallelDofs called for BaseVector, is not parallel" << endl; 
    }
    */
  };

  
  template <class SCAL>
  class NGS_DLL_HEADER S_ParallelBaseVector 
    : virtual public S_BaseVector<SCAL>, 
      virtual public ParallelBaseVector
  {
  protected:
    virtual SCAL InnerProduct (const BaseVector & v2) const;
  };


  template <class SCAL>
  class NGS_DLL_HEADER S_ParallelBaseVectorPtr
    : virtual public S_BaseVectorPtr<SCAL>, 
      virtual public S_ParallelBaseVector<SCAL>
  {
  protected:
    typedef SCAL TSCAL;
    using ParallelBaseVector :: status;
    using ParallelBaseVector :: paralleldofs;

    Table<SCAL> * recvvalues;

  public:
    // S_ParallelBaseVectorPtr (int as, int aes, void * adata) throw();
    S_ParallelBaseVectorPtr (int as, int aes, const ParallelDofs * apd, PARALLEL_STATUS stat) throw();

    virtual ~S_ParallelBaseVectorPtr ();
    virtual void SetParallelDofs (const ParallelDofs * aparalleldofs, const Array<int> * procs=0 );

    virtual void Distribute() const;
    virtual ostream & Print (ostream & ost) const;

    virtual void  IRecvVec ( int dest, MPI_Request & request );
    // virtual void  RecvVec ( int dest );
    virtual void AddRecvValues( int sender );
    virtual BaseVector * CreateVector () const;

    virtual double L2Norm () const;
  };
 



  template <typename T = double>
  class ParallelVVector : public VVector<T>, 
			  public S_ParallelBaseVectorPtr<typename mat_traits<T>::TSCAL>
  {
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { ES = sizeof(T) / sizeof(TSCAL) };

  public:
    explicit ParallelVVector (int as, const ParallelDofs * aparalleldofs,
			      PARALLEL_STATUS astatus = CUMULATED)
      : S_BaseVectorPtr<TSCAL> (as, ES), VVector<T> (as), 
	S_ParallelBaseVectorPtr<TSCAL> (as, ES, aparalleldofs, astatus)
    { ; }

    explicit ParallelVVector (const ParallelDofs * aparalleldofs,
			      PARALLEL_STATUS astatus = CUMULATED)
      : S_BaseVectorPtr<TSCAL> (aparalleldofs->GetNDofLocal(), ES), 
	VVector<T> (aparalleldofs->GetNDofLocal()), 
	S_ParallelBaseVectorPtr<TSCAL> (aparalleldofs->GetNDofLocal(), ES, aparalleldofs, astatus)
    { ; }


    virtual ~ParallelVVector() throw()
    { ; }
  };


  template <typename T = double>
  class ParallelVFlatVector : public VFlatVector<T>,
			      public S_ParallelBaseVectorPtr<typename mat_traits<T>::TSCAL>
  {
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { ES = sizeof(T) / sizeof(TSCAL) };

  public:
    explicit ParallelVFlatVector (int as, T * adata, 
				  ParallelDofs * aparalleldofs, 
				  PARALLEL_STATUS astatus)
    : S_BaseVectorPtr<TSCAL> (as, ES, adata),
      VFlatVector<T> (as, adata),
      S_ParallelBaseVectorPtr<TSCAL> (as, ES, aparalleldofs, astatus)
    { ; }

    explicit ParallelVFlatVector ()
    : S_BaseVectorPtr<TSCAL> (0, ES, NULL),
      S_ParallelBaseVectorPtr<TSCAL> (0, ES, NULL)
    { ; }
      
    virtual ~ParallelVFlatVector() throw()
    { ; }
  };
}

#endif
#endif
