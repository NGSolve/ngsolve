namespace ngla
{

  template <class T> class VFlatVector;
  template <class T> class VVector;


  /**
     The T_BaseVector specifies the type of the vector element
  */
  template <typename T = double> NGS_DLL_HEADER
  class T_BaseVector : public S_BaseVector<typename mat_traits<T>::TSCAL>
  {

  public:
    typedef T TELEM;
    typedef typename mat_traits<T>::TSCAL TSCAL;

    T_BaseVector & operator= (TSCAL s1)
    {
      BaseVector::operator= (s1);
      return *this;
    }

    T_BaseVector & operator= (const T_BaseVector & v2)
    {
      BaseVector::operator= (v2);
      return *this;
    }

    template <typename T2>
    T_BaseVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }



    FlatVector<T> FV () const throw()
    {
      return FlatVector<T> (this->size, this->Memory());
    }

    FlatVector<T> FV () throw()
    {
      return FlatVector<T> (this->size, this->Memory());
    }

    FlatVector<T> FV (const BaseVector & v2) const
    {
      return FlatVector<T> (v2.Size(), v2.Memory());
    }

  
    virtual BaseVector * Range (int begin, int end) const
    {
      return new VFlatVector<T> (end-begin, &FV()[0]+begin);
    }
  
    // create vector, procs is set of processors on which the vector exists
    // default 0 pointer means all procs
    virtual BaseVector * CreateVector ( const Array<int> * procs = 0) const;
  
    virtual ostream & Print (ostream & ost) const;

  };




  /**
     A specific vector based on Vector.
  */
  template <typename T = double>
  class VFlatVector : public T_BaseVector<T>
  {
  protected:     
    // int s;
    T * data;


  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    explicit VFlatVector (int as, T * adata) throw(); 
    explicit VFlatVector () throw();
    virtual ~VFlatVector() throw();

    virtual void * Memory () const throw()
    { return data; }

    void AssignMemory (int as, void * adata)
    {
      this->size = as; data = static_cast<T*> (adata); 
    }

    VFlatVector & operator= (TSCAL s1)
    {
      BaseVector::operator= (s1);
      return *this;
    }

    template <typename T2>
    VFlatVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }


  
#ifdef PARALLEL
    const T * Data() const { return data; }
    T * Data() { return data; }
#endif
  };



  /**
     A specific vector based on Vector.
  */
  template <typename T = double> NGS_DLL_HEADER
  class VVector : public T_BaseVector<T>
  {
  protected:
    DynamicMem<T> data;
  private:

  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    explicit VVector (int as);

    virtual ~VVector() throw();

    virtual void * Memory () const throw()
    { return const_cast<void*> (static_cast<const void*> (&data[0])); }

    void SetSize(int as);

    VVector & operator= (TSCAL s1)
    {
      BaseVector::operator= (s1);
      return *this;
    }

    VVector & operator= (const VVector & v2)
    {
      BaseVector::operator= (v2);
      return *this;
    }

    T & operator() (int i)
    {
      return data[i];
    }
    const T & operator() (int i) const
    {
      return data[i];
    }

    template <typename T2>
    VVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }


  };


#ifdef PARALLEL
  template <typename T = double>
  class ParallelVVector : public VVector<T>
  {
    ngparallel::ParallelDofs * paralleldofs;

    PARALLEL_STATUS status;
    Array<int> sendvector_size;
    Array<int> recvvector_size;
    // Array<T> * recvvalues;
    Table<T> * recvvalues;

  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    explicit ParallelVVector (int as);
    explicit ParallelVVector (int as, ngparallel::ParallelDofs * aparalleldofs ) ;
    explicit ParallelVVector (int as, ngparallel::ParallelDofs * aparalleldofs,
			      PARALLEL_STATUS astatus );

    virtual ~ParallelVVector() throw();
 
    virtual const PARALLEL_STATUS Status () const { return status; }

    virtual void SetStatus ( PARALLEL_STATUS astatus );

    virtual void SetParallelDofs ( ngparallel::ParallelDofs * aparalleldofs, const Array<int> * procs=0 );

    virtual class ngparallel::ParallelDofs * GetParallelDofs () const
    { return paralleldofs; }

    virtual void PrintParallelDofs() const;

    virtual bool IsParallelVector () const
    {
      if ( this->Status() == NOT_PARALLEL ) return false;
      return true;
    }

    /// values from reduceprocs are added up,
    /// vectors in sendtoprocs are set to the cumulated values
    /// default pointer 0 means send to proc 0
    virtual void AllReduce ( Array<int> * reduceprocs, Array<int> * sendtoprocs = 0 ) const;

    virtual void Distribute() const;

    virtual void PrintStatus ( ostream & ost ) const
    {
      if ( this->status == NOT_PARALLEL )
	ost << "NOT PARALLEL" << endl;
      else if ( this->status == DISTRIBUTED )
	ost << "DISTRIBUTED" << endl ;
      else if (this->status == CUMULATED )
	ost << "CUMULATED" << endl ;
    }

    virtual void ISend ( const int dest, MPI_Request & request );
    virtual void Send ( const int dest );

    virtual void  IRecvVec ( const int dest, MPI_Request & request );

    const T & RecvValue( int dest, int i ) const { return (*recvvalues)[dest][i] ; }

    T & RecvValue( int dest, int i ) { return (*recvvalues)[dest][i] ; }

    virtual BaseVector * CreateVector ( const Array<int> * procs = 0) const;

    virtual ostream & Print (ostream & ost) const;

    void AddRecvValues( int sender );

  };



  template <typename T = double>
  class ParallelVFlatVector : public VFlatVector<T>
  {
    ngparallel::ParallelDofs * paralleldofs;

    PARALLEL_STATUS status;

    Array<int> sendvector_size;
    Array<int> recvvector_size;
    //   Array<T> * recvvalues;
    Table<T> * recvvalues;


  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    explicit ParallelVFlatVector (int as, T * adata) ;
    explicit ParallelVFlatVector (int as, T * adata, ngparallel::ParallelDofs * aparalleldofs ) ;
    explicit ParallelVFlatVector (int as, T * adata, ngparallel::ParallelDofs * aparalleldofs, 
				  PARALLEL_STATUS astatus);
    explicit ParallelVFlatVector ();
    virtual ~ParallelVFlatVector() throw();

    virtual const PARALLEL_STATUS Status () const { return status; }

    virtual void SetStatus ( PARALLEL_STATUS astatus );

    virtual void SetParallelDofs ( ngparallel::ParallelDofs * aparalleldofs, const Array<int> * procs=0 );

    virtual class ngparallel::ParallelDofs * GetParallelDofs () const
    { return paralleldofs; }

    virtual void PrintParallelDofs() const;

    virtual bool IsParallelVector () const
    {
      if ( this->Status() == NOT_PARALLEL ) return false;
      return true;
    }

    /// values from reduceprocs are added up,
    /// vectors in sendtoprocs are set to the cumulated values
    /// default pointer 0 means send to proc 0
    virtual void AllReduce ( Array<int> * reduceprocs, Array<int> * sendtoprocs = 0 ) const;

    virtual void Distribute() const;

    virtual void PrintStatus ( ostream & ost ) const
    {
      if ( this->status == NOT_PARALLEL )
	ost << "NOT PARALLEL" << endl;
      else if ( this->status == DISTRIBUTED )
	ost << "DISTRIBUTED" << endl ;
      else if (this->status == CUMULATED )
	ost << "CUMULATED" << endl ;
    }

    virtual void ISend ( const int dest, MPI_Request & request );
    virtual void Send ( const int dest );

    virtual void IRecvVec ( const int dest, MPI_Request & request );

    const T & RecvValue( int dest, int i ) const
    { return (*recvvalues)[dest][i] ; }

    T & RecvValue( int dest, int i ) 
    { return (*recvvalues)[dest][i] ; }

    virtual BaseVector * CreateVector ( const Array<int> * procs = 0) const;

    virtual ostream & Print (ostream & ost) const;

    void AddRecvValues( int sender );
  };
#endif





  /*
    war wohl nur ein Versuch ???
    template <typename T>
    class VFlatVector<FlatVector<T> > : public T_BaseVector<T>
    {
    // int s;
    T * data;
    public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    explicit VFlatVector (int as, T * adata) throw(); 
    explicit VFlatVector () throw();
    virtual ~VFlatVector() throw();

    virtual void * Memory () const throw()
    { return data; }

    void AssignMemory (int as, void * adata)
    {
    this->size = as; data = static_cast<T*> (adata); 
    }

    virtual TempVector Range (int begin, int end)
    {
    return TempVector (new VFlatVector<FlatVector<T> > (end-begin, &FV()[0]+begin));
    }
  
    virtual TempVector Range (int begin, int end) const
    {
    return TempVector (new VFlatVector<FlatVector<T> > (end-begin, &FV()[0]+begin));
    }

    VFlatVector & operator= (TSCAL s1)
    {
    BaseVector::operator= (s1);
    return *this;
    }

    template <typename T2>
    VFlatVector & operator= (const VVecExpr<T2> & v)
    {
    BaseVector::operator= (v);
    return *this;
    }
    };


    template <typename T>
    class VVector<FlatVector<T> > : public T_BaseVector<T>
    {
    MoveableMem<T> data;
    int bs;
    public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    explicit VVector (int as, int abs);

    virtual ~VVector() throw();

    virtual void * Memory () const throw()
    { return const_cast<void*> (static_cast<const void*> (&data[0])); }

    void SetSize(int as, int abs);

    VVector & operator= (TSCAL s1)
    {
    BaseVector::operator= (s1);
    return *this;
    }

    VVector & operator= (const VVector & v2)
    {
    BaseVector::operator= (v2);
    return *this;
    }

    FlatVector<T> operator() (int i)
    {
    return FlatVector<> (bs, &data[i*bs]);
    }

    FlatVector<T> operator() (int i) const 
    {
    return FlatVector<> (bs, const_cast<T*> (&data[i*bs]));
    }

    template <typename T2>
    VVector & operator= (const VVecExpr<T2> & v)
    {
    BaseVector::operator= (v);
    return *this;
    }

    };

  */
}
