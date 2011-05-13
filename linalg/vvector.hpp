namespace ngla
{


  template <class T> class VFlatVector;
  template <class T> class VVector;




  template <typename TSCAL = double> 
  class S_BaseVectorPtr : virtual public S_BaseVector<TSCAL>
  {
  protected:
    TSCAL * pdata;
    int es;
    bool ownmem;
    
  public:
    S_BaseVectorPtr (int as, int aes, void * adata) throw()
    {
      this->size = as;
      es = aes;
      pdata = static_cast<TSCAL*> (adata);
      ownmem = false;
      this->entrysize = es * sizeof(TSCAL) / sizeof(double);
    }

    S_BaseVectorPtr (int as, int aes) 
    {
      this->size = as;
      es = aes;
      pdata = new TSCAL[as*aes];
      ownmem = true;
      this->entrysize = es * sizeof(TSCAL) / sizeof(double);
    }

    void SetSize (int as) 
    {
      if (ownmem) delete [] pdata;
      this->size = as;
      pdata = new TSCAL[as*es];
      ownmem = true;
    }

    void AssignMemory (int as, void * adata)
    {
      this->size = as; 
      this->pdata = static_cast<TSCAL*> (adata); 
    }
    
    virtual ~S_BaseVectorPtr ()
    {
      if (ownmem) delete [] pdata;
    }

    virtual void * Memory () const throw()
    {
      return pdata; 
    }
    
    virtual BaseVector * Range (int begin, int end) const
    {
      return new S_BaseVectorPtr<TSCAL> (end-begin, es, pdata+begin*es);
    }

    virtual BaseVector * CreateVector ( const Array<int> * procs = 0) const
    {
      switch (es)
	{
	case 1: return new VVector<TSCAL> (this->size);
	case 2: return new VVector<Vec<2,TSCAL> > (this->size);
	case 3: return new VVector<Vec<3,TSCAL> > (this->size);
	}
      return new S_BaseVectorPtr<TSCAL> (this->size, es);
    }
  
    virtual ostream & Print (ostream & ost) const
    {
      // return (ost << FV() << endl);
      if (es == 1)
	ost << FlatVector<TSCAL> (this->size, pdata) << endl;
      else
	ost << FlatSysVector<TSCAL> (this->size, es, pdata);
      return ost;
    }
  };






  /**
     A specific vector based on Vector.
  */
  template <typename T = double>
  class VFlatVector :  public S_BaseVectorPtr<typename mat_traits<T>::TSCAL>
  {
  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { ES = sizeof(T) / sizeof(TSCAL) };

    explicit VFlatVector () throw()
      : S_BaseVectorPtr<TSCAL> (0, ES, NULL)
    { ; }

    explicit VFlatVector (int as, T * adata) throw()
      : S_BaseVectorPtr<TSCAL> (as, ES, (void*)adata)
    { ; }

    /*
    void AssignMemory (int as, void * adata)
    {
      this->size = as; 
      this->pdata = static_cast<T*> (adata); 
    }
    */

    FlatVector<T> FV () const throw()
    {
      return FlatVector<T> (this->size, this->pdata);
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






  /**
     A specific vector based on Vector.
  */
  template <typename T = double> 
  class  VVector : public S_BaseVectorPtr<typename mat_traits<T>::TSCAL>
  {
  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { ES = sizeof(T) / sizeof(TSCAL) };

    explicit VVector (int as)
      : S_BaseVectorPtr<TSCAL> (as, ES) 
    { ; }

    virtual ~VVector() { ; }

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

    template <typename T2>
    VVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }

    FlatVector<T> FV () const throw()
    {
      return FlatVector<T> (this->size, this->pdata);
    }

    T & operator() (int i) const
    {
      return static_cast<T*> (static_cast<void*> (this->pdata))[i];
    }
  };






  

 

#ifdef OLD



  template <class T> class VFlatVector;
  template <class T> class VVector;


  /**
     The T_BaseVector specifies the type of the vector element
  */
  template <typename T = double> 
  class NGS_DLL_HEADER T_BaseVector : virtual public S_BaseVector<typename mat_traits<T>::TSCAL>
  {
  protected:
    T * pdata;
    
  public:
    typedef T TELEM;
    typedef typename mat_traits<T>::TSCAL TSCAL;

    T_BaseVector () throw()
    { 
      this->entrysize = sizeof(T) / sizeof(double);
    }

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

    virtual void * Memory () const throw()
    {
      return pdata; 
    }

    FlatVector<T> FV () const throw()
    {
      return FlatVector<T> (this->size, pdata); // this->Memory());
    }

    FlatVector<T> FV (const BaseVector & v2) const
    {
      return FlatVector<T> (v2.Size(), v2.Memory());
    }

    T & operator() (int i) const
    {
      return pdata[i];
    }

    virtual BaseVector * Range (int begin, int end) const
    {
      return new VFlatVector<T> (end-begin, pdata+begin);
    }

    virtual BaseVector * CreateVector ( const Array<int> * procs = 0) const
    {
      return new VVector<T> (this->size);
    }
  
    virtual ostream & Print (ostream & ost) const
    {
      return (ost << FV() << endl);
    }
  };



  /**
     A specific vector based on Vector.
  */
  template <typename T = double>
  class NGS_DLL_HEADER  VFlatVector : public T_BaseVector<T>
  {
  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    explicit VFlatVector () throw()
    { 
      this->size = 0;
      this->pdata = NULL;
    }

    explicit VFlatVector (int as, T * adata) throw()
    { 
      this->size = as;
      this->pdata = adata;
    }

    virtual ~VFlatVector() throw()  
    { ; }

    void AssignMemory (int as, void * adata)
    {
      this->size = as; 
      this->pdata = static_cast<T*> (adata); 
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





  /**
     A specific vector based on Vector.
  */
  template <typename T = double> 
  class NGS_DLL_HEADER VVector : public T_BaseVector<T>
  {
  protected:
    DynamicMem<T> data;

  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    explicit VVector (int as)
      : data(as)
    { 
      this->size = as;
      this->pdata = &data[0]; 
      data.SetName ("VVector");
    }

    virtual ~VVector() throw()
    { 
      ; 
    }

    void SetSize(int as)
    {
      if (this->size == as) return;
      this->size = as;
      data.Alloc (this->size);
      this->pdata = &data[0]; 
    }

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

    template <typename T2>
    VVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }
  };

#endif


}
