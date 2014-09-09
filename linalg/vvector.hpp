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
    
    NGS_DLL_HEADER virtual shared_ptr<BaseVector> Range (int begin, int end) const;
    NGS_DLL_HEADER virtual shared_ptr<BaseVector> Range (IntRange range) const;

    FlatVector<TSCAL> operator() (int i) const
    {
      return FlatVector<TSCAL> (es, pdata+i*es);
    }

    NGS_DLL_HEADER virtual shared_ptr<BaseVector> CreateVector () const;
    /*
    {
      switch (es)
	{
	case 1: return new VVector<TSCAL> (this->size);
	case 2: return new VVector<Vec<2,TSCAL> > (this->size);
	case 3: return new VVector<Vec<3,TSCAL> > (this->size);
	}
      return new S_BaseVectorPtr<TSCAL> (this->size, es);
    }
    */
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


  extern template class S_BaseVectorPtr<double>;
  extern template class S_BaseVectorPtr<Complex>;




  /**
     A specific vector based on Vector.
  */
  template <typename T = double>
  class VFlatVector :  virtual public S_BaseVectorPtr<typename mat_traits<T>::TSCAL>
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
  class  VVector : virtual public S_BaseVectorPtr<typename mat_traits<T>::TSCAL>
  {
  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { ES = sizeof(T) / sizeof(TSCAL) };

    explicit VVector (int as)
      : S_BaseVectorPtr<TSCAL> (as, ES) 
    { ; }

    // virtual ~VVector() { ; }

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



  extern template class VVector<double>;
  extern template class VVector<Complex>;
  
}
