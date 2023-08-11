namespace ngla
{


  template <class T> class VFlatVector;
  template <class T> class VVector;




  template <typename TSCAL = double> 
  class NGS_DLL_HEADER S_BaseVectorPtr : virtual public S_BaseVector<TSCAL>
  {
  protected:
    TSCAL * pdata;
    int es;
    bool ownmem;
    
  public:
    S_BaseVectorPtr (size_t as, int aes, void * adata) throw()
    {
      this->size = as;
      es = aes;
      pdata = static_cast<TSCAL*> (adata);
      ownmem = false;
      this->entrysize = es * sizeof(TSCAL) / sizeof(double);
    }

    S_BaseVectorPtr (size_t as, int aes)
    {
      this->size = as;
      es = aes;
      pdata = new TSCAL[as*aes];
      ownmem = true;
      GetMemoryTracer().Alloc(sizeof(TSCAL) * as * aes);
      this->entrysize = es * sizeof(TSCAL) / sizeof(double);
    }

    void SetSize (size_t as)
    {
      if (ownmem)
        {
          GetMemoryTracer().Free(sizeof(TSCAL) * this->size * es);
          delete [] pdata;
        }
      this->size = as;
      pdata = new TSCAL[as*es];
      ownmem = true;
      GetMemoryTracer().Alloc(sizeof(TSCAL) * as * es);
    }

    void AssignMemory (size_t as, void * adata)
    {
      this->size = as; 
      this->pdata = static_cast<TSCAL*> (adata); 
    }
    
    virtual ~S_BaseVectorPtr ();

    virtual void * Memory () const throw() override
    {
      return pdata; 
    }

    virtual Array<MemoryUsage> GetMemoryUsage () const override;

    
    // virtual AutoVector Range (size_t begin, size_t end) const override;
    virtual AutoVector Range (T_Range<size_t> range) const override;

    template<typename TIND, typename std::enable_if<std::is_integral<TIND>::value, int>::type = 0>
    FlatVector<TSCAL> operator() (TIND i) const
    {
      return FlatVector<TSCAL> (es, pdata+i*es);
    }

    virtual AutoVector CreateVector () const override;
    virtual unique_ptr<MultiVector> CreateMultiVector (size_t cnt) const override;
    
    virtual ostream & Print (ostream & ost) const override;
    using BaseVector::GetMemoryTracer;
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

    explicit VFlatVector (size_t as, T * adata) throw()
      : S_BaseVectorPtr<TSCAL> (as, ES, (void*)adata)
    { ; }

    VFlatVector (FlatVector<T> fv) 
      : S_BaseVectorPtr<TSCAL> (fv.Size(), ES, fv.Data())
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

    FlatVector<T> FV () const throw()
    {
      return FlatVector<T> (this->size, this->pdata);
    }
    
    T & operator() (size_t i) const
    {
      return static_cast<T*> (static_cast<void*> (this->pdata))[i];
    }

  };






  /**
     A specific vector based on Vector.
  */
  template <typename T = double> 
  class VVector : virtual public S_BaseVectorPtr<typename mat_traits<T>::TSCAL>
  {
  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { ES = sizeof(T) / sizeof(TSCAL) };

    explicit VVector (size_t as)
      : S_BaseVectorPtr<TSCAL> (as, ES) 
    { ; }

    explicit VVector (const VVector & v2)
      : S_BaseVectorPtr<TSCAL> (v2.Size(), ES)
    {
      *this = v2;
    }
    
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

    using S_BaseVectorPtr<TSCAL> :: FV;

    FlatVector<T> FV () const throw()
    {
      return FlatVector<T> (this->size, (T*)this->pdata);
    }

    T & operator() (size_t i) const
    {
      return static_cast<T*> (static_cast<void*> (this->pdata))[i];
    }
  };

  extern template class VVector<double>;
  extern template class VVector<Complex>;
}
