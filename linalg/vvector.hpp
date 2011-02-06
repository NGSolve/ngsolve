namespace ngla
{

  template <class T> class VFlatVector;
  template <class T> class VVector;


  /**
     The T_BaseVector specifies the type of the vector element
  */
  template <typename T = double> 
  class NGS_DLL_HEADER T_BaseVector : virtual public S_BaseVector<typename mat_traits<T>::TSCAL>
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
  class NGS_DLL_HEADER  VFlatVector : public T_BaseVector<T>
  {
  protected:     
    T * data;

  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    explicit VFlatVector (int as, T * adata) throw(); 
    explicit VFlatVector () throw();
    virtual ~VFlatVector() throw();

    virtual void * Memory () const throw()
    {
      return data; 
    }

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
  };



  /**
     A specific vector based on Vector.
  */
  template <typename T = double> 
  class NGS_DLL_HEADER VVector : public T_BaseVector<T>
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

    T & operator() (int i) const
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







}
