#ifndef UNIFIED_VECTOR_HPP
#define UNIFIED_VECTOR_HPP


namespace ngla
{

  class UnifiedVector : public S_BaseVector<double>
  {
  protected:
    double * host_data;
    double * dev_data;
  
    mutable bool host_uptodate;
    mutable bool dev_uptodate;
    
  public:
    UnifiedVector () = default;
    UnifiedVector (int asize);
    UnifiedVector (const BaseVector & vec);
    UnifiedVector (const UnifiedVector & vec);
    UnifiedVector (UnifiedVector && vec);

    virtual ~UnifiedVector();

    BaseVector & operator= (double d);
    BaseVector & operator= (const BaseVector & v2);
    UnifiedVector & operator= (const UnifiedVector & v2);    

    template <typename T2>
    UnifiedVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }

    const double & operator [] (const int ind) const;
    double & operator [] (const int ind);

    virtual AutoVector Range (T_Range<size_t> range) const;

    virtual BaseVector & Scale (double scal);
    virtual BaseVector & SetScalar (double scal);
    virtual BaseVector & Set (double scal, const BaseVector & v);
    virtual BaseVector & Add (double scal, const BaseVector & v);

    virtual double InnerProductD (const BaseVector & v2) const;
    virtual double L2Norm () const;

    void UpdateHost () const;
    void UpdateDevice () const;
    void InvalidateHost() const { host_uptodate = false; }
    void InvalidateDevice() const { dev_uptodate = false; }
    bool IsHostUptodate() const { return host_uptodate; }
    bool IsDevUptodate() const { return dev_uptodate; }
    
    virtual ostream & Print (ostream & ost) const;    
    virtual ostream & PrintStatus (ostream & ost) const;
    /* virtual void PrintDevice () const; */
    virtual AutoVector CreateVector () const;

    virtual FlatVector<double> FVDouble () const;
    virtual FlatVector<Complex> FVComplex () const;
    virtual void * Memory() const throw ();

    virtual double* DevData() const
    { 
      return dev_data; 
    }
    
    virtual double* HostData() const
    {
      return host_data;
    }

    friend class DevDMatrix;
    friend class DevSparseMatrix;
    friend class DevJacobiPrecond;
  };

  class UnifiedVectorWrapper : public UnifiedVector
  {
    mutable bool initial_host_uptodate;
    mutable bool initial_dev_uptodate;
    const BaseVector & vec;

  public:
    UnifiedVectorWrapper(const BaseVector & vec_, optional<IntRange> opt_range = nullopt);
    ~UnifiedVectorWrapper();

    using UnifiedVector::operator=;
  };



}
#endif
