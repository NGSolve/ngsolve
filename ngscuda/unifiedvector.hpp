#ifndef UNIFIED_VECTOR_HPP
#define UNIFIED_VECTOR_HPP


namespace ngla
{

  class UnifiedVector : public S_BaseVector<double>
  {
    double * host_data;
    double * dev_data;
    cusparseDnVecDescr_t descr;
  
    mutable bool host_uptodate;
    mutable bool dev_uptodate;
    
  public:
    UnifiedVector (int asize);
    UnifiedVector (const BaseVector & vec);
    // UnifiedVector (UnifiedVector && vec);

    virtual ~UnifiedVector();

    BaseVector & operator= (double d);
    BaseVector & operator= (const BaseVector & v2);

    template <typename T2>
    UnifiedVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }

    const double & operator [] (const int ind) const;
    double & operator [] (const int ind);

    const cusparseDnVecDescr_t& GetDescr() const;

    cusparseDnVecDescr_t& GetDescr();

    virtual BaseVector & Scale (double scal);
    virtual BaseVector & SetScalar (double scal);
    virtual BaseVector & Set (double scal, const BaseVector & v);
    virtual BaseVector & Add (double scal, const BaseVector & v);

    double InnerProduct (const BaseVector & v2, bool conjugate=false) const;

    void UpdateHost () const;
    void UpdateDevice () const;

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
    
    friend class DevDMatrix;
    friend class DevSparseMatrix;
    friend class DevJacobiPrecond;
  };



}
#endif
