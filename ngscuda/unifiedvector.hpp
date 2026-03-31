#ifndef UNIFIED_VECTOR_HPP
#define UNIFIED_VECTOR_HPP


namespace ngla
{

  class UnifiedVector : public S_BaseVector<double>
  {
  protected:
    double * host_data;
    Dev<double> * dev_data;
  
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

    virtual AutoVector Range (T_Range<size_t> range) const override;

    virtual BaseVector & Scale (double scal) override;
    virtual BaseVector & Scale (BaseScalar & scal) override;
    virtual BaseVector & SetScalar (double scal) override;
    virtual BaseVector & Set (double scal, const BaseVector & v) override;
    virtual BaseVector & Add (double scal, const BaseVector & v) override;

    using BaseVector::InnerProduct;
    virtual double InnerProductD (const BaseVector & v2) const override;
    virtual double L2Norm () const override;

    // BaseScalar overrides for GPU-resident scalars:
    virtual void InnerProduct (const BaseVector & v2, BaseScalar & scal, bool conjugate = false) const override;
    virtual BaseVector & Add (BaseScalar & scal, const BaseVector & v) override;
    virtual shared_ptr<BaseScalar> CreateScalar() const override;

    void UpdateHost () const;
    void UpdateDevice () const;
    void InvalidateHost() const { host_uptodate = false; }
    void InvalidateDevice() const { dev_uptodate = false; }
    bool IsHostUptodate() const { return host_uptodate; }
    bool IsDevUptodate() const { return dev_uptodate; }
    
    virtual ostream & Print (ostream & ost) const override;    
    virtual ostream & PrintStatus (ostream & ost) const;
    /* virtual void PrintDevice () const; */
    virtual AutoVector CreateVector () const override;

    virtual FlatVector<double> FVDouble () const override;
    virtual FlatVector<Complex> FVComplex () const override;
    virtual void * Memory() const throw () override;

    virtual Dev<double> * DevData() const
    { 
      return dev_data; 
    }
    
    virtual double* HostData() const
    {
      return host_data;
    }

    // udpate device and invalidate host
    FlatVector<Dev<double>> FVDev() const;

    // udpate device 
    FlatVector<Dev<double>> FVDevRO() const;
    
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



  // -------------------------------------------------------
  // UnifiedScalar: GPU-resident scalar for graph capture
  // Fixed GPU address allows CUDA graph capture.
  // -------------------------------------------------------
  class UnifiedScalar : public BaseScalar
  {
    mutable double host_val = 0.0;
    double* dev_val = nullptr;

  public:
    UnifiedScalar();
    virtual ~UnifiedScalar();

    void Set(double d) override;
    double GetD() const override;
    double* DevPtr() const { return dev_val; }

    template<typename Expr>
    UnifiedScalar& operator=(const Expr& expr)
    {
      double* d = dev_val;
      auto kernel = [d, expr] __device__ (size_t) {
          *d = expr.Evaluate();
      };
      DeviceParallelFor(1, kernel);
      return *this;
    }
  };

}
#endif
