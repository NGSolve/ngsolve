#ifndef FILE_AVECTOR
#define FILE_AVECTOR



namespace ngbla
{

#ifdef NONEALL

  template <typename T> struct FooAVectorType
  {
    typedef FlatVector<T> flattype;
    typedef Vector<T> type;
  };

  template <typename T, ORDERING ORD> struct FooAMatrixType
  {
    typedef FlatMatrix<T> flattype;
    typedef Matrix<T,ORD> type;
  };

  class AFlatVectorD;
  // class AFlatVectorC;
  
  class AFlatMatrixD;
  // class AFlatMatrixC;
  class AVectorD;
  class AMatrixD;
  // class AMatrixC;
  class AMatrixDCol;

  template <> struct FooAVectorType<double>
  {
    typedef AFlatVectorD flattype;
    typedef AVectorD type;
  };

  /*
  template <> struct FooAVectorType<Complex>
  {
    typedef AFlatVectorC flattype;
    // typedef AVectorC type;
  };
  */
  
  template <> struct FooAMatrixType<double,RowMajor>
  {
    typedef AFlatMatrixD flattype;
    typedef AMatrixD type;
  };
  
  template <> struct FooAMatrixType<double,ColMajor>
  {
    // typedef AFlatMatrixD flattype;
    typedef AMatrixDCol type;
  };

  /*
  template <> struct FooAMatrixType<Complex,RowMajor>
  {
    typedef AFlatMatrixC flattype;
    typedef AMatrixC type;
  };
  */

  
  template <typename T = double> 
  using AFlatVector = typename FooAVectorType<T>::flattype;
  template <typename T = double, ORDERING ORD = RowMajor> 
  using AFlatMatrix = typename FooAMatrixType<T,ORD>::flattype;

  template <typename T = double> 
  using AVector = typename FooAVectorType<T>::type;
  template <typename T = double, ORDERING ORD = RowMajor> 
  using AMatrix = typename FooAMatrixType<T,ORD>::type;


  // template <typename T = double> class ABareVector;
  template <typename T = double> class ABareMatrix;
  template <typename T = double> class ASliceMatrix;
  template <typename T = double> class ABareSliceMatrix;




  template <typename T>
  class SIMDExpr : public Expr<T> { };





  /* *************************** SumExpr **************************** */

  template <class TA, class TB> 
  class AVXSumExpr : public SIMDExpr<AVXSumExpr<TA,TB> >
  {
    const TA & a;
    const TB & b;
  public:

    enum { IS_LINEAR = TA::IS_LINEAR && TB::IS_LINEAR };
    enum { IS_LINEAR_VEC = TA::IS_LINEAR_VEC && TB::IS_LINEAR_VEC };
    
    INLINE AVXSumExpr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }

    INLINE auto operator() (size_t i) const { return a(i)+b(i); }
    INLINE auto operator() (size_t i, size_t j) const { return a(i,j)+b(i,j); }
    INLINE auto Get(size_t i) const { return a.Get(i)+b.Get(i); }
    INLINE auto Get(size_t i, size_t j) const { return a.Get(i,j)+b.Get(i,j); } 

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }

    void Dump (ostream & ost) const
    { ost << "("; a.Dump(ost); ost << ") + ("; b.Dump(ost); ost << ")"; }
  };

  template <typename TA, typename TB>
  INLINE AVXSumExpr<TA, TB>
  operator+ (const SIMDExpr<TA> & a, const SIMDExpr<TB> & b)
  {
    return AVXSumExpr<TA, TB> (a.Spec(), b.Spec());
  }


  /* *************************** SubExpr **************************** */

  template <class TA, class TB> 
  class AVXSubExpr : public SIMDExpr<AVXSubExpr<TA,TB> >
  {
    const TA & a;
    const TB & b;
  public:

    enum { IS_LINEAR = TA::IS_LINEAR && TB::IS_LINEAR };
    enum { IS_LINEAR_VEC = TA::IS_LINEAR_VEC && TB::IS_LINEAR_VEC };
    
    INLINE AVXSubExpr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }

    INLINE auto operator() (size_t i) const { return a(i)-b(i); }
    INLINE auto operator() (size_t i, size_t j) const { return a(i,j)-b(i,j); }
    INLINE auto Get(size_t i) const { return a.Get(i)-b.Get(i); }
    INLINE auto Get(size_t i, size_t j) const { return a.Get(i,j)-b.Get(i,j); } 

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }

    void Dump (ostream & ost) const
    { ost << "("; a.Dump(ost); ost << ") - ("; b.Dump(ost); ost << ")"; }
  };

  template <typename TA, typename TB>
  INLINE AVXSubExpr<TA, TB>
  operator- (const SIMDExpr<TA> & a, const SIMDExpr<TB> & b)
  {
    return AVXSubExpr<TA, TB> (a.Spec(), b.Spec());
  }


  
 /* *************************** NegExpr **************************** */

  template <class TA> 
  class AVXNegExpr : public SIMDExpr<AVXNegExpr<TA> >
  {
    const TA & a;
  public:

    enum { IS_LINEAR = TA::IS_LINEAR };
    enum { IS_LINEAR_VEC = TA::IS_LINEAR_VEC };
    
    INLINE AVXNegExpr (const TA & aa) : a(aa) { ; }

    INLINE auto operator() (size_t i) const { return -a(i); }
    INLINE auto operator() (size_t i, size_t j) const { return -a(i,j); }
    INLINE auto Get(size_t i) const { return -a.Get(i); }
    INLINE auto Get(size_t i, size_t j) const { return -a.Get(i,j); } 

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }

    void Dump (ostream & ost) const
    { ost << "-("; a.Dump(ost); ost << ")"; }
  };

  template <typename TA>
  INLINE AVXNegExpr<TA>
  operator- (const SIMDExpr<TA> & a)
  {
    return AVXNegExpr<TA> (a.Spec());
  }



  /* *************************** PW_Mult_Expr **************************** */

  template <class TA, class TB> 
  class AVXPW_Mult_Expr : public SIMDExpr<AVXPW_Mult_Expr<TA,TB> >
  {
    const TA & a;
    const TB & b;
  public:

    enum { IS_LINEAR = TA::IS_LINEAR && TB::IS_LINEAR };
    enum { IS_LINEAR_VEC = TA::IS_LINEAR_VEC && TB::IS_LINEAR_VEC };
    
    INLINE AVXPW_Mult_Expr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }

    INLINE auto operator() (size_t i) const { return a(i)*b(i); }
    INLINE auto operator() (size_t i, size_t j) const { return a(i,j)*b(i,j); }
    auto Get(size_t i) const { return a.Get(i)*b.Get(i); }
    auto Get(size_t i, size_t j) const { return a.Get(i,j)*b.Get(i,j); } 

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }

    void Dump (ostream & ost) const
    { ost << "("; a.Dump(ost); ost << ") + ("; b.Dump(ost); ost << ")"; }
  };

  template <typename TA, typename TB>
  INLINE AVXPW_Mult_Expr<TA, TB>
  pw_mult (const SIMDExpr<TA> & a, const SIMDExpr<TB> & b)
  {
    return AVXPW_Mult_Expr<TA, TB> (a.Spec(), b.Spec());
  }


  /* *************************** ScaleExpr **************************** */

  template <class TA> 
  class AVXScaleExpr : public SIMDExpr<AVXScaleExpr<TA> >
  {
    double s1;
    SIMD<double> s;
    const TA & a;
  public:

    enum { IS_LINEAR = TA::IS_LINEAR };
    enum { IS_LINEAR_VEC = TA::IS_LINEAR_VEC };

    INLINE AVXScaleExpr (double as, const TA & aa) : s1(as), s(as), a(aa) { ; }

    INLINE auto operator() (size_t i) const  { return s1*a(i); }
    INLINE auto operator() (size_t i, size_t j) const  { return s1*a(i,j); }
    INLINE auto Get(size_t i) const { return s * a.Get(i); }
    INLINE auto Get(size_t i, size_t j) const { return s * a.Get(i,j); }

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }

    void Dump (ostream & ost) const
    { ost << s << "*("; a.Dump(ost); ost << ")"; }
  };

  template <typename TA>
  INLINE AVXScaleExpr<TA>
  operator* (double s, const SIMDExpr<TA> & a)
  {
    return AVXScaleExpr<TA> (s, a.Spec());
  }
#endif

#ifdef NONE

  class AFlatVectorD : public SIMDExpr<AFlatVectorD>
  {
  protected:
    size_t size; 
    size_t vsize() const { return (size+SIMD<double>::Size()-1) / SIMD<double>::Size(); }
    SIMD<double> * __restrict data;  
  public:
    AFlatVectorD(size_t asize, LocalHeap & lh)
    {
      size = asize;
      data = lh.Alloc<SIMD<double>> (vsize());
    }

    AFlatVectorD (size_t as, SIMD<double> * adata)
      : size(as), data(adata)
    { ; }

    AFlatVectorD (size_t as, double * adata)
      : size(as), data((SIMD<double>*)(void*)adata)
    { ; }

    enum { IS_LINEAR = true };
    enum { IS_LINEAR_VEC = true };
    
    size_t Size () const { return size; }
    size_t VSize () const { return vsize(); }
    size_t Height () const { return size; }
    size_t Width () const { return 1; }
  
    double & operator() (size_t i) const
    {
      return ((double*)data)[i]; 
    }

    double & operator() (size_t i, size_t j) const
    {
      return ((double*)data)[i]; 
    }

    SIMD<double> & Get(size_t i) const { return data[i]; }

    INLINE const AFlatVectorD & operator= (const AFlatVectorD & v2) const
    {
      for (size_t i = 0; i < vsize(); i++)
        data[i] = v2.data[i];
      return *this;
    }
  
    AFlatVectorD & operator= (double d)
    {
      for (size_t i = 0; i < vsize(); i++)
        data[i] = SIMD<double> (d);
      return *this;
    }

    template<typename TB>
    INLINE const AFlatVectorD & operator= (const Expr<TB> & v) const
    {
      for (size_t i = 0; i < size; i++)
        ((double*)data)[i] = v.Spec()(i);
      return *this;
    }

    template<typename TB>
    INLINE const AFlatVectorD & operator= (const SIMDExpr<TB> & v) const
    {
      // for (size_t i = 0; i < vsize(); i++)
      for (auto i : ngstd::Range(vsize()))       
        data[i] = v.Spec().Get(i);
      return *this;
    }

    template<typename TB>
    INLINE const AFlatVectorD & operator+= (const SIMDExpr<TB> & v) const
    {
      // for (size_t i = 0; i < vsize(); i++)
      for (auto i : ngstd::Range(vsize()))
        data[i] += v.Spec().Get(i);
      return *this;
    }

    FlatVector<> Range (size_t begin, size_t end) const
    {
      return FlatVector<> (end-begin, ((double*)data)+begin);
    }
  };

  class AVectorD : public AFlatVectorD
  {
  public:
    AVectorD (size_t as)
    // : AFlatVectorD (as, (double*) _mm_malloc(sizeof(SIMD<double>) * ((as+SIMD<double>::Size()-1) / SIMD<double>::Size()), 64)) { ; }
      : AFlatVectorD (as, new SIMD<double>[(as+SIMD<double>::Size()-1) / SIMD<double>::Size()]) { ; }
    AVectorD () : AFlatVectorD (0, (double*)nullptr) { ; }
    AVectorD (AVectorD && av2) : AFlatVectorD (0, (double*)nullptr) { Swap(size, av2.size); Swap (data, av2.data); }
    ~AVectorD()
    {
      // _mm_free(data);
      delete data;
    }
    using AFlatVectorD::operator=;
    AVectorD & operator= (AVectorD && av2) { Swap(size, av2.size); Swap (data, av2.data); return *this; }
  };


#endif
  

#ifdef NONE

  class AFlatVectorC : public SIMDExpr<AFlatVectorC>
  {
  protected:
    size_t size; 
    size_t vsize() const { return (size+SIMD<double>::Size()-1) / SIMD<double>::Size(); }
    SIMD<Complex> * __restrict data;  
  public:
    AFlatVectorC(size_t asize, LocalHeap & lh)
    {
      size = asize;
      data = lh.Alloc<SIMD<Complex>> (vsize());
    }

    AFlatVectorC (size_t as, SIMD<Complex> * adata)
      : size(as), data(adata)
    { ; }

    enum { IS_LINEAR = true };
    enum { IS_LINEAR_VEC = true };
    
    size_t Size () const { return size; }
    size_t VSize () const { return vsize(); }
    size_t Height () const { return size; }
    size_t Width () const { return 1; }

    /*
    double & operator() (size_t i) const
    {
      return ((double*)data)[i]; 
    }

    double & operator() (size_t i, size_t j) const
    {
      return ((double*)data)[i]; 
    }
    */
    
    SIMD<Complex> & Get(size_t i) const { return data[i]; }

    INLINE const AFlatVectorC & operator= (const AFlatVectorC & v2) const
    {
      for (size_t i = 0; i < vsize(); i++)
        data[i] = v2.data[i];
      return *this;
    }
  
    AFlatVectorC & operator= (double d)
    {
      for (size_t i = 0; i < vsize(); i++)
        data[i] = SIMD<double> (d);
      return *this;
    }
    
    AFlatVectorC & operator*= (double d)
    {
      for (size_t i = 0; i < vsize(); i++) data[i] *= d;
      return *this;
    }

    /*
    template<typename TB>
    INLINE const AFlatVectorD & operator= (const Expr<TB> & v) const
    {
      for (size_t i = 0; i < size; i++)
        ((double*)data)[i] = v.Spec()(i);
      return *this;
    }

    template<typename TB>
    INLINE const AFlatVectorD & operator= (const SIMDExpr<TB> & v) const
    {
      for (size_t i = 0; i < vsize(); i++)
        data[i] = v.Spec().Get(i);
      return *this;
    }

    template<typename TB>
    INLINE const AFlatVectorD & operator+= (const SIMDExpr<TB> & v) const
    {
      for (size_t i = 0; i < vsize(); i++)
        data[i] += v.Spec().Get(i);
      return *this;
    }
    */
    AFlatVectorC Range (size_t begin, size_t end) const
    {
      return AFlatVectorC (end-begin, data+begin);
    }
  };

  class AVectorC : public AFlatVectorC
  {
  public:
    AVectorC (size_t as)
      : AFlatVectorC (as, new SIMD<Complex>[(as+SIMD<double>::Size()-1) / SIMD<double>::Size()]) { ; }
    AVectorC () : AFlatVectorC (0, (SIMD<Complex>*)nullptr) { ; }
    AVectorC (AVectorC && av2) : AFlatVectorC (0, (SIMD<Complex>*)nullptr) { Swap(size, av2.size); Swap (data, av2.data); }
    ~AVectorC() { delete data; }
    using AFlatVectorC::operator=;
    AVectorC & operator= (AVectorC && av2) { Swap(size, av2.size); Swap (data, av2.data); return *this; }
  };

#endif
  



  

#ifdef NONE
  class AFlatMatrixD : public SIMDExpr<AFlatMatrixD>
  {
  protected:
    size_t h, w;
    SIMD<double> * __restrict data;
  public:
    AFlatMatrixD () = default;
    AFlatMatrixD (const AFlatMatrixD &) = default;
    AFlatMatrixD (size_t ah, size_t aw, LocalHeap & lh)
    {
      h = ah;
      w = aw;
      // data = lh.Alloc<SIMD<double>> (h* ((unsigned(w)+SIMD<double>::Size()-1)/SIMD<double>::Size()));
      data = lh.Alloc<SIMD<double>> (h* ((w+SIMD<double>::Size()-1)/SIMD<double>::Size()));
    }
    
    AFlatMatrixD(size_t ah, size_t aw, double * adata)
      : h(ah), w(aw), data((SIMD<double>*)(void*)adata) { ; } 
    
    AFlatMatrixD(size_t ah, size_t aw, SIMD<double> * mem)
      : h(ah), w(aw), data(mem) { ; } 

    explicit AFlatMatrixD (const FlatMatrix<SIMD<double>> & m2)
      : h(m2.Height()), w(m2.Width()*SIMD<double>::Size()), data(&m2(0,0)) { ; } 
    
    void AssignMemory (size_t ah, size_t aw, SIMD<double> * mem)
    {
      h = ah;
      w = aw;
      data = mem;
    }
    
    enum { IS_LINEAR = false };
    enum { IS_LINEAR_VEC = true };
    
    // size_t Size () const { return h*w; }
    size_t Height () const { return h; }
    size_t Width () const { return w; }
    // unsigned int VWidth() const { return (unsigned(w)+3)/4; }
    size_t VWidth() const { return (w+SIMD<double>::Size()-1)/SIMD<double>::Size(); }
    double & operator() (size_t i) const
    {
      return ((double*)data)[i]; 
    }

    double & operator() (size_t i, size_t j) const
    {
      size_t vw = VWidth(); // (w+3)/4;
      return ((double*)data)[SIMD<double>::Size()*i*vw+j]; 
    }

    SIMD<double> & Get(size_t i) const { return data[i]; }
    SIMD<double> & Get(size_t i, size_t j) const { return data[i*VWidth()+j]; }

    const AFlatMatrixD & operator= (const AFlatMatrixD & m2) const
    {
      size_t els = h*VWidth();
      SIMD<double> * hdata = data;
      SIMD<double> * hdata2 = m2.data;
      for (size_t i = 0; i < els; i++)
        hdata[i] = hdata2[i];
      return *this;
    }
  
    AFlatMatrixD & operator= (double d)
    {
      auto vw = VWidth(); //  (unsigned(w)+3)/4;
      size_t els = h*vw; 
      for (size_t i = 0; i < els; i++)
        data[i] = SIMD<double>(d);
      return *this;
    }

    template<typename TB>
    INLINE const AFlatMatrixD & operator= (const Expr<TB> & v) const
    {
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < w; j++)
          (*this)(i,j) = v.Spec()(i,j);
      return *this;
    }

    template<typename TB>
    INLINE const AFlatMatrixD & operator= (const SIMDExpr<TB> & v) const
    {
      if (TB::IS_LINEAR_VEC)
        for (size_t i = 0; i < h*VWidth(); i++)
          Get(i) = v.Spec().Get(i);
      else
        for (size_t i = 0; i < h; i++)
          for (size_t j = 0; j < VWidth(); j++)
            Get(i,j) = v.Spec().Get(i,j);
      return *this;
    }

    AFlatMatrixD & operator*= (double d)
    {
      size_t vw = VWidth(); //  (w+3)/4;
      for (size_t i = 0; i < h*vw; i++)
        data[i] *= SIMD<double>(d);
      return *this;
    }
  
    operator SliceMatrix<double> () const
    {
      size_t vw = VWidth(); // (w+3)/4;
      return SliceMatrix<double> (h, w, SIMD<double>::Size()*vw, (double*)data);
    }
   
    SliceVector<> Col (size_t c) const
    {
      size_t vw = VWidth(); // (w+3)/4;    
      return SliceVector<> (h, SIMD<double>::Size()*vw, ((double*)data)+c);
    }

    AFlatVector<double> Row (size_t r) const
    {
      // return AFlatVector<double> (w, (double*)&Get(r,0));
      // return AFlatVector<double> (w, &Get(r,0));
      return AFlatVector<double> (w, data+r*VWidth());
    }
  
    AFlatMatrixD Rows(size_t begin, size_t end) const
    {
      size_t vw = VWidth(); 
      return AFlatMatrixD(end-begin, w, data+begin*vw);
    }
    AFlatMatrixD Rows(IntRange r) const
    { return Rows(r.begin(), r.end()); }

    SliceMatrix<> Cols(size_t begin, size_t end) const
    {
      size_t vw = VWidth(); 
      return SliceMatrix<>(h, end-begin, SIMD<double>::Size()*vw, ((double*)data)+begin);
    }
    SliceMatrix<> Cols(IntRange r) const
    { return Cols(r.begin(), r.end()); }
    
    INLINE ABareSliceMatrix<> VCols (size_t begin, size_t end) const;
    operator BareSliceMatrix<SIMD<double>> () const { return BareSliceMatrix<SIMD<double>> (VWidth(), data, DummySize(Height(), VWidth())); }
  };

  class AMatrixD : public AFlatMatrixD
  {
  public:
    AMatrixD (size_t ah, size_t aw)
    // : AFlatMatrixD (ah, aw, (double*)_mm_malloc(sizeof(SIMD<double>)*ah* ((aw+SIMD<double>::Size()-1)/SIMD<double>::Size()), 64)) { ; }
      : AFlatMatrixD (ah, aw, new SIMD<double>[ah* ((aw+SIMD<double>::Size()-1)/SIMD<double>::Size())]) { ; }
    ~AMatrixD ()
    {
      // _mm_free (data);
      delete data;
    }
    using AFlatMatrixD::operator=;
  };
#endif



#ifdef NONE
  class AFlatMatrixDCol
  {
  protected:
    size_t h, w;
    SIMD<double> * __restrict data;
    
  public:
    AFlatMatrixDCol (size_t ah, size_t aw, SIMD<double> * adata)
      : h(ah), w(aw), data(adata) { ; }

    size_t VHeight() const { return (h+SIMD<double>::Size()-1)/SIMD<double>::Size(); }

    double & operator() (size_t i, size_t j) const
    {
      size_t vh = VHeight();
      return ((double*)data)[SIMD<double>::Size()*j*vh+i];
    }
      
    AFlatVector<double> Col (size_t i) const
    {
      return AFlatVector<double> (h, data+i*VHeight());
    }
    
    AFlatMatrixDCol Cols (size_t begin, size_t end) const
    {
      return AFlatMatrixDCol (h, end-begin, data+begin*VHeight());
    }
    AFlatMatrixDCol Cols (IntRange r) const
    { return Cols(r.begin(), r.end()); }

    SliceMatrix<double,ColMajor> Rows (size_t begin, size_t end) const
    {
      return SliceMatrix<double,ColMajor> (end-begin, w, VHeight()*SIMD<double>::Size(), ((double*)data)+begin);
    }
    
    AFlatMatrixDCol & operator= (double scal)
    {
      for (size_t i = 0; i < w*VHeight(); i++)
        data[i] = scal;
      return *this;
    }
  };
  
  // col major
  class AMatrixDCol : public AFlatMatrixDCol
  {
  public:
    AMatrixDCol (size_t ah, size_t aw)
      : AFlatMatrixDCol (ah, aw, (SIMD<double>*)_mm_malloc(sizeof(SIMD<double>)*aw* ((ah+SIMD<double>::Size()-1)/SIMD<double>::Size()), 64)) { ; }
      
    ~AMatrixDCol ()
    {
      _mm_free (data);
    }

    using AFlatMatrixDCol::operator=;
  };
#endif
  




#ifdef NONE
 
  class AFlatMatrixC : public SIMDExpr<AFlatMatrixC>
  {
  protected:
    size_t h, w;
    SIMD<Complex> * __restrict data;
  public:
    AFlatMatrixC () = default;
    AFlatMatrixC (const AFlatMatrixC &) = default;
    AFlatMatrixC (size_t ah, size_t aw, LocalHeap & lh)
    {
      h = ah;
      w = aw;
      data = lh.Alloc<SIMD<Complex>> (h* ((w+SIMD<double>::Size()-1)/SIMD<double>::Size()));
    }
    
    AFlatMatrixC(size_t ah, size_t aw, SIMD<Complex> * mem)
      : h(ah), w(aw), data(mem) { ; } 
    
    void AssignMemory (size_t ah, size_t aw, SIMD<Complex> * mem)
    {
      h = ah;
      w = aw;
      data = mem;
    }
    
    enum { IS_LINEAR = false };
    enum { IS_LINEAR_VEC = true };
    
    // size_t Size () const { return h*w; }
    size_t Height () const { return h; }
    size_t Width () const { return w; }
    // unsigned int VWidth() const { return (unsigned(w)+3)/4; }
    size_t VWidth() const { return (w+SIMD<double>::Size()-1)/SIMD<double>::Size(); }
    /*
    double & operator() (size_t i) const
    {
      return ((double*)data)[i]; 
    }

    double & operator() (size_t i, size_t j) const
    {
      size_t vw = VWidth(); // (w+3)/4;
      return ((double*)data)[SIMD<double>::Size()*i*vw+j]; 
    }
    */

    SIMD<Complex> & Get(size_t i) const { return data[i]; }
    SIMD<Complex> & Get(size_t i, size_t j) const { return data[i*VWidth()+j]; }

    /*
    const AFlatMatrixD & operator= (const AFlatMatrixD & m2) const
    {
      size_t els = h*VWidth();
      SIMD<double> * hdata = data;
      SIMD<double> * hdata2 = m2.data;
      for (size_t i = 0; i < els; i++)
        hdata[i] = hdata2[i];
      return *this;
    }

    AFlatMatrixD & operator= (double d)
    {
      auto vw = VWidth(); //  (unsigned(w)+3)/4;
      size_t els = h*vw; 
      for (size_t i = 0; i < els; i++)
        data[i] = SIMD<double>(d);
      return *this;
    }

    template<typename TB>
    INLINE const AFlatMatrixD & operator= (const Expr<TB> & v) const
    {
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < w; j++)
          (*this)(i,j) = v.Spec()(i,j);
      return *this;
    }

    template<typename TB>
    INLINE const AFlatMatrixD & operator= (const SIMDExpr<TB> & v) const
    {
      if (TB::IS_LINEAR_VEC)
        for (size_t i = 0; i < h*VWidth(); i++)
          Get(i) = v.Spec().Get(i);
      else
        for (size_t i = 0; i < h; i++)
          for (size_t j = 0; j < VWidth(); j++)
            Get(i,j) = v.Spec().Get(i,j);
      return *this;
    }

    AFlatMatrixD & operator*= (double d)
    {
      size_t vw = VWidth(); //  (w+3)/4;
      for (size_t i = 0; i < h*vw; i++)
        data[i] *= SIMD<double>(d);
      return *this;
    }
    */

    /*
    operator SliceMatrix<double> () const
    {
      size_t vw = VWidth(); // (w+3)/4;
      return SliceMatrix<double> (h, w, SIMD<double>::Size()*vw, (double*)data);
    }
   
    SliceVector<> Col (size_t c) const
    {
      size_t vw = VWidth(); // (w+3)/4;    
      return SliceVector<> (h, SIMD<double>::Size()*vw, ((double*)data)+c);
    }
    */
    AFlatVector<Complex> Row (size_t r) const
    {
      return AFlatVector<Complex> (w, &Get(r,0));
    }
  
    AFlatMatrixC Rows(size_t begin, size_t end) const
    {
      size_t vw = VWidth(); 
      return AFlatMatrixC(end-begin, w, data+begin*vw);
    }
    /*
    AFlatMatrixD Rows(IntRange r) const
    { return Rows(r.begin(), r.end()); }

    SliceMatrix<> Cols(size_t begin, size_t end) const
    {
      size_t vw = VWidth(); 
      return SliceMatrix<>(h, end-begin, SIMD<double>::Size()*vw, ((double*)data)+begin);
    }
    SliceMatrix<> Cols(IntRange r) const
    { return Cols(r.begin(), r.end()); }
    
    INLINE ABareSliceMatrix<> VCols (size_t begin, size_t end) const;
    */
    operator BareSliceMatrix<SIMD<Complex>> () const { return BareSliceMatrix<SIMD<Complex>> (VWidth(), data); }    
  };

  class AMatrixC : public AFlatMatrixC
  {
  public:
    AMatrixC (size_t ah, size_t aw)
      : AFlatMatrixC (ah, aw, new SIMD<Complex>[ah* ((aw+SIMD<double>::Size()-1)/SIMD<double>::Size())]) { ; }
    ~AMatrixC ()
    {
      delete data;
    }
    using AFlatMatrixC::operator=;
  };

#endif
  







  
  /*
    template <typename T>
    INLINE void AddABt (AFlatMatrix<T> a, AFlatMatrix<T> b, SliceMatrix<T> c)
    {
    c += a * Trans(b) | Lapack; 
    }
  */

  /*
  INLINE SliceMatrix<double,ColMajor> Trans (const AFlatMatrixD & mat)
  {
    return SliceMatrix<double,ColMajor> (mat.Width(), mat.Height(), SIMD<double>::Size()*mat.VWidth(), &mat(0,0));
  }
  */
  
#ifdef NONE

  template <>
  class ABareVector<double>
  {
    SIMD<double> * __restrict data;
  public:
    ABareVector(SIMD<double> * _data) : data(_data) { ; }
    ABareVector(AFlatVector<double> vec) : data(&vec.Get(0)) { ; } 
    ABareVector(const ABareVector &) = default;

    double & operator() (size_t i) const
    {
      return ((double*)data)[i]; 
    }
    
    double & operator() (size_t i, size_t j) const
    {
      return ((double*)data)[i]; 
    }
    SIMD<double> & Get(size_t i) const { return data[i]; }
    AFlatVector<> AddSize(size_t s) const { return AFlatVector<> (s, data); }
    AFlatVector<> AddVSize(size_t s) const { return AFlatVector<> (s*SIMD<double>::Size(), data); }
  };
#endif
  
  
#ifdef NONE
  template <>
  class ABareVector<Complex>
  {
    SIMD<Complex> * __restrict data;
  public:
    ABareVector(SIMD<Complex> * _data) : data(_data) { ; }
    ABareVector(const ABareVector &) = default;

    SIMD<Complex> & Get(size_t i) const { return data[i]; }
    AFlatVector<Complex> AddSize(size_t s) const { return AFlatVector<Complex> (s, data); }
    AFlatVector<Complex> AddVSize(size_t s) const { return AFlatVector<Complex> (s*SIMD<Complex>::Size(), data); }
  };
#endif



#ifdef NONE
  template <>
  class ABareMatrix<double> : public DummySize
  {
    SIMD<double> * __restrict data;
    size_t dist;   // dist in simds
  public:
    enum { IS_LINEAR = false };
    enum { IS_LINEAR_VEC = true };
    
    ABareMatrix(SIMD<double> * _data, size_t _dist, size_t ah = -1, size_t aw = -1)
      : DummySize(ah, aw), data(_data), dist(_dist) { ; }
    ABareMatrix(AFlatMatrix<double> mat)
      : DummySize(mat.Height(), mat.Width()),
        data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
    ABareMatrix(const ABareMatrix &) = default;
    ABareMatrix & operator= (const ABareMatrix&) = delete;
    
    double & operator() (size_t i, size_t j) const
    { return ((double*)data)[SIMD<double>::Size()*i*dist+j]; }
    size_t Dist() const { return dist; }
    SIMD<double> & Get(size_t i, size_t j) const { return data[i*dist+j]; }
    SIMD<double> & Get(size_t i) const { return data[i]; }
    ABareVector<double> Row(size_t i) const { return ABareVector<double> (data+i*dist); }
    ABareMatrix<double> Rows(size_t first, size_t /* next */) const { return ABareMatrix<double> (data+first*dist, dist); }
    ABareMatrix<double> Rows(IntRange r) const { return Rows(r.First(), r.Next()); } 
    ABareMatrix<double> RowSlice(size_t first, size_t adist) const { return ABareMatrix<double> (data+first*dist, dist*adist); }
    operator BareSliceMatrix<SIMD<double>> () const { return BareSliceMatrix<SIMD<double>> (dist, data, *this); }
  };
  

  template <>
  class ASliceMatrix<double> : public SIMDExpr<ASliceMatrix<double>>
  {
    SIMD<double> * __restrict data;
    size_t dist;   // dist in simds
    size_t h, w;
  public:
    enum { IS_LINEAR = false };
    enum { IS_LINEAR_VEC = false };
    
    ASliceMatrix(size_t ah, size_t aw, size_t _dist, SIMD<double> * _data)
      : data(_data), dist(_dist), h(ah), w(aw) { ; }
    ASliceMatrix(AFlatMatrix<double> mat)
      : data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)), h(mat.Height()), w(mat.Width()) { ; }
    ASliceMatrix(const ASliceMatrix &) = default;

    auto & operator= (const ASliceMatrix<> & m2)
    {
      auto vw = VWidth(); 
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < vw; j++)
          data[i*dist+j] = m2.Get(i,j);
      return *this;
    }

    auto & operator= (const AFlatMatrix<> & m2)
    {
      *this = ASliceMatrix(m2);
      return *this;
    }
    
    auto & operator= (double d)
    {
      auto vw = VWidth(); 
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < vw; j++)
          data[i*dist+j] = SIMD<double>(d);
      return *this;
    }
    
    auto & operator*= (double d)
    {
      auto vw = VWidth(); 
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < vw; j++)
          data[i*dist+j] *= SIMD<double>(d);
      return *this;
    }

    template<typename TB>
    auto & operator= (const Expr<TB> & v) const
    {
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < w; j++)
          (*this)(i,j) = v.Spec()(i,j);
      return *this;
    }
    
    template<typename TB>
    auto & operator= (const SIMDExpr<TB> & v) const
    {
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < VWidth(); j++)
          Get(i,j) = v.Spec().Get(i,j);
      return *this;
    }

    size_t Height() const { return h; }
    size_t Width() const { return w; }
    size_t VWidth() const { return (w+SIMD<double>::Size()-1)/SIMD<double>::Size(); }
    
    double & operator() (size_t i, size_t j) const
    { return ((double*)data)[SIMD<double>::Size()*i*dist+j]; }
    size_t Dist() const { return dist; }
    SIMD<double> & Get(size_t i, size_t j) const { return data[i*dist+j]; }
    SIMD<double> & Get(size_t i) const { return data[i]; }
    AFlatVector<double> Row(size_t i) const { return AFlatVector<double> (w, data+i*dist); }
    ASliceMatrix<double> Rows(size_t first, size_t next) const { return ASliceMatrix<double> (next-first, w, dist, data+first*dist); }
    ASliceMatrix<double> Rows(IntRange r) const { return Rows(r.First(), r.Next()); } 
    // ASliceMatrix<double> RowSlice(size_t first, size_t adist) const { return ABareSliceMatrix<double> (data+first*dist, dist*adist); } 
  };


#endif
  
#ifdef NONE
  template <>
  class ASliceMatrix<Complex> : public SIMDExpr<ASliceMatrix<Complex>>
  {
    SIMD<Complex> * __restrict data;
    size_t dist;   // dist in simds
    size_t h, w;
  public:
    enum { IS_LINEAR = false };
    enum { IS_LINEAR_VEC = false };
    
    ASliceMatrix(size_t ah, size_t aw, size_t _dist, SIMD<Complex> * _data)
      : data(_data), dist(_dist), h(ah), w(aw) { ; }
    ASliceMatrix(AFlatMatrix<Complex> mat)
      : data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)), h(mat.Height()), w(mat.Width()) { ; }
    ASliceMatrix(const ASliceMatrix &) = default;

    auto & operator= (const ASliceMatrix & m2)
    {
      auto vw = VWidth(); 
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < vw; j++)
          data[i*dist+j] = SIMD<Complex> (m2.Get(i,j));
      return *this;
    }

    auto & operator= (const AFlatMatrix<Complex> & m2)
    {
      *this = ASliceMatrix(m2);
      return *this;
    }
    
    auto & operator= (Complex d)
    {
      auto vw = VWidth(); 
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < vw; j++)
          data[i*dist+j] = SIMD<Complex>(d);
      return *this;
    }
    
    auto & operator*= (double d)
    {
      auto vw = VWidth(); 
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < vw; j++)
          data[i*dist+j] *= SIMD<double>(d);
      return *this;
    }
    
    auto & operator*= (Complex d)
    {
      auto vw = VWidth(); 
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < vw; j++)
          data[i*dist+j] *= SIMD<Complex>(d);
      return *this;
    }

    /*
    template<typename TB>
    auto & operator= (const Expr<TB> & v) const
    {
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < w; j++)
          (*this)(i,j) = v.Spec()(i,j);
      return *this;
    }
    */
    template<typename TB>
    auto & operator= (const SIMDExpr<TB> & v) const
    {
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < VWidth(); j++)
          Get(i,j) = v.Spec().Get(i,j);
      return *this;
    }

    size_t Height() const { return h; }
    size_t Width() const { return w; }
    size_t VWidth() const { return (w+SIMD<double>::Size()-1)/SIMD<double>::Size(); }
    
    // double & operator() (size_t i, size_t j) const
    // { return ((double*)data)[SIMD<double>::Size()*i*dist+j]; }
    size_t Dist() const { return dist; }
    SIMD<Complex> & Get(size_t i, size_t j) const { return data[i*dist+j]; }
    SIMD<Complex> & Get(size_t i) const { return data[i]; }
    AFlatVector<Complex> Row(size_t i) const { return AFlatVector<Complex> (w, data+i*dist); }
    ASliceMatrix<Complex> Rows(size_t first, size_t next) const { return ASliceMatrix<Complex> (next-first, w, dist, data+first*dist); }
    ASliceMatrix<Complex> Rows(IntRange r) const { return Rows(r.First(), r.Next()); } 
    // ASliceMatrix<double> RowSlice(size_t first, size_t adist) const { return ABareSliceMatrix<double> (data+first*dist, dist*adist); } 
  };
#endif



#ifdef NONE
  template <>
  class ABareSliceMatrix<double> : public DummySize, public SIMDExpr<ABareSliceMatrix<double>>
  {
    SIMD<double> * __restrict data;
    size_t dist;   // dist in simds
  public:
    enum { IS_LINEAR = false };
    enum { IS_LINEAR_VEC = false };
    
    ABareSliceMatrix(SIMD<double> * _data, size_t _dist, size_t ah = -1, size_t aw = -1)
      : DummySize(ah, aw), data(_data), dist(_dist) { ; }
    ABareSliceMatrix(AFlatMatrix<double> mat)
      : DummySize(mat.Height(), mat.Width()),
        data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
    ABareSliceMatrix(ABareMatrix<double> mat)
      : DummySize(mat), data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
    ABareSliceMatrix(ASliceMatrix<double> mat)
      : DummySize(mat.Height(), mat.Width()),
        data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
    ABareSliceMatrix(FlatMatrix<SIMD<double>> mat)
      : DummySize(-1, -1), data(&mat(0,0)), dist(&mat(1,0)-&mat(0,0)) { ; } 
    ABareSliceMatrix(BareSliceMatrix<SIMD<double>> mat)
      : DummySize(-1, -1), data(&mat(0,0)), dist(mat.Dist()) { ; } 
    
    ABareSliceMatrix(const ABareSliceMatrix &) = default;

    ABareSliceMatrix & operator= (const ABareSliceMatrix&) = delete;
    double & operator() (size_t i, size_t j) const
    { return ((double*)data)[SIMD<double>::Size()*i*dist+j]; }
    size_t Dist() const { return dist; }
    ASliceMatrix<> AddSize(size_t h, size_t w) const
    { return ASliceMatrix<> (h,w,dist,data); }
    ASliceMatrix<> AddVSize(size_t h, size_t vw) const
    { return ASliceMatrix<> (h,SIMD<double>::Size()*vw,dist,data); }
    SIMD<double> & Get(size_t i, size_t j) const { return data[i*dist+j]; }
    SIMD<double> & Get(size_t i) const { return data[i]; }
    ABareVector<double> Row(size_t i) const { return ABareVector<double> (data+i*dist); }
    ABareSliceMatrix<double> Rows(size_t first, size_t /* next */) const { return ABareSliceMatrix<double> (data+first*dist, dist); }
    ABareSliceMatrix<double> Rows(IntRange r) const { return Rows(r.First(), r.Next()); } 
    ABareSliceMatrix<double> RowSlice(size_t first, size_t adist) const { return ABareSliceMatrix<double> (data+first*dist, dist*adist); }
    operator BareSliceMatrix<SIMD<double>> () const { return BareSliceMatrix<SIMD<double>> (dist, data, *this); }
  };
#endif

  /*
  template <>
  class ABareSliceMatrix<Complex> : public DummySize
  {
    SIMD<Complex> * __restrict data;
    size_t dist;   // dist in simds
  public:
    enum { IS_LINEAR = false };
    enum { IS_LINEAR_VEC = false };
    
    ABareSliceMatrix(AFlatMatrix<Complex> mat)
      : DummySize(mat.Height(), mat.Width()),
        data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
    ABareSliceMatrix(ASliceMatrix<Complex> mat)
      : DummySize(mat.Height(), mat.Width()),
        data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
    ABareSliceMatrix(FlatMatrix<SIMD<Complex>> mat)
      : DummySize(-1, -1), data(&mat(0,0)), dist(&mat(1,0)-&mat(0,0)) { ; }     
    ABareSliceMatrix(BareSliceMatrix<SIMD<Complex>> mat)
      : DummySize(-1, -1), data(&mat(0,0)), dist(mat.Dist()) { ; }     
    ABareSliceMatrix(SIMD<Complex> * _data, size_t _dist, size_t ah = -1, size_t aw = -1)
      : DummySize(ah, aw), data(_data), dist(_dist) { ; }

    ABareSliceMatrix & operator= (const ABareSliceMatrix&) = delete;
    // Complex operator() (size_t i, size_t j) const { return ((double*)data)[SIMD<double>::Size()*i*dist+j]; }
    size_t Dist() const { return dist; }
    // ASliceMatrix<> AddSize(size_t h, size_t w) const { return ASliceMatrix<> (h,w,dist,data); }
    ASliceMatrix<Complex> AddVSize(size_t h, size_t vw) const { return ASliceMatrix<Complex> (h,SIMD<double>::Size()*vw,dist,data); }
    SIMD<Complex> & Get (size_t i, size_t j) const { return data[i*dist+j]; }
    SIMD<Complex> & Get (size_t i) const { return data[i]; }
    ABareVector<Complex> Row (size_t i) const { return ABareVector<Complex> (data+i*dist); }
    ABareSliceMatrix<Complex> Rows (size_t first, size_t ) const { return ABareSliceMatrix<Complex> (data+first*dist, dist); }
    ABareSliceMatrix<Complex> Rows (IntRange r) const { return Rows(r.First(), r.Next()); } 
    ABareSliceMatrix<Complex> RowSlice (size_t first, size_t adist) const { return ABareSliceMatrix<Complex> (data+first*dist, dist*adist); }
    operator BareSliceMatrix<SIMD<Complex>> () const { return BareSliceMatrix<SIMD<Complex>> (dist, data); }

  };
*/



  /*
  ABareSliceMatrix<> AFlatMatrixD::VCols (size_t begin, size_t end) const
  { return ABareSliceMatrix<> (data+begin, VWidth(), Height(), Width()); }
  */


  template <typename TA, typename TB>
  void TransposeMatrix(SliceMatrix<TA> a, SliceMatrix<TB> b)
  {
    b = Trans(a);
  }

  /*
  */

  /*
  // c = a * Diag (diag)
  template <typename TA, typename TB, typename TC>
  void MultMatDiagMat(TA a, TB diag, TC c)
  {
    for (int i = 0; i < a.Width(); i++)
      c.Col(i) = diag(i) * a.Col(i);
  }
  */
  
  void TransposeMatrix(SliceMatrix<> a, SliceMatrix<> b);

#if defined(__AVX__) && !defined(__AVX512F__)

  extern void AddABt (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);
  extern void AddABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);

  //extern void AddABtSym (AFlatMatrix<double> a, AFlatMatrix<double> b, BareSliceMatrix<double> c);
  // extern void AddABtSym (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);
  extern void AddABtSym (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);
  extern void AddABtSym (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);

  extern void SubAtB (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c);
  inline void SubABt (SliceMatrix<double,ColMajor> a, SliceMatrix<double,ColMajor> b, SliceMatrix<double,ColMajor> c)
  {
    SubAtB (Trans(b), Trans(a), Trans(c));
  }
    
  // extern void MultMatDiagMat(AFlatMatrixD a, AFlatVectorD diag, AFlatMatrixD c);


#else // __AVX__

  // INLINE void AddABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  // { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }
  
  // INLINE void AddABtSym (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  // { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }

  INLINE void AddABt (SliceMatrix<double> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c)
  { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }
  
  INLINE void AddABtSym (SliceMatrix<double> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c)
  { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }

  INLINE void AddABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c)
  { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }
  
  INLINE void AddABtSym (SliceMatrix<Complex> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c)
  { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }

#endif // __AVX__


  
  template <typename TA, typename TB, typename TC, ORDERING ORD>
  INLINE void SubABt (const TA & a, const TB & b, SliceMatrix<TC,ORD> c)
  {
    c -= a * Trans(b) | Lapack;
    // LapackMultAdd (a, Trans(b), 1.0, c, 1.0);
  }

  
}
#endif // FILE_AVECTOR
