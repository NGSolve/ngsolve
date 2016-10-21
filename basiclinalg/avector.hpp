#ifndef FILE_AVECTOR
#define FILE_AVECTOR



namespace ngbla
{
  
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
  class AFlatMatrixD;
  class AVectorD;
  class AMatrixD;
  class AMatrixDCol;

  template <> struct FooAVectorType<double>
  {
    typedef AFlatVectorD flattype;
    typedef AVectorD type;
  };

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


  template <typename T = double> 
  using AFlatVector = typename FooAVectorType<T>::flattype;
  template <typename T = double, ORDERING ORD = RowMajor> 
  using AFlatMatrix = typename FooAMatrixType<T,ORD>::flattype;

  template <typename T = double> 
  using AVector = typename FooAVectorType<T>::type;
  template <typename T = double, ORDERING ORD = RowMajor> 
  using AMatrix = typename FooAMatrixType<T,ORD>::type;


  template <typename T = double> class ABareVector;
  template <typename T = double> class ABareMatrix;




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

    FlatVector<> Range (size_t begin, size_t end) const
    {
      return FlatVector<> (end-begin, ((double*)data)+begin);
    }
  };

  class AVectorD : public AFlatVectorD
  {
  public:
    AVectorD (size_t as)
      : AFlatVectorD (as, (double*) _mm_malloc(sizeof(SIMD<double>) * ((as+SIMD<double>::Size()-1) / SIMD<double>::Size()), 64)) { ; }
    AVectorD () : AFlatVectorD (0, (double*)nullptr) { ; }
    AVectorD (AVectorD && av2) : AFlatVectorD (0, (double*)nullptr) { Swap(size, av2.size); Swap (data, av2.data); }
    ~AVectorD()
    {
      _mm_free(data);
    }
    using AFlatVectorD::operator=;
    AVectorD & operator= (AVectorD && av2) { Swap(size, av2.size); Swap (data, av2.data); return *this; }
  };

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
    
    void AssignMemory (size_t ah, size_t aw, SIMD<double> * mem)
    {
      h = ah;
      w = aw;
      data = mem;
    }
    
    enum { IS_LINEAR = false };
    enum { IS_LINEAR_VEC = true };
    
    //total size of data (in doubles)
    size_t AllocSize() const { return h*VWidth()*SIMD<double>::Size(); }

    size_t Size () const { return h*w; }
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
      size_t vw = VWidth();
      for (size_t i = 0; i < h*vw; i++)
        data[i] = m2.data[i];
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
      return AFlatVector<double> (w, (double*)&Get(r,0));
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
  };

  class AMatrixD : public AFlatMatrixD
  {
  public:
    AMatrixD (size_t ah, size_t aw)
      : AFlatMatrixD (ah, aw, (double*)_mm_malloc(sizeof(SIMD<double>)*ah* ((aw+SIMD<double>::Size()-1)/SIMD<double>::Size()), 64)) { ; }
    ~AMatrixD ()
    {
      _mm_free (data);
    }
    using AFlatMatrixD::operator=;
  };


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
  
  

  /*
    template <typename T>
    INLINE void AddABt (AFlatMatrix<T> a, AFlatMatrix<T> b, SliceMatrix<T> c)
    {
    c += a * Trans(b) | Lapack; 
    }
  */

  INLINE SliceMatrix<double,ColMajor> Trans (const AFlatMatrixD & mat)
  {
    return SliceMatrix<double,ColMajor> (mat.Width(), mat.Height(), SIMD<double>::Size()*mat.VWidth(), &mat(0,0));
  }





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
  };

  template <>
  class ABareMatrix<double>
  {
    SIMD<double> * __restrict data;
    size_t dist;   // dist in simds
  public:
    ABareMatrix(SIMD<double> * _data, size_t _dist) : data(_data), dist(_dist) { ; }
    ABareMatrix(AFlatMatrix<double> mat) : data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
    ABareMatrix(const ABareMatrix &) = default;

    double & operator() (size_t i, size_t j) const
    {
      return ((double*)data)[SIMD<double>::Size()*i*dist+j]; 
    }
    size_t Dist() const { return dist; }
    SIMD<double> & Get(size_t i, size_t j) const { return data[i*dist+j]; }
    ABareVector<double> Row(size_t i) const { return ABareVector<double> (data+i*dist); }
    ABareMatrix<double> Rows(size_t first, size_t /* next */) const { return ABareMatrix<double> (data+first*dist, dist); }
    ABareMatrix<double> Rows(IntRange r) const { return Rows(r.First(), r.Next()); } 
    ABareMatrix<double> RowSlice(size_t first, size_t adist) const { return ABareMatrix<double> (data+first*dist, dist*adist); } 
  };









  template <typename TA, typename TB>
  void TransposeMatrix(SliceMatrix<TA> a, SliceMatrix<TB> b)
  {
    b = Trans(a);
  }

  template <typename TA, typename TB, typename TC>
  void MultMatMat(SliceMatrix<TA> a, SliceMatrix<TB> b, SliceMatrix<TC> c)
  {
    c = a * b;
  }

  // c = a * Diag (diag)
  template <typename TA, typename TB, typename TC>
  void MultMatDiagMat(TA a, TB diag, TC c)
  {
    for (int i = 0; i < a.Width(); i++)
      c.Col(i) = diag(i) * a.Col(i);
  }


#if defined(__AVX__)

  void TransposeMatrix(SliceMatrix<> a, SliceMatrix<> b);
  extern void MultMatMat(SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c);

  extern void AddABt (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c);
  extern void AddABt (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);
  extern void AddABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);

  extern void AddABtSym (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c);
  extern void AddABtSym (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);
  extern void AddABtSym (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);

  extern void SubABt (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c);
  extern void SubAtB (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c);
  inline void SubABt (SliceMatrix<double,ColMajor> a, SliceMatrix<double,ColMajor> b, SliceMatrix<double,ColMajor> c)
  {
    SubAtB (Trans(b), Trans(a), Trans(c));
  }
    
  extern void MultMatDiagMat(AFlatMatrixD a, AFlatVectorD diag, AFlatMatrixD c);


#else // __AVX__

  // template <typename T>
  // using AFlatMatrix = FlatMatrix<T>;



  /*
    class AFlatVectorD : public FlatVector<double>
    {
    public:
    AFlatVectorD (int as, LocalHeap & lh)
    : FlatVector<double> (as, lh) { ; }
    AFlatVectorD (int as, double * p)
    : FlatVector<double> (as, p) { ; }
    AFlatVectorD (int as, SIMD<double> * p)
    : FlatVector<double> (as, &p->Data()) { ; }
  
    AFlatVectorD & operator= (const AFlatVectorD & v2)  
    { FlatVector<double>::operator= (v2); return *this; }
    AFlatVectorD & operator= (double d)
    { FlatVector<double>::operator= (d); return *this; }
    template<typename TB>
    const AFlatVectorD & operator= (const Expr<TB> & v) const
    { FlatVector<double>::operator= (v); return *this; }
  
    int VSize() const { return size; }
    double & Get(int i) const { return (*this)[i]; }
    double & Get(int i) { return (*this)[i]; }  
    };

    class AFlatMatrixD : public FlatMatrix<double>
    {
    public:
    AFlatMatrixD () = default;  
    AFlatMatrixD (int ah, int aw, LocalHeap & lh)
    : FlatMatrix<double> (ah, aw, lh) { ; }
    AFlatMatrixD (int ah, int aw, double * p)
    : FlatMatrix<double> (ah, aw, p) { ; }
    AFlatMatrixD (int ah, int aw, SIMD<double> * p)
    : FlatMatrix<double> (ah, aw, &p->Data()) { ; }


    AFlatMatrixD & operator= (const AFlatMatrixD & m2)  
    { FlatMatrix<double>::operator= (m2); return *this; }
  
    AFlatMatrixD & operator= (double d)
    { FlatMatrix<double>::operator= (d); return *this; }

    template<typename TB>
    INLINE const AFlatMatrixD & operator= (const Expr<TB> & v) const
    { FlatMatrix<double>::operator= (v); return *this; }


    int VWidth() const { return Width(); }

    double & Get(int i) const { return (*this)(i); }
    double & Get(int i) { return (*this)(i); }  
    double & Get(int i, int j) const { return (*this)(i,j); }
    double & Get(int i, int j) { return (*this)(i,j); }  

    AFlatVectorD Row(int i) const
    { return AFlatVectorD(Width(), &(*this)(i,0)); } 
    AFlatMatrixD Rows(int i, int j) const
    { return AFlatMatrixD(j-i, Width(), &(*this)(i,0)); }
    AFlatMatrixD Rows(IntRange r) const
    { return Rows(r.begin(), r.end()); }
    };




    template <>
    class ABareVector<double>
    {
    double * __restrict data;
    public:
    ABareVector(double * _data) : data(_data) { ; }
    ABareVector(AFlatVector<double> vec) : data(&vec.Get(0)) { ; } 
    ABareVector(const ABareVector &) = default;

    double & operator() (int i) const
    {
    return ((double*)data)[i]; 
    }

    double & operator() (int i, int j) const
    {
    return ((double*)data)[i]; 
    }
    double & Get(int i) const { return data[i]; }
    };

    template <>
    class ABareMatrix<double>
    {
    double * __restrict data;
    int dist;   // dist in simds
    public:
    ABareMatrix(double * _data, int _dist) : data(_data), dist(_dist) { ; }
    ABareMatrix(SIMD<double> * _data, int _dist) : data(&_data->Data()), dist(_dist) { ; }
    ABareMatrix(AFlatMatrix<double> mat) : data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
    ABareMatrix(const ABareMatrix &) = default;

    double & operator() (int i, int j) const
    {
    return ((double*)data)[i*dist+j]; 
    }
    double & Get(int i, int j) const { return data[i*dist+j]; }
    ABareVector<double> Row(int i) const { return ABareVector<double> (data+i*dist); }
    ABareMatrix<double> Rows(int first, int next) const { return ABareMatrix<double> (data+first*dist, dist); }
    ABareMatrix<double> Rows(IntRange r) const { return Rows(r.First(), r.Next()); } 
    ABareMatrix<double> RowSlice(int first, int adist) const { return ABareMatrix<double> (data+first*dist, dist*adist); } 
    };
  */



  template <typename TA, typename TB, typename TC>
  INLINE void AddABt (const TA & a, const TB & b, SliceMatrix<TC> c)
  {
    c += a * Trans(b) | Lapack;
    // LapackMultAdd (a, Trans(b), 1.0, c, 1.0);
  }

  template <typename TA, typename TB, typename TC>
  INLINE void AddABtSym (const TA & a, const TB & b, SliceMatrix<TC> c)
  {
    c += a * Trans(b) | Lapack;
    // LapackMultAdd (a, Trans(b), 1.0, c, 1.0);
  }


#endif // __AVX__

  template <typename TA, typename TB, typename TC, ORDERING ORD>
  INLINE void SubABt (const TA & a, const TB & b, SliceMatrix<TC,ORD> c)
  {
    c -= a * Trans(b) | Lapack;
    // LapackMultAdd (a, Trans(b), 1.0, c, 1.0);
  }

}
#endif // FILE_AVECTOR
