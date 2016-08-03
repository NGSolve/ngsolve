#ifndef FILE_AVECTOR
#define FILE_AVECTOR




template <typename T> struct FooAMatrixType
{
  typedef FlatMatrix<T> type;
};

template <typename T> struct FooAVectorType
{
  typedef FlatVector<T> type;
};


class AFlatVectorD;
class AFlatMatrixD;

template <> struct FooAVectorType<double>
{
  typedef AFlatVectorD type;
};
template <> struct FooAMatrixType<double>
{
  typedef AFlatMatrixD type;
};




template <typename T = double> 
using AFlatVector = typename FooAVectorType<T>::type;
template <typename T = double> 
using AFlatMatrix = typename FooAMatrixType<T>::type;


template <typename T> class ABareVector;
template <typename T> class ABareMatrix;




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
    
  INLINE AVXSumExpr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }

  INLINE auto operator() (int i) const { return a(i)+b(i); }
  INLINE auto operator() (int i, int j) const { return a(i,j)+b(i,j); }
  INLINE auto Get(int i) const { return a.Get(i)+b.Get(i); } 

  INLINE int Height() const { return a.Height(); }
  INLINE int Width() const { return a.Width(); }

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
    
  INLINE AVXPW_Mult_Expr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }

  INLINE auto operator() (int i) const { return a(i)*b(i); }
  INLINE auto operator() (int i, int j) const { return a(i,j)*b(i,j); }
  auto Get(int i) const { return a.Get(i)*b.Get(i); } 

  INLINE int Height() const { return a.Height(); }
  INLINE int Width() const { return a.Width(); }

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

  INLINE AVXScaleExpr (double as, const TA & aa) : s1(as), s(as), a(aa) { ; }

  INLINE auto operator() (int i) const  { return s1*a(i); }
  INLINE auto operator() (int i, int j) const  { return s1*a(i,j); }
  INLINE auto Get(int i) const { return s * a.Get(i); }

  INLINE int Height() const { return a.Height(); }
  INLINE int Width() const { return a.Width(); }

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
  int size; 
  // unsigned vsize() const { return (unsigned(size)+3) / 4; }
  // __m256d * __restrict data;
  unsigned vsize() const { return (unsigned(size)+SIMD<double>::Size()-1) / SIMD<double>::Size(); }
  SIMD<double> * __restrict data;  
public:
  AFlatVectorD(int asize, LocalHeap & lh)
  {
    size = asize;
    data = lh.Alloc<SIMD<double>> (vsize());
  }

  AFlatVectorD (int as, double * adata)
    : size(as), data((SIMD<double>*)(void*)adata)
  { ; }

  enum { IS_LINEAR = true };
  int Size () const { return size; }
  int VSize () const { return vsize(); }
  int Height () const { return size; }
  int Width () const { return 1; }
  
  double & operator() (int i) const
  {
    return ((double*)data)[i]; 
  }

  double & operator() (int i, int j) const
  {
    return ((double*)data)[i]; 
  }

  SIMD<double> & Get(int i) const { return data[i]; }

  INLINE const AFlatVectorD & operator= (const AFlatVectorD & v2) const
  {
    for (int i = 0; i < vsize(); i++)
      data[i] = v2.data[i];
    return *this;
  }
  
  AFlatVectorD & operator= (double d)
  {
    for (int i = 0; i < vsize(); i++)
      data[i] = SIMD<double> (d);
    return *this;
  }

  template<typename TB>
  INLINE const AFlatVectorD & operator= (const Expr<TB> & v) const
  {
    for (int i = 0; i < size; i++)
      ((double*)data)[i] = v.Spec()(i);
    return *this;
  }

  template<typename TB>
  INLINE const AFlatVectorD & operator= (const SIMDExpr<TB> & v) const
  {
    for (int i = 0; i < vsize(); i++)
      data[i] = v.Spec().Get(i);
    return *this;
  }

  template<typename TB>
  INLINE const AFlatVectorD & operator+= (const SIMDExpr<TB> & v) const
  {
    for (auto i : Range(vsize()))
      data[i] += v.Spec().Get(i);
    return *this;
  }
};



class AFlatMatrixD : public SIMDExpr<AFlatMatrixD>
{
  int h, w;
  SIMD<double> * __restrict data;
public:
  AFlatMatrixD () = default;
  AFlatMatrixD (const AFlatMatrixD &) = default;
  AFlatMatrixD(int ah, int aw, LocalHeap & lh)
  {
    h = ah;
    w = aw;
    // int rw = (w+SIMD<double>::Size()-1)&(-SIMD<double>::Size());
    // data = (__m256d*)lh.Alloc<double> (h*rw);
    data = lh.Alloc<SIMD<double>> (h* ((unsigned(w)+SIMD<double>::Size()-1)/SIMD<double>::Size()));
  }

  AFlatMatrixD(int ah, int aw, double * adata)
    : h(ah), w(aw), data((SIMD<double>*)(void*)adata) { ; } 

  AFlatMatrixD(int ah, int aw, SIMD<double> * mem)
    : h(ah), w(aw), data(mem) { ; } 

  void AssignMemory (int ah, int aw, SIMD<double> * mem)
  {
    h = ah;
    w = aw;
    data = mem;
  }
  
  enum { IS_LINEAR = false };
  int Size () const { return h*w; }
  int Height () const { return h; }
  int Width () const { return w; }
  // unsigned int VWidth() const { return (unsigned(w)+3)/4; }
  unsigned int VWidth() const { return (unsigned(w)+SIMD<double>::Size()-1)/SIMD<double>::Size(); }
  double & operator() (int i) const
  {
    return ((double*)data)[i]; 
  }

  double & operator() (int i, int j) const
  {
    int vw = VWidth(); // (w+3)/4;
    return ((double*)data)[SIMD<double>::Size()*i*vw+j]; 
  }

  SIMD<double> & Get(int i) const { return data[i]; }
  SIMD<double> & Get(int i, int j) const { return data[i*VWidth()+j]; }

  const AFlatMatrixD & operator= (const AFlatMatrixD & m2) const
  {
    int vw = VWidth(); // (w+3)/4;
    for (int i = 0; i < h*vw; i++)
      data[i] = m2.data[i];
    return *this;
  }
  
  AFlatMatrixD & operator= (double d)
  {
    auto vw = VWidth(); //  (unsigned(w)+3)/4;
    int els = unsigned(h)*unsigned(vw);
    for (int i = 0; i < els; i++)
      data[i] = SIMD<double>(d);
    return *this;
  }


  template<typename TB>
  INLINE const AFlatMatrixD & operator= (const Expr<TB> & v) const
  {
    for (int i = 0; i < h; i++)
      for (int j = 0; j < w; j++)
        (*this)(i,j) = v.Spec()(i,j);
    return *this;
  }

  AFlatMatrixD & operator*= (double d)
  {
    int vw = VWidth(); //  (w+3)/4;
    for (int i = 0; i < h*vw; i++)
      data[i] *= SIMD<double>(d);
    return *this;
  }
  
  operator SliceMatrix<double> () const
  {
    int vw = VWidth(); // (w+3)/4;
    return SliceMatrix<double> (h, w, SIMD<double>::Size()*vw, (double*)data);
  }
   
  SliceVector<> Col (int c) const
  {
    int vw = VWidth(); // (w+3)/4;    
    return SliceVector<> (h, SIMD<double>::Size()*vw, ((double*)data)+c);
  }

  AFlatVector<double> Row (int r) const
  {
    return AFlatVector<double> (w, (double*)&Get(r,0));
  }
  
  AFlatMatrixD Rows(int begin, int end) const
  {
    int vw = VWidth(); // (w+3)/4;
    return AFlatMatrixD(end-begin, w, data+begin*vw);
  }
  AFlatMatrixD Rows(IntRange r) const
  { return Rows(r.begin(), r.end()); }

  SliceMatrix<> Cols(int begin, int end) const
  {
    int vw = VWidth(); // (w+3)/4;
    return SliceMatrix<>(h, end-begin, SIMD<double>::Size()*vw, ((double*)data)+begin);
  }
  SliceMatrix<> Cols(IntRange r) const
  { return Cols(r.begin(), r.end()); }


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

  double & operator() (int i) const
  {
    return ((double*)data)[i]; 
  }

  double & operator() (int i, int j) const
  {
    return ((double*)data)[i]; 
  }
  SIMD<double> & Get(int i) const { return data[i]; }
};

template <>
class ABareMatrix<double>
{
  SIMD<double> * __restrict data;
  int dist;   // dist in simds
public:
  ABareMatrix(SIMD<double> * _data, int _dist) : data(_data), dist(_dist) { ; }
  ABareMatrix(AFlatMatrix<double> mat) : data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
  ABareMatrix(const ABareMatrix &) = default;

  double & operator() (int i, int j) const
  {
    return ((double*)data)[SIMD<double>::Size()*i*dist+j]; 
  }
  SIMD<double> & Get(int i, int j) const { return data[i*dist+j]; }
  ABareVector<double> Row(int i) const { return ABareVector<double> (data+i*dist); }
  ABareMatrix<double> Rows(int first, int /* next */) const { return ABareMatrix<double> (data+first*dist, dist); }
  ABareMatrix<double> Rows(IntRange r) const { return Rows(r.First(), r.Next()); } 
  ABareMatrix<double> RowSlice(int first, int adist) const { return ABareMatrix<double> (data+first*dist, dist*adist); } 
};










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


#if defined(__AVX2__)


extern void MultMatMat(SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c);

extern void AddABt (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c);
extern void AddABt (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);
extern void AddABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);

extern void AddABtSym (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c);
extern void AddABtSym (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);
extern void AddABtSym (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);


extern void MultMatDiagMat(AFlatMatrixD a, AFlatVectorD diag, AFlatMatrixD c);


#else // __AVX2__

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


#endif // __AVX2__

#endif // FILE_AVECTOR
