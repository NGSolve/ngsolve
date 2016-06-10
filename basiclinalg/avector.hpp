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


template <typename TA, typename TB, typename TC>
void MultMatMat(SliceMatrix<TA> a, SliceMatrix<TB> b, SliceMatrix<TC> c)
{
  c = a * b;
}

template <typename TA, typename TB, typename TC>
void MultMatDiagMat(TA a, TB b, TC c)
{
  for (int i = 0; i < a.Width(); i++)
    c.Col(i) = b(i) * a.Col(i);
}









#if defined(__AVX2__)

// mat-mat product

// b.Width <= 4
INLINE
void MultMatMat4(SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
{
  __m256i mask = _mm256_cmpgt_epi64(_mm256_set1_epi64x(b.Width()),
                                    _mm256_set_epi64x(3, 2, 1, 0));

  /*
  __m256i mask;
  switch (b.Width())
    {
    case 1:
      mask = _mm256_set_epi64x(0,0,0,-1); break;
    case 2:
      mask = _mm256_set_epi64x(0,0,-1,-1); break;
    case 3:
      mask = _mm256_set_epi64x(0,-1,-1,-1); break;
    case 4:
      mask = _mm256_set_epi64x(-1,-1,-1,-1); break;
    }
  */
  unsigned int da = a.Dist();
  int wa = a.Width();
  int r = 0;
  double * bpc = &c(0,0);
  unsigned int dc = c.Dist();
  double * ar = &a(0,0);
  for ( ; r < a.Height()-7; r+=8)
    {
      __m256d sum1 = _mm256_setzero_pd();
      __m256d sum2 = _mm256_setzero_pd();
      __m256d sum3 = _mm256_setzero_pd();
      __m256d sum4 = _mm256_setzero_pd();
      __m256d sum5 = _mm256_setzero_pd();
      __m256d sum6 = _mm256_setzero_pd();
      __m256d sum7 = _mm256_setzero_pd();
      __m256d sum8 = _mm256_setzero_pd();

      for (int j = 0; j < wa; j++)
        {
          __m256d rb =  _mm256_blendv_pd(_mm256_setzero_pd(),
                                         _mm256_loadu_pd(&b(j,0)),
                                         _mm256_castsi256_pd(mask));
          double * arj = ar + j;
          double * arj4 = arj + 4*da;
          sum1 += _mm256_set1_pd(*arj) * rb;
          sum2 += _mm256_set1_pd(*(arj+da)) * rb;
          sum3 += _mm256_set1_pd(*(arj+2*da)) * rb;
          sum4 += _mm256_set1_pd(*(arj+3*da)) * rb;
          sum5 += _mm256_set1_pd(*(arj4)) * rb;
          sum6 += _mm256_set1_pd(*(arj4+da)) * rb;
          sum7 += _mm256_set1_pd(*(arj4+2*da)) * rb;
          sum8 += _mm256_set1_pd(*(arj4+3*da)) * rb;
        }

      _mm256_maskstore_pd(bpc, mask, sum1);
      _mm256_maskstore_pd(bpc+dc, mask, sum2);
      bpc += 2*dc;
      _mm256_maskstore_pd(bpc, mask, sum3);
      _mm256_maskstore_pd(bpc+dc, mask, sum4);
      bpc += 2*dc;
      _mm256_maskstore_pd(bpc, mask, sum5);
      _mm256_maskstore_pd(bpc+dc, mask, sum6);
      bpc += 2*dc;
      _mm256_maskstore_pd(bpc, mask, sum7);
      _mm256_maskstore_pd(bpc+dc, mask, sum8);
      bpc += 2*dc;
      ar += 8*da;
    }

  if (r < a.Height()-3)
    {
      __m256d sum1 = _mm256_setzero_pd();
      __m256d sum2 = _mm256_setzero_pd();
      __m256d sum3 = _mm256_setzero_pd();
      __m256d sum4 = _mm256_setzero_pd();
      
      for (int j = 0; j < wa; j++)
        {
          __m256d rb =  _mm256_blendv_pd(_mm256_setzero_pd(),
                                         _mm256_loadu_pd(&b(j,0)),
                                         _mm256_castsi256_pd(mask));

          double * arj = ar + j;
          sum1 += _mm256_set1_pd(*arj) * rb;
          sum2 += _mm256_set1_pd(*(arj+da)) * rb;
          sum3 += _mm256_set1_pd(*(arj+2*da)) * rb;
          sum4 += _mm256_set1_pd(*(arj+3*da)) * rb;
        }

      _mm256_maskstore_pd(bpc, mask, sum1);
      _mm256_maskstore_pd(bpc+dc, mask, sum2);
      bpc += 2*dc;
      _mm256_maskstore_pd(bpc, mask, sum3);
      _mm256_maskstore_pd(bpc+dc, mask, sum4);
      bpc += 2*dc;
      r += 4;
      ar += 4*da;
    }
  if (r < a.Height()-1)
    {
      __m256d sum1 = _mm256_setzero_pd();
      __m256d sum2 = _mm256_setzero_pd();
      
      for (int j = 0; j < wa; j++)
        {
          __m256d rb =  _mm256_blendv_pd(_mm256_setzero_pd(),
                                         _mm256_loadu_pd(&b(j,0)),
                                         _mm256_castsi256_pd(mask));
          double * arj = ar + j;
          sum1 += _mm256_set1_pd(*arj) * rb;
          sum2 += _mm256_set1_pd(*(arj+da)) * rb;
        }

      _mm256_maskstore_pd(bpc + 0*dc, mask, sum1);
      _mm256_maskstore_pd(bpc + 1*dc, mask, sum2);
      bpc += 2*dc;
      r += 2;
      ar += 2*da;
    }

  if (r < a.Height())
    {
      __m256d sum = _mm256_setzero_pd();
      for (int j = 0; j < wa; j++)
        {
          __m256d rb =  _mm256_loadu_pd(&b(j,0));
          double * arj = ar + j;
          sum += _mm256_set1_pd(*arj) * rb;

        }

      _mm256_maskstore_pd(bpc + 0*dc, mask, sum);
    }
}

// b.Width() = 8
INLINE
void MultMatMat8(SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
{
  unsigned int da = a.Dist();
  int wa = a.Width();
  int r = 0;
  double * bpc = &c(0,0);
  unsigned int dc = c.Dist();
  double * ar = &a(0,0);
  for ( ; r < a.Height()-3; r+=4)
    {
      __m256d sum11 = _mm256_setzero_pd();
      __m256d sum21 = _mm256_setzero_pd();
      __m256d sum31 = _mm256_setzero_pd();
      __m256d sum41 = _mm256_setzero_pd();
      __m256d sum12 = _mm256_setzero_pd();
      __m256d sum22 = _mm256_setzero_pd();
      __m256d sum32 = _mm256_setzero_pd();
      __m256d sum42 = _mm256_setzero_pd();

      for (int j = 0; j < wa; j++)
        {
          __m256d rb1 = _mm256_loadu_pd(&b(j,0));
          __m256d rb2 = _mm256_loadu_pd(&b(j,4));

          double * arj = ar + j;

          __m256d a1 = _mm256_set1_pd(*arj);
          __m256d a2 = _mm256_set1_pd(*(arj+da));
          __m256d a3 = _mm256_set1_pd(*(arj+2*da));
          __m256d a4 = _mm256_set1_pd(*(arj+3*da));

          sum11 += a1 * rb1;
          sum21 += a2 * rb1;
          sum31 += a3 * rb1;
          sum41 += a4 * rb1;
          sum12 += a1 * rb2;
          sum22 += a2 * rb2;
          sum32 += a3 * rb2;
          sum42 += a4 * rb2;
        }

      _mm256_storeu_pd(bpc, sum11);
      _mm256_storeu_pd(bpc+4, sum12);
      bpc += dc;
      _mm256_storeu_pd(bpc, sum21);
      _mm256_storeu_pd(bpc+4, sum22);
      bpc += dc;
      _mm256_storeu_pd(bpc, sum31);
      _mm256_storeu_pd(bpc+4, sum32);
      bpc += dc;
      _mm256_storeu_pd(bpc, sum41);
      _mm256_storeu_pd(bpc+4, sum42);
      bpc += dc;
      ar += 4*da;
    }

  for ( ; r < a.Height()-1; r+=2)
    {
      __m256d sum11 = _mm256_setzero_pd();
      __m256d sum21 = _mm256_setzero_pd();
      __m256d sum12 = _mm256_setzero_pd();
      __m256d sum22 = _mm256_setzero_pd();

      for (int j = 0; j < wa; j++)
        {
          __m256d rb1 = _mm256_loadu_pd(&b(j,0));
          __m256d rb2 = _mm256_loadu_pd(&b(j,4));

          double * arj = ar + j;

          __m256d a1 = _mm256_set1_pd(*arj);
          __m256d a2 = _mm256_set1_pd(*(arj+da));

          sum11 += a1 * rb1;
          sum21 += a2 * rb1;
          sum12 += a1 * rb2;
          sum22 += a2 * rb2;
        }

      _mm256_storeu_pd(bpc, sum11);
      _mm256_storeu_pd(bpc+4, sum12);
      bpc += dc;
      _mm256_storeu_pd(bpc, sum21);
      _mm256_storeu_pd(bpc+4, sum22);
      bpc += dc;
      ar += 2*da;
    }

  for ( ; r < a.Height(); r++)
    {
      __m256d sum11 = _mm256_setzero_pd();
      __m256d sum12 = _mm256_setzero_pd();

      for (int j = 0; j < wa; j++)
        {
          __m256d rb1 = _mm256_loadu_pd(&b(j,0));
          __m256d rb2 = _mm256_loadu_pd(&b(j,4));

          double * arj = ar + j;

          __m256d a1 = _mm256_set1_pd(*arj);
          sum11 += a1 * rb1;
          sum12 += a1 * rb2;
        }

      _mm256_storeu_pd(bpc, sum11);
      _mm256_storeu_pd(bpc+4, sum12);
      bpc += dc;
      ar += da;
    }
}





// c = a * b
static
void MultMatMat(SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
{
  int k = 0;
  for ( ; k < b.Width()-7; k += 8)
    MultMatMat8(a, b.Cols(k,k+8), c.Cols(k,k+8));
  for ( ; k < b.Width(); k += 4)
    {
      int end = min2(b.Width(), k+4);
      MultMatMat4(a, b.Cols(k,end), c.Cols(k,end));
    }
}



/*
template <typename TA, typename TB, typename TC>
INLINE void AddABtSym (const TA & a, const TB & b, SliceMatrix<TC> c)
{
  c += a * Trans(b) | Lapack;
  // LapackMultAdd (a, Trans(b), 1.0, c, 1.0);
}

template <typename T>
INLINE void AddABt (AFlatMatrix<T> a, AFlatMatrix<T> b, SliceMatrix<T> c)
{
  c += a * Trans(b) | Lapack; 
}
*/


// aligned AVX vectors

INLINE ostream & operator<< (ostream & ost, __m256d v4)
{
  double * v = reinterpret_cast<double*>(&v4);
  ost << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3];
  return ost;
}

template <typename T>
class AVXExpr : public Expr<T> { };


/* *************************** SumExpr **************************** */

template <class TA, class TB> 
class AVXSumExpr : public AVXExpr<AVXSumExpr<TA,TB> >
{
  const TA & a;
  const TB & b;
public:

  enum { IS_LINEAR = TA::IS_LINEAR && TB::IS_LINEAR };
    
  INLINE AVXSumExpr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }

  INLINE auto operator() (int i) const -> decltype(a(i)+b(i)) { return a(i)+b(i); }
  INLINE auto operator() (int i, int j) const -> decltype(a(i,j)+b(i,j)) { return a(i,j)+b(i,j); }
  __m256d Get(int i) const { return a.Get(i)+b.Get(i); } 

  INLINE int Height() const { return a.Height(); }
  INLINE int Width() const { return a.Width(); }

  void Dump (ostream & ost) const
  { ost << "("; a.Dump(ost); ost << ") + ("; b.Dump(ost); ost << ")"; }
};

template <typename TA, typename TB>
INLINE AVXSumExpr<TA, TB>
operator+ (const AVXExpr<TA> & a, const AVXExpr<TB> & b)
{
  return AVXSumExpr<TA, TB> (a.Spec(), b.Spec());
}



/* *************************** PW_Mult_Expr **************************** */

template <class TA, class TB> 
class AVXPW_Mult_Expr : public AVXExpr<AVXPW_Mult_Expr<TA,TB> >
{
  const TA & a;
  const TB & b;
public:

  enum { IS_LINEAR = TA::IS_LINEAR && TB::IS_LINEAR };
    
  INLINE AVXPW_Mult_Expr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }

  INLINE auto operator() (int i) const -> decltype(a(i)*b(i)) { return a(i)*b(i); }
  INLINE auto operator() (int i, int j) const -> decltype(a(i,j)*b(i,j)) { return a(i,j)*b(i,j); }
  __m256d Get(int i) const { return a.Get(i)*b.Get(i); } 

  INLINE int Height() const { return a.Height(); }
  INLINE int Width() const { return a.Width(); }

  void Dump (ostream & ost) const
  { ost << "("; a.Dump(ost); ost << ") + ("; b.Dump(ost); ost << ")"; }
};

template <typename TA, typename TB>
INLINE AVXPW_Mult_Expr<TA, TB>
pw_mult (const AVXExpr<TA> & a, const AVXExpr<TB> & b)
{
  return AVXPW_Mult_Expr<TA, TB> (a.Spec(), b.Spec());
}



/* *************************** ScaleExpr **************************** */

template <class TA> 
class AVXScaleExpr : public AVXExpr<AVXScaleExpr<TA> >
{
  double s1;
  __m256d s;
  const TA & a;
public:

  enum { IS_LINEAR = TA::IS_LINEAR };
    
  INLINE AVXScaleExpr (double as, const TA & aa) : s1(as), s(_mm256_set1_pd(as)), a(aa) { ; }

  INLINE auto operator() (int i) const -> decltype(s1*a(i)) { return s1*a(i); }
  INLINE auto operator() (int i, int j) const -> decltype(s1*a(i,j)) { return s1*a(i,j); }
  __m256d Get(int i) const { return s * a.Get(i); }

  INLINE int Height() const { return a.Height(); }
  INLINE int Width() const { return a.Width(); }

  void Dump (ostream & ost) const
  { ost << s << "*("; a.Dump(ost); ost << ")"; }
};

template <typename TA>
INLINE AVXScaleExpr<TA>
operator* (double s, const AVXExpr<TA> & a)
{
  return AVXScaleExpr<TA> (s, a.Spec());
}





class AFlatVectorD : public AVXExpr<AFlatVectorD>
{
  int size; // , vsize;
  unsigned vsize() const { return (unsigned(size)+3) / 4; }
  __m256d * __restrict data;
public:
  AFlatVectorD(int asize, LocalHeap & lh)
  {
    size = asize;
    // vsize = (unsigned(size)+3) / 4;
    data = (__m256d*)lh.Alloc<double> (4*vsize());
  }

  AFlatVectorD (int as, double * adata) : size(as), data((__m256d*)(void*)adata)
  {
    // vsize = (unsigned(size)+3) / 4;    
  }
  

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

  __m256d & Get(int i) const { return data[i]; }

  INLINE const AFlatVectorD & operator= (const AFlatVectorD & v2) const
  {
    for (int i = 0; i < vsize(); i++)
      data[i] = v2.data[i];
    return *this;
  }
  
  AFlatVectorD & operator= (double d)
  {
    for (int i = 0; i < vsize(); i++)
      data[i] = _mm256_set1_pd(d);
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
  INLINE const AFlatVectorD & operator= (const AVXExpr<TB> & v) const
  {
    for (int i = 0; i < vsize(); i++)
      data[i] = v.Spec().Get(i);
    return *this;
    /*
    int i = 0;
    for (; i < vsize-1; i+=2)
      {
        data[i] = v.Spec().Get(i);
        data[i+1] = v.Spec().Get(i+1);
      }
    if (i < vsize)
      {
        data[i] = v.Spec().Get(i);
        i++;
      }
    return *this;
    */
  }


  template<typename TB>
  INLINE const AFlatVectorD & operator+= (const AVXExpr<TB> & v) const
  {
    // for (int i = 0; i < vsize(); i++)
    for (auto i : Range(vsize()))
      data[i] += v.Spec().Get(i);
    return *this;
  }
  
};



class AFlatMatrixD : public AVXExpr<AFlatMatrixD>
{
  int h, w;
  // static int vw;
  __m256d * __restrict data;
public:
  AFlatMatrixD () = default;
  AFlatMatrixD (const AFlatMatrixD &) = default;
  AFlatMatrixD(int ah, int aw, LocalHeap & lh)
  {
    h = ah;
    w = aw;
    // vw = (w+3)/4;
    int rw = (w+3)&(-4);
    data = (__m256d*)lh.Alloc<double> (h*rw);
  }

  AFlatMatrixD(int ah, int aw, double * adata)
    : h(ah), w(aw), data((__m256d*)(void*)adata) { ; } 

  AFlatMatrixD(int ah, int aw, __m256d * mem)
    : h(ah), w(aw), data(mem) { ; } 

  void AssignMemory (int ah, int aw, __m256d * mem)
  {
    h = ah;
    w = aw;
    data = mem;
  }
  
  enum { IS_LINEAR = false };
  int Size () const { return h*w; }
  int Height () const { return h; }
  int Width () const { return w; }
  unsigned int VWidth() const { return (unsigned(w)+3)/4; }
  double & operator() (int i) const
  {
    return ((double*)data)[i]; 
  }

  double & operator() (int i, int j) const
  {
    int vw = VWidth(); // (w+3)/4;
    return ((double*)data)[4*i*vw+j]; 
  }

  __m256d & Get(int i) const { return data[i]; }
  __m256d & Get(int i, int j) const { return data[i*VWidth()+j]; }

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
      data[i] = _mm256_set1_pd(d);
      // _mm256_store_pd((double*)&data[i], _mm256_set1_pd(d));
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
      data[i] *= _mm256_set1_pd(d);
    return *this;
  }
  
  operator SliceMatrix<double> () const
  {
    int vw = VWidth(); // (w+3)/4;
    return SliceMatrix<double> (h, w, 4*vw, (double*)data);
  }
   
  SliceVector<> Col (int c) const
  {
    int vw = VWidth(); // (w+3)/4;    
    return SliceVector<> (h, 4*vw, ((double*)data)+c);
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
    return SliceMatrix<>(h, end-begin, 4*vw, ((double*)data)+begin);
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
  return SliceMatrix<double,ColMajor> (mat.Width(), mat.Height(), 4*mat.VWidth(), &mat(0,0));
}



INLINE __m256d HAdd (__m256d v1, __m256d v2, __m256d v3, __m256d v4)
{
  __m256d hsum1 = _mm256_hadd_pd (v1, v2);
  __m256d hsum2 = _mm256_hadd_pd (v3, v4);
  
  __m256d hsum = 
    _mm256_add_pd (_mm256_insertf128_pd (hsum1, 
                                         _mm256_extractf128_pd (hsum2, 0), 1),
                   _mm256_insertf128_pd (hsum2, 
                                         _mm256_extractf128_pd (hsum1, 1), 0));
  return hsum;
}




extern void AddABt (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c);
extern void AddABt (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);
extern void AddABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);

extern void AddABtSym (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c);
extern void AddABtSym (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);
extern void AddABtSym (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);




// c = a * Diag (d)
static 
void MultMatDiagMat(AFlatMatrixD a, AFlatVectorD diag, AFlatMatrixD c)
{
  /*
  for (int i = 0; i < diag.Size(); i++)
    c.Col(i) = diag(i) * a.Col(i);
  */
  int rest = 4*diag.VSize() - diag.Size();
  int loops = diag.VSize();
  if (rest) loops--;
  
  for (int i = 0; i < c.Height(); i++)
    for (int j = 0; j < loops; j++)
      c.Get(i,j) = a.Get(i,j) * diag.Get(j);

  if (rest)
    {
      __m256i mask = _mm256_cmpgt_epi64(_mm256_set1_epi64x(4-rest),
                                        _mm256_set_epi64x(3, 2, 1, 0));

      __m256d md = _mm256_maskload_pd((double*)&diag.Get(loops), mask);

      for (int i = 0; i < c.Height(); i++)
        {
          __m256d ma = _mm256_maskload_pd((double*)&a.Get(i,loops), mask);
          __m256d prod = md * ma;
          _mm256_maskstore_pd((double*)&c.Get(i,loops), mask, prod);
        }
    }
}



template <>
class ABareVector<double>
{
  __m256d * __restrict data;
public:
  ABareVector(__m256d * _data) : data(_data) { ; }
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
  __m256d & Get(int i) const { return data[i]; }
};

template <>
class ABareMatrix<double>
{
  __m256d * __restrict data;
  int dist;   // dist in simds
public:
  ABareMatrix(__m256d * _data, int _dist) : data(_data), dist(_dist) { ; }
  ABareMatrix(AFlatMatrix<double> mat) : data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
  ABareMatrix(const ABareMatrix &) = default;

  double & operator() (int i, int j) const
  {
    return ((double*)data)[4*i*dist+j]; 
  }
  __m256d & Get(int i, int j) const { return data[i*dist+j]; }
  ABareVector<double> Row(int i) const { return ABareVector<double> (data+i*dist); }
  ABareMatrix<double> Rows(int first, int /* next */) const { return ABareMatrix<double> (data+first*dist, dist); }
  ABareMatrix<double> Rows(IntRange r) const { return Rows(r.First(), r.Next()); } 
  ABareMatrix<double> RowSlice(int first, int adist) const { return ABareMatrix<double> (data+first*dist, dist*adist); } 
};




#else // __AVX2__

// template <typename T>
// using AFlatMatrix = FlatMatrix<T>;

class AFlatVectorD : public FlatVector<double>
{
public:
  AFlatVectorD (int as, LocalHeap & lh)
    : FlatVector<double> (as, lh) { ; }
  AFlatVectorD (int as, double * p)
    : FlatVector<double> (as, p) { ; }
  
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
  ABareMatrix(AFlatMatrix<double> mat) : data(&mat.Get(0,0)), dist(&mat.Get(1,0)-&mat.Get(0,0)) { ; }
  ABareMatrix(const ABareMatrix &) = default;

  double & operator() (int i, int j) const
  {
    return ((double*)data)[i*dist+j]; 
  }
  double & Get(int i, int j) const { return data[i*dist+j]; }
  ABareVector<double> Row(int i) const { return ABareVector<double> (data+i*dist); }
  ABareMatrix<double> Rows(int first, int /* next */) const { return ABareMatrix<double> (data+first*dist, dist); }
  ABareMatrix<double> Rows(IntRange r) const { return Rows(r.First(), r.Next()); } 
  ABareMatrix<double> RowSlice(int first, int adist) const { return ABareMatrix<double> (data+first*dist, dist*adist); } 
};


template <typename TA, typename TB, typename TC>
INLINE void AddABt (const TA & a, const TB & b, SliceMatrix<TC> c)
{
  c += a * Trans(b) | Lapack;
  // LapackMultAdd (a, Trans(b), 1.0, c, 1.0);
}

/*
template <typename T>
INLINE void AddABtSym (AFlatMatrix<T> a, AFlatMatrix<T> b, SliceMatrix<T> c)
*/

template <typename TA, typename TB, typename TC>
INLINE void AddABtSym (const TA & a, const TB & b, SliceMatrix<TC> c)
{
  c += a * Trans(b) | Lapack;
  // LapackMultAdd (a, Trans(b), 1.0, c, 1.0);
}


#endif // __AVX2__

#endif // FILE_AVECTOR
