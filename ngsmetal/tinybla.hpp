string code_tinybla = R"(

#include <metal_stdlib>
using namespace metal;

#ifdef __CUDACC__
#define TB_HD __host__ __device__
#define thread
#else
#define TB_HD
#endif


namespace tinybla {






  template <int S, typename T>
  class HTVec {
    HTVec<S-1,T> tail;
    T head;

  public:
    HTVec () { }
    HTVec (HTVec<S-1,T> _tail, T _head) : tail(_tail), head(_head) {}

    auto Head() const { return head; } 
    auto Tail() const { return tail; } 

    void operator= (T val) {
      tail = val;
      head = val;
    }

    void operator+= (HTVec b) {
      tail += b.Tail();
      head += b.Head();
    }

    template <typename TP>
    void Load(TP data) {
      tail.Load(data);
      head = data[S-1];
    }
    
    template <typename TP>
    void Store(TP data) {
      tail.Store(data);
      data[S-1] = head;
    }
  };



  template <typename T>
  class HTVec<1,T> {
    T head;

  public:
    auto Head() const { return head; } 
    
    HTVec () { }
    HTVec (T _head) : head(_head) {}

    void operator= (T val) {
      head = val;
    }

    void operator+= (HTVec b) {
      head += b.Head();
    }

    template <typename TP>
    void Load(TP data) { head = data[0]; }
    
    template <typename TP>
    void Store(TP data) { data[0] = head; }
  };




  template <int S, typename T>
  inline auto operator+ (HTVec<S,T> a, HTVec<S,T> b) -> HTVec<S,T> {
    if constexpr (S==1)
      return a.Head()+b.Head();
    else
      return { a.Tail()+b.Tail(), a.Head()+b.Head() };
  }

  template <int S, typename T>
  inline auto operator* (float a, HTVec<S,T> b) -> HTVec<S,T> {
    if constexpr (S==1)
      return a*b.Head();
    else
      return { a*b.Tail(), a*b.Head() };
  }
  
  




  template <int H, int W, typename T>
  class HTMat {
    HTMat<H-1,W,T> tail;
    HTVec<W,T> head;

  public:
    HTMat () { }
    HTMat (HTMat<H-1,W,T> _tail, HTVec<W,T> _head) : tail(_tail), head(_head) {}
    template <typename TP>
    HTMat (TP data, unsigned dist) { Load(data, dist); }

    auto Head() const { return head; } 
    auto Tail() const { return tail; } 

    void operator= (T val) {
      tail = val;
      head = val;
    }

    template <typename TP>
    void Load(TP data, unsigned dist) {
      tail.Load(data,dist);
      head.Load(data+(H-1)*dist);
    }
    
    template <typename TP>
    void Store(TP data, unsigned dist) {
      tail.Store(data,dist);
      head.Store(data+(H-1)*dist);
    }

    void operator+= (HTMat b) {
      tail += b.Tail();
      head += b.Head();
    }
  };

  
  template <int W, typename T>
  class HTMat<1,W,T> {
    HTVec<W,T> head;

  public:
    auto Head() const { return head; } 
    
    HTMat () { }
    HTMat (HTVec<W,T> _head) : head(_head) {}
    template <typename TP>
    HTMat (TP data, unsigned dist) { Load(data, dist); }

    void operator= (T val) {
      head = val;
    }


    template <typename TP>
    void Load(TP data, unsigned dist) { head.Load(data); }
    
    template <typename TP>
    void Store(TP data, unsigned dist) { head.Store(data); }


    void operator+= (HTMat b) {
      head += b.Head();
    }
  };

  
  template <int H, int W, typename T>
  auto operator+ (HTMat<H,W,T> a, HTMat<H,W,T> b) -> HTMat<H,W,T> {
    if constexpr (H==1)
      return a.Head()+b.Head();
    else
      return { a.Tail()+b.Tail(), a.Head()+b.Head() };
  }

  // alpha * x + y
  template <int S, typename TA, typename T>
  auto FMA (TA alpha, HTVec<S,T> x, HTVec<S,T> y) -> HTVec<S,T>
  {
    if constexpr (S==1)
      return fma(static_cast<T>(alpha), x.Head(), y.Head());
    else
      return { FMA(alpha, x.Tail(), y.Tail()), fma(static_cast<T>(alpha), x.Head(), y.Head()) };
  } 


  template <int S, int W, typename TA, typename T>
  inline auto operator* (HTVec<S,TA> a, HTMat<S,W,T> b) -> HTVec<W,T> {
    if constexpr (S==1)
      return a.Head()*b.Head();
    else
      // return a.Tail()*b.Tail() + a.Head()*b.Head();
      return FMA (a.Head(),b.Head(), a.Tail()*b.Tail());
  }
  
  template <int H, int K, int W, typename TA, typename T>
  inline auto operator* (HTMat<H,K,TA> a, HTMat<K,W,T> b) -> HTMat<H,W,T> {
    if constexpr (H==1)
      return a.Head()*b;
    else
      return { a.Tail()*b, a.Head()*b };
  }







  template <int H, int K, int W, typename TA, typename T>
  inline auto FMA (HTMat<H,K,TA> a, HTMat<K,W,T> b, HTMat<H,W,T> c) -> HTMat<H,W,T> {
    // return c + a*b;
    if constexpr (K==1)
      return AddOuter(LastCol(a), b.Head(), c);
    else
      return AddOuter(LastCol(a), b.Head(), FMA(RemoveCol(a), b.Tail(), c));
  }

  


  template <int S, typename T>
  auto ColVec (HTVec<S,T> v) -> HTMat<S,1,T>
  {
    HTVec<1,T> last;
    last = v.Head();
    if constexpr (S==1)
      return last;
    else
      return { ColVec(v.Tail()), last };
  }

  template <int H, int W, typename T>
  auto AddCol(HTMat<H,W,T> m, HTVec<H,T> v) -> HTMat<H,W+1,T>
  {
    if constexpr (H==1)
      return HTVec<W+1,T> (m.Head(), v.Head());
    else
      return { AddCol(m.Tail(), v.Tail()), HTVec<W+1,T> (m.Head(), v.Head()) };
  }

  template <int H, int W, typename T>
  auto Trans (HTMat<H,W,T> a) -> HTMat<W,H,T> {
    if constexpr (H==1)
      return ColVec (a.Head()); 
    else
      return AddCol (Trans(a.Tail()), a.Head() );
  }

  template <int H, int W, typename T>
  auto RemoveCol(HTMat<H,W,T> m) -> HTMat<H,W-1,T>
  {
    if constexpr (H==1)
      return m.Head().Tail();
    else
      return { RemoveCol(m.Tail()), m.Head().Tail() };
  }

  template <int H, int W, typename T>
  auto LastCol(HTMat<H,W,T> m) -> HTVec<H,T>
  {
    if constexpr (H==1)
      return m.Head().Head();
    else
      return { LastCol(m.Tail()), m.Head().Head() };
  }


  template <int H, int W, typename T>
  auto Outer (HTVec<H,T> a, HTVec<W,T> b) -> HTMat<H,W,T>
  {
    if constexpr (H==1)
      return a.Head()*b;
    else
      return { Outer(a.Tail(), b), a.Head()*b };
  }

  template <int H, int W, typename TA, typename T>
  auto AddOuter (HTVec<H,TA> a, HTVec<W,T> b, HTMat<H,W,T> m) -> HTMat<H,W,T>
  {
    if constexpr (H==1)
      // return a.Head()*b+m.Head();
      return FMA(a.Head(), b, m.Head());
    else
      // return { AddOuter(a.Tail(), b, m.Tail()), a.Head()*b+m.Head() };
      return { AddOuter(a.Tail(), b, m.Tail()), FMA(a.Head(),b,m.Head()) };
  }






  template <int S, typename T>
  class Vec {
    T data[S];
  public:
    Vec() = default;
    template <typename... Args>
    TB_HD Vec(Args... args) : data{static_cast<T>(args)...} {
      static_assert(sizeof...(args) == S, "wrong number of arguments");
    }
    TB_HD constexpr int Size() const { return S; }
    TB_HD thread T & operator()(int i) { return data[i]; }
    TB_HD T operator()(int i) const  { return data[i]; }
    TB_HD void operator= (T val) { for (int i = 0; i < S; i++) data[i] = val; }
    TB_HD Vec operator+(Vec b) const {
      Vec r;
      for (int i = 0; i < S; i++) r(i) = data[i] + b(i);
      return r;
    }
    TB_HD Vec operator-(Vec b) const {
      Vec r;
      for (int i = 0; i < S; i++) r(i) = data[i] - b(i);
      return r;
    }

    TB_HD Vec operator+=(Vec b) {
      for (int i = 0; i < S; i++) data[i] += b(i);
      return *this;
    }

    template <int FIRST, int NEXT>
    TB_HD constexpr Vec<NEXT - FIRST, T> Range() const {
    Vec<NEXT - FIRST, T> r;
    for (int i = 0; i < NEXT - FIRST; ++i)
      r(i) = data[FIRST + i];
    return r;
    }

    template <int FIRST, int NEXT>
    TB_HD constexpr void SetRange(const Vec<NEXT - FIRST, T> sub) {
    for (int i = 0; i < NEXT - FIRST; ++i)
      data[FIRST + i] = sub(i);
  }
  };

  template <int S, typename T>
  TB_HD Vec<S,T> operator*(T s, Vec<S,T> v) {
    Vec<S,T> r;
    for (int i = 0; i < S; i++) r(i) = s*v(i);
    return r;
  }


  template <int H, int W, typename T>
  class Mat {
    T data[H][W];
  public:
    Mat() = default;
    TB_HD Mat(Vec<H*W,T> vec) {
      for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
          data[i][j] = vec(i*W+j);
    }
    TB_HD constexpr int Height() const { return H; }
    TB_HD constexpr int Width()  const { return W; }
    TB_HD thread T & operator()(int i, int j) { return data[i][j]; }
    TB_HD T operator()(int i, int j) const  { return data[i][j]; }
    TB_HD void operator= (T val) {
      for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
          data[i][j] = val;
    }
    TB_HD operator Vec<H*W,T> () const {
      Vec<H*W,T> r;
      for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
          r(i*W+j) = data[i][j];
      return r;
    }

    TB_HD Mat operator+(Mat b) const {
      Mat r;
      for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
          r(i,j) = data[i][j] + b(i,j);
      return r;
    }

    TB_HD Mat operator+= (Mat b) {
      for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
          data[i][j] += b(i,j);
      return *this;
    }

    TB_HD Mat operator-(Mat b) const {
      Mat r;
      for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
          r(i,j) = data[i][j] - b(i,j);
      return r;
    }
    TB_HD Mat operator*(T s) const {
      Mat r;
      for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
          r(i,j) = data[i][j] * s;
      return r;
    }

    template <int K>
    TB_HD Mat<H,K,T> operator*(Mat<W,K,T> b) const {
      Mat<H,K,T> r;
      for (int i = 0; i < H; i++)
        for (int k = 0; k < K; k++) {
          r(i,k) = 0;
          for (int j = 0; j < W; j++)
            r(i,k) += data[i][j] * b(j,k);
        }
      return r;
    }

    TB_HD Vec<H,T> operator*(Vec<W,T> v) const {
      Vec<H,T> r;
      for (int i = 0; i < H; i++) {
        r(i) = 0;
        for (int j = 0; j < W; j++)
          r(i) += data[i][j] * v(j);
      }
      return r;
    }

  };

  template <int H, int W, typename T>
  TB_HD Mat<H,W,T> operator*(T s, Mat<H,W,T> m) { return m * s; }


  template <int H, int W, typename T>
  TB_HD Mat<W,H,T> Trans(Mat<H,W,T> mat) {
    Mat<W,H,T> r;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        r(j,i) = mat(i,j);
    return r;
  }


  // Cof — square only, S<=3

  template <typename T>
  TB_HD Mat<1,1,T> Cof(Mat<1,1,T> m) {
    Mat<1,1,T> r;
    r(0,0) = T(1);
    return r;
  }

  template <typename T>
  TB_HD Mat<2,2,T> Cof(Mat<2,2,T> m) {
    Mat<2,2,T> r;
    r(0,0) =  m(1,1); r(0,1) = -m(1,0);
    r(1,0) = -m(0,1); r(1,1) =  m(0,0);
    return r;
  }

  template <typename T>
  TB_HD Mat<3,3,T> Cof(Mat<3,3,T> m) {
    Mat<3,3,T> r;
    r(0,0) =  m(1,1)*m(2,2) - m(1,2)*m(2,1);
    r(0,1) = -(m(1,0)*m(2,2) - m(1,2)*m(2,0));
    r(0,2) =  m(1,0)*m(2,1) - m(1,1)*m(2,0);
    r(1,0) = -(m(0,1)*m(2,2) - m(0,2)*m(2,1));
    r(1,1) =  m(0,0)*m(2,2) - m(0,2)*m(2,0);
    r(1,2) = -(m(0,0)*m(2,1) - m(0,1)*m(2,0));
    r(2,0) =  m(0,1)*m(1,2) - m(0,2)*m(1,1);
    r(2,1) = -(m(0,0)*m(1,2) - m(0,2)*m(1,0));
    r(2,2) =  m(0,0)*m(1,1) - m(0,1)*m(1,0);
    return r;
  }


  // Det

  template <typename T>
  TB_HD T Det(Mat<1,1,T> m) { return m(0,0); }

  template <typename T>
  TB_HD T Det(Mat<2,2,T> m) {
    return m(0,0)*m(1,1) - m(0,1)*m(1,0);
  }

  template <typename T>
  TB_HD T Det(Mat<3,3,T> m) {
    return m(0,0)*(m(1,1)*m(2,2) - m(1,2)*m(2,1))
      -m(0,1)*(m(1,0)*m(2,2) - m(1,2)*m(2,0))
      +m(0,2)*(m(1,0)*m(2,1) - m(1,1)*m(2,0));
  }


  // Inv = Cof^T / Det

  template <int S, typename T>
  TB_HD Mat<S,S,T> Inv(Mat<S,S,T> m) {
    return Trans(Cof(m)) * (T(1) / Det(m));
  }


  template <int H, int W, typename T>
  TB_HD auto ToMat(Vec<H*W,T> vec) { return Mat<H,W,T>(vec); }

  enum ORDERING { ColMajor, RowMajor };
  constexpr ORDERING operator! (ORDERING o) { return (o==RowMajor) ? ColMajor : RowMajor; }


  template <ORDERING ORD, typename Tp>
  class BareMatrix
  {
     Tp data;
     int ld;
     using ElementType = remove_addrspace_t<remove_cv_t<remove_reference_t<decltype(*data)>>>;

  public:
     BareMatrix (Tp _data, int _ld) : data(_data), ld(_ld) { }

     int Offset (int r, int c) const {
       if constexpr (ORD==RowMajor) return r*ld+c;
       else return c*ld+r;
     }

     auto operator() (int r, int c) const { return data[Offset(r,c)]; }

     template <int h, int w>
     auto GetTile(int r, int c) const {
      if constexpr (ORD==RowMajor)
        return HTMat<h,w,ElementType> (data+Offset(r,c), ld);
      else
        return Trans(HTMat<w,h,ElementType> (data+Offset(r,c), ld));

/*
       Mat<h,w,float> res;
       for (int i = 0; i < h; i++)
         for (int j = 0; j < w; j++)
           res(i,j) = (*this)(r+i,c+j);
       return res;
*/
     }

     template <int h, int w, typename T2>
     void SetTile(int r, int c, HTMat<h,w,T2> tile) {
      if constexpr (ORD==RowMajor)
        tile.Store(data+Offset(r,c), ld);
      else
        Trans(tile).Store(data+Offset(r,c), ld);

/*
       for (int i = 0; i < h; i++)
         for (int j = 0; j < w; j++)
           data[Offset(r+i,c+j)] = tile(i,j);
*/
     }


     auto SubMatrix (int r, int c) const { return BareMatrix<ORD,Tp> { data+Offset(r,c), ld }; }
     auto Transpose() const { return BareMatrix<!ORD, Tp> { data, ld }; }
     int LD() const { return ld; }
     Tp Data() const { return data; }
     bool IsColMajor() const { return ORD==ColMajor; }
  };

  template <ORDERING ORD, typename Tp>
  inline auto MakeBareMatrix (Tp data, int ld) {
    return BareMatrix<ORD,Tp> (data, ld);
  }


  template <int H, int W, typename T = float>
  class WarpMatrix
  {
    static_assert(H==8);
    static_assert(W==8);
    static constant constexpr int BW = sizeof(float2)/sizeof(T);
    static constant constexpr int BH = 1;
    HTMat<BH,BW,T> myvals; 
  public:
    WarpMatrix() { }

    WarpMatrix(T val) { myvals = val; }

    template <ORDERING ORD, typename Tp>
    WarpMatrix(BareMatrix<ORD,Tp> mat, unsigned tid) {
      unsigned r = MyRow(tid);
      unsigned c = MyCol(tid);
      myvals = mat.template GetTile<BH,BW>(r,c);
    }

    void operator= (T val) { myvals = val; }

    template <int K, ORDERING ORD1, typename Tp1, ORDERING ORD2, typename Tp2>
    void AddMM(BareMatrix<ORD1,Tp1> m1, BareMatrix<ORD2,Tp2> m2, unsigned tid)
    {
      unsigned r = MyRow(tid);
      unsigned c = MyCol(tid);
      // myvals += m1.template GetTile<BH,K>(r,0) * m2.template GetTile<K,BW>(0,c);
      myvals = FMA (m1.template GetTile<BH,K>(r,0), m2.template GetTile<K,BW>(0,c), myvals);
    }

    template <ORDERING ORD, typename Tp>
    void Store(BareMatrix<ORD,Tp> mat, unsigned tid) {
      unsigned r = MyRow(tid);
      unsigned c = MyCol(tid);
      mat.template SetTile<BH,BW>(r,c, myvals);
    }

    unsigned MyCol(unsigned tid) const { return BW*(tid%4); }
    unsigned MyRow(unsigned tid) const { return (tid/4)%8; }
  };

/*
  template <>
  class WarpMatrix<8,8,float>
  {
    typedef float T;
    metal::simdgroup_float8x8 m;
  public:
    WarpMatrix() { }
    WarpMatrix(T val) { m = metal::make_filled_simdgroup_matrix<float, 8, 8>(0.0f); }
    WarpMatrix(threadgroup T * data, int ld, int tid) {
      metal::simdgroup_load(m, data, ld, ulong2(0,0));
    }
    template <typename M>
    WarpMatrix(M mat, int tid) {
      metal::simdgroup_load(m, mat.Data(), mat.LD(), ulong2(0,0), mat.IsColMajor());
    }

    void operator= (T val) { m = metal::make_filled_simdgroup_matrix<float, 8, 8>(0.0f); }

    template <int K, typename M1, typename M2>
    void AddMM(M1 m1, M2 m2, int tid)
    {
      static_assert(K==8);
      metal::simdgroup_float8x8 ma, mb;
      metal::simdgroup_load(ma, m1.Data(), m1.LD(), ulong2(0, 0), m1.IsColMajor());
      metal::simdgroup_load(mb, m2.Data(), m2.LD(), ulong2(0, 0), m2.IsColMajor());
      metal::simdgroup_multiply_accumulate(m, ma, mb, m);
    }

    void Store(threadgroup T * data, int ld, int tid) {
      metal::simdgroup_store(m, data, ld, ulong2(0,0));
    }

    template <typename M>
    void Store(M mat, int tid) {
      metal::simdgroup_store(m, mat.Data(), mat.LD(), ulong2(0,0), mat.IsColMajor());
    }
  };
*/

} // namespace tinybla



)";
