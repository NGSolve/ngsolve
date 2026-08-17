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
    void SetHead (T val) { head = val; }
    template <typename TP>
    void TailLoad(TP ptr) { }

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


  // broadcast inside quad
  template <int S, typename T>
  auto QuadBroadcast (HTVec<S,T> m, ushort lane) -> HTVec<S,T>
  {
    if constexpr (S==1)
      return quad_broadcast(m.Head(), lane);
    else
      return { QuadBroadcast(m.Tail(), lane), quad_broadcast(m.Head(), lane) };
  }

  template <int H, int W, typename T>
  auto QuadBroadcast (HTMat<H,W,T> m, ushort lane) -> HTMat<H,W,T>
  {
    if constexpr (H==1)
      return QuadBroadcast(m.Head(), lane);
    else
      return { QuadBroadcast(m.Tail(), lane), QuadBroadcast(m.Head(), lane) };
  }


  // broadcast across quad WARNING: dysnamic shuffle is SLOW
  template <int S, typename T>
  auto BroadcastQuad (HTVec<S,T> m, unsigned qlane, unsigned from) -> HTVec<S,T>
  {
    if constexpr (S==1)
      return simd_shuffle(m.Head(), 4*from+qlane);
    else
      return { BroadcastQuad(m.Tail(), qlane, from), simd_shuffle(m.Head(), 4*from+qlane) };
  }

  template <int H, int W, typename T>
  auto BroadcastQuad (HTMat<H,W,T> m, unsigned qlane, unsigned from) -> HTMat<H,W,T>
  {
    if constexpr (H==1)
      return BroadcastQuad(m.Head(), qlane, from);
    else
      return { BroadcastQuad(m.Tail(), qlane, from), BroadcastQuad(m.Head(), qlane, from) };
  }


  // broadcast across quad WARNING: dysnamic shuffle is SLOW
  template <ushort mask, int S, typename T>
  auto ShuffleXor (HTVec<S,T> m) -> HTVec<S,T>
  {
    if constexpr (S==1)
      return simd_shuffle_xor(m.Head(), mask);
    else
      return { ShuffleXor<mask>(m.Tail()), simd_shuffle_xor(m.Head(), mask) };
  }

  template <ushort mask, int H, int W, typename T>
  auto ShuffleXor (HTMat<H,W,T> m) -> HTMat<H,W,T>
  {
    if constexpr (H==1)
      return ShuffleXor<mask> (m.Head());
    else
      return { ShuffleXor<mask>(m.Tail()), ShuffleXor<mask>(m.Head()) };
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
     uint ld;
     using ElementType = remove_addrspace_t<remove_cv_t<remove_reference_t<decltype(*data)>>>;

  public:
     BareMatrix (Tp _data, unsigned _ld) : data(_data), ld(_ld) { }

     unsigned Offset (unsigned r, unsigned c) const {
       if constexpr (ORD==RowMajor) return r*ld+c;
       else return c*ld+r;
     }

     auto operator() (unsigned r, unsigned c) const { return data[Offset(r,c)]; }

     template<int N>
     BareMatrix ShiftCols() const {
     if constexpr (ORD==RowMajor) return BareMatrix(data+N, ld);
     else return BareMatrix(data+N*ld, ld);
     }

     template<int N>
     BareMatrix ShiftRows() const {
     if constexpr (ORD==RowMajor) return BareMatrix(data+N*ld, ld);
     else return BareMatrix(data+N, ld);
     }

     template <typename itype> 
     BareMatrix ShiftCols(itype n) const {
     if constexpr (ORD==RowMajor) return BareMatrix(data+n, ld);
     else return BareMatrix(data+n*ld, ld);
     }

     template <typename itype> 
     BareMatrix ShiftRows(itype n) const {
     if constexpr (ORD==RowMajor) return BareMatrix(data+n*ld, ld);
     else return BareMatrix(data+n, ld);
     }


     template <unsigned h, unsigned w>
     auto GetTile() const {
      if constexpr (ORD==RowMajor)
        return HTMat<h,w,ElementType> (data, ld);
      else
        return Trans(HTMat<w,h,ElementType> (data, ld));
    }

     template <unsigned h, unsigned w>
     auto GetTile(unsigned r, unsigned c) const {
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
     void SetTile(unsigned r, unsigned c, HTMat<h,w,T2> tile) {
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


     auto SubMatrix (unsigned r, unsigned c) const { return BareMatrix<ORD,Tp> { data+Offset(r,c), ld }; }
     auto Transpose() const { return BareMatrix<!ORD, Tp> { data, ld }; }
     unsigned LD() const { return ld; }
     Tp Data() const { return data; }
     bool IsColMajor() const { return ORD==ColMajor; }
  };

  template <ORDERING ORD, typename Tp>
  inline auto MakeBareMatrix (Tp data, int ld) {
    return BareMatrix<ORD,Tp> (data, ld);
  }


  template <unsigned H, unsigned W, typename T = float>
  class WarpMatrix
  {
    static_assert(H%8==0);
    static_assert(W%4==0);
    static constant constexpr unsigned BW = W/4; // sizeof(float2)/sizeof(T);
    static constant constexpr unsigned BH = H/8;
    HTMat<BH,BW,T> myvals; 
  public:
    WarpMatrix() { }

    WarpMatrix(T val) { myvals = val; }

    template <ORDERING ORD, typename Tp>
    WarpMatrix(BareMatrix<ORD,Tp> mat, uint tid) {
      myvals = mat.template GetTile<BH,BW>(MyRow(tid),MyCol(tid));
    }

    void operator= (T val) { myvals = val; }

    template <uint K, ORDERING ORD1, typename Tp1, ORDERING ORD2, typename Tp2>
    void AddMM(BareMatrix<ORD1,Tp1> m1, BareMatrix<ORD2,Tp2> m2, uint tid)
    {
      auto r = MyRow(tid);
      auto c = MyCol(tid);
      // myvals += m1.template GetTile<BH,K>(r,0) * m2.template GetTile<K,BW>(0,c);
      // myvals = FMA (m1.template GetTile<BH,K>(r,0), m2.template GetTile<K,BW>(0,c), myvals);


/*
      constexpr int KTILE = 1;
      for (unsigned k = 0; k < K; k+=KTILE)
        myvals = FMA (m1.template GetTile<BH,KTILE>(r,k), m2.template GetTile<KTILE,BW>(k,c), myvals);
*/



      // best for btdtb
      // static_assert (K%4==0, "K must be a multiple of 4");
      constexpr int KTILE = 1;
      unsigned k = 0;
      for ( ; k+4 <= K; k+=4)
        {
          auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + (tid&3));
          for (unsigned k1 = 0; k1 < 4; k1++)
            myvals = FMA (QuadBroadcast(ATileQuad, k1), m2.template GetTile<KTILE,BW>(k+k1,c), myvals);
        }
      for ( ;  k < K; k++)
        myvals = FMA (m1.template GetTile<BH,KTILE>(r,k), m2.template GetTile<KTILE,BW>(k,c), myvals);


/*
      static_assert (K%4==0, "K must be a multiple of 4");
      m1 = m1.ShiftRows(r).ShiftCols(tid&3);
      m2 = m2.ShiftCols(c);
      constexpr int KTILE = 1;
      for (uint k=0 ; k+4 <= K; k+=4)
        {
          auto ATileQuad = m1.template GetTile<BH,KTILE>(0,0);
          m1 = m1.template ShiftCols<4>();
          for (uint k1 = 0; k1 < 4; k1++) {
            myvals = FMA (QuadBroadcast(ATileQuad, k1), m2.template GetTile<KTILE,BW>(0,0), myvals);
            m2 = m2.template ShiftRows<1>();
          }
        }
*/

/*
      static_assert (K%4==0, "K must be a multiple of 4");
      // good MM, but not yet BTDTB     
      // lane info inside simdgroup
      ushort lane  = tid & 31;
      ushort qlane = lane & 3;         // 0..3 in quad
      ushort blane = (lane>>2) & 1;
      constexpr ushort KTILE = 1;

      m1 = m1.ShiftRows(r).ShiftCols(qlane);
      m2 = m2.ShiftRows(blane).ShiftCols(c);

      for (uint k = 0; k < K; k+=4)
        {
          auto ATileQuad = m1.template GetTile<BH,KTILE>();
          m1 = m1.ShiftCols(4*KTILE);
          for (ushort k1 = 0; k1 < 4; k1+=2) {
            auto b0 = m2.template GetTile<1,BW> ();
            m2 = m2.ShiftRows(2);
            auto b1 = ShuffleXor<4> (b0);
            auto lo = blane ? b1 : b0;
            auto hi = blane ? b0 : b1;

            myvals = FMA (QuadBroadcast(ATileQuad, k1), lo, myvals); 
            myvals = FMA (QuadBroadcast(ATileQuad, k1+1), hi, myvals);

          }
        }
*/









/*
      // lane info inside simdgroup
      uint lane  = tid & 31;
       uint qlane = lane & 3;         // 0..3 in quad
      uint blane = (lane>>2) & 3;
      constexpr int KTILE = 1;
      for (unsigned k = 0; k < K; k+=4)
        {

          auto b0 = m2.template GetTile<1,BW> (k+blane, c);
          auto b1 = ShuffleXor<4> (b0);
          auto b2 = ShuffleXor<8> (b0);
          auto b3 = ShuffleXor<8> (b1);

          auto pick_b = [&](uint idx)  {
            auto lo = (idx & 1) ? b1 : b0;
            auto hi = (idx & 1) ? b3 : b2;
            return (idx & 2) ? hi : lo;
          };

          auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + qlane);

#pragma unroll
          for (int k1 = 0; k1 < 4; k1++)
            myvals = FMA (QuadBroadcast(ATileQuad, k1), pick_b(blane^k1), myvals);
        }
*/










/*
     // WIP: use Xor/Select
      constexpr int KTILE = 1;
      unsigned qlane = tid%4;
      unsigned qgroup = (tid/4)%8;

      for (unsigned k = 0; k < K; k+=8)
        {
          auto BTileQuad = m2.template GetTile<KTILE,BW>(k+qgroup,c);

          auto BShuff4 = ShuffleXor<4> (BTileQuad);
          auto Bi = (qgroup&4) ? BTileQuad : BShuff4;

          auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + tid%4);
          for (unsigned k1 = 0; k1 < 4; k1++)
           {
            myvals = FMA (QuadBroadcast(ATileQuad, k1), Bi, myvals);
           }
        }
*/


/*
      unsigned qlane = tid%4;
      unsigned qgroup = (tid/4)%8;

      constexpr unsigned KTILE = 1;
//          auto BTileQuad = m2.template GetTile<KTILE,BW>(qgroup,c);

      for (unsigned k = 0; k < K; k+=8)
        {
          auto BTileQuad = m2.template GetTile<KTILE,BW>(k+qgroup,c);

          auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k+qlane);
          // for (unsigned k1 = 0; k1 < 4; k1++)
          // myvals = FMA (QuadBroadcast(ATileQuad, k1), BroadcastQuad(BTileQuad, qlane, k1), myvals);
          myvals = FMA (QuadBroadcast(ATileQuad, 0), BroadcastQuad(BTileQuad, qlane, 0), myvals);
          myvals = FMA (QuadBroadcast(ATileQuad, 1), BroadcastQuad(BTileQuad, qlane, 1), myvals);
          myvals = FMA (QuadBroadcast(ATileQuad, 2), BroadcastQuad(BTileQuad, qlane, 2), myvals);
          myvals = FMA (QuadBroadcast(ATileQuad, 3), BroadcastQuad(BTileQuad, qlane, 3), myvals);
          auto ATileQuad2 = m1.template GetTile<BH,KTILE>(r,k+4+qlane);
          // for (unsigned k1 = 0; k1 < 4; k1++)
            // myvals = FMA (QuadBroadcast(ATileQuad2, k1), BroadcastQuad(BTileQuad, qlane, 4+k1), myvals);
          myvals = FMA (QuadBroadcast(ATileQuad2, 0), BroadcastQuad(BTileQuad, qlane, 4), myvals);
          myvals = FMA (QuadBroadcast(ATileQuad2, 1), BroadcastQuad(BTileQuad, qlane, 5), myvals);
          myvals = FMA (QuadBroadcast(ATileQuad2, 2), BroadcastQuad(BTileQuad, qlane, 6), myvals);
          myvals = FMA (QuadBroadcast(ATileQuad2, 3), BroadcastQuad(BTileQuad, qlane, 7), myvals);
        }
*/


/*
      // 7 TF float, float2, float4
      constexpr int KTILE = 1;
      for (unsigned k = 0; k < K; k+=4)
        {
          auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + tid%4);
          for (unsigned k1 = 0; k1 < 4; k1++)
            {
              auto myB = m2.template GetTile<KTILE,BW>(k+k1,c);
              for (unsigned runs = 0; runs < 10; runs++)
                myvals = FMA (QuadBroadcast(ATileQuad, k1), myB, myvals);
            } 
        }
*/


/*
// 6.8 TF
      constexpr int KTILE = 1;
      for (int k = 0; k < K; k+=KTILE)
{
        auto atile = m1.template GetTile<BH,KTILE>(r,k);
        auto btile = m2.template GetTile<KTILE,BW>(k,c);
        for (int r = 0; r < 10; r++)
          myvals = FMA (atile, btile, myvals);
}
*/



/*
6 TF
      constexpr int KTILE = 1;
      for (int k = 0; k < K; k+=KTILE)
{
        auto atile = m1.template GetTile<BH,KTILE>(r,k);
        for (int r = 0; r < 10; r++)
          myvals = FMA (atile, m2.template GetTile<KTILE,BW>(k,c), myvals);
}
*/

/*
//  4-5 TF
      constexpr int KTILE = 1;
      for (int k = 0; k < K; k+=KTILE)
{
        auto btile = m2.template GetTile<KTILE,BW>(k,c);
        for (int r = 0; r < 10; r++)
          myvals = FMA (m1.template GetTile<BH,KTILE>(r,k), btile, myvals);
}
*/



/*
      constexpr int KTILE = 1;
      auto curr_A = m1.template GetTile<BH,KTILE>(r, 0);
      auto curr_B = m2.template GetTile<KTILE,BW>(0, c);

      for (int k = 0; k < K; k+=KTILE) {
        auto next_A = m1.template GetTile<BH, KTILE>(r, k + 1);
        auto next_B = m2.template GetTile<KTILE, BW>(k + 1, c);

        myvals = FMA(curr_A, curr_B, myvals); 

        curr_A = next_A;
        curr_B = next_B;
      } 
*/      
    }

    template <ORDERING ORD, typename Tp>
    void Store(BareMatrix<ORD,Tp> mat, unsigned tid) {
      mat.template SetTile<BH,BW>(MyRow(tid),MyCol(tid), myvals);
    }

    auto MyCol(uint tid) const { return BW*(tid&3); }
    auto MyRow(uint tid) const { return BH*((tid>>2)&7); }
  };


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

} // namespace tinybla



)";
