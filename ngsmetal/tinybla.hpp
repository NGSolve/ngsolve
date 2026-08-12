string code_tinybla = R"(

#ifdef __CUDACC__
#define TB_HD __host__ __device__
#define thread
#else
#define TB_HD
#endif


namespace tinybla {

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

  template <ORDERING ORD, typename T = float>
  class BareMatrixShared
  {
     threadgroup T * data;
     int ld;
  public:
     BareMatrixShared (threadgroup T * _data, int _ld) : data(_data), ld(_ld) { }
     T operator() (int r, int c) const
     {
        if constexpr (ORD==RowMajor) return data[r*ld+c];
        else  return data[c*ld+r];
     }

     int Offset (int r, int c) const {
       if constexpr (ORD==RowMajor) return r*ld+c;
       else return c*ld+r;
     }

     template <int h, int w>
     auto GetTile(int r, int c) const {
       Mat<h,w,float> res;
       for (int i = 0; i < h; i++)
         for (int j = 0; j < w; j++)
           res(i,j) = (*this)(r+i,c+j);
       return res;
     }

     auto SubMatrix (int r, int c) const { return BareMatrixShared<ORD,T> { data+Offset(r,c), ld }; }
     auto Transpose() const { return BareMatrixShared<!ORD, T> { data, ld }; }
     int LD() const { return ld; }
     threadgroup T* Data() const { return data; }
     bool IsColMajor() const { return ORD==ColMajor; }
  };

  template <ORDERING ORD, typename T = float>
  class BareMatrixDevice
  {
     device T * data;
     int ld;
  public:
     BareMatrixDevice (device T * _data, int _ld) : data(_data), ld(_ld) { }
     T operator() (int r, int c) const
     {
        if constexpr (ORD==RowMajor) return data[r*ld+c];
        else  return data[c*ld+r];
     }

     int Offset (int r, int c) const {
       if constexpr (ORD==RowMajor) return r*ld+c;
       else return c*ld+r;
     }

     template <int h, int w>
     auto GetTile(int r, int c) const {
       Mat<h,w,float> res;
       for (int i = 0; i < h; i++)
         for (int j = 0; j < w; j++)
           res(i,j) = (*this)(r+i,c+j);
       return res;
     }

     auto SubMatrix (int r, int c) const { return BareMatrixDevice<ORD,T> { data+Offset(r,c), ld }; }
     auto Transpose() const { return BareMatrixShared<!ORD, T> { data, ld }; }
     int LD() const { return ld; }
     device T* Data() const { return data; }
     bool IsColMajor() const { return ORD==ColMajor; }
};



  template <int H, int W, typename T = float>
  class WarpMatrix
  {
    static_assert(H==8);
    static_assert(W==8);
    static constant constexpr int BW = 2;
    static constant constexpr int BH = 1;
    Mat<BH,BW,T> myvals; 
  public:
    WarpMatrix() { }
    WarpMatrix(T val) { myvals = val; }
    WarpMatrix(threadgroup T * data, int ld, int tid) {
      int r = MyRow(tid);
      int c = MyCol(tid);
      for (int i = 0; i < BW; i++)
         myvals(0,i) = data + r*ld + c + i;
    }

    void operator= (T val) { myvals = val; }

    template <int K, typename M1, typename M2>
    void AddMM(M1 m1, M2 m2, int tid)
    {
      int r = MyRow(tid);
      int c = MyCol(tid);
      for (int k = 0; k < K; k+=2)
        myvals += m1.template GetTile<BH,2>(r,k) * m2.template GetTile<2,BW>(k,c);          
    }

    void Store(threadgroup T * data, int ld, int tid) {
      int r = MyRow(tid);
      int c = MyCol(tid);
      for (int i = 0; i < BW; i++)
         data[r*ld + c + i] = myvals(0,i);
    }

    int MyCol(int tid) const { return 2*(tid%4); }
    int MyRow(int tid) const { return (tid/4)%8; }
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

    void operator= (T val) { m = metal::make_filled_simdgroup_matrix<float, 8, 8>(0.0f); }

    template <int K, typename M1, typename M2>
    void AddMM(M1 m1, M2 m2, int tid)
    {
      metal::simdgroup_float8x8 ma, mb;
      metal::simdgroup_load(ma, m1.Data(), m1.LD(), ulong2(0, 0), m1.IsColMajor());
      metal::simdgroup_load(mb, m2.Data(), m2.LD(), ulong2(0, 0), m2.IsColMajor());
      metal::simdgroup_multiply_accumulate(m, ma, mb, m);
    }

    void Store(threadgroup T * data, int ld, int tid) {
      metal::simdgroup_store(m, data, ld, ulong2(0,0));
    }
  };



} // namespace tinybla



)";
