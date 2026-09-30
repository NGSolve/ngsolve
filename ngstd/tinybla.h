/*********************************************************************/
/* File:   tinybla.h                                                 */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   3. Aug. 2026                                              */
/*********************************************************************/

#ifdef __CUDACC__
#define TB_CUDA
#define TB_HD __host__ __device__
#define thread
#define constant
#define threadgroup
#define TB_DEVICE

#elif defined(NGS_GPU_CPU)

#define TB_CPU
#define TB_HD
#define thread
#define constant
#define threadgroup
#define TB_DEVICE
// the lane shuffles of the host reference backend, defined in gpukernel.h
namespace ngs_cpu { unsigned LaneId (); template <unsigned W, typename T> T SubgroupShuffle (T x, unsigned src); }

#else
#include <metal_stdlib>
using namespace metal;
#define TB_METAL
#define TB_HD
#define TB_DEVICE device   // pointer members need an address space

#endif


namespace tinybla {

#ifndef TB_METAL
  // metal supplies these, elsewhere tinybla has to be self-contained
  // (nvrtc cannot include <bla.hpp>). Declared inside the namespace so a
  // global "using namespace std" cannot make them ambiguous.
  typedef unsigned int   uint;
  typedef unsigned short ushort;

  template <typename T> struct tb_remove_reference     { using type = T; };
  template <typename T> struct tb_remove_reference<T&> { using type = T; };
  template <typename T> using remove_reference_t = typename tb_remove_reference<T>::type;

  template <typename T> struct tb_remove_cv                   { using type = T; };
  template <typename T> struct tb_remove_cv<const T>          { using type = T; };
  template <typename T> struct tb_remove_cv<volatile T>       { using type = T; };
  template <typename T> struct tb_remove_cv<const volatile T> { using type = T; };
  template <typename T> using remove_cv_t = typename tb_remove_cv<T>::type;

  template <typename T> using remove_addrspace_t = T;
#endif


  template <uint N>
  struct ic {
      using value_type = uint;
      static constant constexpr uint value = N;
      constexpr operator uint() const { return N; }
      constexpr uint operator()() const { return N; }
  };



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

    template <int DIST, typename TP>
    auto LoadSlice(TP ptr) {
      tail.template LoadSlice<DIST>(ptr);
      head = ptr[DIST*(S-1)];
      return *this;
    }

    template <int DIST, typename TP>
    void StoreSlice(TP ptr) {
      tail.template StoreSlice<DIST>(ptr);
      ptr[DIST*(S-1)] = head;
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
    void Store(TP ptr) { ptr[0] = head; }

    template <int DIST, typename TP>
    auto LoadSlice(TP ptr) { head = ptr[0]; return *this; }

    template <int DIST, typename TP>
    auto StoreSlice(TP ptr) { ptr[0] = head; }
  };


  template <typename T>
  class HTVec<0,T> {
  public:
    HTVec () { }
    void operator= (T val) { }
    void operator+= (HTVec b) { }

    template <typename TP>
    void Load(TP data) { }

    template <typename TP>
    void Store(TP data) { }

    template <int DIST, typename TP>
    auto LoadSlice(TP ptr) { return *this; }

    template <int DIST, typename TP>
    void StoreSlice(TP ptr) { }
  };



  template <int S, typename T>
  inline auto operator+ (HTVec<S,T> a, HTVec<S,T> b) -> HTVec<S,T> {
    if constexpr (S==0)
      return {};
    else if constexpr (S==1)
      return a.Head()+b.Head();
    else
      return { a.Tail()+b.Tail(), a.Head()+b.Head() };
  }

  template <int S, typename T>
  inline auto operator* (float a, HTVec<S,T> b) -> HTVec<S,T> {
    if constexpr (S==0)
      return {};
    else if constexpr (S==1)
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
    template <typename TP, typename TD=uint>
    HTMat (TP data, TD dist) { Load(data, dist); }

    auto Head() const { return head; } 
    auto Tail() const { return tail; } 

    void operator= (T val) {
      tail = val;
      head = val;
    }

    template <uint r>
    HTVec<W,T> GetRow() {
    if constexpr (r == H-1)
      return head;
    else
      return tail.template GetRow<r>();
    }


    template <typename TP>
    void Load(TP data, unsigned dist) {
      tail.Load(data,dist);
      head.Load(data+(H-1)*dist);
    }
    
    template <typename TP, typename TD>
    void Store(TP data, TD dist) {
      tail.Store(data,dist);
      head.Store(data+(H-1)*dist);
    }

    template <int DIST, typename TP, typename TD=uint>
    void LoadSlice(TP data, TD dist) {
      tail.template LoadSlice<DIST>(data,dist);
      head.template LoadSlice<DIST>(data+(H-1)*dist);
    }

    template <int DIST, typename TP, typename TD=uint>
    void StoreSlice(TP data, unsigned dist) {
      tail.template StoreSlice<DIST>(data,dist);
      head.template StoreSlice<DIST>(data+(H-1)*dist);
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


    template <uint r>
    HTVec<W,T> GetRow() {
      static_assert(r==0);
      return head;
    }


    template <typename TP>
    void Load(TP data, unsigned dist) { head.Load(data); }
    
    template <typename TP>
    void Store(TP data, unsigned dist) { head.Store(data); }

    template <int DIST, typename TP>
    void LoadSlice(TP data, unsigned dist) {
      head.template LoadSlice<DIST>(data);
    }

    template <int DIST, typename TP>
    void StoreSlice(TP data, unsigned dist) {
      head.template StoreSlice<DIST>(data);
    }


    void operator+= (HTMat b) {
      head += b.Head();
    }
  };

  template <int W, typename T>
  class HTMat<0,W,T> {
  public:
    HTMat () { }
    template <typename TP, typename TD=uint>
    HTMat (TP data, TD dist) { }

    void operator= (T val) { }
    void operator+= (HTMat b) { }

    template <typename TP>
    void Load(TP data, unsigned dist) { }

    template <typename TP, typename TD>
    void Store(TP data, TD dist) { }

    template <int DIST, typename TP, typename TD=uint>
    void LoadSlice(TP data, TD dist) { }

    template <int DIST, typename TP, typename TD=uint>
    void StoreSlice(TP data, TD dist) { }
  };


  
  template <int H, int W, typename T>
  auto operator+ (HTMat<H,W,T> a, HTMat<H,W,T> b) -> HTMat<H,W,T> {
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
      return a.Head()+b.Head();
    else
      return { a.Tail()+b.Tail(), a.Head()+b.Head() };
  }

  // alpha * x + y
  template <int S, typename TA, typename T>
  auto FMA (TA alpha, HTVec<S,T> x, HTVec<S,T> y) -> HTVec<S,T>
  {
    if constexpr (S==0)
      return {};
    else if constexpr (S==1)
      return fma(alpha, x.Head(), y.Head());
    else
      return { FMA(alpha, x.Tail(), y.Tail()), fma(alpha, x.Head(), y.Head()) };
  } 


  template <int S, int W, typename TA, typename T>
  inline auto operator* (HTVec<S,TA> a, HTMat<S,W,T> b) -> HTVec<W,T> {
    if constexpr (S==0)
      return {};
    else if constexpr (S==1)
      return a.Head()*b.Head();
    else
      // return a.Tail()*b.Tail() + a.Head()*b.Head();
      return FMA (a.Head(),b.Head(), a.Tail()*b.Tail());
  }
  
  template <int H, int K, int W, typename TA, typename T>
  inline auto operator* (HTMat<H,K,TA> a, HTMat<K,W,T> b) -> HTMat<H,W,T> {
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
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
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
      return HTVec<W+1,T> (m.Head(), v.Head());
    else
      return { AddCol(m.Tail(), v.Tail()), HTVec<W+1,T> (m.Head(), v.Head()) };
  }

  template <int H, int W, typename T>
  auto Trans (HTMat<H,W,T> a) -> HTMat<W,H,T> {
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
      return ColVec (a.Head()); 
    else
      return AddCol (Trans(a.Tail()), a.Head() );
  }

  template <int H, int W, typename T>
  auto RemoveCol(HTMat<H,W,T> m) -> HTMat<H,W-1,T>
  {
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
      return m.Head().Tail();
    else
      return { RemoveCol(m.Tail()), m.Head().Tail() };
  }

  template <int H, int W, typename T>
  auto LastCol(HTMat<H,W,T> m) -> HTVec<H,T>
  {
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
      return m.Head().Head();
    else
      return { LastCol(m.Tail()), m.Head().Head() };
  }


  template <int H, int W, typename T>
  auto Outer (HTVec<H,T> a, HTVec<W,T> b) -> HTMat<H,W,T>
  {
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
      return a.Head()*b;
    else
      return { Outer(a.Tail(), b), a.Head()*b };
  }

  template <int H, int W, typename TA, typename T>
  auto AddOuter (HTVec<H,TA> a, HTVec<W,T> b, HTMat<H,W,T> m) -> HTMat<H,W,T>
  {
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
      // return a.Head()*b+m.Head();
      return FMA(a.Head(), b, m.Head());
    else
      // return { AddOuter(a.Tail(), b, m.Tail()), a.Head()*b+m.Head() };
      return { AddOuter(a.Tail(), b, m.Tail()), FMA(a.Head(),b,m.Head()) };
  }


#ifdef TB_CUDA
  // the metal simd/quad primitives on top of cuda warp shuffles
  template <typename T> __device__ T quad_broadcast (T x, unsigned lane)
  { return __shfl_sync (0xffffffff, x, lane, 4); }
  template <typename T> __device__ T simd_shuffle (T x, unsigned lane)
  { return __shfl_sync (0xffffffff, x, lane); }
  template <typename T> __device__ T simd_shuffle_xor (T x, unsigned mask)
  { return __shfl_xor_sync (0xffffffff, x, mask); }
#endif

#ifdef TB_CPU
  // the metal simd/quad primitives on the host reference backend
  template <typename T> T quad_broadcast (T x, unsigned lane)
  { return ngs_cpu::SubgroupShuffle<4> (x, lane); }
  template <typename T> T simd_shuffle (T x, unsigned lane)
  { return ngs_cpu::SubgroupShuffle<32> (x, lane); }
  template <typename T> T simd_shuffle_xor (T x, unsigned mask)
  { return ngs_cpu::SubgroupShuffle<32> (x, (ngs_cpu::LaneId() % 32) ^ mask); }
#endif


  // broadcast inside quad
  template <int S, typename T>
  auto QuadBroadcast (HTVec<S,T> m, ushort lane) -> HTVec<S,T>
  {
    if constexpr (S==0)
      return {};
    else if constexpr (S==1)
      return quad_broadcast(m.Head(), lane);
    else
      return { QuadBroadcast(m.Tail(), lane), quad_broadcast(m.Head(), lane) };
  }

  template <int H, int W, typename T>
  auto QuadBroadcast (HTMat<H,W,T> m, ushort lane) -> HTMat<H,W,T>
  {
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
      return QuadBroadcast(m.Head(), lane);
    else
      return { QuadBroadcast(m.Tail(), lane), QuadBroadcast(m.Head(), lane) };
  }


/*
  // shorter air-code, but no performance difference
  inline auto QuadBroadcast (HTMat<4,1,float> m, ushort lane) -> HTMat<4,1,float>
  {
    float4 f4 { m.Tail().Tail().Tail().Head().Head(),
                m.Tail().Tail().Head().Head(),
                m.Tail().Head().Head(),
                m.Head().Head() };

     auto res = quad_broadcast(f4, lane);
     return HTMat<4,1,float>
            (HTMat<3,1,float>
             (HTMat<2,1,float>
              (HTMat<1,1,float>(HTVec<1,float>(res.x)), HTVec<1,float>(res.y)), HTVec<1,float>(res.z)), HTVec<1,float>(res.w));
  }
*/



  // broadcast across quad WARNING: dysnamic shuffle is SLOW
  template <int S, typename T>
  auto BroadcastQuad (HTVec<S,T> m, unsigned qlane, unsigned from) -> HTVec<S,T>
  {
    if constexpr (S==0)
      return {};
    else if constexpr (S==1)
      return simd_shuffle(m.Head(), 4*from+qlane);
    else
      return { BroadcastQuad(m.Tail(), qlane, from), simd_shuffle(m.Head(), 4*from+qlane) };
  }

  template <int H, int W, typename T>
  auto BroadcastQuad (HTMat<H,W,T> m, unsigned qlane, unsigned from) -> HTMat<H,W,T>
  {
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
      return BroadcastQuad(m.Head(), qlane, from);
    else
      return { BroadcastQuad(m.Tail(), qlane, from), BroadcastQuad(m.Head(), qlane, from) };
  }


  // broadcast across quad WARNING: dysnamic shuffle is SLOW
  template <ushort mask, int S, typename T>
  auto ShuffleXor (HTVec<S,T> m) -> HTVec<S,T>
  {
    if constexpr (S==0)
      return {};
    else if constexpr (S==1)
      return simd_shuffle_xor(m.Head(), mask);
    else
      return { ShuffleXor<mask>(m.Tail()), simd_shuffle_xor(m.Head(), mask) };
  }

  template <ushort mask, int H, int W, typename T>
  auto ShuffleXor (HTMat<H,W,T> m) -> HTMat<H,W,T>
  {
    if constexpr (H==0)
      return {};
    else if constexpr (H==1)
      return ShuffleXor<mask> (m.Head());
    else
      return { ShuffleXor<mask>(m.Tail()), ShuffleXor<mask>(m.Head()) };
  }








  // tensor stuff 



  // Contract:   R_ab = w_i H_iab
  // contracts the Vec index i with a weight vector.
  // covariant (H1-Hessian / HCurl-gradient) correction:  w = F^{-T} uhat = uhat*Finv
  // generic in the element type: works for elements float (-> inner product),
  // HTVec (-> weighted vector sum), HTMat (-> weighted matrix sum)
  template <int S, typename TW, typename TE>
  auto Contract (HTVec<S,TW> w, HTVec<S,TE> H) -> TE
  {
    if constexpr (S==0)
      return {};
    else if constexpr (S==1)
      return w.Head() * H.Head();
    else
      return Contract(w.Tail(), H.Tail()) + w.Head() * H.Head();
  }

  // RowContract:   R_ib = H_iab u_a
  // contracts one matrix index with the bare reference vector u,
  // the Vec index i becomes the row index of the result.
  // Piola (HDiv-gradient) correction, before mixing with F^{-1}:
  //     grad-correction = Finv * RowContract(H, uhat) - Outer(uhat, d)
  // NOTE: implemented as  u * H(i)  ( = (H_i^T u)^T ),
  //       correct only for symmetric H_i !
  template <int S, int D, typename T>
  auto RowContract (HTVec<S, HTMat<D,D,T>> H, HTVec<D,T> u) -> HTMat<S,D,T>
  {
    auto row = u * H.Head();          // row i = (H_i u)^T   (H_i symmetric)
    if constexpr (S==1)
      return HTMat<1,D,T> (row);
    else
      return { RowContract(H.Tail(), u), row };
  }


  // TraceContract:   d_b = (F^{-1})_ai H_iab  =  FinvT_ia H_iab   ( = d^_b ln J )
  // double contraction: Vec index i and matrix index a against F^{-1}.
  // pass FinvT = Trans(Finv), so row i of FinvT (= column i of F^{-1})
  // pairs index-aligned with H(i)  (both Head() = index S-1)
  template <int S, int D, typename T>
  auto TraceContract (HTVec<S, HTMat<D,D,T>> H, HTMat<S,D,T> FinvT) -> HTVec<D,T>
  {
    auto term = FinvT.Head() * H.Head();
    if constexpr (S==1)
      return term;
    else
      return TraceContract(H.Tail(), FinvT.Tail()) + term;
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

  template <int S, typename T>
  TB_HD T InnerProduct (Vec<S,T> a, Vec<S,T> b) {
    T sum = 0;
    for (int i = 0; i < S; i++) sum += a(i)*b(i);
    return sum;
  }

  template <int S, typename T>
  TB_HD T Norm (Vec<S,T> v) { return sqrt(InnerProduct(v, v)); }


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


  // rectangular mapping matrices of manifold elements, as in ngbla: Det is
  // the measure sqrt(det(F^T F)) and Inv the pseudo-inverse (F^T F)^{-1} F^T,
  // so 1/J * F * u, Trace(Inv(F) * dF), ... read the same for volume and
  // surface elements

  template <typename T>
  TB_HD T Det(Mat<3,2,T> m) {       // |F0 x F1|
    T c0 = m(1,0)*m(2,1) - m(2,0)*m(1,1);
    T c1 = m(2,0)*m(0,1) - m(0,0)*m(2,1);
    T c2 = m(0,0)*m(1,1) - m(1,0)*m(0,1);
    return sqrt(c0*c0 + c1*c1 + c2*c2);
  }

  template <typename T>
  TB_HD T Det(Mat<2,1,T> m) { return sqrt(m(0,0)*m(0,0) + m(1,0)*m(1,0)); }

  template <typename T>
  TB_HD T Det(Mat<3,1,T> m) { return sqrt(m(0,0)*m(0,0) + m(1,0)*m(1,0) + m(2,0)*m(2,0)); }

  template <int H, int W, typename T>
  TB_HD Mat<W,H,T> PseudoInv(Mat<H,W,T> m) {
    Mat<W,H,T> mt = Trans(m);
    return Inv(mt * m) * mt;
  }

  template <typename T> TB_HD Mat<2,3,T> Inv(Mat<3,2,T> m) { return PseudoInv(m); }
  template <typename T> TB_HD Mat<1,2,T> Inv(Mat<2,1,T> m) { return PseudoInv(m); }
  template <typename T> TB_HD Mat<1,3,T> Inv(Mat<3,1,T> m) { return PseudoInv(m); }


  template <int H, int W, typename T>
  TB_HD auto ToMat(Vec<H*W,T> vec) { return Mat<H,W,T>(vec); }

  template <int H, int W, typename T>
  TB_HD auto ToVec(Mat<H,W,T> mat) {
    Vec<H*W,T> vec;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        vec(i*W+j) = mat(i,j);
    return vec;
  }

  // --- shape-typed views of strided memory ------------------------------------
  // P is the pointer type: device memory by default, pass e.g.
  // const threadgroup T * for a view into shared memory
  template <int H, int W, typename T, typename P = const TB_DEVICE T *>
  struct DoubleSliceMat {
    P data; int dr, dc;
    TB_HD T operator() (int i, int j) const { return data[i*dr+j*dc]; }
  };
  template <int S, typename T, typename P = const TB_DEVICE T *>
  struct SliceVec {
    P data; int d;
    TB_HD T operator() (int i) const { return data[i*d]; }
  };
  // rank 3, (i, a, b)
  template <int S, int H, int W, typename T, typename P = const TB_DEVICE T *>
  struct SliceTens3 {
    P data; int d0, d1, d2;
    TB_HD T operator() (int i, int a, int b) const { return data[i*d0+a*d1+b*d2]; }
  };

  // rank 3 with a symmetric pair of trailing indices, (i, a, b) = (i, b, a),
  // stored by pair p(a,b) = a*D - a*(a-1)/2 + (b-a) for a <= b (D(D+1)/2 pairs)
  template <int S, int D, typename T, typename P = const TB_DEVICE T *>
  struct SymSliceTens3 {
    P data; int d0, d1;    // node stride, pair stride
    static TB_HD int Pair (int a, int b) { return (a <= b) ? a*D - a*(a-1)/2 + (b-a) : b*D - b*(b-1)/2 + (a-b); }
    TB_HD T operator() (int i, int a, int b) const { return data[i*d0 + Pair(a,b)*d1]; }
  };

  // sum_k a(k) * b(k)
  template <int K, typename T, typename FA, typename FB>
  TB_HD T Contract (FA a, FB b) {
    T sum = T(0);
    for (int k = 0; k < K; k++) sum += a(k) * b(k);
    return sum;
  }

  template <int H, int K, int W, typename T, typename PA, typename PB>
  TB_HD Mat<H,W,T> operator* (DoubleSliceMat<H,K,T,PA> a, DoubleSliceMat<K,W,T,PB> b) {
    Mat<H,W,T> r;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        r(i,j) = Contract<K,T>([=] (int k) { return a(i,k); }, [=] (int k) { return b(k,j); });
    return r;
  }
  template <int H, int K, typename T, typename PA, typename PX>
  TB_HD Vec<H,T> operator* (DoubleSliceMat<H,K,T,PA> a, SliceVec<K,T,PX> x) {
    Vec<H,T> v;
    for (int i = 0; i < H; i++)
      v(i) = Contract<K,T>([=] (int k) { return a(i,k); }, x);
    return v;
  }
  // matrix times rank-3, sliced by the last index:  (A T)(b) = A T(:,:,b).
  // With the coefficients A and the basis second derivatives T(k,a,b) this is
  // the derivative of the mapping matrix, dF(b) = d_b F
  template <int S, int K, int D, typename T, typename PA, typename PT>
  TB_HD Vec<D,Mat<S,D,T>> operator* (DoubleSliceMat<S,K,T,PA> a, SliceTens3<K,D,D,T,PT> t) {
    Vec<D,Mat<S,D,T>> r;
    for (int b = 0; b < D; b++)
      for (int i = 0; i < S; i++)
        for (int p = 0; p < D; p++)
          r(b)(i,p) = Contract<K,T>([=] (int k) { return a(i,k); }, [=] (int k) { return t(k,p,b); });
    return r;
  }

  // the same for symmetric second derivatives: one contraction per pair (a,b)
  // fills dF(b)(i,a) and dF(a)(i,b)
  template <int S, int K, int D, typename T, typename PA, typename PT>
  TB_HD Vec<D,Mat<S,D,T>> operator* (DoubleSliceMat<S,K,T,PA> a, SymSliceTens3<K,D,T,PT> t) {
    Vec<D,Mat<S,D,T>> r;
    for (int p = 0; p < D; p++)
      for (int b = p; b < D; b++)
        for (int i = 0; i < S; i++)
          {
            T v = Contract<K,T>([=] (int k) { return a(i,k); }, [=] (int k) { return t(k,p,b); });
            r(b)(i,p) = v;
            r(p)(i,b) = v;
          }
    return r;
  }

  // --- Vec of matrices ---------------------------------------------------------
  // (A V)(b) = A V(b) is the generic  T * Vec<S,T>  above with T a Mat
  template <int S, int H, int W, typename T>
  TB_HD Vec<S,Mat<W,H,T>> Trans (Vec<S,Mat<H,W,T>> v) {
    Vec<S,Mat<W,H,T>> r;
    for (int b = 0; b < S; b++) r(b) = Trans(v(b));
    return r;
  }
  // column b of the result is V(b) u  (e.g. dF u = d_b (F u) with u frozen)
  template <int D, int S, typename T>
  TB_HD Mat<S,D,T> operator* (Vec<D,Mat<S,D,T>> v, Vec<D,T> u) {
    Mat<S,D,T> r;
    for (int b = 0; b < D; b++)
      {
        Vec<S,T> vu = v(b) * u;
        for (int i = 0; i < S; i++) r(i,b) = vu(i);
      }
    return r;
  }
  // its adjoint: sum_b V(b) (column b of Q),  <Q, V u> = <Trans(V) Q, u>
  template <int D, int S, typename T>
  TB_HD Vec<D,T> operator* (Vec<D,Mat<D,S,T>> v, Mat<S,D,T> q) {
    Vec<D,T> r; r = T(0);
    for (int b = 0; b < D; b++)
      {
        Vec<S,T> qb;
        for (int i = 0; i < S; i++) qb(i) = q(i,b);
        r += v(b) * qb;
      }
    return r;
  }
  template <int D, typename T>
  TB_HD T Trace (Mat<D,D,T> m) {
    T t = T(0);
    for (int i = 0; i < D; i++) t += m(i,i);
    return t;
  }
  template <int S, int D, typename T>
  TB_HD Vec<S,T> Trace (Vec<S,Mat<D,D,T>> v) {
    Vec<S,T> r;
    for (int b = 0; b < S; b++) r(b) = Trace(v(b));
    return r;
  }

  // (a b^T)_ij = a_i b_j
  template <int H, int W, typename T>
  TB_HD Mat<H,W,T> Outer (Vec<H,T> a, Vec<W,T> b) {
    Mat<H,W,T> r;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        r(i,j) = a(i) * b(j);
    return r;
  }

  enum ORDERING { ColMajor, RowMajor };
  constexpr ORDERING Transposed (ORDERING o) { return (o==RowMajor) ? ColMajor : RowMajor; }





  template <ORDERING ORD, typename Tp, typename Tld = uint>
  class BareMatrix
  {
     Tp data;
     Tld ld;
     using ElementType = remove_addrspace_t<remove_cv_t<remove_reference_t<decltype(*data)>>>;

  public:
     BareMatrix (Tp _data, Tld _ld) : data(_data), ld(_ld) { }

     unsigned Offset (unsigned r, unsigned c) const {
       if constexpr (ORD==RowMajor) return r*ld+c;
       else return c*ld+r;
     }

     auto operator() (unsigned r, unsigned c) const { return data[Offset(r,c)]; }
     auto Addr(unsigned r, unsigned c) const { return data + Offset(r,c); }

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


     template <unsigned H, unsigned W, unsigned DIST>
     auto GetTileSlice(unsigned r, ushort qlane) const {
       static_assert(ORD==RowMajor);
       HTMat<H,W,ElementType> res;
       res.template LoadSlice<DIST> (Addr(r, qlane), ld);
       return res;
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


     template <unsigned H, unsigned W, unsigned DIST>
     void SetTileSlice(unsigned r, ushort qlane, HTMat<H,W,ElementType> mat) {
      static_assert(ORD==RowMajor);
      /*
       HTMat<H,W,ElementType> res;
       res.template LoadSlice<DIST> (Addr(r, qlane), ld);
       return res;
*/
       mat.template StoreSlice<DIST> (Addr(r, qlane), ld);
     }




     static constant constexpr bool is_colmajor = (ORD==ColMajor);
     auto SubMatrix (unsigned r, unsigned c) const { return BareMatrix { data+Offset(r,c), ld }; }
     auto Transpose() const { return BareMatrix<Transposed(ORD), Tp,Tld> { data, ld }; }
     auto LD() const { return ld; }
     Tp Data() const { return data; }
     bool IsColMajor() const { return ORD==ColMajor; }
  };

  template <ORDERING ORD, typename Tp, typename Tld>
  inline auto MakeBareMatrix (Tp data, Tld ld) {
    return BareMatrix<ORD,Tp,Tld> (data, ld);
  }

  template <ORDERING ORD, typename Tp, uint W>
  inline auto MakeBareMatrix (Tp data[][W]) {
    return MakeBareMatrix<ORD> (&(data[0][0]), ic<W>());
  }


  /*
    RotMatrix: row-major matrix in group memory with the 16-byte chunks of every
    row rotated by rota*row + rotb*(row/8) chunks. With a row length that is a
    multiple of 8 reals a plain layout puts rows 8 apart on the same banks; the
    rotation spreads them without padding and keeps 16-byte alignment. Measured
    rotations: (1,1) for operands read as column chunks by the quad-broadcast
    AddMM, (2,1) for the mma fragment loads. Tiles taken through GetTile/SetTile
    have contiguous rows inside one chunk: column offsets and widths must divide
    the chunk (4 floats / 2 doubles).
  */
  template <typename Tp, unsigned LD, unsigned ROTA, unsigned ROTB, bool XORSWZ> class RotMatrixT;

  // XORSWZ: xor the chunk index instead of rotating (one bit operation per access, chunk count must be a power of two)
  template <typename Tp, unsigned LD, unsigned ROTA, unsigned ROTB, bool XORSWZ = false>
  class RotMatrix
  {
    Tp data;
    unsigned r0, c0;
  public:
    using ElementType = remove_addrspace_t<remove_cv_t<remove_reference_t<decltype(*data)>>>;
    static constant constexpr unsigned CH = 16/sizeof(ElementType);
    static constant constexpr unsigned NCH = LD/CH;
    static constant constexpr bool is_colmajor = false;
    static_assert(LD % CH == 0, "RotMatrix: row length must be a multiple of the 16-byte chunk");

    RotMatrix (Tp _data, unsigned _r0 = 0, unsigned _c0 = 0) : data(_data), r0(_r0), c0(_c0) { }

    unsigned Offset (unsigned r, unsigned c) const
    {
      r += r0; c += c0;
      unsigned ch;
      if constexpr (XORSWZ)
        {
          static_assert((NCH & (NCH-1)) == 0, "RotMatrix xor swizzle: chunks per row must be a power of two");
          ch = (c/CH) ^ ((ROTA*r + ROTB*(r>>3)) & (NCH-1));
        }
      else
        {
          unsigned rot = (ROTA*r + ROTB*(r>>3)) % NCH;   // per row, hoisted when r is loop invariant
          ch = c/CH + rot;                                // both < NCH: one conditional wrap instead of a modulo
          if (ch >= NCH) ch -= NCH;
        }
      return r*LD + ch*CH + c%CH;
    }
    auto operator() (unsigned r, unsigned c) const { return data[Offset(r,c)]; }
    auto Addr (unsigned r, unsigned c) const { return data + Offset(r,c); }
    RotMatrix SubMatrix (unsigned r, unsigned c) const { return RotMatrix (data, r0+r, c0+c); }
    RotMatrixT<Tp,LD,ROTA,ROTB,XORSWZ> Transpose() const;
    bool IsColMajor() const { return false; }
    Tp Data() const { return data; }

    /*
      tile rows of w elements start at column offsets that are multiples of w
      inside one chunk, so the address is w*sizeof(T) aligned: tell the compiler,
      it then issues one vector load per row instead of w scalar ones
    */
    template <unsigned w>
    auto AlignedAddr (unsigned r, unsigned c) const
    {
      auto p = Addr(r, c);
#ifdef TB_CUDA
      return (decltype(p)) __builtin_assume_aligned (p, (w*sizeof(ElementType) > 16) ? 16 : w*sizeof(ElementType));
#else
      return p;
#endif
    }
    template <unsigned h, unsigned w>
    HTMat<h,w,ElementType> GetTile (unsigned r, unsigned c) const
    {
      HTVec<w,ElementType> v;
      v.Load (AlignedAddr<w>(r+h-1, c));
      if constexpr (h == 1) return HTMat<1,w,ElementType> (v);
      else return HTMat<h,w,ElementType> (GetTile<h-1,w>(r, c), v);
    }
    template <int h, int w, typename T2>
    void SetTile (unsigned r, unsigned c, HTMat<h,w,T2> tile)
    {
      if constexpr (h > 1) SetTile<h-1,w> (r, c, tile.Tail());
      tile.Head().Store (AlignedAddr<w>(r+h-1, c));
    }
  };

  // column-major view of a RotMatrix
  template <typename Tp, unsigned LD, unsigned ROTA, unsigned ROTB, bool XORSWZ>
  class RotMatrixT
  {
    RotMatrix<Tp,LD,ROTA,ROTB,XORSWZ> m;
  public:
    using ElementType = typename RotMatrix<Tp,LD,ROTA,ROTB,XORSWZ>::ElementType;
    static constant constexpr bool is_colmajor = true;
    RotMatrixT (RotMatrix<Tp,LD,ROTA,ROTB,XORSWZ> _m) : m(_m) { }
    auto operator() (unsigned r, unsigned c) const { return m(c,r); }
    auto Addr (unsigned r, unsigned c) const { return m.Addr(c,r); }
    RotMatrixT SubMatrix (unsigned r, unsigned c) const { return RotMatrixT (m.SubMatrix(c,r)); }
    RotMatrix<Tp,LD,ROTA,ROTB,XORSWZ> Transpose() const { return m; }
    bool IsColMajor() const { return true; }
    template <unsigned h, unsigned w>
    auto GetTile (unsigned r, unsigned c) const { return Trans (m.template GetTile<w,h>(c, r)); }
    template <int h, int w, typename T2>
    void SetTile (unsigned r, unsigned c, HTMat<h,w,T2> tile) { m.template SetTile<w,h> (c, r, Trans(tile)); }
  };

  template <typename Tp, unsigned LD, unsigned ROTA, unsigned ROTB, bool XORSWZ>
  RotMatrixT<Tp,LD,ROTA,ROTB,XORSWZ> RotMatrix<Tp,LD,ROTA,ROTB,XORSWZ>::Transpose() const { return RotMatrixT<Tp,LD,ROTA,ROTB,XORSWZ> (*this); }

  template <unsigned LD, unsigned ROTA, unsigned ROTB, bool XORSWZ = false, typename Tp>
  inline auto MakeRotMatrix (Tp data) { return RotMatrix<Tp,LD,ROTA,ROTB,XORSWZ> (data); }





  template <unsigned H, unsigned W, typename T = float>
  class WarpMatrix
  {
    static_assert(H%8==0, "WarpMatrix height must be a multiple of 8");
    static_assert(W%4==0, "WarpMatrix width must be a multiple of 4");
    static constant constexpr unsigned BW = W/4; // sizeof(float2)/sizeof(T);
    static constant constexpr unsigned BH = H/8;
    HTMat<BH,BW,T> myvals; 
  public:
    WarpMatrix() { }

    WarpMatrix(T val) { myvals = val; }

    template <typename M>
    WarpMatrix(M mat, uint tid) {
      myvals = mat.template GetTile<BH,BW>(MyRow(tid),MyCol(tid));
    }

    void operator= (T val) { myvals = val; }

    WarpMatrix (HTMat<BH,BW,T> amyvals) : myvals(amyvals) { } 
    auto GetValues() const { return myvals; }

    template <uint K, typename M1, typename M2>
    void AddMM(M1 m1, M2 m2, uint tid)
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
      if constexpr (!M2::is_colmajor) {
      // if constexpr (true) {
        constexpr int KTILE = 1;
        uint k = 0;
        for ( ; k+4 <= K; k+=4)
          {
            auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + (tid&3));
            for (uint k1 = 0; k1 < 4; k1++)
               myvals = FMA (QuadBroadcast(ATileQuad, k1), m2.template GetTile<KTILE,BW>(k+k1,c), myvals);
          }
        // remainder
        if constexpr (K%4 != 0)
          for ( ;  k < K; k++)
            myvals = FMA (m1.template GetTile<BH,KTILE>(r,k), m2.template GetTile<KTILE,BW>(k,c), myvals);
      } else {

        constexpr int KTILE = 1;
        uint k = 0;
        // auto tmpmyvals = myvals;
        for ( ; k+4 <= K; k+=4)
          {
            auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + (tid&3));
            auto BTileQuad = m2.template GetTile<4*KTILE,BW>(k, c);


            myvals = FMA (QuadBroadcast(ATileQuad, 0), HTMat<1,BW,T> (BTileQuad.template GetRow<0>()), myvals);
            myvals = FMA (QuadBroadcast(ATileQuad, 1), HTMat<1,BW,T> (BTileQuad.template GetRow<1>()), myvals);
            myvals = FMA (QuadBroadcast(ATileQuad, 2), HTMat<1,BW,T> (BTileQuad.template GetRow<2>()), myvals);
            myvals = FMA (QuadBroadcast(ATileQuad, 3), HTMat<1,BW,T> (BTileQuad.template GetRow<3>()), myvals);
/*
            tmpmyvals = FMA (QuadBroadcast(ATileQuad, 0), HTMat<1,BW,float> (BTileQuad.template GetRow<0>()), 
                       FMA (QuadBroadcast(ATileQuad, 1), HTMat<1,BW,T> (BTileQuad.template GetRow<1>()), 
                         FMA (QuadBroadcast(ATileQuad, 2), HTMat<1,BW,T> (BTileQuad.template GetRow<2>()), 
                           FMA (QuadBroadcast(ATileQuad, 3), HTMat<1,BW,T> (BTileQuad.template GetRow<3>()), tmpmyvals))));
*/
            // for (uint k1 = 0; k1 < 4; k1++)
            // myvals = FMA (QuadBroadcast(ATileQuad, k1), m2.template GetTile<KTILE,BW>(k+k1,c), myvals);
          }
        // myvals = tmpmyvals;

        // remainder
        if constexpr (K%4 != 0)
          for ( ;  k < K; k++)
            myvals = FMA (m1.template GetTile<BH,KTILE>(r,k), m2.template GetTile<KTILE,BW>(k,c), myvals);

      }



/*
      static_assert (K%4==0, "K must be a multiple of 4");
      // BTile: load 4 contiguous elements into quad with HTVec.LoadSlice (needs adaption for C)
      // 8x8x8:  2.0 -> 2.5 TF
      constexpr int KTILE = 1;
      uint k = 0;
      m2 = m2.ShiftCols((tid&3));
      for ( ; k+4 <= K; k+=4)
        {
          auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + (tid&3));
          for (uint k1 = 0; k1 < 4; k1++)
            {
              auto BTile = HTMat<1,BW,T> (HTVec<BW,float>().template LoadSlice<4> (m2.Addr(0,0)));
              m2 = m2.template ShiftRows<1>();
              myvals = FMA (QuadBroadcast(ATileQuad, k1), BTile, myvals);
           }
        }
*/

/*
      // BTile: load 4 contiguous elements into quad with HTVec.LoadSize (needs adaption for C)
      // 8x8x8:  2.0 -> 2.5 TF
      // with unrolling -> no difference
      constexpr int KTILE = 1;
      uint k = 0;
      m2 = m2.ShiftCols((tid&3));
      for ( ; k+4 <= K; k+=4)
        {
          auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + (tid&3));

          auto BTile0 = HTMat<1,BW,T> (HTVec<BW,float>().template LoadSlice<4> (m2.Addr(0,0)));
          auto BTile1 = HTMat<1,BW,T> (HTVec<BW,float>().template LoadSlice<4> (m2.Addr(1,0)));
          auto BTile2 = HTMat<1,BW,T> (HTVec<BW,float>().template LoadSlice<4> (m2.Addr(2,0)));
          auto BTile3 = HTMat<1,BW,T> (HTVec<BW,float>().template LoadSlice<4> (m2.Addr(3,0)));
          auto ATile0 = QuadBroadcast(ATileQuad, 0);
          auto ATile1 = QuadBroadcast(ATileQuad, 1);
          auto ATile2 = QuadBroadcast(ATileQuad, 2);
          auto ATile3 = QuadBroadcast(ATileQuad, 3);
          myvals = FMA (ATile0, BTile0, myvals);
          myvals = FMA (ATile1, BTile1, myvals);
          myvals = FMA (ATile2, BTile2, myvals);
          myvals = FMA (ATile3, BTile3, myvals);

          m2 = m2.template ShiftRows<4>();
        }
*/







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


    template <typename M>
    void Store(M mat, unsigned tid) {
      mat.template SetTile<BH,BW>(MyRow(tid),MyCol(tid), myvals);
    }

    static auto MyCol(uint tid) { return BW*(tid&3); }
    static auto MyRow(uint tid) { return BH*((tid>>2)&7); }
  };


  template <uint K, uint H, uint W, typename T, ORDERING ORD1, typename Tp1, typename Tld1,  ORDERING ORD2, typename Tp2, typename Tld2>
  inline auto AddMM(WarpMatrix<H,W,T> m, BareMatrix<ORD1,Tp1,Tld1> m1, BareMatrix<ORD2,Tp2,Tld2> m2, uint tid)
  {
    static_assert(H%8==0, "AddMM height must be a multiple of 8");
    static_assert(W%4==0, "AddMM width must be a multiple of 4");
     constexpr unsigned BW = W/4; // sizeof(float2)/sizeof(T);
     constexpr unsigned BH = H/8;


     auto myvals = m.GetValues();
     auto r = m.MyRow(tid);
     auto c = m.MyCol(tid);

     // tmp.template AddMM<K> (m1, m2, tid);

      if constexpr (ORD2==RowMajor) {
       constexpr int KTILE = 1;
       uint k = 0;
#pragma unroll 2
       for ( ; k+4 <= K; k+=4)
          {
            auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + (tid&3));
            for (uint k1 = 0; k1 < 4; k1++)
               myvals = FMA (QuadBroadcast(ATileQuad, k1), m2.template GetTile<KTILE,BW>(k+k1,c), myvals);
          }
        // remainder
        if constexpr (K%4 != 0)
          for ( ;  k < K; k++)
            myvals = FMA (m1.template GetTile<BH,KTILE>(r,k), m2.template GetTile<KTILE,BW>(k,c), myvals);
      } else {

        constexpr int KTILE = 1;
        uint k = 0;
        for ( ; k+4 <= K; k+=4)
          {
            auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + (tid&3));
            auto BTileQuad = m2.template GetTile<4*KTILE,BW>(k, c);
            myvals = FMA (QuadBroadcast(ATileQuad, 0), HTMat<1,BW,T> (BTileQuad.template GetRow<0>()), myvals);
            myvals = FMA (QuadBroadcast(ATileQuad, 1), HTMat<1,BW,T> (BTileQuad.template GetRow<1>()), myvals);
            myvals = FMA (QuadBroadcast(ATileQuad, 2), HTMat<1,BW,T> (BTileQuad.template GetRow<2>()), myvals);
            myvals = FMA (QuadBroadcast(ATileQuad, 3), HTMat<1,BW,T> (BTileQuad.template GetRow<3>()), myvals);
          }

        // remainder
        if constexpr (K%4 != 0)
          for ( ;  k < K; k++)
            myvals = FMA (m1.template GetTile<BH,KTILE>(r,k), m2.template GetTile<KTILE,BW>(k,c), myvals);
      }

     return WarpMatrix<H,W,float> (myvals);
  }





  // interleaved accumulator entries and B, B-mat must be RowMajor (by now)
  template <unsigned H, unsigned W, typename T = float>
  class WarpMatrixV2
  {
    static_assert(H%8==0, "WarpMatrixV2 height must be a multiple of 8");
    static_assert(W%4==0, "WarpMatrixV2 width must be a multiple of 4");
    static constant constexpr unsigned BW = W/4;
    static constant constexpr unsigned BH = H/8;
    HTMat<BH,BW,T> myvals; 
  public:
    WarpMatrixV2() { }
    WarpMatrixV2(T val) { myvals = val; }

    template <typename Tp, typename Tld>
    WarpMatrixV2(BareMatrix<RowMajor,Tp,Tld> mat, uint tid)
      : myvals { mat.template GetTileSlice<BH,BW,4>(MyRow(tid), tid&3) } { }

    void operator= (T val) { myvals = val; }

    template <uint K, ORDERING ORD1, typename Tp1, typename Tld1, typename Tp2, typename Tld2>
    void AddMM(BareMatrix<ORD1,Tp1,Tld1> m1, BareMatrix<RowMajor,Tp2,Tld2> m2, uint tid)
    {
      auto r = MyRow(tid);

      // static_assert (K%4==0, "K must be a multiple of 4");
      constexpr int KTILE = 1;
      uint k = 0;
      m2 = m2.ShiftCols((tid&3));
      for ( ; k+4 <= K; k+=4)
        {
          auto ATileQuad = m1.template GetTile<BH,KTILE>(r,k + (tid&3));
          for (uint k1 = 0; k1 < 4; k1++)
            {
              // HTMat<1,BW,T> BTile;
              // BTile.template LoadSlice<4> (m2.Addr(0,0), m2.LD());
              auto BTile = m2.template GetTileSlice<KTILE,BW,4> (0,0);
              m2 = m2.template ShiftRows<1>();
              myvals = FMA (QuadBroadcast(ATileQuad, k1), BTile, myvals);
            }
        }
      for ( ;  k < K; k++)
        {
          auto ATile = m1.template GetTile<BH,KTILE>(r,k);
          auto BTile = m2.template GetTileSlice<KTILE,BW,4> (0,0);
          m2 = m2.template ShiftRows<1>();
          myvals = FMA (ATile, BTile, myvals);
        }

    }

    template <typename Tp, typename Tld>
    void Store(BareMatrix<RowMajor,Tp,Tld> mat, unsigned tid) {
      mat.template SetTileSlice<BH,BW,4> (MyRow(tid), tid&3, myvals);
    }

    // auto MyCol(uint tid) const { return BW*(tid&3); }
    auto MyRow(uint tid) const { return BH*((tid>>2)&7); }
  };



  template <unsigned W>
  class WarpMatrixV2<4,W,float>
  {
    typedef float T;
    static_assert(W%8==0);
    static constant constexpr unsigned BW = W/8;
    static constant constexpr unsigned BH = 1;
    HTMat<BH,BW,float> myvals; 
  public:
    WarpMatrixV2() { }
    WarpMatrixV2(float val) { myvals = val; }

    template <typename Tp, typename Tld>
    WarpMatrixV2(BareMatrix<RowMajor,Tp,Tld> mat, uint tid)
      : myvals { mat.template GetTileSlice<BH,BW,8>(MyRow(tid), tid&7) } { }

    void operator= (T val) { myvals = val; }

    template <ORDERING ORD1, typename Tp1, typename Tld1, typename Tp2, typename Tld2>
    void AddMM(uint K, BareMatrix<ORD1,Tp1,Tld1> m1, BareMatrix<RowMajor,Tp2,Tld2> m2, uint tid)
    {
      auto r = MyRow(tid);

      // static_assert (K%4==0, "K must be a multiple of 4");
      constexpr int KTILE = 1;
      m2 = m2.ShiftCols((tid&7));
      for (uint k=0;  k < K; k++)
        {
          auto ATile = m1.template GetTile<BH,KTILE>(r,k);
          auto BTile = m2.template GetTileSlice<KTILE,BW,8> (0,0);
          m2 = m2.template ShiftRows<1>();
          myvals = FMA (ATile, BTile, myvals);
        }

    }

    template <typename Tp, typename Tld>
    void Store(BareMatrix<RowMajor,Tp,Tld> mat, unsigned tid) {
      mat.template SetTileSlice<BH,BW,8> (MyRow(tid), tid&7, myvals);
    }

    // auto MyCol(uint tid) const { return BW*(tid&3); }
    auto MyRow(uint tid) const { return BH*((tid>>3)&3); }
  };




#ifdef TB_METAL
  // simdgroup matrices exist only in metal: float tiles with rows and cols
  // multiples of 8 are blocks of simdgroup_float8x8. A partial specialisation,
  // so it must be guarded - it is parsed even if unused.
  template <unsigned H, unsigned W>
  class WarpMatrix<H,W,float>
  {
    static_assert(H%8==0 && W%8==0, "metal WarpMatrix<float>: rows and cols must be multiples of 8");
    typedef float T;
    static constant constexpr unsigned BH = H/8;
    static constant constexpr unsigned BW = W/8;
    metal::simdgroup_float8x8 m[BH][BW];
  public:
    WarpMatrix() { }
    WarpMatrix(T val) { *this = val; }
    WarpMatrix(threadgroup T * data, int ld, int tid) {
      Load (BareMatrix<RowMajor, threadgroup T*, int> (data, ld));
    }
    template <typename M>
    WarpMatrix(M mat, int tid) { Load(mat); }

    void operator= (T val) {
      for (uint i = 0; i < BH; i++)
        for (uint j = 0; j < BW; j++)
          m[i][j] = metal::make_filled_simdgroup_matrix<float, 8, 8>(val);
    }

    template <typename M>
    void Load(M mat) {
      for (uint i = 0; i < BH; i++)
        for (uint j = 0; j < BW; j++)
          {
            auto sub = mat.SubMatrix(8*i, 8*j);
            metal::simdgroup_load(m[i][j], sub.Data(), sub.LD(), ulong2(0,0), sub.IsColMajor());
          }
    }

    template <int K, typename M1, typename M2>
    void AddMM(M1 m1, M2 m2, int tid)
    {
      static_assert(K%8==0, "metal WarpMatrix::AddMM: K must be a multiple of 8");
      metal::simdgroup_float8x8 ma[BH], mb[BW];
      for (uint k = 0; k < K; k += 8) {
        for (uint i = 0; i < BH; i++) {
          auto sub = m1.SubMatrix(8*i, k);
          metal::simdgroup_load(ma[i], sub.Data(), sub.LD(), ulong2(0,0), sub.IsColMajor());
        }
        for (uint j = 0; j < BW; j++) {
          auto sub = m2.SubMatrix(k, 8*j);
          metal::simdgroup_load(mb[j], sub.Data(), sub.LD(), ulong2(0,0), sub.IsColMajor());
        }
        for (uint i = 0; i < BH; i++)
          for (uint j = 0; j < BW; j++)
            metal::simdgroup_multiply_accumulate(m[i][j], ma[i], mb[j], m[i][j]);
      }
    }

    void Store(threadgroup T * data, int ld, int tid) {
      Store (BareMatrix<RowMajor, threadgroup T*, int> (data, ld), tid);
    }

    template <typename M>
    void Store(M mat, int tid) {
      for (uint i = 0; i < BH; i++)
        for (uint j = 0; j < BW; j++)
          {
            auto sub = mat.SubMatrix(8*i, 8*j);
            metal::simdgroup_store(m[i][j], sub.Data(), sub.LD(), ulong2(0,0), sub.IsColMajor());
          }
    }
  };
#endif // TB_METAL

#if defined(TB_CUDA) && defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 800
  /*
    fp64 tensor cores (sm_80+). Fragment layout per lane (gid = lane/4, tig = lane%4):
    A(gid [+8], k+tig+4s), B(k+tig+4s, gid), C(gid [+8], 2*tig+{0,1}). The accumulator is
    kept as 8x8 blocks c[i][j][2]; the m16 shapes cover two row blocks at once.
    Shape by throughput on H100 (A fragments register-resident): m16n8k16 64 TF, m16n8k8
    62 TF, m8n8k4 32 TF; with A loaded per k-step: L1 34-36 TF, shared 46-50 TF. sm_90 has
    the m16 shapes, sm_80 only m8n8k4. Partial specialisation like the metal simdgroup one.
  */
  template <unsigned H, unsigned W>
  class WarpMatrix<H,W,double>
  {
    static_assert(H%8==0 && W%8==0, "cuda WarpMatrix<double>: rows and cols must be multiples of 8");
    typedef double T;
    static constant constexpr unsigned BH = H/8;
    static constant constexpr unsigned BW = W/8;
    T c[BH][BW][2];

    static __device__ __forceinline__ void MMA8x8x4 (T * c, T a, T b)
    {
      asm volatile("mma.sync.aligned.m8n8k4.row.col.f64.f64.f64.f64 {%0,%1}, {%2}, {%3}, {%0,%1};"
                   : "+d"(c[0]), "+d"(c[1]) : "d"(a), "d"(b));
    }
    static __device__ __forceinline__ void MMA16x8x8 (T * c0, T * c1, const T * a, const T * b)
    {
      asm volatile("mma.sync.aligned.m16n8k8.row.col.f64.f64.f64.f64 {%0,%1,%2,%3}, {%4,%5,%6,%7}, {%8,%9}, {%0,%1,%2,%3};"
                   : "+d"(c0[0]), "+d"(c0[1]), "+d"(c1[0]), "+d"(c1[1])
                   : "d"(a[0]), "d"(a[1]), "d"(a[2]), "d"(a[3]), "d"(b[0]), "d"(b[1]));
    }
    static __device__ __forceinline__ void MMA16x8x16 (T * c0, T * c1, const T * a, const T * b)
    {
      asm volatile("mma.sync.aligned.m16n8k16.row.col.f64.f64.f64.f64 {%0,%1,%2,%3}, {%4,%5,%6,%7,%8,%9,%10,%11}, {%12,%13,%14,%15}, {%0,%1,%2,%3};"
                   : "+d"(c0[0]), "+d"(c0[1]), "+d"(c1[0]), "+d"(c1[1])
                   : "d"(a[0]), "d"(a[1]), "d"(a[2]), "d"(a[3]), "d"(a[4]), "d"(a[5]), "d"(a[6]), "d"(a[7]),
                     "d"(b[0]), "d"(b[1]), "d"(b[2]), "d"(b[3]));
    }
  public:
    WarpMatrix() { }
    WarpMatrix(T val) { *this = val; }
    template <typename M>
    WarpMatrix(M mat, unsigned tid) { Load(mat, tid); }

    void operator= (T val) {
      for (unsigned i = 0; i < BH; i++)
        for (unsigned j = 0; j < BW; j++)
          c[i][j][0] = c[i][j][1] = val;
    }

    template <typename M>
    void Load(M mat, unsigned tid) {
      unsigned gid = (tid&31)>>2, tig = tid&3;
      for (unsigned i = 0; i < BH; i++)
        for (unsigned j = 0; j < BW; j++)
          for (unsigned l = 0; l < 2; l++)
            c[i][j][l] = mat(8*i+gid, 8*j+2*tig+l);
    }

    template <unsigned K, typename M1, typename M2>
    void AddMM(M1 m1, M2 m2, unsigned tid)
    {
      static_assert(K%4==0, "cuda WarpMatrix<double>::AddMM: K must be a multiple of 4");
      unsigned gid = (tid&31)>>2, tig = tid&3;
      /*
        TB_FP64_MMA_K selects the k-depth of the instruction (4: m8n8k4, 8: m16n8k8,
        16: m16n8k16). The m16 shapes reach 62 TF only with register-resident or shared
        A fragments; fed from L1 per k-step they are no faster than m8n8k4 in apply_btdtb
        (measured Sept 2026), so 4 stays the default until the operands are staged.
      */
#ifndef TB_FP64_MMA_K
#define TB_FP64_MMA_K 4
#endif
#if __CUDA_ARCH__ >= 900
      constexpr unsigned KS = (H%16 != 0 || K%TB_FP64_MMA_K != 0) ? 4 : TB_FP64_MMA_K;
#else
      constexpr unsigned KS = 4;
#endif
      // fully unrolled: rolled loops push the accumulator array into local memory
#pragma unroll
      for (unsigned k = 0; k < K; k += KS)
        {
          if constexpr (KS == 4)
            {
              T a[BH], b[BW];
              #pragma unroll
              for (unsigned i = 0; i < BH; i++) a[i] = m1(8*i+gid, k+tig);
              #pragma unroll
              for (unsigned j = 0; j < BW; j++) b[j] = m2(k+tig, 8*j+gid);
              #pragma unroll
              for (unsigned i = 0; i < BH; i++)
                for (unsigned j = 0; j < BW; j++)
                  MMA8x8x4 (c[i][j], a[i], b[j]);
            }
          else
            {
              T a[BH/2][KS/2], b[BW][KS/4];
              #pragma unroll
              for (unsigned i = 0; i < BH/2; i++)
                for (unsigned s = 0; s < KS/4; s++)
                  {
                    a[i][2*s]   = m1(16*i+gid,   k+tig+4*s);
                    a[i][2*s+1] = m1(16*i+gid+8, k+tig+4*s);
                  }
              #pragma unroll
              for (unsigned j = 0; j < BW; j++)
                for (unsigned s = 0; s < KS/4; s++)
                  b[j][s] = m2(k+tig+4*s, 8*j+gid);
              #pragma unroll
              for (unsigned i = 0; i < BH/2; i++)
                for (unsigned j = 0; j < BW; j++)
                  {
                    if constexpr (KS == 16) MMA16x8x16 (c[2*i][j], c[2*i+1][j], a[i], b[j]);
                    else                    MMA16x8x8  (c[2*i][j], c[2*i+1][j], a[i], b[j]);
                  }
            }
        }
    }

    template <typename M>
    void Store(M mat, unsigned tid) {
      unsigned gid = (tid&31)>>2, tig = tid&3;
      for (unsigned i = 0; i < BH; i++)
        for (unsigned j = 0; j < BW; j++)
          for (unsigned l = 0; l < 2; l++)
            *mat.Addr(8*i+gid, 8*j+2*tig+l) = c[i][j][l];
    }
  };
#endif // TB_CUDA fp64 tensor cores


} // namespace tinybla



