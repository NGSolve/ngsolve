#ifndef FILE_NORMALIZED_LEGENDRE
#define FILE_NORMALIZED_LEGENDRE

#include <bla.hpp>

namespace ngsfmm
{
  using namespace ngbla;

  /*
    sqrt( (2n+1) * (n-m)! / (n+m)! ) * P^m_n
  */
  class NormalizedLegendreFunctions
  {
    static Matrix<double,RowMajor> matA, matB;   // y(n) = A x y(n-1) + B y(n-2)
    static Vector<double> sqrt2np1, diagfac;   // sqrt(2n+1), sqrt((2m-1)/(2m))
  public:
    static void Calc(int mmax)
    {
      if (matA.Height() >= mmax+1) return;
      // width padded to SIMD blocks of up to 8 columns
      size_t w = (mmax+8) & ~size_t(7);
      matA.SetSize(mmax+1, w);
      matB.SetSize(mmax+1, w);
      // n <= m: A=0, B=1, the recursion y(n) = y(n-2) carries the seeds to the diagonal
      matA = 0.0;
      matB = 1.0;
      sqrt2np1.SetSize(mmax+1);
      diagfac.SetSize(mmax+1);
      for (int n = 0; n <= mmax; n++)
        sqrt2np1(n) = sqrt(2*n+1.0);
      diagfac(0) = 1;
      for (int m = 1; m <= mmax; m++)
        diagfac(m) = sqrt((2*m-1.0)/(2*m));
      for (int n = 1; n <= mmax; n++)
        {
          matA(n,n-1) = sqrt(2*n-1.0);
          matB(n,n-1) = 0.0;
          for (int m = 0; m+2 <= n; m++)
            {
              matA(n,m) = (2*n-1) / sqrt((n-m+0.0)*(n+m));
              matB(n,m) = -sqrt((n+m-1.0)*(n-m-1.0)) / sqrt((n-m+0.0)*(n+m));
            }
        }
    }
    
  private:
    template <typename T> struct SIMDElem { using type = double; };
    template <typename E, int N> struct SIMDElem<SIMD<E,N>> { using type = E; };

    // blocks of W columns, y(n,m), needs nmax < matA.Height()
    // scalar x: one SIMD<double,W> across the columns
    // SIMD x (over points): W independent T registers, one per column
    template <typename T, typename TY>
    static void EvalSIMD (size_t mmax, size_t nmax, T x, TY y)
    {
      constexpr bool scalar = std::is_same_v<T,double> || std::is_same_v<T,float>;
      using TS = std::conditional_t<scalar, double, T>;
      // fastest on M4: W=8 for scalar, 128 bytes in flight for SIMD x
      constexpr int W = scalar ? 8 : std::clamp(int(128/sizeof(T)), 1, 8);
      static_assert(W <= 8, "matA width padded to 8");
      using TV = std::conditional_t<scalar, SIMD<double,W>, std::array<T,W>>;
      // coefficient type: double, or float for SIMD<float> x
      using TE = typename SIMDElem<T>::type;

      TS xd = x;
      TS u = -sqrt((1-xd)*(1+xd));
      size_t mtop = std::min(mmax, nmax);
      auto scale = sqrt2np1.Data();
      TS dm = 1.0;   // running unscaled diagonal y(m,m)

      for (size_t m0 = 0; m0 <= mtop; m0 += W)
        {
          size_t nl = std::min<size_t>(W, mtop+1-m0);
          TS diag[W];
          for (int k = 0; k < W; k++)
            {
              size_t m = m0+k;
              if (m > 0 && m <= mtop)
                dm = dm*u*TE(diagfac(m));
              diag[k] = (m <= mtop) ? dm : TS(0.0);
            }
          TS sub0 = xd*diag[0]*TE(scale[m0]);
          // seed rows m0, m0+1: lane k gets diag in row m0+k%2, see Calc
          TS s0[W], s1[W];
          for (int k = 0; k < W; k++)
            {
              s0[k] = (k % 2) ? TS(0.0) : diag[k];
              s1[k] = (k % 2) ? diag[k] : TS(0.0);
            }
          s1[0] = sub0;
          TV v2, v1;
          if constexpr (scalar)
            {
              v2 = TV(&s0[0]);
              v1 = TV(&s1[0]);
            }
          else
            for (int k = 0; k < W; k++)
              {
                v2[k] = s0[k];
                v1[k] = s1[k];
              }

          y(m0,m0) = TE(scale[m0])*diag[0];
          if (m0+1 <= nmax)
            {
              y(m0+1,m0) = TE(scale[m0+1])*sub0;
              if (nl > 1) y(m0+1,m0+1) = TE(scale[m0+1])*diag[1];
            }

          for (size_t n = m0+2; n <= nmax; n++)
            {
              if constexpr (scalar)
                {
                  SIMD<double,W> a(&matA(n,m0)), b(&matB(n,m0));
                  SIMD<double,W> v = FMA(a * xd, v1, b * v2);
                  v2 = v1;
                  v1 = v;
                }
              else
                {
                  const double * a = &matA(n,m0);
                  const double * b = &matB(n,m0);
                  for (int k = 0; k < W; k++)
                    {
                      T v = FMA(TE(a[k])*x, v1[k], TE(b[k])*v2[k]);
                      v2[k] = v1[k];
                      v1[k] = v;
                    }
                }
              // lanes above the diagonal carry seeds, not stored
              size_t nst = std::min<size_t>(nl, n-m0+1);
              if constexpr (scalar)
                {
                  SIMD<double,W> s = scale[n] * v1;
                  if constexpr (TY::Ordering() == RowMajor && std::is_same_v<T,double>)
                    {
                      // contiguous row: explicit SIMD store, the loop below becomes memcpy
                      if (nst == W)
                        s.Store(&y(n,m0));
                      else
                        s.Store(&y(n,m0), SIMD<mask64,W>(int64_t(nst)));
                    }
                  else if (nst == W)
                    for (int k = 0; k < W; k++)
                      y(n,m0+k) = s[k];
                  else
                    for (size_t k = 0; k < nst; k++)
                      y(n,m0+k) = s[k];
                }
              else if (nst == W)
                for (int k = 0; k < W; k++)
                  y(n,m0+k) = TE(scale[n])*v1[k];
              else
                for (size_t k = 0; k < nst; k++)
                  y(n,m0+k) = TE(scale[n])*v1[k];
            }
        }
    }

  public:
    // values(n,m) for 0 <= m <= n <= nmax, entries m > n are not touched
    // row-major values give contiguous stores
    template <typename T, typename S>
    static void EvalAll (size_t nmax, T x, S && values)
    {
      size_t mmax = nmax;
      auto y = values.View();

      T u = -sqrt((1-x)*(1+x));
      y(0,0)=1;

      if constexpr (std::is_same_v<T,double> || std::is_same_v<T,float> || requires { T::Size(); })
        if (nmax+1 <= matA.Height())
          {
            EvalSIMD (mmax, nmax, x, y);
            return;
          }

      if (mmax+1 > matA.Height() || nmax > matA.Width())
        {
          for (int m = 0; m <= mmax; m++)
            {
              if (m > 0)
                y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0)/(2*m));
              if (m < nmax)
                y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0);
              for (int n = m+2; n <= nmax; n++)
                y(n,m)=((2*n-1)*x*y(n-1,m) - 
                        sqrt((n+m-1.0)*(n-m-1.0))*y(n-2,m))
                  /sqrt((n-m+0.0)*(n+m));
            }
        }
      else
        {
          for (int m = 0; m <= mmax; m++)
            {
              if (m > 0)
                y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0)/(2*m));
              if (m < nmax)
                {
                  T valold = y(m,m);
                  
                  T val = x*valold*sqrt(2*m+1.0);
                  y(m+1,m)=val;
                  
                  auto coefsA = matA.Col(m);
                  auto coefsB = matB.Col(m);
                  
                  for (int n = m+2; n <= nmax; n++)
                    {
                      T valnew = coefsA(n) * x * val + coefsB(n) * valold;
                      y(n,m) = valnew;
                      valold = val;
                      val = valnew;
                    }
                }
              /*
              for (int n = m+2; n <= nmax; n++)
                y(n,m) = matA(n,m) * x * y(n-1,m) + matB(n,m) * y(n-2,m);
              */
            }
        }
      
      for (int n = 0; n <= nmax; n++)
        y.Row(n).Range(0,n+1) *= sqrt(2*n+1.0);
        /*
        for (int m = 0; m <= n; m++)
          y(n,m)=y(n,m)*sqrt(2*n+1.0);
        */
    }
  };
}

#endif
