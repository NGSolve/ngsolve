#ifndef FILE_CHOLESKY
#define FILE_CHOLESKY

/****************************************************************************/
/* File:   cholesky.hpp                                                     */
/* Author: Joachim Schoeberl                                                */
/* Date:   25. Mar. 2000, 16. June 2002                                     */
/****************************************************************************/

#include "matrix.hpp"

namespace ngbla
{

  /**
     The Cholesky-factorization of a symmetric dense matrix.
     A = L D L^T
  */
  template <class T>
  class FlatCholeskyFactors
  {
  protected:
    /// matrix size
    int n;
    /// left factor
    T * lfact;
    /// inverse diagonal
    T * diag;
  public:
    // typedef typename mat_traits<T>::TV_COL TV;
    /// Factor the matrix A
    FlatCholeskyFactors (const FlatMatrix<T> & a, T * data)
    {
      diag = data;
      Factor (a);
    }

    /// Factor the matrix A
    FlatCholeskyFactors (const FlatMatrix<T> & a, LocalHeap & lh)
    {
      diag = (T*)lh.Alloc(sizeof(T)*RequiredMem(a.Height()));
      Factor (a);
    }

    ///
    NGS_DLL_HEADER void Factor (const FlatMatrix<T> & a);
    /// Multiply with the inverse of A 
    template <typename TV1, typename TV2>
    // NGS_DLL_HEADER void Mult (SliceVector<TV> x, SliceVector<TV> y) const
    void Mult (TV1 && x, TV2 && y) const
    {   
      // TV sum, val;
      // decltype (y(0)) sum, val;


      const T *pj;

      for (int i = 0; i < n; i++)
        y(i) = x(i);
      
      for (int i = 0; i < n; i++)
        {
          auto sum = y(i);
          
          pj = PRow(i);
          for (int j = 0; j < i; ++j)
            sum -= pj[j] * y(j);
          
          y(i) = sum;
        }
      
      for (int i = 0; i < n; i++)
        {
          auto sum = diag[i] * y(i);
          y(i) = sum;
      }
      
      for (int i = n-1; i >= 0; i--)
        {
          pj = PRow(i);
          auto val = y(i);
          for (int j = 0; j < i; ++j)
            y(j) -= pj[j] * val;
        }
    }
    
    /// Print factorization
    NGS_DLL_HEADER ostream & Print (ostream & ost) const;


    /// computes required memory
    static int RequiredMem (int n)
    { return n*(n+1)/2; }

  private:
    /// first element in row
    T * PRow (int i) const { return lfact + (i*(i-1)) / 2; }
  };


  ///  output operator.
  template<typename T>
  inline std::ostream & operator<< (std::ostream & s, const FlatCholeskyFactors<T> & m)
  {
    m.Print (s);
    return s;
  }



  template <class T>
  class CholeskyFactors : public FlatCholeskyFactors<T>
  {
  public:
    /// Factor the matrix A
    CholeskyFactors (const FlatMatrix<T> & a)
      : FlatCholeskyFactors<T> (a, new T[this->RequiredMem(a.Height())])
    { ; }
    /// Delete memory
    ~CholeskyFactors ()
    {
      delete [] this->diag;
    }
  };








  // high performance LDL factorization as used for SparseCholesky
  /*
  template <typename T, ORDERING ORD>
  INLINE void MySubABt (SliceMatrix<T,ORD> a,
                        SliceMatrix<T,ORD> b,
                        SliceMatrix<T,ORD> c)
  {
    static Timer timer1("SparseCholesky::Factor gemm 1", NoTracing);
    static Timer timer2("SparseCholesky::Factor gemm 2", NoTracing);
    static Timer timer3("SparseCholesky::Factor gemm 3", NoTracing);
            
    // if (c.Height() < 10 && c.Width() < 10) //  && a.Width() < 10)
    if (c.Height() < 10 || c.Width() < 10 || a.Width() < 10)
    // if (false)
      {
        // timer1.Start();
        c -= a * Trans(b);
        // timer1.Stop();
        // timer1.AddFlops(c.Height()*c.Width()*a.Width());
      }
    else
      {
        if (c.Height() < 128 && c.Width() < 128)
          // if (true)
          {
            // timer2.Start();
            // c -= a * Trans(b) | Lapack;
            ngbla::SubABt(a,b,c);
            // timer2.Stop();
            // timer2.AddFlops(c.Height()*c.Width()*a.Width());
          }
        else
          {
            timer3.Start();
            int nr = c.Height()/128+1;
            int nc = c.Width()/128+1;
            task_manager -> CreateJob
                         ( [&] (const TaskInfo & ti)
                           {
                             int br = ti.task_nr % nr;
                             int bc = ti.task_nr / nr;
                             auto rowr = Range(c.Height()).Split (br, nr);
                             auto colr = Range(c.Width()).Split (bc, nc);
                             // c.Rows(rowr).Cols(colr) -= a.Rows(rowr) * Trans(b.Rows(colr)) | Lapack;
                             ngbla::SubABt(a.Rows(rowr),b.Rows(colr), c.Rows(rowr).Cols(colr));
                           }, nr*nc);
            timer3.AddFlops(c.Height()*c.Width()*a.Width());
            timer3.Stop();
          }
      }
  }
  */

  
/*
  A   B^t     =  L1        D1  0       L1^t  B1 
  B   C       =  B1 L2      0  D2            L2^t
 */ 

  /*
  template <typename T, ORDERING ORD>
  void MySubADBt (SliceMatrix<T,ORD> A,
                  SliceVector<T> diag,
                  SliceMatrix<T,ORD> B,
                  SliceMatrix<T,ORD> C)
  {
    Matrix<T,ORD> hB(B.Height(), B.Width());
    for (int i = 0; i < hB.Width(); i++)
      hB.Col(i) = diag(i) * B.Col(i);
    C -= A * Trans(hB);
  }
  */

  
  template <typename T, ORDERING ORD>
  INLINE void MySubADBt (SliceMatrix<T,ORD> a,
                         SliceVector<T> diag,
                         SliceMatrix<T,ORD> b,
                         SliceMatrix<T,ORD> c,
                         bool symmetric)
  {
    // static Timer timer1("SparseCholesky::Factor gemm 1", NoTracing);
    // static Timer timer2("SparseCholesky::Factor gemm 2", NoTracing);
    // static Timer timer3("SparseCholesky::Factor gemm 3", NoTracing);

    /*
    // if (c.Height() < 10 && c.Width() < 10) //  && a.Width() < 10)
    if (c.Height() < 10 && c.Width() < 10 && a.Width() < 10)
    // if (false)
      {
        // timer1.Start();
        T hmem[100];
        FlatMatrix<T,ORD> hb(b.Height(), b.Width(), &hmem[0]);
        for (int i = 0; i < hb.Width(); i++)
          hb.Col(i) = diag(i) * b.Col(i);
        c -= a * Trans(hb);
        // timer1.Stop();
        // timer1.AddFlops(c.Height()*c.Width()*a.Width());
      }
    else
    */
      {
        if ( (c.Height() < 128 && c.Width() < 128) ||
             (size_t(c.Height())*c.Width()*a.Width() < 10000) )
          // if (true)
          {
            // timer2.Start();
            ngbla::SubADBt(a,diag,b,c);
            // timer2.Stop();
            // timer2.AddFlops(size_t(c.Height())*c.Width()*a.Width());
          }
        else
          {
            // timer3.Start();
            // int nr = c.Height()/128+1;
            // int nc = c.Width()/128+1;
            
            // constexpr int BH = 96;
            // constexpr int BW = 128;
            // avoid warning for capturing constexpr int, but still needed on MSVC 19.16
            IC<96> BH;
            IC<128> BW;
            int nr = (c.Height()+BH-1) / BH;
            int nc = (c.Width()+BW-1) / BW;
            task_manager -> CreateJob
                           ( [a,b,c,diag,nr,symmetric,BH,BW] (const TaskInfo & ti)
                           {
                             size_t br = ti.task_nr % nr;
                             size_t bc = ti.task_nr / nr;
                             // auto rowr = Range(c.Height()).Split (br, nr);
                             // auto colr = Range(c.Width()).Split (bc, nc);
                             auto rowr = Range(BH*br, min(BH*(br+1), c.Height()));
                             auto colr = Range(BW*bc, min(BW*(bc+1), c.Width()));
                             if (symmetric)
                               if (rowr.Next() <= colr.First())
                                 return; // need only lower half
                             
                             // c.Rows(rowr).Cols(colr) -= a.Rows(rowr) * Trans(b.Rows(colr)) | Lapack;
                             ngbla::SubADBt(a.Rows(rowr),diag, b.Rows(colr), c.Rows(rowr).Cols(colr));
                           }, nr*nc);
            // timer3.AddFlops(size_t(c.Height())*c.Width()*a.Width());
            // timer3.Stop();
          }
      }
  }

  template <typename T, ORDERING ORD>
  INLINE void MySubADBh (SliceMatrix<T,ORD> a,
                         SliceVector<T> diag,
                         SliceMatrix<T,ORD> b,
                         SliceMatrix<T,ORD> c,
                         bool symmetric)
  {
    Matrix<T,ORD> bconj(b.Height(), b.Width());
    bconj = Conj(b);
    MySubADBt (a, diag, SliceMatrix<T,ORD>(bconj), c, false); // symmetric);
  }


  
  
// Solve for B1:   B1 D1 L1^t = B
  template <typename T, ORDERING ORD>
  void CalcLDL_SolveL (SliceMatrix<T,ORD> L, SliceMatrix<T,ORD> B)
  {
    size_t n = L.Height();
    if (n == 1) return;

    if (n >= 2)
      {
        IntRange r1(0,n/2), r2(n/2,n);
        auto L1 = L.Rows(r1).Cols(r1);
        auto L21 = L.Rows(r2).Cols(r1);
        auto L2 = L.Rows(r2).Cols(r2);
        auto B1 = B.Cols(r1);
        auto B2 = B.Cols(r2);
        
        CalcLDL_SolveL(L1, B1);
        MySubADBt (B1, L1.Diag(), L21, B2, false);
        CalcLDL_SolveL(L2, B2);
        return;
      }
    
    static Timer t("LDL - Solve L work", NoTracing);
    t.Start();
    /*
      for (int i = 0; i < L.Height(); i++)
      for (int j = i+1; j < L.Height(); j++)
      for (int k = 0; k < B.Height(); k++)
      B(k,j) -= L(j,i) * B(k,i);
      // B.Col(j) -= L(j,i) * B.Col(i);
      */
    /*
      for (int k = 0; k < B.Height(); k++)
      for (int i = 0; i < L.Height(); i++)
      for (int j = i+1; j < L.Height(); j++)
      B(k,j) -= L(j,i) * B(k,i);
    */
    auto solve_row = [&] (size_t k)
      {
        auto Brow = B.Row(k);
        for (size_t i = 0; i < L.Height(); i++)
          for (size_t j = i+1; j < L.Height(); j++)
            Brow(j) -= L(j,i) * Brow(i);
      };
    if (B.Height() < 1000)
      for (size_t k = 0; k < B.Height(); k++)
        solve_row(k);
    else
      ParallelFor (B.Height(), solve_row);
    
    t.Stop();
  }
  
  // calc new A22-block
  // A2 -= B D B^t
  template <typename T, ORDERING ORD>
  void CalcLDL_A2 (SliceVector<T> diag, SliceMatrix<T,ORD> B, SliceMatrix<T,ORD> A2)
  {
    MySubADBt (B, diag, B, A2, true);
  }

  // calc new A22-block hermitsch
  // A2 -= B D B^h
  template <typename T, ORDERING ORD>
  void CalcLDL_A2H (SliceVector<T> diag, SliceMatrix<T,ORD> B, SliceMatrix<T,ORD> A2)
  {
    MySubADBh (B, diag, B, A2, true);
  }

  

// Calc A = L D L^t
  template <typename T, ORDERING ORD>
  void CalcLDL (SliceMatrix<T,ORD> mat)
  {
    size_t n = mat.Height();
    
    if (n >= 2)
      {
        size_t n1 = n/2;
        auto L1 = mat.Rows(0,n1).Cols(0,n1);
        auto L2 = mat.Rows(n1,n).Cols(n1,n);
        auto B = mat.Rows(n1,n).Cols(0,n1);
        CalcLDL (L1);
        CalcLDL_SolveL (L1,B);
        CalcLDL_A2 (L1.Diag(),B,L2);
        CalcLDL (L2);
        return;
      }

    if (n == 1)
      {
        // auto hm = mat(0,0);
        // CalcInverse (hm, mat(0,0));
        mat(0,0) = Inv(mat(0,0));
        return;
      }

    /*
      not working anymore
    for (size_t i = 0; i < n; i++)
      {
        T dii = mat(i,i);
        T inv_dii;
        CalcInverse (dii, inv_dii);
        for (size_t j = i+1; j < n; j++)
          {
            T hji = mat(j,i);
            T hjiD = hji * inv_dii;
            mat(j,i) = hjiD;
            for (size_t k = i+1; k <= j; k++)
              mat(j,k) -= hji * Trans(mat(k,i));
          }
      }
    */
  }


// Calc A = L D L^h
  template <typename T, ORDERING ORD>
  void CalcLDLH (SliceMatrix<T,ORD> mat)
  {
    size_t n = mat.Height();
    
    if (n >= 2)
      {
        size_t n1 = n/2;
        auto L1 = mat.Rows(0,n1).Cols(0,n1);
        auto L2 = mat.Rows(n1,n).Cols(n1,n);
        auto B = mat.Rows(n1,n).Cols(0,n1);
        CalcLDLH (L1);
        CalcLDL_SolveL (L1,B);
        CalcLDL_A2H (L1.Diag(),B,L2);
        CalcLDLH (L2);
        return;
      }

    if (n == 1)
      {
        // auto hm = mat(0,0);
        // CalcInverse (hm, mat(0,0));
        mat(0,0) = Inv(mat(0,0));
        return;
      }
  }


  
  template <typename T, ORDERING ORD>
  void SolveLDL (SliceMatrix<T,ORD> mat, FlatVector<T> sol)
  {
    size_t n = mat.Height();
    
    for (size_t i = 0; i < n; i++)
      {
        T tmp = mat(i,i)*sol(i);
        for (size_t j = i+1; j < n; j++)
          sol(j) -= mat(j,i) * tmp;
      }
    
    for (size_t i = 0; i < n; i++)
      sol(i) *= mat(i,i);
    
    for (size_t i = n; i--> 0; )
      {
        T hsum{0};
        for (size_t j = i+1; j < n; j++)
          hsum += mat(j,i)*sol(j);
        sol(i) -= mat(i,i) * hsum;
      }
  }

  template <typename T, ORDERING ORD>
  void SolveLDLH (SliceMatrix<T,ORD> mat, FlatVector<T> sol)
  {
    size_t n = mat.Height();
    
    for (size_t i = 0; i < n; i++)
      {
        T tmp = mat(i,i)*sol(i);
        for (size_t j = i+1; j < n; j++)
          sol(j) -= mat(j,i) * tmp;
      }
    
    for (size_t i = 0; i < n; i++)
      sol(i) *= mat(i,i);
    
    for (size_t i = n; i--> 0; )
      {
        T hsum{0};
        for (size_t j = i+1; j < n; j++)
          hsum += Conj(mat(j,i))*sol(j);
        sol(i) -= mat(i,i) * hsum;
      }
  }





  
  // invert lower left matrix
  template <typename T, ORDERING ORD>
  void CalcInverseL (FlatMatrix<T,ORD> mat)
  {
    // M = L^{-1}
    // i>j:   (M L)(i,j) = M_ik L_kj = 0   j <= k <= i
    // M_ij L_jj + M_ik L_kj  = 0      j < k <= i
    // M_ij = -1/{L_jj}  sum j < k <= i:  M_ik L_kj
    
    int n = mat.Height();
    STACK_ARRAY(T, mem, n);
    FlatVector<T> dinv(n, &mem);
    
    for (size_t i = 0; i < n; i++)
      // CalcInverse (mat(i,i), dinv(i));
      dinv(i) = Inv(mat(i,i));
    
    for (int i = n-1; i >= 0; i--)
      {
        mat(i,i) = dinv(i);
        for (int j = i-1; j >= 0; j--)
          {
            T sum(0.0);
            for (int k = j+1; k <= i; k++)
              sum += mat(i,k) * mat(k,j);
            mat(i,j) = -dinv(j) * sum;
          }
      }
  }
  
  // calculate inverse from LDL factorization
  template <typename T, ORDERING ORD>
  void CalcInverseLDL (FlatMatrix<T,ORD> mat)
  {
    size_t n = mat.Height();
    STACK_ARRAY(T,mem, n);
    FlatVector<T> dinv(n, &mem);
    
    for (size_t i = 0; i < n; i++)
      {
        // CalcInverse (mat(i,i), dinv(i));
        dinv(i) = Inv(mat(i,i));
        mat(i,i) = dinv(i);
      }
    
    CalcInverseL (mat);
    
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j <= i; j++)
        {
          T sum = mat(i,j);
          for (int k = i+1; k < n; k++)
            sum += Trans(mat(k,i)) * dinv(k) * mat(k,j);
          mat(i,j) = sum;
        }
    for (size_t i = 0; i < n; i++)
      for (size_t j = i+1; j < n; j++)
        mat(i,j) = mat(j,i);
  }
  




  
  // Solve for B1:   B1 D1 L1^t = B
  template <typename T, ORDERING ORD>
  void CalcLDLNew_SolveL (SliceMatrix<T,ORD> L, SliceMatrix<T,ORD> B)
  {
    size_t n = L.Height();
    if (n <= 1) return;

    IntRange r1(0,n/2), r2(n/2,n);
    auto L1 = L.Rows(r1).Cols(r1);
    auto L21 = L.Rows(r2).Cols(r1);
    auto L2 = L.Rows(r2).Cols(r2);
    auto B1 = B.Cols(r1);
    auto B2 = B.Cols(r2);
    
    CalcLDLNew_SolveL(L1, B1);
    // MySubADBt (B1, L1.Diag(), L21, B2, false);
    B2 -= B1 * Trans(L21);
    CalcLDLNew_SolveL(L2, B2);
  }

  
  // Calc A = L D^{-1} L^t
  template <typename T, ORDERING ORD>
  void CalcLDLNew (SliceMatrix<T,ORD> mat)
  {
    size_t n = mat.Height();
    if (n == 0) return;
    if (n == 1)
      {
        mat(0,0) = Inv(mat(0,0));
        // CalcInverse (mat(0,0));
        return;
      }      
    
    size_t n1 = n/2;
    auto L1 = mat.Rows(0,n1).Cols(0,n1);
    auto L2 = mat.Rows(n1,n).Cols(n1,n);
    auto B = mat.Rows(n1,n).Cols(0,n1);
    CalcLDLNew (L1);
    CalcLDLNew_SolveL (L1,B);
    auto diag = L1.Diag();
    MySubADBt (B, diag, B, L2, true);
    /*
    for (int i = 0; i < B.Height(); i++)
      for (int j = 0; j < B.Height(); j++)
        for (int k = 0; k < B.Width(); k++)
          L2(i,j) -= diag(k)*B(i,k)*B(j,k);
    */
    CalcLDLNew (L2);
    ScaleCols(B, diag);
    /*
    for (int i = 0; i < B.Width(); i++)
      B.Col(i) *= diag(i);
    */
  }



  /*
  template <typename T, ORDERING ORD>
  void SolveLDLNew (SliceMatrix<T,ORD> mat, FlatVector<T> sol)
  {
    size_t n = mat.Height();
    
    for (size_t i = 0; i < n; i++)
      {
        T tmp = sol(i);
        for (size_t j = i+1; j < n; j++)
          sol(j) -= mat(j,i) * tmp;
      }
    
    for (size_t i = 0; i < n; i++)
      sol(i) *= mat(i,i);
    
    for (size_t i = n; i--> 0; )
      {
        T hsum{0};
        for (size_t j = i+1; j < n; j++)
          hsum += mat(j,i)*sol(j);
        sol(i) -= hsum;
      }
  }
  */

  template <typename T, ORDERING ORD>
  void SolveL (SliceMatrix<T,ORD> mat, FlatVector<T> sol)
  {
    size_t n = mat.Height();
    if (n <= 1) return;

    if (n < 32)
      {
        for (size_t i = 0; i < n; i++)
          {
            T tmp = sol(i);
            for (size_t j = i+1; j < n; j++)
              sol(j) -= mat(j,i) * tmp;
          }
        return;
      }

    IntRange r1(0,n/2), r2(n/2,n);
    auto L1 = mat.Rows(r1).Cols(r1);
    auto L21 = mat.Rows(r2).Cols(r1);
    auto L2 = mat.Rows(r2).Cols(r2);
    auto sol1 = sol.Range(r1);
    auto sol2 = sol.Range(r2);

    SolveL (L1, sol1);
    sol2 -= L21 * sol1;
    SolveL (L2, sol2);
  }

  template <typename T, ORDERING ORD>
  void SolveLT (SliceMatrix<T,ORD> mat, FlatVector<T> sol)
  {
    size_t n = mat.Height();

    if (n <= 1) return;

    if (n < 32)
      {
        for (size_t i = n; i--> 0; )
          {
            T hsum{0};
            for (size_t j = i+1; j < n; j++)
              hsum += mat(j,i)*sol(j);
            sol(i) -= hsum;
          }
        return;
      }

    IntRange r1(0,n/2), r2(n/2,n);
    auto L1 = mat.Rows(r1).Cols(r1);
    auto L21 = mat.Rows(r2).Cols(r1);
    auto L2 = mat.Rows(r2).Cols(r2);
    auto sol1 = sol.Range(r1);
    auto sol2 = sol.Range(r2);

    SolveLT (L2, sol2);
    sol1 -= Trans(L21) * sol2;
    SolveLT (L1, sol1);
  }

  template <typename T, ORDERING ORD>
  void SolveLDLNew (SliceMatrix<T,ORD> mat, FlatVector<T> sol)
  {
    size_t n = mat.Height();
    
    SolveL (mat, sol);
    
    auto diag = mat.Diag();
    for (size_t i = 0; i < n; i++)
      sol(i) *= diag(i);
    
    SolveLT (mat, sol);    
  }

}

#endif
