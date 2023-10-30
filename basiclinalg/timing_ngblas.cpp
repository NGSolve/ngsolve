#define COMPILE_NGBLAS

#include <bla.hpp>


namespace ngbla
{

#include "matkernel.hpp"

  

  /**************** timings *********************** */

  template <int SW>
  inline void CopyMatrixIn (size_t h, size_t w,
                     double * ps, size_t dists,
                     SIMD<double,SW> * pd, size_t distd) { ; } 

  template <OPERATION OP, size_t MAXHA>
  inline void REGCALL MultAtB_intern2 (size_t ha, size_t wa, size_t wb,
                                BareSliceMatrix<double> a, BareSliceMatrix<double> b, BareSliceMatrix<double> c) { ; } 


  
  extern void MultUL (SliceMatrix<> A);
  extern void LapackSVD (SliceMatrix<> A,
                         SliceMatrix<double, ColMajor> U,
                         SliceMatrix<double, ColMajor> V);


  string GetTimingHelpString()
  {
    return string(R"raw_string(
Available options timings are:
          -1 .. this help
          0 ... run all timings
          1 ... A = B,   A,B = n*m,   A = aligned, fixed dist
          2 ... A = 0,   A = n*m,     but sliced
          3 ... A = B^t, A = n*m, 
          5 ... y = A*x,   A = n*m
          6 ... y = A^t*x,   A = n*m
          7 ... y += A^t*x(ind),   A = n*m
          10 .. C = A * B,   A=n*m, B=m*k, C=n*k
          11 .. C += A * B,   A=n*m, B=m*k, C=n*k
          // "20 .. C = A * B    A=n*m, B=n*k', C=n*k', k'=round(k), B aligned
          20 .. X = T * X       T=n*n triangular, X=n*m "
          21 .. X = T^-1 * X     T=n*n triangular, X=n*m "
          22 .. T^-1             T=n*n triangular"
          50 .. C += A * B^t,   A=n*k, B=m*k, C=n*m
          51 .. C += A * B^t,   A=n*k, B=m*k, C=n*m,  A,B aligned
          52 .. C = A * B^t,   A=n*k, B=m*k, C=n*m
          60 .. C -= A^t * D B,  A=n*k, B=n*m, C = k*m, D=diag
          61 .. C = A^t B,  A=n*k, B=n*m, C = k*m
          70 .. C += A B^t,  A=n*k, B=m*k, C = n*m, A,B SIMD
	  80 .. (x,y)        inner product, size n
          100.. MultAddKernel  C += A * B,  A=4*n, B=n*3SW
          101.. MultAddKernel  C += A * B,  A=4*n, B=n*3SW, B aligned
          110.. MultAddKernel2  C += A * B,  A=4*n, B=n*m, m multiple of 3*SW
          111.. MultAddKernel2  C += A * B,  A=4*n, B=n*m, m multiple of 3*SW, B aligned
          150.. ScalKernel     C = A * B^t,  A=4*n, B = 3*n
          151.. ScalKernel     C = A * B^t,  A=4*n, B = 3*n\n, A,B aligned
          200.. CalcInverse        A = nxn
          201.. CalcInverse by LU  A = nxn          
          205.. LDL                A = nxn
          210.. CalcInverseLapack  A = nxn
          300.. CalcSVD            A = nxn
)raw_string");
  }


  list<tuple<string,double>> Timing (int what, size_t n, size_t m, size_t k,
                                     bool lapack, bool doubleprec, size_t maxits)
  {
    if (what < 0)
      {

        cout << GetTimingHelpString() << flush;
        /*
        cout << "Available options timings are:\n"
          "-1 .. this help\n"
          "0 ... run all timings\n"
          "1 ... A = B,   A,B = n*m,   A = aligned, fixed dist\n"
          "2 ... A = 0,   A = n*m,     but sliced\n"
          "3 ... A = B^t, A = n*m, \n"
          "5 ... y = A*x,   A = n*m\n"
          "6 ... y = A^t*x,   A = n*m\n"
          "7 ... y += A^t*x(ind),   A = n*m\n"
          "10 .. C = A * B,   A=n*m, B=m*k, C=n*k\n"
          "11 .. C += A * B,   A=n*m, B=m*k, C=n*k\n"
          // "20 .. C = A * B    A=n*m, B=n*k', C=n*k', k'=round(k), B aligned\n"
          "20 .. X = T * X       T=n*n triangular, X=n*m "
          "21 .. X = T^-1 * X     T=n*n triangular, X=n*m "
          "22 .. T^-1             T=n*n triangular"
          "50 .. C += A * B^t,   A=n*k, B=m*k, C=n*m\n"
          "51 .. C += A * B^t,   A=n*k, B=m*k, C=n*m,  A,B aligned\n"
          "52 .. C = A * B^t,   A=n*k, B=m*k, C=n*m\n"
          "53 .. C += A * B^t,   A=n*k, B=m*k, C=n*m\n, complex B"          
          "54 .. C += A * B^t,   A=n*k, B=m*k, C=n*m\n, complex A,B"          
          "60 .. C -= A^t * D B,  A=n*k, B=n*m, C = k*m, D=diag\n"
          "61 .. C = A^t B,  A=n*k, B=n*m, C = k*m\n"
          "70 .. C += A B^t,  A=n*k, B=m*k, C = n*m, A,B SIMD\n"
	  "80 .. (x,y)        inner product, size n\n"
          "100.. MultAddKernel  C += A * B,  A=4*n, B=n*3SW\n"
          "101.. MultAddKernel  C += A * B,  A=4*n, B=n*3SW, B aligned\n"
          "110.. MultAddKernel2  C += A * B,  A=4*n, B=n*m, m multiple of 3*SW\n"
          "111.. MultAddKernel2  C += A * B,  A=4*n, B=n*m, m multiple of 3*SW, B aligned\n"
          "150.. ScalKernel     C = A * B^t,  A=4*n, B = 3*n\n"
          "151.. ScalKernel     C = A * B^t,  A=4*n, B = 3*n\n, A,B aligned\n"
          "200.. CalcInverse        A = nxn\n"
          "201.. CalcInverse by LU  A = nxn\n"          
          "205.. LDL                A = nxn\n"
          "210.. CalcInverseLapack  A = nxn\n"
          "300.. CalcSVD            A = nxn\n"
             << endl;
        */
        return list<tuple<string,double>>();
      }

    list<tuple<string,double>> timings;
    constexpr int SW = SIMD<double>::Size();
    if (what == 0 || what == 1)
      {
        // A = B
        constexpr size_t WA = 128;
        if (m > WA)
          {
            m = WA;
            cout << "max width = " << WA << endl;
          }
        Matrix<> b(n,m);
        STACK_ARRAY(SIMD<double>, mema, n*WA/SIMD<double>::Size());
        FlatMatrix<SIMD<double>> a(n,WA/SIMD<double>::Size(),&mema[0]);
        b = 1;
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("Copy matrix, packed dest");
          t.Start();
          for (size_t j = 0; j < its; j++)
            CopyMatrixIn(n,m, &b(0,0), m, &a(0,0), a.Width());
          t.Stop();
          cout << "Lapack GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Copy matrix, packed dest", 1e-9 * n*m*its / t.GetTime()));
        }
      }


    if (what == 0 || what == 2)
      {
        // A = 0
        Matrix<> a(n,m);
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("Zero matrix, packed dest");
          t.Start();
          for (size_t j = 0; j < its; j++)
            a.Rows(0,n).Cols(0,m) = j;
          t.Stop();
          cout << "Zero matrix GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Zero matrix", 1e-9 * n*m*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 3)
      {
        // A = B^t
        Matrix<> a(n,m), b(m,n);
        b = 1;
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("Matrix Transpose");
          t.Start();
          for (size_t j = 0; j < its; j++)
            TransposeMatrix(b, a);
          t.Stop();
          cout << "Lapack GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Transpose matrix", 1e-9 * tot*its / t.GetTime()));
        }
      }


    
    if (what == 0 || what == 5)
      {
        // y = A*x
        Matrix<> a(n,m);
        Vector<> x(m), y(n);
        a = 1; x = 2;
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          MultMatVec(a,x,y);
          if (L2Norm(a*x-y) > 1e-8)
            throw Exception("MultMatVec is faulty");
          Timer t("y = A*x");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MultMatVec(a,x,y);
          t.Stop();
          cout << "MultMatVec GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatVec", 1e-9 * n*m*its / t.GetTime()));
        }
    #ifdef LAPACK
        {
          Timer t("y = A*x, Lapack");
          t.Start();
          for (size_t j = 0; j < its; j++)
            LapackMultAx (a, x, y);
          t.Stop();
          cout << "MultMatVec Lapack GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatVecLapack", 1e-9 * n*m*its / t.GetTime()));
        }
    #endif // LAPACK
      }

    if (what == 0 || what == 6)
      {
        // y = A*x
        Matrix<> a(n,m);
        Vector<> x(n), y(m);
        a = 1; x = 2;
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("y = A*x");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MultMatTransVec(a,x,y);
          t.Stop();
          cout << "MultMatTransVec GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatVec", 1e-9 * n*m*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 7)
      {
        // y = A*x
        Matrix<> a(n,m);
        Vector<> x(1000), y(m);
        Array<int> index(n);
        for (size_t i = 0; i < n; i++)
          index[i] = (17*i)%1000;
        a = 1; x = 2; y = 0;
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("y = A*x");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MultAddMatTransVecIndirect(1, a,x,y, index);
          t.Stop();
          cout << "MultAddMatTransVecIndirect GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultAddMatVecIndirect", 1e-9 * n*m*its / t.GetTime()));
        }
      }

    
    if (what == 0 || what == 10)
      {
        // C=A*B
        auto doit = [&] (auto a, auto b, auto c) {
          a = 1; b = 2;
          for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
              a(i,j) = sin(i+1) * cos(j);
          for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < k; j++)
              b(i,j) = cos(i+3) * cos(j);
          
          double tot = n*m*k;
          size_t its = 1e10 / tot + 1;
          // MultMatMat(a,b,c);
          if (tot < 1e6)
            {
              c = a * b;
              double err = L2Norm(a*b-c);
              if (err > 1e-8)
                throw Exception("MultMatMat is faulty");
          }
        
          {
            Timer t("C = A*B");
            t.Start();
            if (!lapack)
              for (size_t j = 0; j < its; j++)
                // MultMatMat(a,b,c);
                c = a*b;
            else
              for (size_t j = 0; j < its; j++)
                c = a*b | Lapack;
            t.Stop();
            cout << "MultMatMat GFlops = " << 1e-9 * n*m*k*its / t.GetTime() << endl;
            timings.push_back(make_tuple("MultMatMat", 1e-9 * n*m*k*its / t.GetTime()));
          }
        };

        if (doubleprec)
          {
            Matrix<double> a(n,m), b(m,k), c(n,k);
            doit (a,b,c);
          }
        else
          {
            Matrix<float> a(n,m), b(m,k), c(n,k);
            doit (a,b,c);
          }
        
      }
        

    if (what == 0 || what == 11)
      {
        // C=A*B
        Matrix<> a(n,m), b(m,k), c(n,k);
        a = 1; b = 2;
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < m; j++)
            a(i,j) = sin(i+1) * cos(j);
        for (size_t i = 0; i < m; i++)
          for (size_t j = 0; j < k; j++)
            b(i,j) = cos(i+3) * cos(j);
        c = 0.0;
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        // MultMatMat(a,b,c);
        {
          Timer t("C += A*B");
          t.Start();
          if (!lapack)
            for (size_t j = 0; j < its; j++)
              c += a*b;
          else
            for (size_t j = 0; j < its; j++)
              c += a*b | Lapack;
          t.Stop();
          cout << "Add AB GFlops = " << 1e-9 * n*m*k*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatMat", 1e-9 * n*m*k*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 20)
      {
        Matrix<> a(n,n), b(n,m);
        a = 1; b = 2;
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < n; j++)
            a(i,j) = sin(i+1) * cos(j);
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < m; j++)
            b(i,j) = cos(i+3) * cos(j);
        Matrix<> saveb = b;
        Matrix<double,ColMajor> at = a;
        double tot = n*n*m/2;
        size_t its = 1e10 / tot + 1;

        {
          b = saveb;
          TriangularMult<LowerLeft> (a, b);
          TriangularSolve<LowerLeft> (a, b);
          if (double err = L2Norm(b-saveb); err > 1e-7)
            {
              cout << "diff = " << Truncate(b-saveb) << endl;
              throw Exception("TriangularMult/Solve<LowerLeft> is buggy, err = "+ToString(err));
            }
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularMult<LowerLeft> (a, b);
            }
          t.Stop();
          cout << "TriangularMult<L> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularMult<L>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }
        {
          b = saveb;
          TriangularMult<UpperRight> (a, b);
          TriangularSolve<UpperRight> (a, b);
          if (double err = L2Norm(b-saveb); err > 1e-7)
            throw Exception("TriangularMult/Solve<UpperRight> is buggy, err = " + ToString(err));

          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularMult<UpperRight> (a, b);
            }
          t.Stop();
          cout << "TriangularMult<R> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularMult<R>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }

        {
          b = saveb;          
          TriangularMult<LowerLeft,Normalized> (a, b);
          TriangularSolve<LowerLeft,Normalized> (a, b);
          if (double err = L2Norm(b-saveb); err > 1e-7)
            throw Exception("TriangularMult/Solve<LowerLeft,Normalized> is buggy, err = "+ToString(err));

          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularMult<LowerLeft,Normalized> (a, b);
            }
          t.Stop();
          cout << "TriangularMult<L,N> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularMult<L,N>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }
        {
          b = saveb;          
          TriangularMult<UpperRight,Normalized> (a, b);
          TriangularSolve<UpperRight,Normalized> (a, b);
          if (double err = L2Norm(b-saveb); err > 1e-7)
            throw Exception("TriangularMult/Solve<UpperRight,Normalized> is buggy, err = "+ToString(err));

          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularMult<UpperRight,Normalized> (a, b);
            }
          t.Stop();
          cout << "TriangularMult<R,N> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularMult<R,N>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }
        {
          Timer t("X = Lt * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularMult<LowerLeft> (at, b);
            }
          t.Stop();
          cout << "TriangularMult<Lt> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularMult<Lt>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }
        {
          Timer t("X = Rt * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularMult<UpperRight> (at, b);
            }
          t.Stop();
          cout << "TriangularMult<Rt> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularMult<Rt>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }

      }
    if (what == 0 || what == 21)
      {
        Matrix<> a(n,n), b(n,m);
        a = 1; b = 2;
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < n; j++)
            a(i,j) = sin(i+1) * cos(j);
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < m; j++)
            b(i,j) = cos(i+3) * cos(j);
        Matrix<> saveb = b;
        
        double tot = n*n*m/2;
        size_t its = 1e10 / tot + 1;
        // MultMatMat(a,b,c);
        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularSolve<LowerLeft> (a, b);
            }
          t.Stop();
          cout << "TriangularSolve<L> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularSolve<L>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }
        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularSolve<UpperRight> (a, b);
            }
          t.Stop();
          cout << "TriangularSolve<R> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularSolve<R>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }

        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularSolve<LowerLeft, Normalized> (a, b);
            }
          t.Stop();
          cout << "TriangularSolve<L,N> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularSolve<L,N>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }
        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularSolve<UpperRight, Normalized> (a, b);
            }
          t.Stop();
          cout << "TriangularSolve<R,N> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularSolve<R,N>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }


      }




    if (what == 0 || what == 22)
      {
        // T^{-1}
        Matrix<> a(n,n);
        a = 1; 
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < n; j++)
            a(i,j) = sin(i+1) * cos(j);
        Matrix<> savea = a;
        
        double tot = n*n*n/6;
        size_t its = 1e10 / tot + 1;
        if (its > maxits) its = maxits;
        {
          Timer t("L^-1");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              // a = savea;
              TriangularInvert<LowerLeft> (a);
              TriangularInvert<LowerLeft> (a);
            }
          t.Stop();
          cout << "TriangularInvert<L> GFlops = " << 1e-9 * 2*n*n*n/6*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularInvert<L>", 1e-9 * 2*n*n*n/6*its / t.GetTime()));
        }
        {
          Timer t("R^-1");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              // a = savea;
              TriangularInvert<UpperRight> (a);
              TriangularInvert<UpperRight> (a);
            }
          t.Stop();
          cout << "TriangularInvert<R> GFlops = " << 1e-9 * 2*n*n*n/6*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularInvert<R>", 1e-9 * 2*n*n*n/6*its / t.GetTime()));
        }
        {
          Timer t("L^-1");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              // a = savea;
              TriangularInvert<LowerLeft, Normalized> (a);
              TriangularInvert<LowerLeft, Normalized> (a);
            }
          t.Stop();
          cout << "TriangularInvert<LN> GFlops = " << 1e-9 * 2*n*n*n/6*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularInvert<LN>", 1e-9 * 2*n*n*n/6*its / t.GetTime()));
        }
        {
          Timer t("R^-1");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              // a = savea;
              TriangularInvert<UpperRight, Normalized> (a);
              TriangularInvert<UpperRight, Normalized> (a);
            }
          t.Stop();
          cout << "TriangularInvert<RN> GFlops = " << 1e-9 * 2*n*n*n/6*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularInvert<RN>", 1e-9 * 2*n*n*n/6*its / t.GetTime()));
        }

        {
          Timer t("UL");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              a = savea;
              MultUL (a);
            }
          t.Stop();
          cout << "MultUL GFlops = " << 1e-9 * n*n*n/3*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultUL", 1e-9 * n*n*n/3*its / t.GetTime()));
        }
      }





    
    
    if (what == 0 || what == 50)
      {
        // C=A*B^t
        Matrix<> a(n,k), b(m,k), c(n,m);
        a = 1; b = 2;
        c = 0.0;        
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            // AddABt(a,b,c);
            c += a * Trans(b);
          t.Stop();
          cout << "AddABt GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddABt", 1e-9 * tot *its / t.GetTime()));
        }
      }

    if (what == 0 || what == 51)
      {
        // C=A*B^t
        if (k % SW != 0)
          cout << "k should be a multiple of " << SW << endl;
        size_t ks = k/SW;
        Matrix<SIMD<double>> a(n,ks), b(m,ks);
        Matrix<> c(n,m);
        a = SIMD<double>(1); b = SIMD<double>(2);
        c = 0.0;
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            AddABt(SliceMatrix<double> (a.Height(), SW*a.Width(), SW*a.Width(), (double*)a.Data()),
                   SliceMatrix<double> (b.Height(), SW*b.Width(), SW*b.Width(), (double*)&b(0)),
                   // SliceMatrix<double> (AFlatMatrix<double>(b)),
                   c);
          t.Stop();
          cout << "AddABt GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddABt", 1e-9 * tot *its / t.GetTime()));
        }
      }

    if (what == 0 || what == 52)
      {
        // C=A*B^t
        Matrix<> a(n,k), b(m,k), c(n,m);
        a = 1; b = 2;
        c = 0.0;        
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          if (!lapack)
            for (size_t j = 0; j < its; j++)
              c = a * Trans(b);
          else
            for (size_t j = 0; j < its; j++)
              c = a * Trans(b) | Lapack;
          t.Stop();
          cout << "MultABt GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultABt", 1e-9 * tot *its / t.GetTime()));
        }
      }


    if (what == 0 || what == 53)
      {
        // C=A*B^t
        Matrix<> a(n,k);
        Matrix<Complex> b(m,k), c(n,m);
        a = 1; b = Complex(2.0);
        c = Complex(0.0);
        double tot = 2*n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            AddABt(a,b,c);
          // c += a * Trans(b);
          t.Stop();
          cout << "AddABt GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddABt", 1e-9 * tot *its / t.GetTime()));
        }
      }

    


    if (what == 0 || what == 54)
      {
        // C=A*B^t
        Matrix<Complex> a(n,k);
        Matrix<Complex> b(m,k), c(n,m);
        a = Complex(1); b = Complex(2.0);
        c = Complex(0.0);

        for (int i = 0; i < a.Height()*a.Width(); i++)
          a(i) = Complex(cos(i), sin(i));
        for (int i = 0; i < b.Height()*b.Width(); i++)
          b(i) = Complex(cos(2*i), sin(2*i));


        /*
        AddABt(a,b,c);
        cout << "c = " << endl << c << endl;
        // c += a*Trans(b);
        
        c = a*Trans(b);
        cout << "cex = " << endl << c << endl;
        // cout << "norm(c) = " << L2Norm(c) << endl;
        */
        AddABt(a,b,c);
        c -= a*Trans(b);        
        cout << "norm(c) = " << L2Norm(c) << endl;        
        
        double tot = 4*n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            AddABt(a,b,c);
            // c += a * Trans(b) | Lapack;
          t.Stop();
          cout << "AddABt GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddABt", 1e-9 * tot *its / t.GetTime()));
        }
      }

    
    if (what == 0 || what == 60)
      {
        // C=A*B^t
        Matrix<> a(n,k), b(n,m), c(k,m);
        Vector<> d(n);
        a = 1, b = 1, d = 2;
        c = 0.0;
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C -= A^t*D*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            SubAtDB(a, d, b, c);
          t.Stop();
          cout << "AddAtDB GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddAtDB", 1e-9 * tot *its / t.GetTime()));
        }
      }

    if (what == 0 || what == 61)
      {
        // C=A^t*B
        Matrix<> a(n,k), b(n,m), c(k,m);
        for (size_t i = 0; i < a.Height(); i++)
          for (size_t j = 0; j < a.Width(); j++)
            a(i,j) = sin(i+1) * cos(j);
        for (size_t i = 0; i < b.Height(); i++)
          for (size_t j = 0; j < b.Width(); j++)
            b(i,j) = cos(i+3) * cos(j);
        
        c = 0.0;
        // MultAtB (a,b,c);
        c = Trans(a)*b;
        if (n <= 1000)
          {
            double err = L2Norm(Trans(a)*b-c);
            if (err > 1e-8)
              throw Exception("MultAtB is faulty");
          }
        double tot = n*m*k;
        size_t its = 5e7 / tot + 1;
        {
          Timer t("C = A^t*B");
          t.Start();
          if (!lapack)
            for (size_t j = 0; j < its; j++)
              c = Trans(a)*b;
          else
            for (size_t j = 0; j < its; j++)
              c = Trans(a)*b | Lapack;
            
          // MultAtB(a, b, c);
          t.Stop();
          cout << "MultAtB GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultAtB", 1e-9 * tot *its / t.GetTime()));
        }
        /*
        {
          Timer t("C = A^t*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            c = Trans(a)*b;
          t.Stop();
          cout << "C=A^t*B GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("A^T*B", 1e-9 * tot *its / t.GetTime()));
        }
        */
        /*
        {
          Timer t("C = A^t*B, block block");
          constexpr size_t BS = 96;
          if (n >= BS && m >= BS && k >= BS)
            {
              tot = BS*BS*BS;
              its = 5e7 / tot+1;
              auto ba = a.Rows(BS).Cols(BS);
              // auto bb = b.Rows(BS).Cols(BS);
              auto bc = c.Rows(BS).Cols(BS);
              // Matrix ba(BS,BS);
              Matrix bb(BS,BS);
              t.Start();
              for (size_t j = 0; j < its; j++)
                MultAtB_intern2<SET,BS> (ba.Height(), ba.Width(), bb.Width(), ba, bb, bc);
              t.Stop();
              cout << "MultAtB - block GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
              timings.push_back(make_tuple("MultAtB - block", 1e-9 * tot *its / t.GetTime()));
            }
        }
        */
      }
    
    if (what == 0 || what == 70)
      {
        // C=A*B^t
        if (k % SW != 0)
          cout << "k should be a multiple of " << SW << endl;
        size_t ks = k/SW;
        Matrix<SIMD<double>> a(n,ks), b(m,ks);
        Matrix<> c(n,m);
        a = SIMD<double>(1); b = SIMD<double>(2);
        c = 0.0;
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C += A*Bt, sym");
          t.Start();
          for (size_t j = 0; j < its; j++)
            AddABtSym(a, b, c);
          t.Stop();
          cout << "AddABt, sym GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddABt", 1e-9 * tot *its / t.GetTime()));
        }
      }
    
    if (what == 0 || what == 80)
      {
	Vector x(n), y(n);
        double tot = n;
        size_t its = 1e8 / tot + 1;
	x = 1; y = 1;
        {
          Timer t("InnerProduct");
	  double sum = 0.0;
          t.Start();
          for (size_t j = 0; j < its; j++)
            sum += InnerProduct(x,y);
          t.Stop();
	  if (sum == 0) cout << "xx" << endl; // print only in exceptional case
          cout << "InnerProduct GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("InnerProduct", 1e-9 * tot *its / t.GetTime()));
        }
      }

    
    if (what == 0 || what == 100)
      {
        // C=A*B
        Matrix<> a(4,n), b(n,3*SW), c(4,3*SW);
        a = 1; b = 2; c = 0;
        double tot = n*4*3*SW;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MatKernelMultAB<4,3,ADD>(n,a.Data(), a.Width(), &b(0), b.Width(), &c(0), c.Width());
          t.Stop();
          cout << "MatKernelAddAB 3x4 = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB", 1e-9 * tot*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 101)
      {
        // C=A*B
        Matrix<> a(4,n), c(4,3*SW);
        Matrix<SIMD<double>> b(n, 3);
        a = 1; b = SIMD<double>(2); c = 0;
        double tot = n*4*3*SW;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MatKernelMultAB<4,3,ADD>(n,a.Data(), a.Width(), &b(0), b.Width(), &c(0), c.Width());
          t.Stop();
          cout << "MatKernelAddAB 3x4, algined GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB aligned", 1e-9 * tot*its / t.GetTime()));
        }
      }
    
    if (what == 0 || what == 102)
      {
        // C=A*B
        constexpr int W = 8;
        Matrix<> a(2,n), b(n,W*SW), c(2,W*SW);
        a = 1; b = 2; c = 0;
        double tot = a.Width()*a.Height()*b.Width();
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MatKernelMultAB<2,W,ADD>(n,a.Data(), a.Width(), &b(0), b.Width(), &c(0), c.Width());
            // MatKernel2AddAB<1,SET>(b.Height(), b.Width(), a.Data(), a.Dist(), &b(0), b.Dist(), &c(0), c.Dist());
          t.Stop();
          cout << "MatKernelAddAB 1x4 = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB", 1e-9 * tot*its / t.GetTime()));
        }
      }


    if (what == 0 || what == 110)
      {
        // C=A*B
        if (m % (3*SW) != 0)
          cout << "m should be a multiple of 3*SIMD::Size" << endl;
        Matrix<> a(4,n), b(n,m), c(4,m);
        a = 1; b = 2; c = 0;
        double tot = n*4*m;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            for (size_t i = 0; i+3*SW <= m; i += 3*SW)
              MatKernelMultAB<4,3,ADD>(n,a.Data(), a.Width(), &b(i), b.Width(), &c(i), c.Width());
          t.Stop();
          cout << "MatKernel2AddAB 3x4 = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB", 1e-9 * tot*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 111)
      {
        // C=A*B
        if (m % (3*SW) != 0)
          cout << "m should be a multiple of 3*SIMD::Size" << endl;
        Matrix<> a(4,n), c(4,m);
        Matrix<SIMD<double>> b(n, m/SW);
        a = 1; b = SIMD<double>(2); c = 0;
        double tot = n*4*m;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            for (size_t i = 0; i+3*SW <= m; i += 3*SW)            
              MatKernelMultAB<4,3,ADD>(n,a.Data(), a.Width(), &b(i/SW), b.Width(), &c(i), c.Width());
          t.Stop();
          cout << "MatKernel2AddAB 3x4, algined GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB aligned", 1e-9 * tot*its / t.GetTime()));
        }
      }




    if (what == 0 || what == 150)
      {
        // C=A*B
        Matrix<> a(4,n), b(4,n), c(3,4);
        a = 1; b = 2; c = 0;
        double tot = n*4*3;
        size_t its = 1e10 / tot + 1;
        SIMD<double,4> sum(0);
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              auto res = MatKernelScalAB<3,4>(n,a.Data(), a.Width(), &b(0), b.Width());
              sum += get<0>(res) + get<1>(res) + get<2>(res);
            }
          t.Stop();
          cout << sum;
          cout << "MatKernelScalAB 4x3 = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelScalAB 4x3", 1e-9 * tot*its / t.GetTime()));
        }
      }
    

    if (what == 0 || what == 151)
      {
        // C=A*B
        Matrix<SIMD<double>> a(4,n), b(4,n);
        Matrix<> c(3,4);
        a = SIMD<double>(1); b = SIMD<double>(2); c = 0;
        double tot = n*4*3*SW;
        size_t its = 1e10 / tot + 1;
        SIMD<double,4> sum(0);
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              auto res = MatKernelScalAB<3,4>(n,a.Data(), a.Width(), &b(0), b.Width());
              sum += get<0>(res) + get<1>(res) + get<2>(res);
            }
          t.Stop();
          cout << sum;
          cout << "MatKernelScalAB, simd 4x3 = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelScalAB, simd 4x3", 1e-9 * tot*its / t.GetTime()));
        }
      }
    


    if (what == 0 || what == 200)
      {
        // CalcInverse
        Matrix<> a(n,n);
        a = 1;
        a.Diag() = 10000;
        double tot = n*n*n;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("Inv(A)");
          t.Start();
          for (size_t j = 0; j < its; j++)
            CalcInverse(a, INVERSE_LIB::INV_NGBLA);
          t.Stop();
          cout << "Inv(A) GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Inv(A)", 1e-9 * tot *its / t.GetTime()));
        }
      }

    if (what == 0 || what == 201)
      {
        // CalcInverse
        Matrix<> a(n,n);
        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++)
            a(i,j) = cos(i+j);
        // a = 1;
        // a.Diag() = 1.1;
        double tot = n*n*n;
        size_t its = 1e10 / tot + 1;
        if (its > maxits) its = maxits;
        {
          Timer t("Inv(A)");
          t.Start();
          for (size_t j = 0; j < its; j++)
            CalcInverse(a, INVERSE_LIB::INV_NGBLA_LU);
          t.Stop();
          cout << "Inv(A) GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Inv(A)", 1e-9 * tot *its / t.GetTime()));
        }

        
        {
          Timer t("CalcLU");
          Array<int> p(n);
          Matrix<> ha = a;
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              ha = a;
              CalcLU(ha, p);
            }
          t.Stop();
          cout << "CalcLU GFlops = " << 1e-9 * tot/3*its / t.GetTime() << endl;
          timings.push_back(make_tuple("CalcLU", 1e-9 * tot/3 *its / t.GetTime()));
        }

        {
          Timer t("InvFromLU");
          Array<int> p(n);
          CalcLU(a, p);          
          Matrix<> ha = a;
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              ha = a;
              InverseFromLU(ha, p);
            }
          t.Stop();
          cout << "InvFromLU GFlops = " << 1e-9 * tot*2/3*its / t.GetTime() << endl;
          timings.push_back(make_tuple("InvFromLU", 1e-9 * tot*2/3 *its / t.GetTime()));
        }



      }

    

    if (what == 0 || what == 205)
      {
        // CalcInverse
        Matrix<double,ColMajor> a(n,n);
        a = 1;
        Trans(a).Diag() = 10000;
        double tot = n*n*n;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("Inv(A)");
          t.Start();
          for (size_t j = 0; j < its; j++)
            CalcLDL (SliceMatrix<double,ColMajor> (a));
          t.Stop();
          cout << "Inv(A) GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Inv(A)", 1e-9 * tot *its / t.GetTime()));
        }
      }

    
     if (what == 0 || what == 210)
      {
        // CalcInverse
        Matrix<> a(n,n);
        a = 1;
        a.Diag() = 10000;
        double tot = n*n*n;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("Inv(A)");
          t.Start();
          for (size_t j = 0; j < its; j++)
            LapackInverse(a);
          t.Stop();
          cout << "LapackInv(A) GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("LapackInv(A)", 1e-9 * tot *its / t.GetTime()));
        }
      }

  #ifdef LAPACK
     if (what == 0 || what == 211)
      {
        // CalcInverse
        Matrix<> a(n,n);
        a = 1;
        a.Diag() = 10000;
        double tot = n*n*n;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("Inv(A)");
          t.Start();
          for (size_t j = 0; j < its; j++)
            LapackInverseSPD(a);
          t.Stop();
          cout << "LapackInv(A) GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("LapackInv(A)", 1e-9 * tot *its / t.GetTime()));
        }
      }
  #endif // LAPACK
     
     if (what == 0 || what == 300)
      {
        // CalcSVD
        Matrix<> a(n,n);
        Matrix<double,ColMajor> U(n,n), V(n,n);
        for (int i = 0; i < a.Height(); i++)
          for (int j = 0; j < a.Width(); j++)
            a(i,j) = rand() / double(RAND_MAX);
        
        Matrix aorig = a;
        double tot = 5 * n*n*n;
        size_t its = 1e10 / tot + 1;
        its = min(maxits, its);
        if (!lapack)
          {
            Timer t("CalcSVD");
            t.Start();
            for (size_t j = 0; j < its; j++)
              {
                a = aorig;
                CalcSVD(a, U, V);
            }
            t.Stop();
            cout << "CalcSVD GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
            timings.push_back(make_tuple("CalcSVD", 1e-9 * tot *its / t.GetTime()));
          }
        else
          {
          #ifdef LAPACK
            Timer t("LapackSVD");
            t.Start();
            for (size_t j = 0; j < its; j++)
              {
                a = aorig;
                LapackSVD(a, U, V);
              }
            t.Stop();
            cout << "LapackSVD GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
            timings.push_back(make_tuple("LapackSVD", 1e-9 * tot *its / t.GetTime()));
          #endif // LAPACK
          }
      }

    
    return timings;
  }


}



