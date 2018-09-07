#include "catch.hpp"
#include <bla.hpp>
using namespace ngbla;

void SetRandom (SliceMatrix<> mat)
{
  for (int i = 0; i < mat.Height(); i++)
    for (int j = 0; j < mat.Width(); j++)
      mat(i,j) = sin(2+3*i+5*j);
}

void SetRandom (SliceVector<> vec)
{
  for (int i = 0; i < vec.Size(); i++)
    vec(i) = sin(3+3*i);
}

void SetRandom (SliceMatrix<Complex> mat)
{
  for (int i = 0; i < mat.Height(); i++)
    for (int j = 0; j < mat.Width(); j++)
      mat(i,j) = Complex(4+ sin(3*i+5*j), cos(2*i-j) );
}

void SetRandom (SliceVector<Complex> vec)
{
  for (int i = 0; i < vec.Size(); i++)
    vec(i) = Complex(5+ sin(3*i), cos(2*i) );
}


TEST_CASE ("MatVec", "[ngblas]") {
    for (int n = 1; n < 20; n++) {
        SECTION ("n = "+to_string(n)) {
            for (int m = 1; m < 20; m++) {
                SECTION ("m = "+to_string(m)) {
                  Matrix<> a(n,m);
                  Vector<> x(m), y(n);
                  SetRandom(a);
                  SetRandom(x);
                  y = a*x;
                  double err = L2Norm (a*x-y);
                  CHECK(err < 1e-13);
                }
            }
        }
    }
}

TEST_CASE ("MatTransVec", "[ngblas]") {
    for (int n = 1; n < 20; n++) {
        SECTION ("n = "+to_string(n)) {
            for (int m = 1; m < 20; m++) {
                SECTION ("m = "+to_string(m)) {
                  Matrix<> a(n,m);
                  Vector<> x(n), y(m);
                  SetRandom(a);
                  SetRandom(x);
                  y = Trans(a)*x;
                  // MultMatTransVec(a, x, y);
                  double err = L2Norm (Trans(a)*x-y);
                  CHECK(err < 1e-13);
                }
            }
        }
    }
}






TEST_CASE ("MultMatMat", "[ngblas]") {
    for (int n = 1; n < 20; n++) {
        SECTION ("n = "+to_string(n)) {
            for (int m = 1; m < 20; m++) {
                SECTION ("m = "+to_string(m)) {
                    for (int k = 1; k < 20; k++) {
                        SECTION ("k = "+to_string(k)) {
                            Matrix<> a(n,m), b(m,k), c(n,k);
                            SetRandom(a);
                            SetRandom(b);
                            MultMatMat (a, b, c);
                            double err = L2Norm (a*b-c);
                            CHECK(err < 1e-13);
                        }
                    }
                }
            }
        }
    }
}


TEST_CASE ("MinusMultAB", "[ngblas]") {
    for (int n = 1; n < 20; n++) {
        SECTION ("n = "+to_string(n)) {
            for (int m = 1; m < 20; m++) {
                SECTION ("m = "+to_string(m)) {
                    for (int k = 1; k < 20; k++) {
                        SECTION ("k = "+to_string(k)) {
                            Matrix<> a(n,m), b(m,k), c(n,k);
                            SetRandom(a);
                            SetRandom(b);
                            MinusMultAB (a, b, c);
                            double err = L2Norm (a*b+c);
                            CHECK(err < 1e-13);
                        }
                    }
                }
            }
        }
    }
}


TEST_CASE ("AddAB", "[ngblas]") {
    for (int n = 1; n < 20; n++) {
        SECTION ("n = "+to_string(n)) {
            for (int m = 1; m < 20; m++) {
                SECTION ("m = "+to_string(m)) {
                    for (int k = 1; k < 20; k++) {
                        SECTION ("k = "+to_string(k)) {
                          Matrix<> a(n,m), b(m,k), c(n,k), c2(n,k);
                            SetRandom(a);
                            SetRandom(b);
                            SetRandom(c);
                            c2 = c;
                            AddAB (a, b, c);
                            double err = L2Norm (c2+a*b-c);
                            CHECK(err < 1e-13);
                        }
                    }
                }
            }
        }
    }
}


TEST_CASE ("SubAB", "[ngblas]") {
    for (int n = 1; n < 20; n++) {
        SECTION ("n = "+to_string(n)) {
            for (int m = 1; m < 20; m++) {
                SECTION ("m = "+to_string(m)) {
                    for (int k = 1; k < 20; k++) {
                        SECTION ("k = "+to_string(k)) {
                          Matrix<> a(n,m), b(m,k), c(n,k), c2(n,k);
                            SetRandom(a);
                            SetRandom(b);
                            SetRandom(c);
                            c2 = c;
                            SubAB (a, b, c);
                            double err = L2Norm (c2-a*b-c);
                            CHECK(err < 1e-13);
                        }
                    }
                }
            }
        }
    }
}




TEST_CASE ("SubAtDB", "[ngblas]") {
    for (int n = 1; n < 20; n++) {
        SECTION ("n = "+to_string(n)) {
            for (int m = 1; m < 20; m++) {
                SECTION ("m = "+to_string(m)) {
                    for (int k = 1; k < 20; k++) {
                        SECTION ("k = "+to_string(k)) {
                            Matrix<Complex> a(n,m), b(n,k), c(m,k);
                            Vector<Complex> diag(n);

                            SetRandom(a);
                            SetRandom(b);
                            SetRandom(c);
                            SetRandom(diag);

                            Matrix<Complex> c2 = c;

                            SubAtDB (a, diag, b, c);

                            for (int i = 0; i < diag.Size(); i++)
                                b.Row(i) *= diag(i);
                            c2 -= Trans(a) * b;

                            double err = L2Norm(c-c2);
                            CHECK(err < 1e-10);
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE ("SIMD<double>", "[simd]") {
    constexpr size_t N = SIMD<double>::Size();
    double src[N];
    double dst[N];
    for (auto i : Range(N)) {
        src[i] = i+1;
        dst[i] = 0.0;
    }

    SECTION ("Mask load") {
        for (auto k : Range(N+1)) {
            SIMD<double> simd(src,k);
            for (auto i : Range(N)) {
                CHECK(simd[i] == ( i<k? src[i] : 0.0 ));
            }
        }
    }

    SECTION ("Mask store") {
        for (auto k : Range(N+1)) {
            SIMD<double> simd(src);
            simd.Store(dst, k);
            for (auto i : Range(N)) {
                CHECK(dst[i] == ( i<k? src[i] : 0.0 ));
            }
        }
    }
}

TEST_CASE ("SIMD<Complex>", "[simd]") {
    constexpr size_t N = SIMD<Complex>::Size();
    Complex src[N];
    Complex dst[N];
    for (auto i : Range(N)) {
        src[i].real(i+1);
        src[i].imag(N*(i+1));
        dst[i] = 0.0;
    }

    SECTION ("Mask load/store") {
        for (auto k : Range(N+1)) {
            SIMD<Complex> simd;
            simd.LoadFast(src,k);
            simd.StoreFast(dst,k);
            for (auto i : Range(N)) {
                CHECK(dst[i].real() == ( i<k? src[i].real() : 0.0 ));
                CHECK(dst[i].imag() == ( i<k? src[i].imag() : 0.0 ));
            }
        }
    }
}

TEST_CASE ("Vec", "[double]") {
    Vec<1,double> v2{42};
    CHECK(v2[0] == 42);

    constexpr int N = 5;

    Vec<N, double> v0{42.42};
    for (int i : Range(N))
      CHECK(v0[i] == 42.42);

    double vals [5] = {1.2,3.4,5.6,7.8,9.1};
    Vec<N, double> v1{vals[0], vals[1], vals[2], vals[3], vals[4]};
    for (int i : Range(N))
      CHECK(v1[i] == vals[i]);
}
