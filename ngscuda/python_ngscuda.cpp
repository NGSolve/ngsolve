#include <core/python_ngcore.hpp>

#include "../ngstd/python_ngstd.hpp"
#include "cuda_linalg.hpp"

// TODO: always use ngs_cuda?
using namespace ngbla;
using namespace ngla;
using namespace ngs_cuda;

namespace ngla {
  extern void InitApplyIntegrationPoints ();
}

PYBIND11_MODULE(ngscuda, m) {

  InitCUDA(1);
  InitCuLinalg();
  InitApplyIntegrationPoints();

  m.def("InitCuLinalg", &InitCuLinalg, "Initializing cublas and cusparse.");
  
  py::class_<UnifiedVector, BaseVector, shared_ptr<UnifiedVector>>
    (m, "UnifiedVector", "UnifiedVector for CUDA applications", py::multiple_inheritance())
    
    .def(py::init([] (int size)
                  { 
                    return make_shared<UnifiedVector>(size); 
                  }))
    .def(py::init([] (const BaseVector &vec) 
                  {
                    return make_shared<UnifiedVector>(vec);
                  }))
    .def(py::init([] (py::array_t<double> bvec)
                  {
                    auto vec = bvec.template unchecked<1>();
                    shared_ptr<UnifiedVector> uv = make_shared<UnifiedVector>(vec.size());
                    FlatVector<double> fv = uv->FVDouble();
                    for (size_t i = 0; i < vec.size(); i++)
                      {
                        fv(i) = vec(i);
                      }
                    return uv;
                  }))

    /*
      JS: should be inherited from BaseVector ??
  //  TODO: extend for splicing (define UnifiedVector.Range?)
  //    not that important. maybe delete
  .def("__getitem__", [] (UnifiedVector & self, int ind)
      {
        if (ind < 0)
          ind += self.Size();
        if (ind < 0 || ind >= self.Size())
          py::index_error();
        return py::cast(self[ind]);
      }, py::arg("ind"), "Return value at given position")
  .def("__setitem__", [] (UnifiedVector & self, int ind, double z)
      { 
        if (ind < 0)
          ind += self.Size();
        if (ind < 0 || ind >= self.Size())
          py::index_error();
        self[ind] = z; 
      })
    */
    

    /*
  .def("UpdateHost", [] (UnifiedVector &self) -> void
        { self.UpdateHost(); }) 
  .def("UpdateDevice", [] (UnifiedVector &self) -> void
        { self.UpdateDevice(); });
    */
    .def("UpdateHost", &UnifiedVector::UpdateHost)
    .def("UpdateDevice", &UnifiedVector::UpdateDevice)
    ;

  /* .def("__str__", [] (UnifiedVector & self) { return ToString<UnifiedVector>(self); } ) */
  /* .def("__repr__", [] (UnifiedVector & self) { return "unfiedvector"; } ) */

  /* .def("__len__", [] (UnifiedVector & self) {return self.Size(); }) */

  /* .def("Scale", [] (UnifiedVector & self, double d) */
  /*       { */
  /*         self.Scale(d); */
  /*       }); */

  /* TODO: */  
  /* .def("__add__", [] (shared_ptr<UnifiedVector> e1, shared_ptr<UnifiedVector> e2) */
  /*       { */
  /*         cerr << "aunivec, vec add" << endl; */
  /*         return make_shared<BaseVector>((*e1) + (*e2)); */
  /*       }) */
  /* .def("__neg__", [] (shared_ptr<UnifiedVector> e1) */
  /*       { */
  /*         cerr << "univec neg" << endl; */
  /*         return make_shared<BaseVector>(-(*e1)); */
  /*       }) */
  /* .def("__sub__", [] (shared_ptr<UnifiedVector> e1, shared_ptr<UnifiedVector> e2) */
  /*       { */
  /*         cerr << "univec, base sub" << endl; */
  /*         return (*e1) - (*e2); */
  /*       }) */
  /* .def("__mul__", [] (shared_ptr<UnifiedVector> x, shared_ptr<DevSparseMatrix> a) */
  /*       { */
  /*         cerr << "univec mat prod" << endl; */
  /*         return (*x) * (*a); */
  /*       }) */
  /* .def("__rmul__", [] (double x, const UnifiedVector& v) */
  /*       { */
  /*         cerr << "scal, mat prod" << endl; */
  /*         return x * v; */
  /*       }); */

  /* .def("Set") */
  /* .def("__add__", (UnifiedVector & v1, BaseVector & v2, py::object s) -> void */
  /*     { */
  /*       return v1 + v2; */
  /*     }) */
  /* .def("Scale", [](UnifiedVector & self, double scal) -> shared_ptr<UnifiedVector> */
  /*     { cerr << "scaling by " << scal << endl; return self.Scale(scal); }) */

  // TODO: probably useless... use ToGPU instead?
  //
  //   alternative: BaesVector::ToGPU
  /* .def("Assign", [] (UnifiedVector &self, BaseVector &v2) -> void */
  /*       { */
  /*         self = v2; */
  /*         return; */
  /*       }) */

  /* .def("Assign", [] (UnifiedVector &self, py::object d) -> void */
  /*       { */
  /*         self.SetScalar(py::extract<double>(d)()); */
  /*         return; */
  /*       }) */

  // TODO:
  /* m.def("InnerProdcut", [] (UnifiedVector &a, UnifiedVector &b)) */

  py::class_<DevMatrix, BaseMatrix, shared_ptr<DevMatrix>>
    (m, "DevBaseMatrix", "device matrix for CUDA applications");

  py::class_<DevSparseMatrix, DevMatrix, shared_ptr<DevSparseMatrix>>
    (m, "DevSparseMatrix", "DevSparseMatrix for CUDA applications")

    .def(py::init ( [] (SparseMatrix<double>& mat) -> shared_ptr<DevSparseMatrix>
          {
            return make_shared<DevSparseMatrix>(mat);
          }
          ));
    /* .def("Mult", [] (DevSparseMatrix &self, UnifiedVector &x, UnifiedVector &y) */
    /*     { */
    /*       self.Mult(x, y); */
    /*     }) */
    /* .def("Mult", [] (DevSparseMatrix &self, double s, UnifiedVector &x, UnifiedVector &y) */
    /*     { */
    /*       self.MultAdd(s, x, y); */
    /*     }); */

  // TODO:
  // make operators work instead of Add, Mult, ..
  py::class_<DevDMatrix, DevMatrix, shared_ptr<DevDMatrix>>
    (m, "DevDMatrix", "dense device matrix for CUDA applications")
    
    .def(py::init ( [] (const Matrix<>& mat)
          {
            return make_shared<DevDMatrix>(mat);
          }))
    .def("Add", [] (DevDMatrix& self, DevDMatrix& b)
          {
            self.Add(b);
          })
    /* .def("Mult", [] (DevDMatrix& self, DevDMatrix& b, DevDMatrix& c) */
    /*       { */
    /*         self.Mult(b, c); */
    /*       }) */
    .def("Scale", [] (DevDMatrix& self, double d)
          {
            self.Scale(d);
          })
    .def("SetZero", [] (DevDMatrix& self)
          {
            self.SetZero();
          });
  m.def("MatMult", [] (const DevDMatrix& mata, const DevDMatrix& matb)
          {
            shared_ptr<DevDMatrix> matptr = MatMult(mata, matb);
            cerr << *matptr << endl;
            return matptr;
          });
    
    /* .def(py::init ( [] (BaseMatrix & mat) -> shared_ptr<BaseMatrix> */

  m.def("CreateDevMatrix", [] (BaseMatrix &mat)
          {
            return CreateDevMatrix(mat);
          });

  m.def("CreateDevMatrix", [] (Matrix<> &mat)
          {
            return CreateDevMatrix(mat);
          });

  py::class_<DevJacobiPrecond, DevSparseMatrix, shared_ptr<DevJacobiPrecond>>
    (m, "DevJacobiPrecond", "Jacobi Preconditioner working on device");

  m.def("CreateDevSmoother", [] (SparseMatrix<double> &mat, shared_ptr<BitArray> ba)
          {
            shared_ptr<DevJacobiPrecond> res_ptr = make_shared<DevJacobiPrecond>(mat, ba);
            return res_ptr;
          },
          py::arg("mat"),
          py::arg("freedofs") = shared_ptr<BitArray>());


  
  
  
  py::class_<Matrix<Dev<double>>> (m, "DevMatrix")
    .def(py::init<FlatMatrix<double>>())
    .def("D2H", &Matrix<Dev<double>>::D2H)
    .def("__matmul__", [&](const Matrix<Dev<double>> & a, const Matrix<Dev<double>> & b)
         {
           Matrix<Dev<double>> c(a.Height(), b.Width());
           MultMatMat (a, b, c);
           return c;
         })
    
    .def("__timing__", []()
         {
           for (int n = 100; n < 5000; n *= 2)
             {
               cout << "n = " << n << endl;
               Matrix a(n,n), b(n,n);
               a = 1; b = 2;
               Matrix<Dev<double>> deva(a);
               Matrix<Dev<double>> devb(b);
               Matrix<Dev<double>> devc(n, n);
               int runs = 1e11 / (double(n)*double(n)*double(n)) + 1;
               Timer t("matmat");
               t.Start();
               cudaDeviceSynchronize();
               for (int i = 0; i < runs; i++)
                 MultMatMat(deva, devb, devc);
               cudaDeviceSynchronize();
               t.Stop();
               t.AddFlops (double(runs)*n*n*n);
               cout << "time = " << t.GetTime() << " MFlops = " << t.GetMFlops() << endl;
             }
         })
    ;
}

