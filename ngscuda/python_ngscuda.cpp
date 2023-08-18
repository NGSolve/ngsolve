#include <core/python_ngcore.hpp>

#include "../ngstd/python_ngstd.hpp"
#include "cuda_linalg.hpp"
#include "cuda_profiler.hpp"

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

    .def("UpdateHost", &UnifiedVector::UpdateHost)
    .def("UpdateDevice", &UnifiedVector::UpdateDevice)
    ;


  py::class_<DevMatrix, BaseMatrix, shared_ptr<DevMatrix>>
    (m, "DevBaseMatrix", "device matrix for CUDA applications");

  py::class_<DevSparseMatrix, DevMatrix, shared_ptr<DevSparseMatrix>>
    (m, "DevSparseMatrix", "DevSparseMatrix for CUDA applications")

    .def(py::init ( [] (SparseMatrix<double>& mat) -> shared_ptr<DevSparseMatrix>
          {
            return make_shared<DevSparseMatrix>(mat);
          }
          ));

#ifdef NONE
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
#endif   
    /* .def(py::init ( [] (BaseMatrix & mat) -> shared_ptr<BaseMatrix> */

  m.def("CreateDevMatrix", [] (BaseMatrix &mat)
          {
            return CreateDevMatrix(mat);
          });
/*
  m.def("CreateDevMatrix", [] (Matrix<> &mat)
          {
            return CreateDevMatrix(mat);
          });
*/
    
    /*
  py::class_<DevJacobiPrecond, DevSparseMatrix, shared_ptr<DevJacobiPrecond>>
    (m, "DevJacobiPrecond", "Jacobi Preconditioner working on device");

  m.def("CreateDevSmoother", [] (SparseMatrix<double> &mat, shared_ptr<BitArray> ba)
          {
            shared_ptr<DevJacobiPrecond> res_ptr = make_shared<DevJacobiPrecond>(mat, ba);
            return res_ptr;
          },
          py::arg("mat"),
          py::arg("freedofs") = shared_ptr<BitArray>());
*/

  m.def("__time_tracer__", TimeProfiler);
  m.def("SetCudaTimer", CudaRegionTimer::SetCudaTimer);
  
  
  
  py::class_<Matrix<Dev<double>>> (m, "DevMatrix")
    .def(py::init<FlatMatrix<double>>())
    .def("D2H", [&](const Matrix<Dev<double>> & a) 
         { return a.D2H(); }) 
    .def("__matmul__", [&](const Matrix<Dev<double>> & a, const Matrix<Dev<double>> & b)
         {
           Matrix<Dev<double>> c(a.Height(), b.Width());
           MultMatMat (a, b, c);
           return c;
         })
    
    .def("__timing__", []()
         {
            if (false)
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
             
            for (int n = 5; n <= 30; n++)
             {
            int nummats = 1000000;
            int math = n, matw = n;
            Array<double> hmatdata(nummats*math*matw);
            Array<double> hvecxdata(nummats*matw);
            Array<double> hvecydata(nummats*math);
            hvecxdata = 0.5;
            hmatdata = 2;
             
            Array<Dev<double>> matdata(hmatdata), vecxdata(hvecxdata), vecydata(hvecydata);
             
            Array<MatVecData> hmvdata(nummats);
             
            for (int i = 0; i < nummats; i++)
            {
                new (&hmvdata[i].mat) SliceMatrix<Dev<double>>(math, matw, matw, matdata.Data()+i*math*matw);
                // new (&hmvdata[i].x) BareVector(vecxdata.Data()+i*matw);
                // new (&hmvdata[i].y) BareVector(vecydata.Data()+i*math);
                hmvdata[i].offsetx = i*matw;
                hmvdata[i].offsety = (i/4*4)*math;
            }
            
            Array<Dev<MatVecData>> mvdata(hmvdata);
            FlatVector<Dev<double>> vx(vecxdata.Size(), vecxdata.Data());
            FlatVector<Dev<double>> vy(vecydata.Size(), vecydata.Data());
            SetScalar (0, vy);
                /*
            Dev<double> * pdev = Dev<double>::Malloc(1);
            pdev->H2D(0.0);
            vy = *pdev;
            */
                
            Timer t("manymatvec");
             size_t runs = 1;
               cudaDeviceSynchronize();
               t.Start();
               for (int i = 0; i < runs; i++)
                 ManyMatVec(mvdata, vx, vy);
               cudaDeviceSynchronize();
               t.Stop();
               t.AddFlops (double(runs)*nummats*math*matw);
               cout << "manymatvec, vecy = " << vecydata[0].D2H() << "," << vecydata[math].D2H()
                   << ", time = " << t.GetTime() << " MFlops = " << t.GetMFlops() << endl;
             }
         })
    ;
}

