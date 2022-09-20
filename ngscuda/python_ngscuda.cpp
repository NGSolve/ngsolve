#include <core/python_ngcore.hpp>

#include "../ngstd/python_ngstd.hpp"
#include "cuda_linalg.hpp"
#include "cuda_ngstd.hpp"

/*
 * TODO:
 *   always use ngs_cuda?
 */
using namespace ngla;
using namespace ngs_cuda;

PYBIND11_MODULE(ngscuda, m) {

  InitCUDA(1);

  py::class_<UnifiedVector, BaseVector, shared_ptr<UnifiedVector>>
    (m, "UnifiedVector", "UnifiedVector for CUDA applications")
                 
  .def(py::init([] (int asize) -> shared_ptr<UnifiedVector>
        { 
          return make_shared<UnifiedVector>(asize); 
        }))
  .def(py::init([] (const UnifiedVector &vec) -> shared_ptr<UnifiedVector>
        {
          return make_shared<UnifiedVector>(vec);
        }))
  .def(py::init([] (const BaseVector &vec) 
        {
          return make_shared<UnifiedVector>(vec);
        }))
  .def(py::init([] (py::array_t<double> bvec)
        {
          auto vec = bvec.template unchecked<1>();
          shared_ptr<UnifiedVector> v = make_shared<UnifiedVector>(vec.size());
          for (size_t i = 0; i < vec.size(); i++)
          {
            (*v)[i] = vec(i);
          }
          return v;
        }))
  /* .def(py::init([] (UnifiedVector &vec) -> shared_ptr<UnifiedVector> */
  /*       { */
  /*         return make_shared<UnifiedVector>(dynamic_cast<BaseVector&>(vec)); */
  /*       })) */
  /* .def(py::init([] (UnifiedVector &vec) -> shared_ptr<UnifiedVector> */
  /*       { */
  /*         cerr << "i am here." << endl; */
  /*         throw Exception("TODO"); */
  /*         /1* return make_shared<UnifiedVector>(vec); *1/ */
  /*       })) */
  /* .def(py::init([] (DynamicVectorExpression expr) -> shared_ptr<UnifiedVector> */
  /*       { */
  /*         cerr << "expr." << endl; */
  /*         throw Exception("TODO"); */
  /*       })) */
  .def("CreateVector", [] (UnifiedVector & self, bool copy)
         {
           auto newvec = self.CreateVector();
           if (copy) newvec = self;
           return shared_ptr<BaseVector>(newvec);
         }, py::arg("copy")=false,
         "creates a new vector of same type, contents is undefined if copy is false")
   

  /* .def("__len__", &UnifiedVector::Size) */

//  TODO: extend for splicing (define UnifiedVector.Range?)
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

  /* .def("PrintDevice", [] (UnifiedVector & self) */
  /*       { */
  /*         self.PrintDevice(); */
  /*       }) */

  // TODO: fix?
  .def("Add", [] (UnifiedVector & self, UnifiedVector & v2, py::object s) -> void
        {
          self.Add (py::extract<double>(s)(), v2);
          return;
        })

  // TODO: add complex / conjugate
  .def("InnerProduct", [] (UnifiedVector & self, UnifiedVector & v2) -> double
      {
        /* cerr << "innerproduct as method" << endl; */
        cerr << "InnerProduct sizes: " << self.Size() << " " << v2.Size() << endl; 
        return self.InnerProduct(v2);
      })
  .def("InnerProduct", [] (UnifiedVector & self, BaseVector & v2) -> double
      {
        /* cerr << "innerproduct as method" << endl; */
        cerr << "InnerProduct sizes: " << self.Size() << " " << v2.Size() << endl; 
        return self.InnerProduct(v2);
      })
  .def("UpdateHost", [] (UnifiedVector &self) -> void
        { self.UpdateHost(); }) 
  .def("UpdateDevice", [] (UnifiedVector &self) -> void
        { self.UpdateDevice(); })

  .def("__str__", [] (UnifiedVector & self) { return ToString<UnifiedVector>(self); } )
  .def("__repr__", [] (UnifiedVector & self) { return "unfiedvector"; } )

  .def("__len__", [] (UnifiedVector & self) {return self.Size(); })

  .def("Scale", [] (UnifiedVector & self, double d)
        {
          self.Scale(d);
        });

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

  py::class_<DevSparseMatrix, BaseMatrix, shared_ptr<DevSparseMatrix>>
    (m, "DevSparseMatrix", "DevSparseMatrix for CUDA applications")

    .def(py::init ( [] (SparseMatrix<double>& mat) -> shared_ptr<DevSparseMatrix>
          {
            return make_shared<DevSparseMatrix>(mat);
          }
          ))
    .def("Mult", [] (DevSparseMatrix &self, UnifiedVector &x, UnifiedVector &y)
        {
          self.Mult(x, y);
        })
    .def("Mult", [] (DevSparseMatrix &self, double s, UnifiedVector &x, UnifiedVector &y)
        {
          self.MultAdd(s, x, y);
        });
    /* .def("__mul__", [] (const DevSparseMatrix & a, const DevSparseMatrix & b) */
    /*       { */
    /*         cerr << "mat, mat prod" << endl; */
    /*         return a * b; */
    /*       }) */
    /* .def("__add__", [] (const DevSparseMatrix & a, const DevSparseMatrix & b) */
    /*       { */
    /*         cerr << "mat, mat add" << endl; */
    /*         return a + b; */
    /*       }) */
    /* .def("__sub__", [] (const DevSparseMatrix & a, const DevSparseMatrix & b) */
    /*       { */
    /*         cerr << "mat, mat sub" << endl; */
    /*         return a - b; */
    /*       }); */

  py::class_<DevMatrix, BaseMatrix, shared_ptr<DevMatrix>>
    (m, "DevMatrix", "device matrix for CUDA applications");

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
    .def("Mult", [] (DevDMatrix& self, DevDMatrix& b, DevDMatrix& c)
          {
            self.Mult(b, c);
          })
    .def("Scale", [] (DevDMatrix& self, double d)
          {
            self.Scale(d);
          })
    .def("SetZero", [] (DevDMatrix& self)
          {
            self.SetZero();
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

}

