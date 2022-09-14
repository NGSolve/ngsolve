#include <core/python_ngcore.hpp>

#include "../ngstd/python_ngstd.hpp"
#include "cuda_linalg.hpp"
#include "cuda_ngstd.hpp"

/*
 * TODO:
 * 	always use ngs_cuda?
 */
using namespace ngla;
using namespace ngs_cuda;

PYBIND11_MODULE(ngscuda, m) {
	/*
	 * TODO: 
	 * 	-) matrix to DeviceMatrix Conversion
	 * 	-) mat-rix-vector product
	 * 	-) Is it possible to use the usual functions (f.e. iterative method
	 * 		using jacobi) without redifining only using f.e. operator+
	 */

	InitCUDA();

  py::class_<UnifiedVector, BaseVector, shared_ptr<UnifiedVector>>
    (m, "UnifiedVector", "UnifiedVector for CUDA applications")
                 
	.def(py::init([] (int asize) -> shared_ptr<UnifiedVector>
				{ 
					return make_shared<UnifiedVector>(asize); 
				}))
	.def(py::init([] (BaseVector &vec) -> shared_ptr<UnifiedVector>
				{
					return make_shared<UnifiedVector>(vec);
				}))
	.def("CreateVector", [] (UnifiedVector & self, bool copy)
         {
           auto newvec = self.CreateVector();
           if (copy) newvec = self;
           return shared_ptr<BaseVector>(newvec);
         }, py::arg("copy")=false,
         "creates a new vector of same type, contents is undefined if copy is false")
   

	.def("__len__", &UnifiedVector::Size)

//	TODO: extend for splicing (define UnifiedVector.Range?)
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

	.def("PrintDevice", [] (UnifiedVector & self)
				{
					self.PrintDevice();
				})

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
				return self.InnerProduct(v2);
			})

	.def("UpdateHost", [] (UnifiedVector &self) -> void
				{ self.UpdateHost(); }) 
	.def("UpdateDevice", [] (UnifiedVector &self) -> void
				{ self.UpdateDevice(); })

	// TODO: somewhere must be a cpy Host2Device... check it
  .def("__str__", [](UnifiedVector &self) { return ToString<UnifiedVector>(self); } )
  .def("__repr__", [](UnifiedVector &self) { return "unfiedvector"; } );

	/* TODO: */	
	/* .def("Set") */
	/* .def("__add__", (UnifiedVector & v1, BaseVector & v2, py::object s) -> void */
	/* 		{ */
	/* 			return v1 + v2; */
	/* 		}) */
	/* .def("Scale", [](UnifiedVector & self, double scal) -> shared_ptr<UnifiedVector> */
	/* 		{ cerr << "scaling by " << scal << endl; return self.Scale(scal); }) */

	// TODO: probably useless... use ToGPU instead?
	//
	// 	alternative: BaesVector::ToGPU
	/* .def("Assign", [] (UnifiedVector &self, BaseVector &v2) -> void */
	/* 			{ */
	/* 				self = v2; */
	/* 				return; */
	/* 			}) */

	/* .def("Assign", [] (UnifiedVector &self, py::object d) -> void */
	/* 			{ */
	/* 				self.SetScalar(py::extract<double>(d)()); */
	/* 				return; */
	/* 			}) */

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

<<<<<<< HEAD
	py::class_<DevMatrix, BaseMatrix, shared_ptr<DevMatrix>>
		(m, "DevMatrix", "DevMatrix for CUDA applications");
		
		/* .def(py::init ( [] (BaseMatrix & mat) -> shared_ptr<BaseMatrix> */
	m.def("CreateDevMatrix", [] (BaseMatrix &mat)
					{
						/* return make_shared<DevMatrix>(mat); */
						return CreateDevMatrix(mat);
					});

	/* m.def("CreateDevMatrix", [] (BaseMatrix & mat) */
	/* 		{ */
	/* 			return CreateDevMatrix(mat); */
	/* 		}); */
=======
>>>>>>> 6d28ada61482d85ff3a9bc60e9d1ab0e8a5ca709
}

