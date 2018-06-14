#ifdef NGS_PYTHON
#include <la.hpp>
#include <parallelngs.hpp>
#include "../ngstd/python_ngstd.hpp"
using namespace ngla;


template<typename T>
void ExportSparseMatrix(py::module m)
{
  py::class_<SparseMatrix<T>, shared_ptr<SparseMatrix<T>>, BaseSparseMatrix, S_BaseMatrix<typename mat_traits<T>::TSCAL>>
    (m, (string("SparseMatrix") + typeid(T).name()).c_str(),
     "a sparse matrix in CSR storage")
    .def("__getitem__",
         [](const SparseMatrix<T> & self, py::tuple t)
         {
           size_t row = t[0].cast<size_t>();
           size_t col = t[1].cast<size_t>();
           return self(row,col);
         })
    .def("__setitem__",
         [](SparseMatrix<T> & self, py::tuple t, T value)
         {
           size_t row = t[0].cast<size_t>();
           size_t col = t[1].cast<size_t>();
           self(row,col) = value;
         })

    .def("COO", [] (SparseMatrix<T> * sp) -> py::object
         {
           size_t nze = sp->NZE();
           Array<int> ri(nze), ci(nze);
           Vector<T> vals(nze);
           for (size_t i = 0, ii = 0; i < sp->Height(); i++)
             {
               FlatArray<int> ind = sp->GetRowIndices(i);
               FlatVector<T> rv = sp->GetRowValues(i);
               for (int j = 0; j < ind.Size(); j++, ii++)
                 {
                   ri[ii] = i;
                   ci[ii] = ind[j];
                   vals[ii] = rv[j];
                 }
             }
           /*
           t2.Start();
           // still copies, we don't understand why
           py::object pyri = py::cast(std::move(ri));
           py::object pyci = py::cast(std::move(ci));
           py::object pyvals = py::cast(std::move(vals));
           t2.Stop();
           return py::make_tuple (pyri, pyci, pyvals);
           */
           // moves the arrays
           return py::make_tuple (move(ri), move(ci), move(vals));
         })
    
    .def("CRS", [] (SparseMatrix<T> * sp) -> py::object
         {
           FlatArray<int> colind(sp->NZE(), sp->GetRowIndices(0).Addr(0));
           FlatVector<T> values(sp->NZE(), sp->GetRowValues(0).Addr(0));
           FlatArray<size_t> first = sp->GetFirstArray();
           return py::make_tuple (colind, values, first); 
         })
    
    .def_static("CreateFromCOO",
                [] (py::list indi, py::list indj, py::list values, size_t h, size_t w)
                {
                  auto cindi = makeCArray<int>(indi);
                  auto cindj = makeCArray<int>(indj);
                  auto cvalues = makeCArray<double>(values);
                  return SparseMatrix<double>::CreateFromCOO (cindi,cindj,cvalues, h,w);
                })
    
    .def("CreateTranspose", [] (const SparseMatrix<double> & sp)
         { return TransposeMatrix (sp); })

    .def("__matmul__", [] (const SparseMatrix<double> & a, const SparseMatrix<double> & b)
         { return MatMult(a,b); })
    .def("__matmul__", [](shared_ptr<SparseMatrix<double>> a, shared_ptr<BaseMatrix> mb)
         ->shared_ptr<BaseMatrix> { return make_shared<ProductMatrix> (a, mb); })
    ;

  py::class_<SparseMatrixSymmetric<T>, shared_ptr<SparseMatrixSymmetric<T>>, SparseMatrix<T>>
    (m, (string("SparseMatrixSymmetric") + typeid(T).name()).c_str());
}

void NGS_DLL_HEADER ExportNgla(py::module &m) {

  py::enum_<PARALLEL_STATUS>(m, "PARALLEL_STATUS", "enum of possible parallel ")
    .value("DISTRIBUTED", DISTRIBUTED)
    .value("CUMULATED", CUMULATED)
    .value("NOT_PARALLEL", NOT_PARALLEL)
    .export_values()
    ;
    
  py::class_<ParallelDofs, shared_ptr<ParallelDofs>> (m, "ParallelDofs")
#ifdef PARALLEL
    .def("SubSet", [](const ParallelDofs & self, shared_ptr<BitArray> take_dofs) { 
        return self.SubSet(take_dofs); })
#endif
    .def_property_readonly ("ndoflocal", [](const ParallelDofs & self) 
			    { return self.GetNDofLocal(); },
                            "number of degrees of freedom")

    .def_property_readonly ("ndofglobal",
                            [](const ParallelDofs & self) 
			    { return self.GetNDofGlobal(); },    
                            "number of global degrees of freedom")
    .def("ExchangeProcs", [] (const ParallelDofs & self)
         { return self.GetDistantProcs(); } )
    .def("Dof2Proc", [] (const ParallelDofs & self, int dof)
         { return self.GetDistantProcs(dof); })
    .def("Proc2Dof", [] (const ParallelDofs & self, int proc)
         { return self.GetExchangeDofs(proc); })
    ;

    m.def("CreateVVector",
          [] (size_t s, bool is_complex, int es) -> shared_ptr<BaseVector>
          { return CreateBaseVector(s,is_complex, es); },
          "size"_a, "complex"_a=false, "entrysize"_a=1);
    
  py::class_<BaseVector, shared_ptr<BaseVector>>(m, "BaseVector",
        py::dynamic_attr() // add dynamic attributes
      )
    .def(py::init([] (size_t s, bool is_complex, int es) -> shared_ptr<BaseVector>
                  { return CreateBaseVector(s,is_complex, es); }),
                  "size"_a, "complex"_a=false, "entrysize"_a=1)
#ifdef PARALLEL
    .def_property_readonly("local_vec", [](BaseVector & self) -> shared_ptr<BaseVector> {
	auto * pv = dynamic_cast_ParallelBaseVector (&self);
	if (pv==NULL) throw Exception("Only ParallelVectors have a local Vector!");
	auto rv = pv->GetLocalVector();
	return rv;
      } )
#else
    .def_property_readonly("local_vec", [](BaseVector & self) -> shared_ptr<BaseVector> {
	return self;
      } )
#endif
    .def(py::pickle([] (const BaseVector& bv)
                    {
                      MemoryView mv((void*) &bv.FVDouble()[0], sizeof(double) * bv.FVDouble().Size());
                      return py::make_tuple(bv.Size(),bv.IsComplex(),bv.EntrySize(),mv);
                    },
                    [] (py::tuple state) -> shared_ptr<BaseVector>
                    {
                      auto mv = state[3].cast<MemoryView>();
                      shared_ptr<BaseVector> bv;
                      if (state[1].cast<bool>())
                        {
                          // create basevector with owning pointer and afterwards assign it to mem
                          auto bptr = make_shared<S_BaseVectorPtr<Complex>>(0, state[2].cast<size_t>());
                          bptr->AssignMemory(state[0].cast<size_t>(), mv.Ptr());
                          return bptr;
                        }
                      else
                        {
                          // create basevector with owning pointer and afterwards assign it to mem
                          auto bptr = make_shared<S_BaseVectorPtr<double>>(0, state[2].cast<size_t>());
                          bptr->AssignMemory(state[0].cast<size_t>(), mv.Ptr());
                          return bptr;
                        }
                    }
                    ))
    .def("__str__", [](BaseVector &self) { return ToString<BaseVector>(self); } )
    .def("__repr__", [](BaseVector &self) { return "basevector"; } )
    .def_property_readonly("size", py::cpp_function( [] (BaseVector &self) { return self.Size(); } ) )
    .def("__len__", [] (BaseVector &self) { return self.Size(); })
    .def("CreateVector", [] (BaseVector & self)
         { return shared_ptr<BaseVector>(self.CreateVector()); },
         "creates a new vector of same type, contents is undefined")
    
    .def("Copy", [] (BaseVector & self)
         {
           auto hv = shared_ptr<BaseVector>(self.CreateVector());
           *hv = self;
           return hv;
         },
         "creates a new vector of same type, copy contents")
    
    .def("Assign",[](BaseVector & self, BaseVector & v2, py::object s)->void
                                   { 
                                     if ( py::extract<double>(s).check() )
                                       {
                                         self.Set (py::extract<double>(s)(), v2);
                                         return;
                                       }
                                     if ( py::extract<Complex>(s).check() )
                                       {
                                         self.Set (py::extract<Complex>(s)(), v2);
                                         return;
                                       }
                                     throw Exception ("BaseVector::Assign called with non-scalar type");
                                   })
    .def("Add",[](BaseVector & self, BaseVector & v2, py::object s)->void
                                   { 
                                     if ( py::extract<double>(s).check() )
                                       {
                                         self.Add (py::extract<double>(s)(), v2);
                                         return;
                                       }
                                     if ( py::extract<Complex>(s).check() )
                                       {
                                         self.Add (py::extract<Complex>(s)(), v2);
                                         return;
                                       }
                                     throw Exception ("BaseVector::Assign called with non-scalar type");
                                   })


    // TODO
//     .add_property("expr", py::object(expr_namespace["VecExpr"]) )
//     .add_property("data", py::object(expr_namespace["VecExpr"]), py::object(expr_namespace["expr_data"] ))
//     .def("__add__" , py::object(expr_namespace["expr_add"]) )
//     .def("__sub__" , py::object(expr_namespace["expr_sub"]) )
//     .def("__rmul__" , py::object(expr_namespace["expr_rmul"]) )
    .def("__getitem__",
          [](BaseVector & self,  int ind )
           {
             if (ind < 0 || ind >= self.Size()) 
               throw py::index_error();
             int entrysize = self.EntrySize();
             if( self.IsComplex() ) entrysize/=2;
             if(entrysize == 1)
             {
                 if( !self.IsComplex() )
                     return py::cast(self.FVDouble()[ind]);
                 else
                     return py::cast(self.FVComplex()[ind]);
             }
             else
             {
                 // return FlatVector<T>
                 if( self.IsComplex() )
                   return py::cast(self.SV<Complex>()(ind));
                 else
                   return py::cast(self.SV<double>()(ind));
             }
           } )
    .def("__getitem__", [](BaseVector & self,  py::slice inds )
      {
          size_t start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if (step != 1)
            throw Exception ("slices with non-unit distance not allowed");
          return shared_ptr<BaseVector>(self.Range(start, start+n));
      } )
    .def("__setitem__", [](BaseVector & self,  int ind, double d )
      {
          self.Range(ind,ind+1) = d;
      } )
    .def("__setitem__", [](BaseVector & self,  int ind, Complex z )
      {
        self.Range(ind,ind+1) = z;
      } )
    .def("__setitem__", [](BaseVector & self,  py::slice inds, double d )
      {
          size_t start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if (step != 1)
            throw Exception ("slices with non-unit distance not allowed");          
	  if(n==self.Size()) {
	    self.SetScalar(d);
	    return;
	  }
          self.Range(start,start+n) = d;
      } )
    .def("__setitem__", [](BaseVector & self,  py::slice inds, Complex z )
      {
          size_t start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if (step != 1)
            throw Exception ("slices with non-unit distance not allowed");          
          self.Range(start,start+n) = z;
      } )
    .def("__setitem__", [](BaseVector & self, py::slice inds, shared_ptr<BaseVector> v )
      {
        size_t start, step, n;
        InitSlice( inds, self.Size(), start, step, n );
        if (step != 1)
          throw Exception ("slices with non-unit distance not allowed");        
        self.Range(start, start+n) = *v;
      } )
    .def("__setitem__", [](BaseVector & self,  int ind, FlatVector<double> & v )
      {
          if( self.IsComplex() )
            self.SV<Complex>()(ind) = v;
          else
            self.SV<double>()(ind) = v;
      } )
    .def("__setitem__", [](BaseVector & self,  int ind, FlatVector<Complex> & v )
      {
          if( self.IsComplex() )
            self.SV<Complex>()(ind) = v;
          else
            throw py::index_error("cannot assign complex values to real vector");
      } )
    .def("__iadd__", [](BaseVector & self,  BaseVector & other) -> BaseVector& { self += other; return self;})
    .def("__isub__", [](BaseVector & self,  BaseVector & other) -> BaseVector& { self -= other; return self;})
    .def("__imul__", [](BaseVector & self,  double scal) -> BaseVector& { self *= scal; return self;})
    .def("__imul__", [](BaseVector & self,  Complex scal) -> BaseVector& { self *= scal; return self;})
    .def("__itruediv__", [](BaseVector & self,  double scal) -> BaseVector& { self /= scal; return self;})
    .def("__itruediv__", [](BaseVector & self,  Complex scal) -> BaseVector& { self /= scal; return self;})
    .def("InnerProduct", [](BaseVector & self, BaseVector & other, bool conjugate)
                                          {
                                            if (self.IsComplex())
                                              {
                                                if (conjugate)
                                                  return py::cast (S_InnerProduct<ComplexConjugate> (self, other));
                                                else
                                                  return py::cast (S_InnerProduct<Complex> (self, other));
                                              }
                                            else
                                              return py::cast (InnerProduct (self, other));
                                          },
         "InnerProduct", py::arg("other"), py::arg("conjugate")=py::cast(true)         
         )
    .def("Norm",  [](BaseVector & self) { return self.L2Norm(); })
    .def("Range", [](BaseVector & self, int from, int to) -> shared_ptr<BaseVector>
                                   {
                                     return shared_ptr<BaseVector>(self.Range(from,to));
                                   })
    .def("FV", [] (BaseVector & self)
                                {
                                  if (!self.IsComplex())
                                    return py::cast(self.FVDouble());
                                  else
                                    return py::cast(self.FVComplex());
                                })
    .def("Distribute", [] (BaseVector & self) { self.Distribute(); } ) 
    .def("Cumulate", [] (BaseVector & self) { self.Cumulate(); } ) 
    .def("GetParallelStatus", [] (BaseVector & self) { return self.GetParallelStatus(); } )
    .def("SetParallelStatus", [] (BaseVector & self, PARALLEL_STATUS stat) { self.SetParallelStatus(stat); });

  // m.def("InnerProduct",[](BaseVector & v1, BaseVector & v2)->double { return InnerProduct(v1,v2); })
  m.def ("InnerProduct",
           [] (py::object x, py::object y) -> py::object
                            { return py::handle(x.attr("InnerProduct")) (y); });
  ;
  

  py::class_<BlockVector, BaseVector, shared_ptr<BlockVector>> (m, "BlockVector")
    .def(py::init<> ([] (vector<shared_ptr<BaseVector>> vecs)
                     {
                       Array<shared_ptr<BaseVector>> v2;
                       for (auto v : vecs) v2 += v;
                       return make_shared<BlockVector> (v2);
                     }))
    
    .def("__getitem__", [](BlockVector & self, int ind) { return self[ind]; })
    ;




  
  typedef BaseMatrix BM;
  // typedef BaseVector BV;

    class BaseMatrixTrampoline : public BaseMatrix {
    public:
      using BaseMatrix::BaseMatrix;
      BaseMatrixTrampoline() : BaseMatrix()
          {
            static_assert( sizeof(BaseMatrix)==sizeof(BaseMatrixTrampoline), "slkdf");
          }

      bool IsComplex() const override { 
        PYBIND11_OVERLOAD_PURE(
            bool, /* Return type */
            BaseMatrix,      /* Parent class */
            IsComplex,          /* Name of function */
            );
      }
      
      int VHeight() const override { 
        PYBIND11_OVERLOAD_PURE(
            int, /* Return type */
            BaseMatrix,      /* Parent class */
            Height,          /* Name of function */
            );
      }

      int VWidth() const override { 
        PYBIND11_OVERLOAD_PURE(
            int, /* Return type */
            BaseMatrix,      /* Parent class */
            Width,          /* Name of function */
            );
      }
      
      void Mult (const BaseVector & x, BaseVector & y) const override {
        pybind11::gil_scoped_acquire gil;
        pybind11::function overload = pybind11::get_overload(this, "Mult");
        if (overload) {
	  const AutoVector * avecx = dynamic_cast<const AutoVector*>(&x);
          auto sx = shared_ptr<BaseVector>(const_cast<BaseVector*>((avecx!=NULL)?&(**avecx):&x),
					   NOOP_Deleter);
	  const AutoVector * avecy = dynamic_cast<const AutoVector*>(&y);
          auto sy = shared_ptr<BaseVector>(const_cast<BaseVector*>((avecy!=NULL)?&(**avecy):&y),
					   NOOP_Deleter);
          overload(sx,sy);
        }
        else
          BaseMatrix::Mult(x,y);
      }
      void MultAdd (double s, const BaseVector & x, BaseVector & y) const override {
        pybind11::gil_scoped_acquire gil;
        pybind11::function overload = pybind11::get_overload(this, "MultAdd");
        if (overload) {
	  const AutoVector * avecx = dynamic_cast<const AutoVector*>(&x);
          auto sx = shared_ptr<BaseVector>(const_cast<BaseVector*>((avecx!=NULL)?&(**avecx):&x),
					   NOOP_Deleter);
	  const AutoVector * avecy = dynamic_cast<const AutoVector*>(&y);
          auto sy = shared_ptr<BaseVector>(const_cast<BaseVector*>((avecy!=NULL)?&(**avecy):&y),
					   NOOP_Deleter);
          overload(s, sx,sy);
        }
        else
          BaseMatrix::MultAdd(s, x, y);
      }
      void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override {
        pybind11::gil_scoped_acquire gil;
        pybind11::function overload = pybind11::get_overload(this, "MultTransAdd");
        if (overload) {
	  const AutoVector * avecx = dynamic_cast<const AutoVector*>(&x);
          auto sx = shared_ptr<BaseVector>(const_cast<BaseVector*>((avecx!=NULL)?&(**avecx):&x),
					   NOOP_Deleter);
	  const AutoVector * avecy = dynamic_cast<const AutoVector*>(&y);
          auto sy = shared_ptr<BaseVector>(const_cast<BaseVector*>((avecy!=NULL)?&(**avecy):&y),
					   NOOP_Deleter);
          overload(s, sx,sy);
        }
        else
          BaseMatrix::MultTransAdd(s, x, y);
      }


      void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override {
        pybind11::gil_scoped_acquire gil;
        pybind11::function overload = pybind11::get_overload(this, "MultAdd");
        if (overload) {
	  const AutoVector * avecx = dynamic_cast<const AutoVector*>(&x);
          auto sx = shared_ptr<BaseVector>(const_cast<BaseVector*>((avecx!=NULL)?&(**avecx):&x),
					   NOOP_Deleter);
	  const AutoVector * avecy = dynamic_cast<const AutoVector*>(&y);
          auto sy = shared_ptr<BaseVector>(const_cast<BaseVector*>((avecy!=NULL)?&(**avecy):&y),
					   NOOP_Deleter);
          overload(s, sx,sy);
        }
        else
          BaseMatrix::MultAdd(s, x, y);
      }
      void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override {
        pybind11::gil_scoped_acquire gil;
        pybind11::function overload = pybind11::get_overload(this, "MultTransAdd");
        if (overload) {
	  const AutoVector * avecx = dynamic_cast<const AutoVector*>(&x);
          auto sx = shared_ptr<BaseVector>(const_cast<BaseVector*>((avecx!=NULL)?&(**avecx):&x),
					   NOOP_Deleter);
	  const AutoVector * avecy = dynamic_cast<const AutoVector*>(&y);
          auto sy = shared_ptr<BaseVector>(const_cast<BaseVector*>((avecy!=NULL)?&(**avecy):&y),
					   NOOP_Deleter);
          overload(s, sx,sy);
        }
        else
          BaseMatrix::MultTransAdd(s, x, y);
      }


    };
  


  py::class_<BaseMatrix, shared_ptr<BaseMatrix>, BaseMatrixTrampoline>(m, "BaseMatrix")
    /*
    .def("__init__", [](BaseMatrix *instance) { 
        new (instance) BaseMatrixTrampoline(); }
        )
    */
    .def(py::init<> ([]() { return new BaseMatrixTrampoline(); }))
    .def("__str__", [](BaseMatrix &self) { return ToString<BaseMatrix>(self); } )
    .def_property_readonly("height", [] ( BaseMatrix & self)
        { return self.Height(); } )
    .def_property_readonly("width", [] ( BaseMatrix & self)
        { return self.Width(); } )
    .def_property_readonly("nze", [] ( BaseMatrix & self)
                           { return self.NZE(); }, "number of non-zero elements")
    // .def("CreateMatrix", &BaseMatrix::CreateMatrix)
    .def("CreateMatrix", [] ( BaseMatrix & self)
        { return self.CreateMatrix(); } )

    
    .def("CreateRowVector", [] ( BaseMatrix & self)
        { return shared_ptr<BaseVector>(self.CreateRowVector()); } )
    .def("CreateColVector", [] ( BaseMatrix & self)
        { return shared_ptr<BaseVector>(self.CreateColVector()); } )
    
    .def("AsVector", [] (BM & m)
                                      {
                                        return shared_ptr<BaseVector> (&m.AsVector(), NOOP_Deleter);
                                      })

    .def("Mult",         [](BaseMatrix &m, BaseVector &x, BaseVector &y) { m.Mult(x, y); }, py::call_guard<py::gil_scoped_release>())
    .def("MultAdd",      [](BaseMatrix &m, double s, BaseVector &x, BaseVector &y) { m.MultAdd (s, x, y); }, py::call_guard<py::gil_scoped_release>())
    .def("MultTrans",    [](BaseMatrix &m, double s, BaseVector &x, BaseVector &y) { y=0; m.MultTransAdd (1.0, x, y); }, py::call_guard<py::gil_scoped_release>())
    .def("MultTransAdd",  [](BaseMatrix &m, double s, BaseVector &x, BaseVector &y) { m.MultTransAdd (s, x, y); }, py::call_guard<py::gil_scoped_release>())
    .def("MultScale",    [](BaseMatrix &m, double s, BaseVector &x, BaseVector &y)
          {
              m.Mult (x,y);
              if(s!=1.0)
                  y *= s;
          } , py::call_guard<py::gil_scoped_release>())
    .def("MultAdd",      [](BaseMatrix &m, Complex s, BaseVector &x, BaseVector &y) { m.MultAdd (s, x, y); }, py::call_guard<py::gil_scoped_release>())
    .def("MultTrans",    [](BaseMatrix &m, Complex s, BaseVector &x, BaseVector &y) { y=0; m.MultTransAdd (1.0, x, y); }, py::call_guard<py::gil_scoped_release>())
    .def("MultTransAdd",  [](BaseMatrix &m, Complex s, BaseVector &x, BaseVector &y) { m.MultTransAdd (s, x, y); }, py::call_guard<py::gil_scoped_release>())
    .def("MultScale",    [](BaseMatrix &m, Complex s, BaseVector &x, BaseVector &y)
          {
              m.Mult (x,y);
              if(s!=1.0)
                  y *= s;
          }, py::call_guard<py::gil_scoped_release>() )

    .def("__iadd__", [] (BM &m, BM &m2) { 
        m.AsVector()+=m2.AsVector();
    }, py::call_guard<py::gil_scoped_release>())

    .def("GetInverseType", [](BM & m)
                                            {
                                              return GetInverseName( m.GetInverseType());
                                            })

    .def("Inverse", [](BM &m, shared_ptr<BitArray> freedofs, string inverse)
                                     { 
                                       if (inverse != "") m.SetInverseType(inverse);
                                       return m.InverseMatrix(freedofs);
                                     }
         ,"Inverse", py::arg("freedofs")=nullptr, py::arg("inverse")=py::str(""), 
         docu_string(R"raw_string(Calculate inverse of sparse matrix
Parameters

freedofs : BitArray
  If set, invert only the rows/columns the matrix defined by the bit array, otherwise invert the whole matrix

inverse : string
  Solver to use, allowed values are:
    sparsecholesky - internal solver of NGSolve for symmetric matrices
    umfpack        - solver by Suitesparse/UMFPACK (if NGSolve was configured with USE_UMFPACK=ON)
    pardiso        - PARDISO, either provided by libpardiso (USE_PARDISO=ON) or Intel MKL (USE_MKL=ON).
                     If neither Pardiso nor Intel MKL was linked at compile-time, NGSolve will look
                     for libmkl_rt in LD_LIBRARY_PATH (Unix) or PATH (Windows) at run-time.
)raw_string"), py::call_guard<py::gil_scoped_release>())
    // .def("Inverse", [](BM &m)  { return m.InverseMatrix(); })

    .def_property_readonly("T", [](shared_ptr<BM> m)->shared_ptr<BaseMatrix> { return make_shared<Transpose> (m); })
    .def("__matmul__", [](shared_ptr<BM> ma, shared_ptr<BM> mb)->shared_ptr<BaseMatrix>
         { return make_shared<ProductMatrix> (ma, mb); })
    .def("__add__", [](shared_ptr<BM> ma, shared_ptr<BM> mb)->shared_ptr<BaseMatrix>
         { return make_shared<SumMatrix> (ma, mb, 1, 1); })
    .def("__sub__", [](shared_ptr<BM> ma, shared_ptr<BM> mb)->shared_ptr<BaseMatrix>
         { return make_shared<SumMatrix> (ma, mb, 1, -1); })
    .def("__rmul__", [](shared_ptr<BM> ma, double a)->shared_ptr<BaseMatrix>
         { return make_shared<VScaleMatrix<double>> (ma, a); })
    .def("__rmul__", [](shared_ptr<BM> ma, Complex a)->shared_ptr<BaseMatrix>
         { return make_shared<VScaleMatrix<Complex>> (ma, a); })
    .def("Update", [](BM &m) { m.Update(); }, py::call_guard<py::gil_scoped_release>())
    ;

  py::class_<BaseSparseMatrix, shared_ptr<BaseSparseMatrix>, BaseMatrix>
    (m, "BaseSparseMatrix", "sparse matrix of any type")
    
    .def("CreateSmoother", [](BaseSparseMatrix & m, shared_ptr<BitArray> ba) 
         { return m.CreateJacobiPrecond(ba); },
         py::arg("freedofs") = shared_ptr<BitArray>())
    
    .def("CreateBlockSmoother", [](BaseSparseMatrix & m, py::object blocks)
         {
           size_t size = py::len(blocks);
           
           Array<int> cnt(size);
           size_t i = 0;
           for (auto block : blocks)
             cnt[i++] = py::len(block);
           
           i = 0;
           Table<int> blocktable(cnt);
           for (auto block : blocks)
             {
               auto row = blocktable[i++];
               size_t j = 0;
               for (auto val : block)
                 row[j++] = val.cast<int>();
             }

           auto pre = m.CreateBlockJacobiPrecond (make_shared<Table<int>> (move(blocktable)));
           return pre;
         })
     ;

  py::class_<S_BaseMatrix<double>, shared_ptr<S_BaseMatrix<double>>, BaseMatrix>
    (m, "S_BaseMatrixD", "base sparse matrix");
  py::class_<S_BaseMatrix<Complex>, shared_ptr<S_BaseMatrix<Complex>>, BaseMatrix>
    (m, "S_BaseMatrixC", "base sparse matrix");


  py::class_<BlockMatrix, BaseMatrix, shared_ptr<BlockMatrix>> (m, "BlockMatrix")
    .def(py::init<> ([] (vector<vector<shared_ptr<BaseMatrix>>> mats)
                     {
                       Array<Array<shared_ptr<BaseMatrix>>> m2;
                       for (auto mrow : mats)
                         {
                           Array<shared_ptr<BaseMatrix>> mrow2;
                           for (auto m : mrow) mrow2 += m;
                           m2 += mrow2;
                         }
                       return make_shared<BlockMatrix> (m2);
                     }))
    
    // .def("__getitem__", [](BlockMatrix & self, int row, int col) { return self(row,rol); })
    ;


  
#ifdef PARALLEL
  py::class_<ParallelMatrix, shared_ptr<ParallelMatrix>, BaseMatrix>
    (m, "ParallelMatrix", "MPI-distributed matrix")
    .def(py::init<shared_ptr<BaseMatrix>, shared_ptr<ParallelDofs>>())
    .def_property_readonly("local_mat", [](ParallelMatrix & mat) { return mat.GetMatrix(); })
    ;

  py::class_<FETI_Jump_Matrix, shared_ptr<FETI_Jump_Matrix>, BaseMatrix>
    (m, "FETI_Jump", "B-matrix of the FETI-system")
    .def(py::init<shared_ptr<ParallelDofs>>())
    ;
#endif
  

  
  ExportSparseMatrix<double>(m);
  ExportSparseMatrix<Complex>(m);
  ExportSparseMatrix<Mat<2,2,double>>(m);
  ExportSparseMatrix<Mat<2,2,Complex>>(m);
  ExportSparseMatrix<Mat<3,3,double>>(m);
  ExportSparseMatrix<Mat<3,3,Complex>>(m);
  
  py::class_<BaseBlockJacobiPrecond, shared_ptr<BaseBlockJacobiPrecond>, BaseMatrix>
    (m, "BlockSmoother",
     "block Jacobi and block Gauss-Seidel smoothing")
    .def("Smooth", &BaseBlockJacobiPrecond::GSSmooth,
         py::arg("x"), py::arg("b"), py::arg("steps"),
         "performs steps block-Gauss-Seidel iterations for the linear system A x = b")
    .def("SmoothBack", &BaseBlockJacobiPrecond::GSSmoothBack,
         py::arg("x"), py::arg("b"), py::arg("steps"),
         "performs steps block-Gauss-Seidel iterations for the linear system A x = b in reverse order")
    ;

  py::class_<BaseJacobiPrecond, shared_ptr<BaseJacobiPrecond>, BaseMatrix>
    (m, "Smoother",
     "Jacobi and Gauss-Seidel smoothing")
    .def("Smooth", [&](BaseJacobiPrecond & jac, BaseVector & x, BaseVector & b)
         { jac.GSSmooth (x, b); },
         py::arg("x"), py::arg("b"),
         "performs steps Gauss-Seidel iterations for the linear system A x = b")
    .def("SmoothBack", &BaseJacobiPrecond::GSSmoothBack,
         py::arg("x"), py::arg("b"),
         "performs steps Gauss-Seidel iterations for the linear system A x = b in reverse order")
    ;

  py::class_<SparseFactorization, shared_ptr<SparseFactorization>, BaseMatrix>
    (m, "SparseFactorization")
    .def("Smooth", [] (SparseFactorization & self, BaseVector & u, BaseVector & y)
         {
           self.Smooth (u, y /* this is not needed */, y);
         }, "perform smoothing step (needs non-symmetric storage so symmetric sparse matrix)")
    ;

  py::class_<SparseCholesky<double>, shared_ptr<SparseCholesky<double>>, SparseFactorization> (m, "SparseCholesky_d");
  py::class_<SparseCholesky<Complex>, shared_ptr<SparseCholesky<Complex>>, SparseFactorization> (m, "SparseCholesky_c");
  
  py::class_<Projector, shared_ptr<Projector>, BaseMatrix> (m, "Projector")
    .def(py::init<shared_ptr<BitArray>,bool>());
    ;

    py::class_<ngla::IdentityMatrix, shared_ptr<ngla::IdentityMatrix>, BaseMatrix> (m, "IdentityMatrix")
    .def(py::init<>())
    ;

  py::class_<KrylovSpaceSolver, shared_ptr<KrylovSpaceSolver>, BaseMatrix> (m, "KrylovSpaceSolver")
    .def("GetSteps", &KrylovSpaceSolver::GetSteps)
    ;

  m.def("CGSolver", [](const BaseMatrix & mat, const BaseMatrix & pre,
                                          bool iscomplex, bool printrates, 
                                          double precision, int maxsteps)
                                       {
                                         KrylovSpaceSolver * solver;
                                         if(mat.IsComplex()) iscomplex = true;
                                         
                                         if (iscomplex)
                                           solver = new CGSolver<Complex> (mat, pre);
                                         else
                                           solver = new CGSolver<double> (mat, pre);
                                         solver->SetPrecision(precision);
                                         solver->SetMaxSteps(maxsteps);
                                         solver->SetPrintRates (printrates);
                                         return shared_ptr<KrylovSpaceSolver>(solver);
                                       },
          "CG Solver", py::arg("mat"), py::arg("pre"), py::arg("complex") = false, py::arg("printrates")=true,
           py::arg("precision")=1e-8, py::arg("maxsteps")=200
          )
    ;

  m.def("GMRESSolver", [](const BaseMatrix & mat, const BaseMatrix & pre,
                                           bool printrates, 
                                           double precision, int maxsteps)
                                        {
                                          KrylovSpaceSolver * solver;
                                          if (!mat.IsComplex())
                                            solver = new GMRESSolver<double> (mat, pre);
                                          else
                                            solver = new GMRESSolver<Complex> (mat, pre);                                            
                                          solver->SetPrecision(precision);
                                          solver->SetMaxSteps(maxsteps);
                                          solver->SetPrintRates (printrates);
                                          return shared_ptr<KrylovSpaceSolver>(solver);
                                        },
          "GMRES Solver", py::arg("mat"), py::arg("pre"), py::arg("printrates")=true,
           py::arg("precision")=1e-8, py::arg("maxsteps")=200
          )
    ;

  m.def("TestPC", [](const BaseMatrix & mat, const BaseMatrix & pre) {
      EigenSystem eigen(mat, pre);
      eigen.Calc();
      cout << IM(1) << " Min Eigenvalue : "  << eigen.EigenValue(1) << endl; 
      cout << IM(1) << " Max Eigenvalue : " << eigen.MaxEigenValue() << endl; 
      cout << IM(1) << " Condition   " << eigen.MaxEigenValue()/eigen.EigenValue(1) << endl; 
      return;
    },
    py::arg("mat"), py::arg("pre"));
  
  py::class_<QMRSolver<double>, shared_ptr<QMRSolver<double>>, BaseMatrix> (m, "QMRSolverD")
    ;
  py::class_<QMRSolver<Complex>, shared_ptr<QMRSolver<Complex>>, BaseMatrix> (m, "QMRSolverC")
    ;

  m.def("QMRSolver", [](const BaseMatrix & mat, const BaseMatrix & pre,
                                           bool printrates, 
                                           double precision, int maxsteps)
                                        {
                                          KrylovSpaceSolver * solver;
                                          if (!mat.IsComplex())
                                            solver = new QMRSolver<double> (mat, pre);
                                          else
                                            solver = new QMRSolver<Complex> (mat, pre);                                            
                                          solver->SetPrecision(precision);
                                          solver->SetMaxSteps(maxsteps);
                                          solver->SetPrintRates (printrates);
                                          return shared_ptr<KrylovSpaceSolver>(solver);
                                        },
          "QMR Solver", py::arg("mat"), py::arg("pre"), py::arg("printrates")=true,
           py::arg("precision")=1e-8, py::arg("maxsteps")=200
          )
    ;
  
  m.def("ArnoldiSolver", [](BaseMatrix & mata, BaseMatrix & matm, shared_ptr<BitArray> freedofs,
                            py::list vecs, py::object bpshift)
        {
          if (py::len(vecs) > mata.Height())
            throw Exception ("number of eigenvectors to compute "+ToString(py::len(vecs))
                             + " is greater than matrix dimension "
                             + ToString(mata.Height()));
          if (mata.IsComplex())
            {
              Arnoldi<Complex> arnoldi (mata, matm, freedofs);
              Complex shift = 0.0;
              shift = py::cast<Complex>(bpshift);
              // cout << "shift = " << shift << endl;
              arnoldi.SetShift (shift);
              
              int nev = py::len(vecs);
              // cout << "num vecs: " << nev << endl;
              Array<shared_ptr<BaseVector>> evecs(nev);
                                                  
              Array<Complex> lam(nev);
              arnoldi.Calc (2*nev+1, lam, nev, evecs, 0);
              
              for (int i = 0; i < nev; i++)
                vecs[i].cast<BaseVector&>() = *evecs[i];
              
              Vector<Complex> vlam(nev);
              for (int i = 0; i < nev; i++)
                vlam(i) = lam[i];
              return vlam;
            }
          else
            {
              Arnoldi<double> arnoldi (mata, matm, freedofs);
              double shift = py::cast<double>(bpshift);
              // cout << "shift = " << shift << endl;
              arnoldi.SetShift (shift);
              
              int nev = py::len(vecs);
              // cout << "num vecs: " << nev << endl;
              Array<shared_ptr<BaseVector>> evecs(nev);
              
              Array<Complex> lam(nev);
              arnoldi.Calc (2*nev+1, lam, nev, evecs, 0);
              
              for (int i = 0; i < nev; i++)
                vecs[i].cast<BaseVector&>() = *evecs[i];
              
              Vector<Complex> vlam(nev);
              for (int i = 0; i < nev; i++)
                vlam(i) = lam[i];
              return vlam;
            }
        },
        "Arnoldi Solver", py::arg("mata"), py::arg("matm"), py::arg("freedofs"), py::arg("vecs"), py::arg("shift")=DummyArgument()
        )
    ;
  
  

  m.def("DoArchive" , [](shared_ptr<Archive> & arch, BaseMatrix & mat)
                                         { cout << "output basematrix" << endl;
                                           mat.DoArchive(*arch); return arch; });
                                           
}



PYBIND11_MODULE(libngla, m) {
  m.attr("__name__") = "la";
  ExportNgla(m);
}




#endif // NGS_PYTHON
