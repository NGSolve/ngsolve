#ifdef NGS_PYTHON
#include <la.hpp>
#include <parallelngs.hpp>
#include "../ngstd/python_ngstd.hpp"
using namespace ngla;
// include netgen-header to get access to PyMPI
#include <myadt.hpp>

#include "python_linalg.hpp"

template <typename T>
Array<T> ArrayFromVector (const std::vector<T> & vec)
{
  Array<T> a (vec.size());
  for (int i = 0; i < a.Size(); i++)
    a[i] = vec[i];
  return a;
}


class PyLinearOperator : public BaseMatrix
{
protected:
  py::object pyop;
  size_t h, w;
  bool is_complex;
public:
  PyLinearOperator (py::object apyop)
    : pyop(apyop)
  {
    py::object shape = pyop.attr("shape");
    h = py::cast<size_t> (shape.attr("__getitem__")(0));
    w = py::cast<size_t> (shape.attr("__getitem__")(1));
    // TODO: check for complex
    // auto dtype = py::cast<string> (pyop.attr("dtype"));
    // cout << "dtype = " << dtype << endl;
    // py::print (pyop.attr("dtype"));
    is_complex=false;
  }

  bool IsComplex() const override { return is_complex; }
  int VHeight() const override { return h; }
  int VWidth() const override { return w; }
  AutoVector CreateRowVector () const override { return CreateBaseVector(w, is_complex, 1); }
  AutoVector CreateColVector () const override { return CreateBaseVector(h, is_complex, 1); }

  void Mult (const BaseVector & x, BaseVector & y) const override
  {
    shared_ptr<BaseVector> spx(const_cast<BaseVector*>(&x), &NOOP_Deleter);
    py::object pyy = pyop * py::cast(spx);
    auto pya = py::array_t<double> (pyy);
    auto vec = pya. template unchecked<1>();
    VFlatVector<const double> vv(vec.size(), &vec(0));
    y = vv;
  }
};


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
         }, py::arg("pos"), "Return value at given position")
    .def("__setitem__",
         [](SparseMatrix<T> & self, py::tuple t, T value)
         {
           size_t row = t[0].cast<size_t>();
           size_t col = t[1].cast<size_t>();
           self(row,col) = value;
         }, py::arg("pos"), py::arg("value"), "Set value at given position")

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
    
    .def("CSR", [] (shared_ptr<SparseMatrix<T>> sp) -> py::object
         {
           // FlatArray<int> colind(sp->NZE(), sp->GetRowIndices().Addr(0));
           FlatArray<int> colind = sp->GetColIndices();
           // FlatVector<T> values(sp->NZE(), sp->GetRowValues().Addr(0));
           FlatVector<T> values = sp->GetValues();
           typedef typename mat_traits<T>::TSCAL TSCAL;
           FlatVector<TSCAL> svalues (values.Size()*sizeof(T)/sizeof(TSCAL), (TSCAL*)(void*)values.Addr(0));
           FlatArray<size_t> first = sp->GetFirstArray();
           if (sp->NZE() != colind.Size() || sp->NZE() != values.Size())
             {
               cout << "sizes don't match:" << endl
                    << "nze = " << sp->NZE() << endl
                    << "val.size = " << values.Size() << endl
                    << "colind.size = " << colind.Size() << endl;
             }
           return py::make_tuple (svalues, colind, first); 
         },
         py::return_value_policy::reference_internal)

    .def_property_readonly("entrysizes", [](shared_ptr<SparseMatrix<T>> self)
                           { return self->EntrySizes(); })
    
    .def_static("CreateFromCOO",
                [] (py::list indi, py::list indj, py::list values, size_t h, size_t w)
                {
                  auto cindi = makeCArray<int>(indi);
                  auto cindj = makeCArray<int>(indj);
                  auto cvalues = makeCArray<double>(values);
                  return SparseMatrix<double>::CreateFromCOO (cindi,cindj,cvalues, h,w);
                }, py::arg("indi"), py::arg("indj"), py::arg("values"), py::arg("h"), py::arg("w"))

    .def_static("CreateFromElmat",
                [] (py::list coldnums, py::list rowdnums, py::list elmats, size_t h, size_t w)
                {
                  auto cdnums1 = makeCTable<int>(coldnums);
                  auto rdnums1 = makeCTable<int>(rowdnums);
                  auto sparsemat = make_shared<SparseMatrix<double>>(h, w, cdnums1, rdnums1, false);
                  sparsemat->SetZero();
                  auto cdnums = makeCTable<int>(coldnums);
                  auto rdnums = makeCTable<int>(rowdnums);
                  for (int i = 0; i < py::len(elmats); i++)
                    {
                      const Matrix<double> & m = py::cast<const Matrix<double>&> (elmats[i]);
                      sparsemat -> AddElementMatrix(cdnums[i], rdnums[i], m, false);
                    }
                  return sparsemat;
                  // auto cvalues = makeCArray<double>(values);
                  // return SparseMatrix<double>::CreateFromCOO (cindi,cindj,cvalues, h,w);
                }, py::arg("col_ind"), py::arg("row_ind"), py::arg("matrices"), py::arg("h"), py::arg("w"))
    
    .def("CreateTranspose", [] (const SparseMatrix<T> & sp)
         { return sp.CreateTranspose (); }, "Return transposed matrix")

    .def("__matmul__", [] (const SparseMatrix<double> & a, const SparseMatrix<double> & b)
         { return MatMult(a,b); }, py::arg("mat"))
    .def("__matmul__", [](shared_ptr<SparseMatrix<T>> a, shared_ptr<BaseMatrix> mb)
         ->shared_ptr<BaseMatrix> { return make_shared<ProductMatrix> (a, mb); }, py::arg("mat"))
    ;

  py::class_<SparseMatrixSymmetric<T>, shared_ptr<SparseMatrixSymmetric<T>>, SparseMatrix<T>>
    (m, (string("SparseMatrixSymmetric") + typeid(T).name()).c_str());
}

void NGS_DLL_HEADER ExportNgla(py::module &m) {

  py::enum_<PARALLEL_STATUS>(m, "PARALLEL_STATUS", "enum of possible parallel statuses")
    .value("DISTRIBUTED", DISTRIBUTED)
    .value("CUMULATED", CUMULATED)
    .value("NOT_PARALLEL", NOT_PARALLEL)
    ;


  py::class_<ParallelDofs, shared_ptr<ParallelDofs>> (m, "ParallelDofs")
#ifdef PARALLEL
    .def("SubSet", [](const ParallelDofs & self, shared_ptr<BitArray> take_dofs) { 
        return self.SubSet(take_dofs); }, py::arg("dofs"))
    .def(py::init([](py::object procs, NgMPI_Comm comm) {
	  size_t n = py::len(procs);
	  TableCreator<int> ct(n);
	  while (!ct.Done()) {
	    size_t rn = 0;
	    for (auto row:procs) {
	      for (auto v:row)
		ct.Add(rn,v.cast<int>()); 
	      rn++;
	    }
	    ct++;
	  }
	  return new ParallelDofs(comm, ct.MoveTable());
	}), py::arg("dist_procs"), py::arg("comm"))
#endif
    .def_property_readonly ("comm", [](const ParallelDofs & self) { return self.GetCommunicator(); })
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
         { return self.GetDistantProcs(dof); }, py::arg("dof"))
    .def("Proc2Dof", [] (const ParallelDofs & self, int proc)
         { return self.GetExchangeDofs(proc); }, py::arg("proc"))
    .def("EnumerateGlobally", [] (shared_ptr<ParallelDofs> pardofs, shared_ptr<BitArray> freedofs) {
        Array<int> globnum;
        int num_glob_dofs;
        pardofs->EnumerateGlobally (freedofs, globnum, num_glob_dofs);
        return tuple ( py::cast(globnum), py::cast(num_glob_dofs) );
      }, py::arg("freedofs")=nullptr)
    .def("MasterDofs", &ParallelDofs::MasterDofs)
    .def_property_readonly("entrysize", [](shared_ptr<ParallelDofs> self)  { return self->GetEntrySize(); })
    ;

  py::class_<DofRange, IntRange> (m, "DofRange")
    ;

    m.def("CreateVVector",
          [] (size_t s, bool is_complex, int es) -> shared_ptr<BaseVector>
          { return CreateBaseVector(s,is_complex, es); },
          "size"_a, "complex"_a=false, "entrysize"_a=1);

    m.def("CreateParallelVector",
          [] (shared_ptr<ParallelDofs> pardofs, PARALLEL_STATUS status) -> shared_ptr<BaseVector>
          {
            return CreateParallelVector(pardofs, status);
            /*
#ifdef PARALLEL
	    if(pardofs->IsComplex())
	      return make_shared<S_ParallelBaseVectorPtr<Complex>> (pardofs->GetNDofLocal(), pardofs->GetEntrySize(), pardofs, DISTRIBUTED);
	    else
	      return make_shared<S_ParallelBaseVectorPtr<double>> (pardofs->GetNDofLocal(), pardofs->GetEntrySize(), pardofs, DISTRIBUTED);
#else
	    if(pardofs->IsComplex())
	      return make_shared<VVector<Complex>>(pardofs->GetNDofLocal());
	    else
	      return make_shared<VVector<double>>(pardofs->GetNDofLocal());
#endif
            */
	  },
          py::arg("pardofs"), py::arg("status"));
     
    
  py::class_<BaseVector, shared_ptr<BaseVector>>(m, "BaseVector",
        py::dynamic_attr() // add dynamic attributes
      )
    .def(py::init([] (size_t s, bool is_complex, int es) -> shared_ptr<BaseVector>
                  { return CreateBaseVector(s,is_complex, es); }),
         "size"_a, "complex"_a=false, "entrysize"_a=1)
    .def(py::init([] (DynamicVectorExpression expr)
                  { cout << IM(5) << "experimental: vector from expression" << endl;
                    return shared_ptr<BaseVector> (expr.Evaluate()); }))
    .def(py::init([] (py::array_t<double> bvec)
                  { 
                    auto vec = bvec. template unchecked<1>();
                    shared_ptr<BaseVector> bv = CreateBaseVector(vec.size(), false, 1);
                    FlatVector<double> fv = bv->FV<double>();
                    for (size_t i = 0; i < fv.Size(); i++)
                      fv(i) = vec(i);
                    return bv;
                  }))
    .def_property_readonly("local_vec", [](shared_ptr<BaseVector> self) -> shared_ptr<BaseVector> {
#ifdef PARALLEL
	auto pv = dynamic_cast_ParallelBaseVector (self.get());
	return (pv==nullptr) ? self : pv->GetLocalVector();
#else
	return self;
#endif
      } )
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
    .def_property_readonly("is_complex", &BaseVector::IsComplex)

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
                                   }, py::arg("vec"), py::arg("value"))
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
                                   }, py::arg("vec"), py::arg("value"))


    // TODO
//     .add_property("expr", py::object(expr_namespace["VecExpr"]) )
//     .add_property("data", py::object(expr_namespace["VecExpr"]), py::object(expr_namespace["expr_data"] ))
//     .def("__add__" , py::object(expr_namespace["expr_add"]) )
//     .def("__sub__" , py::object(expr_namespace["expr_sub"]) )
//     .def("__rmul__" , py::object(expr_namespace["expr_rmul"]) )
    .def("__getitem__",
          [](BaseVector & self,  int ind )
           {
             if (ind < 0) ind+=self.Size();
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
           }, py::arg("ind"), "Return value at given position" )
    .def("__getitem__", [](BaseVector & self,  py::slice inds )
      {
          size_t start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if (step != 1)
            throw Exception ("slices with non-unit distance not allowed");
          return shared_ptr<BaseVector>(self.Range(start, start+n));
      }, py::arg("inds"), "Return values at given position" )
    .def("__getitem__", [](BaseVector& self, IntRange range)
         {
           return shared_ptr<BaseVector>(self.Range(range));
         })
    .def("__getitem__", [](BaseVector& self, DofRange range)
         {
           return shared_ptr<BaseVector>(self.Range(range));
         })
    .def("__setitem__", [](BaseVector & self,  int ind, double d )
      {
          if (ind < 0) ind+=self.Size();
          if (ind < 0 || ind >= self.Size())
            throw py::index_error();
          self.Range(ind,ind+1) = d;
      }, py::arg("ind"), py::arg("value"), "Set value at given position" )
    .def("__setitem__", [](BaseVector & self,  int ind, Complex z )
      {
          if (ind < 0) ind+=self.Size();
          if (ind < 0 || ind >= self.Size())
            throw py::index_error();
        self.Range(ind,ind+1) = z;
      }, py::arg("ind"), py::arg("value"), "Set value at given position" )
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
      }, py::arg("inds"), py::arg("value"), "Set value at given positions" )
    .def("__setitem__", [](BaseVector & self,  py::slice inds, Complex z )
      {
          size_t start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if (step != 1)
            throw Exception ("slices with non-unit distance not allowed");          
          self.Range(start,start+n) = z;
      }, py::arg("inds"), py::arg("value"), "Set value at given positions" )
    .def("__setitem__", [](BaseVector & self,  IntRange range, double d )
      {
        self.Range(range) = d;
      }, py::arg("range"), py::arg("value"), "Set value for range of indices" )
    .def("__setitem__", [](BaseVector & self,  IntRange range, Complex z )
      {
        self.Range(range) = z;
      }, py::arg("range"), py::arg("value"), "Set value for range of indices" )
    .def("__setitem__", [](BaseVector & self, IntRange range, shared_ptr<BaseVector> v )
      {
        self.Range(range) = *v;
      }, py::arg("range"), py::arg("vec") )
    .def("__setitem__", [](BaseVector & self,  int ind, FlatVector<double> & v )
      {
          if (ind < 0) ind+=self.Size();
          if (ind < 0 || ind >= self.Size())
            throw py::index_error();
          if( self.IsComplex() )
            self.SV<Complex>()(ind) = v;
          else
            self.SV<double>()(ind) = v;
      }, py::arg("ind"), py::arg("vec") )
    .def("__setitem__", [](BaseVector & self,  int ind, FlatVector<Complex> & v )
      {
          if (ind < 0) ind+=self.Size();
          if (ind < 0 || ind >= self.Size())
            throw py::index_error();
          if( self.IsComplex() )
            self.SV<Complex>()(ind) = v;
          else
            throw py::index_error("cannot assign complex values to real vector");
      }, py::arg("ind"), py::arg("vec") )
    .def("__iadd__", [](BaseVector & self,  BaseVector & other) -> BaseVector& { self += other; return self;}, py::arg("vec"))
    .def("__isub__", [](BaseVector & self,  BaseVector & other) -> BaseVector& { self -= other; return self;}, py::arg("vec"))
    .def("__imul__", [](BaseVector & self,  double scal) -> BaseVector& { self *= scal; return self;}, py::arg("value"))
    .def("__imul__", [](BaseVector & self,  Complex scal) -> BaseVector& { self *= scal; return self;}, py::arg("value"))
    .def("__itruediv__", [](BaseVector & self,  double scal) -> BaseVector& { self /= scal; return self;}, py::arg("value"))
    .def("__itruediv__", [](BaseVector & self,  Complex scal) -> BaseVector& { self /= scal; return self;}, py::arg("value"))

    .def_property("data",
                  [](shared_ptr<BaseVector> self)
                  { return self; },
                  [](shared_ptr<BaseVector> self, DynamicVectorExpression v2)
                  { v2.AssignTo(1, *self); })
    
    .def("__add__", [] (shared_ptr<BaseVector> a, DynamicVectorExpression b)
         { return a+b; })
    .def("__iadd__", [] (shared_ptr<BaseVector> a, DynamicVectorExpression b) 
         { b.AddTo(1, *a); return a; })
    
    .def("__sub__", [] (shared_ptr<BaseVector> a, DynamicVectorExpression b)
         { return a-b; })
    .def("__isub__", [] (shared_ptr<BaseVector> a, DynamicVectorExpression b) 
         { b.AddTo(-1, *a); return a; })

    .def("__neg__", [] (shared_ptr<BaseVector> a) { return (-1.0)*a; })
    .def("__rmul__", [] (shared_ptr<BaseVector> a, double scal) { return scal*a; })
    .def("__rmul__", [] (shared_ptr<BaseVector> a, Complex scal) { return scal*a; })
    
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
                                          }, py::arg("other"), py::arg("conjugate")=py::cast(true), "Computes (complex) InnerProduct"         
         )
    .def("Norm",  [](BaseVector & self) { return self.L2Norm(); }, "Calculate Norm")
    .def("Range", [](BaseVector & self, int from, int to) -> shared_ptr<BaseVector>
                                   {
                                     return shared_ptr<BaseVector>(self.Range(from,to));
                                   }, py::arg("from"), py::arg("to"), "Return values from given range")
    .def("FV", [] (BaseVector & self)
                                {
                                  if (!self.IsComplex())
                                    return py::cast(self.FVDouble());
                                  else
                                    return py::cast(self.FVComplex());
                                })
    .def("Reshape", [] (BaseVector & self, size_t w)
         {
           size_t h = self.Size()/w;
           return FlatMatrix<> (h, w, &self.FVDouble()(0));
         }, py::arg("width"))
    .def("SetRandom", [] (BaseVector & self, optional<unsigned int> seed)
         {
           if (seed.has_value())
             srand(seed.value_or(0));
           self.SetRandom();
         }, py::arg("seed") = optional<unsigned int>())
    .def("Distribute", [] (BaseVector & self) { self.Distribute(); } ) 
    .def("Cumulate", [] (BaseVector & self) { self.Cumulate(); } ) 
    .def("GetParallelStatus", [] (BaseVector & self) { return self.GetParallelStatus(); } )
    .def("SetParallelStatus", [] (BaseVector & self, PARALLEL_STATUS stat) { self.SetParallelStatus(stat); }, py::arg("stat"));

  // m.def("InnerProduct",[](BaseVector & v1, BaseVector & v2)->double { return InnerProduct(v1,v2); })
  m.def ("InnerProduct",
         [] (py::object x, py::object y, py::kwargs kw) -> py::object
         { return py::handle(x.attr("InnerProduct")) (y, **kw); }, py::arg("x"), py::arg("y"), "Computes InnerProduct of given objects");
  ;
  

  py::class_<BlockVector, BaseVector, shared_ptr<BlockVector>> (m, "BlockVector")
    .def(py::init<> ([] (vector<shared_ptr<BaseVector>> vecs)
                     {
                       Array<shared_ptr<BaseVector>> v2;
                       for (auto v : vecs) v2 += v;
                       return make_shared<BlockVector> (v2);
                     }), py::arg("vecs"), "Makes BlockVector by given array of vectors")
    
    .def("__getitem__", [](BlockVector & self, int ind) { return self[ind]; }, py::arg("ind"), "Return block at given position")
    .def_property_readonly ("nblocks", [](const BlockVector & self) 
			    { return self.NBlocks(); },
                            "number of blocks in BlockVector")
    ;

  py::class_<MultiVectorExpr, shared_ptr<MultiVectorExpr> >(m, "MultiVectorExpr")
    .def("Scale", [](shared_ptr<MultiVectorExpr> e, Vector<double> v)
         -> shared_ptr<MultiVectorExpr>
         { return make_shared<ScaledMultiVectorExpr<double>> (e, v); })
    .def("__add__", [](shared_ptr<MultiVectorExpr> e1, shared_ptr<MultiVectorExpr> e2)
         { return e1+e2; })
    .def("__neg__", [](shared_ptr<MultiVectorExpr> e1)
         { return -e1; })
    .def("__sub__", [](shared_ptr<MultiVectorExpr> e1, shared_ptr<MultiVectorExpr> e2)
         { return e1-e2; })
    ;

  py::class_<MultiVector, MultiVectorExpr, shared_ptr<MultiVector>> (m, "MultiVector")
    // .def(py::init<shared_ptr<BaseVector>,size_t>([] ))
    .def(py::init<>([] (shared_ptr<BaseVector> bv, size_t cnt) { return bv->CreateMultiVector(cnt); } ))
    .def(py::init<size_t,size_t,bool>())
    .def("__len__", &MultiVector::Size)
    // .def("__getitem__", &MultiVector::operator[])
    .def("__getitem__",
         [](MultiVector & self, int ind )
         {
           if (ind < 0) ind+=self.Size();
           if (ind < 0 || ind >= self.Size()) 
             throw py::index_error();
           return self[ind];
         })
    .def("__getitem__", [](MultiVector & self, py::slice inds) {
        size_t start, step, n;
        InitSlice( inds, self.Size(), start, step, n );
        if (step != 1)
          throw Exception ("slices with non-unit distance not allowed");
        return shared_ptr<MultiVector>(self.Range(IntRange(start, start+n)));
      })
    /*
    .def("__getitem__", [](MultiVector & self, std::vector<int> inds) {
        return shared_ptr<MultiVector>(self.SubSet(ArrayFromVector(inds)));
      })
    */
    .def("__getitem__", [](MultiVector & self, Array<int> inds) {
        return shared_ptr<MultiVector>(self.SubSet(inds));
      })
    .def("__setitem__", [](MultiVector & x, int nr, DynamicVectorExpression & expr)
        {
          if( !x.IsComplex() )
              expr.AssignTo (1, *x[nr]);
          else
              expr.AssignTo (Complex(1), *x[nr]);
        })
    .def("__setitem__", [](MultiVector & x, int nr, double val)
        {
          *x[nr] = val;
        })
    .def("__setitem__", [](MultiVector & x, int nr, Complex val)
        {
          *x[nr] = val;
        })
    
    .def("__setitem__", [](MultiVector & self, py::slice inds, MultiVector & y) {
          size_t start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if (step != 1)
            throw Exception ("slices with non-unit distance not allowed");
          *self.Range(IntRange(start, start+n)) = y;
      })
    .def("__setitem__", [](MultiVector & self, py::slice inds, const MultiVectorExpr & v2) {
          size_t start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if (step != 1)
            throw Exception ("slices with non-unit distance not allowed");
          auto selfr = self.Range(IntRange(start, start+n));
          Vector<> ones(n); ones = 1;
          v2.AssignTo (ones, *selfr);
      })
    .def("__setitem__", [](MultiVector & self, std::vector<int> inds, MultiVector & y) {
        return *self.SubSet(ArrayFromVector(inds));
      })
    .def("__setitem__", [](MultiVector & self, std::vector<int> inds, const MultiVectorExpr & v2) {
        auto selfr = self.SubSet(ArrayFromVector(inds));
        *selfr = v2; })

    
    .def_property("data",
                  [](shared_ptr<MultiVector> self)
                  { return self; },
                  [](shared_ptr<MultiVector> self, const MultiVectorExpr & v2)
                  { *self = v2; })
    .def("Expand", &MultiVector::Extend, "deprecated, use Extend instead")
    .def("Extend", &MultiVector::Extend)
    .def("Append", &MultiVector::Append)
    .def("AppendOrthogonalize", &MultiVector::AppendOrthogonalize,
         py::arg("vec"), py::arg("ipmat")=nullptr,
         "assumes that existing vectors are orthogonal, and orthogonalize new vector against existing vectors")
    .def("Replace", [](MultiVector & x, int ind, shared_ptr<BaseVector> v2)
         {
           x.Replace(ind, v2);
         }, py::arg("ind"), py::arg("v2"))
    .def("Replace", [](MultiVector & x, std::vector<int> inds, MultiVector & v2)
         {
           for (size_t i = 0; i < inds.size(); i++)
             x.Replace(inds[i], v2[i]);
         }, py::arg("inds"), py::arg("mv2"))
    .def("InnerProduct", [](MultiVector & x, MultiVector & y, bool conjugate)
        { 
          if( !x.IsComplex() )
            return py::cast(x.InnerProductD(y));
          else
            return py::cast(x.InnerProductC(y, conjugate));
        }, py::arg("other"), py::arg("conjugate")=py::cast(true))
    .def("InnerProduct", [](MultiVector & x, MultiVectorExpr & y, bool conjugate)
        { 
          if( !x.IsComplex() )
            return py::cast(x.InnerProductD(y));
          else
            return py::cast(x.InnerProductC(y, conjugate));
        }, py::arg("other"), py::arg("conjugate")=py::cast(true))
    /*
    .def("Orthogonalize", [](MultiVector & x, BaseMatrix * ipmat)
         {
           x.Orthogonalize(ipmat);
         },py::arg("ipmat")=nullptr)
    */
    .def("Orthogonalize",[](MultiVector & x, BaseMatrix * ipmat)
         {
           if (!x.IsComplex())
             return py::cast(x.T_Orthogonalize<double>(ipmat));
           else
             return py::cast(x.T_Orthogonalize<Complex>(ipmat));
         }, py::arg("ipmat")=nullptr,
         "Orthogonalize vectors by modified Gram-Schmidt, returns R-factor of QR decomposition (only ipmat version, for the moment)")
    .def("__mul__", [](shared_ptr<MultiVector> x, Vector<double> a) 
         { // cout << "in double __mul__" << endl;
           return DynamicVectorExpression(make_shared<MultiVecAxpyExpr<double>>(a, x)); })
    .def("__mul__", [](shared_ptr<MultiVector> x, Vector<Complex> a) 
         { // cout << "in complex __mul__" << endl;
           return DynamicVectorExpression(make_shared<MultiVecAxpyExpr<Complex>>(a, x)); })

    .def("__mul__", [](shared_ptr<MultiVector> x, Matrix<double> a)
         -> shared_ptr<MultiVectorExpr>
         {
           return make_shared<MultiVecMatrixExpr<double>>(a, x);
         })
    .def("__mul__", [](shared_ptr<MultiVector> x, Matrix<Complex> a)
         -> shared_ptr<MultiVectorExpr>
         { return make_shared<MultiVecMatrixExpr<Complex>>(a, x); })
    ;
    /*
    // not taken, thus moved to BaseMatrix
    .def("__rmul__", [](shared_ptr<MultiVector> x, shared_ptr<BaseMatrix> mat)
         -> shared_ptr<MultiVectorExpr>
         {
           return make_shared<MatMultiVecExpr> (mat, x); 
         })
    */
  

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
      
      AutoVector CreateRowVector () const override {
        py::gil_scoped_acquire gil;
        pybind11::function overload = pybind11::get_overload(this, "CreateRowVector");
        if (overload) {
          auto vec = py::cast<shared_ptr<BaseVector>> (overload());
          return vec->CreateVector();
        }
        throw Exception ("CreateRowVector not overloaded from python");        
        // python can only create shared_ptr<BaseVector>, create another unique_ptr<vector> from C++
#ifdef NONE
        PYBIND11_OVERLOAD_PURE(
           shared_ptr<BaseVector>, /* Return type */
           BaseMatrix,             /* Parent class */
           CreateRowVector,        /* Name of function */
           );
#endif
      }

      AutoVector CreateColVector () const override {
        py::gil_scoped_acquire gil;
        pybind11::function overload = pybind11::get_overload(this, "CreateColVector");
        if (overload) {
          auto vec = py::cast<shared_ptr<BaseVector>> (overload());
          return vec->CreateVector();
        }
        throw Exception ("CreateColVector not overloaded from python");        
#ifdef NONE        
        PYBIND11_OVERLOAD_PURE(
           shared_ptr<BaseVector>, /* Return type */
           BaseMatrix,             /* Parent class */
           CreateColVector,        /* Name of function */
           );
#endif
      }

#ifdef NONE                
      AutoVector CreateVector () const override {
        PYBIND11_OVERLOAD_PURE(
           shared_ptr<BaseVector>, /* Return type */
           BaseMatrix,             /* Parent class */
           CreateVector,           /* Name of function */
           );
      }
#endif

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

      void MultTrans (const BaseVector & x, BaseVector & y) const override {
        pybind11::gil_scoped_acquire gil;
        pybind11::function overload = pybind11::get_overload(this, "MultTrans");
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
          BaseMatrix::MultTrans(x,y);
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
    .def(py::init<> ())
    .def(py::init<>([] (py::object pyob)
                    { return make_shared<PyLinearOperator> (pyob); }))
    .def("__str__", [](BaseMatrix &self) { return ToString<BaseMatrix>(self); } )
    .def_property_readonly("height", [] ( BaseMatrix & self)
                           { return self.Height(); }, "Height of the matrix" )
    .def_property_readonly("width", [] ( BaseMatrix & self)
                           { return self.Width(); }, "Width of the matrix" )
    .def_property_readonly("nze", [] ( BaseMatrix & self)
                           { return self.NZE(); }, "number of non-zero elements")
    .def_property_readonly("local_mat", [](shared_ptr<BaseMatrix> & mat) { return mat; })
    // .def("CreateMatrix", &BaseMatrix::CreateMatrix)

    .def("GetOperatorInfo", [] (BaseMatrix & self)
         {
           stringstream str;
           self.PrintOperatorInfo(str);
           return str.str();
         })
    
    .def("CreateMatrix", [] ( BaseMatrix & self)
         { return self.CreateMatrix(); }, "Create matrix of same dimension and same sparsestructure" )

    
    .def("CreateRowVector", [] ( BaseMatrix & self)
        { return shared_ptr<BaseVector>(self.CreateRowVector()); } )
    .def("CreateColVector", [] ( BaseMatrix & self)
        { return shared_ptr<BaseVector>(self.CreateColVector()); } )
    
    .def("AsVector", [] (BM & m)
                                      {
                                        return shared_ptr<BaseVector> (&m.AsVector(), NOOP_Deleter);
                                      }, "Interprets the matrix values as a vector")

    .def("Mult",         [](BaseMatrix &m, BaseVector &x, BaseVector &y) { m.Mult(x, y); }, py::call_guard<py::gil_scoped_release>(), py::arg("x"), py::arg("y"))
    .def("MultAdd",      [](BaseMatrix &m, double s, BaseVector &x, BaseVector &y) { m.MultAdd (s, x, y); }, py::arg("value"), py::arg("x"), py::arg("y"), py::call_guard<py::gil_scoped_release>())
    .def("MultTrans",    [](BaseMatrix &m, double s, BaseVector &x, BaseVector &y) { y=0; m.MultTransAdd (1.0, x, y); }, py::arg("value"), py::arg("x"), py::arg("y"), py::call_guard<py::gil_scoped_release>())
    .def("MultTransAdd",  [](BaseMatrix &m, double s, BaseVector &x, BaseVector &y) { m.MultTransAdd (s, x, y); }, py::arg("value"), py::arg("x"), py::arg("y"), py::call_guard<py::gil_scoped_release>())
    .def("MultScale",    [](BaseMatrix &m, double s, BaseVector &x, BaseVector &y)
          {
              m.Mult (x,y);
              if(s!=1.0)
                  y *= s;
          }, py::arg("value"), py::arg("x"), py::arg("y") , py::call_guard<py::gil_scoped_release>())
    .def("MultAdd",      [](BaseMatrix &m, Complex s, BaseVector &x, BaseVector &y) { m.MultAdd (s, x, y); }, py::arg("value"), py::arg("x"), py::arg("y"), py::call_guard<py::gil_scoped_release>())
    .def("MultTrans",    [](BaseMatrix &m, Complex s, BaseVector &x, BaseVector &y) { y=0; m.MultTransAdd (1.0, x, y); }, py::arg("value"), py::arg("x"), py::arg("y"), py::call_guard<py::gil_scoped_release>())
    .def("MultTransAdd",  [](BaseMatrix &m, Complex s, BaseVector &x, BaseVector &y) { m.MultTransAdd (s, x, y); }, py::arg("value"), py::arg("x"), py::arg("y"), py::call_guard<py::gil_scoped_release>())
    .def("MultScale",    [](BaseMatrix &m, Complex s, BaseVector &x, BaseVector &y)
          {
              m.Mult (x,y);
              if(s!=1.0)
                  y *= s;
          }, py::arg("value"), py::arg("x"), py::arg("y"), py::call_guard<py::gil_scoped_release>() )

    .def("__iadd__", [] (BM &m, BM &m2) { 
        m.AsVector()+=m2.AsVector();
      }, py::arg("mat"), py::call_guard<py::gil_scoped_release>())

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
Parameters:

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

    .def_property_readonly("T", [](shared_ptr<BM> m)->shared_ptr<BaseMatrix>
                           { return TransposeOperator (m); }, "Return transpose of matrix")
    .def_property_readonly("H", [](shared_ptr<BM> m)->shared_ptr<BaseMatrix> { return make_shared<ConjTrans> (m); }, "Return conjugate transpose of matrix (WIP, only partially supported)")
    /*
    .def("__matmul__", [](shared_ptr<BM> ma, shared_ptr<EmbeddingTranspose> mb)->shared_ptr<BaseMatrix>
         { return make_shared<EmbeddedTransposeMatrix> (mb->Width(), mb->GetRange(), ma); }, py::arg("mat"))
    .def("__matmul__", [](shared_ptr<BM> ma, shared_ptr<BM> mb)->shared_ptr<BaseMatrix>
         { return make_shared<ProductMatrix> (ma, mb); }, py::arg("mat"))
    */
    .def("__matmul__", [](shared_ptr<BM> ma, shared_ptr<BM> mb)->shared_ptr<BaseMatrix>
         { return ComposeOperators (ma, mb); }, py::arg("mat"))
    
    .def("__add__", [](shared_ptr<BM> ma, shared_ptr<BM> mb)->shared_ptr<BaseMatrix>
         { return make_shared<SumMatrix> (ma, mb, 1, 1); }, py::arg("mat"))
    .def("__radd__", [](shared_ptr<BM> ma, int i) {
        if (i != 0) throw Exception("can only add integer 0 to BaseMatrix (for Python sum(list))");
        return ma; })
    .def("__sub__", [](shared_ptr<BM> ma, shared_ptr<BM> mb)->shared_ptr<BaseMatrix>
         { return make_shared<SumMatrix> (ma, mb, 1, -1); }, py::arg("mat"))
    .def("__rmul__", [](shared_ptr<BM> ma, double a)->shared_ptr<BaseMatrix>
         { return make_shared<VScaleMatrix<double>> (ma, a); }, py::arg("value"))
    .def("__rmul__", [](shared_ptr<BM> ma, Complex a)->shared_ptr<BaseMatrix>
         { return make_shared<VScaleMatrix<Complex>> (ma, a); }, py::arg("value"))
    .def("__neg__", [](shared_ptr<BM> ma)->shared_ptr<BaseMatrix>
         { return make_shared<VScaleMatrix<double>> (ma, -1); })
    .def("__mul__", [](shared_ptr<BaseMatrix> m, shared_ptr<BaseVector> v)
         { return DynamicVectorExpression(make_shared<DynamicMatVecExpression>(m,v)); })
    // TODO: solve complex problem
    .def("__mul__", [](shared_ptr<BaseMatrix> mat, shared_ptr<MultiVector> x)
         -> shared_ptr<MultiVectorExpr>
         {
           return make_shared<MatMultiVecExpr> (mat, x); 
         })
             
    
    .def("Update", [](BM &m) { m.Update(); }, py::call_guard<py::gil_scoped_release>(), "Update matrix")
    ;

  py::class_<BaseSparseMatrix, shared_ptr<BaseSparseMatrix>, BaseMatrix>
    (m, "BaseSparseMatrix", "sparse matrix of any type")
    
    .def("CreateSmoother", [](BaseSparseMatrix & m, shared_ptr<BitArray> ba) 
         { return m.CreateJacobiPrecond(ba); }, py::call_guard<py::gil_scoped_release>(),
         py::arg("freedofs") = shared_ptr<BitArray>())
    
    .def("CreateBlockSmoother", [](BaseSparseMatrix & m, py::object blocks, bool parallel)
         {
           shared_ptr<Table<int>> blocktable;
           {
             py::gil_scoped_acquire aq;
             size_t size = py::len(blocks);
           
             Array<int> cnt(size);
             size_t i = 0;
             for (auto block : blocks)
               cnt[i++] = py::len(block);
           
             i = 0;
             blocktable = make_shared<Table<int>>(cnt);
             for (auto block : blocks)
               {
                 auto row = (*blocktable)[i++];
                 size_t j = 0;
                 for (auto val : block)
                   row[j++] = val.cast<int>();
               }
           }
           return m.CreateBlockJacobiPrecond (blocktable, nullptr, parallel);
         }, py::call_guard<py::gil_scoped_release>(), py::arg("blocks"), py::arg("parallel")=false)
     ;
  
  py::class_<S_BaseMatrix<double>, shared_ptr<S_BaseMatrix<double>>, BaseMatrix>
    (m, "S_BaseMatrixD", "base sparse matrix");
  py::class_<S_BaseMatrix<Complex>, shared_ptr<S_BaseMatrix<Complex>>, BaseMatrix>
    (m, "S_BaseMatrixC", "base sparse matrix");

  py::class_<CumulationOperator, shared_ptr<CumulationOperator>, BaseMatrix> (m, "CumulationOperator")
    .def(py::init<shared_ptr<ParallelDofs>>())
    ;

  py::class_<LoggingMatrix, shared_ptr<LoggingMatrix>, BaseMatrix> (m, "LoggingMatrix")
    .def(py::init<shared_ptr<BaseMatrix>,string,string,optional<NgMPI_Comm>>(),
         py::arg("mat"), py::arg("label"), py::arg("logfile")="stdout", py::arg("comm")=std::nullopt)
    ;
  
  py::class_<ConstantElementByElementMatrix, shared_ptr<ConstantElementByElementMatrix>, BaseMatrix>
    (m, "ConstEBEMatrix")
    .def(py::init<> ([] (size_t h, size_t w, Matrix<> mat,
                         py::list pycdofs, py::list pyrdofs)
                     {
                       /*
                       size_t n = py::len(pyrdofs);
                       Array<int> entrysize(n);
                       for (size_t i = 0; i < n; i++)
                         entrysize[i] = py::len(pyrdofs[i]);
                       Table<int> rdofs(entrysize);
                       for (size_t i = 0; i < n; i++)
                         {
                           const py::object & obj = pyrdofs[i];
                           rdofs[i] = makeCArray<int> (obj);
                         }

                       for (size_t i = 0; i < n; i++)
                         entrysize[i] = py::len(pycdofs[i]);
                       Table<int> cdofs(entrysize);
                       for (size_t i = 0; i < n; i++)
                         {
                           const py::object & obj = pycdofs[i];
                           cdofs[i] = makeCArray<int> (obj);
                         }
                       */
                       auto rdofs = makeCTable<int> (pyrdofs);
                       auto cdofs = makeCTable<int> (pycdofs);
                       
                       return make_shared<ConstantElementByElementMatrix> (h, w, mat,
                                                                           std::move(cdofs), std::move(rdofs));
                     }),
         py::arg("h"), py::arg("w"), py::arg("matrix"),
         py::arg("col_ind"), py::arg("row_ind"))
    ;

  m.def("ChebyshevIteration", [](shared_ptr<BaseMatrix> mat, shared_ptr<BaseMatrix> pre,
				 int steps, double lambda_min, double lambda_max)
	-> shared_ptr<BaseMatrix> {
	  auto cheb = make_shared<ChebyshevIteration>(*mat, *pre, steps);
	  cheb->SetBounds(lambda_min, lambda_max);
	  return cheb;
	}, py::arg("mat") = nullptr, py::arg("pre") = nullptr,
	py::arg("steps") = 3, py::arg("lam_min") = 1, py::arg("lam_max") = 1);
  
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
                     }), py::arg("mats"), "Make BlockMatrix with given array of matrices")
    .def("__getitem__", [](BlockMatrix & self, py::tuple inds) { 
        if (py::len(inds) != 2)
          throw Exception ("BlockMatrix needs two indices to access block");

        int row = inds[0].cast<int>();
        int col = inds[1].cast<int>();
        return self(row,col); 
      }, py::arg("inds"), "Return value at given position")
    .def_property_readonly("row_nblocks", [](BlockMatrix & mat) { return mat.BlockRows(); })
    .def_property_readonly("col_nblocks", [](BlockMatrix & mat) { return mat.BlockCols(); })
    ;

  py::class_<DynamicVectorExpression> (m, "DynamicVectorExpression")
    .def(py::init<shared_ptr<BaseVector>>())
    .def(py::init([] (py::array_t<double> bvec)
                  {
                    auto vec = bvec. template unchecked<1>();
                    shared_ptr<BaseVector> bv = make_shared<VFlatVector<const double>> (vec.size(), &vec(0));
                    return DynamicVectorExpression(bv);
                  }), py::keep_alive<1,2>())
    .def(py::self+py::self)
    .def(py::self-py::self)
    .def("__neg__", [] (DynamicVectorExpression a) { return (-1.0)*a; })    
    .def(double()*py::self)
    .def("__rmul__", [] (DynamicVectorExpression a, Complex scal) { return scal*a; })

    /*
      // crashing, why ? 
    .def("InnerProduct", [](DynamicVectorExpression a, py::args la, py::kwargs kw)
         { auto evalv = a.Evaluate();
           reeturn py::cast(evalv).attr("InnerProduct")(la, kw);
         })
    */
    // only the real case, for the moment ..
    .def("InnerProduct", [](DynamicVectorExpression v1, BaseVector & v2)
         { auto v = v1.Evaluate();
           return InnerProduct<double> (v, v2); })

    .def("Evaluate", [](DynamicVectorExpression expr)
         {
           return shared_ptr<BaseVector> (expr.Evaluate());
         }, "create vector and evaluate expression into it")
  ;

  py::implicitly_convertible<BaseVector, DynamicVectorExpression>();
  py::implicitly_convertible<DynamicVectorExpression, BaseVector>();
  py::implicitly_convertible<py::array_t<double>, DynamicVectorExpression>();
  py::implicitly_convertible<py::array_t<double>, BaseVector>();
  
#ifndef PARALLEL

  m.def("ParallelMatrix", [](py::object mat, py::object row_pardofs, py::object col_pardofs, py::object op) {
      throw Exception("Sorry, ParallelMatrix only available in MPI version!");
    }, py::arg("mat")=py::none(), py::arg("row_pardofs")=py::none(), py::arg("col_pardofs")=py::none(), py::arg("op")=py::none());
  m.def("ParallelMatrix", [](py::object mat, py::object pardofs, py::object op) {
      throw Exception("Sorry, ParallelMatrix only available in MPI version!");
    }, py::arg("mat")=py::none(), py::arg("pardofs")=py::none(), py::arg("op")=py::none());
    
#else
  
  auto parmat = py::class_<ParallelMatrix, shared_ptr<ParallelMatrix>, BaseMatrix>
    (m, "ParallelMatrix", "MPI-distributed matrix");

  py::enum_<PARALLEL_OP>(parmat, "PARALLEL_OP", "enum of possible parallel ops")
    .value("C2C", C2C)
    .value("C2D", C2D)
    .value("D2C", D2C)
    .value("D2D", D2D)
    .export_values()
    ;

  parmat.def(py::init<shared_ptr<BaseMatrix>, shared_ptr<ParallelDofs>, PARALLEL_OP>(),
	     py::arg("mat"), py::arg("pardofs"), py::arg("op")=C2D)
    .def(py::init<shared_ptr<BaseMatrix>, shared_ptr<ParallelDofs>, shared_ptr<ParallelDofs>, PARALLEL_OP>(),
	 py::arg("mat"), py::arg("row_pardofs"), py::arg("col_pardofs"), py::arg("op")=C2D)
    .def_property_readonly("row_pardofs", [](ParallelMatrix & mat) { return mat.GetRowParallelDofs(); })
    .def_property_readonly("col_pardofs", [](ParallelMatrix & mat) { return mat.GetColParallelDofs(); })
    .def_property_readonly("local_mat", [](ParallelMatrix & mat) { return mat.GetMatrix(); })
    .def_property_readonly("op_type", [](ParallelMatrix & mat) { return mat.GetOpType(); })
    ;


  py::class_<FETI_Jump_Matrix, shared_ptr<FETI_Jump_Matrix>, BaseMatrix>
    (m, "FETI_Jump", "B-matrix of the FETI-system")
    .def(py::init<shared_ptr<ParallelDofs>>(),
	 py::arg("pardofs"))
    .def(py::init<shared_ptr<ParallelDofs>, shared_ptr<ParallelDofs>>(),
	 py::arg("pardofs"), py::arg("u_pardofs"))
    .def_property_readonly("row_pardofs", [](FETI_Jump_Matrix & mat) { return mat.GetRowParallelDofs(); })
    .def_property_readonly("col_pardofs", [](FETI_Jump_Matrix & mat) { return mat.GetColParallelDofs(); })
    ;
#endif
  

  
  ExportSparseMatrix<double>(m);
  ExportSparseMatrix<Complex>(m);
  ExportSparseMatrix<Mat<2,2,double>>(m);
  ExportSparseMatrix<Mat<2,2,Complex>>(m);
  ExportSparseMatrix<Mat<3,3,double>>(m);
  ExportSparseMatrix<Mat<3,3,Complex>>(m);


  py::class_<SparseMatrixDynamic<double>, shared_ptr<SparseMatrixDynamic<double>>, BaseMatrix>
    (m, "SparseMatrixDynamic")
    .def(py::init([] (const BaseMatrix & mat) -> shared_ptr<SparseMatrixDynamic<double>>
                  {
                    if (auto ptr = dynamic_cast<const SparseMatrixTM<double>*> (&mat); ptr)
                      return make_shared<SparseMatrixDynamic<double>> (*ptr);
                    if (auto ptr = dynamic_cast<const SparseMatrixTM<Mat<2,2>>*> (&mat); ptr)
                      return make_shared<SparseMatrixDynamic<double>> (*ptr);
                    if (auto ptr = dynamic_cast<const SparseMatrixTM<Mat<3,3>>*> (&mat); ptr)
                      return make_shared<SparseMatrixDynamic<double>> (*ptr);
#if MAX_SYS_DIM >= 4
                    if (auto ptr = dynamic_cast<const SparseMatrixTM<Mat<4,4>>*> (&mat); ptr)
                      return make_shared<SparseMatrixDynamic<double>> (*ptr);
#endif
#if MAX_SYS_DIM >= 5                    
                    if (auto ptr = dynamic_cast<const SparseMatrixTM<Mat<5,5>>*> (&mat); ptr)
                      return make_shared<SparseMatrixDynamic<double>> (*ptr);
#endif
#if MAX_SYS_DIM >= 6                    
                    if (auto ptr = dynamic_cast<const SparseMatrixTM<Mat<6,6>>*> (&mat); ptr)
                      return make_shared<SparseMatrixDynamic<double>> (*ptr);
#endif                    
                    return nullptr;
                  })
         )
         ;


  py::class_<SparseMatrixVariableBlocks<double>, shared_ptr<SparseMatrixVariableBlocks<double>>, BaseMatrix>
    (m, "SparseMatrixVariableBlocks")
    .def(py::init([] (const BaseMatrix & mat)
                  {
                    if (auto ptr = dynamic_cast<const SparseMatrixTM<double>*> (&mat); ptr)
                      return make_shared<SparseMatrixVariableBlocks<double>> (*ptr);
                    throw Exception("cannot create SparseMatrixVariableBlocks");
                  }))
    ;

  
  py::class_<BaseBlockJacobiPrecond, shared_ptr<BaseBlockJacobiPrecond>, BaseMatrix>
    (m, "BlockSmoother",
     "block Jacobi and block Gauss-Seidel smoothing")
    .def("Smooth", &BaseBlockJacobiPrecond::GSSmooth, py::call_guard<py::gil_scoped_release>(),
         py::arg("x"), py::arg("b"), py::arg("steps")=1,
         "performs steps block-Gauss-Seidel iterations for the linear system A x = b")
    .def("SmoothBack", &BaseBlockJacobiPrecond::GSSmoothBack,
         py::arg("x"), py::arg("b"), py::arg("steps")=1, py::call_guard<py::gil_scoped_release>(),
         "performs steps block-Gauss-Seidel iterations for the linear system A x = b in reverse order")
    ;

  py::class_<BaseJacobiPrecond, shared_ptr<BaseJacobiPrecond>, BaseMatrix>
    (m, "Smoother",
     "Jacobi and Gauss-Seidel smoothing")
    .def("Smooth", [&](BaseJacobiPrecond & jac, BaseVector & x, BaseVector & b)
         { jac.GSSmooth (x, b); }, py::call_guard<py::gil_scoped_release>(),
         py::arg("x"), py::arg("b"),
         "performs one step Gauss-Seidel iteration for the linear system A x = b")
    .def("SmoothBack", &BaseJacobiPrecond::GSSmoothBack,
         py::arg("x"), py::arg("b"), py::call_guard<py::gil_scoped_release>(),
         "performs one step Gauss-Seidel iteration for the linear system A x = b in reverse order")
    ;

  py::class_<SparseFactorization, shared_ptr<SparseFactorization>, BaseMatrix>
    (m, "SparseFactorization")
    .def("Smooth", [] (SparseFactorization & self, BaseVector & u, BaseVector & y)
         {
           self.Smooth (u, y /* this is not needed */, y);
         }, py::call_guard<py::gil_scoped_release>(),
         "perform smoothing step (needs non-symmetric storage so symmetric sparse matrix)")
    ;

  py::class_<SparseCholesky<double>, shared_ptr<SparseCholesky<double>>, SparseFactorization> (m, "SparseCholesky_d");
  py::class_<SparseCholesky<Complex>, shared_ptr<SparseCholesky<Complex>>, SparseFactorization> (m, "SparseCholesky_c");
  
  py::class_<Projector, shared_ptr<Projector>, BaseMatrix> (m, "Projector")
    .def(py::init<shared_ptr<BitArray>,bool>(),
         py::arg("mask"), py::arg("range"),
         "Linear operator projecting to true/false bits of BitArray mask, depending on argument range")
    .def("Project", &Projector::Project, "project vector inline")
    ;
  
  py::class_<ngla::IdentityMatrix, shared_ptr<ngla::IdentityMatrix>, BaseMatrix> (m, "IdentityMatrix")
    .def(py::init<>())
    .def(py::init<size_t, bool>(),
         py::arg("size"), py::arg("complex")=false)
    ;

  py::class_<Real2ComplexMatrix<double,Complex>, shared_ptr<Real2ComplexMatrix<double,Complex>>,
             BaseMatrix> (m, "Real2ComplexMatrix")
    .def(py::init<shared_ptr<BaseMatrix>>())
    ;
  
  py::class_<PermutationMatrix, shared_ptr<PermutationMatrix>, BaseMatrix> (m, "PermutationMatrix")
    .def(py::init([](size_t w, std::vector<size_t> ind)
                  {
                    Array<size_t> inda(ind.size());
                    for (size_t i = 0; i < inda.Size(); i++)
                      inda[i] = ind[i];
                    return make_shared<PermutationMatrix> (w, move(inda)); 
                  }),
         py::arg("w"), py::arg("ind"))
    ;
  
  py::class_<Embedding, shared_ptr<Embedding>, BaseMatrix> (m, "Embedding")
    .def(py::init<size_t, IntRange, bool>(),
         py::arg("height"), py::arg("range"), py::arg("complex")=false,
         "Linear operator embedding a shorter vector into a longer vector")
    /*
    .def_property_readonly("T", [](shared_ptr<Embedding> m)->shared_ptr<EmbeddingTranspose>
                           { return make_shared<EmbeddingTranspose> (m->Height(), m->GetRange()); }, "Return transpose of matrix")
    */
    /*
    .def("__matmul__", [](shared_ptr<Embedding> ma, shared_ptr<BM> mb)->shared_ptr<BaseMatrix>
         { return make_shared<EmbeddedMatrix> (ma->Height(), ma->GetRange(), mb); }, py::arg("mat"))
    */
    .def("__matmul__", [](shared_ptr<Embedding> ma, shared_ptr<BM> mb)
         { return ComposeOperators(ma, mb); }, py::arg("mat"))
    
    ;
  
  py::class_<EmbeddingTranspose, shared_ptr<EmbeddingTranspose>, BaseMatrix> (m, "EmbeddingTranspose")
    /*
    .def("__rmatmul__", [](shared_ptr<EmbeddingTranspose> ma, shared_ptr<BM> mb)->shared_ptr<BaseMatrix>
         { return make_shared<EmbeddedTransposeMatrix> (ma->Width(), ma->GetRange(), mb); }, py::arg("mat"))
    */
    .def("__rmatmul__", [](shared_ptr<EmbeddingTranspose> ma, shared_ptr<BM> mb)
         { return ComposeOperators(mb, ma); }, py::arg("mat"))
    
    ;
    
  py::class_<KrylovSpaceSolver, shared_ptr<KrylovSpaceSolver>, BaseMatrix> (m, "KrylovSpaceSolver")
    .def("GetSteps", &KrylovSpaceSolver::GetSteps)
    .def_property("tol", &KrylovSpaceSolver::GetPrecision, &KrylovSpaceSolver::SetPrecision)
    .def_property("maxsteps", &KrylovSpaceSolver::GetMaxSteps, &KrylovSpaceSolver::SetMaxSteps)
    .def("SetAbsolutePrecision", &KrylovSpaceSolver::SetAbsolutePrecision)
    ;

  m.def("CGSolver", [](shared_ptr<BaseMatrix> mat, shared_ptr<BaseMatrix> pre,
                       bool iscomplex, bool printrates,
                       double precision, int maxsteps, bool conjugate)
                                       {
                                         shared_ptr<KrylovSpaceSolver> solver;
                                         if(mat->IsComplex()) iscomplex = true;
                                         
                                         if (iscomplex)
                                           {
                                             if(conjugate)
                                               solver = make_shared<CGSolver<ComplexConjugate>>(mat, pre);
                                             else
                                               solver = make_shared<CGSolver<Complex>> (mat, pre);
                                           }
                                         else
                                           solver = make_shared<CGSolver<double>> (mat, pre);
                                         solver->SetPrecision(precision);
                                         solver->SetMaxSteps(maxsteps);
                                         solver->SetPrintRates (printrates);
                                         return solver;
                                       },
           py::arg("mat"), py::arg("pre"), py::arg("complex") = false, py::arg("printrates")=true,
        py::arg("precision")=1e-8, py::arg("maxsteps")=200, py::arg("conjugate")=false,
        docu_string(R"raw_string(
A CG Solver.

Parameters:

mat : ngsolve.la.BaseMatrix
  input matrix 

pre : ngsolve.la.BaseMatrix
  input preconditioner matrix

complex : bool
  input complex, if not set it is deduced from matrix type

printrates : bool
  input printrates

precision : float
  input requested precision. CGSolver stops if precision is reached.

maxsteps : int
  input maximal steps. CGSolver stops after this steps.

)raw_string"))
    ;

  m.def("GMRESSolver", [](shared_ptr<BaseMatrix> mat, shared_ptr<BaseMatrix> pre,
                          bool printrates, 
                          double precision, int maxsteps)
        {
          shared_ptr<KrylovSpaceSolver> solver;
          if (!mat->IsComplex())
            solver = make_shared<GMRESSolver<double>> (mat, pre);
          else
            solver = make_shared<GMRESSolver<Complex>> (mat, pre);                                            
          solver->SetPrecision(precision);
          solver->SetMaxSteps(maxsteps);
          solver->SetPrintRates (printrates);
          return shared_ptr<KrylovSpaceSolver>(solver);
        },
        py::arg("mat"), py::arg("pre"), py::arg("printrates")=true,
        py::arg("precision")=1e-8, py::arg("maxsteps")=200, docu_string(R"raw_string(
A General Minimal Residuum (GMRES) Solver.

Parameters:

mat : ngsolve.la.BaseMatrix
  input matrix 

pre : ngsolve.la.BaseMatrix
  input preconditioner matrix

printrates : bool
  input printrates

precision : float
  input requested precision. GMRESSolver stops if precision is reached.

maxsteps : int
  input maximal steps. GMRESSolver stops after this steps.

)raw_string"))
    ;

  m.def("EigenValues_Preconditioner", [](const BaseMatrix & mat, const BaseMatrix & pre, double tol) {
      EigenSystem eigen(mat, pre);
      eigen.SetPrecision(tol);
      eigen.Calc();
      Vector<double> ev(eigen.NumEigenValues());
      for (size_t i = 0; i < ev.Size(); i++)
        ev[i] = eigen.EigenValue(i+1);
      return ev;
    },
    py::arg("mat"), py::arg("pre"), py::arg("tol")=1e-10,
    "Calculate eigenvalues of pre * mat, where pre and mat are positive definite matrices.\n"
    "The typical usecase of this function is to calculate the condition number of a preconditioner."
    "It uses the Lanczos algorithm and bisection for the tridiagonal matrix"
    );
  
  py::class_<QMRSolver<double>, shared_ptr<QMRSolver<double>>, BaseMatrix> (m, "QMRSolverD")
    ;
  py::class_<QMRSolver<Complex>, shared_ptr<QMRSolver<Complex>>, BaseMatrix> (m, "QMRSolverC")
    ;

  m.def("QMRSolver", [](shared_ptr<BaseMatrix> mat, shared_ptr<BaseMatrix> pre,
                        bool printrates, 
                        double precision, int maxsteps)
        {
          shared_ptr<KrylovSpaceSolver> solver;
          if (!mat->IsComplex())
            solver = make_shared<QMRSolver<double>> (mat, pre);
          else
            solver = make_shared<QMRSolver<Complex>> (mat, pre);                                            
          solver->SetPrecision(precision);
          solver->SetMaxSteps(maxsteps);
          solver->SetPrintRates (printrates);
          return solver;
        },
        py::arg("mat"), py::arg("pre"), py::arg("printrates")=true,
        py::arg("precision")=1e-8, py::arg("maxsteps")=200, docu_string(R"raw_string(
A Quasi Minimal Residuum (QMR) Solver.

Parameters:

mat : ngsolve.la.BaseMatrix
  input matrix 

pre : ngsolve.la.BaseMatrix
  input preconditioner matrix

printrates : bool
  input printrates

precision : float
  input requested precision. QMRSolver stops if precision is reached.

maxsteps : int
  input maximal steps. QMRSolver stops after this steps.

)raw_string"))
    ;
  
  m.def("ArnoldiSolver", [](shared_ptr<BaseMatrix> mata, shared_ptr<BaseMatrix> matm,
                            shared_ptr<BitArray> freedofs,
                            py::list vecs, Complex shift)
        {
          int nev;
          {
            py::gil_scoped_acquire acq;
            if (py::len(vecs) > mata->Height())
              throw Exception ("number of eigenvectors to compute "+ToString(py::len(vecs))
                               + " is greater than matrix dimension "
                               + ToString(mata->Height()));
            nev = py::len(vecs);
          }
          if (mata->IsComplex())
            {
              Arnoldi<Complex> arnoldi (mata, matm, freedofs);
              arnoldi.SetShift (shift);
              
              Array<shared_ptr<BaseVector>> evecs(nev);
                                                  
              Array<Complex> lam(nev);
              arnoldi.Calc (2*nev+1, lam, nev, evecs, 0);

              {
                py::gil_scoped_acquire acq;
                for (int i = 0; i < nev; i++)
                  vecs[i].cast<BaseVector&>() = *evecs[i];
              }
              
              Vector<Complex> vlam(nev);
              for (int i = 0; i < nev; i++)
                vlam(i) = lam[i];
              return vlam;
            }
          else
            {
              Arnoldi<double> arnoldi (mata, matm, freedofs);
              if (shift.imag())
                throw Exception("Only real shifts allowed for real arnoldi");
              arnoldi.SetShift (shift.real());
              
              Array<shared_ptr<BaseVector>> evecs(nev);
              
              Array<Complex> lam(nev);
              arnoldi.Calc (2*nev+1, lam, nev, evecs, 0);

              {
                py::gil_scoped_acquire acq;
                for (int i = 0; i < nev; i++)
                  vecs[i].cast<BaseVector&>() = *evecs[i];
              }
              
              Vector<Complex> vlam(nev);
              for (int i = 0; i < nev; i++)
                vlam(i) = lam[i];
              return vlam;
            }
        },
          py::arg("mata"), py::arg("matm"), py::arg("freedofs"), py::arg("vecs"), py::arg("shift")=DummyArgument(),
        py::call_guard<py::gil_scoped_release>(),
        docu_string(R"raw_string(
Shift-and-invert Arnoldi eigenvalue solver

Solves the generalized linear EVP A*u = M*lam*u using an Arnoldi iteration for the 
shifted EVP (A-shift*M)^(-1)*M*u = lam*u with a Krylow space of dimension 2*len(vecs)+1.
len(vecs) eigenpairs with the closest eigenvalues to the shift are returned.

Parameters:

mata : ngsolve.la.BaseMatrix
  matrix A

matm : ngsolve.la.BaseMatrix
  matrix M

freedofs : nsolve.ngstd.BitArray
  correct degrees of freedom

vecs : list
  list of BaseVectors for writing eigenvectors

shift : object
  complex or real shift
)raw_string"));
  
  

  m.def("DoArchive" , [](shared_ptr<Archive> & arch, BaseMatrix & mat)
                                         { cout << "output basematrix" << endl;
                                           mat.DoArchive(*arch); return arch; });
                                           
}


#endif // NGS_PYTHON
