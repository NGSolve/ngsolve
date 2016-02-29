#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <boost/python/slice.hpp>
#include <la.hpp>
using namespace ngla;


static void InitSlice( const bp::slice &inds, int len, int &start, int &step, int &n ) {
    bp::object indices = inds.attr("indices")(len);
    start = 0;
    step = 1;
    n = 0;

    try {
        start = bp::extract<int>(indices[0]);
        int stop  = bp::extract<int>(indices[1]);
        step  = bp::extract<int>(indices[2]);
        n = (stop-start+step-1) / step;
    }
    catch (bp::error_already_set const &) {
        cout << "Error in InitSlice(slice,...): " << endl;
        PyErr_Print();
    }

}




void NGS_DLL_HEADER ExportNgla() {
    std::string nested_name = "la";
    if( bp::scope() )
         nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".la");
    
    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << IM(1) << "exporting la as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("la") = module ;

    bp::scope ngla_scope(module);

    bp::object expr_module = bp::import("ngsolve.__expr");
    bp::object expr_namespace = expr_module.attr("__dict__");




    struct BaseVector_pickle_suite : bp::pickle_suite
    {
      static
      bp::tuple getinitargs(const BaseVector & v)
      {
        return bp::make_tuple(v.Size(), v.IsComplex()); 
      }

      static
      bp::tuple getstate(bp::object obj)
      {
        /*
        const BaseVector & vec = bp::extract<BaseVector const&>(obj)();
        bp::list data;
        for (int i : Range(vec))
          data.append (vec.FV<double>()(i));
        */
        auto vec = bp::extract<BaseVector const&>(obj)().FV<double>();
        bp::list data;
        for (int i : Range(vec))
          data.append (vec(i));
        return bp::make_tuple (obj.attr("__dict__"), data);
      }
    
      static
      void setstate(bp::object obj, bp::tuple state)
      {
        bp::dict d = bp::extract<bp::dict>(obj.attr("__dict__"))();
        d.update(state[0]);
        bp::list data = bp::extract<bp::list>(state[1]);
        auto vec = bp::extract<BaseVector const&>(obj)().FV<double>();
        for (int i : Range(vec))
          vec(i) = bp::extract<double> (data[i]);
        /*
        BaseVector & vec = bp::extract<BaseVector&>(obj)();
        for (int i : Range(vec))
          vec.FV<double>()(i) = bp::extract<double> (data[i]);
        */
      }

    static bool getstate_manages_dict() { return true; }
  };
    
    
  
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<BaseVector>);
  bp::class_<BaseVector, shared_ptr<BaseVector>, boost::noncopyable>("BaseVector", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](int size, bool is_complex) -> shared_ptr<BaseVector>
                           {
                             if (is_complex)
                               return make_shared<VVector<Complex>> (size);
                             else
                               return make_shared<VVector<double>> (size);                               
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("size"), bp::arg("complex")=false)
          ))
    .def_pickle(BaseVector_pickle_suite())
    .def("__ngsid__", FunctionPointer( [] ( BaseVector & self)
        { return reinterpret_cast<std::uintptr_t>(&self); } ) )
    
    .def("__str__", &ToString<BaseVector>)
    .add_property("size", &BaseVector::Size)
    .def("__len__", &BaseVector::Size)
    .def("CreateVector", FunctionPointer( [] ( BaseVector & self)
        { return shared_ptr<BaseVector>(self.CreateVector()); } ))

    /*
    .def("Assign", FunctionPointer([](BaseVector & self, BaseVector & v2, double s)->void { self.Set(s, v2); }))
    .def("Add", FunctionPointer([](BaseVector & self, BaseVector & v2, double s)->void { self.Add(s, v2); }))
    .def("Assign", FunctionPointer([](BaseVector & self, BaseVector & v2, Complex s)->void { self.Set(s, v2); }))
    .def("Add", FunctionPointer([](BaseVector & self, BaseVector & v2, Complex s)->void { self.Add(s, v2); }))
    */
    .def("Assign", FunctionPointer([](BaseVector & self, BaseVector & v2, bp::object s)->void 
                                   { 
                                     if ( bp::extract<double>(s).check() )
                                       {
                                         self.Set (bp::extract<double>(s)(), v2);
                                         return;
                                       }
                                     if ( bp::extract<Complex>(s).check() )
                                       {
                                         self.Set (bp::extract<Complex>(s)(), v2);
                                         return;
                                       }
                                     throw Exception ("BaseVector::Assign called with non-scalar type");
                                   }))
    .def("Add", FunctionPointer([](BaseVector & self, BaseVector & v2, bp::object s)->void 
                                   { 
                                     if ( bp::extract<double>(s).check() )
                                       {
                                         self.Add (bp::extract<double>(s)(), v2);
                                         return;
                                       }
                                     if ( bp::extract<Complex>(s).check() )
                                       {
                                         self.Add (bp::extract<Complex>(s)(), v2);
                                         return;
                                       }
                                     throw Exception ("BaseVector::Assign called with non-scalar type");
                                   }))


    .add_property("expr", bp::object(expr_namespace["VecExpr"]) )
    .add_property("data", bp::object(expr_namespace["VecExpr"]), bp::object(expr_namespace["expr_data"] ))
    .def("__add__" , bp::object(expr_namespace["expr_add"]) )
    .def("__sub__" , bp::object(expr_namespace["expr_sub"]) )
    .def("__rmul__" , bp::object(expr_namespace["expr_rmul"]) )
    .def("__getitem__", FunctionPointer
         ( [](BaseVector & self,  int ind )
           {
             if (ind < 0 || ind >= self.Size()) 
               bp::exec("raise IndexError()\n");
             if( self.IsComplex() )
               return bp::object(self.FVComplex()[ind]);
             else
               return bp::object(self.FVDouble()[ind]);
           } ))
    .def("__getitem__", FunctionPointer( [](BaseVector & self,  bp::slice inds )
      {
          int start, step, n;
          InitSlice( inds, self.Size(), start, step, n );

          if( self.IsComplex() )
            {
              if( step == 1 )
                return bp::object(self.FVComplex().Range(start, start+n));
              else
                return bp::object(self.FVComplex().Slice(start, step).Range(n));
            }
          else
            {
              if( step == 1 )
                return bp::object(self.FVDouble().Range(start, start+n));
              else
                return bp::object(self.FVDouble().Slice(start, step).Range(n));
            }
      } ))
    .def("__setitem__", FunctionPointer( [](BaseVector & self,  int ind, Complex z )
      {
          self.FVComplex()[ind] = z;
      } ))
    .def("__setitem__", FunctionPointer( [](BaseVector & self,  int ind, double d )
      {
          self.FVDouble()[ind] = d;
      } ))
    .def("__setitem__", FunctionPointer( [](BaseVector & self,  bp::slice inds, Complex z )
      {
          int start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if( step == 1 )
            self.FVComplex().Range(start,start+n) = z;
          else
            self.FVComplex().Slice(start,step).Range(n) = z;
      } ))
    .def("__setitem__", FunctionPointer( [](BaseVector & self,  bp::slice inds, double d )
      {
          int start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if( step == 1 )
            self.FVDouble().Range(start,start+n) = d;
          else
            self.FVDouble().Slice(start,step).Range(n) = d;
      } ))
    .def("__setitem__", FunctionPointer( [](BaseVector & self,  bp::slice inds, FlatVector<Complex> & v )
      {
          int start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if( step == 1 )
            self.FVComplex().Range(start,start+n) = v;
          else
            self.FVComplex().Slice(start,step).Range(n) = v;
      } ))
    .def("__setitem__", FunctionPointer( [](BaseVector & self,  bp::slice inds, FlatVector<double> & v )
      {
          int start, step, n;
          InitSlice( inds, self.Size(), start, step, n );
          if( step == 1 )
            self.FVDouble().Range(start,start+n) = v;
          else
            self.FVDouble().Slice(start,step).Range(n) = v;
      } ))
    .def(bp::self+=bp::self)
    .def(bp::self-=bp::self)
    .def(bp::self*=double())
    .def("InnerProduct", FunctionPointer( [](BaseVector & self, BaseVector & other)
                                          {
                                            if (self.IsComplex())
                                              return bp::object (S_InnerProduct<ComplexConjugate> (self, other));
                                            else
                                              return bp::object (InnerProduct (self, other));
                                          }))
    .def("Norm", FunctionPointer ( [](BaseVector & self) { return self.L2Norm(); }));
  ;       

  // bp::def("InnerProduct", FunctionPointer([](BaseVector & v1, BaseVector & v2)->double { return InnerProduct(v1,v2); }))
  bp::def ("InnerProduct",
           FunctionPointer( [] (bp::object x, bp::object y) -> bp::object
                            { return x.attr("InnerProduct") (y); }));
  ;
  

  typedef BaseMatrix BM;
  typedef BaseVector BV;

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<BaseMatrix>);
  bp::class_<BaseMatrix, shared_ptr<BaseMatrix>, boost::noncopyable>("BaseMatrix", bp::no_init)
    .def("__str__", &ToString<BaseMatrix>)
    .add_property("height", &BaseMatrix::Height)
    .add_property("width", &BaseMatrix::Width)

    .def("CreateMatrix", &BaseMatrix::CreateMatrix)

    .def("CreateRowVector", FunctionPointer( [] ( BaseMatrix & self)
        { return shared_ptr<BaseVector>(self.CreateRowVector()); } ))
    .def("CreateColVector", FunctionPointer( [] ( BaseMatrix & self)
        { return shared_ptr<BaseVector>(self.CreateColVector()); } ))

    .def("AsVector", FunctionPointer( [] (BM & m)
                                      {
                                        return shared_ptr<BaseVector> (&m.AsVector(), NOOP_Deleter);
                                      }))
    .def("COO", FunctionPointer( [] (BM & m) -> bp::object
                                 {
                                   SparseMatrix<double> * sp = dynamic_cast<SparseMatrix<double>*> (&m);
                                   if (sp)
                                     {
                                       Array<int> ri, ci;
                                       Array<double> vals;
                                       for (int i = 0; i < sp->Height(); i++)
                                         {
                                           FlatArray<int> ind = sp->GetRowIndices(i);
                                           FlatVector<double> rv = sp->GetRowValues(i);
                                           for (int j = 0; j < ind.Size(); j++)
                                             {
                                               ri.Append (i);
                                               ci.Append (ind[j]);
                                               vals.Append (rv[j]);
                                             }
                                         }

                                       bp::list pyri (ri);
                                       bp::list pyci (ci);
                                       bp::list pyvals (vals);
                                       return bp::make_tuple (pyri, pyci, pyvals);
                                     }
				   throw Exception ("COO needs real-valued sparse matrix");
                                 }))

    .def("Mult",        FunctionPointer( [](BM &m, BV &x, BV &y, double s) { m.Mult (x,y); y *= s; }) )
    .def("MultAdd",     FunctionPointer( [](BM &m, BV &x, BV &y, double s) { m.MultAdd (s, x, y); }))
    // .def("MultTrans",   FunctionPointer( [](BM &m, BV &x, BV &y, double s) { y  = s*Trans(m)*x; }) )
    // .def("MultTransAdd",FunctionPointer( [](BM &m, BV &x, BV &y, double s) { y += s*Trans(m)*x; }) )

    .add_property("expr", bp::object(expr_namespace["MatExpr"]) )
    .def("__mul__" , bp::object(expr_namespace["expr_mul"]) )
    .def("__rmul__" , bp::object(expr_namespace["expr_rmul"]) )

    .def("__iadd__", FunctionPointer( [] (BM &m, BM &m2) { 
        m.AsVector()+=m2.AsVector();
    }))

    .def("Inverse", FunctionPointer( [](BM &m, BitArray * freedofs, string inverse)
                                     ->shared_ptr<BaseMatrix>
                                     { 
                                       if (inverse != "") m.SetInverseType(inverse);
                                       return m.InverseMatrix(freedofs);
                                     }),
         (bp::arg("self"), bp::arg("freedofs"), bp::arg("inverse")=""))
    .def("Inverse", FunctionPointer( [](BM &m)->shared_ptr<BaseMatrix>
                                     { return m.InverseMatrix(); }))
    .def("Transpose", FunctionPointer( [](BM &m)->shared_ptr<BaseMatrix>
                                       { return make_shared<Transpose> (m); }))
    // bp::return_value_policy<bp::manage_new_object>())
    ;



  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<CGSolver<double>>);
  bp::class_<CGSolver<double>, shared_ptr<CGSolver<double>>,bp::bases<BaseMatrix>,boost::noncopyable> ("CGSolverD", bp::no_init)
    .def("GetSteps", &CGSolver<double>::GetSteps)
    ;

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<CGSolver<Complex>>);
  bp::class_<CGSolver<Complex>, shared_ptr<CGSolver<Complex>>,bp::bases<BaseMatrix>,boost::noncopyable> ("CGSolverC", bp::no_init)
    .def("GetSteps", &CGSolver<Complex>::GetSteps)
    ;

  bp::def("CGSolver", FunctionPointer ([](const BaseMatrix & mat, const BaseMatrix & pre,
                                          bool iscomplex, bool printrates, 
                                          double precision, int maxsteps) -> BaseMatrix *
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
                                         return solver;
                                       }),
          (bp::arg("mat"), bp::arg("pre"), bp::arg("complex") = false, bp::arg("printrates")=true,
           bp::arg("precision")=1e-8, bp::arg("maxsteps")=200),
          bp::return_value_policy<bp::manage_new_object>()
          )
    ;

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<QMRSolver<double>>);
  bp::class_<QMRSolver<double>, shared_ptr<QMRSolver<double>>,bp::bases<BaseMatrix>,boost::noncopyable> ("QMRSolverD", bp::no_init)
    ;
  bp::def("QMRSolver", FunctionPointer ([](const BaseMatrix & mat, const BaseMatrix & pre,
                                           bool printrates, 
                                           double precision, int maxsteps) -> BaseMatrix *
                                        {
                                          KrylovSpaceSolver * solver;
                                          solver = new QMRSolver<double> (mat, pre);
                                          solver->SetPrecision(precision);
                                          solver->SetMaxSteps(maxsteps);
                                          solver->SetPrintRates (printrates);
                                          return solver;
                                        }),
          (bp::arg("mat"), bp::arg("pre"), bp::arg("printrates")=true,
           bp::arg("precision")=1e-8, bp::arg("maxsteps")=200),
          bp::return_value_policy<bp::manage_new_object>()
          )
    ;
  
  bp::def("ArnoldiSolver", FunctionPointer ([](BaseMatrix & mata, BaseMatrix & matm, const BitArray & freedofs,
                                               bp::list vecs, bp::object bpshift)
                                            {
                                              if (mata.IsComplex())
                                                {
                                                  Arnoldi<Complex> arnoldi (mata, matm, &freedofs);
                                                  Complex shift = 0.0;
                                                  if (bp::extract<Complex>(bpshift).check())
                                                    shift = bp::extract<Complex>(bpshift)();
                                                  cout << "shift = " << shift << endl;
                                                  arnoldi.SetShift (shift);
                                                  
                                                  int nev = bp::len(vecs);
                                                  cout << "num vecs: " << nev << endl;
                                                  Array<shared_ptr<BaseVector>> evecs(nev);
                                                  
                                                  Array<Complex> lam(nev);
                                                  arnoldi.Calc (2*nev+1, lam, nev, evecs, 0);
                                            
                                                  for (int i = 0; i < nev; i++)
                                                    * bp::extract<shared_ptr<BaseVector>>(vecs[i])() = *evecs[i];

                                                  Vector<Complex> vlam(nev);
                                                  for (int i = 0; i < nev; i++)
                                                    vlam(i) = lam[i];
                                                  return vlam;
                                                }
                                              
                                              cout << "real Arnoldi not supported" << endl;
                                              Vector<Complex> lam(5);
                                              return lam;
                                            }),
          (bp::arg("mata"), bp::arg("matm"), bp::arg("freedofs"), bp::arg("vecs"), bp::arg("shift")=bp::object())
          // bp::return_value_policy<bp::manage_new_object>()
          )
    ;

  

  bp::def("DoArchive" , FunctionPointer( [](shared_ptr<Archive> & arch, BaseMatrix & mat) 
                                         { cout << "output basematrix" << endl;
                                           mat.DoArchive(*arch); return arch; }));
                                           
}



void ExportNgbla();

BOOST_PYTHON_MODULE(libngla) {
  // ExportNgbla();
  ExportNgla();
}






#endif // NGS_PYTHON
