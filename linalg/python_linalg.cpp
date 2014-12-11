#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <boost/python/slice.hpp>
#include <la.hpp>
using namespace ngla;





void ExportNgla() {
    std::string nested_name = "la";
    if( bp::scope() )
         nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".la");
    
    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << "exporting la as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("la") = module ;

    bp::scope ngla_scope(module);

    bp::object expr_module = bp::import("ngsolve.__expr");
    bp::object expr_namespace = expr_module.attr("__dict__");


  
  bp::class_<BaseVector, shared_ptr<BaseVector>, boost::noncopyable>("BaseVector", bp::no_init)
    .def("__str__", &ToString<BaseVector>)
    .add_property("size", &BaseVector::Size)
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
    .def("__getitem__", FunctionPointer( [](BaseVector & self,  bp::slice inds ) {
        return self.FVDouble();
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
    ;       

  // bp::def("InnerProduct", FunctionPointer([](BaseVector & v1, BaseVector & v2)->double { return InnerProduct(v1,v2); }))
  bp::def ("InnerProduct",
           FunctionPointer( [] (bp::object x, bp::object y) -> bp::object
                            { return x.attr("InnerProduct") (y); }));
  ;
  

  typedef BaseMatrix BM;
  typedef BaseVector BV;

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
				   throw Exception ("COO needs sparse matrix");
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

    .def("Inverse", FunctionPointer( [](BM &m, BitArray & freedofs)->shared_ptr<BaseMatrix>
                                     { return m.InverseMatrix(&freedofs); }))
    // bp::return_value_policy<bp::manage_new_object>(),
    // (bp::arg("self"), bp::arg("freedofs")))
    .def("Inverse", FunctionPointer( [](BM &m)->shared_ptr<BaseMatrix>
                                     { return m.InverseMatrix(); }))
    // bp::return_value_policy<bp::manage_new_object>())
    ;



  bp::class_<CGSolver<double>, shared_ptr<CGSolver<double>>,bp::bases<BaseMatrix>,boost::noncopyable> ("CGSolverD", bp::no_init)
    // .def(bp::init<const BaseMatrix &, const BaseMatrix &>())
    ;
  bp::class_<CGSolver<Complex>, shared_ptr<CGSolver<Complex>>,bp::bases<BaseMatrix>,boost::noncopyable> ("CGSolverC", bp::no_init)
    // .def(bp::init<const BaseMatrix &, const BaseMatrix &>())
    ;

  bp::def("CGSolver", FunctionPointer ([](const BaseMatrix & mat, const BaseMatrix & pre,
                                          bool iscomplex, bool printrates) -> BaseMatrix *
                                       {
                                         KrylovSpaceSolver * solver;
                                         if (iscomplex)
                                           solver = new CGSolver<Complex> (mat, pre);
                                         else
                                           solver = new CGSolver<double> (mat, pre);                               
                                         solver->SetPrintRates (printrates);
                                         return solver;
                                       }),
          (bp::arg("mat"), bp::arg("pre"), bp::arg("complex") = false, bp::arg("printrates")=true),
          bp::return_value_policy<bp::manage_new_object>()
          )
    ;
}



void ExportNgbla();

BOOST_PYTHON_MODULE(libngla) {
  // ExportNgbla();
  ExportNgla();
}






#endif // NGS_PYTHON
