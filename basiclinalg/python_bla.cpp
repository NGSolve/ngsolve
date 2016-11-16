#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <bla.hpp>

using namespace ngbla;

template <typename T> char * getPythonFormatString();
template <> char * getPythonFormatString<float>() { static char f[2]="f"; return f; }
template <> char * getPythonFormatString<double>() { static char f[2]="d"; return f; }
template <> char * getPythonFormatString<Complex>() { static char f[2]="D"; return f; }


template <typename T, int DIM>
struct PyBufferProtocol {
    typedef typename T::TSCAL TSCAL;

    template <class Tclass>
    void visit(Tclass& c) const {
//         c.def("NumPy",
//             [] (py::object & self)
//               {
//                 try {
//                   T& fv = py::cast<T&>(self);
//                   auto numpy = py::import("numpy");
//                   auto frombuffer = numpy.attr("frombuffer");
//                   if (DIM==1)
//                     return frombuffer(self);
//                   if (DIM==2)
//                     return frombuffer(self).attr("reshape")(fv.Height(),fv.Width());
//                 } catch (py::error_already_set const &) {
//                   PyErr_Print();
//                 }
//                 return py::object();
//               }  );

        // update exported type so it implements the buffer protocol
//         const py::converter::registration& fb_reg(py::converter::registry::lookup(py::type_id<T>())); TODO
//         PyTypeObject* fb_type = fb_reg.get_class_object();
// 
//         // buffer protocol
//         static PyBufferProcs buffer_functions = {
//             PyBufferProtocol<T, DIM>::getbuffer,/* bf_getbuffer */
//             PyBufferProtocol<T, DIM>::releasebuffer/* bf_releasebuffer */
//         };
// 
//         // register buffer protocol functions in PyType object
//         fb_type->tp_as_buffer = &buffer_functions;
    }

//     static int getbuffer(PyObject *exporter, Py_buffer *view, int flags) {
//         auto b = py::cast<T&>(exporter);
// 
// //         if (!b.check()) {
// //             PyErr_SetString(PyExc_BufferError, "Invalid exporter instance");
// //             view->obj = NULL;
// //             return -1;
// //         }
// 
//         T& array = b();
// 
//         if (view == NULL)
//           return 0;
// 
//         view->obj = exporter;
//         Py_INCREF(view->obj);
//         view->buf = &array(0);
//         view->len = sizeof(TSCAL) * array.Width() * array.Height();
//         view->readonly = PyBUF_WRITABLE;
//         view->itemsize = sizeof(TSCAL);
//         view->format = getPythonFormatString<TSCAL>();
//         view->ndim = DIM;
//         static_assert(DIM==1 || DIM==2, "invalid dimension for buffer protocol export!");
//         view->shape = NULL;
//         if(DIM==2)
//           {
//             view->shape = new Py_ssize_t[2];
//             view->shape[0] = array.Width();
//             view->shape[1] = array.Height();
//           }
//         view->strides = NULL;
//         view->suboffsets = NULL;
//         view->internal = NULL;
//         return 0;
//     }
// 
//     static void releasebuffer(PyObject *exporter, Py_buffer *view) {
//       auto b = py::cast<T&>(exporter);
// 
// //         if (!b.check()) {
// //             PyErr_SetString(PyExc_BufferError, "Invalid buffer exporter instance");
// //             return;
// //         }
// 
//         if (view->shape)
//           {
//             delete view->shape;
//             view->shape = NULL;
//           }
//     }
};

template <typename T, typename TNEW = T, typename TCLASS = py::class_<T> >
void PyVecAccess( py::module &m, TCLASS &c )
{
        typedef typename T::TSCAL TSCAL;
        c.def("__getitem__", [](T &self, py::slice inds )-> TNEW {
            size_t start, step, n;
            InitSlice( inds, self.Size(), start, step, n );
            TNEW res(n);
            for (int i=0; i<n; i++, start+=step)
                res[i] = self[start];
            return res;
            }  );
        c.def("__getitem__", [](T &v, py::list ind )-> TNEW {
                int n = py::len(ind);
                TNEW res(n);
                for (int i=0; i<n; i++) {
                    res[i] = v[ ind[i].cast<int>() ];
                }
                return res;
            }  );
        c.def("__setitem__", [](T &self, py::slice inds, const T & rv ) {
            size_t start, step, n;
            InitSlice( inds, self.Size(), start, step, n );
            for (int i=0; i<n; i++, start+=step)
                self[start] = rv[i];
            }  );
        c.def("__setitem__", [](T &self, py::slice inds, TSCAL val ) {
            size_t start, step, n;
            InitSlice( inds, self.Size(), start, step, n );
            for (int i=0; i<n; i++, start+=step)
                self[start] = val;
            }  );
        c.def("__add__" , [](T &self, T &v) { return TNEW(self+v); } );
        c.def("__sub__" , [](T &self, T &v) { return TNEW(self-v); } );
        c.def("__mul__" , [](T &self, TSCAL s) { return TNEW(s*self); } );
        c.def("__rmul__" , [](T &self, TSCAL s) { return TNEW(s*self); } );
        c.def("InnerProduct",  [](T & x, T & y) { return InnerProduct (x, y); });
        c.def("Norm",  [](T & x) { return L2Norm(x); });
}



template <typename TMAT, typename TNEW=TMAT, typename TCLASS = py::class_<TMAT> >
void PyMatAccess( TCLASS &c )
{
        // TODO: correct typedefs
        typedef typename TMAT::TSCAL TSCAL;
        typedef Vector<TSCAL> TROW;
        typedef Vector<TSCAL> TCOL;

        struct PyMatAccessHelper {
          static py::object GetTuple( TMAT & self, py::tuple t) {
            py::object rows = t[0];
            py::object cols = t[1];

            // First element of tuple is of type int
            int row_ind = rows.cast<int>();
            //         if(row_ind.check()) {
            py::object row = py::cast( self.Row(row_ind) );
            py::object f = row.attr("__getitem__");
            return f(cols);
            //         }

            // Second element of tuple is of type int
            auto col_ind = cols.cast<int>();
            //         if(col_ind.check()) {
            py::object col = py::cast( TROW(self.Col(col_ind)) );
            py::object f1 = col.attr("__getitem__");
            return f1(rows);
            //         }

            // Both elements are slices
            try {
              auto row_slice = rows.cast<py::slice>();
              auto col_slice = cols.cast<py::slice>();
              return py::cast(ColGetSlice(RowGetSlice(self, row_slice()), col_slice()));
            } catch (py::error_already_set const &) {
              cerr << "Invalid Matrix access!" << endl;
              PyErr_Print();
            }
            return py::object();
          }

          static void SetTupleVec( TMAT & self, py::tuple t, const FlatVector<TSCAL> &v) {
            py::object rows = t[0];
            py::object cols = t[1];

            // First element of tuple is of type int
            int row_ind = rows.cast<int>();
            //         if(row_ind.check()) {
            py::object row = py::cast( self.Row(row_ind) );
            py::object f = row.attr("__setitem__");
            f(cols, v);
            return;
            //         }

            // Second element of tuple is of type int
            auto col_ind = cols.cast<int>();
            //         if(col_ind.check()) {
            auto row_slice = rows.cast<py::slice>();
            auto col = self.Col(col_ind);
            size_t start, step, n;
            InitSlice( row_slice, self.Height(), start, step, n );
            for (int i=0; i<n; i++, start+=step)
              col[start] = v[i];
            return;
            //         }

            // One of the indices has to be of type int
            cerr << "Invalid Matrix access!" << endl;
          }

          static void SetTupleScal( TMAT & self, py::tuple t, TSCAL val) {
            py::object rows = t[0];
            py::object cols = t[1];

            // First element of tuple is of type int
            int row_ind = rows.cast<int>();
            //         if(row_ind.check()) {
            py::object row = py::cast( self.Row(row_ind) );
            py::object f = row.attr("__setitem__");
            f(cols, val);
            return;
            //         }

            // Second element of tuple is of type int
            int col_ind = cols.cast<int> ();
            //         if(col_ind.check()) {
            py::slice row_slice = rows.cast<py::slice> ();
            auto col = self.Col(col_ind);
            size_t start, step, n;
            InitSlice( row_slice, self.Height(), start, step, n );
            for (int i=0; i<n; i++, start+=step)
              col[start] = val;
            return;
            //         }

            // Both elements are slices
            try {
              py::slice row_slice = rows.cast<py::slice> ();
              size_t start, step, n;
              InitSlice( row_slice(), self.Height(), start, step, n );
              for (int i=0; i<n; i++, start+=step) {
                py::object row = py::cast(self.Row(start));
                py::object f = row.attr("__setitem__");
                f(cols, val);
                return;
              }
            } catch (py::error_already_set const &) {
              cerr << "Invalid Matrix access!" << endl;
              PyErr_Print();
            }
          }

          static void SetTuple( TMAT & self, py::tuple t, const TMAT & rmat) {
            py::object rows = t[0];
            py::object cols = t[1];

            // Both elements have to be slices
            try {
              auto row_slice = rows.cast<py::slice> ();
              auto col_slice = cols.cast<py::slice> ();
              size_t start, step, n;
              InitSlice( row_slice(), self.Height(), start, step, n );
              for (int i=0; i<n; i++, start+=step) {
                py::object row = py::cast(self.Row(start));
                py::object f = row.attr("__setitem__");
                f(self, cols, rmat.Row(i));
              }
            } catch (py::error_already_set const &) {
              cerr << "Invalid Matrix access!" << endl;
              PyErr_Print();
            }
          }

          static TROW RowGetInt( TMAT & self,  int ind )  {
            return self.Row(ind);
          }

          static TNEW RowGetSlice( TMAT & self,  py::slice inds ) {
            size_t start, step, n;
            InitSlice( inds, self.Height(), start, step, n );
            TNEW res(n, self.Width());
            for (int i=0; i<n; i++, start+=step)
              res.Row(i) = self.Row(start);
            return res;
          }

          static void RowSetInt( TMAT & self,  int ind, const TROW &r ) {
            self.Row(ind) = r;
          }

          static void RowSetIntScal( TMAT & self,  int ind, TSCAL r ) {
            self.Row(ind) = r;
          }

          static void RowSetSlice( TMAT & self,  py::slice inds, const TMAT &r ) {
            size_t start, step, n;
            InitSlice( inds, self.Height(), start, step, n );
            for (int i=0; i<n; i++, start+=step)
              self.Row(start) = r.Row(i);
          }

          static void RowSetSliceScal( TMAT & self,  py::slice inds, TSCAL r ) {
            size_t start, step, n;
            InitSlice( inds, self.Height(), start, step, n );
            for (int i=0; i<n; i++, start+=step)
              self.Row(start) = r;
          }

          static Vector<TSCAL> ColGetInt( TMAT & self,  int ind )  {
            return Vector<TSCAL>(self.Col(ind));
          }

          static TNEW ColGetSlice( const TMAT & self,  py::slice inds ) {
            size_t start, step, n;
            InitSlice( inds, self.Width(), start, step, n );
            TNEW res(self.Height(),n);
            for (int i=0; i<n; i++, start+=step)
              res.Col(i) = self.Col(start);
            return res;
          }

          static void ColSetInt( TMAT & self,  int ind, const TCOL &r ) {
            self.Col(ind) = r;
          }

          static void ColSetIntScal( TMAT & self,  int ind, TSCAL r ) {
            self.Col(ind) = r;
          }

          static void ColSetSlice( TMAT & self,  py::slice inds, const TMAT &r ) {
            size_t start, step, n;
            InitSlice( inds, self.Width(), start, step, n );
            for (int i=0; i<n; i++, start+=step)
              self.Col(start) = r.Col(i);
          }

          static void ColSetSliceScal( TMAT & self,  py::slice inds, TSCAL r ) {
            size_t start, step, n;
            InitSlice( inds, self.Width(), start, step, n );
            for (int i=0; i<n; i++, start+=step)
              self.Col(start) = r;
          }
        };
        c.def("__getitem__", &PyMatAccessHelper::GetTuple);
        c.def("__getitem__", &PyMatAccessHelper::RowGetInt);
        c.def("__getitem__", &PyMatAccessHelper::RowGetSlice);

        c.def("__setitem__", &PyMatAccessHelper::SetTuple);
        c.def("__setitem__", &PyMatAccessHelper::SetTupleScal);
        c.def("__setitem__", &PyMatAccessHelper::SetTupleVec);
        c.def("__setitem__", &PyMatAccessHelper::RowSetInt);
        c.def("__setitem__", &PyMatAccessHelper::RowSetIntScal);
        c.def("__setitem__", &PyMatAccessHelper::RowSetSlice);
        c.def("__setitem__", &PyMatAccessHelper::RowSetSliceScal);

        c.def_property("diag",
                py::cpp_function([](TMAT &self) { return Vector<TSCAL>(self.Diag()); }),
                py::cpp_function([](TMAT &self, const FlatVector<TSCAL> &v) { self.Diag() = v; }));
        c.def("__add__" , [](TMAT &self, TMAT &m) { return TNEW(self+m); } );
        c.def("__sub__" , [](TMAT &self, TMAT &m) { return TNEW(self-m); } );
        c.def("__mul__" , [](TMAT &self, TMAT &m) { return TNEW(self*m); } );
        c.def("__mul__" , [](TMAT &self, FlatVector<TSCAL> &v) { return Vector<TSCAL>(self*v); } );
        c.def("__mul__" , [](TMAT &self, TSCAL s) { return TNEW(s*self); } );
        c.def("__rmul__" , [](TMAT &self, TSCAL s) { return TNEW(s*self); } );
        c.def("Height", &TMAT::Height );
        c.def("Width", &TMAT::Width );
        c.def_property_readonly("h", py::cpp_function(&TMAT::Height ));
        c.def_property_readonly("w", py::cpp_function(&TMAT::Width ));
        c.def_property_readonly("T", py::cpp_function([](TMAT &self) { return TNEW(Trans(self)); } ) );
        c.def_property_readonly("A", py::cpp_function([](TMAT &self) { return Vector<TSCAL>(FlatVector<TSCAL>( self.Width()* self.Height(), &self(0,0)) ); } ) );
        c.def("__len__", []( TMAT& self) { return self.Height();}  );
}


template <typename TVEC, typename TNEW, typename TSCAL>
auto ExportVector(py::module &m, const char * name ) -> py::class_<TVEC>
  {
    auto c = py::class_<TVEC >(m, name);
    PyDefVector<TVEC, TSCAL>(m, c);
//     PyBufferProtocol<TVEC, 1>(m, c);
    PyVecAccess< TVEC, TNEW >(m, c);
    PyDefToString<TVEC >(m, c);
    c.def(py::self+=py::self);
    c.def(py::self-=py::self);
    c.def(py::self*=TSCAL());
    return c;
  }

void NGS_DLL_HEADER ExportNgbla(py::module & m) {

    ///////////////////////////////////////////////////////////////////////////////////////
    // Vector types
    typedef FlatVector<double> FVD;
    typedef FlatVector<Complex> FVC;
    typedef SliceVector<double> SVD;
    typedef SliceVector<Complex> SVC;
    typedef Vector<double> VD;
    typedef Vector<Complex> VC;

    ExportVector< FVD, VD, double>(m, "FlatVectorD")
        .def(py::init<int, double *>())
        .def("Range",    static_cast</* const */ FVD (FVD::*)(size_t,size_t) const> (&FVD::Range ) )
      ;
    ExportVector< FVC, VC, Complex>(m, "FlatVectorC")
        .def(py::self*=double())
        .def(py::init<int, Complex *>())
        .def("Range",    static_cast</* const */ FVC (FVC::*)(size_t,size_t) const> (&FVC::Range ) )
        ;

    ExportVector< SVD, VD, double>(m, "SliceVectorD")
        .def("Range",    static_cast<const SVD (SVD::*)(size_t,size_t) const> (&SVD::Range ) )
        ;
    ExportVector< SVC, VC, Complex>(m, "SliceVectorC")
        .def("Range",    static_cast<const SVC (SVC::*)(size_t,size_t) const> (&SVC::Range ) )
        .def(py::self*=double())
        ;

    py::class_<VD, FVD>(m, "VectorD")
        .def(py::init<int>())
        ;

    py::class_<VC, FVC >(m, "VectorC")
        .def(py::init<int>())
        ;

    m.def("Vector",
            [] (int n, bool is_complex) {
                if(is_complex) return py::cast(Vector<Complex>(n));
                else return py::cast(Vector<double>(n));
                },
            py::arg("length"),
            py::arg("complex")=false
           );
    m.def("Vector",
            [] (py::list values) {
                               Vector<> v(len(values));
                               for (int i = 0; i < v.Size(); i++)
                                 v(i) = values[i].cast<double>();
                               return v;
                             },
            py::arg("vals")
           );
    m.def("Vector",
            [] (py::tuple values) ->py::object {
                               bool is_double = true;
//                                for (int i = 0; i < len(values); i++)
//                                  is_double &= values[i].cast<double>();  // TODO
                               if (is_double)
                                 {
                                   Vector<> v(len(values));
                                   for (int i = 0; i < v.Size(); i++)
                                     v(i) = values[i].cast<double>();
                                   return py::cast(v);
                                 }
                               bool is_complex = true;
//                                for (int i = 0; i < len(values); i++)
//                                  is_complex &= values[i].cast<Complex>(); // TODO
                               if (is_complex)
                                 {
                                   Vector<Complex> v(len(values));
                                   for (int i = 0; i < v.Size(); i++)
                                     v(i) = values[i].cast<Complex>();
                                   return py::cast(v);
                                 }
                               throw Exception("cannot make a vector from tuple");
              },
            py::arg("vals")
           );


    ///////////////////////////////////////////////////////////////////////////////////////
    // Matrix types
    typedef FlatMatrix<double> FMD;
    py::class_<FlatMatrix<double> > class_FMD(m, "FlatMatrixD");
        PyMatAccess<FMD, Matrix<double> >(class_FMD);
        PyDefToString<FMD>(m, class_FMD);
        class_FMD.def(py::self+=py::self);
        class_FMD.def(py::self-=py::self);
        class_FMD.def(py::self*=double());
//         .def(PyBufferProtocol<FMD, 2>())
        class_FMD.def("Inverse", [](FMD & self, FMD & inv) {
	    CalcInverse(self,inv); return;
	  });


    typedef FlatMatrix<Complex> FMC;
    py::class_<FlatMatrix<Complex> >(m, "FlatMatrixC")
//         .def(PyDefToString<FMC>())
//         .def(PyMatAccess<FMC, Matrix<Complex> >())
        .def(py::self+=py::self)
        .def(py::self-=py::self)
        .def(py::self*=Complex())
//         .def(PyBufferProtocol<FMC, 2>())
        .def_property("diag",
                py::cpp_function([](const FMC &self) { return Vector<Complex>(self.Diag()); }),
                py::cpp_function([](FMC &self, const FVC &v) { self.Diag() = v; }))
        .def("__add__" , [](FMC &self, FMD &m) { return Matrix<Complex>(self+m); } )
        .def("__sub__" , [](FMC &self, FMD &m) { return Matrix<Complex>(self-m); } )
        .def("__mul__" , [](FMC &self, FMD &m) { return Matrix<Complex>(self*m); } )
        .def("__radd__" , [](FMC &self, FMD &m) { return Matrix<Complex>(self+m); } )
        .def("__rsub__" , [](FMC &self, FMD &m) { return Matrix<Complex>(self-m); } )
        .def("__rmul__" , [](FMC &self, FMD &m) { return Matrix<Complex>(m*self); } )
        .def("__mul__" , [](FMC &self, FVD &v) { return Vector<Complex>(self*v); } )
        .def("__mul__" , [](FMC &self, double s) { return Matrix<Complex>(s*self); } )
        .def("__rmul__" , [](FMC &self, double s) { return Matrix<Complex>(s*self); } )
        .def("Height", &FMC::Height )
        .def("Width", &FMC::Width )
        .def("__len__", []( FMC& self) { return self.Height();}  )
        .def_property_readonly("h", py::cpp_function(&FMC::Height ))
        .def_property_readonly("w", py::cpp_function(&FMC::Width ))
        .def_property_readonly("A", py::cpp_function([](FMC &self) { return Vector<Complex>(FlatVector<Complex>( self.Width()* self.Height(), &self(0,0) )); }  ))
        .def_property_readonly("T", py::cpp_function([](FMC &self) { return Matrix<Complex>(Trans(self)); } ) )
        .def_property_readonly("C", py::cpp_function([](FMC &self) { 
            Matrix<Complex> result( self.Height(), self.Width() );
            for (int i=0; i<self.Height(); i++)
                for (int j=0; j<self.Width(); j++) 
                    result(i,j) = Conj(self(i,j));
            return result;
            } ) )
        .def_property_readonly("H", py::cpp_function([](FMC &self) { 
            Matrix<Complex> result( self.Width(), self.Height() );
            for (int i=0; i<self.Height(); i++)
                for (int j=0; j<self.Width(); j++) 
                    result(j,i) = Conj(self(i,j));
            return result;
            } ) )
        ;

    py::class_<Matrix<double>, FMD>(m, "MatrixD")
        .def(py::init<int, int>())
        ;

    py::class_<Matrix<Complex>, FMC >(m, "MatrixC")
        .def(py::init<int, int>())
        ;

    m.def("Matrix",
            [] (int h, int w, bool is_complex) {
                if(is_complex) return py::cast(Matrix<Complex>(h,w));
                else return py::cast(Matrix<double>(h,w));
                },
            py::arg("height"), 
            py::arg("width"), 
            py::arg("complex")=false
           );

    m.def("InnerProduct",
             [] (py::object x, py::object y) -> py::object
                              { return py::object(x.attr("InnerProduct")) (y); });

    m.def("Norm",
             [] (py::object x) -> py::object
                              { return py::object(x.attr("Norm")) (); });


    m.def("CheckPerformance",
             [] (int n, int m, int k)
                              {
                                Matrix<> a(n,k), b(m,k), c(n,m);
                                a = 1; b = 2;
                                double tot = double(n)*k*m;
                                c = a * Trans(b) | Lapack; // warmup
                                int its = 1e10 / tot + 1;
                                {
                                  Timer t("matmat");
                                  t.Start();
                                  for (int j = 0; j < its; j++)
                                    c = a * Trans(b) | Lapack;
                                  t.Stop();
                                  cout << "Lapack GFlops = " << 1e-9 * n*k*m*its / t.GetTime() << endl;
                                }

                                /*
                                {
                                  Timer t("matmat2");
                                  t.Start();
                                  c = a * b;
                                  t.Stop();
                                  cout << "without Lapack GFlops = " << 1e-9 * n*n*n / t.GetTime() << endl;
                                }
                                */

                                {
                                Timer t2("matmat - par");
                                c = a * Trans(b) | Lapack; // warmup
                                t2.Start();
                                RunWithTaskManager
                                  ([=] ()
                                   {
                                     ParallelFor (Range(8), [=] (int nr)
                                                  {
                                                    Matrix<> a(n,k), b(m,k), c(n,m);
                                                    a = 1; b = 2;
                                                    for (int j = 0; j < its; j++)
                                                      c = a * Trans(b) | Lapack;
                                                  });
                                   }
                                   );
                                t2.Stop();
                                cout << "Task-manager Lapack GFlops = " << 8 * 1e-9 * n*k*m*its / t2.GetTime() << endl;
                                }
                              });
             }

PYBIND11_PLUGIN(libngbla) {
  py::module m("bla", "pybind bla");
  ExportNgbla(m);
  return m.ptr();
}

#endif // NGS_PYTHON
