#ifdef NGS_PYTHON
#include <pybind11/numpy.h>
#include "../ngstd/python_ngstd.hpp"
#include <bla.hpp>

using namespace ngbla;

template<typename T, typename TCLASS>
void PyDefVecBuffer( TCLASS & c )
{
    typedef typename T::TSCAL TSCAL;
    c.def_buffer([](T &self) -> py::buffer_info {
        return py::buffer_info(
            &self(0),                                     /* Pointer to buffer */
            sizeof(TSCAL),                                /* Size of one scalar */
            py::format_descriptor<TSCAL>::format(),       /* Python struct-style format descriptor */
            1,                                            /* Number of dimensions */
            { self.Size() },                              /* Buffer dimensions */
            { sizeof(TSCAL)*(self.Addr(1)-self.Addr(0)) } /* Strides (in bytes) for each index */
        );
    });
    c.def("NumPy", [] (py::object & self) {
        // T& fv = py::cast<T&>(self);
        auto numpy = py::module::import("numpy");
        auto frombuffer = numpy.attr("frombuffer");
        return frombuffer(self, py::detail::npy_format_descriptor<TSCAL>::dtype());
      }, "Return NumPy object");
}

template<typename T, typename TCLASS>
void PyDefMatBuffer( TCLASS & c )
{
    typedef typename T::TSCAL TSCAL;
    c.def_buffer([](T &self) -> py::buffer_info {
        return py::buffer_info(
            &self(0),                                     /* Pointer to buffer */
            sizeof(TSCAL),                                /* Size of one scalar */
            py::format_descriptor<TSCAL>::format(),       /* Python struct-style format descriptor */
            2,                                            /* Number of dimensions */
            { self.Height(), self.Width() },              /* Buffer dimensions */
            { sizeof(TSCAL)*self.Width(), sizeof(TSCAL) } /* Strides (in bytes) for each index */
        );
    });
    c.def("NumPy", [] (py::object & self) {
        T& fv = py::cast<T&>(self);
        auto numpy = py::module::import("numpy");
        auto frombuffer = numpy.attr("frombuffer");
        return frombuffer(self, py::detail::npy_format_descriptor<TSCAL>::dtype()).attr("reshape")(fv.Height(),fv.Width());
      }, "Return NumPy object");
}

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
          }, py::arg("inds"), "Return values at given positions"  );
        c.def("__getitem__", [](T &v, py::list ind )-> TNEW {
                int n = py::len(ind);
                TNEW res(n);
                for (int i=0; i<n; i++) {
                    res[i] = v[ ind[i].cast<int>() ];
                }
                return res;
          }, py::arg("ind"), "Return values at given positions"  );
        c.def("__setitem__", [](T &self, py::slice inds, const T & rv ) {
            size_t start, step, n;
            InitSlice( inds, self.Size(), start, step, n );
            for (int i=0; i<n; i++, start+=step)
                self[start] = rv[i];
          }, py::arg("inds"), py::arg("rv"), "Set values at given positions" );
        c.def("__setitem__", [](T &self, py::slice inds, TSCAL val ) {
            size_t start, step, n;
            InitSlice( inds, self.Size(), start, step, n );
            for (int i=0; i<n; i++, start+=step)
                self[start] = val;
          }, py::arg("inds"), py::arg("value"), "Set value at given positions" );
        c.def("__add__" , [](T &self, T &v) { return TNEW(self+v); }, py::arg("vec") );
        c.def("__sub__" , [](T &self, T &v) { return TNEW(self-v); }, py::arg("vec") );
        c.def("__mul__" , [](T &self, TSCAL s) { return TNEW(s*self); }, py::arg("value") );
        c.def("__rmul__" , [](T &self, TSCAL s) { return TNEW(s*self); }, py::arg("value") );
        c.def("__neg__" , [](T &self) { return TNEW(-self); });        
        c.def("InnerProduct",  [](T & x, T & y) { return InnerProduct (x, y); }, py::arg("y"), "Returns InnerProduct with other object");
        c.def("Norm",  [](T & x) { return L2Norm(x); }, "Returns L2-norm");
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
            if(py::isinstance<py::int_>(rows)) {
              py::object row = py::cast( self.Row(rows.cast<int>()) );
              return row.attr("__getitem__")(cols);
            }

            // Second element of tuple is of type int
            if(py::isinstance<py::int_>(cols)) {
              py::object col = py::cast( TROW(self.Col(cols.cast<int>())) );
              return col.attr("__getitem__")(rows);
            }

            // Both elements are slices
            try {
              auto row_slice = rows.cast<py::slice>();
              auto col_slice = cols.cast<py::slice>();
              return py::cast(ColGetSlice(RowGetSlice(self, row_slice), col_slice));
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
            if(py::isinstance<py::int_>(rows)) {
              py::object row = py::cast( self.Row(rows.cast<int>()) );
              row.attr("__setitem__")(cols, v);
              return;
            }

            // Second element of tuple is of type int
            if(py::isinstance<py::int_>(cols)) {
              auto row_slice = rows.cast<py::slice>();
              auto col = self.Col(cols.cast<int>());
              size_t start, step, n;
              InitSlice( row_slice, self.Height(), start, step, n );
              for (int i=0; i<n; i++, start+=step)
                col[start] = v[i];
              return;
            }

            // One of the indices has to be of type int
            cerr << "Invalid Matrix access!" << endl;
          }

          static void SetTupleScal( TMAT & self, py::tuple t, TSCAL val) {
            py::object rows = t[0];
            py::object cols = t[1];

            // First element of tuple is of type int
            if(py::isinstance<py::int_>(rows)) {
              py::object row = py::cast( self.Row(rows.cast<int>()) );
              row.attr("__setitem__")(cols,val);
              return;
            }

            // Second element of tuple is of type int
            if(py::isinstance<py::int_>(cols)) {
              auto row_slice = rows.cast<py::slice>();
              auto col = self.Col(cols.cast<int>());
              size_t start, step, n;
              InitSlice( row_slice, self.Height(), start, step, n );
              for (int i=0; i<n; i++, start+=step)
                col[start] = val;
              return;
            }

            // Both elements are slices
            try {
              py::slice row_slice = rows.cast<py::slice> ();
              size_t start, step, n;
              InitSlice( row_slice, self.Height(), start, step, n );
              for (int i=0; i<n; i++, start+=step) {
                py::object row = py::cast(self.Row(start));
                row.attr("__setitem__")(cols,val);
              }
              return;
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
              /*
              size_t start, step, n;
              InitSlice( row_slice, self.Height(), start, step, n );
              for (int i=0; i<n; i++, start+=step) {
                py::object row = py::cast(self.Row(start));
                py::object f = row.attr("__setitem__");
                f(row, cols, rmat.Row(i));
              }
              */
              size_t rstart, rstep, rn;              
              size_t cstart, cstep, cn;
              InitSlice( row_slice, self.Height(), rstart, rstep, rn );
              InitSlice( col_slice, self.Width(), cstart, cstep, cn );
              for (int i = 0, ii=rstart; i < rn; i++, ii+=rstep)
                for (int j = 0, jj = cstart; j < cn; j++, jj+=cstep)
                  self(ii,jj) = rmat(i,j);

            } catch (py::error_already_set const &) {
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
        c.def("__add__" , [](TMAT &self, TMAT &m) { return TNEW(self+m); }, py::arg("mat") );
        c.def("__sub__" , [](TMAT &self, TMAT &m) { return TNEW(self-m); }, py::arg("mat") );
        c.def("__mul__" , [](TMAT &self, TMAT &m) { return TNEW(self*m); }, py::arg("mat") );
        c.def("__mul__" , [](TMAT &self, FlatVector<TSCAL> &v) { return Vector<TSCAL>(self*v); }, py::arg("vec") );
        c.def("__mul__" , [](TMAT &self, TSCAL s) { return TNEW(s*self); }, py::arg("values") );
        c.def("__rmul__" , [](TMAT &self, TSCAL s) { return TNEW(s*self); }, py::arg("value") );
        c.def("__neg__" , [](TMAT &self) { return TNEW(-self); });
        c.def("Height", &TMAT::Height, "Return height of matrix" );
        c.def("Width", &TMAT::Width, "Return width of matrix" );
        c.def_property_readonly("h", py::cpp_function(&TMAT::Height ), "Height of the matrix");
        c.def_property_readonly("w", py::cpp_function(&TMAT::Width ), "Width of the matrix");
        c.def_property_readonly("T", py::cpp_function([](TMAT &self) { return TNEW(Trans(self)); } ), "return transpose of matrix" );
        c.def_property_readonly("A", py::cpp_function([](TMAT &self) { return Vector<TSCAL>(FlatVector<TSCAL>( self.Width()* self.Height(), &self(0,0)) ); } ), "Returns matrix as vector" );
        c.def("__len__", []( TMAT& self) { return self.Height();}, "Return height of matrix"  );
}


template <typename TVEC, typename TNEW, typename TSCAL>
auto ExportVector(py::module &m, const char * name ) -> py::class_<TVEC>
  {
    auto c = py::class_<TVEC >(m, name, py::buffer_protocol());
    PyDefVector<TVEC, TSCAL>(m, c);
    PyVecAccess< TVEC, TNEW >(m, c);
    PyDefVecBuffer<TVEC>(c);
    c.def(py::self+=py::self);
    c.def(py::self-=py::self);
    c.def(py::self*=TSCAL());
    c.def("__str__", &ToString<TVEC>);
    c.def("__repr__", &ToString<TVEC>);
    return c;
  }

template <typename T, typename TSCAL, typename TPYCLASS>
void ExportImmediateOperators(TPYCLASS &c)
  {
      // "return self;" is important here!
      c.def("__iadd__", [] (T &self, T &rhs) { self+=rhs; return self; });
      c.def("__isub__", [] (T &self, T &rhs) { self-=rhs; return self; });
      c.def("__imul__", [] (T &self, TSCAL &rhs) { self*=rhs; return self; });
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
        .def(py::init<size_t, double *>())
        .def("Range",    static_cast</* const */ FVD (FVD::*)(size_t,size_t) const> (&FVD::Range ) )
      ;
    ExportVector< FVC, VC, Complex>(m, "FlatVectorC")
        .def(py::self*=double())
        .def(py::init<size_t, Complex *>())
        .def("Range",    static_cast</* const */ FVC (FVC::*)(size_t,size_t) const> (&FVC::Range ) )
        ;

    ExportVector< SVD, VD, double>(m, "SliceVectorD")
        .def("Range",    static_cast<const SVD (SVD::*)(size_t,size_t) const> (&SVD::Range ) )
        ;
    ExportVector< SVC, VC, Complex>(m, "SliceVectorC")
        .def("Range",    static_cast<const SVC (SVC::*)(size_t,size_t) const> (&SVC::Range ) )
        .def(py::self*=double())
        ;

    py::class_<VD, FVD> cvd(m, "VectorD", py::buffer_protocol());
    cvd.def(py::init( [] (int n) { return new VD(n); }));
    PyDefVecBuffer<VD>(cvd);
    ExportImmediateOperators<VD, double>(cvd);

    py::class_<VC, FVC > cvc(m, "VectorC", py::buffer_protocol());
    cvc.def(py::init( [] (int n) { return new VC(n); }));
    PyDefVecBuffer<VC>(cvc);
    ExportImmediateOperators<VC, Complex>(cvc);

    m.def("Vector",
            [] (int n, bool is_complex) {
                if(is_complex) {
                  return py::cast(Vector<Complex>(n));
                }
                else return py::cast(Vector<double>(n));
                },
            py::arg("length"),
          py::arg("complex")=false, docu_string(R"raw_string(

Parameters:

length : int
  input length

complex : bool
  input complex values
)raw_string")
           );
//     m.def("Vector",
//             [] (py::list values) {
//                                Vector<> v(len(values));
//                                for (int i = 0; i < v.Size(); i++)
//                                  v(i) = values[i].cast<double>();
//                                return v;
//                              },
//             py::arg("vals"), docu_string(R"raw_string(

// Parameters:

// vals : list
//        input list of values
// )raw_string")
//            );

    m.def("Vector", [] (const std::vector<double> & values)
      {
        Vector<double> v(values.size());
        for (auto i : Range(values.size()))
          v[i] = values[i];
      });

    m.def("Vector", [] (const std::vector<Complex> & values)
      {
        Vector<Complex> v(values.size());
        for (auto i : Range(values.size()))
          v[i] = values[i];
      });

/*
    m.def("Vector",
            [] (py::list values) -> py::object {
                              py::object tmp = values[0];
                              // TODO: How to do that better? What is the complex equivalent to py::_int?
                              // to be able to use something like py::isinstance<py::int_>(tmp)
                              py::object type = py::eval("complex"); 
                              if (py::isinstance(tmp, type)){
                                Vector<Complex> v(len(values));
                                for (int i = 0; i < v.Size(); i++)
                                v(i) = values[i].cast<Complex>();
                                return py::cast(v);
                              } 
                              Vector<> v(len(values));
                              for (int i = 0; i < v.Size(); i++)
                              v(i) = values[i].cast<double>();
                              return py::cast(v);

                             },
            py::arg("vals"), docu_string(R"raw_string(

Parameters:

vals : list
       input list of values
)raw_string")
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
            py::arg("vals"), docu_string(R"raw_string(
Parameters:

vals : tuple
       input tuple of values

)raw_string")
           );
*/

    py::class_<Vec<1>> v1(m, "Vec1D");
    PyVecAccess<Vec<1>>(m, v1);
    PyDefROBracketOperator<Vec<1>, double>(m, v1);

    py::class_<Vec<2>> v2(m, "Vec2D");
    PyVecAccess<Vec<2>>(m, v2);
    PyDefROBracketOperator<Vec<2>, double>(m, v2);

    py::class_<Vec<3>> v3(m, "Vec3D");
    PyVecAccess<Vec<3>>(m, v3);
    PyDefROBracketOperator<Vec<3>, double>(m, v3);

    ///////////////////////////////////////////////////////////////////////////////////////
    // Matrix types
    typedef FlatMatrix<double> FMD;
    py::class_<FlatMatrix<double> > class_FMD(m, "FlatMatrixD", py::buffer_protocol());
        PyMatAccess<FMD, Matrix<double> >(class_FMD);
        class_FMD.def(py::self+=py::self);
        class_FMD.def(py::self-=py::self);
        class_FMD.def(py::self*=double());
        class_FMD.def("Inverse", [](FMD & self, FMD & inv) {
	    CalcInverse(self,inv); return;
	  });
        class_FMD.def_property_readonly("I", py::cpp_function([](FMD &self) { return Inv(self); } ) );        
        class_FMD.def("__str__", &ToString<FMD>);
        class_FMD.def("__repr__", &ToString<FMD>);
        PyDefMatBuffer<FMD>(class_FMD);


    typedef FlatMatrix<Complex> FMC;
    auto class_FMC = py::class_<FlatMatrix<Complex> > (m, "FlatMatrixC", py::buffer_protocol());
        PyMatAccess<FMC, Matrix<Complex> >(class_FMC);
        class_FMC.def("__str__", &ToString<FMC>);
        class_FMC.def("__repr__", &ToString<FMC>);
        class_FMC.def(py::self+=py::self)
        .def(py::self-=py::self)
        .def(py::self*=Complex())
        .def_property("diag",
                py::cpp_function([](const FMC &self) { return Vector<Complex>(self.Diag()); }),
                py::cpp_function([](FMC &self, const FVC &v) { self.Diag() = v; }))
          .def("__add__" , [](FMC &self, FMD &m) { return Matrix<Complex>(self+m); }, py::arg("mat") )
        .def("__sub__" , [](FMC &self, FMD &m) { return Matrix<Complex>(self-m); }, py::arg("mat") )
        .def("__mul__" , [](FMC &self, FMD &m) { return Matrix<Complex>(self*m); }, py::arg("mat") )
        .def("__radd__" , [](FMC &self, FMD &m) { return Matrix<Complex>(self+m); }, py::arg("mat") )
        .def("__rsub__" , [](FMC &self, FMD &m) { return Matrix<Complex>(self-m); }, py::arg("mat") )
        .def("__rmul__" , [](FMC &self, FMD &m) { return Matrix<Complex>(m*self); }, py::arg("mat") )
          .def("__mul__" , [](FMC &self, FVD &v) { return Vector<Complex>(self*v); }, py::arg("vec") )
          .def("__mul__" , [](FMC &self, double s) { return Matrix<Complex>(s*self); }, py::arg("value") )
          .def("__rmul__" , [](FMC &self, double s) { return Matrix<Complex>(s*self); }, py::arg("value") )
          .def("Height", &FMC::Height, "Returns height of the matrix" )
          .def("Width", &FMC::Width, "Returns width of the matrix" )
        .def("__len__", []( FMC& self) { return self.Height();}  )
          .def_property_readonly("h", py::cpp_function(&FMC::Height ), "Height of the matrix")
          .def_property_readonly("w", py::cpp_function(&FMC::Width ), "Width of the matrix")
          .def_property_readonly("A", py::cpp_function([](FMC &self) { return Vector<Complex>(FlatVector<Complex>( self.Width()* self.Height(), &self(0,0) )); }  ), "Returns matrix as vector")
          .def_property_readonly("T", py::cpp_function([](FMC &self) { return Matrix<Complex>(Trans(self)); } ), "Return transpose of matrix" )
        .def_property_readonly("C", py::cpp_function([](FMC &self) { 
            Matrix<Complex> result( self.Height(), self.Width() );
            for (int i=0; i<self.Height(); i++)
                for (int j=0; j<self.Width(); j++) 
                    result(i,j) = Conj(self(i,j));
            return result;
            } ), "Return conjugate matrix" )
        .def_property_readonly("H", py::cpp_function([](FMC &self) { 
            Matrix<Complex> result( self.Width(), self.Height() );
            for (int i=0; i<self.Height(); i++)
                for (int j=0; j<self.Width(); j++) 
                    result(j,i) = Conj(self(i,j));
            return result;
            } ), "Return conjugate and transposed matrix" )
          .def_property_readonly("I", py::cpp_function([](FMC & self) { return Inv(self); }))
        ;
    PyDefMatBuffer<FMC>(class_FMC);

    auto class_MD = py::class_<Matrix<double>, FMD>(m, "MatrixD", py::buffer_protocol())
      .def(py::init( [] (int n, int m) { return new Matrix<double>(n, m); }), py::arg("n"), py::arg("m"), "Makes matrix of dimension n x m")
        ;
    PyDefMatBuffer<Matrix<>>(class_MD);
    ExportImmediateOperators<Matrix<double>, double>(class_MD);

    auto class_MC = py::class_<Matrix<Complex>, FMC >(m, "MatrixC", py::buffer_protocol())
      .def(py::init( [] (int n, int m) { return new Matrix<Complex>(n, m); }), py::arg("n"), py::arg("m"), "Makes matrix of dimension n x m")
        ;
    PyDefMatBuffer<Matrix<Complex>>(class_MC);
    ExportImmediateOperators<Matrix<Complex>, Complex>(class_MC);

    auto class_Mat2D = py::class_<Mat<2,2,double>>(m,"Mat2D", py::buffer_protocol());
    PyDefMatBuffer<Mat<2,2,double>>(class_Mat2D);
    class_Mat2D.def("__getitem__", [](Mat<2,2,double> self, py::tuple i)
                    { return self(i[0].cast<size_t>(),i[1].cast<size_t>()); });
    auto class_Mat2C = py::class_<Mat<2,2,Complex>>(m,"Mat2C", py::buffer_protocol());
    PyDefMatBuffer<Mat<2,2,Complex>>(class_Mat2C);
    class_Mat2C.def("__getitem__", [](Mat<2,2,Complex> self, py::tuple i)
                    { return self(i[0].cast<size_t>(),i[1].cast<size_t>()); });
    auto class_Mat3D = py::class_<Mat<3,3,double>>(m,"Mat3D", py::buffer_protocol());
    PyDefMatBuffer<Mat<3,3,double>>(class_Mat3D);
    class_Mat3D.def("__getitem__", [](Mat<3,3,double> self, py::tuple i)
                    { return self(i[0].cast<size_t>(),i[1].cast<size_t>()); });
    auto class_Mat3C = py::class_<Mat<3,3,Complex>>(m,"Mat3C", py::buffer_protocol());
    PyDefMatBuffer<Mat<3,3,Complex>>(class_Mat3C);
    class_Mat3C.def("__getitem__", [](Mat<3,3,Complex> self, py::tuple i)
                    { return self(i[0].cast<size_t>(),i[1].cast<size_t>()); });

    m.def("Matrix",
            [] (int h, int w, bool is_complex) {
                if(is_complex) return py::cast(Matrix<Complex>(h,w));
                else return py::cast(Matrix<double>(h,w));
                },
            py::arg("height"), 
            py::arg("width"), 
          py::arg("complex")=false, docu_string(R"raw_string(
Creates a matrix of given height and width.

Parameters:

height : int
  input height

width : int
  input width

complex : bool
  input complex values

)raw_string")
           );

    m.def("InnerProduct",
             [] (py::object x, py::object y) -> py::object
          { return py::object(x.attr("InnerProduct")) (y); }, py::arg("x"), py::arg("y"), "Compute InnerProduct");

    m.def("Norm",
             [] (py::object x) -> py::object
          { return py::object(x.attr("Norm")) (); }, py::arg("x"),"Compute Norm");

    m.def("__timing__", &ngbla::Timing, py::arg("what"), py::arg("n"), py::arg("m"), py::arg("k"), py::arg("lapack")=false);
    m.def("CheckPerformance",
             [] (size_t n, size_t m, size_t k)
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


                                {
                                  Timer t("own AddABt");
                                  t.Start();
                                  c = 0.0;
                                  for (int j = 0; j < its; j++)
                                    // c = a * Trans(b) | Lapack;
                                    AddABt (a, b, c);
                                  t.Stop();
                                  cout << "own AddABt GFlops = " << 1e-9 * n*k*m*its / t.GetTime() << endl;
                                }


                                {
                                  Timer t("own MultMatMat");
                                  t.Start();
                                  c = 0.0;
                                  Matrix<> bt = Trans(b);
                                  for (int j = 0; j < its; j++)
                                    // c = a * Trans(b) | Lapack;
                                    MultMatMat (a, bt, c);
                                  t.Stop();
                                  cout << "own MultMatMat GFlops = " << 1e-9 * n*k*m*its / t.GetTime() << endl;
                                }

                                {
                                  Timer t("own AlignedMultMatMat");
                                  t.Start();
                                  Matrix<SIMD<double>> bt(k,m);
                                  Matrix<SIMD<double>> c(n,m);
                                  c = SIMD<double>(0.0);
                                  for (int j = 0; j < its; j++)
                                    // c = a * Trans(b) | Lapack;
                                    MultMatMat (a, bt, c);
                                  t.Stop();
                                  cout << "own AlignedMultMatMat GFlops = " << SIMD<double>::Size()*1e-9 * n*k*m*its / t.GetTime() << endl;
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


                                { // Lapack - Inverse
                                  Matrix<> a(n,n);
                                  a = 1e-5;
                                  for (size_t i : Range(n)) a(i,i) = 1;

                                  size_t ops = n*n*n;
                                  size_t runs = 1e10/ops+1;
                                  LapackInverse (a);
                                  
                                  Timer t("inverse");
                                  t.Start();
                                  for (size_t j = 0; j < runs; j++)
                                    LapackInverse (a);
                                  t.Stop();
                                  cout << "LapackInverse GFlops = " << 1e-9 * ops*runs / t.GetTime() << endl;
                                }

                                { // Lapack - Inverse
                                  Matrix<> a(n,n);
                                  a = 1e-5;
                                  for (size_t i : Range(n)) a(i,i) = 1;

                                  size_t ops = n*n*n;
                                  size_t runs = 1e10/ops+1;
                                  LapackInverseSPD (a);
                                  
                                  Timer t("inverse");
                                  t.Start();
                                  for (size_t j = 0; j < runs; j++)
                                    LapackInverse (a);
                                  t.Stop();
                                  cout << "LapackInverse GFlops = " << 1e-9 * ops*runs / t.GetTime() << endl;
                                }



                                { // LDL 
                                  Matrix<double,ColMajor> a(n,n), ah(n,n);
                                  ah = 1e-5;
                                  for (size_t i : Range(n)) ah(i,i) = 1;

                                  size_t ops = n*n*n/6;
                                  size_t runs = 1e10/ops+1;
                                  a = ah;
                                  CalcLDL (a.Rows(0,n).Cols(0,n));  
                                  
                                  Timer t("LDL");
                                  t.Start();
                                  for (size_t j = 0; j < runs; j++)
                                    {
                                      a = ah;
                                      CalcLDL (a.Rows(0,n).Cols(0,n));
                                    }
                                  t.Stop();
                                  cout << "CalcLDL GFlops = " << 1e-9 * ops*runs / t.GetTime() << endl;
                                }

                                
                              }, py::arg("n"), py::arg("m"), py::arg("k"));
             }

#endif // NGS_PYTHON
