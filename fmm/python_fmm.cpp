#ifdef NGS_PYTHON

#include <python_ngstd.hpp>
#include "fmmoperator.hpp"
#include <python_linalg.hpp>

#include "mptools.hpp"
#include "fmmkernels.hpp"


using namespace ngsfmm;

namespace
{
  void CheckSHIndex(int order, tuple<int,int> nm)
  {
    auto [n, m] = nm;
    if (n < 0 || n > order || m < -n || m > n)
      throw py::index_error("Sphericalharmonics index (n,m)=(" + ToString(n) + "," + ToString(m) + ") out of range for order " + ToString(order));
  }

  static constexpr const char * pointfmm_doc = R"raw_string(
FMM matrix-vector product on point clouds, without mesh or space:
  (op @ q)[j] = sum_i K(x_j, y_i) q[i],  y_i = sources[i], x_j = targets[j], r = |x_j - y_i|.

Arguments:
  sources : (N,3) float array.
  targets : (M,3) float array, or None for targets = sources. Coincident
    pairs x_j == y_i contribute zero, so targets=None drops the self term.
  kernel : "laplace": G = 1/(4 pi r), real, no kappa.
    "helmholtz": G = exp(i kappa r)/(4 pi r), complex, kappa real or complex, nonzero.
  source : "charges": K = G.
    "dipoles": K = dG/dn_y = <x-y,n_i> (1 - i kappa r) exp(i kappa r)/(4 pi r^3),
    n_i = source_normals[i], used as given (not normalised).
    "curl": the field of the currents J_i = q[i] (3-vectors), only with
    kernel="helmholtz", components=3 and target="potential":
      (op @ q)[j] = sum_i curl_x (G J_i) = sum_i grad_x G x J_i
                  = sum_i (i kappa r - 1) exp(i kappa r)/(4 pi r^3) (x_j - y_i) x J_i.
  target : "potential": K as above. "gradient": grad_x K, the gradient with
    respect to the target point.
  source_normals : (N,3) float array, required by and only allowed for "dipoles".
  components : 1 or 3, default 1. Densities per source. With 3, "charges" and
    "dipoles" apply K to three independent densities q[:,c] at once (one tree,
    the dipoles of a source share n_i); "curl" needs 3, the current vector.

op @ q and op(q) take the densities (real for "laplace", real or complex for
"helmholtz") and return a new array:
  components=1: q of shape (N,) -> (M,) for "potential", (M,3) for "gradient"
    with u[j,k] = d/dx_k u(x_j).
  components=3: q of shape (N,3), q[i,c] = component c at source i -> (M,3) for
    "potential" and "curl", (M,3,3) for "gradient" with u[j,c,k] = d/dx_k of
    component c at x_j. A flat q of length 3N (q[3*i+c]) gives the flat result,
    the ravel() of the shaped one.
op.mat is the operator as BaseMatrix on the flat vectors (the ravel() of these
arrays: gradient entry 3*j+k for components=1, 9*j+3*c+k for components=3), and
op.shape = (op.mat.height, op.mat.width) = (M*components, N*components), the
height times 3 for "gradient".

Keyword arguments (the numbers as int, float or numpy scalars):
  fmm_maxdirect : int, default 100. A box with fmm_maxdirect or more sources (targets) is subdivided, down to fmm_maxlevel.
  fmm_minorder : int, default 20. Base expansion order.
  fmm_order_factor : float, default 2.0. order = fmm_minorder + fmm_order_factor*abs(kappa)*r.
  fmm_separation : float, default 2.0. Box-box admissibility factor of the M2L translations.
    Values below sqrt(3) ~ 1.73 treat boxes touching at an edge or corner as separated and
    lose accuracy.
  fmm_eval_separation : float, default 3.0. The same for single-point evaluation.
  fmm_split_kr : float, default 20.0. A leaf with abs(kappa)*r >= fmm_split_kr keeps subdividing.
  fmm_maxlevel : int, default 20. Maximum FMM tree level.
  device : str, default "auto". "auto": the GPU only after `import ngsolve.gpu` and
    if it has float64 or fp32=True, else native. "none": native. "gpu": the GPU, else
    an error; a GPU without float64 (Apple Metal) needs fp32=True. "cpu": the device
    kernels on the host. The device kernels run in float with fp32=True, else in double.
    The GPU layer is not thread-safe: other threads must not use it
    while a device operator is applied.
  fp32 : bool (or numpy.bool_), default False. True: the FMM computes in single
    precision, the multipole and local expansions (Complex32), their translations with
    the Bessel and Hankel functions, and the direct sums between nearby points. The
    points stay double, the densities and results float64/complex128 numpy arrays as
    above, op.mat a BaseMatrix on double/complex vectors. Single precision gives about 6 to 7
    significant digits.
    With kernel="laplace" (run as Helmholtz with kappa 1e-16) single precision needs a point
    cloud of extent about 1e-2 or more: in smaller clouds the error can grow, down to NaN
    below about 1e-6. Rescale the coordinates of such clouds.
use_fmm is not accepted.
)raw_string";

  static constexpr const char * addcharges_doc = R"raw_string(
Add charges c_i at the points y_i, the same as AddCharge(points[i], charges[i]) for all i.
Their potential is
  u(x) = sum_i c_i G(x - y_i),  G(z) = exp(i kappa |z|) / (4 pi |z|).
For the Laplace kernel 1/(4 pi |z|) construct the expansion with a tiny kappa such as
1e-16; kappa = 0 gives NaN after Calc().

  points : (n,3) real array, one point per row, inside the cube with the center of the
    expansion and half side r (points outside make the expansions inaccurate).
  charges : (n,) real or complex array, or one number for all points.

Add all sources before Calc(). The GIL is released while adding: do not run Calc() or
Evaluate() on this expansion in another thread meanwhile.
)raw_string";

  static constexpr const char * adddipoles_doc = R"raw_string(
Add dipoles with directions d_i and strengths s_i at the points y_i, the same as
AddDipole(points[i], directions[i], strengths[i]) for all i. Their potential is the
derivative along d_i with respect to the evaluation point x,
  u(x) = sum_i s_i <d_i, grad_x G(x - y_i)>
       = sum_i s_i <x - y_i, d_i> (i kappa r_i - 1) exp(i kappa r_i) / (4 pi r_i^3),
r_i = |x - y_i|, G as in AddCharges; for Laplace (tiny kappa) -s_i <x - y_i, d_i> / (4 pi r_i^3).

  points : (n,3) real array, one point per row, inside the cube as in AddCharges.
  directions : (n,3) real array, the d_i, used as given (not normalised).
  strengths : (n,) real or complex array, the s_i, or one number for all points.

Add all sources before Calc(). The GIL is released as in AddCharges.
)raw_string";

  static constexpr const char * evaluate_doc = R"raw_string(
The potential u of the sources (see AddCharges, AddDipoles) at points.

  points : (m,3) real array, one point per row -> complex array of shape (m,);
    a single point (3 numbers) -> complex number.

Before Calc() u is summed directly over all sources, after Calc() the multipole
expansions replace the sums over the sources far from the point. A source at the
evaluation point contributes 0. The GIL is released and the points are evaluated in
parallel inside a TaskManager; Evaluate may run in several threads at once, but not
together with adding sources or Calc().
)raw_string";

  static constexpr const char * evaluate3_doc = R"raw_string(
The field u of the currents (see AddCurrent) at points: a current j on the segment from
sp to ep contributes
  u(x) = j int_sp^ep grad_x G(x - y) x dy,  G(z) = exp(i kappa |z|) / (4 pi |z|),
by the midpoint rule with num points on every piece the FMM boxes cut the segment into.

  points : (m,3) real array, one point per row -> complex array of shape (m,3);
    a single point (3 numbers) -> complex array of shape (3,).

Before Calc() u is summed directly, after Calc() the multipole expansions replace the
sums over the currents far from the point. A quadrature point at the evaluation point
contributes 0. The GIL is released and the points are evaluated in parallel inside a
TaskManager; Evaluate may run in several threads at once, but not together with
AddCurrent or Calc().
)raw_string";

  // an FMM_Operator on point clouds and the numpy shapes of its vectors
  struct PointFMM
  {
    shared_ptr<BaseMatrix> mat;
    bool is_complex;
    size_t nsrc, ntgt;
    int components;
    int target_comps;
  };

  Array<Vec<3>> PointsFromNumpy (py::handle obj, const string & name)
  {
    auto arr = py::array::ensure(obj);
    if (!arr || arr.dtype().kind() == 'c' || arr.ndim() != 2 || arr.shape(1) != 3)
      {
        string msg = name + " must be a real (n,3) array";
        if (arr && arr.ndim() == 2 && arr.shape(0) == 3 && arr.shape(1) != 3)
          msg += ", one point per row; got shape (3," + ToString(arr.shape(1)) + "), pass its transpose";
        throw py::value_error(msg);
      }
    py::array_t<double, py::array::c_style | py::array::forcecast> darr(arr);
    auto a = darr.unchecked<2>();
    Array<Vec<3>> pts(a.shape(0));
    for (size_t i = 0; i < pts.Size(); i++)
      {
        pts[i] = Vec<3>(a(i,0), a(i,1), a(i,2));
        if (!std::isfinite(a(i,0)) || !std::isfinite(a(i,1)) || !std::isfinite(a(i,2)))
          throw py::value_error(name + " must be finite, got a non-finite coordinate in row " + ToString(i));
      }
    return pts;
  }

  // n values from a number (broadcast) or a 1-D array of length n
  Array<Complex> ValuesFromNumpy (py::handle obj, size_t n, const string & name)
  {
    auto arr = py::array::ensure(obj);
    if (!arr || string("biufc").find(arr.dtype().kind()) == string::npos)
      throw py::value_error(name + " must be real or complex numbers");
    if (arr.ndim() > 1 || (arr.ndim() == 1 && size_t(arr.shape(0)) != n))
      throw py::value_error(name + " must be a number or a 1-D array of length " + ToString(n)
                            + " (one value per point)");
    py::array_t<Complex, py::array::c_style | py::array::forcecast> carr(arr);
    Array<Complex> vals(n);
    if (carr.ndim() == 0)
      vals = *carr.data();
    else
      for (size_t i = 0; i < n; i++)
        vals[i] = carr.data()[i];
    return vals;
  }

  template <typename T>
  py::object EvaluateExpansion (const SingularMLExpansion<T> & self, py::handle points, const string & name)
  {
    auto arr = py::array::ensure(points);
    bool single = arr && arr.ndim() == 1;
    if (single && arr.shape(0) != 3)
      throw py::value_error(name + ": a single point must have 3 coordinates, got " + ToString(arr.shape(0)));
    Array<Vec<3>> pts = PointsFromNumpy(single ? py::object(arr.reshape({1, 3})) : py::reinterpret_borrow<py::object>(points),
                                        name + ": points");
    Array<T> vals(pts.Size());
    {
      py::gil_scoped_release rel;
      ParallelFor (pts.Size(), [&](size_t i) { vals[i] = self.Evaluate(pts[i]); }, TasksPerThread(4));
    }
    if constexpr (std::is_same_v<T,Complex>)
      {
        if (single)
          return py::cast(vals[0]);
        py::array_t<Complex> res(vals.Size());
        for (size_t i = 0; i < vals.Size(); i++)
          res.mutable_data()[i] = vals[i];
        return res;
      }
    else
      {
        std::vector<py::ssize_t> shape { py::ssize_t(vals.Size()), 3 };
        if (single)
          shape.erase(shape.begin());
        py::array_t<Complex> res(shape);
        for (size_t i = 0; i < vals.Size(); i++)
          for (int k = 0; k < 3; k++)
            res.mutable_data()[3*i+k] = vals[i](k);
        return res;
      }
  }

  FMM_Parameters PointFMMParameters (const py::kwargs & kwargs, bool fp32)
  {
    FMM_Parameters params;
    auto number = [&kwargs] (const char * name, double value) -> double
    {
      if (!kwargs.contains(name))
        return value;
      py::object obj = kwargs[name];
      if (PyBool_Check(obj.ptr()) || !py::isinstance(obj, py::module_::import("numbers").attr("Real")))
        throw py::type_error(string("FMM: ") + name + " must be an int or float, got "
                             + py::cast<string>(py::type::handle_of(obj).attr("__name__")));
      try
        {
          return double(py::float_(obj));
        }
      catch (py::error_already_set & e)
        {
          if (!e.matches(PyExc_OverflowError))
            throw;
          throw py::value_error(string("FMM: ") + name + " is out of range");
        }
    };
    auto number_in = [&number, &kwargs] (const char * name, double value, double lo, bool strict, const char * what) -> double
    {
      double v = number(name, value);
      if (!kwargs.contains(name))
        return v;
      if (!std::isfinite(v) || (strict ? v <= lo : v < lo))
        throw py::value_error(string("FMM: ") + name + " must be finite" + what + ", got "
                              + py::cast<string>(py::str(kwargs[name])));
      return v;
    };
    auto integer = [&number, &kwargs] (const char * name, int value) -> int
    {
      double v = number(name, value);
      if (!(v > double(std::numeric_limits<int>::min()) - 1 && v < double(std::numeric_limits<int>::max()) + 1))
        throw py::value_error(string("FMM: ") + name + " = " + py::cast<string>(py::str(kwargs[name]))
                              + " does not fit an int");
      return int(v);
    };
    params.maxdirect = integer("fmm_maxdirect", params.maxdirect);
    if (params.maxdirect < 2)
      throw py::value_error("FMM: fmm_maxdirect must be >= 2, got " + py::cast<string>(py::str(kwargs["fmm_maxdirect"])));
    params.minorder = integer("fmm_minorder", params.minorder);
    if (params.minorder < 0)
      throw py::value_error("FMM: fmm_minorder must be >= 0, got " + py::cast<string>(py::str(kwargs["fmm_minorder"])));
    params.order_factor = number_in("fmm_order_factor", params.order_factor, 0, false, " and >= 0");
    params.separation = number_in("fmm_separation", params.separation, 1, true, " and > 1");
    params.eval_separation = number_in("fmm_eval_separation", params.eval_separation, 0, true, " and > 0");
    params.split_kr = number_in("fmm_split_kr", params.split_kr, 0, true, " and > 0");
    params.maxlevel = integer("fmm_maxlevel", params.maxlevel);

    string device = "auto";
    if (kwargs.contains("device"))
      {
        py::object obj = kwargs["device"];
        if (!py::isinstance<py::str>(obj))
          throw py::type_error("FMM: device must be a str ('auto', 'gpu', 'cpu', 'none'), got "
                               + py::cast<string>(py::type::handle_of(obj).attr("__name__")));
        device = py::cast<string>(obj);
        if (device != "auto" && device != "gpu" && device != "cpu" && device != "none")
          throw py::value_error("FMM: device must be one of 'auto', 'gpu', 'cpu', 'none', got '" + device + "'");
      }
    if (device == "auto" && !py::module_::import("sys").attr("modules").contains("ngsolve.gpu"))
      device = "none";
    params.device = ResolveFMMDevice(device, fp32, /*use_fmm*/ true);
    // the device translations need expansion orders >= 1 (order = minorder + order_factor kappa r)
    if (params.minorder < 1)
      params.device = FMM_Device::None;
    return params;
  }

  template <int COMPS, typename TVAL, typename TFMM, typename TKAPPA>
  void MakePointFMMInterfaces (TKAPPA kappa, bool dipoles, bool curl, bool gradient,
                               shared_ptr<BaseFMMInterface> & src, shared_ptr<BaseFMMInterface> & tgt)
  {
    if constexpr (COMPS == 3 && IsComplex<TVAL>())
      {
        if (curl)
          {
            // also the target, as in MaxwellDLKernel, so that MultTrans is the transpose
            src = make_shared<MaxwellCurlDipoles<TKAPPA,TFMM>> (kappa);
            tgt = make_shared<MaxwellCurlDipoles<TKAPPA,TFMM>> (kappa);
            return;
          }
      }
    if (dipoles)
      src = make_shared<Dipoles<COMPS,TVAL,TFMM,TKAPPA>> (kappa);
    else
      src = make_shared<Charges<COMPS,TVAL,TFMM,TKAPPA>> (kappa);
    if (gradient)
      tgt = make_shared<GradientEval<COMPS,TVAL,TFMM,TKAPPA>> (kappa);
    else
      tgt = make_shared<Charges<COMPS,TVAL,TFMM,TKAPPA>> (kappa);
  }

  template <typename TVAL, typename TKAPPA>
  shared_ptr<BaseMatrix> MakePointFMMOperator (const string & name, TKAPPA kappa, int comps, bool fp32, bool dipoles,
                                               bool curl, bool gradient, Array<Vec<3>> xpts, Array<Vec<3>> ypts,
                                               Array<Vec<3>> xnv, const FMM_Parameters & params)
  {
    using TVAL32 = std::conditional_t<IsComplex<TVAL>(), Complex32, float>;
    shared_ptr<BaseFMMInterface> src, tgt;
    if (comps == 3 && fp32)
      MakePointFMMInterfaces<3,TVAL,TVAL32> (kappa, dipoles, curl, gradient, src, tgt);
    else if (comps == 3)
      MakePointFMMInterfaces<3,TVAL,TVAL> (kappa, dipoles, curl, gradient, src, tgt);
    else if (fp32)
      MakePointFMMInterfaces<1,TVAL,TVAL32> (kappa, dipoles, curl, gradient, src, tgt);
    else
      MakePointFMMInterfaces<1,TVAL,TVAL> (kappa, dipoles, curl, gradient, src, tgt);
    Array<Vec<3>> ynv(ypts.Size());
    ynv = Vec<3>(0.0);
    return make_shared<FMM_Operator<TVAL>> (name, IVec<2>(gradient ? 3*comps : comps, comps), src, tgt,
                                            std::move(xpts), std::move(ypts), std::move(xnv), std::move(ynv), params);
  }

  template <typename T>
  py::array ApplyPointFMM (const PointFMM & self, const py::array & q, const std::vector<py::ssize_t> & shape)
  {
    py::array_t<T, py::array::c_style | py::array::forcecast> tq(q);
    const T * qv = tq.data();
    VVector<T> x(self.nsrc*self.components), y(self.ntgt*self.target_comps*self.components);
    for (size_t i = 0; i < x.Size(); i++)
      x(i) = qv[i];
    {
      py::gil_scoped_release rel;
      self.mat->Mult(x, y);
    }
    py::array_t<T> res(shape);
    T * data = res.mutable_data();
    auto fy = y.FV();
    for (size_t i = 0; i < fy.Size(); i++)
      data[i] = fy(i);
    return res;
  }

  py::array PointFMMApplyNumpy (const PointFMM & self, py::handle q)
  {
    auto arr = py::array::ensure(q);
    if (arr && string("biufc").find(arr.dtype().kind()) == string::npos)
      throw py::value_error("FMM: the density must be real or complex numbers");
    std::vector<py::ssize_t> shape { py::ssize_t(self.ntgt) };
    if (self.components == 1)
      {
        if (!arr || arr.ndim() != 1 || size_t(arr.shape(0)) != self.nsrc)
          throw py::value_error("FMM: the density must be a 1-D array of length " + ToString(self.nsrc)
                                + " (one value per source)");
        if (self.target_comps > 1)
          shape.push_back(self.target_comps);
      }
    else
      {
        bool shaped = arr && arr.ndim() == 2 && size_t(arr.shape(0)) == self.nsrc && arr.shape(1) == 3;
        bool flat = arr && arr.ndim() == 1 && size_t(arr.shape(0)) == 3*self.nsrc;
        if (!shaped && !flat)
          {
            string msg = "FMM: with components=3 the density must be an (" + ToString(self.nsrc)
              + ",3) array or a 1-D array of length " + ToString(3*self.nsrc) + " (entry 3*i+c = component c of source i)";
            if (arr && arr.ndim() == 2 && arr.shape(0) == 3 && size_t(arr.shape(1)) == self.nsrc)
              msg += "; got shape (3," + ToString(self.nsrc) + "), pass its transpose";
            throw py::value_error(msg);
          }
        if (flat)
          shape[0] *= 3*self.target_comps;
        else
          {
            shape.push_back(3);
            if (self.target_comps > 1)
              shape.push_back(self.target_comps);
          }
      }
    if (arr.dtype().kind() == 'c' && !self.is_complex)
      throw py::value_error("FMM: complex density for the real kernel 'laplace'; "
                            "apply the operator to the real and imaginary parts separately");
    if (self.is_complex)
      return ApplyPointFMM<Complex> (self, arr, shape);
    return ApplyPointFMM<double> (self, arr, shape);
  }
}



void NGS_DLL_HEADER ExportNgsfmm(py::module & m)
{
  py::class_<SphericalHarmonics<Complex>> (m, "Sphericalharmonics")
    .def_property_readonly("order", [](SphericalHarmonics<Complex>& self) { return self.Order(); })
    .def("__setitem__", [](SphericalHarmonics<Complex>& self, tuple<int,int> nm, Complex val)
    { CheckSHIndex(self.Order(), nm); self.Coef(get<0>(nm), get<1>(nm)) = val; })
    .def("__getitem__", [](SphericalHarmonics<Complex>& self, tuple<int,int> nm)
    { CheckSHIndex(self.Order(), nm); return self.Coef(get<0>(nm), get<1>(nm)); })
    .def_property_readonly("coefs",
                           py::cpp_function([](SphericalHarmonics<Complex>& self) { return self.Coefs(); },
                                            py::keep_alive<0,1>()),
                           "coefficient vector")
    .def("RotateZ", [](SphericalHarmonics<Complex>& self, double alpha) { self.RotateZ(alpha); })
    .def("RotateY", [](SphericalHarmonics<Complex>& self, double alpha) { self.RotateY(alpha); })
    .def("FlipZ", [](SphericalHarmonics<Complex>& self) { self.FlipZ(); })    
    ;

  py::class_<SingularMLExpansion<Complex>, shared_ptr<SingularMLExpansion<Complex>>> (m, "SingularMLExpansion")
    .def(py::init<Vec<3>,double,double>())
    .def("AddCharge", &SingularMLExpansion<Complex>::AddCharge)
    .def("AddDipole", &SingularMLExpansion<Complex>::AddDipole)
    .def("AddCharges", [](SingularMLExpansion<Complex> & self, py::handle points, py::handle charges)
    {
      Array<Vec<3>> pts = PointsFromNumpy(points, "SingularMLExpansion.AddCharges: points");
      Array<Complex> c = ValuesFromNumpy(charges, pts.Size(), "SingularMLExpansion.AddCharges: charges");
      py::gil_scoped_release rel;
      for (size_t i = 0; i < pts.Size(); i++)
        self.AddCharge(pts[i], c[i]);
    }, py::arg("points"), py::arg("charges"), addcharges_doc)
    .def("AddDipoles", [](SingularMLExpansion<Complex> & self, py::handle points, py::handle directions,
                          py::handle strengths)
    {
      Array<Vec<3>> pts = PointsFromNumpy(points, "SingularMLExpansion.AddDipoles: points");
      Array<Vec<3>> dirs = PointsFromNumpy(directions, "SingularMLExpansion.AddDipoles: directions");
      if (dirs.Size() != pts.Size())
        throw py::value_error("SingularMLExpansion.AddDipoles: directions must have the shape of points, ("
                              + ToString(pts.Size()) + ",3)");
      Array<Complex> c = ValuesFromNumpy(strengths, pts.Size(), "SingularMLExpansion.AddDipoles: strengths");
      py::gil_scoped_release rel;
      for (size_t i = 0; i < pts.Size(); i++)
        self.AddDipole(pts[i], dirs[i], c[i]);
    }, py::arg("points"), py::arg("directions"), py::arg("strengths"), adddipoles_doc)
    .def("Evaluate", [](const SingularMLExpansion<Complex> & self, py::handle points)
    { return EvaluateExpansion(self, points, "SingularMLExpansion.Evaluate"); }, py::arg("points"), evaluate_doc)
    .def("Calc", &SingularMLExpansion<Complex>::CalcMP)
    .def("Norm", &SingularMLExpansion<Complex>::Norm)    
    .def("__str__", [](SingularMLExpansion<Complex>& mlmp) { return ToString<>(mlmp); })
    .def("Print", [](const SingularMLExpansion<Complex> &self) {
        py::scoped_ostream_redirect stream(
            std::cout,
            py::module_::import("sys").attr("stdout")
        );
        self.Print(std::cout);
    })
    ;

  py::class_<SingularMLExpansion<Vec<3,Complex>>, shared_ptr<SingularMLExpansion<Vec<3,Complex>>>> (m, "SingularMLExpansion3")
    .def(py::init<Vec<3>,double,double>())
    .def("AddCurrent", &SingularMLExpansion<Vec<3,Complex>>::AddCurrent, py::arg("sp"), py::arg("ep"), py::arg("j"), py::arg("num")=100)
    .def("Evaluate", [](const SingularMLExpansion<Vec<3,Complex>> & self, py::handle points)
    { return EvaluateExpansion(self, points, "SingularMLExpansion3.Evaluate"); }, py::arg("points"), evaluate3_doc)
    .def("Calc", &SingularMLExpansion<Vec<3,Complex>>::CalcMP)
    ;

  
  py::class_<RegularMLExpansion<Complex>, shared_ptr<RegularMLExpansion<Complex>>> (m, "RegularMLExpansion")
    .def("Norm", &RegularMLExpansion<Complex>::Norm)    
    .def("__str__", [](RegularMLExpansion<Complex>& mlmp) { return ToString<>(mlmp); })
    .def("Print", [](const RegularMLExpansion<Complex> &self) {
        py::scoped_ostream_redirect stream(
            std::cout,
            py::module_::import("sys").attr("stdout")
        );
        self.Print(std::cout);
    })
    ;

  py::class_<PointFMM, shared_ptr<PointFMM>> (m, "FMM", pointfmm_doc)
    .def(py::init([](py::object sources, py::object targets, string kernel,
                     optional<std::variant<double,Complex>> kappa, string source, string target,
                     py::object source_normals, py::object components, py::kwargs kwargs)
    {
      const char * fmm_options[] = { "fmm_maxdirect", "fmm_minorder", "fmm_order_factor", "fmm_separation",
                                     "fmm_eval_separation", "fmm_split_kr", "fmm_maxlevel", "device", "fp32" };
      for (auto item : kwargs)
        {
          string key = py::cast<string>(item.first);
          if (key == "use_fmm")
            throw py::type_error("FMM: 'use_fmm' is not supported, the operator always uses the multipole expansions");
          if (std::find(std::begin(fmm_options), std::end(fmm_options), key) == std::end(fmm_options))
            throw py::type_error("FMM: unexpected keyword argument '" + key + "'");
        }

      bool fp32 = false;
      if (kwargs.contains("fp32"))
        {
          py::object obj = kwargs["fp32"];
          if (!PyBool_Check(obj.ptr()) && !py::isinstance(obj, py::module_::import("numpy").attr("bool_")))
            throw py::type_error("FMM: fp32 must be True or False, got " + py::cast<string>(py::repr(obj)));
          fp32 = py::cast<bool>(obj);
        }

      if (kernel != "laplace" && kernel != "helmholtz")
        throw py::value_error("FMM: kernel must be 'laplace' or 'helmholtz', got '" + kernel + "'");
      if (source != "charges" && source != "dipoles" && source != "curl")
        throw py::value_error("FMM: source must be 'charges', 'dipoles' or 'curl', got '" + source + "'");
      if (target != "potential" && target != "gradient")
        throw py::value_error("FMM: target must be 'potential' or 'gradient', got '" + target + "'");
      bool laplace = kernel == "laplace", dipoles = source == "dipoles", curl = source == "curl",
        gradient = target == "gradient";

      bool integral = !PyBool_Check(components.ptr())
        && py::isinstance(components, py::module_::import("numbers").attr("Integral"));
      int comps = integral && components.equal(py::int_(1)) ? 1 : integral && components.equal(py::int_(3)) ? 3 : 0;
      if (!comps)
        throw py::value_error("FMM: components must be 1 or 3, got " + py::cast<string>(py::repr(components)));
      if (curl && (laplace || comps != 3 || gradient))
        throw py::value_error("FMM: source='curl' needs kernel='helmholtz', components=3 and target='potential', got kernel='"
                              + kernel + "', components=" + ToString(comps) + ", target='" + target + "'");

      bool zero_kappa = !kappa || std::visit([](auto k) { return k == 0.0; }, *kappa);
      if (laplace && !zero_kappa)
        throw py::value_error("FMM: kernel='laplace' takes no kappa, use kernel='helmholtz'");
      if (!laplace && zero_kappa)
        throw py::value_error("FMM: kernel='helmholtz' needs a nonzero kappa, use kernel='laplace' for kappa=0");

      Array<Vec<3>> xpts = PointsFromNumpy(sources, "FMM: sources");
      Array<Vec<3>> ypts = targets.is_none() ? Array<Vec<3>>(xpts) : PointsFromNumpy(targets, "FMM: targets");
      Array<Vec<3>> xnv(xpts.Size());
      if (dipoles)
        {
          if (source_normals.is_none())
            throw py::value_error("FMM: source='dipoles' needs source_normals");
          xnv = PointsFromNumpy(source_normals, "FMM: source_normals");
          if (xnv.Size() != xpts.Size())
            throw py::value_error("FMM: source_normals must have the shape of sources");
        }
      else
        {
          if (!source_normals.is_none())
            throw py::value_error("FMM: source_normals are only used with source='dipoles'");
          xnv = Vec<3>(0.0);
        }

      FMM_Parameters params = PointFMMParameters(kwargs, fp32);
      string name = "FMM " + kernel + " " + source + " " + target;

      auto self = make_shared<PointFMM>();
      self->is_complex = !laplace;
      self->nsrc = xpts.Size();
      self->ntgt = ypts.Size();
      self->components = comps;
      self->target_comps = gradient ? 3 : 1;
      if (laplace)
        self->mat = MakePointFMMOperator<double,double> (name, 1e-16, comps, fp32, dipoles, curl, gradient, std::move(xpts),
                                                        std::move(ypts), std::move(xnv), params);
      else if (std::holds_alternative<double>(*kappa))
        self->mat = MakePointFMMOperator<Complex,double> (name, std::get<double>(*kappa), comps, fp32, dipoles, curl,
                                                         gradient, std::move(xpts), std::move(ypts), std::move(xnv), params);
      else
        self->mat = MakePointFMMOperator<Complex,Complex> (name, std::get<Complex>(*kappa), comps, fp32, dipoles, curl,
                                                          gradient, std::move(xpts), std::move(ypts), std::move(xnv), params);
      return self;
    }), py::arg("sources"), py::arg("targets")=py::none(), py::arg("kernel")="laplace",
        py::arg("kappa")=py::none(), py::arg("source")="charges", py::arg("target")="potential",
        py::arg("source_normals")=py::none(), py::arg("components")=1)
    .def("__matmul__", &PointFMMApplyNumpy, py::arg("q"))
    .def("__call__", &PointFMMApplyNumpy, py::arg("q"))
    .def_property_readonly("mat", [](const PointFMM & self) { return self.mat; })
    .def_property_readonly("shape", [](const PointFMM & self)
                           { return py::make_tuple(self.mat->Height(), self.mat->Width()); })
    ;
}

#endif // NGS_PYTHON
