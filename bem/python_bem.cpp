#ifdef NGS_PYTHON
#include <regex>
#include <variant>

#include <pybind11/warnings.h>

#include "../ngstd/python_ngstd.hpp"
#include "python_comp.hpp"

#include "mptools.hpp"
#include "potentialtools.hpp"
#include "ngbem.hpp"
#include "mp_coefficient.hpp"


using namespace ngsbem;

namespace
{
  void WarnDeprecated(const string & oldname, const string & replacement)
  {
    string message = oldname + " is deprecated; use " + replacement + " instead.";
    py::warnings::warn(message.c_str(), PyExc_FutureWarning, 1);
  }

  template <class KernelReal, class KernelComplex>
  shared_ptr<BasePotentialOperator>
  MakePotentialFromVariantKappa(shared_ptr<ProxyFunction> proxy,
                                VorB source_vb,
                                optional<Region> definedon,
                                shared_ptr<DifferentialOperator> eval,
                                IntOp_Parameters ioparams,
                                int intorder,
                                const std::variant<double, Complex> & kappa)
  {
    if (std::holds_alternative<double>(kappa))
      return make_shared<PotentialOperator<KernelReal>>(proxy, source_vb, definedon, eval,
                                                        KernelReal(std::get<double>(kappa)),
                                                        ioparams, intorder);

    return make_shared<PotentialOperator<KernelComplex>>(proxy, source_vb, definedon, eval,
                                                         KernelComplex(std::get<Complex>(kappa)),
                                                         ioparams, intorder);
  }

  inline py::dict FMMInfoToDict (const FMMOperatorInfo & info)
  {
    double abs_kappa = std::visit([](auto k) { return std::abs(k); }, info.kappa);
    bool has_wavelength = abs_kappa > 1e-14;
    double wavelength = has_wavelength ? 2.0 * M_PI / abs_kappa : 0.0;
    auto maybe = [] (bool condition, auto value) -> py::object {
      return condition ? py::cast(value) : py::none();
    };
    auto mean_or_none = [] (double sum, size_t count) -> py::object {
      return count > 0 ? py::cast(sum / double(count)) : py::none();
    };

    auto tree_to_dict = [&](const FMMTreeStats & s, double bbox_radius) {
      py::dict d;
      d["depth"] = s.max_level;
      d["num_nodes"] = s.num_nodes;
      d["num_leaves"] = s.num_leaves;
      d["active_leaves"] = s.active_leaves;
      d["active_leaves_fraction"] = (s.num_leaves > 0) ? double(s.active_leaves) / double(s.num_leaves) : 0.0;
      d["bbox_radius"] = bbox_radius;
      d["diameter_in_wavelengths"] = maybe(has_wavelength, 2.0 * bbox_radius / wavelength);
      py::list npl;
      for (auto n : s.nodes_per_level) npl.append(n);
      d["nodes_per_level"] = npl;
      d["leaf_size_min"] = maybe(s.num_leaves > 0, s.leaf_size_min);
      d["leaf_size_max"] = maybe(s.num_leaves > 0, s.leaf_size_max);
      d["leaf_size_mean"] = mean_or_none(s.leaf_size_sum, s.num_leaves);
      d["order_min"] = maybe(s.num_allocated_multipoles > 0, s.order_min);
      d["order_max"] = maybe(s.num_allocated_multipoles > 0, s.order_max);
      d["order_mean"] = mean_or_none(s.order_sum, s.num_allocated_multipoles);
      d["num_allocated_multipoles"] = s.num_allocated_multipoles;
      d["total_coefficients"] = s.total_coefficients;
      d["multipole_mb"] = s.multipole_bytes / 1.0e6;
      return d;
    };

    auto params_to_dict = [](const FMM_Parameters & p) {
      py::dict d;
      d["fmm_maxdirect"] = p.maxdirect;
      d["fmm_minorder"] = p.minorder;
      d["fmm_order_factor"] = p.order_factor;
      d["fmm_separation"] = p.separation;
      d["fmm_eval_separation"] = p.eval_separation;
      d["fmm_split_kr"] = p.split_kr;
      d["fmm_maxlevel"] = p.maxlevel;
      return d;
    };

    py::dict d;
    d["kernel_name"] = info.kernel_name;
    d["source_size"] = info.source_size;
    d["target_size"] = info.target_size;
    d["source_dofs"] = info.source_dofs;
    d["target_dofs"] = info.target_dofs;
    if (std::holds_alternative<Complex>(info.kappa))
      d["kappa"] = std::get<Complex>(info.kappa);
    else
      d["kappa"] = std::get<double>(info.kappa);
    d["wavelength"] = maybe(has_wavelength, wavelength);
    d["total_num_nodes"] = info.total_num_nodes;
    d["total_num_leaves"] = info.total_num_leaves;
    d["total_multipole_coefficients"] = info.total_multipole_coefficients;
    d["multipole_memory_mb"] = info.total_memory_bytes / 1.0e6;
    d["num_s2r"] = info.num_s2r;
    d["num_direct_evaluations"] = info.num_direct_evaluations;
    d["direct_fallback_fraction"] = info.direct_fallback_fraction;
    d["nearfield_nze"] = info.nearfield_nze;
    d["nearfield_fraction"] = info.nearfield_fraction;
    d["source"] = tree_to_dict(info.source_tree, info.source_bbox_radius);
    d["target"] = tree_to_dict(info.target_tree, info.target_bbox_radius);
    d["parameters"] = params_to_dict(info.parameters);
    return d;
  }

  static constexpr const char * bem_operator_kwargs_doc = R"raw_string(
Keyword arguments:
  use_fmm : bool, default True
    Enable FMM acceleration.
  fmm_maxdirect : int, default 100
    Maximum number of direct source or target contributions in an FMM leaf.
  fmm_minorder : int, default 20
    Base spherical expansion order.
  fmm_order_factor : float, default 2.0
    Frequency-dependent order factor: order = fmm_minorder + fmm_order_factor*abs(kappa)*r.
  fmm_separation : float, default 2.0
    Box-box admissibility factor for M2L translations: a source/target box pair
    is treated as far if dist(centers) > fmm_separation*(r_source+r_target).
    Raise for more accurate translations at the cost of more near-field work.
  fmm_eval_separation : float, default 3.0
    Same as fmm_separation but for single-point evaluation;
  fmm_split_kr : float, default 5.0
    A leaf with abs(kappa)*r >= fmm_split_kr keeps subdividing even when it
    holds fewer than fmm_maxdirect points.
  fmm_maxlevel : int, default 20
    Maximum FMM tree level.
)raw_string";
}



void NGS_DLL_HEADER ExportNgsbem(py::module &m)
{

  // ****************************   Multipole stuff *********************************
  
  py::class_<SphericalHarmonics<Complex>> (m, "Sphericalharmonics")
    .def_property_readonly("order", [](SphericalHarmonics<Complex>& self) { return self.Order(); })
    .def("__setitem__", [](SphericalHarmonics<Complex>& self, tuple<int,int> nm, Complex val)
    { self.Coef(get<0>(nm), get<1>(nm)) = val; })
    .def("__getitem__", [](SphericalHarmonics<Complex>& self, tuple<int,int> nm)
    { return self.Coef(get<0>(nm), get<1>(nm)); })
    .def_property_readonly("coefs",
                           [](SphericalHarmonics<Complex>& self) { return self.Coefs(); },
                           "coefficient vector")
    .def("RotateZ", [](SphericalHarmonics<Complex>& self, double alpha) { self.RotateZ(alpha); })
    .def("RotateY", [](SphericalHarmonics<Complex>& self, double alpha) { self.RotateY(alpha); })
    .def("FlipZ", [](SphericalHarmonics<Complex>& self) { self.FlipZ(); })    
    ;

  py::class_<SingularMLExpansion<Complex>, shared_ptr<SingularMLExpansion<Complex>>> (m, "SingularMLExpansion")
    .def(py::init<Vec<3>,double,double>())
    .def("AddCharge", &SingularMLExpansion<Complex>::AddCharge)
    .def("AddDipole", &SingularMLExpansion<Complex>::AddDipole)
    .def("AddChargeDensity", [](SingularMLExpansion<Complex> & mp, shared_ptr<CoefficientFunction> charge,
                                ngcomp::Region reg) { AddChargeDensity(mp,charge,reg); })
    
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
    .def("AddCurrentDensity", [](SingularMLExpansion<Vec<3,Complex>> & mp, shared_ptr<CoefficientFunction> current,
                                 ngcomp::Region reg) { AddCurrentDensity(mp,current,reg); })
    
    .def("Calc", &SingularMLExpansion<Vec<3,Complex>>::CalcMP)
    // .def("Norm", &SingularMLExpansion<Complex>::Norm)    
    // .def("__str__", [](SingularMLExpansion<Complex>& mlmp) { return ToString<>(mlmp); })
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

  
  py::class_<SphericalHarmonicsCF, shared_ptr<SphericalHarmonicsCF>, CoefficientFunction> (m, "SphericalHarmonicsCF")
    .def(py::init<int>())
    .def_property_readonly("sh", [](SphericalHarmonicsCF& self) -> SphericalHarmonics<Complex>& { return self.SH(); })
    ;

  py::class_<SphericalExpansionCF<Regular>, shared_ptr<SphericalExpansionCF<Regular>>, CoefficientFunction> (m, "RegularExpansionCF")
    .def(py::init<int,double,Vec<3>,double>(), py::arg("order"), py::arg("kappa"),  py::arg("center"), py::arg("rad")=1.0)    
    .def_property_readonly("sh", [](SphericalExpansionCF<Regular>& self) -> SphericalHarmonics<Complex>& { return self.SH(); })
    .def("AddPlaneWave", [](SphericalExpansionCF<Regular>& self, Vec<3> d, Complex c) { self.MP().AddPlaneWave(d, c); })
    .def("In2Out", [](SphericalExpansionCF<Regular>& self, SphericalExpansionCF<Singular>& other, double r) { self.MP().In2Out(other.MP(), r); },
         py::arg("sing"), py::arg("rad"))
    .def("ShiftZ", [](SphericalExpansionCF<Regular>& self, double z, SphericalExpansionCF<Regular> & target) { self.ShiftZ(z, target.MP()); })
    .def("Transform", [](SphericalExpansionCF<Regular>& self, SphericalExpansionCF<Regular> & target) { self.Transform(target); })
    .def("TransformAdd", [](SphericalExpansionCF<Regular>& self, SphericalExpansionCF<Regular> & target) { self.TransformAdd(target); })
    .def("Spectrum", [](SphericalExpansionCF<Regular>& self, bool scaled) { return self.MP().Spectrum(scaled); }, py::arg("scaled"))    
    ;

  py::class_<SphericalExpansionCF<Singular>, shared_ptr<SphericalExpansionCF<Singular>>, CoefficientFunction> (m, "SingularExpansionCF")
    .def(py::init<int,double,Vec<3>,double>(), py::arg("order"), py::arg("kappa"),  py::arg("center"), py::arg("rad")=1.0)
    .def_property_readonly("sh", [](SphericalExpansionCF<Singular>& self) -> SphericalHarmonics<Complex>& { return self.SH(); })
    .def("AddCharge", [](SphericalExpansionCF<Singular>& self, Vec<3> x, Complex c) { self.MP().AddCharge(x, c); })
    .def("AddDipole", [](SphericalExpansionCF<Singular>& self, Vec<3> x, Vec<3> d, Complex c) { self.MP().AddDipole(x, d, c); })
    .def("AddPlaneWave", [](SphericalExpansionCF<Singular>& self, Vec<3> d, Complex c) { self.MP().AddPlaneWave(d, c); })    
    .def("ShiftZ", [](SphericalExpansionCF<Singular>& self, double z, SphericalExpansionCF<Regular> & target) { self.ShiftZ(z, target.MP()); })
    .def("ShiftZ", [](SphericalExpansionCF<Singular>& self, double z, SphericalExpansionCF<Singular> & target) { self.ShiftZ(z, target.MP()); })        
    .def("Transform", [](SphericalExpansionCF<Singular>& self, SphericalExpansionCF<Regular> & target) { self.Transform(target); })
    .def("Transform", [](SphericalExpansionCF<Singular>& self, SphericalExpansionCF<Singular> & target) { self.Transform(target); })
    .def("TransformAdd", [](SphericalExpansionCF<Singular>& self, SphericalExpansionCF<Regular> & target) { self.TransformAdd(target); })
    .def("TransformAdd", [](SphericalExpansionCF<Singular>& self, SphericalExpansionCF<Singular> & target) { self.TransformAdd(target); })
    .def("Spectrum", [](SphericalExpansionCF<Singular>& self, bool scaled) { return self.MP().Spectrum(scaled); }, py::arg("scaled"))
    ;

  py::class_<SphericalExpansionCF<Singular,Vec<3,Complex>>, shared_ptr<SphericalExpansionCF<Singular,Vec<3,Complex>>>, CoefficientFunction> (m, "BiotSavartCF")
    .def(py::init<int,double,Vec<3>,double>(), py::arg("order"), py::arg("kappa"),  py::arg("center"), py::arg("rad")=1.0)
    .def("AddCurrent", [](SphericalExpansionCF<Singular,Vec<3,Complex>>& self, Vec<3> sp, Vec<3> ep, Complex j, int num)
    { self.MP().AddCurrent(sp, ep, j, num); },
      py::arg("sp"), py::arg("ep"), py::arg("j"), py::arg("num")=100
      )
    ;

   // from potentialtools



  
  py::class_<SingularMLExpansionCF<Complex>, shared_ptr<SingularMLExpansionCF<Complex>>, CoefficientFunction> (m, "SingularMLExpansionCF")
    .def(py::init<Vec<3>, double, double>(), py::arg("center"), py::arg("r"), py::arg("kappa"))
    .def_property_readonly("expansion", [](SingularMLExpansionCF<Complex>& self) { return self.MLExpansion(); })
    .def("CreateRegularExpansion", &SingularMLExpansionCF<Complex>::CreateRegularExpansion, py::arg("center"), py::arg("r"))
    ;
  py::class_<RegularMLExpansionCF<Complex>, shared_ptr<RegularMLExpansionCF<Complex>>, CoefficientFunction> (m, "RegularMLExpansionCF")
    .def(py::init<shared_ptr<SingularMLExpansionCF<Complex>>,Vec<3>, double>(), py::arg("singularexpansion"), py::arg("center"), py::arg("r"))
    .def_property_readonly("expansion", [](RegularMLExpansionCF<Complex>& self) { return self.MLExpansion(); })
    ;


  py::class_<SingularMLExpansionCF<Vec<3,Complex>>, shared_ptr<SingularMLExpansionCF<Vec<3,Complex>>>, CoefficientFunction> (m, "BiotSavartSingularMLCF")
    .def(py::init<Vec<3>, double, double>(), py::arg("center"), py::arg("r"), py::arg("kappa"))
    .def("CreateRegularExpansion", &SingularMLExpansionCF<Vec<3,Complex>>::CreateRegularExpansion, py::arg("center"), py::arg("r"))    
    .def_property_readonly("expansion", [](SingularMLExpansionCF<Vec<3,Complex>>& self) { return self.MLExpansion(); })
    ;
  py::class_<RegularMLExpansionCF<Vec<3,Complex>>, shared_ptr<RegularMLExpansionCF<Vec<3,Complex>>>, CoefficientFunction> (m, "BiotSavartRegularMLCF")
    .def(py::init<shared_ptr<SingularMLExpansionCF<Vec<3,Complex>>>,Vec<3>, double>(), py::arg("mp"), py::arg("center"), py::arg("r"))
    .def_property_readonly("expansion", [](RegularMLExpansionCF<Vec<3,Complex>>& self) { return self.MLExpansion(); })
    ;

  

   

  // ************************** Potential and integral operators *******************************************

  py::class_<IntegralOperator,shared_ptr<IntegralOperator>> (m, "IntegralOperator")
    .def_property_readonly("mat", &IntegralOperator::GetMatrix)
    .def("NearFieldMatrix", &IntegralOperator::GetNearFieldMatrix)
    .def("GetPotential", &IntegralOperator::GetPotential,
         py::arg("gf"), py::arg("intorder")=nullopt, py::arg("nearfield_experimental")=false)

    .def("GetFMMInfo", [](shared_ptr<IntegralOperator> iop) {
      return FMMInfoToDict(iop->GetFMMInfo());
    })

    .def("CalcSubMatrix", [](shared_ptr<IntegralOperator> iop,
                             py::array_t<DofId> rowids, py::array_t<DofId> colids) {

      FlatArray<DofId> rowidsa(rowids.size(), rowids.mutable_data());
      FlatArray<DofId> colidsa(colids.size(), colids.mutable_data());      

      py::gil_scoped_release rel;
      LocalHeapMem<1000000> lh("CalcSubMatrix lh"); // TODO: form LocalHeapProvider in python_comp.cpp
      // iop->CalcSubMatrix(rowidsa, colidsa, lh);  // call it twice, for timing
      return iop->CalcSubMatrix(rowidsa, colidsa, lh);
    }, py::arg("rowids"), py::arg("colids"))

    .def("CalcSubMatrixCapsule", [](std::shared_ptr<IntegralOperator> iop)
    {
      using backend_callback_t = void(*)(int, int,
                                         const int*, const int*,
                                         double*, void*);

      struct Backend {
        backend_callback_t callback;
        void* ctx;
      };

      struct Context {
        std::shared_ptr<IntegralOperator> iop;
      };

      Context* ctx = new Context{iop};

      auto callback = [](int n, int m,
                         const int* rowid,
                         const int* colid,
                         double* data,
                         void* vctx)
      {
        Context* ctx = static_cast<Context*>(vctx);
        auto& iop = *ctx->iop;

        LocalHeapMem<1000000> lh("CalcSubMatrix lh");

        FlatArray<DofId> rowidsa(n, const_cast<int*>(rowid));
        FlatArray<DofId> colidsa(m, const_cast<int*>(colid));

        auto mat = iop.CalcSubMatrix(rowidsa, colidsa, lh);

        if(std::holds_alternative<Matrix<double>>(mat))
          {
            auto& matd = std::get<Matrix<double>>(mat);
            for (int i = 0; i < n; i++)
              for (int j = 0; j < m; j++)
                data[i + j*n] = matd(i, j);
          }
        else
          {
            auto& matc = std::get<Matrix<Complex>>(mat);
            for (int i = 0; i < n; i++)
              for (int j = 0; j < m; j++)
                {
                  data[2*(i + j*n)] = matc(i, j).real();
                  data[2*(i + j*n)+1] = matc(i, j).imag();
                }
          }
      };

      Backend* backend = new Backend{
        callback,
        ctx,
      };

      return py::capsule(backend, "backend", [](void* p)
      {
        auto* backend = static_cast<Backend*>(p);
        if (backend->ctx)
          delete static_cast<Context*>(backend->ctx);
        delete backend;
      });
    })
    
    .def("__add__", [](shared_ptr<IntegralOperator> a, shared_ptr<IntegralOperator> b)
    {
      return AddIntegralOperators(a, b);
    })
    .def("__sub__", [](shared_ptr<IntegralOperator> a, shared_ptr<IntegralOperator> b)
    {
      return AddIntegralOperators(a, ScaleIntegralOperator(b, -1.0));
    })
    .def("__rmul__", [](shared_ptr<IntegralOperator> a, double fac)
    {
      return ScaleIntegralOperator(a, fac);
    })
    .def("__rmul__", [](shared_ptr<IntegralOperator> a, Complex fac)
    {
      return ScaleIntegralOperator(a, fac);
    })
    .def("__neg__", [](shared_ptr<IntegralOperator> a)
    {
      return ScaleIntegralOperator(a, -1.0);
    })
    ;

  m.def("SingleLayerPotentialOperator", [](shared_ptr<FESpace> space, int intorder) -> shared_ptr<IntegralOperator>
  {
    WarnDeprecated("SingleLayerPotentialOperator", "LaplaceSL");
    return make_unique<GenericIntegralOperator<LaplaceSLKernel<3>>>(space, space, nullopt, nullopt,
                                                                    space->GetEvaluator(BND), space->GetEvaluator(BND),
                                                                    LaplaceSLKernel<3>(), intorder);
    
  }, py::arg("space"), py::arg("intorder")=3);

  m.def("SingleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space,
                                           optional<Region> trial_definedon, optional<Region> test_definedon,
                                           int intorder) -> shared_ptr<IntegralOperator>
  {
    WarnDeprecated("SingleLayerPotentialOperator", "LaplaceSL");
    return make_unique<GenericIntegralOperator<LaplaceSLKernel<3>>>(trial_space, test_space,
                                                                    trial_definedon, test_definedon,
                                                                    trial_space -> GetEvaluator(BND),
                                                                    test_space -> GetEvaluator(BND), LaplaceSLKernel<3>(), intorder);
    
  }, py::arg("trial_space"), py::arg("test_space"),
        py::arg("trial_definedon")=nullopt, py::arg("test_definedon")=nullopt,
        py::arg("intorder")=3);
  
  

  m.def("DoubleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space,
                                           optional<Region> trial_definedon, optional<Region> test_definedon,
                                           int intorder) -> shared_ptr<IntegralOperator>
  {
    WarnDeprecated("DoubleLayerPotentialOperator", "LaplaceDL");
    return make_unique<GenericIntegralOperator<LaplaceDLKernel<3>>>(trial_space, test_space,
                                                                    trial_definedon, test_definedon,
                                                                    trial_space -> GetEvaluator(BND),
                                                                    test_space -> GetEvaluator(BND),
                                                                    LaplaceDLKernel<3>(), intorder);    
  }, py::arg("trial_space"), py::arg("test_space"),
        py::arg("trial_definedon")=nullopt, py::arg("test_definedon")=nullopt,
        py::arg("intorder")=3);


  m.def("HypersingularOperator", [](shared_ptr<FESpace> space, optional<Region> definedon,
                                    int intorder) -> shared_ptr<IntegralOperator>
  {
    WarnDeprecated("HypersingularOperator", "LaplaceSL");
    return make_unique<GenericIntegralOperator<LaplaceSLKernel<3,3>>>(space, space, definedon, definedon,
                                                                    make_shared<T_DifferentialOperator<DiffOpBoundaryRot>>(),
                                                                    make_shared<T_DifferentialOperator<DiffOpBoundaryRot>>(),
                                                                    LaplaceSLKernel<3,3>(), intorder);
    
  }, py::arg("space"), py::arg("definedon")=nullopt,
        py::arg("intorder")=3);
  
  

  m.def("HelmholtzSingleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space, double kappa,
                                                    int intorder) -> shared_ptr<IntegralOperator>
  {
    WarnDeprecated("HelmholtzSingleLayerPotentialOperator", "HelmholtzSL");
    return make_unique<GenericIntegralOperator<HelmholtzSLKernel<3>>>(trial_space, test_space, nullopt, nullopt,
                                                                      trial_space -> GetEvaluator(BND),
                                                                      test_space -> GetEvaluator(BND),
                                                                      HelmholtzSLKernel<3>(kappa), intorder);
    
  }, py::arg("trial_space"), py::arg("test_space")=nullptr, py::arg("kappa"), py::arg("intorder")=3);



  m.def("HelmholtzDoubleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space, double kappa,
                                                    int intorder) -> shared_ptr<IntegralOperator>
  {
    WarnDeprecated("HelmholtzDoubleLayerPotentialOperator", "HelmholtzDL");
    return make_unique<GenericIntegralOperator<HelmholtzDLKernel<3>>>(trial_space, test_space, nullopt, nullopt,
                                                                      trial_space -> GetEvaluator(BND),
                                                                      test_space -> GetEvaluator(BND),
                                                                      HelmholtzDLKernel<3>(kappa), intorder);
    
  }, py::arg("trial_space"), py::arg("test_space")=nullptr, py::arg("kappa"), py::arg("intorder")=3);


  m.def("HelmholtzCombinedFieldOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space,
                                             optional<Region> trial_definedon, optional<Region> test_definedon,
                                             double kappa,
                                             int intorder) -> shared_ptr<IntegralOperator>
  {
    WarnDeprecated("HelmholtzCombinedFieldOperator", "HelmholtzCF");
    return make_unique<GenericIntegralOperator<CombinedFieldKernel<3>>>(trial_space, test_space, trial_definedon, test_definedon,
                                                                        trial_space -> GetEvaluator(BND),
                                                                        test_space -> GetEvaluator(BND),
                                                                        CombinedFieldKernel<3>(kappa), intorder);
    
  }, py::arg("trial_space"), py::arg("test_space")=nullptr,
        py::arg("trial_definedon")=nullopt, py::arg("test_definedon")=nullopt,        
        py::arg("kappa"), py::arg("intorder")=3);


    m.def("HelmholtzHypersingularOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space, double kappa,  int intorder) -> shared_ptr<IntegralOperator>
  {
    return make_unique<GenericIntegralOperator<HelmholtzHSKernel<3>>>(trial_space, test_space, nullopt, nullopt,
                                                                      make_shared<T_DifferentialOperator<DiffOpHelmholtz>>(),
                                                                      make_shared<T_DifferentialOperator<DiffOpHelmholtz>>(),
                                                                      HelmholtzHSKernel<3>(kappa), intorder);
    
  }, py::arg("trial_space"), py::arg("test_space")=nullptr, py::arg("kappa"), py::arg("intorder")=3);


  m.def("MaxwellSingleLayerPotentialOperator", [](shared_ptr<FESpace> space, double kappa, optional<Region> definedon,
                                                  int intorder) -> shared_ptr<IntegralOperator>
  {
    WarnDeprecated("MaxwellSingleLayerPotentialOperator", "HelmholtzSL");
    return make_unique<GenericIntegralOperator<MaxwellSLKernel<3>>>(space, space, definedon, definedon,
                                                                    make_shared<T_DifferentialOperator<DiffOpMaxwellNew>>(),
                                                                    make_shared<T_DifferentialOperator<DiffOpMaxwellNew>>(), 
                                                                    MaxwellSLKernel<3>(kappa), intorder);
    
  }, py::arg("space"), py::arg("kappa"), py::arg("definedon")=nullopt,
        py::arg("intorder")=3);

  
  m.def("MaxwellSingleLayerPotentialOperatorCurl", [](shared_ptr<FESpace> space, double kappa, optional<Region> definedon,
                                                      int intorder) -> shared_ptr<IntegralOperator>
  {
    WarnDeprecated("MaxwellSingleLayerPotentialOperatorCurl", "HelmholtzSL");
    return make_unique<GenericIntegralOperator<MaxwellSLKernel<3>>>(space, space, definedon, definedon,
                                                                    make_shared<T_DifferentialOperator<DiffOpMaxwell>>(),
                                                                    make_shared<T_DifferentialOperator<DiffOpMaxwell>>(), 
                                                                    MaxwellSLKernel<3>(kappa), intorder);
    
  }, py::arg("space"), py::arg("kappa"), py::arg("definedon")=nullopt,
        py::arg("intorder")=3);
  
  

  
  m.def("MaxwellDoubleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space,
                                                  double kappa, 
                                                  optional<Region> trial_definedon, optional<Region> test_definedon,
                                                  int intorder) -> shared_ptr<IntegralOperator>
  {
    WarnDeprecated("MaxwellDoubleLayerPotentialOperator", "MaxwellDL");
    return make_unique<GenericIntegralOperator<MaxwellDLKernel<3>>>(trial_space, test_space,
                                                                    trial_definedon, test_definedon,
                                                                    make_shared<T_DifferentialOperator<DiffOpRotatedTrace>>(),
                                                                    test_space->GetEvaluator(BND),
                                                                    MaxwellDLKernel<3>(kappa), intorder);
  }, py::arg("trial_space"), py::arg("test_space"), py::arg("kappa"),
        py::arg("trial_definedon")=nullopt, py::arg("test_definedon")=nullopt,        
        py::arg("intorder")=3);
  





  // ******************** Potential operators ************************************


  
  py::class_<BasePotentialCF, CoefficientFunction, shared_ptr<BasePotentialCF>> (m, "PotentialCF")
    .def("BuildLocalExpansion", [](shared_ptr<BasePotentialCF> potcf, const Region & region)
    {
      potcf->BuildLocalExpansion(region);
      return potcf;
    })
    ;

  py::class_<BasePotentialOperator, shared_ptr<BasePotentialOperator>> (m, "PotentialOperator")
    .def("Operator", [](shared_ptr<BasePotentialOperator> pot, string name) {
        return pot->MakeDiffBasePotential(name);
    })
    .def("__mul__", [](shared_ptr<BasePotentialOperator> pot, shared_ptr<CoefficientFunction> test_proxy) {
      return BasePotentialOperatorAndTest (pot, test_proxy); //  { pot, dynamic_pointer_cast<ProxyFunction>(test_proxy) };
    })
    .def("__rmul__", [](shared_ptr<BasePotentialOperator> pot, double fac) {
      return SumOfPotentialOperators(Scalar(fac), pot); 
    })
    .def("__rmul__", [](shared_ptr<BasePotentialOperator> pot, Complex fac) {
      return SumOfPotentialOperators(Scalar(fac), pot); 
    })
    
    .def("__call__", [](shared_ptr<BasePotentialOperator> pot, shared_ptr<GridFunction> gf) {
      return pot->MakePotentialCF(gf);
    })
    .def("__call__", [](shared_ptr<BasePotentialOperator> pot, shared_ptr<GridFunction> gf, const Region & region) {
      auto cf = pot->MakePotentialCF(gf);
      cf->BuildLocalExpansion(region);
      return cf;
    })
    .def("__str__", [](BasePotentialOperator & self) { return ToString(self); })
    ;
  
  py::class_<BasePotentialOperatorAndTest> (m, "BasePotentialOperatorAndTest")
    .def("__mul__", [](BasePotentialOperatorAndTest pottest, DifferentialSymbol dx)
    {
      return pottest.MakeIntegralOperator(dx); 
    })
    ;

  py::class_<SumOfPotentialOperatorsAndTest> (m, "SumOfPotentialOperatorsAndTest")
    .def("__mul__", [](SumOfPotentialOperatorsAndTest pottest, DifferentialSymbol dx)
    {
      return pottest.MakeIntegralOperator(dx);
    })
    ;


  py::class_<SumOfPotentialOperators> (m, "SumOfPotentialOperators")
    .def(py::self+py::self)
    .def(py::self-py::self)
    .def(-py::self)
    .def("__mul__", [](SumOfPotentialOperators sumpot, shared_ptr<CoefficientFunction> test_proxy)
    {
      return SumOfPotentialOperatorsAndTest(sumpot.Summands(), test_proxy);
    })
    .def("__call__", [](SumOfPotentialOperators sumpot, shared_ptr<GridFunction> gf)
    {
      return sumpot.MakePotentialCF(gf);
    })
    .def("__call__", [](SumOfPotentialOperators sumpot, shared_ptr<GridFunction> gf,
                        const Region & region)
    {
      return sumpot.MakePotentialCF(gf, region);
    })
    ;

  m.def("LaplaceSL", [&](shared_ptr<SumOfIntegrals> potential, py::kwargs kwargs) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    
    auto proxy = GetProxyWithFactor(igl->cf, true, igl->dx.vb);

    auto fes = proxy->GetFESpace();
    int fesorder = GetFESOrder (proxy);

    auto flags = CreateFlagsFromKwArgs(kwargs); 
    IntOp_Parameters ioparams(flags);
    // cout << ioflags << endl;
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    switch (proxy->Dimension())
      {
      case 1:
        if (fes->IsComplex())
          return make_shared<PotentialOperator<LaplaceSLKernel<3,1,Complex>>>
            (proxy, igl->dx.vb, definedon, proxy->Evaluator(),
             LaplaceSLKernel<3,1,Complex>{}, ioparams,
             fesorder+igl->dx.bonus_intorder);
        return make_shared<PotentialOperator<LaplaceSLKernel<3>>> (proxy, igl->dx.vb, definedon, proxy->Evaluator(),
                                                                   LaplaceSLKernel<3>{}, ioparams, 
                                                                   fesorder+igl->dx.bonus_intorder);
      case 3:
        if (fes->IsComplex())
          return make_shared<PotentialOperator<LaplaceSLKernel<3,3,Complex>>>
            (proxy, igl->dx.vb, definedon, proxy->Evaluator(),
             LaplaceSLKernel<3,3,Complex>{}, ioparams,
             fesorder+igl->dx.bonus_intorder);
        return make_shared<PotentialOperator<LaplaceSLKernel<3,3>>> (proxy, igl->dx.vb, definedon, proxy->Evaluator(),
                                                                     LaplaceSLKernel<3,3>{}, ioparams,
                                                                     fesorder+igl->dx.bonus_intorder);
      default:
        ;
      }
    throw Exception("only dim=1 and dim=3 LaplaceSL are supported");
  }, py::arg("potential"), docu_string(bem_operator_kwargs_doc));



  m.def("LaplaceDL", [](shared_ptr<SumOfIntegrals> potential, py::kwargs kwargs) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    if (igl->dx.vb != BND) throw Exception("need boundary integral");
    auto proxy = dynamic_pointer_cast<ProxyFunction>(igl->cf);
    auto fes = proxy->GetFESpace();
    /*
    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();
    while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
      {
        tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
        tmpeval = compeval->BaseDiffOp();
      }
    */
    int fesorder = GetFESOrder (proxy);
    auto flags = CreateFlagsFromKwArgs(kwargs);
    IntOp_Parameters ioparams(flags);
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));
    if (proxy->Dimension() == 1)
      {
        if (fes->IsComplex())
          return make_shared<PotentialOperator<LaplaceDLKernel<3,1,Complex>>> (proxy, igl->dx.vb, definedon, proxy->Evaluator(),
                                                                              LaplaceDLKernel<3,1,Complex>{}, ioparams, fesorder+igl->dx.bonus_intorder);
        return make_shared<PotentialOperator<LaplaceDLKernel<3>>> (proxy, igl->dx.vb, definedon, proxy->Evaluator(),
                                                                  LaplaceDLKernel<3>{}, ioparams, fesorder+igl->dx.bonus_intorder);
      }
    if (proxy->Dimension() == 3)
      {
        if (fes->IsComplex())
          return make_shared<PotentialOperator<LaplaceDLKernel<3,3,Complex>>> (proxy, igl->dx.vb, definedon, proxy->Evaluator(),
                                                                              LaplaceDLKernel<3,3,Complex>{}, ioparams, fesorder+igl->dx.bonus_intorder);
        return make_shared<PotentialOperator<LaplaceDLKernel<3,3>>> (proxy, igl->dx.vb, definedon, proxy->Evaluator(),
                                                                    LaplaceDLKernel<3,3>{}, ioparams, fesorder+igl->dx.bonus_intorder);
      }
    throw Exception("only dim=1 and dim=3 LaplaceDL are supported");
  }, py::arg("potential"), docu_string(bem_operator_kwargs_doc));

  m.def("HelmholtzSL", [](shared_ptr<SumOfIntegrals> potential, std::variant<double, Complex> kappa, py::kwargs kwargs) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    /*
    auto proxy = dynamic_pointer_cast<ProxyFunction>(igl->cf);
    auto fes = proxy->GetFESpace();
    */
    auto proxy = GetProxyWithFactor(igl->cf, true, igl->dx.vb);
    auto fes = proxy->GetFESpace();

    
    /*
    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();
    while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
      {
        tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
        tmpeval = compeval->BaseDiffOp();
      }
    */
    int fesorder = GetFESOrder (proxy);
    auto flags = CreateFlagsFromKwArgs(kwargs);
    IntOp_Parameters ioparams(flags);
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 3)
    {
      return std::visit( [&] (auto val) -> shared_ptr<BasePotentialOperator> {
        return make_shared<PotentialOperator<HelmholtzSLKernel<3,3,decltype(val)>>>
            (proxy, igl->dx.vb, definedon, proxy->Evaluator(), val, ioparams, fesorder+igl->dx.bonus_intorder);
        }, kappa); 
    }
    else if (proxy->Dimension() == 1)
      {
      return std::visit( [&] (auto val) -> shared_ptr<BasePotentialOperator> {
        return make_shared<PotentialOperator<HelmholtzSLKernel<3,1,decltype(val)>>>
            (proxy, igl->dx.vb, definedon, proxy->Evaluator(), val, ioparams, fesorder+igl->dx.bonus_intorder);
        }, kappa); 
    }
    else
      throw Exception("only dim=1 and dim=3 HelmholtzSL are supported");
  }, py::arg("potential"), py::arg("kappa"), docu_string(bem_operator_kwargs_doc));

  m.def("HelmholtzDL", [](shared_ptr<SumOfIntegrals> potential, std::variant<double, Complex> kappa, py::kwargs kwargs) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    if (igl->dx.vb != BND) throw Exception("need boundary integral");

    
    // auto proxy = dynamic_pointer_cast<ProxyFunction>(igl->cf);
    // auto fes = proxy->GetFESpace();
    auto proxy = GetProxyWithFactor(igl->cf, true);
    auto fes = proxy->GetFESpace();

    /*
    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();
    while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
      {
        tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
        tmpeval = compeval->BaseDiffOp();
      }
    */
    int fesorder = GetFESOrder (proxy);    
    auto flags = CreateFlagsFromKwArgs(kwargs);
    IntOp_Parameters ioparams(flags);
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 3)
    {
      return std::visit( [&] (auto val) -> shared_ptr<BasePotentialOperator> {
        return make_shared<PotentialOperator<HelmholtzDLKernel<3,3,decltype(val)>>>
            (proxy, igl->dx.vb, definedon, proxy->Evaluator(), val, ioparams, fesorder+igl->dx.bonus_intorder);
        }, kappa); 
    }
    if (proxy->Dimension() == 1)
    {
      return std::visit( [&] (auto val) -> shared_ptr<BasePotentialOperator> {
        return make_shared<PotentialOperator<HelmholtzDLKernel<3,1,decltype(val)>>>
            (proxy, igl->dx.vb, definedon, proxy->Evaluator(), val, ioparams, fesorder+igl->dx.bonus_intorder);
        }, kappa); 
    }
    else
      throw Exception("only dim=1 and dim=3 HelmholtzDL are supported");
  }, py::arg("potential"), py::arg("kappa"), docu_string(bem_operator_kwargs_doc));



  m.def("HelmholtzCF", [](shared_ptr<SumOfIntegrals> potential, std::variant<double, Complex> kappa, py::kwargs kwargs) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    if (igl->dx.vb != BND) throw Exception("need boundary integral");
    auto proxy = dynamic_pointer_cast<ProxyFunction>(igl->cf);
    auto fes = proxy->GetFESpace();

    /*
    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();
    while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
      {
        tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
        tmpeval = compeval->BaseDiffOp();
      }
    */
    int fesorder = GetFESOrder (proxy);

    auto flags = CreateFlagsFromKwArgs(kwargs); 
    IntOp_Parameters ioparams(flags);
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 1)
      return std::visit( [&] (auto val) -> shared_ptr<BasePotentialOperator> {
        return make_shared<PotentialOperator<CombinedFieldKernel<3,1,decltype(val)>>>
            (proxy, igl->dx.vb, definedon, proxy->Evaluator(), val, ioparams, fesorder+igl->dx.bonus_intorder);
        }, kappa);
    if (proxy->Dimension() == 3)
      return std::visit( [&] (auto val) -> shared_ptr<BasePotentialOperator> {
        return make_shared<PotentialOperator<CombinedFieldKernel<3,3,decltype(val)>>>
            (proxy, igl->dx.vb, definedon, proxy->Evaluator(), val, ioparams, fesorder+igl->dx.bonus_intorder);
        }, kappa); 
    throw Exception("only dim=1 HelmholtzCF is supported");
  }, py::arg("potential"), py::arg("kappa"), docu_string(bem_operator_kwargs_doc));

  m.def("MaxwellDL", [](shared_ptr<SumOfIntegrals> potential, std::variant<double, Complex> kappa, py::kwargs kwargs) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    // if (igl->dx.vb != BND) throw Exception("need boundary integral");

    auto proxy = GetProxyWithFactor(igl->cf, true);
    auto fes = proxy->GetFESpace();

    int fesorder = GetFESOrder (proxy);    
    auto flags = CreateFlagsFromKwArgs(kwargs);
    IntOp_Parameters ioparams(flags);
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 3)
    {
      return std::visit( [&] (auto val) -> shared_ptr<BasePotentialOperator> {
        return make_shared<PotentialOperator<MaxwellDLKernel<3,decltype(val)>>>
            (proxy, igl->dx.vb, definedon, proxy->Evaluator(), val, ioparams, fesorder+igl->dx.bonus_intorder);
        }, kappa); 
    }
    else
      throw Exception("only dim=3 MaxwellDL are supported");
  }, py::arg("potential"), py::arg("kappa"), docu_string(bem_operator_kwargs_doc));


  m.def("LameSL", [](shared_ptr<SumOfIntegrals> potential, double E, double nu, py::kwargs kwargs) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    if (igl->dx.vb != BND) throw Exception("need boundary integral");
    
    // auto [proxy,factor] = GetProxyAndFactor(igl->cf, true);
    auto proxy = GetProxyWithFactor(igl->cf, true);
    
    auto fes = proxy->GetFESpace();
    auto flags = CreateFlagsFromKwArgs(kwargs);
    IntOp_Parameters ioparams(flags);

    /*
    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();
    while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
      {
        tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
        tmpeval = compeval->BaseDiffOp();
      }
    */
    int fesorder = GetFESOrder (proxy);
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 3)
      return make_shared<PotentialOperator<LameSLKernel<3>>> (proxy, igl->dx.vb, definedon, proxy->Evaluator(),
                                                              LameSLKernel<3>{E,nu}, ioparams, fesorder /* tmpfes->GetOrder()*/ +igl->dx.bonus_intorder);

    throw Exception("only dim=3 LameSL is supported");
    }, py::arg("term"), py::arg("E"), py::arg("nu"), docu_string(bem_operator_kwargs_doc));



  m.def("GetDofCoordinates", [](shared_ptr<FESpace> fes)
  {
    auto ma = fes->GetMeshAccess();
    Matrix pnts(fes->GetNDof(), 3);
    
    pnts.Col(0) = 1;
    pnts.Col(1) = 2;
    pnts.Col(2) = 3;

    Array<DofId> dofs;
    for (auto node : ma->Nodes(NT_VERTEX))
      {
        fes->GetDofNrs(node, dofs);
        Vec<3> p = ma->GetPoint<3>(node.GetNr());
        for (auto d : dofs)
          pnts.Row(d) = p;
      }

    for (auto node : ma->Nodes(NT_FACE))
      {
        fes->GetDofNrs(node, dofs);
        auto vertices = ma->GetFacePNums(node.GetNr());
        Vec<3> p(0,0,0);
        for (auto v : vertices)
          p += ma->GetPoint<3>(v);
        p /= vertices.Size();
        for (auto d : dofs)
          pnts.Row(d) = p;
      }

    for (auto node : ma->Nodes(NT_EDGE))
      {
        fes->GetDofNrs(node, dofs);
        auto vertices = ma->GetEdgePNums(node.GetNr());
        Vec<3> p(0,0,0);
        for (auto v : vertices)
          p += ma->GetPoint<3>(v);
        p /= vertices.Size();
        for (auto d : dofs)
          pnts.Row(d) = p;
      }
    
    return pnts;
  });
}

#endif // NGS_PYTHON
