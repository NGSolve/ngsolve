#ifdef NGS_PYTHON
#include <regex>

#include "../ngstd/python_ngstd.hpp"
#include "python_comp.hpp"

#include "mptools.hpp"
#include "potentialtools.hpp"
#include "ngbem.hpp"
#include "mp_coefficient.hpp"


using namespace ngsbem;



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

  py::class_<SingularMLMultiPole<Complex>, shared_ptr<SingularMLMultiPole<Complex>>> (m, "SingularMLMP")
    .def(py::init<Vec<3>,double,double>())
    .def("AddCharge", &SingularMLMultiPole<Complex>::AddCharge)
    .def("AddDipole", &SingularMLMultiPole<Complex>::AddDipole)
    .def("AddChargeDensity", [](SingularMLMultiPole<Complex> & mp, shared_ptr<CoefficientFunction> charge,
                                ngcomp::Region reg) { AddChargeDensity(mp,charge,reg); })
    
    .def("Calc", &SingularMLMultiPole<Complex>::CalcMP)
    .def("Norm", &SingularMLMultiPole<Complex>::Norm)    
    .def("__str__", [](SingularMLMultiPole<Complex>& mlmp) { return ToString<>(mlmp); })
    .def("Print", [](const SingularMLMultiPole<Complex> &self) {
        py::scoped_ostream_redirect stream(
            std::cout,
            py::module_::import("sys").attr("stdout")
        );
        self.Print(std::cout);
    })
    ;

  py::class_<SingularMLMultiPole<Vec<3,Complex>>, shared_ptr<SingularMLMultiPole<Vec<3,Complex>>>> (m, "SingularMLMP3")
    .def(py::init<Vec<3>,double,double>())
    .def("AddCurrent", &SingularMLMultiPole<Vec<3,Complex>>::AddCurrent, py::arg("sp"), py::arg("ep"), py::arg("j"), py::arg("num")=100)
    .def("AddCurrentDensity", [](SingularMLMultiPole<Vec<3,Complex>> & mp, shared_ptr<CoefficientFunction> current,
                                 ngcomp::Region reg) { AddCurrentDensity(mp,current,reg); })
    
    .def("Calc", &SingularMLMultiPole<Vec<3,Complex>>::CalcMP)
    // .def("Norm", &SingularMLMultiPole<Complex>::Norm)    
    // .def("__str__", [](SingularMLMultiPole<Complex>& mlmp) { return ToString<>(mlmp); })
    ;

  
  py::class_<RegularMLMultiPole<Complex>, shared_ptr<RegularMLMultiPole<Complex>>> (m, "RegularMLMP")
    .def("Norm", &RegularMLMultiPole<Complex>::Norm)    
    .def("__str__", [](RegularMLMultiPole<Complex>& mlmp) { return ToString<>(mlmp); })
    .def("Print", [](const RegularMLMultiPole<Complex> &self) {
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

  py::class_<MultiPoleCF<MPRegular>, shared_ptr<MultiPoleCF<MPRegular>>, CoefficientFunction> (m, "RegularMultiPoleCF")
    .def(py::init<int,double,Vec<3>,double>(), py::arg("order"), py::arg("kappa"),  py::arg("center"), py::arg("rad")=1.0)    
    .def_property_readonly("sh", [](MultiPoleCF<MPRegular>& self) -> SphericalHarmonics<Complex>& { return self.SH(); })
    .def("ShiftZ", [](MultiPoleCF<MPRegular>& self, double z, MultiPoleCF<MPRegular> & target) { self.ShiftZ(z, target.MP()); })
    .def("Transform", [](MultiPoleCF<MPRegular>& self, MultiPoleCF<MPRegular> & target) { self.Transform(target); })
    .def("TransformAdd", [](MultiPoleCF<MPRegular>& self, MultiPoleCF<MPRegular> & target) { self.TransformAdd(target); })
    .def("Spectrum", [](MultiPoleCF<MPRegular>& self, bool scaled) { return self.MP().Spectrum(scaled); }, py::arg("scaled"))    
    ;

  py::class_<MultiPoleCF<MPSingular>, shared_ptr<MultiPoleCF<MPSingular>>, CoefficientFunction> (m, "SingularMultiPoleCF")
    .def(py::init<int,double,Vec<3>,double>(), py::arg("order"), py::arg("kappa"),  py::arg("center"), py::arg("rad")=1.0)
    .def_property_readonly("sh", [](MultiPoleCF<MPSingular>& self) -> SphericalHarmonics<Complex>& { return self.SH(); })
    .def("AddCharge", [](MultiPoleCF<MPSingular>& self, Vec<3> x, Complex c) { self.MP().AddCharge(x, c); })
    .def("AddDipole", [](MultiPoleCF<MPSingular>& self, Vec<3> x, Vec<3> d, Complex c) { self.MP().AddDipole(x, d, c); })
    .def("ShiftZ", [](MultiPoleCF<MPSingular>& self, double z, MultiPoleCF<MPRegular> & target) { self.ShiftZ(z, target.MP()); })
    .def("ShiftZ", [](MultiPoleCF<MPSingular>& self, double z, MultiPoleCF<MPSingular> & target) { self.ShiftZ(z, target.MP()); })        
    .def("Transform", [](MultiPoleCF<MPSingular>& self, MultiPoleCF<MPRegular> & target) { self.Transform(target); })
    .def("Transform", [](MultiPoleCF<MPSingular>& self, MultiPoleCF<MPSingular> & target) { self.Transform(target); })
    .def("TransformAdd", [](MultiPoleCF<MPSingular>& self, MultiPoleCF<MPRegular> & target) { self.TransformAdd(target); })
    .def("TransformAdd", [](MultiPoleCF<MPSingular>& self, MultiPoleCF<MPSingular> & target) { self.TransformAdd(target); })
    .def("Spectrum", [](MultiPoleCF<MPSingular>& self, bool scaled) { return self.MP().Spectrum(scaled); }, py::arg("scaled"))
    ;

  py::class_<MultiPoleCF<MPSingular,Vec<3,Complex>>, shared_ptr<MultiPoleCF<MPSingular,Vec<3,Complex>>>, CoefficientFunction> (m, "BiotSavartCF")
    .def(py::init<int,double,Vec<3>,double>(), py::arg("order"), py::arg("kappa"),  py::arg("center"), py::arg("rad")=1.0)
    .def("AddCurrent", [](MultiPoleCF<MPSingular,Vec<3,Complex>>& self, Vec<3> sp, Vec<3> ep, Complex j, int num)
    { self.MP().AddCurrent(sp, ep, j, num); },
      py::arg("sp"), py::arg("ep"), py::arg("j"), py::arg("num")=100
      )
    ;

   // from potentialtools



  
  py::class_<SingularMLMultiPoleCF<Complex>, shared_ptr<SingularMLMultiPoleCF<Complex>>, CoefficientFunction> (m, "SingularMLMultiPoleCF")
    .def(py::init<Vec<3>, double, double>(), py::arg("center"), py::arg("r"), py::arg("kappa"))
    .def_property_readonly("mlmp", [](SingularMLMultiPoleCF<Complex>& self) { return self.MLMP(); })
    .def("CreateRegularExpansion", &SingularMLMultiPoleCF<Complex>::CreateRegularExpansion, py::arg("center"), py::arg("r"))
    ;
  py::class_<RegularMLMultiPoleCF<Complex>, shared_ptr<RegularMLMultiPoleCF<Complex>>, CoefficientFunction> (m, "RegularMLMultiPoleCF")
    .def(py::init<shared_ptr<SingularMLMultiPoleCF<Complex>>,Vec<3>, double>(), py::arg("mp"), py::arg("center"), py::arg("r"))
    .def_property_readonly("mlmp", [](RegularMLMultiPoleCF<Complex>& self) { return self.MLMP(); })
    ;


  py::class_<SingularMLMultiPoleCF<Vec<3,Complex>>, shared_ptr<SingularMLMultiPoleCF<Vec<3,Complex>>>, CoefficientFunction> (m, "BiotSavartSingularMLCF")
    .def(py::init<Vec<3>, double, double>(), py::arg("center"), py::arg("r"), py::arg("kappa"))
    .def("CreateRegularExpansion", &SingularMLMultiPoleCF<Vec<3,Complex>>::CreateRegularExpansion, py::arg("center"), py::arg("r"))    
    .def_property_readonly("mlmp", [](SingularMLMultiPoleCF<Vec<3,Complex>>& self) { return self.MLMP(); })
    ;
  py::class_<RegularMLMultiPoleCF<Vec<3,Complex>>, shared_ptr<RegularMLMultiPoleCF<Vec<3,Complex>>>, CoefficientFunction> (m, "BiotSavartRegularMLCF")
    .def(py::init<shared_ptr<SingularMLMultiPoleCF<Vec<3,Complex>>>,Vec<3>, double>(), py::arg("mp"), py::arg("center"), py::arg("r"))
    .def_property_readonly("mlmp", [](RegularMLMultiPoleCF<Vec<3,Complex>>& self) { return self.MLMP(); })
    ;

  

   

  // ************************** Potential and integral operators *******************************************



  py::class_<IntegralOperator,shared_ptr<IntegralOperator>> (m, "IntegralOperator")
    .def_property_readonly("mat", &IntegralOperator::GetMatrix)
    .def("GetPotential", &IntegralOperator::GetPotential,
         py::arg("gf"), py::arg("intorder")=nullopt, py::arg("nearfield_experimental")=false)
    ;
  
  m.def("SingleLayerPotentialOperator", [](shared_ptr<FESpace> space, int intorder) -> shared_ptr<IntegralOperator>
  {
    return make_unique<GenericIntegralOperator<LaplaceSLKernel<3>>>(space, space, nullopt, nullopt,
                                                                    space->GetEvaluator(BND), space->GetEvaluator(BND),
                                                                    LaplaceSLKernel<3>(), intorder);
    
  }, py::arg("space"), py::arg("intorder")=3);

  m.def("SingleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space,
                                           optional<Region> trial_definedon, optional<Region> test_definedon,
                                           int intorder) -> shared_ptr<IntegralOperator>
  {
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
    return make_unique<GenericIntegralOperator<LaplaceHSKernel<3>>>(space, space, definedon, definedon,
                                                                    make_shared<T_DifferentialOperator<DiffOpBoundaryRot>>(),
                                                                    make_shared<T_DifferentialOperator<DiffOpBoundaryRot>>(), 
                                                                    LaplaceHSKernel<3>(), intorder);
    
  }, py::arg("space"), py::arg("definedon")=nullopt,
        py::arg("intorder")=3);
  
  

  m.def("HelmholtzSingleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space, double kappa,
                                                    int intorder) -> shared_ptr<IntegralOperator>
  {
    return make_unique<GenericIntegralOperator<HelmholtzSLKernel<3>>>(trial_space, test_space, nullopt, nullopt,
                                                                      trial_space -> GetEvaluator(BND),
                                                                      test_space -> GetEvaluator(BND),
                                                                      HelmholtzSLKernel<3>(kappa), intorder);
    
  }, py::arg("trial_space"), py::arg("test_space")=nullptr, py::arg("kappa"), py::arg("intorder")=3);



  m.def("HelmholtzDoubleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space, double kappa,
                                                    int intorder) -> shared_ptr<IntegralOperator>
  {
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
    return make_unique<GenericIntegralOperator<MaxwellSLKernel<3>>>(space, space, definedon, definedon,
                                                                    make_shared<T_DifferentialOperator<DiffOpMaxwellNew>>(),
                                                                    make_shared<T_DifferentialOperator<DiffOpMaxwellNew>>(), 
                                                                    MaxwellSLKernel<3>(kappa), intorder);
    
  }, py::arg("space"), py::arg("kappa"), py::arg("definedon")=nullopt,
        py::arg("intorder")=3);

  
  m.def("MaxwellSingleLayerPotentialOperatorCurl", [](shared_ptr<FESpace> space, double kappa, optional<Region> definedon,
                                                      int intorder) -> shared_ptr<IntegralOperator>
  {
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
    .def("BuildLocalExpansion", [](BasePotentialCF & potcf, const Region & region) { potcf.BuildLocalExpansion(region); })
    ;

  py::class_<BasePotentialOperator, shared_ptr<BasePotentialOperator>> (m, "PotentialOperator")
    .def("__mul__", [](shared_ptr<BasePotentialOperator> pot, shared_ptr<CoefficientFunction> test_proxy) {
      return BasePotentialOperatorAndTest (pot, test_proxy); //  { pot, dynamic_pointer_cast<ProxyFunction>(test_proxy) };
    })
    .def("__call__", [](shared_ptr<BasePotentialOperator> pot, shared_ptr<GridFunction> gf) {
      return pot->MakePotentialCF(gf);
    })
    ;
  
  py::class_<BasePotentialOperatorAndTest> (m, "BasePotentialOperatorAndTest")
    .def("__mul__", [](BasePotentialOperatorAndTest pottest, DifferentialSymbol dx)
    {
      return pottest.MakeIntegralOperator(dx); 
    })
    ;
  
  m.def("LaplaceSL", [](shared_ptr<SumOfIntegrals> potential) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    if (igl->dx.vb != BND) throw Exception("need boundary integral");
    
    auto [proxy,factor] = GetProxyAndFactor(igl->cf, true);
    
    auto fes = proxy->GetFESpace();
    
    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();
    while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
      {
        tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
        tmpeval = compeval->BaseDiffOp();
      }
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 1)
      return make_shared<PotentialOperator<LaplaceSLKernel<3>>> (proxy, factor, definedon, proxy->Evaluator(),
                                                                 LaplaceSLKernel<3>{}, tmpfes->GetOrder()+igl->dx.bonus_intorder);
    if (proxy->Dimension() == 3)
      return make_shared<PotentialOperator<LaplaceHSKernel<3>>> (proxy, factor, definedon, proxy->Evaluator(),
                                                                 LaplaceHSKernel<3>{}, tmpfes->GetOrder()+igl->dx.bonus_intorder);

    throw Exception("only dim=1 and dim=3 LaplaceSL are supported");
  });



  m.def("LaplaceDL", [](shared_ptr<SumOfIntegrals> potential) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    if (igl->dx.vb != BND) throw Exception("need boundary integral");
    auto proxy = dynamic_pointer_cast<ProxyFunction>(igl->cf);
    auto fes = proxy->GetFESpace();

    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();
    while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
      {
        tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
        tmpeval = compeval->BaseDiffOp();
      }
        
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));
    return make_shared<PotentialOperator<LaplaceDLKernel<3>>> (proxy, nullptr, definedon, proxy->Evaluator(),
                                                               LaplaceDLKernel<3>{}, tmpfes->GetOrder()+igl->dx.bonus_intorder);
  });

  m.def("HelmholtzSL", [](shared_ptr<SumOfIntegrals> potential, double kappa) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    if (igl->dx.vb != BND) throw Exception("need boundary integral");
    auto proxy = dynamic_pointer_cast<ProxyFunction>(igl->cf);
    auto fes = proxy->GetFESpace();

    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();
    while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
      {
        tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
        tmpeval = compeval->BaseDiffOp();
      }
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 3)
    {
      return make_shared<PotentialOperator<HelmholtzSLVecKernel<3>>> (proxy, nullptr, definedon, proxy->Evaluator(),
                                                                      HelmholtzSLVecKernel<3>(kappa), tmpfes->GetOrder()+igl->dx.bonus_intorder);
    }
    else if (proxy->Dimension() == 1)
      {
      HelmholtzSLKernel<3> kernel(kappa);
      return make_shared<PotentialOperator<HelmholtzSLKernel<3>>> (proxy, nullptr, definedon, proxy->Evaluator(),
                                                                   kernel, tmpfes->GetOrder()+igl->dx.bonus_intorder);
    }
    else
      throw Exception("only dim=1 and dim=3 HelmholtzSL are supported");
  });

  m.def("HelmholtzDL", [](shared_ptr<SumOfIntegrals> potential, double kappa) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    if (igl->dx.vb != BND) throw Exception("need boundary integral");
    auto proxy = dynamic_pointer_cast<ProxyFunction>(igl->cf);
    auto fes = proxy->GetFESpace();

    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();
    while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
      {
        tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
        tmpeval = compeval->BaseDiffOp();
      }
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 3)
    {
      MaxwellDLKernel<3> kernel(kappa);
      return make_shared<PotentialOperator<MaxwellDLKernel<3>>> (proxy, nullptr, definedon, proxy->Evaluator(),
                                                                 kernel, tmpfes->GetOrder()+igl->dx.bonus_intorder);
    }
    else if (proxy->Dimension() == 1)
    {
      HelmholtzDLKernel<3> kernel(kappa);
      return make_shared<PotentialOperator<HelmholtzDLKernel<3>>> (proxy, nullptr, definedon, proxy->Evaluator(),
                                                                   kernel, tmpfes->GetOrder()+igl->dx.bonus_intorder);
    }
    else
      throw Exception("only dim=1 and dim=3 HelmholtzDL are supported");
  });



  m.def("HelmholtzCF", [](shared_ptr<SumOfIntegrals> potential, double kappa) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    if (igl->dx.vb != BND) throw Exception("need boundary integral");
    auto proxy = dynamic_pointer_cast<ProxyFunction>(igl->cf);
    auto fes = proxy->GetFESpace();

    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();
    while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
      {
        tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
        tmpeval = compeval->BaseDiffOp();
      }
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 1)
      return make_shared<PotentialOperator<CombinedFieldKernel<3>>> (proxy, nullptr, definedon, proxy->Evaluator(),
                                                                     CombinedFieldKernel<3>(kappa), tmpfes->GetOrder()+igl->dx.bonus_intorder);
    throw Exception("only dim=1 and dim=3 HelmholtzDL are supported");
  });

  
}

#endif // NGS_PYTHON
