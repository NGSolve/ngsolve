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
    return make_unique<GenericIntegralOperator<LaplaceSLKernel<3,3>>>(space, space, definedon, definedon,
                                                                    make_shared<T_DifferentialOperator<DiffOpBoundaryRot>>(),
                                                                    make_shared<T_DifferentialOperator<DiffOpBoundaryRot>>(), 
                                                                    LaplaceSLKernel<3,3>(), intorder);
    
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
    .def("BuildLocalExpansion", [](shared_ptr<BasePotentialCF> potcf, const Region & region)
    {
      potcf->BuildLocalExpansion(region);
      return potcf;
    })
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
    
    auto proxy = GetProxyWithFactor(igl->cf, true);

    auto fes = proxy->GetFESpace();
    int fesorder = GetFESOrder (proxy);
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    switch (proxy->Dimension())
      {
      case 1:
        return make_shared<PotentialOperator<LaplaceSLKernel<3>>> (proxy, definedon, proxy->Evaluator(),
                                                                   LaplaceSLKernel<3>{}, fesorder+igl->dx.bonus_intorder);
      case 3:
        return make_shared<PotentialOperator<LaplaceSLKernel<3,3>>> (proxy, definedon, proxy->Evaluator(),
                                                                     LaplaceSLKernel<3,3>{}, fesorder+igl->dx.bonus_intorder);
      default:
        ;
      }
    throw Exception("only dim=1 and dim=3 LaplaceSL are supported");
  });



  m.def("LaplaceDL", [](shared_ptr<SumOfIntegrals> potential) -> shared_ptr<BasePotentialOperator> {
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
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));
    if (proxy->Dimension() == 1)
      return make_shared<PotentialOperator<LaplaceDLKernel<3>>> (proxy, definedon, proxy->Evaluator(),
                                                                LaplaceDLKernel<3>{}, fesorder+igl->dx.bonus_intorder);
    if (proxy->Dimension() == 3)
      return make_shared<PotentialOperator<LaplaceDLKernel<3,3>>> (proxy, definedon, proxy->Evaluator(),
                                                                LaplaceDLKernel<3,3>{}, fesorder+igl->dx.bonus_intorder);
    throw Exception("only dim=1 and dim=3 LaplaceDL are supported");
  });

  m.def("HelmholtzSL", [](shared_ptr<SumOfIntegrals> potential, double kappa) -> shared_ptr<BasePotentialOperator> {
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
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 3)
    {
      return make_shared<PotentialOperator<HelmholtzSLVecKernel<3>>> (proxy, definedon, proxy->Evaluator(),
                                                                      HelmholtzSLVecKernel<3>(kappa), fesorder+igl->dx.bonus_intorder);
    }
    else if (proxy->Dimension() == 1)
      {
      HelmholtzSLKernel<3> kernel(kappa);
      return make_shared<PotentialOperator<HelmholtzSLKernel<3>>> (proxy, definedon, proxy->Evaluator(),
                                                                   kernel, fesorder+igl->dx.bonus_intorder);
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
    {
      MaxwellDLKernel<3> kernel(kappa);
      return make_shared<PotentialOperator<MaxwellDLKernel<3>>> (proxy, definedon, proxy->Evaluator(),
                                                                 kernel, fesorder+igl->dx.bonus_intorder);
    }
    else if (proxy->Dimension() == 1)
    {
      HelmholtzDLKernel<3> kernel(kappa);
      return make_shared<PotentialOperator<HelmholtzDLKernel<3>>> (proxy, definedon, proxy->Evaluator(),
                                                                   kernel, fesorder+igl->dx.bonus_intorder);
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

    if (proxy->Dimension() == 1)
      return make_shared<PotentialOperator<CombinedFieldKernel<3>>> (proxy, definedon, proxy->Evaluator(),
                                                                     CombinedFieldKernel<3>(kappa), fesorder+igl->dx.bonus_intorder);
    throw Exception("only dim=1 HelmholtzCF is supported");
  });



  m.def("LameSL", [](shared_ptr<SumOfIntegrals> potential, double E, double nu) -> shared_ptr<BasePotentialOperator> {
    if (potential->icfs.Size()!=1) throw Exception("need one integral");
    auto igl = potential->icfs[0];
    if (igl->dx.vb != BND) throw Exception("need boundary integral");
    
    // auto [proxy,factor] = GetProxyAndFactor(igl->cf, true);
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
    
    optional<Region> definedon;
    if (igl->dx.definedon)
      definedon = Region(fes->GetMeshAccess(), igl->dx.vb, get<1> (*(igl->dx.definedon)));

    if (proxy->Dimension() == 3)
      return make_shared<PotentialOperator<LameSLKernel<3>>> (proxy, definedon, proxy->Evaluator(),
                                                              LameSLKernel<3>{E,nu}, fesorder /* tmpfes->GetOrder()*/ +igl->dx.bonus_intorder);

    throw Exception("only dim=3 LameSL is supported");
  });



  
}

#endif // NGS_PYTHON
