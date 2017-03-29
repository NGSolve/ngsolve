
#include "catch.hpp"
#include <fem.hpp>

using namespace ngfem;

#define TEST_OPERATOR_COEFFICIENTFUNCTION(CF)   \
  SECTION(#CF)                                  \
  {                                             \
    TestCoefficientFunction(CF);                \
  }


void TestCoefficientFunction(shared_ptr<CoefficientFunction> cf)
{
  LocalHeap lh(100000, "lh");
  IntegrationRule ir(ET_TET, 3);
  FE_ElementTransformation<3,3> trafo(ET_TET);
  MappedIntegrationRule<3,3> mir(ir, trafo, lh);

  SIMD_IntegrationRule simd_ir(ET_TET, 3);
  SIMD_MappedIntegrationRule<3,3> simd_mir(simd_ir, trafo, lh);
  shared_ptr<CoefficientFunction> c_cf_f = Compile(cf,false);
  shared_ptr<CoefficientFunction> c_cf_t = Compile(cf,false);

  // we compute directional derivative by the first proxy ...
  Array<ProxyFunction*> proxies;
  cf->TraverseTree
    ([&] (CoefficientFunction & nodecf)
     {
       auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
       if (proxy && !proxies.Contains(proxy))
         proxies.Append (proxy);
     });

  ProxyUserData ud(proxies.Size(), lh);
  trafo.userdata = &ud;
  ScalarDummyFE<ET_TET> dummy_fe;
  Vector<> dummy_elvec;
  ud.fel = &dummy_fe;
  ud.elx = &dummy_elvec;
  ud.lh = &lh;
  // we give some value to the proxy functions
  for (ProxyFunction * proxy : proxies)
    {
      ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);
      ud.GetAMemory(proxy) = 100;
      ud.GetMemory(proxy) = 100;
    }
  ud.trialfunction = nullptr;
  if (proxies.Size())
    ud.trialfunction = proxies[0];
  ud.trial_comp = 0;
  ud.testfunction = nullptr;
  ud.test_comp = 0;


  Matrix<> vals(ir.Size(),cf->Dimension());
  cf->Evaluate(mir,vals);
  AMatrix<> simd_vals (cf->Dimension(), simd_ir.GetNIP());
  Matrix<> c_vals(ir.Size(), cf->Dimension());

  SECTION ("Evaluate")
    {
      SECTION("SIMD")
        {
          simd_vals = 0;
          cf->Evaluate(simd_mir, simd_vals);
          CHECK(L2Norm(vals-Trans(simd_vals)) < 1e-12);
        }

      SECTION("Compile")
        {
          c_vals = 0;
          c_cf_f->Evaluate(mir,c_vals);
          CHECK(L2Norm(vals-c_vals) < 1e-12);
        }

      SECTION("TrueCompile")
        {
          c_vals = 0;
          c_cf_t->Evaluate(mir,c_vals);
          CHECK(L2Norm(vals-c_vals) < 1e-12);
        }

      try
        {
          SECTION("SIMD")
            {
              simd_vals = 0;
              cf->Evaluate(simd_mir, simd_vals);
              CHECK(L2Norm(vals-Trans(simd_vals)) < 1e-12);

              SECTION("Compile")
                {
                  simd_vals = 0;
                  c_cf_f->Evaluate(simd_mir, simd_vals);
                  CHECK(L2Norm(vals-Trans(simd_vals)) < 1e-12);
                }

              SECTION("TrueCompile")
                {
                  simd_vals = 0;
                  c_cf_t->Evaluate(simd_mir, simd_vals);
                  CHECK(L2Norm(vals-Trans(simd_vals)) < 1e-12);
                }
            }
        }
      catch (Exception ex)
        {
          CHECK_THROWS_AS(throw ex, ExceptionNOSIMD);
          WARN("SIMD not implemented");
        }
    }

  SECTION("Derivative")
    {
      Matrix<> derivs(ir.Size(), cf->Dimension());
      cf->EvaluateDeriv(mir,vals,derivs);
      AMatrix<> simd_derivs(cf->Dimension(),simd_ir.GetNIP());
      Matrix<> c_derivs(ir.Size(), cf->Dimension());

      SECTION("Compile")
        {
          c_vals = 0;
          c_derivs = 0;
          c_cf_f->EvaluateDeriv(mir,c_vals, c_derivs);
          CHECK(L2Norm(vals-c_vals) < 1e-12);
          CHECK(L2Norm(derivs-c_derivs) < 1e-12);
        }

      SECTION("TrueCompile")
        {
          c_vals = 0;
          c_derivs = 0;
          c_cf_t->EvaluateDeriv(mir,c_vals, c_derivs);
          CHECK(L2Norm(vals-c_vals) < 1e-12);
          CHECK(L2Norm(derivs-c_derivs) < 1e-12);
        }

      try
        {
          SECTION("SIMD")
            {
              simd_vals = 0;
              simd_derivs = 0;
              cf->EvaluateDeriv(simd_mir, simd_vals, simd_derivs);
              CHECK(L2Norm(vals-Trans(simd_vals)) < 1e-12);
              CHECK(L2Norm(derivs-Trans(simd_derivs)) < 1e-12);
              SECTION("Compile SIMD")
                {
                  simd_vals = 0;
                  simd_derivs = 0;
                  c_cf_f->EvaluateDeriv(simd_mir, simd_vals, simd_derivs);
                  CHECK(L2Norm(vals-Trans(simd_vals)) < 1e-12);
                  CHECK(L2Norm(derivs-Trans(simd_derivs)) < 1e-12);
                }

              SECTION("TrueCompile SIMD")
                {
                  simd_vals = 0;
                  simd_derivs = 0;
                  c_cf_t->EvaluateDeriv(simd_mir, simd_vals, simd_derivs);
                  CHECK(L2Norm(vals-Trans(simd_vals)) < 1e-12);
                  CHECK(L2Norm(derivs - Trans(simd_derivs)) < 1e-12);
                }
            }
        }
      catch (Exception ex)
        {
          CHECK_THROWS_WITH(throw ex, Catch::Matchers::Contains("EvaluateDeriv(simd) not overloaded"));
          WARN("SIMD not implemented");
        }
    }
}


TEST_CASE ("CoefficientFunctions","[fem][coefficient]")
{
  printmessage_importance = 1;
  SECTION("ConstantCoefficientFunction")
    {
      shared_ptr<CoefficientFunction> cf = make_shared<ConstantCoefficientFunction>(1);
      TestCoefficientFunction(cf);
    }

  auto a = make_shared<ConstantCoefficientFunction> (1);
  auto b = make_shared<ConstantCoefficientFunction> (2);
  auto y = MakeCoordinateCoefficientFunction(1);
  
  TEST_OPERATOR_COEFFICIENTFUNCTION(a+b+y);

  auto diffop = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
  auto u = make_shared<ProxyFunction>(false, false, diffop, nullptr, nullptr, nullptr, nullptr, nullptr);

  TEST_OPERATOR_COEFFICIENTFUNCTION(3*u*u);
  TEST_OPERATOR_COEFFICIENTFUNCTION(MakeVectorialCoefficientFunction({u,u*u}));
}
