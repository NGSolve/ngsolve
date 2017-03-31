
#include "catch.hpp"
#include <fem.hpp>

using namespace ngfem;
  
#define MERGE_(a,b)  a##b
#define LABEL_(a) MERGE_(mynamespace, a)
#define UNIQUE_NAME LABEL_(__COUNTER__)

constexpr double tolerance = 1e-8;

#define TEST_OPERATOR_COEFFICIENTFUNCTION(CF)                 \
  namespace UNIQUE_NAME {                                          \
    auto c_cf_f = Compile(CF,false);                               \
    auto c_cf_t = Compile(CF,true);                                \
    TEST_CASE(#CF) {                                               \
      TestCoefficientFunction(CF,c_cf_f,c_cf_t);                   \
    }                                                              \
  }

void TestCoefficientFunction(shared_ptr<CoefficientFunction> cf,
                             shared_ptr<CoefficientFunction> c_cf_f,
                             shared_ptr<CoefficientFunction> c_cf_t)
{
  LocalHeap lh(100000, "lh");
  IntegrationRule ir(ET_TET, 3);
  FE_ElementTransformation<3,3> trafo(ET_TET);
  MappedIntegrationRule<3,3> mir(ir, trafo, lh);

  SIMD_IntegrationRule simd_ir(ET_TET, 3);
  SIMD_MappedIntegrationRule<3,3> simd_mir(simd_ir, trafo, lh);

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
  Matrix<> derivs(ir.Size(), cf->Dimension());
  cf->EvaluateDeriv(mir,vals,derivs);
  AMatrix<> simd_derivs(cf->Dimension(),simd_ir.GetNIP());
  Matrix<> c_derivs(ir.Size(), cf->Dimension());

  SECTION ("Evaluate")
    {
      SECTION("Compile")
        {
          c_vals = 0;
          c_cf_f->Evaluate(mir,c_vals);
          CHECK(L2Norm(vals-c_vals) < tolerance);
        }

      SECTION("TrueCompile")
        {
          c_vals = 0;
          c_cf_t->Evaluate(mir,c_vals);
          CHECK(L2Norm(vals-c_vals) < tolerance);
        }

      try
        {
          SECTION("SIMD")
            {
              simd_vals = 0;
              cf->Evaluate(simd_mir, simd_vals);
              CHECK(L2Norm(vals-Trans(simd_vals)) < tolerance);

              SECTION("Compile")
                {
                  simd_vals = 0;
                  c_cf_f->Evaluate(simd_mir, simd_vals);
                  CHECK(L2Norm(vals-Trans(simd_vals)) < tolerance);
                }

              SECTION("TrueCompile")
                {
                  simd_vals = 0;
                  c_cf_t->Evaluate(simd_mir, simd_vals);
                  CHECK(L2Norm(vals-Trans(simd_vals)) < tolerance);
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

      SECTION("Compile")
        {
          c_vals = 0;
          c_derivs = 0;
          c_cf_f->EvaluateDeriv(mir,c_vals, c_derivs);
          CHECK(L2Norm(vals-c_vals) < tolerance);
          CHECK(L2Norm(derivs-c_derivs) < tolerance);
        }

      SECTION("TrueCompile")
        {
          c_vals = 0;
          c_derivs = 0;
          c_cf_t->EvaluateDeriv(mir,c_vals, c_derivs);
          CHECK(L2Norm(vals-c_vals) < tolerance);
          CHECK(L2Norm(derivs-c_derivs) < tolerance);
        }

      try
        {
          SECTION("SIMD")
            {
              SECTION("no-compile SIMD")
                {
                  simd_vals = 0;
                  simd_derivs = 0;
                  cf->EvaluateDeriv(simd_mir, simd_vals, simd_derivs);
                  CHECK(L2Norm(vals-Trans(simd_vals)) < tolerance);
                  CHECK(L2Norm(derivs-Trans(simd_derivs)) < tolerance);
                }
              SECTION("Compile SIMD")
                {
                  simd_vals = 0;
                  simd_derivs = 0;
                  c_cf_f->EvaluateDeriv(simd_mir, simd_vals, simd_derivs);
                  CHECK(L2Norm(vals-Trans(simd_vals)) < tolerance);
                  CHECK(L2Norm(derivs-Trans(simd_derivs)) < tolerance);
                }

              SECTION("TrueCompile SIMD")
                {
                  simd_vals = 0;
                  simd_derivs = 0;
                  c_cf_t->EvaluateDeriv(simd_mir, simd_vals, simd_derivs);
                  CHECK(L2Norm(vals-Trans(simd_vals)) < tolerance);
                  CHECK(L2Norm(derivs - Trans(simd_derivs)) < tolerance);
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


int set_printmessage_level = [](){
    printmessage_importance = 1;
    return 0;
}();

auto a = make_shared<ConstantCoefficientFunction> (1);
auto b = make_shared<ConstantCoefficientFunction> (2);
auto x = MakeCoordinateCoefficientFunction(0);
auto y = MakeCoordinateCoefficientFunction(1);
auto z = MakeCoordinateCoefficientFunction(2);
auto p = make_shared<ParameterCoefficientFunction>(42.39);

TEST_OPERATOR_COEFFICIENTFUNCTION(p);
TEST_OPERATOR_COEFFICIENTFUNCTION(a*x+y/(b+z));

auto diffop = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
auto u = make_shared<ProxyFunction>(false, false, diffop, nullptr, nullptr, nullptr, nullptr, nullptr);

TEST_OPERATOR_COEFFICIENTFUNCTION(3*u*u);
TEST_OPERATOR_COEFFICIENTFUNCTION(MakeVectorialCoefficientFunction({u,u*u}));

auto xyz = MakeVectorialCoefficientFunction({x,u,z});
auto mat = MakeVectorialCoefficientFunction({a,b,x,y,z,b,u,a,b});
int set_dimension = [](auto mat) {
mat->SetDimensions (Array<int>({3,3}));
return 1; }(mat);
TEST_OPERATOR_COEFFICIENTFUNCTION(xyz);
// TEST_OPERATOR_COEFFICIENTFUNCTION(NormCF(xyz));
TEST_OPERATOR_COEFFICIENTFUNCTION(InnerProduct(xyz,xyz));
TEST_OPERATOR_COEFFICIENTFUNCTION(MakeComponentCoefficientFunction(xyz, 1));
TEST_OPERATOR_COEFFICIENTFUNCTION(TransposeCF(mat));
TEST_OPERATOR_COEFFICIENTFUNCTION(mat*xyz);

auto longvec = MakeVectorialCoefficientFunction({x,y,z,x,y,z,z,x,y,x,x,x});
// TEST_OPERATOR_COEFFICIENTFUNCTION(InnerProduct(longvec, longvec));
