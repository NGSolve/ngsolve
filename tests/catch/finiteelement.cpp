
#include <catch.hpp>
#include <fem.hpp>

using namespace ngfem;

// visual studio wants this...
using ngfem::ELEMENT_TYPE;

// workaround for visual compiler bug ---------------------------
template<ELEMENT_TYPE ET>
using H1 = H1HighOrderFE<
  ET,
  H1HighOrderFE_Shape<ET>,
  T_ScalarFiniteElement<
    H1HighOrderFE_Shape<ET>,
    ET,
    ScalarFiniteElement<ET_trait<ET>::DIM>
    >
  >;

static_assert(is_same<H1HighOrderFE<ET_TRIG>,H1<ET_TRIG>>::value,"Default arguments for H1HoFE have changed, update catch unit tests finiteelement.cpp");

template<ELEMENT_TYPE ET>
using L2 = L2HighOrderFE<
  ET,
  L2HighOrderFE_Shape<ET>,
  T_ScalarFiniteElement<
    L2HighOrderFE_Shape<ET>,
    ET,
    DGFiniteElement<ET>
    >
  >;

static_assert(is_same<L2HighOrderFE<ET_TRIG>,L2<ET_TRIG>>::value,"Default arguments for L2HoFE have changed, update catch unit tests finiteelement.cpp");

// end workaround -------------------------------------------------

// Run f for ELEMENT_TYPE ET1 (used as the end of variadic template ForET)
template<ELEMENT_TYPE ET1, typename FUNC>
void ForET(FUNC f)
{
  SECTION ("ELEMENT_TYPE = " + std::string(ElementTopology::GetElementName(ET1)),"")
    {
      f(ET_trait<ET1>());
    }
}

// Run function f for all ELEMENT_TYPES specified in template argument
template<ELEMENT_TYPE ET1, ELEMENT_TYPE ET2, ELEMENT_TYPE ... ET_REST, typename FUNC>
void ForET(FUNC f)
{
  SECTION ("ELEMENT_TYPE = " + std::string(ElementTopology::GetElementName(ET1)),"")
    {
      f(ET_trait<ET1>());
    }
  ForET<ET2,ET_REST...>(f);
}

// run f for all ELEMENT_TYPEs
template<typename FUNC>
void ForET (FUNC f)
{
  ForET<ET_POINT,ET_SEGM,ET_TRIG,ET_QUAD,ET_TET,ET_PRISM,ET_PYRAMID,ET_HEX>(f);
}

template<typename T>
void TestFiniteElement(T fel)
{
          IntegrationRule ir(fel.ElementType(),fel.Order());

          SECTION ("Evaluate", "[evaluate]")
            {
              Vector<> coefs(fel.GetNDof()), values(ir.Size());
              for (auto i : Range(coefs.Size()))
                coefs[i] = i*0.5;

              fel.Evaluate(ir,coefs,values);

              SECTION ("SIMD", "[SIMD]")
                {
                  SIMD_IntegrationRule simdir(fel.ElementType(),fel.Order());
                  Vector<SIMD<double>> simd_values(simdir.Size());

                  try
                    {
                      fel.Evaluate(simdir,coefs,simd_values);
                      FlatVector<double> values_ref(values.Size(),(double*)&simd_values[0]);
                      SECTION("SIMD correctness", "[SIMD]")
                        {
                          CHECK(L2Norm(values-values_ref) < 1e-10);
                        }
                    }
                  catch(const Exception & ex)
                    {
                      CHECK(typeid(ex) == typeid(ExceptionNOSIMD));
                      WARN("SIMD not implemented");
                    }
                }
            }
}

TEST_CASE ("FiniteElement", "[fem][finiteelement]")
{
  SECTION ("H1HighOrder", "[h1]")
    {
      //
      ForET([&](auto ET) {
          for (auto order : Range(1,8)) {
            SECTION ("order = " + std::to_string(order),"")
              {
                TestFiniteElement(H1<ET.ElementType()>(order));
              }
          }
        });
    }

  SECTION ("L2HighOrder", "[l2]")
    {
      ForET([&](auto ET) {
          for (auto order : Range(8)) {
            SECTION ("order = " + std::to_string(order),"")
              {
                TestFiniteElement(L2<ET.ElementType()>(order));
              }
          }
        });
    }
}
