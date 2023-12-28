/*********************************************************************/
/* File:   hdivdivfe.cpp                                             */
/* Author: Astrid Pechstein, Joachim Schoeberl                       */
/* Date:   orig 2006, redesign Dec 2016                              */
/*********************************************************************/

#define FILE_HDIVDIVFE_CPP

// #include <fem.hpp>
#include "scalarfe.hpp"
#include "hdivdivfe.hpp"

namespace ngfem
{
  template<int D>
  list<tuple<string,double>> HDivDivFiniteElement<D> :: Timing () const
  {
    list<tuple<string,double>>timings;
    IntegrationRule ir(ElementType(), 2*Order());
    SIMD_IntegrationRule simdir(ElementType(), 2*Order());
    Vector<> coefs(GetNDof());
    Matrix<> shape(GetNDof(),(D+1)*D/2);
    Matrix<> divshape(GetNDof(),D);
    Vector<> values(ir.Size());
    Matrix<> dvalues(ir.Size(), D);
    Matrix<SIMD<double>> simd_shapes(D*D*GetNDof(), simdir.Size());
    FE_ElementTransformation<D,D> trafo(ElementType());
    static LocalHeap lh (10000000, "FE - Timing");
    HeapReset hr(lh);
    // auto & mir = trafo(ir, lh);
    auto & mir = trafo(ir, lh);    
    auto & simdmir = trafo(simdir, lh);

    coefs = 1;
    
    double maxtime = 0.5;
    double time;

    constexpr size_t steps = 1000;
    time = RunTiming([&]() {
        for (size_t i = 0; i < steps; i++)
          for (size_t j = 0; j < ir.Size(); j++)
            this -> CalcShape(ir[j], shape);
      });
    timings.push_back(make_tuple("CalcShape", time/steps*1e9/(D*(D+1)/2*GetNDof()*ir.Size())));

    time = RunTiming([&]() {
        for (size_t i = 0; i < steps; i++)
          for (size_t j = 0; j < ir.Size(); j++)
            this -> CalcDivShape(ir[j], divshape);
      });
    timings.push_back(make_tuple("CalcDivShape", time/steps*1e9/(D*GetNDof()*ir.Size())));

    time = RunTiming([&]() {
        for (size_t i = 0; i < steps; i++)
          for (size_t j = 0; j < ir.Size(); j++)
            this -> CalcMappedDivShape(mir[j], divshape);
      });
    timings.push_back(make_tuple("CalcMappedDivShape", time/steps*1e9/(D*GetNDof()*ir.Size())));

    
    time = RunTiming([&]() {
        for (size_t i = 0; i < steps; i++)
          this -> CalcMappedShape_Matrix(simdmir, simd_shapes);
      }, maxtime);
    timings.push_back(make_tuple("CalcShape (SIMD)", time/steps*1e9/(D*D*GetNDof()*simdir.GetNIP())));
    cout << "simd_shape mem = " << simd_shapes.Height()*simd_shapes.Width()*sizeof(SIMD<double>) << endl;
    return timings;
  }

  template class HDivDivFiniteElement<2>;
  template class HDivDivFiniteElement<3>;

}
