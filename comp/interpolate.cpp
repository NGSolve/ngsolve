/**********************************************************************/
/* File:   interpolate.cpp                                            */
/* Author: L Kogler, M Neunteufel, J Schoeberl                        */
/* Date:   June 2020                                                  */
/**********************************************************************/

/* 
   Interpolation of CoefficientFunctions using
   dual shapes
*/

#include <comp.hpp>
#include <variant>
#include "../fem/integratorcf.hpp"


namespace ngcomp
{
  
  class InterpolationCoefficientFunction : public T_CoefficientFunction<InterpolationCoefficientFunction>
  {
  protected:
    shared_ptr<CoefficientFunction> func;
    shared_ptr<FESpace> fes;


    Array<shared_ptr<BilinearFormIntegrator>> bli;
    Array<shared_ptr<BilinearFormIntegrator>> single_bli;
    shared_ptr<CoefficientFunction> dual_diffop;
  public:
    InterpolationCoefficientFunction (shared_ptr<CoefficientFunction> f, shared_ptr<FESpace> afes ) :
      T_CoefficientFunction<InterpolationCoefficientFunction>(f->Dimension(), f->IsComplex()),
      func(f), fes(afes)
    {
      this->SetDimensions (func->Dimensions());
      this->elementwise_constant = func->ElementwiseConstant();


      // copied from Set (dualshapes)

      
      VorB vb = VOL;    // for the moment only 
      /** Trial-Proxy **/
      // if (!fes->GetEvaluator(vb))
      // throw Exception(fes->GetClassName()+string(" does not have an evaluator for ")+ToString(vb)+string("!"));
      
      auto single_evaluator =  fes->GetEvaluator(vb);
      if (dynamic_pointer_cast<BlockDifferentialOperator>(single_evaluator))
        single_evaluator = dynamic_pointer_cast<BlockDifferentialOperator>(single_evaluator)->BaseDiffOp();
      auto trial = make_shared<ProxyFunction>(fes, false, false, single_evaluator,
                                              nullptr, nullptr, nullptr, nullptr, nullptr);

      /** Test-Proxy (dual) **/
      auto dual_evaluator = fes->GetAdditionalEvaluators()["dual"];
      for (VorB avb = VOL; avb < vb; avb++) {
        dual_evaluator = dual_evaluator->GetTrace();
        if ( dual_evaluator == nullptr )
          { throw Exception(fes->GetClassName() + string(" has no dual trace operator for vb = ") + \
                            to_string(avb) + string(" -> ") + to_string(avb + 1) + string("!")); }
      }
      if (dynamic_pointer_cast<BlockDifferentialOperator>(dual_evaluator))
        dual_evaluator = dynamic_pointer_cast<BlockDifferentialOperator>(dual_evaluator)->BaseDiffOp();
      auto dual = make_shared<ProxyFunction>(fes, true, false, dual_evaluator,
                                             nullptr, nullptr, nullptr, nullptr, nullptr);

      dual_diffop = dual;

      
      
      for (auto element_vb : fes->GetDualShapeNodes(vb))
        {
          shared_ptr<CoefficientFunction> dual_trial;
          if (dual -> Dimension() == 1)
            { dual_trial = dual * trial; }
          else
            { dual_trial = InnerProduct(dual, trial); }
          auto bfi = make_shared<SymbolicBilinearFormIntegrator> (dual_trial, vb, element_vb);
	  bfi->SetSimdEvaluate(false);  // dual are not simded, yet
	  bli.Append(bfi);
	  if (auto block_bfi = dynamic_pointer_cast<BlockBilinearFormIntegrator> (bfi)) {
	    auto sbfi = block_bfi->BlockPtr();
	    sbfi->SetSimdEvaluate(false);
	    single_bli.Append(sbfi);
	  }
	  else
	    { single_bli.Append(bfi); }
	}
    }


    
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
    {
      func->Evaluate(ir, values);


#ifdef FIRSTDRAFT
      LocalHeapMem<100000> lh("interpolate");
    
    
    const ElementTransformation & trafo = ir.GetTransformation();
    const MeshAccess & ma = *static_cast<const MeshAccess*> (trafo.GetMesh());
    ElementId ei = trafo.GetElementId();
    auto & fel = fes->GetFE(ei, lh);

    for (auto el_vb : fes->GetDualShapeNodes(trafo.VB()))
      {
        Facet2ElementTrafo f2el (fel.ElementType(), el_vb);
        for (int locfnr : Range(f2el.GetNFacets()))
          {
            SIMD_IntegrationRule irfacet(f2el.FacetType(locfnr), 2 * fel.Order());
            auto & irvol = f2el(locfnr, irfacet, lh);
            auto & mir = trafo(irvol, lh);
            mir.ComputeNormalsAndMeasure(fel.ElementType(), locfnr);
            // do_ir(mir);

            FlatMatrix<T,ORD> mflux(dim, irfacet.Size(), lh);
            func->Evaluate (mir, mflux);
            for (size_t j : Range(mir))
              mflux.Col(j) *= mir[j].GetWeight();
            diffop->AddTrans (fel, mir, mflux, rhs);
          }
      }
    
    /** Calc Element Matrix **/
    FlatMatrix<double> elmat(fel.GetNDof(), lh);
    elmat = 0.0;
    bool symmetric_so_far = true;
    for (auto sbfi : single_bli)
      sbfi->CalcElementMatrixAdd (fel, eltrans, elmat, symmetric_so_far, lh);
    
    /** Invert Element Matrix and Solve for RHS **/
    CalcInverse(elmat); // Not Symmetric !
    
#endif
    
    }

  template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,
                     BareSliceMatrix<T,ORD> values) const
  {
    T_Evaluate (ir, values);

    
  }

  Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ func }); }

  void PrintReport (ostream & ost) const override
    {
      ost << "InterpolationCF(";
      func->PrintReport(ost);
      ost << ")";
    }

  string GetDescription() const override
    {
      return "InterpolationCF";
    }

  void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> nonzero) const override
    {
      func->NonZeroPattern(ud, nonzero);
    }

  void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
    {
      func->NonZeroPattern(ud, input, values);
    }

  void TraverseTree (const function<void(CoefficientFunction&)> & func_) override
    {
      func->TraverseTree(func_);
      func_(*this);
    }
  };
  
  
  shared_ptr<CoefficientFunction> InterpolateCF (shared_ptr<CoefficientFunction> func, shared_ptr<FESpace> space)
  {
    cout << "called interpolate" << endl;
    cout << "func = " << endl;
    func->PrintReport(cout);
    return make_shared<InterpolationCoefficientFunction> (func, space);
  }

}
