#include <comp.hpp>

namespace ngcomp
{

  shared_ptr<BaseMatrix> ConvertOperator (shared_ptr<CompoundFESpace> comp_space,
					  int inda, int indb,
					  DifferentialOperator * diffop,
					  VorB vb, const Region * reg, LocalHeap & lh,
					  bool localop = false, bool parmat = true, bool use_simd = true)
  {
    /** Use compound space for element-coloring **/
    shared_ptr<FESpace> spacea = (*comp_space)[inda];
    shared_ptr<FESpace> spaceb = (*comp_space)[indb];

    if ( parmat && (spacea->IsParallel() != spaceb->IsParallel()) )
      { throw Exception("Cannot form ConvertOperator between a parallel and a local space!"); }

    if ( parmat && spacea->IsParallel() && spaceb->IsParallel() ) {
      MPI_Comm comma = spacea->GetParallelDofs()->GetCommunicator(), commb = spaceb->GetParallelDofs()->GetCommunicator();
      if (comma != commb)
	{ throw Execption("Cannot convert between spaces defined on different Communicators!"); }
    }
      

    /** Solves gfb = diffop(gfa), where gfb is a function from spaceb and gfa one from spacea **/

    /** Proxies and Integrators **/

    shared_ptr<ProxyFunction> trial_a;
    if ( !diffop ) // Probably most of the time
      { trial_a = spacea->GetEvaluator(vb); }
    else {
      trial_a = make_shared<ProxyFunction>(space_a, false, false, diffop->shared_from_this(),
					   nullptr, nullptr, nullptr, nullptr, nullptr);
    }

    shared_ptr<ProxyFunction> trial_b = spaceb->GetEvaluator(vb);

    if (!spaceb->GetAdditionalEvaluators().Used("dual"))
      throw Exception(string("Dual diffop does not exist for ") + spaceb->GetClassName() + string("!"));
    auto dual_evaluator = spaceb->GetAdditionalEvaluators()["dual"];
    for (VorB avb = VOL; avb < vb; avb++) {
      dual_evaluator = dual_evaluator->GetTrace();
      if ( dual_evaluator == nullptr )
	{ throw Exception(spaceb->GetClassName() + string(" has no dual trace operator for vb = ") + to_string(avb) + string("!")); }
    }
    dual_b = make_shared<ProxyFunction>(space_b, true, false, dual_evaluator,
					nullptr, nullptr, nullptr, nullptr, nullptr);
    
    if ( trial_a->Dimension() != trial_b->Dimension() )
      { throw Exception(string("Cannot convert from ") + space_a->GetClassName() + string(" to ") + space_b->GetClassName() +
			string(" - dimensions mismatch: ") + to_string(trial_a->Dimension()) +
			string(" != ") + to_string(trial_b->Dimension()) + string("!")); }

    Array<shared_ptr<BilinearFormIntegrator>> ab_bfis, bb_bfis;
    for (auto element_vb : spaceb->GetDualShapeNodes(vb)) {
      shared_ptr<CoefficientFunction> ab, bb;
      if (trial_b->Dimension() == 1) {
	ab = trial_a * dual_b;
	bb = trial_b * dual_b;
      }
      else {
	ab = InnerProduct(trial_a, dual_b);
	bb = InnerProduct(trial_b, dual_b);
      }
      auto ab_bfi = make_shared<BilinearFormIntegrator>(ab, vb, element_vb);
      ab_bfi->SetSimdEvaluate(use_simd);
      ab_bfis.Append(ab_bfi);
      auto bb_bfi = make_shared<BilinearFormIntegrator>(bb, vb, element_vb);
      bb_bfi->SetSimdEvaluate(use_simd);
      bb_bfis.Append(bb_bfi);
    }
    
    /** Utility **/
    auto it_els = [&](auto lam) {
      IterateElements(*comp_space, vb, lh,
		      [&](FESpace::Element ei, Localheap & lh)
      {
	if ( (!spacea->DefinedOn(vb, ei.GetIndex())) || (!spaceb->DefinedOn(vb, ei.GetIndex())) )
	  { return; }
	if ( reg && !reg->Mask().Test(ei.GetIndex()) )
	  { return; }
	lam(ei);
      });
    };

    /** Create Matrix Graph **/
    TableCreator<int> cgraph(spaceb->GetNDof());
    Array<DofId> dnums_a(100), dnums_b(100);
    for (; !cgraph.Done(); cgraph++ )
      it_els([&](auto & ei) {
	  ElementId eid(ei);
	  spacea->GetDofNrs(eid, dnums_a, ANY_DOF); // get rid of UNUSED DOFs
	  spaceb->GetDofNrs(eid, dnums_b, ANY_DOF); // get rid of UNUSED DOFs
	  for (auto db : dnums_b)
	    if (IsRegularDof(db))
	      for (auto da : dnums_a)
		if (IsRegularDof(da))
		  { cgraph.Add(db, da); }
	});

    /** Alloc Matrix **/
    auto graph = cgraph.MoveTable();
    auto spmat = make_shared<SparseMatrix<double>> (graph, spacea->GetNDof());
    spmat->AsVector() = 0;

    /** Fill Matrix **/

    Array<int> cntb(spmat->Height()); cntb = 0;

    auto fill_loop = [&]() {
	it_els([&](auto & ei) {
	    /** Get Finite Elements / Element Transformation **/
	    const FiniteElement & fel = ei.GetFE();
	    const auto & comp_fel = static_cast<const CompoundFiniteElement&>(fel);
	    const FiniteElement & fela = comp_fel[inda];
	    const FiniteElement & felb = comp_fel[indb];
	    MixedFiniteElement felab (fela, felb);
	    const ElementTransformation & eltrans = ei.GetTrafo();

	    /** Calc mass/dual and mixed mass/dual matrices **/
	    FlatMatrix<SCAL> bamat(felb.GetNDof(), fela.GetNDof(), lh); abmat = 0.0;
	    FlatMatrix<SCAL> bbmat(felb.GetNDof(), felb.GetNDof(), lh); bbmat = 0.0;
	    for (auto bfi : ab_bfis)
	      { bfi->CalcElementMatrixAdd(felab, eltrans, bamat, lh); }
	    for (auto bfi : bb_bfis)
	      { bfi->CalcElementMatrixAdd(felb, eltrans, bbmat, lh); }

	    /** Calc Elmat **/
	    CalcInverse(abmat); // NOT symmetric!!
	    FlatMatrix<SCAL> elmat(felb.GetNDof(), fela.GetNDof(), lh);
	    elmat = bbmat * bamat;

	    /** Write into SparseMatrix **/
	    spacea->GetDofNrs(eid, dnums_a, ANY_DOF); // get rid of UNUSED DOFs
	    spaceb->GetDofNrs(eid, dnums_b, ANY_DOF); // get rid of UNUSED DOFs
	    spmat->AddElementMatrix(dnums_b, dnums_a, elmat, false); // should be taken care of by fespace coloring

	    for (auto dnum : dnums_b)
	      if (IsRegularDof(dnum))
		{ cntb[dnum]++; }
	  });
    };
    if (use_simd) {
      try { fill_loop(); }
      catch (ExceptionNOSIMD e) { /** Turn off SIMD, try again! **/
	for (auto bfi : ab_bfis)
	  { bfi->SetSimdEvaluate(false); }
	for (auto bfi : bb_bfis)
	  { bfi->SetSimdEvaluate(false); }
	fill_loop();
      }
    } // use_simd
    else
      { fill_loop(); }

#ifdef PARALLEL
    if (spaceb->IsParallel() && !localop)
      { AllReduceDofData (cntb, MPI_SUM, spaceb->GetParallelDofs()); }
#endif

    for (auto rownr : Range(spmat->Height()))
      if (cnta[rownr] > 1) {
	double fac = 1.0 / double(cnta[rownr]);
	auto rvs = spmat->GetRowValues(rownr);
	for (auto & v : rvs)
	  { v &= fac; }
      }

    shared_ptr<BaseMatrix> op = spmat;

    if ( parmat && spacea->IsParallel() ) {
      op = make_shared<ParallelMatrix> (op, spacea->GetParallelDofs(), spaceb->GetParallelDofs(),
					localop ? PARALLEL_OP::C2C : PARALLEL_OP:C2D;);
    }

    return op;
  } // ConvertOperator


  shared_ptr<BaseMatrix> ConvertOperator (shared_ptr<FESpace> spacea, shared_ptr<FESpace> spaceb,
					  DifferentialOperator * diffop,
					  VorB vb, const Region * reg, LocalHeap & lh,
					  bool localop, bool parmat, bool use_simd)
  {
    Array<shared_ptr<FESpace>> spaces ( { spacea, spaceb } );
    auto cs = make_shared<CompoundFESpace>(spaces);
    return ConvertOperator (cs, 0, 1, diffop, vb, reg, lh, localop, parmat, use_simd);
  } // ConvertOperator

} // namespace ngcomp
