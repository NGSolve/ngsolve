#include <comp.hpp>

namespace ngcomp
{

  template<class SCAL>
  shared_ptr<BaseMatrix> ConvertOperator (shared_ptr<CompoundFESpace> comp_space,
					  int inda, int indb,
					  shared_ptr<DifferentialOperator> diffop,
					  VorB vb, const Region * reg, LocalHeap & clh,
					  bool localop = false, bool parmat = true, bool use_simd = true)
  {
    /** Solves gfb = diffop(gfa), where gfb is a function from space_b and gfa one from space_a **/

    /** Compound space is used for element-coloring **/
    shared_ptr<FESpace> space_a = (*comp_space)[inda];
    shared_ptr<FESpace> space_b = (*comp_space)[indb];

    if ( parmat && (space_a->IsParallel() != space_b->IsParallel()) )
      { throw Exception("Cannot form ConvertOperator between a parallel and a local space!"); }

    if ( parmat && space_a->IsParallel() && space_b->IsParallel() ) {
      MPI_Comm comma = space_a->GetParallelDofs()->GetCommunicator(), commb = space_b->GetParallelDofs()->GetCommunicator();
      if (comma != commb)
	{ throw Exception("Cannot convert between spaces defined on different Communicators!"); }
    }

    /** Proxies and Integrators **/

    shared_ptr<ProxyFunction> trial_a;
    if ( !diffop ) { // Probably most of the time
      trial_a = make_shared<ProxyFunction>(space_a, false, false, space_a->GetEvaluator(vb),
					   nullptr, nullptr, nullptr, nullptr, nullptr);
    }
    else {
      trial_a = make_shared<ProxyFunction>(space_a, false, false, diffop,
					   nullptr, nullptr, nullptr, nullptr, nullptr);
    }

    auto trial_b = make_shared<ProxyFunction>(space_a, false, false, space_b->GetEvaluator(vb),
					      nullptr, nullptr, nullptr, nullptr, nullptr);

    if (!space_b->GetAdditionalEvaluators().Used("dual"))
      throw Exception(string("Dual diffop does not exist for ") + space_b->GetClassName() + string("!"));
    auto dual_evaluator = space_b->GetAdditionalEvaluators()["dual"];
    for (VorB avb = VOL; avb < vb; avb++) {
      dual_evaluator = dual_evaluator->GetTrace();
      if ( dual_evaluator == nullptr )
	{ throw Exception(space_b->GetClassName() + string(" has no dual trace operator for vb = ") + to_string(avb) + string("!")); }
    }
    auto dual_b = make_shared<ProxyFunction>(space_b, true, false, dual_evaluator,
					     nullptr, nullptr, nullptr, nullptr, nullptr);
    
    if ( trial_a->Dimension() != trial_b->Dimension() )
      { throw Exception(string("Cannot convert from ") + space_a->GetClassName() + string(" to ") + space_b->GetClassName() +
			string(" - dimensions mismatch: ") + to_string(trial_a->Dimension()) +
			string(" != ") + to_string(trial_b->Dimension()) + string("!")); }

    Array<shared_ptr<BilinearFormIntegrator>> ab_bfis, bb_bfis;
    for (auto element_vb : space_b->GetDualShapeNodes(vb)) {
      shared_ptr<CoefficientFunction> ab, bb;
      if (trial_b->Dimension() == 1) {
	ab = trial_a * dual_b;
	bb = trial_b * dual_b;
      }
      else {
	ab = InnerProduct(trial_a, dual_b);
	bb = InnerProduct(trial_b, dual_b);
      }
      auto ab_bfi = make_shared<SymbolicBilinearFormIntegrator>(ab, vb, element_vb);
      ab_bfi->SetSimdEvaluate(use_simd);
      ab_bfis.Append(ab_bfi);
      auto bb_bfi = make_shared<SymbolicBilinearFormIntegrator>(bb, vb, element_vb);
      bb_bfi->SetSimdEvaluate(use_simd);
      bb_bfis.Append(bb_bfi);
    }
    
    // cout << " ab_bfis: " << endl << ab_bfis << endl;
    // cout << " bb_bfis: " << endl << bb_bfis << endl;

    /** Utility **/
    auto it_els = [&](auto lam) {
      // cout << " called it_els !" << endl;
      // cout << " VB " << vb << endl;
      // cout << " ma NE " << comp_space->GetMeshAccess()->GetNE(vb) << endl;
      IterateElements(*comp_space, vb, clh,
		      [&](FESpace::Element ei, LocalHeap & lh)
      {
	// cout << " use element " << ElementId(ei) << " ? " << endl;
	if ( (!space_a->DefinedOn(vb, ei.GetIndex())) || (!space_b->DefinedOn(vb, ei.GetIndex())) )
	  { return; }
	if ( reg && !reg->Mask().Test(ei.GetIndex()) )
	  { return; }
	// cout << " use element " << ElementId(ei) << " ! " << endl;
	lam(ei, lh);
      });
    };

    /** Create Matrix Graph **/
    int dima = space_a->GetDimension(), dimb = space_b->GetDimension();

    if ( (dima != 1) || (dimb != 1) )
      { throw Exception("multidim space convert TODO!!"); }

    TableCreator<int> cgraph(space_b->GetNDof() * dimb);
    Array<DofId> dnums_a(100), dnums_b(100);
    for (; !cgraph.Done(); cgraph++ )
      it_els([&](FESpace::Element & ei, LocalHeap & lh) {
	  ElementId eid(ei);
	  space_a->GetDofNrs(eid, dnums_a, ANY_DOF); // get rid of UNUSED DOFs
	  space_b->GetDofNrs(eid, dnums_b, ANY_DOF); // get rid of UNUSED DOFs
	  for (auto db : dnums_b)
	    if (IsRegularDof(db))
	      for (auto da : dnums_a)
		if (IsRegularDof(da)) {
		  for (auto lb : Range(dimb)) {
		    auto x = dimb * db + lb;
		    for (auto la : Range(dima))
		      { cgraph.Add(x, dima * da + la); }
		  }
		}
	});

    /** Alloc Matrix **/
    auto graph = cgraph.MoveTable();
    cout << " have graph: " << endl << graph << endl;
    Array<int> perow(space_b->GetNDof() * dimb);
    for (auto k : Range(perow))
      { perow[k] = graph[k].Size(); }
    auto spmat = make_shared<SparseMatrix<SCAL>> (perow, space_a->GetNDof() * dima);
    for (auto k : Range(perow)) {
      QuickSort(graph[k]);
      spmat->GetRowIndices(k) = graph[k];
    }
    spmat->AsVector() = 0;

    /** Fill Matrix **/
    Array<int> cnt_b(spmat->Height()); cnt_b = 0;
    auto fill_lam = [&](FESpace::Element & ei, LocalHeap & lh) {
      /** Get Finite Elements / Element Transformation **/
      const FiniteElement & fel = ei.GetFE();
      const auto & comp_fel = static_cast<const CompoundFiniteElement&>(fel);
      const FiniteElement & fela = comp_fel[inda];
      const FiniteElement & felb = comp_fel[indb];
      MixedFiniteElement felab (fela, felb);
      const ElementTransformation & eltrans = ei.GetTrafo();

      /** Calc mass/dual and mixed mass/dual matrices **/
      FlatMatrix<SCAL> bamat(felb.GetNDof() * dima, fela.GetNDof() * dima, lh); bamat = 0.0;
      FlatMatrix<SCAL> bbmat(felb.GetNDof() * dimb, felb.GetNDof() * dimb, lh); bbmat = 0.0;
      for (auto bfi : ab_bfis)
	{ bfi->CalcElementMatrixAdd(felab, eltrans, bamat, lh); }
      for (auto bfi : bb_bfis)
	{ bfi->CalcElementMatrixAdd(felb, eltrans, bbmat, lh); }

      /** Calc Elmat **/
      // cout << " bamat " << endl << bamat << endl;
      // cout << " bbmat " << endl << bbmat << endl;
      CalcInverse(bbmat); // NOT symmetric!!
      // cout << " bbmat inv " << endl << bbmat << endl;
      FlatMatrix<SCAL> elmat(felb.GetNDof(), fela.GetNDof(), lh);
      elmat = bbmat * bamat;
      // cout << " elmat " << endl << elmat << endl;

      /** Write into SparseMatrix **/
      space_a->GetDofNrs(ei, dnums_a, ANY_DOF); // get rid of UNUSED DOFs
      space_b->GetDofNrs(ei, dnums_b, ANY_DOF); // get rid of UNUSED DOFs

      // cout << " add elmat, rownums " << endl << dnums_b << endl;
      // cout << " add elmat, colnums " << endl << dnums_a << endl;

      spmat->AddElementMatrix(dnums_b, dnums_a, elmat, false); // should be taken care of by fespace coloring

      for (auto dnum : dnums_b)
	if (IsRegularDof(dnum))
	  { cnt_b[dnum]++; }
    };

    if (use_simd) {
      it_els([&](FESpace::Element & ei, LocalHeap & lh) {
	  try { fill_lam(ei, lh); }
	  catch (ExceptionNOSIMD e) { /** Turn off SIMD, try again! **/
	    for (auto bfi : ab_bfis)
	      { bfi->SetSimdEvaluate(false); }
	    for (auto bfi : bb_bfis)
	      { bfi->SetSimdEvaluate(false); }
	    fill_lam(ei, lh);
	  }
	});
    }
    else // use_simd
      { it_els(fill_lam); }

#ifdef PARALLEL
    if (space_b->IsParallel() && !localop)
      { AllReduceDofData (cnt_b, MPI_SUM, space_b->GetParallelDofs()); }
#endif

    for (auto rownr : Range(spmat->Height()))
      if (cnt_b[rownr] > 1) {
	double fac = 1.0 / double(cnt_b[rownr]);
	auto rvs = spmat->GetRowValues(rownr);
	for (auto & v : rvs)
	  { v *= fac; }
      }

    shared_ptr<BaseMatrix> op = spmat;

    if ( parmat && space_a->IsParallel() ) {
      op = make_shared<ParallelMatrix> (op, space_a->GetParallelDofs(), space_b->GetParallelDofs(),
					localop ? PARALLEL_OP::C2C : PARALLEL_OP::C2D);
    }

    return op;
  } // ConvertOperator


  shared_ptr<BaseMatrix> ConvertOperator (shared_ptr<FESpace> space_a, shared_ptr<FESpace> space_b, VorB vb, LocalHeap & lh,
					  shared_ptr<DifferentialOperator> diffop, const Region * reg,
					  bool localop, bool parmat, bool use_simd)
  {
    if ( space_a->IsComplex() != space_b->IsComplex() ) // b complex and a real could work in principle (?)
      { throw Exception("Cannot convert between complex and non-complex space!"); }

    Flags flags;
    auto cs = make_shared<CompoundFESpace>(space_a->GetMeshAccess(), flags);
    cs->AddSpace(space_a);
    cs->AddSpace(space_b);
    cs->Update();
    cs->FinalizeUpdate();

    shared_ptr<BaseMatrix> op;
    if (space_b->IsComplex())
      { op = ConvertOperator<Complex> (cs, 0, 1, diffop, vb, reg, lh, localop, parmat, use_simd); }
    else
      { op = ConvertOperator<double> (cs, 0, 1, diffop, vb, reg, lh, localop, parmat, use_simd); }

    return op;
  } // ConvertOperator

} // namespace ngcomp
