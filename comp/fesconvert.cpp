#include <comp.hpp>

namespace ngcomp
{

  template<class SCAL>
  shared_ptr<BaseMatrix> ConvertOperator (shared_ptr<FESpace> space_a, shared_ptr<FESpace> space_b,
					  // int inda, int indb,
					  shared_ptr<DifferentialOperator> diffop,
					  VorB vb, const Region * reg, LocalHeap & clh,
					  bool localop = false, bool parmat = true, bool use_simd = true)
  {
    /** Solves gfb = diffop(gfa), where gfb is a function from space_b and gfa one from space_a **/

    // /** Compound space is used for element-coloring **/
    // shared_ptr<FESpace> space_a = (*comp_space)[inda];
    // shared_ptr<FESpace> space_b = (*comp_space)[indb];
    auto ma = space_b->GetMeshAccess();

    auto ngmesh = ma->GetNetgenMesh();
    int curve_order = ngmesh->GetCurvedElements().GetOrder();

    int bonus_int_order = (curve_order > 1) ? ma->GetDimension() - 1 : 0;

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
	{ throw Exception(space_b->GetClassName() + string(" has no dual trace operator for vb = ")
			  + to_string(avb) + string(" -> ") + to_string(avb + 1) + string("!")); }
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
      if ( element_vb != VOL ) { // SEEMS to only be a problem for BND (and maybe BBND)
	ab_bfi->SetBonusIntegrationOrder(bonus_int_order);
	bb_bfi->SetBonusIntegrationOrder(bonus_int_order);
      }
    }

    /** Utility **/

    /** Create Matrix Graph **/
    int dima = space_a->GetDimension(), dimb = space_b->GetDimension();

    // if ( (dima != 1) || (dimb != 1) )
      // { throw Exception("multidim space convert TODO!!"); }

    int nthreads = (task_manager == nullptr) ? 1 : task_manager->GetNumThreads();
    Array<int> tmdsa(nthreads), tmdsb(nthreads);
    tmdsa = 0; tmdsb = 0;
    TableCreator<int> crnrs(ma->GetNE(vb)), ccnrs(ma->GetNE(vb));
    for (; !crnrs.Done(); crnrs++, ccnrs++ )
      ParallelForRange
	(ma->GetNE(vb), [&](IntRange r)
	 {
	   Array<DofId> dnums_a(100), dnums_b(100);
	   int tid = (task_manager == nullptr) ? 0 : task_manager->GetThreadId();
	   int& maxdsa = tmdsa[tid], maxdsb = tmdsb[tid];
	   for (auto i : r) {
	     ElementId eid(vb, i);
	     Ngs_Element el = ma->GetElement(eid);
	     if ( (!space_a->DefinedOn(vb, el.GetIndex())) || (!space_b->DefinedOn(vb, el.GetIndex())) )
	       { return; }
	     if ( reg && !reg->Mask().Test(el.GetIndex()) )
	       { return; }
	     space_a->GetDofNrs(eid, dnums_a, ANY_DOF); // get rid of UNUSED DOFs
	     maxdsa = max2(maxdsa, int(dnums_a.Size()));
	     for (auto da : dnums_a)
	       if (IsRegularDof(da)) {
		 auto ba = dima * da;
		 for (auto l : Range(dima))
		   { ccnrs.Add(i, ba + l); }
	       }
	     space_b->GetDofNrs(eid, dnums_b, ANY_DOF); // get rid of UNUSED DOFs
	     maxdsb = max2(maxdsb, int(dnums_b.Size()));
	     for (auto db : dnums_b)
	       if (IsRegularDof(db)) {
		 auto bb = dimb * db;
		 for (auto l : Range(dimb))
		   { crnrs.Add(i, bb + l); }
	       }
	   }
	 });
    int maxds_a = 0, maxds_b = 0;
    for (auto k : Range(tmdsa)) {
      maxds_a = max2(maxds_a, tmdsa[k]);
      maxds_b = max2(maxds_b, tmdsb[k]);
    }

    /** Alloc Matrix **/
    Table<int> rnrs = crnrs.MoveTable(), cnrs = ccnrs.MoveTable();
    MatrixGraph graph (space_b->GetNDof() * dimb, space_a->GetNDof() * dimb, rnrs, cnrs, false);

    auto spmat = make_shared<SparseMatrix<SCAL>> (graph, true);
    spmat->AsVector() = 0;

    /** Fill Matrix - with lambdas, we can use the same code for SIMD and non-SIMD cases **/
    Array<int> cnt_b(space_b->GetNDof()); cnt_b = 0;
    auto it_els = [&](auto lam) {
      IterateElements(*space_b, vb, clh, // !! use space_b element coloring !!
		      [&](FESpace::Element fei, LocalHeap & lh)
      {
	if ( (!space_a->DefinedOn(vb, fei.GetIndex())) || (!space_b->DefinedOn(vb, fei.GetIndex())) )
	  { return; }
	if ( reg && !reg->Mask().Test(fei.GetIndex()) )
	  { return; }
	lam(fei, lh);
      });
    };
    auto fill_lam = [&](FESpace::Element & fei, LocalHeap & lh) {
      /** Get Finite Elements / Element Transformation **/
      const ElementTransformation & eltrans = fei.GetTrafo();
      ElementId ei(fei);
      const FiniteElement & fela = space_a->GetFE(ei, lh);
      const FiniteElement & felb = fei.GetFE();
      MixedFiniteElement felab (fela, felb);

      Array<int> dnums_a(maxds_a, lh), dnums_b(maxds_b, lh);
      space_a->GetDofNrs(ei, dnums_a, ANY_DOF); // get rid of UNUSED DOFs [[ ACTUALLY, NO, also need right elmat dnrs! ]]
      space_b->GetDofNrs(ei, dnums_b, ANY_DOF); // get rid of UNUSED DOFs

      if (dnums_b.Size() == 0) // (compressed space)
	{ return; }

      /** Calc mass/dual and mixed mass/dual matrices **/
      FlatMatrix<SCAL> bamat(felb.GetNDof() * dimb, fela.GetNDof() * dima, lh); bamat = 0.0;
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

      FlatMatrix<SCAL> elmat(felb.GetNDof() * dimb, fela.GetNDof() * dima, lh);
      elmat = bbmat * bamat;

      // cout << " elmat " << endl << elmat << endl;

      // cout << " add elmat, dofs B " << endl << dnums_b << endl;
      // cout << " add elmat, dofs A " << endl << dnums_a << endl;

      /** Write into SparseMatrix **/
      if(  (dima == 1) && (dimb == 1) )
	{ spmat->AddElementMatrix(dnums_b, dnums_a, elmat, false); } // should be taken care of by fespace coloring
      else
	{
	  FlatArray<int> mddb(dnums_b.Size() * dimb, lh);
	  for (auto k : Range(dnums_b))
	    for (auto l : Range(dimb))
	      { mddb[k * dimb + l] = dimb * dnums_b[k] + l; }
	  FlatArray<int> mdda(dnums_a.Size() * dima, lh);
	  for (auto k : Range(dnums_a))
	    for (auto l : Range(dima))
	      { mdda[k * dima + l] = dima * dnums_a[k] + l; }
	  spmat->AddElementMatrix(mddb, mdda, elmat, false);
	}

      for (auto dnum : dnums_b) // space_b element coloring!!
	if (IsRegularDof(dnum))
	  { cnt_b[dnum]++; }
    };

    if (use_simd) {
      it_els([&](FESpace::Element & ei, LocalHeap & lh) {
	  try { fill_lam(ei, lh); }
	  catch (ExceptionNOSIMD e) { /** Turn off SIMD and continue **/
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

    for (auto dofnr : Range(size_t(spmat->Height()/dimb))) {
      if (cnt_b[dofnr] > 1) {
	double fac = 1.0 / double(cnt_b[dofnr]);
	for (auto rownr : Range(dimb * dofnr, dimb * (dofnr+1))) {
	  auto rvs = spmat->GetRowValues(rownr);
	  for (auto & v : rvs)
	    { v *= fac; }
	}
      }
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

    // Flags flags;
    // auto cs = make_shared<CompoundFESpace>(space_a->GetMeshAccess(), flags);
    // cs->AddSpace(space_a);
    // cs->AddSpace(space_b);
    // cs->Update();
    // cs->FinalizeUpdate();

    shared_ptr<BaseMatrix> op;
    if (space_b->IsComplex())
      { op = ConvertOperator<Complex> (space_a, space_b, diffop, vb, reg, lh, localop, parmat, use_simd); }
    else
      { op = ConvertOperator<double> (space_a, space_b, diffop, vb, reg, lh, localop, parmat, use_simd); }

    return op;
  } // ConvertOperator

} // namespace ngcomp
