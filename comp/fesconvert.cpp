#include <comp.hpp>

namespace ngcomp
{

  template<int BSA, int BSB, typename SCAL> struct TM_TRAIT { typedef Mat<BSA, BSB, SCAL> type; };
  template<> struct TM_TRAIT<1, 1, double> { typedef double type; };
  template<> struct TM_TRAIT<1, 1, Complex> { typedef Complex type; };

  template<int BSA, int BSB, typename SCAL> using TM = typename TM_TRAIT<BSA, BSB, SCAL>::type;
  template<int BSA, int BSB, typename SCAL> using TSPM = SparseMatrix<TM<BSA, BSB, SCAL>>;

  template<class SCAL, int DIMA, int DIMB>
  shared_ptr<BaseMatrix> ConvertOperator (shared_ptr<FESpace> space_a, shared_ptr<FESpace> space_b,
					  // int inda, int indb,
					  shared_ptr<DifferentialOperator> diffop, shared_ptr<CoefficientFunction> trial_cf,
					  VorB vb, const Region * reg, LocalHeap & clh,
					  shared_ptr<BitArray> range_dofs = nullptr,
					  bool localop = false, bool parmat = true, bool use_simd = true,
					  int bonus_intorder_ab = 0, int bonus_intorder_bb = 0)
  {
    static Timer t ("ConvertOperator");
    RegionTimer regt(t);
    static Timer tass("CovnertOperator - assemble");
    static Timer telmat("CovnertOperator - elmat");
    
    /** Solves gfb = diffop(gfa), where gfb is a function from space_b and gfa one from space_a **/

    /** Other Sparse-Matrix templates are not instantiated **/
    static_assert( (DIMA <= MAX_SYS_DIM) && (DIMB <= MAX_SYS_DIM) && ( (DIMA == DIMB) || (DIMA == 1) || (DIMB == 1) ) );

    // /** Compound space is used for element-coloring **/
    // shared_ptr<FESpace> space_a = (*comp_space)[inda];
    // shared_ptr<FESpace> space_b = (*comp_space)[indb];
    auto ma = space_b->GetMeshAccess();

    if ( parmat && (space_a->IsParallel() != space_b->IsParallel()) )
      { throw Exception("Cannot form ConvertOperator between a parallel and a local space!"); }

    if ( parmat && space_a->IsParallel() && space_b->IsParallel() ) {
      MPI_Comm comma = space_a->GetParallelDofs()->GetCommunicator(), commb = space_b->GetParallelDofs()->GetCommunicator();
      if (comma != commb)
	{ throw Exception("Cannot convert between spaces defined on different Communicators!"); }
    }

    /** Proxies and Integrators **/

    shared_ptr<CoefficientFunction> trial_a;
    if (trial_cf != nullptr)
      { trial_a = trial_cf; }
    else if ( diffop != nullptr ) {
      trial_a = make_shared<ProxyFunction>(space_a, false, false, diffop,
					   nullptr, nullptr, nullptr, nullptr, nullptr);
    }
    else { // Probably most of the time
      trial_a = make_shared<ProxyFunction>(space_a, false, false, space_a->GetEvaluator(vb),
					   nullptr, nullptr, nullptr, nullptr, nullptr);
    }

    auto trial_b = make_shared<ProxyFunction>(space_b, false, false, space_b->GetEvaluator(vb),
					      nullptr, nullptr, nullptr, nullptr, nullptr);

    if (!space_b->GetAdditionalEvaluators().Used("dual"))
      throw Exception(string("Dual diffop does not exist for ") + space_b->GetClassName() + string("!"));
    auto dual_evaluator = space_b->GetAdditionalEvaluators()["dual"];
    for (VorB avb = dual_evaluator->VB(); avb < vb; avb++) {
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
      ab_bfi->SetBonusIntegrationOrder(bonus_intorder_ab);
      ab_bfis.Append(ab_bfi);
      auto bb_bfi = make_shared<SymbolicBilinearFormIntegrator>(bb, vb, element_vb);
      bb_bfi->SetSimdEvaluate(use_simd);
      bb_bfis.Append(bb_bfi);
      bb_bfi->SetBonusIntegrationOrder(bonus_intorder_bb);
    }

    /** Utility **/

    /** Create Matrix Graph **/
    int dima = space_a->GetDimension(), dimb = space_b->GetDimension();

    if ( (dima != DIMA) || (dimb != DIMB) )
      { throw Exception("Dimensons and template parameters do not match in ConvertOperator!"); }

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
	       { continue; }
	     if ( reg && !reg->Mask().Test(el.GetIndex()) )
	       { continue; }
	     space_a->GetDofNrs(eid, dnums_a, ANY_DOF); // get rid of UNUSED DOFs
	     maxdsa = max2(maxdsa, int(dnums_a.Size()));
	     for (auto da : dnums_a)
	       if (IsRegularDof(da))
		 { ccnrs.Add(i, da); }
	     space_b->GetDofNrs(eid, dnums_b, ANY_DOF); // get rid of UNUSED DOFs
	     maxdsb = max2(maxdsb, int(dnums_b.Size()));
	     for (auto db : dnums_b)
	       if (IsRegularDof(db) && (!range_dofs || range_dofs->Test(db)) )
		 { crnrs.Add(i, db); }
	   }
	 });
    int maxds_a = 0, maxds_b = 0;
    for (auto k : Range(tmdsa)) {
      maxds_a = max2(maxds_a, tmdsa[k]);
      maxds_b = max2(maxds_b, tmdsb[k]);
    }


    /** Alloc Matrix and Fill Matrix - with lambdas, we can use the same code for SIMD and non-SIMD cases **/
    Table<int> rnrs = crnrs.MoveTable(), cnrs = ccnrs.MoveTable();
    MatrixGraph graph (space_b->GetNDof(), space_a->GetNDof(), rnrs, cnrs, false);


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

    shared_ptr<BaseSparseMatrix> bspmat;
    /** This does not work with gcc 7 **/
    // Switch<MAX_SYS_DIM>(dima-1, [&](auto DAMO) { // note: Switch<N> is [0..N-1], dim is [1..N]
    // 	constexpr int DIMA1 = DAMO + 1;
    // 	Switch<MAX_SYS_DIM>(dimb-1, [&](auto DBMO) {
    // 	    constexpr int DIMA = DIMA1; // otherwise gcc complains about captured value not being a constexpression
    // 	    constexpr int DIMB = DBMO + 1;
    // 	    if constexpr( (DIMA == DIMB) ||
    // 			  ( (DIMA == 1) || (DIMB == 1) ) ) { // so we do not try to instantiate SparseMatrix<Mat<2,3,SCAL>> etc.
    auto spmat = make_shared<TSPM<DIMB, DIMA, SCAL>>(std::move(graph));
    spmat->AsVector() = 0;
    bspmat = spmat;

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
      bool symmetric_so_far = true;
      for (auto bfi : ab_bfis)
	{ bfi->CalcElementMatrixAdd(felab, eltrans, bamat, symmetric_so_far, lh); }
      for (auto bfi : bb_bfis)
	{ bfi->CalcElementMatrixAdd(felb, eltrans, bbmat, symmetric_so_far, lh); }

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

      /** Set rownrs we do not want to write to -1 **/
      if (range_dofs != nullptr)
	for (auto & dnum : dnums_b)
	  if ( IsRegularDof(dnum) && !range_dofs->Test(dnum))
	    { dnum = -1; }

      /** Write into SparseMatrix - conversion to block-entries happens in AddElementMatrix, neat! **/
      spmat->AddElementMatrix(dnums_b, dnums_a, elmat, false); // should be taken care of by fespace coloring

      for (auto dnum : dnums_b) // space_b element coloring!!
	if (IsRegularDof(dnum))
	  { cnt_b[dnum]++; }
    };

    tass.Start();
    if (use_simd) {
      it_els([&](FESpace::Element & ei, LocalHeap & lh) {
	  try { fill_lam(ei, lh); }
	  catch (const ExceptionNOSIMD& e) { /** Turn off SIMD and continue **/
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
    tass.Stop();
#ifdef PARALLEL
    if (space_b->IsParallel() && !localop)
      { AllReduceDofData (cnt_b, MPI_SUM, space_b->GetParallelDofs()); }
#endif

    for (auto dofnr : Range(spmat->Height())) {
      if ( (cnt_b[dofnr] > 1) && ( !range_dofs || range_dofs->Test(dofnr) ) ) {
	double fac = 1.0 / double(cnt_b[dofnr]);
	auto rvs = spmat->GetRowValues(dofnr);
	for (auto & v : rvs)
	  { v *= fac; }
      }
    }
      // 	      } // (DIMA == 1) || (DIMB == 1)
      // 	  });
      // });

    shared_ptr<BaseMatrix> op = bspmat;

    if ( parmat && space_a->IsParallel() ) {
      op = make_shared<ParallelMatrix> (op, space_a->GetParallelDofs(), space_b->GetParallelDofs(),
					localop ? PARALLEL_OP::C2C : PARALLEL_OP::C2D);
    }

    return op;
  } // ConvertOperator



  template<class SCAL>
  shared_ptr<BaseMatrix> ConvertOperatorGF (shared_ptr<FESpace> space_a, shared_ptr<FESpace> space_b,
					  // int inda, int indb,
					  shared_ptr<DifferentialOperator> diffop, shared_ptr<CoefficientFunction> trial_cf,
					  VorB vb, const Region * reg, LocalHeap & lh,
					  shared_ptr<BitArray> range_dofs = nullptr,
					  bool localop = false, bool parmat = true, bool use_simd = true,
					  int bonus_intorder_ab = 0, int bonus_intorder_bb = 0)
  {

    static Timer t ("ConvertOperatorGF");
    RegionTimer regt(t);

    auto ma = space_b->GetMeshAccess();

    if ( parmat && (space_a->IsParallel() != space_b->IsParallel()) )
      { throw Exception("Cannot form ConvertOperator between a parallel and a local space!"); }

    if ( parmat && space_a->IsParallel() && space_b->IsParallel() ) {
      MPI_Comm comma = space_a->GetParallelDofs()->GetCommunicator(), commb = space_b->GetParallelDofs()->GetCommunicator();
      if (comma != commb)
	{ throw Exception("Cannot convert between spaces defined on different Communicators!"); }
    }

    /** Proxies and Integrators **/

    shared_ptr<CoefficientFunction> trial_a;
    if (trial_cf != nullptr)
      { trial_a = trial_cf; }
    else if ( diffop != nullptr ) {
      trial_a = make_shared<ProxyFunction>(space_a, false, false, diffop,
					   nullptr, nullptr, nullptr, nullptr, nullptr);
    }
    else { // Probably most of the time
      trial_a = make_shared<ProxyFunction>(space_a, false, false, space_a->GetEvaluator(vb),
					   nullptr, nullptr, nullptr, nullptr, nullptr);
    }

    auto trial_b = make_shared<ProxyFunction>(space_b, false, false, space_b->GetEvaluator(vb),
					      nullptr, nullptr, nullptr, nullptr, nullptr);

    if (!space_b->GetAdditionalEvaluators().Used("dual"))
      throw Exception(string("Dual diffop does not exist for ") + space_b->GetClassName() + string("!"));
    auto dual_evaluator = space_b->GetAdditionalEvaluators()["dual"];
    for (VorB avb = dual_evaluator->VB(); avb < vb; avb++) {
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
      ab_bfi->SetBonusIntegrationOrder(bonus_intorder_ab);
      ab_bfis.Append(ab_bfi);
      auto bb_bfi = make_shared<SymbolicBilinearFormIntegrator>(bb, vb, element_vb);
      bb_bfi->SetSimdEvaluate(use_simd);
      bb_bfis.Append(bb_bfi);
      bb_bfi->SetBonusIntegrationOrder(bonus_intorder_bb);
    }
 
    /** Create Matrix Graph **/
    int dima = space_a->GetDimension(), dimb = space_b->GetDimension();

    /** element equivalence classes **/
    // LocalHeap lh(10000000);

    /** #of equivalence classes per element type **/
    Array<short> et_eqc_cnt(ET_HEX+1); et_eqc_cnt = 0;
    ma->IterateElements
      (vb, lh, [&] (auto el, LocalHeap & llh) {
	short anum = 1 +
	  SwitchET (el.GetType(),
		    [&] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
	et_eqc_cnt[int(el.GetType())] = max(et_eqc_cnt[int(el.GetType())], anum);
      });
    /** prefix-sum so we have a unique nr for every element-type/element-eqc pair **/
    Array<short> et_firsti(ET_HEX+2); et_firsti = 0;
    for (auto k : Range(et_eqc_cnt))
      { et_firsti[k+1] = et_eqc_cnt[k] + et_firsti[k]; }
    /** actual class nrs **/
    Array<short> classnr(ma->GetNE(vb));
    ma->IterateElements
      (vb, lh, [&] (auto el, LocalHeap & llh) 
      {
        if ( (space_a->DefinedOn(vb, el.GetIndex())) && (space_b->DefinedOn(vb, el.GetIndex())) )
        {
            classnr[el.Nr()] =
            SwitchET
            (el.GetType(),
            [&] (auto et) { return et_firsti[int(et)] + ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
        }
        else
        {classnr[el.Nr()] = -1;}
      });
    
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        { if (classnr[i] != -1){creator.Add (classnr[i], i);} } //{ creator.Add (classnr[i], i); }
    Table<size_t> table = creator.MoveTable();

    /** assemble element matrix for every equivalence class **/
    Array<int> cnt_b(space_b->GetNDof()); cnt_b = 0;

    shared_ptr<BaseMatrix> op;

    auto simd_guard = [&](auto some_lam) {
      if (use_simd) {
	try { some_lam(); }
	catch (const ExceptionNOSIMD& e) { /** Turn off SIMD and continue **/
	  for (auto bfi : ab_bfis)
	    { bfi->SetSimdEvaluate(false); }
	  for (auto bfi : bb_bfis)
	    { bfi->SetSimdEvaluate(false); }
	  some_lam();
	}
      }
      else
	{ some_lam(); }
    };

    for (auto elclass_inds : table) {

      HeapReset hr(lh);

      if (elclass_inds.Size() == 0) continue;

      ElementId ei(vb, elclass_inds[0]);

      auto & eltrans = ma->GetTrafo(ei, lh);

      auto & fela = space_a->GetFE (ei, lh); int nda = fela.GetNDof();
      auto & felb = space_b->GetFE (ei, lh); int ndb = felb.GetNDof();
      MixedFiniteElement felab(fela, felb);

      FlatMatrix<SCAL> bamat(ndb*dimb, nda*dima, lh), bbmat(ndb*dimb, ndb*dimb, lh),
	elmat(ndb*dimb, nda*dima, lh);

      simd_guard([&]() {
	  bamat = 0.0; bbmat = 0.0;
	  bool symmetric_so_far = true; // will be set to false here
	  for (auto bfi : ab_bfis)
	    { bfi->CalcElementMatrixAdd(felab, eltrans, bamat, symmetric_so_far, lh); }
	  symmetric_so_far = true; // will probably be set to false
	  for (auto bfi : bb_bfis)
	    { bfi->CalcElementMatrixAdd(felb, eltrans, bbmat, symmetric_so_far, lh); }
	});

      CalcInverse(bbmat);
      elmat = bbmat * bamat;

      Table<DofId> adofs(elclass_inds.Size(), nda*dima),
	bdofs(elclass_inds.Size(), ndb*dimb);
      Array<DofId> dnumsa, dnumsb;
      for (auto i : Range(elclass_inds)) {
	ElementId ei(vb, elclass_inds[i]);
	space_a->GetDofNrs(ei, dnumsa);
	space_b->GetDofNrs(ei, dnumsb);
	for (auto d : dnumsb)
	  { cnt_b[d]++; }
	int c = 0;
	for (auto l : Range(dnumsa))
	  for (auto ll : Range(dima))
	    { adofs[i][c++] = dima * dnumsa[l] + ll; }
	c = 0;
	for (auto l : Range(dnumsb))
	  for (auto ll : Range(dimb))
	    { bdofs[i][c++] = dimb * dnumsb[l] + ll; }
      }

      auto mat = make_shared<ConstantElementByElementMatrix>
	(space_b->GetNDof(), space_a->GetNDof(),
	 elmat, std::move(bdofs), std::move(adofs));

      if (op != nullptr)
	{ op = make_shared<SumMatrix>(op, mat); }
      else
	{ op = mat; }
    }

    if (op == nullptr) {
      // dummy op for empty FESpace
      Table<DofId> xdofs(0, 0), ydofs(0, 0);
      Matrix<> mat(0,0);
      op = make_shared<ConstantElementByElementMatrix>
	(space_b->GetNDof(), space_a->GetNDof(),
	 mat, std::move(ydofs), std::move(xdofs));
    }

#ifdef PARALLEL
    if (space_b->IsParallel() && !localop)
      { AllReduceDofData (cnt_b, MPI_SUM, space_b->GetParallelDofs()); }
#endif

    bool multiple = false;
    for (auto c : cnt_b)
      if (c > 1)
	{ multiple = true; break; }

    if (multiple || range_dofs)
      {
        VVector<> scaling(cnt_b.Size()*dimb);
        for (auto dofnr : Range(cnt_b))
	  if ( range_dofs && !range_dofs->Test(dofnr) )
	    for (auto ll : Range(dimb))
	      { scaling(dimb*dofnr+ll) = 0.0; }
	  else
	    for (auto ll : Range(dimb))
	      { scaling(dimb*dofnr+ll) = cnt_b[dofnr] ? 1.0/cnt_b[dofnr] : 0.0; }
	auto diagmat = make_shared<DiagonalMatrix<>> (scaling);
        op = make_shared<ProductMatrix> (diagmat, op);
      }

    if ( parmat && space_a->IsParallel() ) {
      op = make_shared<ParallelMatrix> (op, space_a->GetParallelDofs(), space_b->GetParallelDofs(),
					localop ? PARALLEL_OP::C2C : PARALLEL_OP::C2D);
    }

    return op;
  } // ConvertOperatorGF



  shared_ptr<BaseMatrix> ConvertOperator (shared_ptr<FESpace> space_a, shared_ptr<FESpace> space_b, VorB vb, LocalHeap & lh,
					  shared_ptr<DifferentialOperator> diffop, shared_ptr<CoefficientFunction> trial_cf, const Region * reg,
					  shared_ptr<BitArray> range_dofs, bool localop, bool parmat, bool use_simd,
					  int bonus_intorder_ab, int bonus_intorder_bb, bool geom_free)
  {
    if ( space_a->IsComplex() != space_b->IsComplex() ) // b complex and a real could work in principle (?)
      { throw Exception("Cannot convert between complex and non-complex space!"); }

    shared_ptr<BaseMatrix> op;

    // This is a workaround because nested Switch does not work with gcc 7
    int dima = space_a->GetDimension();
    if (geom_free) {
      if (space_b->IsComplex())
	{ op = nullptr; throw Exception("ConstEBE not availaible for scalar type Complex!"); }
      else
	{ op = ConvertOperatorGF<double> (space_a, space_b, diffop, trial_cf, vb, reg, lh, range_dofs, localop, parmat, use_simd,
					  bonus_intorder_ab, bonus_intorder_bb); }
    }
    else { // geom_free
      Switch<MAX_SYS_DIM>(dima - 1, [&](auto DAMO) {
	  constexpr int DIMA = DAMO + 1;
	  int dimb = space_b->GetDimension();
	  if (dimb == 1) {
	    if (space_b->IsComplex())
	      { op = ConvertOperator<Complex, DIMA, 1> (space_a, space_b, diffop, trial_cf, vb, reg, lh, range_dofs, localop, parmat, use_simd,
							bonus_intorder_ab, bonus_intorder_bb); }
	    else
	      { op = ConvertOperator<double, DIMA, 1> (space_a, space_b, diffop, trial_cf, vb, reg, lh, range_dofs, localop, parmat, use_simd,
						       bonus_intorder_ab, bonus_intorder_bb); }
	  }
	  else if (dimb == dima) {
	    if (space_b->IsComplex())
	      { op = ConvertOperator<Complex, DIMA, DIMA> (space_a, space_b, diffop, trial_cf, vb, reg, lh, range_dofs, localop, parmat, use_simd,
							   bonus_intorder_ab, bonus_intorder_bb); }
	    else
	      { op = ConvertOperator<double, DIMA, DIMA> (space_a, space_b, diffop, trial_cf, vb, reg, lh, range_dofs, localop, parmat, use_simd,
							  bonus_intorder_ab, bonus_intorder_bb); }
	  }
	});
    } // geom_free

    return op;
  } // ConvertOperator

} // namespace ngcomp
