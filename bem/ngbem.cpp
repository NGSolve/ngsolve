#include "../comp/meshaccess.hpp"
#include "../comp/gridfunction.hpp"

#include "intrules_SauterSchwab.hpp"
#include "ngbem.hpp"
#include "fmmoperator.hpp"
#include "../linalg/elementbyelement.hpp"
#include "../linalg/diagonalmatrix.hpp"




namespace ngsbem
{


  IntOp_Parameters :: IntOp_Parameters (const Flags & flags)
  {
    auto use_fmm_flag = flags.GetDefineFlagX("use_fmm");
    if (use_fmm_flag.IsTrue()) use_fmm = true;
    if (use_fmm_flag.IsFalse()) use_fmm = false;

    fmm_maxdirect = int(flags.GetNumFlag("fmm_maxdirect", fmm_maxdirect));
    fmm_minorder = int(flags.GetNumFlag("fmm_minorder", fmm_minorder));
    fmm_order_factor = flags.GetNumFlag("fmm_order_factor", fmm_order_factor);
    fmm_separation = flags.GetNumFlag("fmm_separation", fmm_separation);
    fmm_eval_separation = flags.GetNumFlag("fmm_eval_separation", fmm_eval_separation);
    fmm_split_kr = flags.GetNumFlag("fmm_split_kr", fmm_split_kr);
    fmm_maxlevel = int(flags.GetNumFlag("fmm_maxlevel", fmm_maxlevel));
  }


  IntegralOperator ::
  IntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                   optional<Region> _trial_definedon, optional<Region> _test_definedon,
                   shared_ptr<DifferentialOperator> _trial_evaluator, // shared_ptr<CoefficientFunction> _trial_factor,
                   shared_ptr<DifferentialOperator> _test_evaluator, // shared_ptr<CoefficientFunction> _test_factor,
                   int _intorder, const IntOp_Parameters & _io_params,
                   VorB _trial_vb, VorB _test_vb)
    : trial_space(_trial_space), test_space(_test_space),
      trial_definedon(_trial_definedon), test_definedon(_test_definedon),
      trial_evaluator(_trial_evaluator), // trial_factor(_trial_factor),
      test_evaluator(_test_evaluator), // test_factor(_test_factor),
      trial_vb(_trial_vb), test_vb(_test_vb),
      intorder(_intorder), io_params(_io_params)
  {
    if (!test_space)
      test_space = trial_space;

    /*
    if (trial_factor)
      trial_evaluator = make_shared<DifferentialOperatorWithFactor>(trial_evaluator, trial_factor);

    if (test_factor)
      test_evaluator = make_shared<DifferentialOperatorWithFactor>(test_evaluator, test_factor);
    */
  }


  /*
  template <typename KERNEL>
  GenericIntegralOperator<KERNEL> ::
  GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                          optional<Region> _definedon_trial, optional<Region> _definedon_test,
                          shared_ptr<DifferentialOperator> _trial_evaluator, shared_ptr<CoefficientFunction> _trial_factor,
                          shared_ptr<DifferentialOperator> _test_evaluator, shared_ptr<CoefficientFunction> _test_factor,
                          KERNEL _kernel,
                          int _intorder)
  : IntegralOperator(_trial_space, _test_space, _definedon_trial, _definedon_test,
                     _trial_evaluator, _trial_factor,
                     _test_evaluator, _test_factor, _intorder), kernel(_kernel)
  {
    LocalHeap lh(100000000);

    matrix = this->CreateMatrixFMM(lh);
  }
  */


  template <typename KERNEL>
  GenericIntegralOperator<KERNEL> ::
  GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                          optional<Region> _definedon_trial, optional<Region> _definedon_test,
                          shared_ptr<DifferentialOperator> _trial_evaluator,
                          shared_ptr<DifferentialOperator> _test_evaluator,
                          KERNEL _kernel,
                          int _intorder, const IntOp_Parameters & _io_params,
                          VorB _trial_vb, VorB _test_vb)
  : IntegralOperator(_trial_space, _test_space, _definedon_trial, _definedon_test,
                     _trial_evaluator, _test_evaluator, _intorder, _io_params,
                     _trial_vb, _test_vb),
    kernel(_kernel)
  {
    if (trial_vb == BND && test_vb == BND)
      {
        tie(identic_panel_x, identic_panel_y, identic_panel_weight) = IdenticPanelIntegrationRule(intorder);
        tie(identic_panel_quad_x, identic_panel_quad_y, identic_panel_quad_weight) = IdenticPanelQuadIntegrationRule(intorder);

        tie(common_vertex_x, common_vertex_y, common_vertex_weight) = CommonVertexIntegrationRule(intorder);
        tie(common_vertex_quad_x, common_vertex_quad_y, common_vertex_quad_weight) = CommonVertexQuadIntegrationRule(intorder);
        tie(common_vertex_quadtrig_x, common_vertex_quadtrig_y, common_vertex_quadtrig_weight) = CommonVertexQuadTrigIntegrationRule(intorder);

        tie(common_edge_x, common_edge_y, common_edge_weight) = CommonEdgeIntegrationRule(intorder);
        tie(common_edge_quad_x, common_edge_quad_y, common_edge_quad_weight) = CommonEdgeQuadIntegrationRule(intorder);
        tie(common_edge_quadtrig_x, common_edge_quadtrig_y, common_edge_quadtrig_weight) = CommonEdgeQuadTrigIntegrationRule(intorder);
      }
    else if (trial_vb == VOL && test_vb == VOL)
      {
        tie(identic_tet_x, identic_tet_y, identic_tet_weight) = IdenticTetrahedronIntegrationRule(intorder);
        tie(common_face_tet_x, common_face_tet_y, common_face_tet_weight) = CommonFaceTetrahedronIntegrationRule(intorder);
        tie(common_edge_tet_x, common_edge_tet_y, common_edge_tet_weight) = CommonEdgeTetrahedronIntegrationRule(intorder);
        tie(common_vertex_tet_x, common_vertex_tet_y, common_vertex_tet_weight) = CommonVertexTetrahedronIntegrationRule(intorder);
      }
  }



  template <typename KERNEL>
  shared_ptr<BaseMatrix> GenericIntegralOperator<KERNEL> ::
  CreateMatrixFMM(LocalHeap & lh) const
  {
    static Timer tall("ngbem fmm setup"); RegionTimer r(tall);
    if (trial_vb != test_vb)
      throw Exception("BEM assembly for mixed source and test domains is not implemented");
    if (trial_vb != BND && trial_vb != VOL)
      throw Exception("BEM assembly supports boundary and volume domains only");

    Array<Vec<3>> xpts, ypts, xnv, ynv;
    IntegrationRule ir_trig(ET_TRIG, intorder);
    IntegrationRule ir_quad(ET_QUAD, intorder);
    IntegrationRule ir_tet(ET_TET, intorder);
    auto get_ir = [&] (ELEMENT_TYPE et) -> const IntegrationRule &
    {
      if (et == ET_TRIG) return ir_trig;
      if (et == ET_QUAD) return ir_quad;
      if (et == ET_TET) return ir_tet;
      throw Exception("BEM volume assembly currently supports tetrahedra only");
    };

    auto trial_mesh = trial_space->GetMeshAccess();
    auto test_mesh = test_space->GetMeshAccess();

    Array<int> compress_trial_els(trial_mesh->GetNE(trial_vb));
    Array<int> compress_test_els(test_mesh->GetNE(test_vb));
    compress_trial_els = -1;
    compress_test_els = -1;

    int cnt = 0;
    for (auto el : trial_mesh->Elements(trial_vb))
      if (trial_space->DefinedOn(el))
        if (!trial_definedon || (*trial_definedon).Mask().Test(trial_mesh->GetElIndex(el)))
          {
            HeapReset hr(lh);
            auto & trafo = trial_mesh->GetTrafo(el, lh);
            auto & ir = get_ir(trafo.GetElementType());
            auto & mir = trafo(ir, lh);
            for (auto & mip : mir)
              {
                xpts.Append(mip.GetPoint());
                if (trial_vb == BND)
                  xnv.Append(static_cast<const MappedIntegrationPoint<2,3>&>(mip).GetNV());
                else
                  xnv.Append(Vec<3>(0.0));
              }
            compress_trial_els[el.Nr()] = cnt++;
          }

    cnt = 0;
    for (auto el : test_mesh->Elements(test_vb))
      if (test_space->DefinedOn(el))
        if (!test_definedon || (*test_definedon).Mask().Test(test_mesh->GetElIndex(el)))
          {
            HeapReset hr(lh);
            auto & trafo = test_mesh->GetTrafo(el, lh);
            auto & ir = get_ir(trafo.GetElementType());
            auto & mir = trafo(ir, lh);
            for (auto & mip : mir)
              {
                ypts.Append(mip.GetPoint());
                if (test_vb == BND)
                  ynv.Append(static_cast<const MappedIntegrationPoint<2,3>&>(mip).GetNV());
                else
                  ynv.Append(Vec<3>(0.0));
              }
            compress_test_els[el.Nr()] = cnt++;
          }

    auto create_eval = [&](const FESpace & fes,
                           const Array<int> & compress_els,
                           const DifferentialOperator & evaluator,
                           VorB vb)
    {
      auto mesh = fes.GetMeshAccess();
      Array<short> classnr(mesh->GetNE(vb));
      mesh->IterateElements
        (vb, lh, [&] (auto el, LocalHeap & llh)
        {
          if (el.GetType() == ET_QUAD)
            classnr[el.Nr()] = 8 /* max trig calsses */  +  ET_trait<ET_QUAD>::GetClassNr(el.Vertices());
          else
            classnr[el.Nr()] =
              SwitchET<ET_SEGM, ET_TRIG, ET_TET>
              (el.GetType(),
               [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
        });

      TableCreator<size_t> creator;
      for ( ; !creator.Done(); creator++)
        for (auto i : Range(classnr))
          if (compress_els[i] != -1)
            creator.Add (classnr[i], i);
      Table<size_t> table = creator.MoveTable();

      // count integration points per compressed element
      int ncomp = 0;
      for (auto nr : compress_els) if (nr!=-1) ncomp++;
      Array<int> first_ip_nr(ncomp);
      int total_npts = 0;
      for (auto el : mesh->Elements(vb))
        if (compress_els[el.Nr()] != -1)
          {
            auto & ir = get_ir(el.GetType());
            first_ip_nr[compress_els[el.Nr()]] = total_npts;
            total_npts += ir.Size();
          }

      shared_ptr<BaseMatrix> evalx;

      for (auto elclass_inds : table)
        {
          if (elclass_inds.Size() == 0) continue;
          ElementId ei(vb, elclass_inds[0]);
          auto & felx = fes.GetFE (ei, lh);
          auto & ir = get_ir(mesh->GetElType(ei));

          int dim = evaluator.DimRef();
          Matrix<double,ColMajor> bmat_(dim*ir.Size(), felx.GetNDof());
          bmat_ = 0.0;

          // IntRange r1 = evaluator.UsedDofs(felx);   // TODO

          for (int i : Range(ir.Size()))
            evaluator.CalcMatrix(felx, ir[i], bmat_.Rows(dim*i, dim*(i+1)), lh);

          Matrix bmat = bmat_;

          Table<DofId> xdofsin(elclass_inds.Size(), felx.GetNDof());
          Table<DofId> xdofsout(elclass_inds.Size(), bmat.Height());

          Array<DofId> dnumsx, dnumsy;
          for (auto i : Range(elclass_inds))
            {
              ElementId ei(vb, elclass_inds[i]);
              fes.GetDofNrs(ei, dnumsx);
              xdofsin[i] = dnumsx;

              // for (int j = 0; j < dim*ir.Size(); j++)
              //   xdofsout[i][j] = compress_els[elclass_inds[i]]*(dim*ir.Size())+j;
              for (int j = 0; j < dim*ir.Size(); j++)
                xdofsout[i][j] = first_ip_nr[compress_els[elclass_inds[i]]]*dim+j;
            }

          auto part_evalx = make_shared<ConstantElementByElementMatrix<typename KERNEL::value_type>>
            // (mesh->GetNE(BND)*ir.Size()*dim, fes.GetNDof(),
            (total_npts*dim, fes.GetNDof(),
             bmat, std::move(xdofsout), std::move(xdofsin));

          if (evalx)
            evalx = evalx + part_evalx;
          else
            evalx = part_evalx;
        }

      // Tensor<3, typename KERNEL::value_type> weights(cnt*ir.Size(),
      //                                                evaluator.Dim(), evaluator.DimRef());
      Tensor<3, typename KERNEL::value_type> weights(total_npts,
                                                     evaluator.Dim(), evaluator.DimRef());
      Matrix<double> transformation(evaluator.Dim(), evaluator.DimRef());

      for (auto el : mesh->Elements(vb))
        if (compress_els[el.Nr()] != -1)
          {
            HeapReset hr(lh);
            auto & trafo = mesh->GetTrafo(el, lh);
            // auto & mir = trafo(ir, lh);
            auto & ir = get_ir(trafo.GetElementType());
            auto & mir = trafo(ir, lh);
            for (auto j : Range(mir.Size()))
              {
                evaluator.CalcTransformationMatrix(mir[j], transformation, lh);
                // weights(compress_els[el.Nr()]*mir.Size()+j,STAR,STAR) =
                weights(first_ip_nr[compress_els[el.Nr()]]+j,STAR,STAR) =
                  mir[j].GetWeight()*transformation;
              }
          }
      auto diagmat = make_shared<BlockDiagonalMatrix<typename KERNEL::value_type>>(std::move(weights));

      return diagmat*evalx;
    };

    auto evalx = create_eval(*trial_space, compress_trial_els, *trial_evaluator, trial_vb);
    auto evaly = create_eval(*test_space, compress_test_els, *test_evaluator, test_vb);
    auto fmmop = make_shared<FMM_Operator<KERNEL>> (kernel, std::move(xpts), std::move(ypts),
                                                    std::move(xnv), std::move(ynv), io_params);


    if (trial_mesh != test_mesh)
      return TransposeOperator(evaly) * fmmop * evalx;


    // **************   nearfield operator *****************

    static Timer tfind("ngbem fmm find near");
    static Timer tsetupgraph("ngbem fmm setup matrixgraph");
    static Timer tassemble("ngbem fmm assemble nearfield");
    static Timer tassSS("ngbem fmm assemble Sauter-Schwab");
    static Timer tasscorr("ngbem fmm assemble correction");

    tfind.Start();
    Array<tuple<size_t, size_t>> pairs;

    Array<size_t> other;
    for (ElementId ei : trial_mesh->Elements(trial_vb))
      if (trial_space->DefinedOn(ei))
        if (!trial_definedon || (*trial_definedon).Mask().Test(trial_mesh->GetElIndex(ei)))
        {
          other.SetSize0();
          for (auto v : trial_mesh->GetElement(ei).Vertices())
            for (auto ej : trial_mesh->GetVertexElements(v,trial_vb))
              if (test_space->DefinedOn(ElementId(test_vb,ej)))
                if (!test_definedon || (*test_definedon).Mask().Test(test_mesh->GetElIndex(ElementId(test_vb,ej))))
                  {
                    if (!other.Contains(ej))
                      {
                        other.Append (ej );
                        pairs.Append ( { ei.Nr(), ej });
                      }
                  }
        }


    if (!io_params.UseFMM())
      {
        pairs.SetSize0();
        for (ElementId ei : trial_mesh->Elements(trial_vb))
          if (trial_space->DefinedOn(ei))
            if (!trial_definedon || (*trial_definedon).Mask().Test(trial_mesh->GetElIndex(ei)))
              {
                for (ElementId ej : test_mesh->Elements(test_vb))
                  if (test_space->DefinedOn(ej))
                    if (!test_definedon || (*test_definedon).Mask().Test(test_mesh->GetElIndex(ej)))
                      pairs.Append ( { ei.Nr(), ej.Nr() });
              }
      }


    /*
    pairs.SetSize0();
    for (ElementId ei : trial_mesh->Elements(BND))
      for (ElementId ej : trial_mesh->Elements(BND))
        pairs.Append ( { ei.Nr(), ej.Nr() });
    */

    tfind.Stop();
    tsetupgraph.Start();
    TableCreator<int> trial_elements_creator(pairs.Size()), test_elements_creator(pairs.Size());

    for ( ; !trial_elements_creator.Done(); trial_elements_creator++, test_elements_creator++)
      {
        Array<DofId> dnums;
        for (auto i : Range(pairs))
          {
            auto [ei_trial, ei_test] = pairs[i];
            trial_space->GetDofNrs( { trial_vb, ei_trial }, dnums);
            trial_elements_creator.Add (i, dnums);

            test_space->GetDofNrs( { test_vb, ei_test }, dnums);
            test_elements_creator.Add (i, dnums);
          }
      }

    Table<int> trial_elements = trial_elements_creator.MoveTable();
    Table<int> test_elements = test_elements_creator.MoveTable();

    auto nearfield_correction =
      make_shared<SparseMatrix<value_type>> (test_space->GetNDof(), trial_space->GetNDof(),
                                             test_elements, trial_elements, false);
    tsetupgraph.Stop();
    tassemble.Start();
    nearfield_correction->SetZero();



    TableCreator<int> create_nbels(trial_mesh->GetNE(trial_vb));
    for ( ; !create_nbels.Done(); create_nbels++)
      for (auto i : Range(pairs))
        create_nbels.Add (get<0>(pairs[i]), get<1>(pairs[i]));
    Table<int> nbels = create_nbels.MoveTable();

    trial_mesh->IterateElements
      (trial_vb, lh, [&](auto ei_trial, LocalHeap & lh)
      {
        for (auto nrtest : nbels[ei_trial.Nr()])
          {
            HeapReset hr(lh);
            ElementId ei_test(test_vb, nrtest);

            auto & trial_trafo = trial_mesh -> GetTrafo(ei_trial, lh);
            auto & test_trafo = test_mesh -> GetTrafo(ei_test, lh);
            auto & trial_fel = trial_space->GetFE(ei_trial, lh);
            auto & test_fel = test_space->GetFE(ei_test, lh);

            Array<DofId> trial_dnums(trial_fel.GetNDof(), lh);
            Array<DofId> test_dnums(test_fel.GetNDof(), lh);

            trial_space->GetDofNrs (ei_trial, trial_dnums);
            test_space->GetDofNrs (ei_test, test_dnums);

            FlatMatrix<value_type> elmat(test_dnums.Size(), trial_dnums.Size(), lh);
            // tassSS.Start();
            CalcElementMatrix (elmat, ei_trial, ei_test, lh);
            // tassSS.Stop();
            // tasscorr.Start();

            // subtract the regular quadrature already represented by the FMM
            auto & ir_trial = get_ir(trial_trafo.GetElementType());
            auto & ir_test = get_ir(test_trafo.GetElementType());

            auto & trial_mir = trial_trafo(ir_trial, lh);
            auto & test_mir = test_trafo(ir_test, lh);

            // FlatMatrix<> shapesi(test_fel.GetNDof(), test_evaluator->Dim()*ir.Size(), lh);
            // FlatMatrix<> shapesj(trial_fel.GetNDof(), trial_evaluator->Dim()*ir.Size(), lh);
            FlatMatrix<> shapesi(test_fel.GetNDof(), test_evaluator->Dim()*ir_test.Size(), lh);
            FlatMatrix<> shapesj(trial_fel.GetNDof(), trial_evaluator->Dim()*ir_trial.Size(), lh);

            IntRange test_range = test_evaluator->UsedDofs(test_fel);
            IntRange trial_range = trial_evaluator->UsedDofs(trial_fel);

            test_evaluator -> CalcMatrix(test_fel, test_mir, Trans(shapesi), lh);
            trial_evaluator-> CalcMatrix(trial_fel, trial_mir, Trans(shapesj), lh);

            for (auto term : kernel.terms)
              {
                HeapReset hr(lh);
                // FlatMatrix<value_type> kernel_ixiy(ir.Size(), ir.Size(), lh);
                // for (int ix = 0; ix < ir.Size(); ix++)
                //   for (int iy = 0; iy < ir.Size(); iy++)
                FlatMatrix<value_type> kernel_ixiy(ir_test.Size(), ir_trial.Size(), lh);
                for (int ix = 0; ix < ir_test.Size(); ix++)
                  {
                    for (int iy = 0; iy < ir_trial.Size(); iy++)
                      {
                        Vec<3> x = test_mir[ix].GetPoint();
                        Vec<3> y = trial_mir[iy].GetPoint();

                        Vec<3> nx(0.0), ny(0.0);
                        if (test_vb == BND)
                          {
                            nx = static_cast<const MappedIntegrationPoint<2,3>&>(test_mir[ix]).GetNV();
                            ny = static_cast<const MappedIntegrationPoint<2,3>&>(trial_mir[iy]).GetNV();
                          }
                        value_type kernel_ = 0.0;
                        if (L2Norm2(x-y) > 0)
                          kernel_ = kernel.Evaluate(x, y, nx, ny)(term.kernel_comp);

                        double fac = test_mir[ix].GetWeight()*trial_mir[iy].GetWeight();
                        kernel_ixiy(ix, iy) = term.fac*fac*kernel_;
                      }
                  }

                // FlatMatrix<value_type> kernel_shapesj(ir.Size(), trial_fel.GetNDof(), lh);
                // FlatMatrix<> shapesi1(test_fel.GetNDof(), ir.Size(), lh);
                // FlatMatrix<> shapesj1(trial_fel.GetNDof(), ir.Size(), lh);
                FlatMatrix<value_type> kernel_shapesj(ir_test.Size(), trial_fel.GetNDof(), lh);
                FlatMatrix<> shapesi1(test_fel.GetNDof(), ir_test.Size(), lh);
                FlatMatrix<> shapesj1(trial_fel.GetNDof(), ir_trial.Size(), lh);

                // for (int j = 0; j < ir.Size(); j++)
                //   {
                //     shapesi1.Col(j) = shapesi.Col(test_evaluator->Dim()*j+term.test_comp);
                //     shapesj1.Col(j) = shapesj.Col(trial_evaluator->Dim()*j+term.trial_comp);
                //   }
                for (int j = 0; j < ir_test.Size(); j++)
                  shapesi1.Col(j) = shapesi.Col(test_evaluator->Dim()*j+term.test_comp);
                for (int j = 0; j < ir_trial.Size(); j++)
                  shapesj1.Col(j) = shapesj.Col(trial_evaluator->Dim()*j+term.trial_comp);

                kernel_shapesj = kernel_ixiy * Trans(shapesj1);

                if (io_params.UseFMM())
                  elmat.Rows(test_range).Cols(trial_range) -= shapesi1.Rows(test_range) * kernel_shapesj.Cols(trial_range);
              }
            // tasscorr.Stop();
            nearfield_correction -> AddElementMatrix (test_dnums, trial_dnums, elmat, true);
          }
      });


    tassemble.Stop();
    if (io_params.UseFMM())
      return TransposeOperator(evaly) * fmmop * evalx + nearfield_correction;
    else
      return nearfield_correction;
  }





  template <typename KERNEL>
  shared_ptr<BaseMatrix>
  GenericIntegralOperator<KERNEL> :: CreateNearFieldMatrix(LocalHeap & lh) const
  {
    if (trial_vb != test_vb)
      throw Exception("BEM assembly for mixed source and test domains is not implemented");

    auto trial_mesh = trial_space->GetMeshAccess();
    auto test_mesh = test_space->GetMeshAccess();


    Array<tuple<size_t, size_t>> pairs;

    Array<size_t> other;
    for (ElementId ei : trial_mesh->Elements(trial_vb))
      if (trial_space->DefinedOn(ei))
        if (!trial_definedon || (*trial_definedon).Mask().Test(trial_mesh->GetElIndex(ei)))
        {
          other.SetSize0();
          for (auto v : trial_mesh->GetElement(ei).Vertices())
            for (auto ej : trial_mesh->GetVertexElements(v,trial_vb))
              if (test_space->DefinedOn(ElementId(test_vb,ej)))
                if (!test_definedon || (*test_definedon).Mask().Test(test_mesh->GetElIndex(ElementId(test_vb,ej))))
                  {
                    if (!other.Contains(ej))
                      {
                        other.Append (ej );
                        pairs.Append ( { ei.Nr(), ej });
                      }
                  }
        }

    TableCreator<int> trial_elements_creator(pairs.Size()), test_elements_creator(pairs.Size());

    for ( ; !trial_elements_creator.Done(); trial_elements_creator++, test_elements_creator++)
      {
        Array<DofId> dnums;
        for (auto i : Range(pairs))
          {
            auto [ei_trial, ei_test] = pairs[i];
            trial_space->GetDofNrs( { trial_vb, ei_trial }, dnums);
            trial_elements_creator.Add (i, dnums);

            test_space->GetDofNrs( { test_vb, ei_test }, dnums);
            test_elements_creator.Add (i, dnums);
          }
      }

    Table<int> trial_elements = trial_elements_creator.MoveTable();
    Table<int> test_elements = test_elements_creator.MoveTable();

    auto nearfield =
      make_shared<SparseMatrix<value_type>> (test_space->GetNDof(), trial_space->GetNDof(),
                                             test_elements, trial_elements, false);

    nearfield->SetZero();

    TableCreator<int> create_nbels(trial_mesh->GetNE(trial_vb));
    for ( ; !create_nbels.Done(); create_nbels++)
      for (auto i : Range(pairs))
        create_nbels.Add (get<0>(pairs[i]), get<1>(pairs[i]));
    Table<int> nbels = create_nbels.MoveTable();

    trial_mesh->IterateElements
      (trial_vb, lh, [&](auto ei_trial, LocalHeap & lh)
      {
        for (auto nrtest : nbels[ei_trial.Nr()])
          {
            HeapReset hr(lh);
            ElementId ei_test(test_vb, nrtest);

            // auto & trial_trafo = trial_mesh -> GetTrafo(ei_trial, lh);
            // auto & test_trafo = test_mesh -> GetTrafo(ei_test, lh);
            auto & trial_fel = trial_space->GetFE(ei_trial, lh);
            auto & test_fel = test_space->GetFE(ei_test, lh);

            Array<DofId> trial_dnums(trial_fel.GetNDof(), lh);
            Array<DofId> test_dnums(test_fel.GetNDof(), lh);

            trial_space->GetDofNrs (ei_trial, trial_dnums);
            test_space->GetDofNrs (ei_test, test_dnums);

            FlatMatrix<value_type> elmat(test_dnums.Size(), trial_dnums.Size(), lh);
            CalcElementMatrix (elmat, ei_trial, ei_test, lh);

            /*
            MappedIntegrationRule<2,3> trial_mir(ir, trial_trafo, lh);
            MappedIntegrationRule<2,3> test_mir(ir, test_trafo, lh);

            FlatMatrix<> shapesi(test_fel.GetNDof(), test_evaluator->Dim()*ir.Size(), lh);
            FlatMatrix<> shapesj(trial_fel.GetNDof(), trial_evaluator->Dim()*ir.Size(), lh);

            IntRange test_range = test_evaluator->UsedDofs(test_fel);
            IntRange trial_range = trial_evaluator->UsedDofs(trial_fel);

            test_evaluator -> CalcMatrix(test_fel, test_mir, Trans(shapesi), lh);
            trial_evaluator-> CalcMatrix(trial_fel, trial_mir, Trans(shapesj), lh);
            */

            nearfield -> AddElementMatrix (test_dnums, trial_dnums, elmat, true);
          }
      });

    return nearfield;
  }


  template <typename KERNEL>
  void GenericIntegralOperator<KERNEL> ::
  CalcElementMatrix(FlatMatrix<value_type> matrix,
                    ElementId ei_trial, ElementId ei_test,
                    LocalHeap &lh) const
  {
    auto mesh = this->trial_space->GetMeshAccess();
    auto mesh2 = this->test_space->GetMeshAccess();

    // static Timer tall("ngbem - elementmatrix " + KERNEL::Name());
    // RegionTimer reg(tall);

    // static Timer t1("ngbem - elementmatrix, part1  " + KERNEL::Name());

    // t1.Start();
    matrix = 0.;
    // t1.Stop();

    Vec<3> x,y,nx,ny;
    typedef decltype(kernel.Evaluate (x,y,nx,ny)) KERNEL_COMPS_T;

    if ((trial_vb == VOL) != (test_vb == VOL))
      throw Exception("mixed boundary-volume element matrices are not implemented");

    if (trial_vb == VOL)
      {
        auto verti = mesh2->GetElement(ei_test).Vertices();
        auto vertj = mesh->GetElement(ei_trial).Vertices();

        FiniteElement & feli = test_space->GetFE(ei_test, lh);
        FiniteElement & felj = trial_space->GetFE(ei_trial, lh);
        ElementTransformation & trafoi = mesh2->GetTrafo(ei_test, lh);
        ElementTransformation & trafoj = mesh->GetTrafo(ei_trial, lh);

        auto Integrate6D = [&] (const IntegrationRule & irx,
                                const IntegrationRule & iry,
                                const FiniteElement & feli,
                                const FiniteElement & felj,
                                const ElementTransformation & trafoi,
                                const ElementTransformation & trafoj,
                                FlatMatrix<value_type> elmat,
                                LocalHeap & lh)
        {
          HeapReset hr(lh);
          SIMD_IntegrationRule simd_irx(irx);
          SIMD_IntegrationRule simd_iry(iry);

          SIMD_MappedIntegrationRule<3,3> mirx(simd_irx, trafoi, lh);
          SIMD_MappedIntegrationRule<3,3> miry(simd_iry, trafoj, lh);

          FlatMatrix<SIMD<double>> mshapesi(feli.GetNDof()*test_evaluator->Dim(), mirx.Size(), lh);
          FlatMatrix<SIMD<value_type>> mshapesi_kern(feli.GetNDof(), mirx.Size(), lh);
          FlatMatrix<SIMD<double>> mshapesj(felj.GetNDof()*trial_evaluator->Dim(), miry.Size(), lh);

          IntRange test_range = test_evaluator->UsedDofs(feli);
          IntRange trial_range = trial_evaluator->UsedDofs(felj);
          test_evaluator->CalcMatrix(feli, mirx, mshapesi);
          trial_evaluator->CalcMatrix(felj, miry, mshapesj);

          FlatVector<Vec<KERNEL_COMPS_T::SIZE, SIMD<value_type>>> kernel_values(mirx.Size(), lh);
          Vec<3,SIMD<double>> zero(0.0);
          for (int k2 = 0; k2 < mirx.Size(); k2++)
            kernel_values(k2) = mirx[k2].GetMeasure()*miry[k2].GetMeasure()
              * simd_irx[k2].Weight()
              * kernel.Evaluate(mirx[k2].Point(), miry[k2].Point(), zero, zero);

          for (auto term : kernel.terms)
            {
              auto mshapesi_comp = mshapesi.RowSlice(term.test_comp, test_evaluator->Dim());
              for (int k2 = 0; k2 < mirx.Size(); k2++)
                mshapesi_kern.Col(k2) = term.fac*kernel_values(k2)(term.kernel_comp)
                  * mshapesi_comp.Col(k2);

              AddABt(mshapesi_kern.Rows(test_range),
                     mshapesj.RowSlice(term.trial_comp, trial_evaluator->Dim()).AddSize(felj.GetNDof(), miry.Size()).Rows(trial_range),
                     elmat.Rows(test_range).Cols(trial_range));
            }
        };

        if (feli.ElementType() != ET_TET || felj.ElementType() != ET_TET)
          throw Exception("singular volume BEM assembly currently supports tetrahedra only");

        Array<IVec<2>> common;
        for (int i = 0; i < verti.Size(); i++)
          for (int j = 0; j < vertj.Size(); j++)
            if (verti[i] == vertj[j])
              common.Append(IVec<2>(i,j));

        if (common.Size() == 0)
          {
            IntegrationRule iri(ET_TET, intorder);
            IntegrationRule irj(ET_TET, intorder);
            MappedIntegrationRule<3,3> mirx(iri, trafoi, lh);
            MappedIntegrationRule<3,3> miry(irj, trafoj, lh);

            FlatMatrix<> shapesi(feli.GetNDof(), test_evaluator->Dim()*iri.Size(), lh);
            FlatMatrix<> shapesj(felj.GetNDof(), trial_evaluator->Dim()*irj.Size(), lh);
            shapesi = 0.0;
            shapesj = 0.0;
            test_evaluator->CalcMatrix(feli, mirx, Trans(shapesi), lh);
            trial_evaluator->CalcMatrix(felj, miry, Trans(shapesj), lh);

            Vec<3> zero(0.0);
            for (auto term : kernel.terms)
              {
                HeapReset hr(lh);
                FlatMatrix<value_type> kernel_ixiy(iri.Size(), irj.Size(), lh);
                for (int ix = 0; ix < iri.Size(); ix++)
                  for (int iy = 0; iy < irj.Size(); iy++)
                    kernel_ixiy(ix,iy) = term.fac * mirx[ix].GetWeight() * miry[iy].GetWeight()
                      * kernel.Evaluate(mirx[ix].GetPoint(), miry[iy].GetPoint(), zero, zero)(term.kernel_comp);

                FlatMatrix<value_type> kernel_shapesj(iri.Size(), felj.GetNDof(), lh);
                FlatMatrix<> shapesi1(feli.GetNDof(), iri.Size(), lh);
                FlatMatrix<> shapesj1(felj.GetNDof(), irj.Size(), lh);
                for (int j = 0; j < iri.Size(); j++)
                  shapesi1.Col(j) = shapesi.Col(test_evaluator->Dim()*j+term.test_comp);
                for (int j = 0; j < irj.Size(); j++)
                  shapesj1.Col(j) = shapesj.Col(trial_evaluator->Dim()*j+term.trial_comp);

                kernel_shapesj = kernel_ixiy * Trans(shapesj1);
                matrix += shapesi1 * kernel_shapesj;
              }
            return;
          }

        const Array<Vec<3>> * quad_x = nullptr;
        const Array<Vec<3>> * quad_y = nullptr;
        const Array<double> * quad_weight = nullptr;
        Array<int> permi, permj;

        switch (common.Size())
          {
          case 4:
            quad_x = &identic_tet_x;
            quad_y = &identic_tet_y;
            quad_weight = &identic_tet_weight;
            break;
          case 3:
            quad_x = &common_face_tet_x;
            quad_y = &common_face_tet_y;
            quad_weight = &common_face_tet_weight;
            break;
          case 2:
            quad_x = &common_edge_tet_x;
            quad_y = &common_edge_tet_y;
            quad_weight = &common_edge_tet_weight;
            break;
          case 1:
            quad_x = &common_vertex_tet_x;
            quad_y = &common_vertex_tet_y;
            quad_weight = &common_vertex_tet_weight;
            break;
          default:
            throw Exception("invalid tetrahedron-pair topology");
          }

        for (auto pair : common)
          {
            permi.Append(pair[0]);
            permj.Append(pair[1]);
          }
        for (int nr = 0; nr < 4; nr++)
          {
            if (!permi.Contains(nr)) permi.Append(nr);
            if (!permj.Contains(nr)) permj.Append(nr);
          }

        auto PermutePoint = [] (Vec<3> p, const Array<int> & perm)
        {
          Vec<4> lambda(p[0], p[1], p[2], 1-p[0]-p[1]-p[2]);
          Vec<4> mapped;
          for (int i = 0; i < 4; i++)
            mapped[perm[i]] = lambda[i];
          return Vec<3>(mapped[0], mapped[1], mapped[2]);
        };

        constexpr int BS = 128;
        for (int k = 0; k < quad_weight->Size(); k += BS)
          {
            int num = std::min(size_t(BS), quad_weight->Size()-k);
            HeapReset hr(lh);
            IntegrationRule irx(num, lh), iry(num, lh);
            for (int k2 = 0; k2 < num; k2++)
              {
                Vec<3> px = PermutePoint((*quad_x)[k+k2], permi);
                Vec<3> py = PermutePoint((*quad_y)[k+k2], permj);
                irx[k2] = IntegrationPoint(px[0], px[1], px[2], (*quad_weight)[k+k2]);
                iry[k2] = IntegrationPoint(py[0], py[1], py[2], 0.0);
              }
            Integrate6D(irx, iry, feli, felj, trafoi, trafoj, matrix, lh);
          }
        return;
      }

    // common code for same panel, common edge, common vertex
    auto Integrate4D = [&] (const IntegrationRule & irx,
                            const IntegrationRule & iry,
                            const FiniteElement & feli,
                            const FiniteElement & felj,
                            const ElementTransformation & trafoi,
                            const ElementTransformation & trafoj,
                            FlatMatrix<value_type> elmat,
                            LocalHeap & lh)
    {
      HeapReset hr(lh);
      SIMD_IntegrationRule simd_irx(irx);
      SIMD_IntegrationRule simd_iry(iry);

      SIMD_MappedIntegrationRule<2,3> mirx(simd_irx, trafoi, lh);
      SIMD_MappedIntegrationRule<2,3> miry(simd_iry, trafoj, lh);

      FlatMatrix<SIMD<double>> mshapesi(feli.GetNDof()*test_evaluator->Dim(), mirx.Size(), lh);
      FlatMatrix<SIMD<value_type>> mshapesi_kern(feli.GetNDof(), mirx.Size(), lh);
      FlatMatrix<SIMD<double>> mshapesj(felj.GetNDof()*trial_evaluator->Dim(), miry.Size(), lh);

      // mshapesi = 0.;
      // mshapesj = 0.;

      IntRange test_range = test_evaluator->UsedDofs(feli);
      IntRange trial_range = trial_evaluator->UsedDofs(felj);
      test_evaluator->CalcMatrix(feli, mirx, mshapesi);  // only used are set for compound fe !!!
      trial_evaluator->CalcMatrix(felj, miry, mshapesj);

      FlatVector<Vec<KERNEL_COMPS_T::SIZE, SIMD<value_type>>> kernel_values(mirx.Size(), lh);
      for (int k2 = 0; k2 < mirx.Size(); k2++)
        {
          Vec<3,SIMD<double>> x = mirx[k2].Point();
          Vec<3,SIMD<double>> y = miry[k2].Point();
          Vec<3,SIMD<double>> nx = mirx[k2].GetNV();
          Vec<3,SIMD<double>> ny = miry[k2].GetNV();
          kernel_values(k2) = mirx[k2].GetMeasure()*miry[k2].GetMeasure()*simd_irx[k2].Weight() *
            kernel.Evaluate(x, y, nx, ny);
        }
      for (auto term : kernel.terms)
        {
          auto mshapesi_comp = mshapesi.RowSlice(term.test_comp, test_evaluator->Dim());
          for (int k2 = 0; k2 < mirx.Size(); k2++)
            {
              SIMD<value_type> kernel_ = kernel_values(k2)(term.kernel_comp);
              mshapesi_kern.Col(k2) = term.fac*kernel_ * mshapesi_comp.Col(k2);
            }

          AddABt (mshapesi_kern.Rows(test_range),
                  mshapesj.RowSlice(term.trial_comp, trial_evaluator->Dim()).AddSize(felj.GetNDof(), miry.Size()).Rows(trial_range),
                  elmat.Rows(test_range).Cols(trial_range));
        }
    };


    auto Integrate4DMapped = [&] (FlatArray<Vec<2>> quad_x, FlatArray<Vec<2>> quad_y, FlatArray<double> weight,
                                  Vec<2> p0x, Mat<2,2> Tx, Vec<2> p0y, Mat<2,2> Ty,
                                  const FiniteElement & feli,
                                  const FiniteElement & felj,
                                  const ElementTransformation & trafoi,
                                  const ElementTransformation & trafoj,
                                  FlatMatrix<value_type> elmat,
                                  LocalHeap & lh)
    {
      constexpr int BS = 128;
      for (int k = 0; k < weight.Size(); k+=BS)
        {
          int num = std::min(size_t(BS), weight.Size()-k);

          HeapReset hr(lh);

          IntegrationRule irx(num, lh);
          IntegrationRule iry(num, lh);

          for (int k2 = 0; k2 < num; k2++)
            {
              Vec<2> px = p0x + Tx * quad_x[k+k2];
              Vec<2> py = p0y + Ty * quad_y[k+k2];

              irx[k2] = IntegrationPoint(px(0), px(1), 0, weight[k+k2]);
              iry[k2] = IntegrationPoint(py(0), py(1), 0, 0);
            }

          Integrate4D (irx, iry, feli, felj, trafoi, trafoj, elmat, lh);
        }
    };




    auto verti = mesh2->GetElement(ei_test).Vertices();
    auto vertj = mesh->GetElement(ei_trial).Vertices();

    FiniteElement &feli = test_space->GetFE(ei_test, lh);
    FiniteElement &felj = trial_space->GetFE(ei_trial, lh);

    ElementTransformation &trafoi = mesh2->GetTrafo(ei_test, lh);
    ElementTransformation &trafoj = mesh->GetTrafo(ei_trial, lh);

    Array<DofId> dnumsi, dnumsj;
    test_space->GetDofNrs(ei_test, dnumsi); // mapping to global dof
    trial_space->GetDofNrs(ei_trial, dnumsj);

    FlatMatrix<> shapei(feli.GetNDof(), test_evaluator->Dim(), lh);
    FlatMatrix<> shapej(felj.GetNDof(), trial_evaluator->Dim(), lh);

    matrix = 0.0;

    int n_common_vertices = 0;
    for (auto vi : verti)
      if (vertj.Contains(vi))
        n_common_vertices++;

    /*
    // treat quad-quad and quad-trig as disjoint
    if ((verti.Size()==4) != (vertj.Size()==4))
      if (n_common_vertices == 1)
        n_common_vertices = 0;
    */


    switch (n_common_vertices)
      {
      case 4: //identical panel quad
        {
          constexpr int BS = 128;
          for (int k = 0; k < identic_panel_quad_weight.Size(); k+=BS)
            {
              int num = std::min(size_t(BS), identic_panel_quad_weight.Size()-k);

              HeapReset hr(lh);

              IntegrationRule irx(num, lh);
              IntegrationRule iry(num, lh);

              for (int k2 = 0; k2 < num; k2++)
                {
                  Vec<2> xk = identic_panel_quad_x[k+k2];
                  Vec<2> yk = identic_panel_quad_y[k+k2];

                  irx[k2] = IntegrationPoint(xk(0), xk(1), 0,
                                             identic_panel_quad_weight[k+k2]);
                  iry[k2] = IntegrationPoint(yk(0), yk(1), 0, 0);
                }

              Integrate4D (irx, iry, feli, felj, trafoi, trafoj, matrix, lh);
            }


          break;
        }

      case 3: //identical panel trig
        {
          constexpr int BS = 128;
          for (int k = 0; k < identic_panel_weight.Size(); k+=BS)
            {
              int num = std::min(size_t(BS), identic_panel_weight.Size()-k);

              HeapReset hr(lh);

              IntegrationRule irx(num, lh);
              IntegrationRule iry(num, lh);

              for (int k2 = 0; k2 < num; k2++)
                {
                  Vec<2> xk = identic_panel_x[k+k2];
                  Vec<2> yk = identic_panel_y[k+k2];

                  irx[k2] = IntegrationPoint(xk(0), xk(1), 0,
                                             identic_panel_weight[k+k2]);
                  iry[k2] = IntegrationPoint(yk(0), yk(1), 0, 0);
                }

              Integrate4D (irx, iry, feli, felj, trafoi, trafoj, matrix, lh);
            }

          break;
        }


      case 2: //common edge
        {
          // RegionTimer reg(t_common_edge);

          //FlatArray<const EDGE> edgesx(ElementTopology::GetNEdges(feli.ElementType()),  ElementTopology::GetEdges (feli.ElementType()));    // 0 1 | 1 2 | 2 0
          // FlatArray<const EDGE> edgesy(ElementTopology::GetNEdges(felj.ElementType()),  ElementTopology::GetEdges (felj.ElementType()));    // 0 1 | 1 2 | 2 0

          auto edgesx = ElementTopology::GetEdges(feli.ElementType());
          auto edgesy = ElementTopology::GetEdges(felj.ElementType());

          int cex=0, cey=0;
          bool same_orientation=false;
          for (int cx = 0; cx < edgesx.Size(); cx++)
            for (int cy = 0; cy < edgesy.Size(); cy++)
              {
                IVec<2> ex (verti[edgesx[cx][0]], verti[edgesx[cx][1]]);
                IVec<2> ey (vertj[edgesy[cy][0]], vertj[edgesy[cy][1]]);
                bool same = ex==ey;
                if (ex.Sort() == ey.Sort())
                  {
                    cex = cx;  // -> "common" edge number element i
                    cey = cy;  // -> "common" edge number element j
                    same_orientation = same;
                    break;
                  }
              }


          Vec<2> trigverts[3] = { { 1,0 }, { 0, 1}, { 0, 0 } };
          Vec<2> quadverts[4] = { { 0,0 }, { 1, 0}, { 1, 1 }, { 0,1}  };
          Vec<2> gradverts[4] = { { -1,-1 }, { 1, -1}, { 1, 1 }, { -1,1}  };


          Vec<2> p0x, p0y;
          Mat<2,2> Tx, Ty;

          if (edgesx.Size() == 3)
            {
              p0x = trigverts[edgesx[cex][0]];
              Tx.Col(0) = trigverts[edgesx[cex][1]] - p0x;
              Tx.Col(1) = trigverts[3-edgesx[cex][0]-edgesx[cex][1]] - p0x;
            }
          else
            {
              p0x = quadverts[edgesx[cex][0]];
              Tx.Col(0) = quadverts[edgesx[cex][1]]-p0x;
              Tx.Col(1) = -0.5*(gradverts[edgesx[cex][0]]+gradverts[edgesx[cex][1]]);
            }

          if (edgesy.Size() == 3)
            {
              if (same_orientation)
                {
                  p0y = trigverts[edgesy[cey][0]];
                  Ty.Col(0) = trigverts[edgesy[cey][1]] - p0y;
                  Ty.Col(1) = trigverts[3-edgesy[cey][0]-edgesy[cey][1]] - p0y;
                }
              else
                {
                  p0y = trigverts[edgesy[cey][1]];
                  Ty.Col(0) = trigverts[edgesy[cey][0]] - p0y;
                  Ty.Col(1) = trigverts[3-edgesy[cey][0]-edgesy[cey][1]] - p0y;
                }
            }
          else
            {
              if (same_orientation)
                {
                  p0y = quadverts[edgesy[cey][0]];
                  Ty.Col(0) = quadverts[edgesy[cey][1]]-p0y;
                  Ty.Col(1) = -0.5*(gradverts[edgesy[cey][0]]+gradverts[edgesy[cey][1]]);
                }
              else
                {
                  p0y = quadverts[edgesy[cey][1]];
                  Ty.Col(0) = quadverts[edgesy[cey][0]]-p0y;
                  Ty.Col(1) = -0.5*(gradverts[edgesy[cey][0]]+gradverts[edgesy[cey][1]]);
                }
            }


          if (edgesx.Size()==3 && edgesy.Size()==3)
            {
              Integrate4DMapped (common_edge_x, common_edge_y, common_edge_weight,
                                 p0x, Tx, p0y, Ty, feli, felj, trafoi, trafoj, matrix, lh);
              /*
              constexpr int BS = 128;
              for (int k = 0; k < common_edge_weight.Size(); k+=BS)
                {
                  int num = std::min(size_t(BS), common_edge_weight.Size()-k);

                  HeapReset hr(lh);

                  IntegrationRule irx(num, lh);
                  IntegrationRule iry(num, lh);

                  for (int k2 = 0; k2 < num; k2++)
                    {
                      Vec<2> px = p0x + Tx * common_edge_x[k+k2];
                      Vec<2> py = p0y + Ty * common_edge_y[k+k2];

                      irx[k2] = IntegrationPoint(px(0), px(1), 0, common_edge_weight[k+k2]);
                      iry[k2] = IntegrationPoint(py(0), py(1), 0, 0);
                    }

                  Integrate4D (irx, iry, feli, felj, trafoi, trafoj, matrix, lh);
                }
              */
            }

          else if (edgesx.Size()==4 && edgesy.Size()==3)
            {
              constexpr int BS = 128;
              for (int k = 0; k < common_edge_quadtrig_weight.Size(); k+=BS)
                {
                  int num = std::min(size_t(BS), common_edge_quadtrig_weight.Size()-k);

                  HeapReset hr(lh);

                  IntegrationRule irx(num, lh);
                  IntegrationRule iry(num, lh);

                  for (int k2 = 0; k2 < num; k2++)
                    {
                      Vec<2> px = p0x + Tx * common_edge_quadtrig_x[k+k2];
                      Vec<2> py = p0y + Ty * common_edge_quadtrig_y[k+k2];

                      irx[k2] = IntegrationPoint(px(0), px(1), 0, common_edge_quadtrig_weight[k+k2]);
                      iry[k2] = IntegrationPoint(py(0), py(1), 0, 0);
                    }

                  Integrate4D (irx, iry, feli, felj, trafoi, trafoj, matrix, lh);
                }
            }

          else if (edgesx.Size()==3 && edgesy.Size()==4)
            {
              constexpr int BS = 128;
              for (int k = 0; k < common_edge_quadtrig_weight.Size(); k+=BS)
                {
                  int num = std::min(size_t(BS), common_edge_quadtrig_weight.Size()-k);

                  HeapReset hr(lh);

                  IntegrationRule irx(num, lh);
                  IntegrationRule iry(num, lh);

                  for (int k2 = 0; k2 < num; k2++)
                    {
                      Vec<2> px = p0x + Tx * common_edge_quadtrig_y[k+k2];
                      Vec<2> py = p0y + Ty * common_edge_quadtrig_x[k+k2];

                      irx[k2] = IntegrationPoint(px(0), px(1), 0, common_edge_quadtrig_weight[k+k2]);
                      iry[k2] = IntegrationPoint(py(0), py(1), 0, 0);
                    }

                  Integrate4D (irx, iry, feli, felj, trafoi, trafoj, matrix, lh);
                }
            }



          else if (edgesx.Size()==4 && edgesy.Size()==4)
            {
              Integrate4DMapped (common_edge_quad_x, common_edge_quad_y, common_edge_quad_weight,
                                 p0x, Tx, p0y, Ty, feli, felj, trafoi, trafoj, matrix, lh);
              /*
              constexpr int BS = 128;
              for (int k = 0; k < common_edge_quad_weight.Size(); k+=BS)
                {
                  int num = std::min(size_t(BS), common_edge_quad_weight.Size()-k);

                  HeapReset hr(lh);

                  IntegrationRule irx(num, lh);
                  IntegrationRule iry(num, lh);

                  for (int k2 = 0; k2 < num; k2++)
                    {
                      Vec<2> px = p0x + Tx * common_edge_quad_x[k+k2];
                      Vec<2> py = p0y + Ty * common_edge_quad_y[k+k2];

                      irx[k2] = IntegrationPoint(px(0), px(1), 0, common_edge_quad_weight[k+k2]);
                      iry[k2] = IntegrationPoint(py(0), py(1), 0, 0);
                    }

                  Integrate4D (irx, iry, feli, felj, trafoi, trafoj, matrix, lh);
                }
              */
            }

          break;
        }




      case 1: //common vertex
        {
          // RegionTimer reg(t_common_vertex);

          int cvx=-1, cvy=-1;
          for (int cx = 0; cx < verti.Size(); cx++)
            for (int cy = 0; cy < vertj.Size(); cy++)
              {
                if (verti[cx] == vertj[cy])
                  {
                    cvx = cx;
                    cvy = cy;
                    break;
                  }
              }



          Vec<2> trigverts[3] = { { 1,0 }, { 0, 1}, { 0, 0 } };
          Vec<2> quadverts[4] = { { 0,0 }, { 1, 0}, { 1, 1 }, { 0,1 } };


          Vec<2> p0x, p0y;
          Mat<2,2> Tx, Ty;

          if (verti.Size() == 3)
            {
              p0x = trigverts[cvx];
              Tx.Col(0) = trigverts[(cvx+1)%3] - p0x;
              Tx.Col(1) = trigverts[(cvx+2)%3] - p0x;
            }
          else
            {
              p0x = quadverts[cvx];
              Tx.Col(0) = quadverts[(cvx+1)%4] - p0x;
              Tx.Col(1) = quadverts[(cvx+3)%4] - p0x;
            }

          if (vertj.Size() == 3)
            {
              p0y = trigverts[cvy];
              Ty.Col(0) = trigverts[(cvy+1)%3] - p0y;
              Ty.Col(1) = trigverts[(cvy+2)%3] - p0y;
            }
          else
            {
              p0y = quadverts[cvy];
              Ty.Col(0) = quadverts[(cvy+1)%4] - p0y;
              Ty.Col(1) = quadverts[(cvy+3)%4] - p0y;
            }

          if (verti.Size() == 3 && vertj.Size() == 3)
            {
              Integrate4DMapped (common_vertex_x, common_vertex_y, common_vertex_weight,
                                 p0x, Tx, p0y, Ty, feli, felj, trafoi, trafoj, matrix, lh);
            }
          else if (verti.Size() == 4 && vertj.Size() == 3)
            {
              Integrate4DMapped (common_vertex_quadtrig_x, common_vertex_quadtrig_y, common_vertex_quadtrig_weight,
                                 p0x, Tx, p0y, Ty, feli, felj, trafoi, trafoj, matrix, lh);
            }
          else if (verti.Size() == 3 && vertj.Size() == 4)
            {
              Integrate4DMapped (common_vertex_quadtrig_y, common_vertex_quadtrig_x, common_vertex_quadtrig_weight,
                                 p0x, Tx, p0y, Ty, feli, felj, trafoi, trafoj, matrix, lh);
            }
          else if (verti.Size() == 4 && vertj.Size() == 4)
            {
              Integrate4DMapped (common_vertex_quad_x, common_vertex_quad_y, common_vertex_quad_weight,
                                 p0x, Tx, p0y, Ty, feli, felj, trafoi, trafoj, matrix, lh);
            }


          break;
        }

      case 0: //disjoint panels
        {
          // RegionTimer r(t_disjoint);
          // shapes+geom out of loop, matrix multiplication
          // IntegrationRule irtrig(ET_TRIG, intorder);
          IntegrationRule iri(feli.ElementType(), intorder);
          IntegrationRule irj(felj.ElementType(), intorder);

          MappedIntegrationRule<2,3> mirx(iri, trafoi, lh);
          MappedIntegrationRule<2,3> miry(irj, trafoj, lh);

          FlatMatrix<> shapesi(feli.GetNDof(), test_evaluator->Dim()*iri.Size(), lh);
          FlatMatrix<> shapesj(felj.GetNDof(), trial_evaluator->Dim()*irj.Size(), lh);
          shapesi = 0.0;
          shapesj = 0.0;
          test_evaluator -> CalcMatrix(feli, mirx, Trans(shapesi), lh);
          trial_evaluator-> CalcMatrix(felj, miry, Trans(shapesj), lh);

          for (auto term : kernel.terms)
            {
              HeapReset hr(lh);
              FlatMatrix<value_type> kernel_ixiy(iri.Size(), irj.Size(), lh);
              for (int ix = 0; ix < iri.Size(); ix++)
                {
                  for (int iy = 0; iy < irj.Size(); iy++)
                    {
                      Vec<3> x = mirx[ix].GetPoint();
                      Vec<3> y = miry[iy].GetPoint();

                      Vec<3> nx = mirx[ix].GetNV();
                      Vec<3> ny = miry[iy].GetNV();
                      value_type kernel_ = kernel.Evaluate(x, y, nx, ny)(term.kernel_comp);

                      double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
                      kernel_ixiy(ix, iy) = term.fac*fac*kernel_;
                    }
                }


              FlatMatrix<value_type> kernel_shapesj(iri.Size(), felj.GetNDof(), lh);
              FlatMatrix<> shapesi1(feli.GetNDof(), iri.Size(), lh);
              FlatMatrix<> shapesj1(felj.GetNDof(), irj.Size(), lh);

              for (int j = 0; j < iri.Size(); j++)
                shapesi1.Col(j) = shapesi.Col(test_evaluator->Dim()*j+term.test_comp);
              for (int j = 0; j < irj.Size(); j++)
                shapesj1.Col(j) = shapesj.Col(trial_evaluator->Dim()*j+term.trial_comp);

              kernel_shapesj = kernel_ixiy * Trans(shapesj1);
              matrix += shapesi1 * kernel_shapesj;
            }
          break;
        }
      default:
        throw Exception ("not possible2 ");
      }
  }



  template <typename KERNEL>
  std::variant<Matrix<double>, Matrix<Complex>> GenericIntegralOperator<KERNEL> ::
  CalcSubMatrix (FlatArray<DofId> target_ids, FlatArray<DofId> source_ids, LocalHeap &lh) const
  {
    auto nearfield = dynamic_pointer_cast<SparseMatrix<value_type>>(GetNearFieldMatrix());
    if (!nearfield)
      throw Exception("GenericIntegralOperator::CalcSubMatrix expects sparse nearfield matrix");

    Matrix<value_type> mat(target_ids.Size(), source_ids.Size());
    mat = value_type(0.0);

    for (int i = 0; i < mat.Height(); i++)
      for (int j = 0; j < mat.Width(); j++)
        if (auto pos = nearfield->GetPositionTest(target_ids[i], source_ids[j]); pos != -1)
          mat(i,j) = nearfield->GetValues()[pos];


    auto CreateDof2El = [&] (const FESpace & fes, VorB vb)
    {
      TableCreator<tuple<ElementId,int>,DofId> creator;
      Array<DofId> dofs;
      for ( ; !creator.Done(); creator++)
        for (auto ei : fes.Elements(vb))
          {
            fes.GetDofNrs(ei, dofs);
            for (auto [nr,d] : Enumerate(dofs))
              creator.Add (d, { ei, nr });
          }
      return make_unique<Table<tuple<ElementId,int>,DofId>> (creator.MoveTable());
    };

    static std::mutex lock;
    if (!trial_dof2el)
      {
        lock_guard<mutex> guard(lock);
        if (!trial_dof2el)
          {
            trial_dof2el = CreateDof2El(*trial_space, trial_vb);
            test_dof2el = CreateDof2El(*test_space, test_vb);
          }
      }

    Array<tuple<ElementId,int,int>> target_map;
    Array<tuple<ElementId,int,int>> source_map;
    Array<ElementId> target_els, source_els; // unique ElementIds

    for (auto i : Range(target_ids))
      for (auto [elid,nr] : (*test_dof2el)[target_ids[i]])
        target_map.Append( { elid, i, nr} );

    for (auto [elid,i,nr] : target_map)
      if (!target_els.Contains(elid))
        target_els.Append(elid);

    TableCreator<tuple<int,int>> creator_target(target_els.Size());
    for ( ; !creator_target.Done(); creator_target++)
      for (auto [elid,i,nr] : target_map)
        creator_target.Add (target_els.Pos(elid), tuple<int,int>{ i,nr });

    auto target_el2dof = creator_target.MoveTable();


    for (auto i : Range(source_ids))
      for (auto [elid,nr] : (*trial_dof2el)[source_ids[i]])
        source_map.Append( { elid, i, nr} );

    for (auto [elid,i,nr] : source_map)
      if (!source_els.Contains(elid))
        source_els.Append(elid);

    TableCreator<tuple<int,int>> creator_source(source_els.Size());
    for ( ; !creator_source.Done(); creator_source++)
      for (auto [elid,i,nr] : source_map)
        creator_source.Add (source_els.Pos(elid), tuple<int,int>{ i,nr });

    auto source_el2dof = creator_source.MoveTable();

    auto source_ma = trial_space->GetMeshAccess();
    auto target_ma = test_space->GetMeshAccess();

    for (auto [selnr, source_el] : Enumerate(source_els))
      for (auto [telnr, target_el] : Enumerate(target_els))
        {
          HeapReset hr(lh);
          const FiniteElement & sfel = trial_space->GetFE(source_el, lh);
          const FiniteElement & tfel = test_space->GetFE(target_el, lh);

          bool common_verts = false;

          if (source_ma == target_ma)
            {
              auto sverts = source_ma->GetElement(source_el).Vertices();
              auto tverts = target_ma->GetElement(target_el).Vertices();

              for (auto sv : sverts)
                for (auto tv : tverts)
                  if (sv==tv) common_verts = true;
            }
          if (common_verts) continue;

          FlatMatrix<value_type> elmat(tfel.GetNDof(), sfel.GetNDof(),lh);
          CalcElementMatrix(elmat, source_el, target_el, lh);

          for (auto [tdofind,tnr] : target_el2dof[telnr])
            for (auto [sdofind,snr] : source_el2dof[selnr])
              mat(tdofind, sdofind) += elmat(tnr, snr);
        }

    return mat;
  }


  // ********************************* Potential **********************************************






  template <typename KERNEL>
  shared_ptr<BasePotentialCF> GenericIntegralOperator<KERNEL> ::
  GetPotential(shared_ptr<GridFunction> gf, optional<int> io, bool nearfield_experimental) const
  {
    return  make_shared<PotentialCF<KERNEL>> (gf, trial_vb, trial_definedon, trial_evaluator,
                                              kernel, io.value_or(intorder), nearfield_experimental, io_params);
  }


  template class GenericIntegralOperator<LaplaceSLKernel<3>>;
  template class GenericIntegralOperator<LaplaceSLKernel<3,3>>;
  template class GenericIntegralOperator<LaplaceSLKernel<3,1,Complex>>;
  template class GenericIntegralOperator<LaplaceSLKernel<3,3,Complex>>;
  template class GenericIntegralOperator<LaplaceDLKernel<3>>;
  template class GenericIntegralOperator<LaplaceDLKernel<3,3>>;
  template class GenericIntegralOperator<LaplaceDLKernel<3,1,Complex>>;
  template class GenericIntegralOperator<LaplaceDLKernel<3,3,Complex>>;
  template class GenericIntegralOperator<LameSLKernel<3>>;

  template class GenericIntegralOperator<HelmholtzSLKernel<3>>;
  template class GenericIntegralOperator<HelmholtzSLKernel<3,3>>;
  template class GenericIntegralOperator<HelmholtzSLKernel<3,1,Complex>>;
  template class GenericIntegralOperator<HelmholtzSLKernel<3,3,Complex>>;
  template class GenericIntegralOperator<HelmholtzDLKernel<3>>;
  template class GenericIntegralOperator<HelmholtzDLKernel<3,3>>;
  template class GenericIntegralOperator<HelmholtzDLKernel<3,1,Complex>>;
  template class GenericIntegralOperator<HelmholtzDLKernel<3,3,Complex>>;
  template class GenericIntegralOperator<HelmholtzHSKernel<3>>;

  template class GenericIntegralOperator<CombinedFieldKernel<3>>;
  template class GenericIntegralOperator<CombinedFieldKernel<3,3>>;
  template class GenericIntegralOperator<CombinedFieldKernel<3,1,Complex>>;
  template class GenericIntegralOperator<CombinedFieldKernel<3,3,Complex>>;

  template class GenericIntegralOperator<MaxwellSLKernel<3>>;
  template class GenericIntegralOperator<MaxwellDLKernel<3>>;
  template class GenericIntegralOperator<MaxwellDLKernel<3,Complex>>;

  template class GenericIntegralOperator<DiffLaplaceSLKernel<3>>;
  template class GenericIntegralOperator<DiffLaplaceSLKernel<3,3>>;
  template class GenericIntegralOperator<DiffLaplaceSLKernel<3,1,Complex>>;
  template class GenericIntegralOperator<DiffLaplaceSLKernel<3,3,Complex>>;
  template class GenericIntegralOperator<DiffHelmholtzSLKernel<3>>;
  template class GenericIntegralOperator<DiffHelmholtzSLKernel<3,3>>;
  template class GenericIntegralOperator<DiffHelmholtzSLKernel<3,1,Complex>>;
  template class GenericIntegralOperator<DiffHelmholtzSLKernel<3,3,Complex>>;
}
