#include <solve.hpp>        // everything from ngsolve
#include <cmath>

#include "intrules_SauterSchwab.hpp"
#include "ngbem.hpp"
// #include "hmat.hpp"
#include "fmmoperator.hpp"


namespace ngsbem
{

  
  class DifferentialOperatorWithFactor : public DifferentialOperator
  {
    shared_ptr<DifferentialOperator> diffop;
    shared_ptr<CoefficientFunction> factor;

  public:
    DifferentialOperatorWithFactor (shared_ptr<DifferentialOperator> adiffop,
                                    shared_ptr<CoefficientFunction> afactor)
      : DifferentialOperator(afactor->Dimensions()[0], 1, adiffop->VB(), adiffop->DiffOrder()),
        diffop(adiffop), factor(afactor)
    { ; }

    int DimRef() const override { return diffop->DimRef(); }



    void CalcMatrix (const FiniteElement & fel,
                     const BaseMappedIntegrationPoint & mip,
                     BareSliceMatrix<double,ColMajor> mat, 
                     LocalHeap & lh) const override
    {
      FlatMatrix<double,ColMajor> hmat(diffop->Dim(), fel.GetNDof(), lh);
      diffop -> CalcMatrix (fel, mip, hmat, lh);

      auto dims = factor->Dimensions();      
      FlatMatrix<double> factorx(dims[0], dims[1], lh);
      factor->Evaluate (mip, factorx.AsVector());

      mat = factorx * hmat;
    }
    
    void CalcMatrix (const FiniteElement & fel,
                     const SIMD_BaseMappedIntegrationRule & mir,
                     BareSliceMatrix<SIMD<double>> mat) const override
    {
      // *testout << "CalcMatrix SIMD" << endl;      
      Matrix<SIMD<double>> hmat (fel.GetNDof()*diffop->Dim(), mir.Size());
      diffop -> CalcMatrix (fel, mir, hmat);

      Matrix<SIMD<double>> fac(factor->Dimension(), mir.Size());
      factor -> Evaluate (mir, fac);

      auto dims = factor -> Dimensions();

      mat.Rows(fel.GetNDof()*dims[1]).Cols(mir.Size()) = SIMD<double>(0.0);
      
      for (size_t i = 0; i < mir.Size(); i++)
        for (size_t j = 0; j < dims[0]; j++)
          for (size_t k = 0; k < dims[1]; k++)
            mat.Col(i).Slice(j,dims[0]) += fac(j*dims[1]+k, i) * hmat.Col(i).Slice(k, dims[1]);
    }
    

    void CalcMatrix (const FiniteElement & fel,
                     const IntegrationPoint & ip,
                     BareSliceMatrix<double,ColMajor> mat,
                     LocalHeap & lh) const override
    {
      diffop -> CalcMatrix(fel, ip, mat, lh);
    }

    void CalcTransformationMatrix (const BaseMappedIntegrationPoint & mip,
                                   SliceMatrix<double> trans,
                                   LocalHeap & lh) const override
    {
      HeapReset hr(lh);
      auto dims = factor->Dimensions();
      
      FlatMatrix<double> factorx(dims[0], dims[1], lh);
      factor->Evaluate (mip, factorx.AsVector());

      FlatMatrix<double> basetrans(diffop->Dim(), diffop->DimRef(), lh);
      diffop -> CalcTransformationMatrix(mip, basetrans, lh);

      trans = factorx * basetrans;
    }
    
  };


  
  IntegralOperator ::
  IntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                   optional<Region> _trial_definedon, optional<Region> _test_definedon,
                   shared_ptr<DifferentialOperator> _trial_evaluator, shared_ptr<CoefficientFunction> _trial_factor,
                   shared_ptr<DifferentialOperator> _test_evaluator, shared_ptr<CoefficientFunction> _test_factor,
                   int _intorder)
    : trial_space(_trial_space), test_space(_test_space),
      trial_definedon(_trial_definedon), test_definedon(_test_definedon),
      trial_evaluator(_trial_evaluator), trial_factor(_trial_factor),
      test_evaluator(_test_evaluator), test_factor(_test_factor),
      intorder(_intorder)
  {
    if (!test_space)
      test_space = trial_space;

    if (trial_factor)
      trial_evaluator = make_shared<DifferentialOperatorWithFactor>(trial_evaluator, trial_factor);

    if (test_factor)
      test_evaluator = make_shared<DifferentialOperatorWithFactor>(test_evaluator, test_factor);
    
    tie(identic_panel_x, identic_panel_y, identic_panel_weight) = IdenticPanelIntegrationRule(intorder);
    tie(common_vertex_x, common_vertex_y, common_vertex_weight) = CommonVertexIntegrationRule(intorder);
    tie(common_edge_x, common_edge_y, common_edge_weight) = CommonEdgeIntegrationRule(intorder);
  }



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


  

  template <typename KERNEL>
  shared_ptr<BaseMatrix> GenericIntegralOperator<KERNEL> ::
  CreateMatrixFMM(LocalHeap & lh) const
  {
    static Timer tall("ngbem fmm setup"); RegionTimer r(tall);
    Array<Vec<3>> xpts, ypts, xnv, ynv;
    IntegrationRule ir(ET_TRIG, intorder);
    auto trial_mesh = trial_space->GetMeshAccess();
    auto test_mesh = test_space->GetMeshAccess();

    Array<int> compress_trial_els(trial_mesh->GetNE(BND));
    Array<int> compress_test_els(test_mesh->GetNE(BND));
    compress_trial_els = -1;
    compress_test_els = -1;
    
    int cnt = 0;
    for (auto el : trial_mesh->Elements(BND))
      if (trial_space->DefinedOn(el))
        if (!trial_definedon || (*trial_definedon).Mask().Test(trial_mesh->GetElIndex(el)))      
          {
            HeapReset hr(lh);
            auto & trafo = trial_mesh->GetTrafo(el, lh);
            auto & mir = static_cast<MappedIntegrationRule<2,3>&>(trafo(ir, lh));
            for (auto & mip : mir)
              {
                xpts.Append(mip.GetPoint());
                xnv.Append(mip.GetNV());
              }
            compress_trial_els[el.Nr()] = cnt++;
          }
    
    cnt = 0;
    for (auto el : test_mesh->Elements(BND))
      if (test_space->DefinedOn(el))      
        if (!test_definedon || (*test_definedon).Mask().Test(test_mesh->GetElIndex(el)))            
          {
            HeapReset hr(lh);
            auto & trafo = test_mesh->GetTrafo(el, lh);
            auto & mir = static_cast<MappedIntegrationRule<2,3>&>(trafo(ir, lh));        
            for (auto & mip : mir)
              {
                ypts.Append(mip.GetPoint());
                ynv.Append(mip.GetNV());        
              }
            compress_test_els[el.Nr()] = cnt++;
          }
    
    auto create_eval = [&](const FESpace & fes,
                           const Array<int> & compress_els,
                           const DifferentialOperator & evaluator)
    {
      auto mesh = fes.GetMeshAccess();
      Array<short> classnr(mesh->GetNE(BND));
      mesh->IterateElements
        (BND, lh, [&] (auto el, LocalHeap & llh)
        {
          if (el.GetType() == ET_QUAD)
            classnr[el.Nr()] = -1;
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

      shared_ptr<BaseMatrix> evalx;
      
      for (auto elclass_inds : table)
        {
          if (elclass_inds.Size() == 0) continue;
          ElementId ei(BND, elclass_inds[0]);
          auto & felx = fes.GetFE (ei, lh);

          int dim = evaluator.DimRef();
          Matrix<double,ColMajor> bmat_(dim*ir.Size(), felx.GetNDof());
          bmat_ = 0.0;
          for (int i : Range(ir.Size()))
            evaluator.CalcMatrix(felx, ir[i], bmat_.Rows(dim*i, dim*(i+1)), lh);
          
          Matrix bmat = bmat_;
          
          Table<DofId> xdofsin(elclass_inds.Size(), felx.GetNDof());
          Table<DofId> xdofsout(elclass_inds.Size(), bmat.Height());

          Array<DofId> dnumsx, dnumsy;
          for (auto i : Range(elclass_inds))
            {
              ElementId ei(BND, elclass_inds[i]);
              fes.GetDofNrs(ei, dnumsx);
              xdofsin[i] = dnumsx;
              
              for (int j = 0; j < dim*ir.Size(); j++)
                xdofsout[i][j] = compress_els[elclass_inds[i]]*(dim*ir.Size())+j;
            }

          auto part_evalx = make_shared<ConstantElementByElementMatrix<typename KERNEL::value_type>>
            (mesh->GetNE(BND)*ir.Size()*dim, fes.GetNDof(),
             bmat, std::move(xdofsout), std::move(xdofsin));
          
          if (evalx)
            evalx = evalx + part_evalx;
          else
            evalx = part_evalx;
        }

      int cnt = 0;
      for (auto nr : compress_els) if (nr!=-1) cnt++;

      Tensor<3, typename KERNEL::value_type> weights(cnt*ir.Size(),
                                                     evaluator.Dim(), evaluator.DimRef());
      Matrix<double> transformation(evaluator.Dim(), evaluator.DimRef()); 
      
      for (auto el : mesh->Elements(BND))
        if (compress_els[el.Nr()] != -1)
          {
            HeapReset hr(lh);
            auto & trafo = mesh->GetTrafo(el, lh);
            auto & mir = trafo(ir, lh);
            for (auto j : Range(mir.Size()))
              {
                evaluator.CalcTransformationMatrix(mir[j], transformation, lh);
                weights(compress_els[el.Nr()]*mir.Size()+j,STAR,STAR) =
                  mir[j].GetWeight()*transformation;
              }
          }

      auto diagmat = make_shared<BlockDiagonalMatrix<typename KERNEL::value_type>>(std::move(weights));
      
      return diagmat*evalx;
    };

    auto evalx = create_eval(*trial_space, compress_trial_els, *trial_evaluator);
    auto evaly = create_eval(*test_space, compress_test_els, *test_evaluator);    
    auto fmmop = make_shared<FMM_Operator<KERNEL>> (kernel, std::move(xpts), std::move(ypts),
                                                    std::move(xnv), std::move(ynv));


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
    for (ElementId ei : trial_mesh->Elements(BND))
      if (trial_space->DefinedOn(ei))      
        if (!trial_definedon || (*trial_definedon).Mask().Test(trial_mesh->GetElIndex(ei)))     
        {
          other.SetSize0();
          for (auto v : trial_mesh->GetElement(ei).Vertices())
            for (auto ej : trial_mesh->GetVertexElements(v,BND))
              if (test_space->DefinedOn(ElementId(BND,ej)))                    
                if (!test_definedon || (*test_definedon).Mask().Test(test_mesh->GetElIndex(ElementId(BND,ej))))
                  {
                    if (!other.Contains(ej))
                      {
                        other.Append (ej );
                        pairs.Append ( { ei.Nr(), ej });
                      }
                  }
        }
    
    tfind.Stop();
    tsetupgraph.Start();
    TableCreator<int> trial_elements_creator(pairs.Size()), test_elements_creator(pairs.Size());

    for ( ; !trial_elements_creator.Done(); trial_elements_creator++, test_elements_creator++)
      {
        Array<DofId> dnums;
        for (auto i : Range(pairs))
          {
            auto [ei_trial, ei_test] = pairs[i];
            trial_space->GetDofNrs( { BND, ei_trial }, dnums);
            trial_elements_creator.Add (i, dnums);

            test_space->GetDofNrs( { BND, ei_test }, dnums);
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



    TableCreator<int> create_nbels(trial_mesh->GetNE(BND));
    for ( ; !create_nbels.Done(); create_nbels++)    
      for (auto i : Range(pairs))
        create_nbels.Add (get<0>(pairs[i]), get<1>(pairs[i]));
    Table<int> nbels = create_nbels.MoveTable();

    trial_mesh->IterateElements
      (BND, lh, [&](auto ei_trial, LocalHeap & lh)
      {
        for (auto nrtest : nbels[ei_trial.Nr()])
          {
            ElementId ei_test(BND, nrtest);
            
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
            
            // subtract terms from fmm:
            
            MappedIntegrationRule<2,3> trial_mir(ir, trial_trafo, lh);
            MappedIntegrationRule<2,3> test_mir(ir, test_trafo, lh);
            
            FlatMatrix<> shapesi(test_fel.GetNDof(), test_evaluator->Dim()*ir.Size(), lh);
            FlatMatrix<> shapesj(trial_fel.GetNDof(), trial_evaluator->Dim()*ir.Size(), lh);
            shapesi = 0.;
            shapesj = 0.;
            test_evaluator -> CalcMatrix(test_fel, test_mir, Trans(shapesi), lh);
            trial_evaluator-> CalcMatrix(trial_fel, trial_mir, Trans(shapesj), lh);
            
            for (auto term : kernel.terms)
              {
                HeapReset hr(lh);
                FlatMatrix<value_type> kernel_ixiy(ir.Size(), ir.Size(), lh);
                for (int ix = 0; ix < ir.Size(); ix++)
                  {
                    for (int iy = 0; iy < ir.Size(); iy++)
                      {
                        Vec<3> x = test_mir[ix].GetPoint();
                        Vec<3> y = trial_mir[iy].GetPoint();
                        
                        Vec<3> nx = test_mir[ix].GetNV();
                        Vec<3> ny = trial_mir[iy].GetNV();
                        value_type kernel_ = 0.0;
                        if (L2Norm2(x-y) > 0)
                          kernel_ = kernel.Evaluate(x, y, nx, ny)(term.kernel_comp);
                        
                        double fac = test_mir[ix].GetWeight()*trial_mir[iy].GetWeight();
                        kernel_ixiy(ix, iy) = term.fac*fac*kernel_;
                      }
                  }
                
                FlatMatrix<value_type> kernel_shapesj(ir.Size(), trial_fel.GetNDof(), lh);
                FlatMatrix<> shapesi1(test_fel.GetNDof(), ir.Size(), lh);
                FlatMatrix<> shapesj1(trial_fel.GetNDof(), ir.Size(), lh);
                
                for (int j = 0; j < ir.Size(); j++)
                  {
                    shapesi1.Col(j) = shapesi.Col(test_evaluator->Dim()*j+term.test_comp);
                    shapesj1.Col(j) = shapesj.Col(trial_evaluator->Dim()*j+term.trial_comp);
                  }
                kernel_shapesj = kernel_ixiy * Trans(shapesj1);
                elmat -= shapesi1 * kernel_shapesj;
              }
            // tasscorr.Stop();        
            nearfield_correction -> AddElementMatrix (test_dnums, trial_dnums, elmat, true);
          }
      });

    
    tassemble.Stop();
    return TransposeOperator(evaly) * fmmop * evalx + nearfield_correction;
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
    IntegrationRule irtrig(ET_TRIG, intorder);
    matrix = 0.; 
    // t1.Stop();

    Vec<3> x,y,nx,ny;
    typedef decltype(kernel.Evaluate (x,y,nx,ny)) KERNEL_COMPS_T;

    
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

      mshapesi = 0.;
      mshapesj = 0.;
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
          
          AddABt (mshapesi_kern, 
                  mshapesj.RowSlice(term.trial_comp, trial_evaluator->Dim()).AddSize(felj.GetNDof(), miry.Size()),
                  elmat);
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
    
    switch (n_common_vertices)
      {
      case 3: //identical panel
        {
          // RegionTimer reg(t_identic);    
          
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
          
          const EDGE * edges = ElementTopology::GetEdges (ET_TRIG); // 0 1 | 1 2 | 2 0 
          int cex, cey;
          for (int cx = 0; cx < 3; cx++)
            for (int cy = 0; cy < 3; cy++)
              {
                IVec<2> ex (verti[edges[cx][0]], verti[edges[cx][1]]); 
                IVec<2> ey (vertj[edges[cy][0]], vertj[edges[cy][1]]); 
                if (ex.Sort() == ey.Sort()) 
                  {
                    cex = cx;  // -> "common" edge number triangle i
                    cey = cy;  // -> "common" edge number triangle j
                    break;
                  }
              }
          int vpermx[3] = { edges[cex][0], edges[cex][1], -1 }; // common edge gets first
          vpermx[2] = 3-vpermx[0]-vpermx[1]; 
          int vpermy[3] = { edges[cey][1], edges[cey][0], -1 }; // common edge gets first
          vpermy[2] = 3-vpermy[0]-vpermy[1];
          
          constexpr int BS = 128;
          for (int k = 0; k < common_edge_weight.Size(); k+=BS)
            {
              int num = std::min(size_t(BS), common_edge_weight.Size()-k);
              
              HeapReset hr(lh);
              
              IntegrationRule irx(num, lh);
              IntegrationRule iry(num, lh);
              
              for (int k2 = 0; k2 < num; k2++)
                {
                  Vec<2> xk = common_edge_x[k+k2];
                  Vec<2> yk = common_edge_y[k+k2];
                  
                  Vec<3> lamx (1-xk(0)-xk(1), xk(0), xk(1) );
                  Vec<3> lamy (1-yk(0)-yk(1), yk(0), yk(1) );
                  
                  Vec<3> plamx, plamy;
                  for (int i = 0; i < 3; i++)
                    {
                      plamx(vpermx[i]) = lamx(i);
                      plamy(vpermy[i]) = lamy(i);
                    }
                  
                  irx[k2] = IntegrationPoint(plamx(0), plamx(1), 0, common_edge_weight[k+k2]);
                  iry[k2] = IntegrationPoint(plamy(0), plamy(1), 0, 0);
                }
              
              Integrate4D (irx, iry, feli, felj, trafoi, trafoj, matrix, lh);
            }
          
          break;
        }
        
      case 1: //common vertex
        {
          // RegionTimer reg(t_common_vertex);    
          
          int cvx=-1, cvy=-1;
          for (int cx = 0; cx < 3; cx++)
            for (int cy = 0; cy < 3; cy++)
              {
                if (verti[cx] == vertj[cy])
                  {
                    cvx = cx;
                    cvy = cy;
                    break;
                  }
              }

          int vpermx[3] = { cvx, (cvx+1)%3, (cvx+2)%3 };
          vpermx[2] = 3-vpermx[0]-vpermx[1];
          int vpermy[3] = { cvy, (cvy+1)%3, (cvy+2)%3 };
          vpermy[2] = 3-vpermy[0]-vpermy[1];
          
          // vectorized version:
          constexpr int BS = 128;
          for (int k = 0; k < common_vertex_weight.Size(); k+=BS)
            {
              int num = std::min(size_t(BS), common_vertex_weight.Size()-k);
              
              HeapReset hr(lh);
              
              IntegrationRule irx(num, lh);
              IntegrationRule iry(num, lh);
              
              for (int k2 = 0; k2 < num; k2++)
                {
                  Vec<2> xk = common_vertex_x[k+k2];
                  Vec<2> yk = common_vertex_y[k+k2];
                  
                  Vec<3> lamx (1-xk(0)-xk(1), xk(0), xk(1) );
                  Vec<3> lamy (1-yk(0)-yk(1), yk(0), yk(1) );
                  
                  Vec<3> plamx, plamy;
                  for (int i = 0; i < 3; i++)
                    {
                      plamx(vpermx[i]) = lamx(i);
                      plamy(vpermy[i]) = lamy(i);
                    }
                  
                  irx[k2] = IntegrationPoint(plamx(0), plamx(1), 0, common_vertex_weight[k+k2]);
                  iry[k2] = IntegrationPoint(plamy(0), plamy(1), 0, 0);
                }
              
              Integrate4D (irx, iry, feli, felj, trafoi, trafoj, matrix, lh);
            }
          
          break;
        }
        
      case 0: //disjoint panels
        {
          // RegionTimer r(t_disjoint);    
          
          // shapes+geom out of loop, matrix multiplication
          MappedIntegrationRule<2,3> mirx(irtrig, trafoi, lh);
          MappedIntegrationRule<2,3> miry(irtrig, trafoj, lh);
          
          FlatMatrix<> shapesi(feli.GetNDof(), test_evaluator->Dim()*irtrig.Size(), lh);
          FlatMatrix<> shapesj(felj.GetNDof(), trial_evaluator->Dim()*irtrig.Size(), lh);
          shapesi = 0.0;
          shapesj = 0.0;
          test_evaluator -> CalcMatrix(feli, mirx, Trans(shapesi), lh);
          trial_evaluator-> CalcMatrix(felj, miry, Trans(shapesj), lh);
          
          for (auto term : kernel.terms)
            {
              HeapReset hr(lh);
              FlatMatrix<value_type> kernel_ixiy(irtrig.Size(), irtrig.Size(), lh);
              for (int ix = 0; ix < irtrig.Size(); ix++)
                {
                  for (int iy = 0; iy < irtrig.Size(); iy++)
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
              
              
              FlatMatrix<value_type> kernel_shapesj(irtrig.Size(), felj.GetNDof(), lh);
              FlatMatrix<> shapesi1(feli.GetNDof(), irtrig.Size(), lh);
              FlatMatrix<> shapesj1(felj.GetNDof(), irtrig.Size(), lh);
              
              for (int j = 0; j < irtrig.Size(); j++)
                {
                  shapesi1.Col(j) = shapesi.Col(test_evaluator->Dim()*j+term.test_comp);
                        shapesj1.Col(j) = shapesj.Col(trial_evaluator->Dim()*j+term.trial_comp);
                }
              kernel_shapesj = kernel_ixiy * Trans(shapesj1);
              matrix += shapesi1 * kernel_shapesj;
            }
          break;
        }
      default:
        throw Exception ("not possible2 ");
      }
  }







  // ********************************* Potential ********************************************** 
  



  
  
  template <typename KERNEL>
  shared_ptr<BasePotentialCF> GenericIntegralOperator<KERNEL> ::
  GetPotential(shared_ptr<GridFunction> gf, optional<int> io, bool nearfield_experimental) const
  {
    return  make_shared<PotentialCF<KERNEL>> (gf, trial_definedon, trial_evaluator,
                                              kernel, io.value_or(intorder), nearfield_experimental);
  }


  template <typename KERNEL>
  PotentialCF<KERNEL> ::
  PotentialCF (shared_ptr<GridFunction> _gf,
               optional<Region> _definedon,                   
               shared_ptr<DifferentialOperator> _evaluator,
               KERNEL _kernel, int _intorder, bool _nearfield)
    : BasePotentialCF(_gf, _definedon, _evaluator, std::is_same<typename KERNEL::value_type,Complex>()),
      kernel(_kernel), intorder(_intorder), nearfield(_nearfield)
  {
    ;
  }



  template <typename KERNEL>
  void PotentialCF<KERNEL> ::
  BuildLocalExpansion(const Region & reg)
  {
    LocalHeapMem<100000> lh("PotentialCF::BuildLocalExpansion");

    auto space = this->gf->GetFESpace();
    auto mesh = space->GetMeshAccess();

    
    Vec<3> smax(-1e99, -1e99, -1e99);
    Vec<3> smin(1e99, 1e99, 1e99);
    
    for (size_t i = 0; i < mesh->GetNSE(); i++)
      {
        HeapReset hr(lh);
        ElementId ei(BND, i);
        if (!space->DefinedOn(ei)) continue;

        // const FiniteElement &fel = space->GetFE(ei, lh);
        const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);
        IntegrationRule ir(trafo.GetElementType(), intorder);
        MappedIntegrationRule<2,3> miry(ir, trafo, lh);

        for (int k = 0; k < miry.Size(); k++)
          for (int j = 0; j < 3; j++)
            {
              smin(j) = min(smin(j), miry[k].GetPoint()(j));
              smax(j) = max(smax(j), miry[k].GetPoint()(j));
            }
      }

    Vec<3> cs = 0.5*(smin+smax);
    double rs = MaxNorm(smax-smin);

    // cout << "cs = " << cs << ", rs = " << rs << endl;
    auto singmp = kernel.CreateMultipoleExpansion(cs, rs);
    
    typedef typename KERNEL::value_type T;
    for (size_t i = 0; i < mesh->GetNSE(); i++)
      {
        HeapReset hr(lh);
        ElementId ei(BND, i);

        if (!space->DefinedOn(ei)) continue;
          
        const FiniteElement &fel = space->GetFE(ei, lh);
        const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);
        
        Array<DofId> dnums(fel.GetNDof(), lh);
        space->GetDofNrs(ei, dnums);
        FlatVector<T> elvec(fel.GetNDof(), lh);
        gf->GetElementVector(dnums, elvec);
        
        IntegrationRule ir(fel.ElementType(), intorder);
        MappedIntegrationRule<2,3> miry(ir, trafo, lh);
        FlatMatrix<T> vals(miry.Size(), evaluator->Dim(), lh);
        
        evaluator->Apply (fel, miry, elvec, vals, lh);
        
        // add vals to multipole ...
        for (int j = 0; j < miry.Size(); j++)
          {
            vals.Row(j) *= miry[j].GetWeight();
            kernel.AddSource (*singmp, miry[j].GetPoint(), miry[j].GetNV(), make_BareSliceVector(vals.Row(j)));
          }
      }

    singmp->CalcMP();

    Vec<3> tmax(-1e99, -1e99, -1e99);
    Vec<3> tmin(1e99, 1e99, 1e99);
    Array<Vec<3>> tpoints;
    auto tmesh = reg.Mesh();
    for (auto el : reg.GetElements())
      {      
        HeapReset hr(lh);
        
        const ElementTransformation &trafo = tmesh->GetTrafo(el, lh);
        IntegrationRule ir(trafo.GetElementType(), intorder);
        MappedIntegrationRule<2,3> miry(ir, trafo, lh);

        for (int k = 0; k < miry.Size(); k++)
          {
            tpoints.Append (Vec<3>(miry[k].GetPoint()));
            for (int j = 0; j < 3; j++)
              {
                tmin(j) = min(tmin(j), miry[k].GetPoint()(j));
                tmax(j) = max(tmax(j), miry[k].GetPoint()(j));
              }
          }
      }

    Vec<3> ct = 0.5*(tmin+tmax);
    double rt = MaxNorm(tmax-tmin);
    // cout << "ct = " << ct << ", rt = " << rt << endl;
    
    local_expansion = kernel.CreateLocalExpansion(ct, rt); 
    
    ParallelFor (tpoints.Size(), [&](int i){
      local_expansion->AddTarget(tpoints[i]);
    });

    local_expansion->CalcMP(singmp, false);
  }



  

  // minimize the function f(x) = 1/2 x a x + b x + c
  // returns (x, f(x))
  tuple<double,double> MinimizeOnSegm (double a, double b, double c)
  {
    if (a > 0)
      {
        double x = -b/a;
        if (x > 0 && x < 1)
          return { x, 0.5*a*x*x + b*x + c };
      }
    double val0 = c;
    double val1 = 0.5*a+b+c;
    if (val0 < val1)
      return { 0, val0 };
    else
      return { 1, val1 };
  }
                         
  // 1/2 x^T A x + b x + c
  Vec<2> MinimizeOnTrig (Mat<2,2> a, Vec<2> b, double c)
  {
    if (a(0,0) > 0 && Det(a) > 0)
      {
        Vec<2> sol = -Inv(a)*b;
        if (sol(0) > 0 && sol(1) > 0 && sol(0)+sol(1) < 1)
          return sol;
      }

    
    auto [x0,val0] = MinimizeOnSegm(a(0,0), b(0), c);
    auto [x1,val1] = MinimizeOnSegm(a(1,1), b(1), c);
    Vec<2> p{1,0}, d{-1,1};
    Vec<2> g = a*p + b;
    auto [x2,val2] = MinimizeOnSegm(InnerProduct(a*d,d), InnerProduct(g,d), 0.5*a(0,0)+b(0)+c);  // (1,0) + (-1,1)*s

    if (val0 < val1 && val0 < val2)
      return { x0, 0 };
    if (val1 < val2)
      return { 0, x1 };
    return p+x2*d;
  }

  
  IntegrationRule GetIntegrationRule(Vec<3> x, const ElementTransformation & trafo, int intorder, bool nearfield)
  {
    if (!nearfield || trafo.GetElementType() != ET_TRIG)
      return IntegrationRule(trafo.GetElementType(), intorder);


    IntegrationPoint ip(1.0/3, 1.0/3);
    MappedIntegrationPoint<2,3>  mip(ip, trafo);
    double elsize = L2Norm(mip.GetJacobian());
    double dist = L2Norm(x-mip.GetPoint());

    if (dist < elsize)
      {
        // use SQP to find projection of x onto (curved) triangle
        
        IntegrationPoint ip(1./3, 1./3);
        for (int j = 0; j < 5; j++) // SQP steps
          {
            MappedIntegrationPoint<2,3>  mip(ip, trafo);
            // dist = || x - (mip+Jac*(uv-ip) + 1/2*Hesse(uv-ip, uv-ip)) ||
            Mat<3,2> jac = mip.GetJacobian();
            Vec<2> ipvec { ip(0), ip(1) };
            auto Hesse = mip.CalcHesse();
            Vec<3,Vec<2>> Hesseip
              {
                Hesse[0]*ipvec,
                Hesse[1]*ipvec,
                Hesse[2]*ipvec
              };
            Vec<3> Hesseipip
              {
                InnerProduct(Hesseip(0), ipvec),
                InnerProduct(Hesseip(1), ipvec),
                InnerProduct(Hesseip(2), ipvec)
              };
            Mat<3,2> jacphip = jac;
            jacphip.Row(0) -= Hesseip(0);
            jacphip.Row(1) -= Hesseip(1);
            jacphip.Row(2) -= Hesseip(2);
            Mat<2,2> a = Trans(jac)*jac;
            Vec<2> b = -Trans(jacphip) * (x-mip.GetPoint()+jac*ipvec + 0.5*Hesseipip);
            Vec<2> uv = MinimizeOnTrig(a, b, 0);
            ip = IntegrationPoint(uv(0), uv(1));
          }

        
        // generate Duffy integration rules on split triangles
        IntegrationRule irsegm(ET_SEGM, intorder);
        IntegrationRule irtrig(trafo.GetElementType(), intorder);            
        IntegrationRule ir;

        Vec<2> corners[] = { Vec<2>(0,0), Vec<2>(1,0), Vec<2>(0,1) };
        for (int j = 0; j < 3; j++)
          {
            Vec<2> v0 = corners[j];
            Vec<2> v1 = corners[(j+1)%3];
            Vec<2> v2 { ip(0), ip(1) };
            Mat<2,2> sides;
            sides.Col(0) = v0-v2;
            sides.Col(1) = v1-v2;
            double factor = Det(sides);

            if (fabs(factor) > 1e-12)
              for (auto ip1 : irsegm)
                for (auto ip2 : irsegm)
                  {
                    Vec<2> ipxy = v0 + ip1(0)*(1-ip2(0))*(v1-v0) + ip2(0)*(v2-v0);
                    ir.AddIntegrationPoint (IntegrationPoint(ipxy(0), ipxy(1), 0,
                                                             ip1.Weight()*ip2.Weight()*(1-ip2(0))*factor));
                  }
          }
        return ir;
      }
    return IntegrationRule(trafo.GetElementType(), intorder);    
  }


  
  template <typename KERNEL> template <typename T>
  void PotentialCF<KERNEL> :: T_Evaluate(const BaseMappedIntegrationPoint & mip,
                                         FlatVector<T> result) const
  {
    static Timer t("ngbem evaluate potential (ip)"); RegionTimer reg(t);
    LocalHeapMem<100000> lh("Potential::Eval");
    auto space = this->gf->GetFESpace();
    auto mesh = space->GetMeshAccess();

    Vector<SIMD<T>> simd_result(Dimension());
    simd_result = SIMD<T>(0.0);
    auto & mip23 = dynamic_cast<const MappedIntegrationPoint<2,3>&>(mip);
    if constexpr (std::is_same<typename KERNEL::value_type,T>())
      for (size_t i = 0; i < mesh->GetNSE(); i++)
        {
          HeapReset hr(lh);
          ElementId ei(BND, i);
          if (!space->DefinedOn(ei)) continue;
          if (definedon &&  !(*definedon).Mask().Test(mesh->GetElIndex(ei))) continue;
            
          const FiniteElement &fel = space->GetFE(ei, lh);
          const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);
          
          Array<DofId> dnums(fel.GetNDof(), lh);
          space->GetDofNrs(ei, dnums);
          FlatVector<T> elvec(fel.GetNDof(), lh);
          gf->GetElementVector(dnums, elvec);
          
          // IntegrationRule ir(fel.ElementType(), intorder);
          IntegrationRule ir = GetIntegrationRule(mip23.GetPoint(), trafo, intorder, nearfield);
          
          SIMD_IntegrationRule simd_ir(ir);

          static constexpr int bs = 64;
          for (int k = 0; k < simd_ir.Size(); k += bs)
            {
              HeapReset hr(lh);
              auto simd_ir_range = simd_ir.Range(k, min(simd_ir.Size(), size_t(k+bs)));
              SIMD_MappedIntegrationRule<2,3> miry(simd_ir_range, trafo, lh);
              FlatMatrix<SIMD<T>> vals(evaluator->Dim(), miry.Size(), lh);
              
              evaluator->Apply (fel, miry, elvec, vals);
              for (int iy = 0; iy < miry.Size(); iy++)
                {
                  Vec<3,SIMD<double>> x = mip23.GetPoint();
                  Vec<3,SIMD<double>> nx = mip23.GetNV();
                  
                  Vec<3,SIMD<double>> y = miry[iy].GetPoint();
                  Vec<3,SIMD<double>> ny = miry[iy].GetNV();
                  
                  for (auto term : kernel.terms)
                    {
                      auto kernel_ = kernel.Evaluate(x, y, nx, ny)(term.kernel_comp);
                      simd_result(term.test_comp) += miry[iy].GetWeight()*kernel_ * vals(term.trial_comp,iy);
                    }
                }
            }
        }
    for (int i = 0; i < Dimension(); i++)
      result(i) = HSum(simd_result(i));
  }


  template <typename KERNEL> template <typename T>
  void PotentialCF<KERNEL> :: T_Evaluate(const BaseMappedIntegrationRule & bmir,
                                         BareSliceMatrix<T> result) const
  {
    
    if constexpr (std::is_same<typename KERNEL::value_type,T>())
      if (local_expansion)
        {
          auto & mir = dynamic_cast<const MappedIntegrationRule<2,3>&>(bmir);        
          for (int j = 0; j < mir.Size(); j++)
            kernel.EvaluateMP (*local_expansion, Vec<3>(mir[j].GetPoint()), Vec<3>(mir[j].GetNV()), make_BareSliceVector(result.Row(j)));
          return;
        }
    

    
    if (nearfield)
      {
        for (int i = 0; i < bmir.Size(); i++)
          T_Evaluate(bmir[i], result.Row(i).Range(0,Dimension()));
        return;
      }
    
    try
      {
        static Timer t("ngbem evaluate potential (ip)"); RegionTimer reg(t);
        LocalHeapMem<100000> lh("Potential::Eval");
        auto space = this->gf->GetFESpace();
        auto mesh = space->GetMeshAccess();

        auto & mirx = dynamic_cast<const MappedIntegrationRule<2,3>&>(bmir);
        
        Matrix<SIMD<T>> simd_result(Dimension(), mirx.Size());
        simd_result = SIMD<T>(0.0);
        if constexpr (std::is_same<typename KERNEL::value_type,T>())
          for (size_t i = 0; i < mesh->GetNSE(); i++)
            {
              HeapReset hr(lh);
              ElementId ei(BND, i);
              if (!space->DefinedOn(ei)) continue;
              
              const FiniteElement &fel = space->GetFE(ei, lh);
              const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);

              Array<DofId> dnums(fel.GetNDof(), lh);
              space->GetDofNrs(ei, dnums);
              FlatVector<T> elvec(fel.GetNDof(), lh);
              gf->GetElementVector(dnums, elvec);
              
              IntegrationRule ir(fel.ElementType(), intorder);
              SIMD_IntegrationRule simd_ir(ir);
              SIMD_MappedIntegrationRule<2,3> miry(simd_ir, trafo, lh);
              FlatMatrix<SIMD<T>> vals(evaluator->Dim(), miry.Size(), lh);
              
              evaluator->Apply (fel, miry, elvec, vals);
              for (int ix = 0; ix < mirx.Size(); ix++)
                for (int iy = 0; iy < miry.Size(); iy++)
                  {
                    Vec<3,SIMD<double>> x = mirx[ix].GetPoint();
                    Vec<3,SIMD<double>> y = miry[iy].GetPoint();
                    
                    Vec<3,SIMD<double>> nx = mirx[ix].GetNV();
                    Vec<3,SIMD<double>> ny = miry[iy].GetNV();
                    
                    for (auto term : kernel.terms)
                      {
                        auto kernel_ = kernel.Evaluate(x, y, nx, ny)(term.kernel_comp);
                        simd_result(term.test_comp, ix) += miry[iy].GetWeight()*kernel_ * vals(term.trial_comp,iy);
                      }
                  }
            }
        for (int i = 0; i < Dimension(); i++)
          for (int j = 0; j < mirx.Size(); j++)
            result(j, i) = HSum(simd_result(i,j));
      }
    catch (ExceptionNOSIMD & e)
      {
        e.Append ("\nin PotentialCF::Evaluate(mir)");
        throw e;
      }
  }

  
  
  template <typename KERNEL> template <typename T>
  void PotentialCF<KERNEL> :: T_Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
                                         BareSliceMatrix<SIMD<T>> result) const
  {
    static Timer t("ngbem evaluate potential (ir-simd)"); RegionTimer reg(t);
    throw ExceptionNOSIMD ("PotentialCF::Evaluate (SIMD) not available");
    
    result.AddSize(Dimension(), ir.Size()) = SIMD<T>(0.0);
    return;
  }




  template class PotentialCF<LaplaceSLKernel<3>>;
  template class PotentialCF<LaplaceDLKernel<3>>;
  template class PotentialCF<LaplaceHSKernel<3>>;  
  template class PotentialCF<HelmholtzSLKernel<3>>;
  template class PotentialCF<HelmholtzSLVecKernel<3>>;
  template class PotentialCF<HelmholtzDLKernel<3>>;
  template class PotentialCF<HelmholtzHSKernel<3>>;
  template class PotentialCF<CombinedFieldKernel<3>>;  
  template class PotentialCF<MaxwellSLKernel<3>>;
  template class PotentialCF<MaxwellDLKernel<3>>;

  
  
  template class GenericIntegralOperator<LaplaceSLKernel<3>>;
  template class GenericIntegralOperator<LaplaceDLKernel<3>>;
  template class GenericIntegralOperator<LaplaceHSKernel<3>>;
  
  template class GenericIntegralOperator<HelmholtzSLKernel<3>>;
  template class GenericIntegralOperator<HelmholtzSLVecKernel<3>>;
  template class GenericIntegralOperator<HelmholtzDLKernel<3>>;
  template class GenericIntegralOperator<HelmholtzHSKernel<3>>;
  
  template class GenericIntegralOperator<CombinedFieldKernel<3>>;

  template class GenericIntegralOperator<MaxwellSLKernel<3>>;
  template class GenericIntegralOperator<MaxwellDLKernel<3>>;    
}
