#include <solve.hpp>        // everything from ngsolve
#include <cmath>

#include "intrules_SauterSchwab.hpp"
#include "ngbem.hpp"
// #include "hmat.hpp"
#include "fmmoperator.hpp"


namespace ngsbem
{


  
  template <typename T>  
  IntegralOperator<T> ::
  IntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                   optional<Region> _trial_definedon, optional<Region> _test_definedon,
                   BEMParameters _param)
    : trial_space(_trial_space), test_space(_test_space),
      trial_definedon(_trial_definedon), test_definedon(_test_definedon),
      param(_param)
  {
    if (!test_space)
      test_space = trial_space;

    /*
    trial_ct = make_shared<ClusterTree>(trial_space, param.leafsize);
    if (trial_space == test_space)
      test_ct = trial_ct; // the same
    else
      test_ct = make_shared<ClusterTree>(test_space, param.leafsize);
    */

    auto mesh = trial_space->GetMeshAccess(); // trialspace
    auto mesh2 = test_space->GetMeshAccess(); // testspace

    if (trial_definedon)
      cout << "trial is definedon: " << (*trial_definedon).Mask() << endl;
    if (test_definedon)
      cout << "test is definedon: " << (*test_definedon).Mask() << endl;
    
    // setup global-2-boundary mappings;
    BitArray bnddofs(trial_space->GetNDof());
    bnddofs.Clear();
    for (int i = 0; i < mesh->GetNSE(); i++)
      {
        ElementId ei(BND, i);
        if (!trial_definedon || (*trial_definedon).Mask().Test(mesh->GetElIndex(ei)))
          {
            Array<DofId> dnums;
            trial_space->GetDofNrs(ei, dnums);
            for (auto d : dnums)
              bnddofs.SetBit(d);
          }
      }

        
    mapglob2bnd.SetSize(trial_space->GetNDof());
    mapglob2bnd = -1;
    for (int i = 0; i < trial_space->GetNDof(); i++)
      if (bnddofs.Test(i))
	{
	  mapglob2bnd[i] = mapbnd2glob.Size();
	  mapbnd2glob.Append(i);
	}

    // run through surface elements and add elem to related trial dofs
    // Table-creator creates table with one big block of memory,
    // avoids memory fragmentation 
    TableCreator<int> creator;
    Array<DofId> dnumsi;
    for ( ; !creator.Done(); creator++)
      for (int i = 0; i < mesh->GetNSE(); i++)
        {
          trial_space->GetDofNrs( ElementId(BND, i), dnumsi); 
          for (auto d : dnumsi)
            creator.Add (d, i);
        }
    elems4dof = creator.MoveTable();

    BitArray bnddofs2(test_space->GetNDof());
    bnddofs2.Clear();
    for (int i = 0; i < mesh2->GetNSE(); i++)
      {
        ElementId ei(BND, i);        
        if (!test_definedon || (*test_definedon).Mask().Test(mesh->GetElIndex(ei)))
          {
            Array<DofId> dnums;
            test_space->GetDofNrs(ElementId(BND,i), dnums);
            for (auto d : dnums)
              bnddofs2.SetBit(d);
          }
      }
    
    mapglob2bnd2.SetSize(test_space->GetNDof());
    mapglob2bnd2 = -1;
    for (int i = 0; i < test_space->GetNDof(); i++)
      if (bnddofs2.Test(i))
	{
	  mapglob2bnd2[i] = mapbnd2glob2.Size();
	  mapbnd2glob2.Append(i);
	}

    // run through surface elements and add elem to related test dofs
    TableCreator<int> creator2;
    for ( ; !creator2.Done(); creator2++)
      for (int i = 0; i < mesh2->GetNSE(); i++)
        {
          test_space->GetDofNrs( ElementId(BND, i), dnumsi); 
          for (auto d : dnumsi)
            creator2.Add (d, i);
        }
    elems4dof2 = creator2.MoveTable();
  }


  template <typename KERNEL>
  shared_ptr<BaseMatrix> GenericIntegralOperator<KERNEL> ::
  CreateMatrixFMM(LocalHeap & lh) const
  {
    static Timer tall("ngbem fmm setup"); RegionTimer r(tall);
    Array<Vec<3>> xpts, ypts, xnv, ynv;
    IntegrationRule ir(ET_TRIG, param.intorder);
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
          classnr[el.Nr()] = 
            SwitchET<ET_SEGM, ET_TRIG,ET_TET>
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


      /*
      VVector<typename KERNEL::value_type> weights(cnt*ir.Size());
      for (auto el : mesh->Elements(BND))
        if (compress_els[el.Nr()] != -1)
          {
            HeapReset hr(lh);
            auto & trafo = mesh->GetTrafo(el, lh);
            auto & mir = trafo(ir, lh);
            for (auto j : Range(mir.Size()))
              {
                Mat<1,1> transformation;  // todo: general dimensions
                evaluator.CalcTransformationMatrix(mir[j], transformation, lh);
                weights(compress_els[el.Nr()]*mir.Size()+j) = mir[j].GetWeight()*transformation(0,0);
              }
          }
      auto diagmat = make_shared<DiagonalMatrix<typename KERNEL::value_type>>(std::move(weights));
      */

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


    /*
    for (auto i : Range(pairs))
      {
        HeapReset hr(lh);
        ElementId ei_trial(BND,get<0> (pairs[i]));
        ElementId ei_test(BND,get<1> (pairs[i]));

        auto & trial_trafo = trial_mesh -> GetTrafo(ei_trial, lh);
        auto & test_trafo = test_mesh -> GetTrafo(ei_test, lh);
        auto & trial_fel = trial_space->GetFE(ei_trial, lh);
        auto & test_fel = test_space->GetFE(ei_test, lh);

        Array<DofId> trial_dnums(trial_fel.GetNDof(), lh);
        Array<DofId> test_dnums(test_fel.GetNDof(), lh);
        
        trial_space->GetDofNrs (ei_trial, trial_dnums);        
        test_space->GetDofNrs (ei_test, test_dnums);

        FlatMatrix<value_type> elmat(test_dnums.Size(), trial_dnums.Size(), lh);
        tassSS.Start();
        CalcElementMatrix (elmat, ei_trial, ei_test, lh);
        tassSS.Stop();
        tasscorr.Start();        

        // subtract terms from fmm:
        
        MappedIntegrationRule<2,3> trial_mir(ir, trial_trafo, lh);
        MappedIntegrationRule<2,3> test_mir(ir, test_trafo, lh);
          
        FlatMatrix<> shapesi(test_fel.GetNDof(), test_evaluator->Dim()*ir.Size(), lh);
        FlatMatrix<> shapesj(trial_fel.GetNDof(), trial_evaluator->Dim()*ir.Size(), lh);
          
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
        tasscorr.Stop();        
        nearfield_correction -> AddElementMatrix (test_dnums, trial_dnums, elmat);
      }
    */

    TableCreator<int> create_nbels;
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
            tassSS.Start();
            CalcElementMatrix (elmat, ei_trial, ei_test, lh);
            tassSS.Stop();
            tasscorr.Start();        
            
            // subtract terms from fmm:
            
            MappedIntegrationRule<2,3> trial_mir(ir, trial_trafo, lh);
            MappedIntegrationRule<2,3> test_mir(ir, test_trafo, lh);
            
            FlatMatrix<> shapesi(test_fel.GetNDof(), test_evaluator->Dim()*ir.Size(), lh);
            FlatMatrix<> shapesj(trial_fel.GetNDof(), trial_evaluator->Dim()*ir.Size(), lh);
            
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
            tasscorr.Stop();        
            nearfield_correction -> AddElementMatrix (test_dnums, trial_dnums, elmat);
          }
      });

    
    tassemble.Stop();
    return TransposeOperator(evaly) * fmmop * evalx + nearfield_correction;
  }
  


  /*
  template <typename T>    
  void IntegralOperator<T> :: CalcHMatrix(HMatrix<T> & hmatrix, LocalHeap & clh, struct BEMParameters &param) const
  {
    static Timer t("ngbem - BaseClass::CalcHMatrix"); RegionTimer reg(t);    
    auto & matList = hmatrix.GetMatList();

    ParallelForRange (matList.Size(), [&](IntRange r)
    {
      LocalHeap lh = clh.Split();
      for (int k : r)
        {
          HeapReset hr(lh);
          BEMBlock<T> & block = matList[k];
          auto trialdofs = block.GetTrialDofs();
          auto testdofs = block.GetTestDofs();
          if(block.IsNearField())
            {
              // Compute dense block
              Matrix<T> mat_near(testdofs.Size(), trialdofs.Size());
              CalcBlockMatrix(mat_near, trialdofs, testdofs, lh);
              block.SetMat(make_unique<BaseMatrixFromMatrix<T>>(std::move(mat_near)));
            }
          else
            {
              // Compute low-rank block
              try
                {
                  block.SetMat(CalcFarFieldBlock(trialdofs, testdofs, lh));
                }
              catch (netgen::NgException & e)
                {
                  // cout << "not seperated, size = " << testdofs.Size() << " x " << trialdofs.Size() << endl;
                  Matrix<T> mat_near(testdofs.Size(), trialdofs.Size());
                  CalcBlockMatrix(mat_near, trialdofs, testdofs, lh);
                  block.SetMat(make_unique<BaseMatrixFromMatrix<T>>(std::move(mat_near)));
                }
            }
        }
    }, TasksPerThread(4));
  }


  template <typename T>
  void StochasticTSVD1 (MatrixView<T> A, MatrixView<T> U, MatrixView<T> V, VectorView<> S)
  {
    for (int i = 0; i < V.Height(); i++)
      for (int j = 0; j < V.Width(); j++)
        V(i,j) = double (rand()) / RAND_MAX;

    Matrix<T> AVt = A * Trans(V);
    Matrix<T> tmp(V.Height(), V.Height());
    LapackSVD (AVt, U, tmp, S, false);    
    Matrix<T> UtA = Trans(U) * A;
    LapackSVD (UtA, tmp, V, S, false);
    U = Matrix<T> (U * tmp);
  }

  void StochasticTSVD1 (MatrixView<Complex> A, MatrixView<Complex> U, MatrixView<Complex> V, VectorView<> S)
  {
    for (int i = 0; i < V.Height(); i++)
      for (int j = 0; j < V.Width(); j++)
        V(i,j) = double (rand()) / RAND_MAX;

    // Matrix<Complex> AVt = A * Conj(Trans(V));
    Matrix<Complex> AVt = A * Trans(V) | Lapack;
    
    Matrix<Complex> tmp(V.Height(), V.Height());
    LapackSVD (AVt, U, tmp, S, false);    
    // Matrix<Complex> UtA = Conj(Trans(U)) * A;
    Matrix<Complex> UtA = Matrix<Complex>(Conj(Trans(U))) * A | Lapack;
    LapackSVD (UtA, tmp, V, S, false);
    U = Matrix<Complex> (U * tmp | Lapack);
  }

  
  template <typename T>
  size_t StochasticTSVD (MatrixView<T> A, MatrixView<T> U, MatrixView<T> V, VectorView<> S, double eps)
  {
    static Timer tsvd("ngbem - StochasticTSVD"); 
    RegionTimer reg(tsvd);
    
    int rank = 5;
    int p = min(A.Height(), A.Width());
    while (rank < p)
      {
        StochasticTSVD1 (A, U.Cols(rank), V.Rows(rank), S.Range(rank));
        if (S[rank-1] < eps)
          {
            for (int j = 1; j < p; j++)
              if (S[j] < eps)
                return j-1;
          }

        rank = int (rank*1.5);
      }
    
    LapackSVD(A, U, V, S, false);
    for (int j = 1; j < p; j++)
      if (S[j] < eps)
        return j-1;
    return p;
  }
  */


  /*
  template <typename T>    
  unique_ptr<LowRankMatrix<T>> IntegralOperator<T> ::
  CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, LocalHeap &lh) const
  {
    static Timer t("ngbem - IntegralOperator::CalcFarFieldBlock"); RegionTimer reg(t);
    int m = testdofs.Size();
    int n = trialdofs.Size();
    int p = min(n, m);
    
    Matrix<T> A(m, n);
    CalcBlockMatrix(A, trialdofs, testdofs, lh);

    // Calculate SVD for A^\top = V S U^\top
    Matrix<T, ColMajor> V(n, p), Ut(p, m);
    Vector<> S(p);

    int k = StochasticTSVD<T> (A, Trans(Ut), Trans(V), S, param.eps);

    // Low-rank approximation from truncated svd
    
    Matrix<T> U_trc(m, k), Vt_trc(k, n);
    for (size_t j = 0; j < k; j++)
      U_trc.Col(j) = sqrt(S(j)) * Ut.Row(j);

    for (size_t i = 0; i < k; i++)    
      Vt_trc.Row(i) = sqrt(S(i)) * V.Col(i);

    return make_unique<LowRankMatrix<T>> (Matrix<T>(U_trc), Matrix<T>(Vt_trc));
  }
  */
  
  

  /*
  template <typename KERNEL>
  GenericIntegralOperator<KERNEL> ::
  GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                          KERNEL _kernel,
                          BEMParameters _param)
    : IntegralOperator<value_type>(_trial_space, _test_space, _param), kernel(_kernel)
  {
    trial_evaluator = trial_space -> GetEvaluator(BND);
    test_evaluator = test_space -> GetEvaluator(BND);
    hmatrix =
      make_shared<HMatrix<value_type>>(trial_ct, test_ct, 
                                       param.eta, trial_space->GetNDof(), test_space->GetNDof());
    
    LocalHeap lh(100000000);
    this->CalcHMatrix(*hmatrix, lh, param);
	
    if (param.testhmatrix)
      {
        Matrix<value_type> dense(test_space->GetNDof(), trial_space->GetNDof());
        CalcBlockMatrix(dense, mapbnd2glob, mapbnd2glob2, lh);            
        cout << "dense: " << dense.Height() << " x " << dense.Width() << endl;
        
        // compute all its blocks
        HeapReset hr(lh);    
        
        // Estimate spectral norm error
        Vector<value_type> x(trial_space->GetNDof()), y(test_space->GetNDof());
	VFlatVector<value_type> x_base(x), y_base(y);
	x = 1.;

	for (int i = 0; i < 30; i++) {
	  x *= 1. / L2Norm(x);
	  y = dense * x;	  
	  y_base -= (*hmatrix) * x_base;
	  y = Conj(y);
	  x = Trans(dense) * y;
	  hmatrix->MultTransAdd(-1., y_base, x_base);
	  x = Conj(x);
	}

	// Get compression rate
	auto & matList = hmatrix->GetMatList();
	size_t nf_c = 0.;
	for (int i = 0; i < matList.Size(); i++) {
	  BEMBlock<value_type> & block = matList[i];
	  if (block.IsNearField()) {
	    nf_c += block.GetTrialDofs().Size() * block.GetTestDofs().Size();
	  } else {
	    nf_c += block.GetMat()->NZE();
	  }
	}
										  	    
	cout << "compression rate " << nf_c / ((double) test_space->GetNDof() * trial_space->GetNDof()) << "  2-norm error " << sqrt(L2Norm(x)) << endl;
      }
  }
  */


  template <typename KERNEL>
  GenericIntegralOperator<KERNEL> ::
  GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                          optional<Region> _definedon_trial, optional<Region> _definedon_test,                          
                          shared_ptr<DifferentialOperator> _trial_evaluator, 
                          shared_ptr<DifferentialOperator> _test_evaluator, 
                          KERNEL _kernel,
                          BEMParameters _param)
  : IntegralOperator<value_type>(_trial_space, _test_space, _definedon_trial, _definedon_test, _param), kernel(_kernel),
    trial_evaluator(_trial_evaluator), test_evaluator(_test_evaluator)
  {
    LocalHeap lh(100000000);

    tie(identic_panel_x, identic_panel_y, identic_panel_weight) = IdenticPanelIntegrationRule(param.intorder);
    tie(common_vertex_x, common_vertex_y, common_vertex_weight) = CommonVertexIntegrationRule(param.intorder);
    tie(common_edge_x, common_edge_y, common_edge_weight) = CommonEdgeIntegrationRule(param.intorder);
    
    if (param.method == "fmm")
      {
        matrix = this->CreateMatrixFMM(lh);
        return;
      }

    /*
    auto hmatrix =
      make_shared<HMatrix<value_type>>(trial_ct, test_ct, 
                                       param.eta, trial_space->GetNDof(), test_space->GetNDof());
    this->matrix = hmatrix;
    this->CalcHMatrix(*hmatrix, lh, param);
    */

    /*
    if (param.testhmatrix)
      {
        Matrix<value_type> dense(test_space->GetNDof(), trial_space->GetNDof());
        CalcBlockMatrix(dense, mapbnd2glob, mapbnd2glob2, lh);            
        cout << "dense: " << dense.Height() << " x " << dense.Width() << endl;
        
        // compute all its blocks
        HeapReset hr(lh);    
        
        // Estimate spectral norm error
        Vector<value_type> x(trial_space->GetNDof()), y(test_space->GetNDof());
	VFlatVector<value_type> x_base(x), y_base(y);
	x = 1.;

	for (int i = 0; i < 30; i++) {
	  x *= 1. / L2Norm(x);
	  y = dense * x;	  
	  y_base -= (*hmatrix) * x_base;
	  y = Conj(y);
	  x = Trans(dense) * y;
	  hmatrix->MultTransAdd(-1., y_base, x_base);
	  x = Conj(x);
	}

	// Get compression rate
	auto & matList = hmatrix->GetMatList();
	size_t nf_c = 0.;
	for (int i = 0; i < matList.Size(); i++) {
	  BEMBlock<value_type> & block = matList[i];
	  if (block.IsNearField()) {
	    nf_c += block.GetTrialDofs().Size() * block.GetTestDofs().Size();
	  } else {
	    nf_c += block.GetMat()->NZE();
	  }
	}
										  	    
	cout << "compression rate " << nf_c / ((double) test_space->GetNDof() * trial_space->GetNDof()) << "  2-norm error " << sqrt(L2Norm(x)) << endl;
      }
    */
  }

  
  template <typename KERNEL>
  void GenericIntegralOperator<KERNEL> ::
  CalcBlockMatrix(FlatMatrix<value_type> matrix, FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, 
                  LocalHeap &lh) const
  {
    auto mesh = this->trial_space->GetMeshAccess();  
    auto mesh2 = this->test_space->GetMeshAccess();  
    
    static Timer tall("ngbem - all " + KERNEL::Name());
    static Timer tloops("ngbem - loops " + KERNEL::Name());
    static Timer t_identic("ngbem identic panel " + KERNEL::Name());
    static Timer t_common_vertex("ngbem common vertex " + KERNEL::Name());

    static Timer t_common_edge("ngbem common edge " + KERNEL::Name());
    static Timer t_disjoint("ngbem disjoint " + KERNEL::Name());

    RegionTimer reg(tall);

    IntegrationRule irtrig(ET_TRIG, param.intorder);
    
    auto [ identic_panel_x, identic_panel_y, identic_panel_weight ] =
      IdenticPanelIntegrationRule(param.intorder);

    auto [ common_vertex_x, common_vertex_y, common_vertex_weight ] =
      CommonVertexIntegrationRule(param.intorder);
    
    auto [ common_edge_x, common_edge_y, common_edge_weight ] =
      CommonEdgeIntegrationRule(param.intorder);

    matrix = 0; 

    // auto evaluator = trial_space->GetEvaluator(BND);
    // auto evaluator2 = test_space->GetEvaluator(BND);

    Array<int> tmp, tmp2;
    Array<int> patchi, patchj;
    Array<int> trialdofsinv(trial_space->GetNDof()); 
    Array<int> testdofsinv(test_space->GetNDof());

    trialdofsinv = -1;
    testdofsinv = -1;
  
    for (int i = 0; i < testdofs.Size(); i++)
      {
	testdofsinv[testdofs[i]] = i;
	tmp2.Append(elems4dof2[ testdofs[i] ]);
      }
    QuickSort( tmp2 );
    for (int i = 0; i < tmp2.Size(); i++)
      {
	patchi.Append(tmp2[i]);
	int tmpi = tmp2[i];
	while (i < tmp2.Size() && tmp2[i] == tmpi)
	  i++;
	i--;
      }
    
    for (int j = 0; j < trialdofs.Size(); j++)
      {
	trialdofsinv[trialdofs[j]] = j;
	tmp.Append(elems4dof[ trialdofs[j] ]);
      }
    QuickSort( tmp );
    for (int j = 0; j < tmp.Size(); j++)
      {
	patchj.Append(tmp[j]);
	int tmpj = tmp[j];
	while (j < tmp.Size() && tmp[j] == tmpj)
	  j++;
	j--;
      }

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
      
      test_evaluator->CalcMatrix(feli, mirx, mshapesi);
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

      
    RegionTimer regloops(tloops);    
    for (int i = 0; i < patchi.Size(); i++) // test
      for (int j = 0; j < patchj.Size(); j++) // trial
	{
	  HeapReset hr(lh);
	  ElementId ei(BND, patchi[i]);
	  ElementId ej(BND, patchj[j]);
          
	  auto verti = mesh2->GetElement(ei).Vertices();
	  auto vertj = mesh->GetElement(ej).Vertices();          
            
	  FiniteElement &feli = test_space->GetFE(ei, lh);
	  FiniteElement &felj = trial_space->GetFE(ej, lh);
              
	  ElementTransformation &trafoi = mesh2->GetTrafo(ei, lh);
	  ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
              
	  Array<DofId> dnumsi, dnumsj;
	  test_space->GetDofNrs(ei, dnumsi); // mapping to global dof
	  trial_space->GetDofNrs(ej, dnumsj);
        
	  FlatMatrix<> shapei(feli.GetNDof(), test_evaluator->Dim(), lh);
	  FlatMatrix<> shapej(felj.GetNDof(), trial_evaluator->Dim(), lh);
          
	  FlatMatrix<value_type> elmat(feli.GetNDof(), felj.GetNDof(), lh); 
	  elmat = 0;
          
	  int n_common_vertices = 0;
	  for (auto vi : verti)
	    if (vertj.Contains(vi))
	      n_common_vertices++;


	  switch (n_common_vertices)
	    {
	    case 3: //identical panel
	      {
                RegionTimer reg(t_identic);    
                      
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

                    Integrate4D (irx, iry, feli, felj, trafoi, trafoj, elmat, lh);
                  }

                
		break;
	      }
	    case 2: //common edge
	      {
                RegionTimer reg(t_common_edge);    
          
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

                    Integrate4D (irx, iry, feli, felj, trafoi, trafoj, elmat, lh);
                  }
                
		break;
	      }

	    case 1: //common vertex
	      {
                RegionTimer reg(t_common_vertex);    

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

                    Integrate4D (irx, iry, feli, felj, trafoi, trafoj, elmat, lh);
                  }

		break;
	      }

	    case 0: //disjoint panels
	      {
                RegionTimer r(t_disjoint);    
                  
		// shapes+geom out of loop, matrix multiplication
		MappedIntegrationRule<2,3> mirx(irtrig, trafoi, lh);
		MappedIntegrationRule<2,3> miry(irtrig, trafoj, lh);

		FlatMatrix<> shapesi(feli.GetNDof(), test_evaluator->Dim()*irtrig.Size(), lh);
		FlatMatrix<> shapesj(felj.GetNDof(), trial_evaluator->Dim()*irtrig.Size(), lh);
                
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
                    elmat += shapesi1 * kernel_shapesj;
                  }
		break;
	      }
	    default:
	      throw Exception ("not possible");
	    }

	  for (int ii = 0; ii < dnumsi.Size(); ii++) // test
	    for (int jj = 0; jj < dnumsj.Size(); jj++) // trial
	      if(trialdofsinv[dnumsj[jj]] != -1 && testdofsinv[dnumsi[ii]] != -1)
		matrix(testdofsinv[dnumsi[ii]], trialdofsinv[dnumsj[jj]]) += elmat(ii, jj);
	}
  }


  template <typename KERNEL>
  void GenericIntegralOperator<KERNEL> ::
  CalcElementMatrix(FlatMatrix<value_type> matrix,
                    ElementId ei_trial, ElementId ei_test,
                    LocalHeap &lh) const
  {
    auto mesh = this->trial_space->GetMeshAccess();  
    auto mesh2 = this->test_space->GetMeshAccess();  
    
    static Timer tall("ngbem - elementmatrix " + KERNEL::Name());
    RegionTimer reg(tall);

    static Timer t1("ngbem - elementmatrix, part1  " + KERNEL::Name());

    t1.Start();
    IntegrationRule irtrig(ET_TRIG, param.intorder);
    /*
    auto [ identic_panel_x, identic_panel_y, identic_panel_weight ] =
      IdenticPanelIntegrationRule(param.intorder);

    auto [ common_vertex_x, common_vertex_y, common_vertex_weight ] =
      CommonVertexIntegrationRule(param.intorder);
    
    auto [ common_edge_x, common_edge_y, common_edge_weight ] =
      CommonEdgeIntegrationRule(param.intorder);
    */
    matrix = 0; 
    t1.Stop();

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
      
      test_evaluator->CalcMatrix(feli, mirx, mshapesi);
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

  

#ifdef USE_HMATRIX
  template <typename KERNEL>
  unique_ptr<LowRankMatrix<typename KERNEL::value_type>> GenericIntegralOperator<KERNEL> ::    
  CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
                    LocalHeap &lh) const 
  
  {
    // if (trial_evaluator->Dim() > 1)
    // throw Exception("ACA not supported for vectorial evaluators");              
    
    auto mesh = this->trial_space->GetMeshAccess();  
    auto mesh2 = this->test_space->GetMeshAccess();  
    
    static Timer tall("ngbem FarFieldBlock " + KERNEL::Name());
    static Timer tACA("ngbem FarFieldBlock - ACA " + KERNEL::Name());
    static Timer tkernel("ngbem FarFieldBlock - kernel " + KERNEL::Name());
    static Timer tnorm("ngbem FarFieldBlock - norm " + KERNEL::Name());
    static Timer tmatvec1("ngbem FarFieldBlock - matvec1 " + KERNEL::Name());              
    static Timer tmatvec2("ngbem FarFieldBlock - matvec2 " + KERNEL::Name());              
    RegionTimer reg(tall);

    IntegrationRule irtrig(ET_TRIG, param.intorder);
    SIMD_IntegrationRule simd_irtrig(irtrig);

    Array<int> tmp, tmp2;
    Array<int> patchi, patchj;
    Array<int> trialdofsinv(trial_space->GetNDof()); 
    Array<int> testdofsinv(test_space->GetNDof());

    trialdofsinv = -1;
    testdofsinv = -1;
  
    for (int i = 0; i < testdofs.Size(); i++)
      {
	testdofsinv[testdofs[i]] = i;
	tmp2.Append(elems4dof2[ testdofs[i] ]);
      }
    QuickSort( tmp2 );
    for (int i = 0; i < tmp2.Size(); i++)
      {
	patchi.Append(tmp2[i]);
	int tmpi = tmp2[i];
	while (i < tmp2.Size() && tmp2[i] == tmpi)
	  i++;
	i--;
      }
    
    for (int j = 0; j < trialdofs.Size(); j++)
      {
	trialdofsinv[trialdofs[j]] = j;
	tmp.Append(elems4dof[ trialdofs[j] ]);
      }
    QuickSort( tmp );
    for (int j = 0; j < tmp.Size(); j++)
      {
	patchj.Append(tmp[j]);
	int tmpj = tmp[j];
	while (j < tmp.Size() && tmp[j] == tmpj)
	  j++;
	j--;
      }



    // new code
    Array<Vec<3>> xi, yj, nxi, nyj;  // i..test, j... trial
    Array<double> wxi, wyj;
    
    BitArray test_vertices(mesh->GetNV());
    test_vertices.Clear();

    // test patches
    for (int i = 0; i < patchi.Size(); i++)
      {
        HeapReset hr(lh);
        ElementId ei(BND, patchi[i]);
          
        for (auto v : mesh2->GetElement(ei).Vertices())
          test_vertices.SetBit(v);
        
        ElementTransformation &trafoi = mesh2->GetTrafo(ei, lh);
        MappedIntegrationRule<2,3> mirx(irtrig, trafoi, lh);
        
        for (int ix = 0; ix < irtrig.Size(); ix++)
          {
            xi.Append (mirx[ix].GetPoint());
            nxi.Append (mirx[ix].GetNV());
            wxi.Append (mirx[ix].GetWeight());
          }
      }

    // trial patches
    for (int j = 0; j < patchj.Size(); j++)
      {
        HeapReset hr(lh);
        ElementId ej(BND, patchj[j]);

        if (mesh == mesh2)
          for (auto v : mesh->GetElement(ej).Vertices())
            if (test_vertices.Test(v))
              throw Exception("far field block must not have common vertices");                     
        ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
        MappedIntegrationRule<2,3> miry(irtrig, trafoj, lh);
            
        for (int iy = 0; iy < irtrig.Size(); iy++)
          {
            yj.Append (miry[iy].GetPoint());
            nyj.Append (miry[iy].GetNV());
            wyj.Append (miry[iy].GetWeight());
          }
      }

    tACA.Start();
    /*
    Matrix<value_type> kernel_matrix(xi.Size(), yj.Size());
    for (int i = 0; i < xi.Size(); i++)
      for (int j = 0; j < yj.Size(); j++)
        kernel_matrix(i,j) = kernel.Evaluate(xi[i], yj[j], nxi[i], nyj[j])(0,0);

    size_t p = 200;
    Matrix<value_type> Umax(kernel_matrix.Height(), p), Vmax(p, kernel_matrix.Width());
    Vector<> S(p);

    int k = StochasticTSVD<value_type> (kernel_matrix, Umax, Vmax, S, param.eps);


    for (size_t j = 0; j < k; j++)
      {
        Umax.Col(j) *= sqrt(S(j));
        Vmax.Row(j) *= sqrt(S(j));
      }
    
    */

    size_t p = min(xi.Size(), yj.Size());
    //int rank = p;
    auto GetRow = [&](int i, SliceVector<value_type> row, int comp)
    {
      RegionTimer reg(tkernel);
      tkernel.AddFlops (yj.Size());
      for (int j = 0; j < yj.Size(); j++)
        row(j) = kernel.Evaluate(xi[i], yj[j], nxi[i], nyj[j])(comp);
    };
    auto GetCol = [&](int j, SliceVector<value_type> col, int comp)
    {
      RegionTimer reg(tkernel);
      tkernel.AddFlops (xi.Size());
      for (int i = 0; i < xi.Size(); i++)
        col(i) = kernel.Evaluate(xi[i], yj[j], nxi[i], nyj[j])(comp);
    };


    Array<Matrix<value_type>> Us;
    Array<Matrix<value_type>> Vs;

    size_t num_kernel_comps = 0;
    for (auto term : kernel.terms)
      num_kernel_comps = std::max(num_kernel_comps, term.kernel_comp+1);
    

    for (size_t comp = 0; comp < num_kernel_comps; comp++)
      {
        Matrix<value_type> Umax(xi.Size(), 20);
        Matrix<value_type> Vmax(20, yj.Size());
    
        size_t ik = 0, jk = 0, ikm1 = 0, jkm1 = yj.Size() + 1, rank = p;
        
        // // for quasi random sequence of pivot indices
        // size_t primes[] = { 71, 73, 79, 83, 89, 97 };
        // size_t prime;
        // for (auto tp : primes)
        //   {
        //     if (xi.Size()%tp != 0)
        //       {
        //         prime = tp;
        //         break;
        //       }
        //   }
        // // ACA compression 
        // for (size_t k = 0; k < p; k++)
        //   {
        //     // int ik = k;  // what else ?
        //     size_t ik = (k*prime)%xi.Size(); 
        
        //     GetRow(ik, Vmax.Row(k));
        //     Vmax.Row(k) -= Trans(Vmax.Rows(0,k)) * Umax.Row(ik).Range(0,k);
        
        //     double err = L2Norm(Vmax.Row(k));
        //     // cout << "Norm vk = " << err << endl;
        //     if (err < param.eps)
        //       {
        //         rank = k;
        //         break;
        //       }
        
        //     int jmax = 0;
        //     for (int j = 0; j < Vmax.Width(); j++)
        //       if (fabs (Vmax(k,j)) > fabs(Vmax(k,jmax)))
        //         jmax = j;
        //     Vmax.Row(k) *= 1.0 / Vmax(k,jmax);
        
        //     GetCol(jmax, Umax.Col(k));
        //     Umax.Col(k) -= Umax.Cols(0,k) * Vmax.Col(jmax).Range(0,k);
        //   }

        
        // Scale eps appropriately (see Bebendorf, Hierarchical Matrices  p. 126 & 135
        //double eps = 2. / 3. * param.eps / sqrt(xi.Size() * yj.Size());
	double eps = param.eps;
        // The Frobenius norm squared of the approximant U * V^H
        double norm2 = 0., norm2_V = 0., norm2_U = 0.;
        
        for (size_t k = 0; k < p; k++) {
          // Get the ik-th row
          if (Umax.Width() == k)
            {
              Umax.ChangeSize(Umax.Height(), Umax.Width()+20);
              Vmax.ChangeSize(Vmax.Height()+20, Vmax.Width());
            }
          
          GetRow(ik, Vmax.Row(k), comp);
          tmatvec1.Start();          
          Vmax.Row(k) -= Trans(Vmax.Rows(0,k)) * Umax.Row(ik).Range(0,k);
          tmatvec1.Stop();
          tmatvec1.AddFlops (k * Vmax.Width());
           
          // Find the new column pivot position jk in the new row
          double vkj = 0.;
          for (int j = 0; j < Vmax.Width(); j++)
            if (fabs (Vmax(k, j)) > vkj && j != jkm1) {
	      vkj = fabs (Vmax(k, j));
	      jk = j;
	    }
          
          // If the pivot element is close to zero, exit
          if (vkj == 0.) {
            rank = k;
            break;
          }
          
          // Scale with inverse of the pivot entry at (ik, jk)
          Vmax.Row(k) *= sqrt(vkj) / Vmax(k, jk);
          
          // Get the jk-th column
          GetCol(jk, Umax.Col(k), comp);
          tmatvec2.Start();
          Umax.Col(k) -= Umax.Cols(0,k) * Vmax.Col(jk).Range(0,k);
	  Umax.Col(k) *= 1. / sqrt(vkj);
          tmatvec2.Stop();
          tmatvec2.AddFlops (k * Umax.Height());          
          
          // Find the new row pivot position ik in the new column
          double uik = 0.;
          for (int i = 0; i < Umax.Height(); i++)
            if (fabs (Umax(i, k)) > uik && i != ikm1) {
	      uik = fabs (Umax(i, k));
	      ik = i;
	    }

          tnorm.Start();

          // Update the Frobenius norm
          double norm_k = L2Norm(Vmax.Row(k)) * L2Norm(Umax.Col(k));
          norm2 += norm_k * norm_k;
          
	  // how it was:
          // for (int l = 0; l < k; l++)
          //   norm2 += 2. * std::real(InnerProduct(Vmax.Row(k), Vmax.Row(l)) *
          //                           InnerProduct(Umax.Col(k), Umax.Col(l)));

          // the correct complex inner product:
          for (int l = 0; l < k; l++)
            norm2 += 2. * std::real(Conj(InnerProduct(Vmax.Row(k), Conj(Vmax.Row(l)))) *
                                    InnerProduct(Umax.Col(k), Conj(Umax.Col(l))));

	  // // Estimate the spectral norm instead of Frobenius norm (likely overestimating)
	  // double norm_Vk = L2Norm(Vmax.Row(k)), norm_Uk = L2Norm(Umax.Col(k));	  
	  // norm2_V += norm_Vk * norm_Vk;
	  // norm2_U += norm_Uk * norm_Uk;
	  
          tnorm.Stop();
          //  *testout << "k = " << k << "norm2 = " << norm2 << endl;
          // *testout << "norm2 = " << norm2 << " =?= " << L2Norm2(Vmax.Rows(k+1))*L2Norm2(Umax.Cols(k+1)) << endl;
           
          // New pivots become old pivots
          ikm1 = ik;
          jkm1 = jk;

	  // // Stop if the updates are separately relatively small
	  // if (norm_Vk < sqrt(eps) * sqrt(norm2_V) && norm_Uk < sqrt(eps) * sqrt(norm2_U)) {
	  //   rank = k + 1;
	  //   break;
	  // }
          // Stop if the new update is relatively small, see Bebendorf pp. 141-142
          if (norm_k < eps * sqrt(norm2)) {
	    rank = k + 1;
	    break;
	  }
        }

        Us += Matrix<value_type> (Umax.Cols(0,rank));
        Vs += Matrix<value_type> (Vmax.Rows(0,rank));

        for (int i = 0; i < Us.Last().Height(); i++)
          Us.Last().Row(i) *= wxi[i];
        for (int j = 0; j < Vs.Last().Width(); j++)
          Vs.Last().Col(j) *= wyj[j];
      }

    Array<IntRange> ranges;
    size_t total_rank = 0;
    for (auto term : kernel.terms)
      {
        ranges += IntRange(total_rank, total_rank+Us[term.kernel_comp].Width());
        total_rank += Us[term.kernel_comp].Width();
      }
    
    // *testout << "rank = " << rank << endl;
    // size_t k = rank;
    tACA.Stop();
    
    // auto U = Umax.Cols(0,k);
    // auto V = Vmax.Rows(0,k);
    // cout << "k = " << k << ", err = " << L2Norm(help-U*V) << endl;


    Matrix<value_type> U2(testdofs.Size(), total_rank); 
    Matrix<value_type> V2(total_rank, trialdofs.Size());
    U2 = value_type(0.0);
    V2 = value_type(0.0);
    
    int cnt = 0;
    for (int i = 0; i < patchi.Size(); i++) // test
      {
        HeapReset hr(lh);
        ElementId ei(BND, patchi[i]);
            
        FiniteElement &feli = test_space->GetFE(ei, lh);
        ElementTransformation &trafoi = mesh2->GetTrafo(ei, lh);
              
        Array<DofId> dnumsi;
        test_space->GetDofNrs(ei, dnumsi);

        int dim = test_evaluator->Dim();
        SIMD_MappedIntegrationRule<2,3> mirx(simd_irtrig, trafoi, lh);
        FlatMatrix<SIMD<double>> shapesi(dim*feli.GetNDof(), simd_irtrig.Size(), lh);
        SliceMatrix<double> dshapesi(dim*feli.GetNDof(), simd_irtrig.GetNIP(), simd_irtrig.Size()*SIMD<double>::Size(),
                                     (double*)shapesi.Data());
        
        test_evaluator -> CalcMatrix(feli, mirx, shapesi);

        /*
        Matrix<value_type> tmp1 = dshapesi * U.Rows(cnt, cnt+irtrig.Size());        
        auto tmp = tmp1.Reshape(feli.GetNDof(), dim*U.Width());
                                
        cnt += irtrig.Size();
        for (int ii = 0; ii < dnumsi.Size(); ii++) // test
          if (testdofsinv[dnumsi[ii]] != -1)
            U2.Row(testdofsinv[dnumsi[ii]]) += tmp.Row(ii);
        */

        // for (auto term : kernel.terms)
        for (int t = 0; t < kernel.terms.Size(); t++)
          {
            auto term = kernel.terms[t];
            auto U2cols = U2.Cols(ranges[t]); // term.test_comp*k, (term.test_comp+1)*k);
            Matrix<value_type> tmp = dshapesi.RowSlice(term.test_comp, dim).AddSize(feli.GetNDof(), simd_irtrig.GetNIP())
              * Us[term.kernel_comp].Rows(cnt, cnt+irtrig.Size());        
            for (int ii = 0; ii < dnumsi.Size(); ii++) // test
              if (testdofsinv[dnumsi[ii]] != -1)
                U2cols.Row(testdofsinv[dnumsi[ii]]) += term.fac * tmp.Row(ii);
          }
        cnt += irtrig.Size();
        
      }

    cnt = 0;
    for (int j = 0; j < patchj.Size(); j++) // test
      {
        HeapReset hr(lh);
        ElementId ej(BND, patchj[j]);
            
        FiniteElement &felj = trial_space->GetFE(ej, lh);
        ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
        
        Array<DofId> dnumsj;
        trial_space->GetDofNrs(ej, dnumsj);


        int dim = trial_evaluator->Dim();        
        SIMD_MappedIntegrationRule<2,3> miry(simd_irtrig, trafoj, lh);
        FlatMatrix<SIMD<double>> shapesj(dim*felj.GetNDof(), simd_irtrig.Size(), lh);
        SliceMatrix<double> dshapesj(dim*felj.GetNDof(), simd_irtrig.GetNIP(), simd_irtrig.Size()*SIMD<double>::Size(),
                                     (double*)shapesj.Data());
        
        trial_evaluator -> CalcMatrix(felj, miry, shapesj);
        /*
        Matrix<value_type> tmp1 = dshapesj * Trans(V).Rows(cnt, cnt+irtrig.Size());        
        auto tmp = tmp1.Reshape(felj.GetNDof(), dim*V.Height());
        
        cnt += irtrig.Size();
        for (int jj = 0; jj < dnumsj.Size(); jj++) // trial
          if(trialdofsinv[dnumsj[jj]] != -1)
            V2.Col(trialdofsinv[dnumsj[jj]]) += tmp.Row(jj);
        */

        for (int t = 0; t < kernel.terms.Size(); t++)
          {
            auto term = kernel.terms[t];
            auto V2rows = V2.Rows(ranges[t]); // term.trial_comp*k, (term.trial_comp+1)*k);
            Matrix<value_type> tmp = dshapesj.RowSlice(term.trial_comp, dim).AddSize(felj.GetNDof(), simd_irtrig.GetNIP())
              * Trans(Vs[term.kernel_comp]).Rows(cnt, cnt+irtrig.Size());        
            for (int jj = 0; jj < dnumsj.Size(); jj++) // trial
              if (trialdofsinv[dnumsj[jj]] != -1)
                V2rows.Col(trialdofsinv[dnumsj[jj]]) += tmp.Row(jj);
          }
        cnt += irtrig.Size();
        
      }

    return make_unique<LowRankMatrix<value_type>> (std::move(U2), std::move(V2));
  }
#endif
  
  
  template <typename KERNEL>
  shared_ptr<CoefficientFunction> GenericIntegralOperator<KERNEL> ::
  GetPotential(shared_ptr<GridFunction> gf) const
  {
    return  make_shared<PotentialCF<KERNEL>> (gf, trial_definedon, trial_evaluator,
                                              kernel, param);
  }


  template <typename KERNEL>
  PotentialCF<KERNEL> ::
  PotentialCF (shared_ptr<GridFunction> _gf,
               optional<Region> _definedon,                   
               shared_ptr<DifferentialOperator> _evaluator,
               KERNEL _kernel, BEMParameters _param)
    : CoefficientFunctionNoDerivative (_evaluator->Dim(), std::is_same<typename KERNEL::value_type,Complex>()),
      gf(_gf), definedon(_definedon), evaluator(_evaluator), kernel(_kernel), param(_param)
  {
    ;
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
          IntegrationRule ir(fel.ElementType(), param.intorder);
          SIMD_IntegrationRule simd_ir(ir);
          SIMD_MappedIntegrationRule<2,3> miry(simd_ir, trafo, lh);
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
    for (int i = 0; i < Dimension(); i++)
      result(i) = HSum(simd_result(i));
  }


  template <typename KERNEL> template <typename T>
  void PotentialCF<KERNEL> :: T_Evaluate(const BaseMappedIntegrationRule & bmir,
                                         BareSliceMatrix<T> result) const
  {
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
              
              IntegrationRule ir(fel.ElementType(), param.intorder);
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
    // cout << "Potential - eval (simd) called" << endl;
    throw ExceptionNOSIMD ("PotentialCF::Evaluate (SIMD) not available");
    
    result.AddSize(Dimension(), ir.Size()) = SIMD<T>(0.0);
    return;
    /*
    LocalHeapMem<100000> lh("Potential::Eval");
    auto space = this->gf->GetFESpace();
    auto mesh = space->GetMeshAccess();
    
    result = T(0.0);
    auto & mip23 = dynamic_cast<const MappedIntegrationPoint<2,3>&>(mip);
    if constexpr (std::is_same<typename KERNEL::value_type,T>())
      for (size_t i = 0; i < mesh->GetNSE(); i++)
        {
          HeapReset hr(lh);
          ElementId ei(BND, i);
          
          const FiniteElement &fel = space->GetFE(ei, lh);
          const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);
          
          Array<DofId> dnums(fel.GetNDof(), lh);
          space->GetDofNrs(ei, dnums);
          FlatVector<T> elvec(fel.GetNDof(), lh);
          gf->GetElementVector(dnums, elvec);
          
          IntegrationRule ir(fel.ElementType(), param.intorder);
          MappedIntegrationRule<2,3> mir(ir, trafo, lh);
          FlatMatrix<T> vals(ir.Size(), evaluator->Dim(), lh);
          
          evaluator->Apply (fel, mir, elvec, vals, lh);
          
          for (int j = 0; j < ir.Size(); j++)
            {
              Vec<3> x = mir[j].GetPoint();
              Vec<3> y = mip23.GetPoint();
              
              Vec<3> nx = mir[j].GetNV();
              Vec<3> ny = mip23.GetNV();
              
              for (auto term : kernel.terms)
                {
                  auto kernel_ = kernel.Evaluate(x, y, nx, ny)(term.kernel_comp);
                  result(term.test_comp) += mir[j].GetWeight()*kernel_ * vals(j, term.trial_comp);
                }
            }
        }
    */
  }
  
  
  template class GenericIntegralOperator<LaplaceSLKernel<3>>;
  template class GenericIntegralOperator<LaplaceDLKernel<3>>;
  template class GenericIntegralOperator<LaplaceHSKernel<3>>;
  
  template class GenericIntegralOperator<HelmholtzSLKernel<3>>;
  template class GenericIntegralOperator<HelmholtzDLKernel<3>>;
  template class GenericIntegralOperator<HelmholtzHSKernel<3>>;
  
  template class GenericIntegralOperator<CombinedFieldKernel<3>>;

  template class GenericIntegralOperator<MaxwellSLKernel<3>>;
  template class GenericIntegralOperator<MaxwellDLKernel<3>>;    
}
