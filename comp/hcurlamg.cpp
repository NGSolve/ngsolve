
#include <core/array.hpp>
#include "h1amg.hpp"
#include "hcurlamg.hpp"

#include <special_matrix.hpp>
#include <jacobi.hpp>

namespace ngcomp
{

  // some sparse matrix utilities -> should go to linalg
  
  shared_ptr<SparseMatrixTM<double>> ToSparseMatrix( const Embedding & e )
  {
    auto r = e.GetRange();
    auto h = e.Height();
    auto w = e.Width();

    Array<int> cnt(h);
    cnt = 0;
    cnt.Range(r) = 1;

    shared_ptr<SparseMatrixTM<double>> sp_mat;
    if(e.IsComplex())
        sp_mat = make_shared<SparseMatrix<double, Complex, Complex>>( cnt, w );
    else
        sp_mat = make_shared<SparseMatrix<double, double, double>>( cnt, w );

    auto & mat = *sp_mat;

    for(auto i : Range(r.Size()))
        mat(r[i], i) = 1.0;

    return sp_mat;
  }

  shared_ptr<SparseMatrixTM<double>> ToSparseMatrix( const EmbeddingTranspose & e )
  {
    auto r = e.GetRange();
    auto h = e.Height();
    auto w = e.Width();

    Array<int> cnt(h);
    cnt = 1;

    shared_ptr<SparseMatrixTM<double>> sp_mat;
    if(e.IsComplex())
        sp_mat = make_shared<SparseMatrix<double, Complex, Complex>>( cnt, w );
    else
        sp_mat = make_shared<SparseMatrix<double, double, double>>( cnt, w );

    auto & mat = *sp_mat;

    for(auto i : Range(r.Size()))
        mat(i, r[i]) = 1.0;

    return sp_mat;
  }


  template<typename TM, typename TROW, typename TCOL>
  shared_ptr<SparseMatrix<TM, TROW, TCOL>> MatAdd( const SparseMatrixTM<TM> & a_, const SparseMatrixTM<TM> & b_ )
    {
      auto & a = *dynamic_cast<const SparseMatrix<TM,TROW,TCOL>*>(&a_);
      auto & b = *dynamic_cast<const SparseMatrix<TM,TROW,TCOL>*>(&b_);

      auto w = a.Width();
      auto h = a.Height();

      NETGEN_CHECK_RANGE(b.Width(), w, w+1);
      NETGEN_CHECK_RANGE(b.Height(), h, h+1);

      auto firsta = a.GetFirstArray();
      auto firstb = b.GetFirstArray();
      auto cola = a.GetColIndices();
      auto colb = b.GetColIndices();

      Array<int> cnt(h);
      cnt = 0;

      for(auto row : Range(h))
        {
          auto cnta = firsta[row+1]-firsta[row];
          auto cntb = firstb[row+1]-firstb[row];

          if(cnta==0 || cntb==0)
            {
              cnt[row] = cnta+cntb;
              continue;
            }
          throw Exception("merge not implemented yet");
        }

      auto sp_mat = make_shared<SparseMatrix<TM, TROW, TCOL>>(cnt, w);
      auto & mat = *sp_mat;

      auto vala = a.GetValues();
      auto valb = b.GetValues();
      for(auto row : Range(h))
        {
          for(auto i : Range(firsta[row], firsta[row+1]))
              mat(row, cola[i]) = vala[i];

          for(auto i : Range(firstb[row], firstb[row+1]))
              mat(row, colb[i]) = valb[i];
        }

      return sp_mat;
    }

  template<typename MAT, typename SCAL>
  void CalcSchurComplement(const MAT a, FlatMatrix<SCAL> s,
                           const BitArray& used, LocalHeap& lh)
  {
    if (s.Height() == 0) return;
    if (s.Height() == a.Height())
      {
        s = a;
        return;
      }

    HeapReset hr(lh);

    int n = a.Height();
    Array<int> used_dofs(n, lh);
    Array<int> unused_dofs(n, lh);
    used_dofs.SetSize(0);
    unused_dofs.SetSize(0);
    for (int i = 0; i < n; i++)
      if (used[i])
        used_dofs.Append(i);
      else
        unused_dofs.Append(i);

    s = a.Rows(used_dofs).Cols(used_dofs);
    FlatMatrix<SCAL> b1 = a.Rows(unused_dofs).Cols(used_dofs) | lh;
    FlatMatrix<SCAL> b2 = a.Rows(used_dofs).Cols(unused_dofs) | lh;
    FlatMatrix<SCAL> c = a.Rows(unused_dofs).Cols(unused_dofs) | lh;
    FlatMatrix<SCAL> hb1 (b1.Height(), b1.Width(), lh);

    if (n > 10)
      {
        LapackInverse (c);
        hb1 = c * b1 | Lapack;
        s -= b2 * hb1 | Lapack;
      }
    else
      {
        CalcInverse (c);
        hb1 = c * b1;
        s -= b2 * hb1;
      }
  }





  /* *********************** HCurlAMG_Matrix  **************************** */


  
  template<typename SCAL>
  class HCurlAMG_Matrix : public BaseMatrix
  {
  protected:
    // parameters:
    HCurlAMG_Parameters param;
    bool node_on_each_level;
    bool need_vertex_prolongation = false;

    size_t size;
    int nv = 0;
    
    shared_ptr<SparseMatrixTM<SCAL>> mat;
    shared_ptr<BaseJacobiPrecond> smoother;
    shared_ptr<BaseBlockJacobiPrecond> blocksmoother;
    shared_ptr<SparseMatrixTM<double>> prolongation, restriction;
    shared_ptr<SparseMatrixTM<double>> vert_prolongation, vert_restriction;
    shared_ptr<SparseMatrixTM<double>> gradient, trans_gradient;
    shared_ptr<BaseMatrix> node_h1;
    shared_ptr<BaseMatrix> coarse_precond;
    
  public:
    HCurlAMG_Matrix (HCurlAMG_Parameters aparam,
                     bool _node_on_each_level = false
                     // bool _use_smoothed_prolongation = true
                     // int _coarsenings_per_level = 3,
                     // bool _blockjacobi_smoother = true
                     )
      
      : param(aparam),
        node_on_each_level(_node_on_each_level)
    {}

    HCurlAMG_Matrix(shared_ptr<SparseMatrixTM<SCAL>> _mat,
                    shared_ptr<BitArray> freedofs,
                    FlatArray<IVec<3>> f2e,
                    FlatArray<IVec<2>> e2v,
                    FlatArray<double> edge_weights,
                    FlatArray<double> face_weights,
                    size_t level,
                    HCurlAMG_Parameters aparam)
      : HCurlAMG_Matrix(aparam)
    {
      nv = 0;
      for (auto verts : e2v)
        nv = max3(nv, verts[0], verts[1]);
      nv++;
      Init(_mat, freedofs, f2e, e2v, edge_weights, face_weights, level);
    }

    void Init(shared_ptr<SparseMatrixTM<SCAL>> _mat,
              shared_ptr<BitArray> freedofs,
              FlatArray<IVec<3>> f2e,
              FlatArray<IVec<2>> e2v,
              FlatArray<double> edge_weights,
              FlatArray<double> face_weights,
              size_t level);

    
    int VHeight() const override { return size; }
    int VWidth() const override { return size; }
    bool IsComplex() const override { return is_same<SCAL, Complex>(); }
    AutoVector CreateRowVector() const override { return mat->CreateColVector(); }
    AutoVector CreateColVector() const override { return mat->CreateRowVector(); }

    void Mult(const BaseVector& f, BaseVector& u) const override;
    
  protected:

    Array<double> CalcEdgeCollapseWeights(FlatArray<IVec<3>> f2e,
                                          FlatArray<IVec<2>> e2v,
                                          FlatArray<double> edge_weights,
                                          FlatArray<double> face_weights,
                                          FlatTable<int> e2f) const;
    struct AMGInfo
    {
      Array<IVec<3>> f2e;
      Array<IVec<2>> e2v;
      Array<double> edge_weights;
      Array<double> face_weights;
      shared_ptr<BitArray> freedofs;
      shared_ptr<SparseMatrixTM<double>> prolongation;
      shared_ptr<SparseMatrixTM<double>> vert_prolongation;
    };
    virtual void BuildCoarseMat(const AMGInfo& cinfo, int level);

    virtual shared_ptr<BitArray>
    GetHCurlFreeDofs(shared_ptr<BitArray> freedofs) const
    {
      return freedofs;
    }

    virtual shared_ptr<BitArray>
    CreateCoarseFreedofs(shared_ptr<BitArray> freedofs,
                         int nce, int ncv,
                         FlatArray<int> edge_map,
                         FlatArray<int> vert_map,
                         FlatArray<IVec<2>> e2v,
                         FlatArray<IVec<2>> ce2v) const;

    virtual shared_ptr<BitArray> GetH1FreeDofs(FlatArray<IVec<2>> e2v,
                                               FlatTable<int> e2f,
                                               shared_ptr<BitArray> freedofs) const;

    AMGInfo CalcCoarsening(FlatArray<double> coll_weights,
                           shared_ptr<BitArray> freedofs,
                           FlatArray<IVec<3>> f2e,
                           FlatArray<IVec<2>> e2v,
                           FlatTable<int> e2f,
                           FlatArray<double> edge_weights,
                           FlatArray<double> face_weights,
                           int nv,
                           int level,
                           int coarsening) const;
  };






  
  template<typename SCAL>
  void HCurlAMG_Matrix<SCAL> ::
  Init (shared_ptr<SparseMatrixTM<SCAL>> _mat,
        shared_ptr<BitArray> freedofs,
        FlatArray<IVec<3>> f2e,
        FlatArray<IVec<2>> e2v,
        FlatArray<double> edge_weights,
        FlatArray<double> face_weights,
        size_t level)
  {
    static Timer timer("HCurlAMG_Matrix"); RegionTimer rt(timer);

    mat = _mat;
    size = mat->Height();
    
    auto ne = edge_weights.Size();
    auto nf = f2e.Size();

    if (param.verbose > 1)
      cout << IM(0) << "init level " << level << ", matsize = " << size << ", nedge = " << ne << ", nface = " << nf << endl;
    
    TableCreator<int> e2f_creator(ne);
    for(;!e2f_creator.Done(); e2f_creator++)
      for(auto f : Range(nf))
        for(int j = 0; j < 3; j++)
          e2f_creator.Add(f2e[f][j], f);
    Table<int> e2f = e2f_creator.MoveTable();

    int coarse_nv = nv;
    AMGInfo cinfo;
    Table<int> ce2f;

    for(auto coarsening : Range(param.coarsenings_per_level))
      {
        if(coarsening == 0)
          {
            auto coll_weights = CalcEdgeCollapseWeights(f2e, e2v, edge_weights,
                                                        face_weights, e2f);
            
            cinfo = CalcCoarsening(coll_weights, freedofs,
                                   f2e, e2v, e2f, edge_weights, face_weights,
                                   nv, level, coarsening);
            prolongation = cinfo.prolongation;
            if(need_vertex_prolongation)
              vert_prolongation = cinfo.vert_prolongation;
            auto nce = cinfo.e2v.Size();
            
            // build smoother
            if(param.block_smoother)
              {
                TableCreator<int> smoothing_blocks_creator(nce);
                for(; !smoothing_blocks_creator.Done(); smoothing_blocks_creator++)
                  for(auto i : Range(ne))
                    if(freedofs->Test(i))
                      for(auto j : prolongation->GetRowIndices(i))
                        smoothing_blocks_creator.Add(j, i);
                auto smoothing_blocks = make_shared<Table<int>>(smoothing_blocks_creator.MoveTable());
                blocksmoother = mat->CreateBlockJacobiPrecond(smoothing_blocks);
              }
            else
              {
                smoother = mat->CreateJacobiPrecond(GetHCurlFreeDofs(freedofs));
              }
          }
        else
          {
            auto coll_weights = CalcEdgeCollapseWeights(cinfo.f2e, cinfo.e2v,
                                                        cinfo.edge_weights,
                                                        cinfo.face_weights, ce2f);

            cinfo = std::move(CalcCoarsening(coll_weights, cinfo.freedofs,
                                             cinfo.f2e, cinfo.e2v, ce2f,
                                             cinfo.edge_weights, cinfo.face_weights,
                                             coarse_nv, level, coarsening));
            prolongation = MatMult(*prolongation, *cinfo.prolongation);
            if(need_vertex_prolongation)
              vert_prolongation = MatMult(*vert_prolongation,
                                          *cinfo.vert_prolongation);
          }

        if ( (cinfo.e2v.Size() < param.max_coarse) ||
             (level >= param.max_level) || 
             (cinfo.e2v.Size() == ne))
          break;
        if(coarsening < 2)
          {
            TableCreator<int> e2f_creator(ne);
            for(;!e2f_creator.Done(); e2f_creator++)
              for(auto f : Range(cinfo.f2e.Size()))
                for(int j = 0; j < 3; j++)
                  e2f_creator.Add(cinfo.f2e[f][j], f);
            ce2f = e2f_creator.MoveTable();
            coarse_nv = 0;
            for(const auto& verts : cinfo.e2v)
              coarse_nv = max3(coarse_nv, verts[0], verts[1]);
            coarse_nv++;
          }
      }

    // smoothed prolongation
    Array<int> nne(ne);
    if(param.use_smoothed_prolongation)
      {
        static Timer tsmprol("smoothed prolongation"); RegionTimer rsmprol(tsmprol);
        nne = 0;
        for(auto i : Range(ne))
          if(e2f[i].Size() > 0)
            nne[i] = 1 + 2 * e2f[i].Size();
        auto smoothprol = make_shared<SparseMatrix<double, SCAL, SCAL>>(nne, ne);
        double alpha = 0.5;
        
        // for(auto i : Range(ne))
        ParallelForRange(ne, [&](IntRange r)
        {
          Array<int> row_indices;
          for (auto i : r)
            {
              if(e2f[i].Size() == 0) continue; // unused dofs
              row_indices.SetSize0();
              row_indices.Append(i);
              for(auto face : e2f[i])
                for(auto edge : f2e[face])
                  if(edge != i)
                    row_indices.Append(edge);
              QuickSort(row_indices);

              for(auto j : row_indices)
                (*smoothprol)(i,j) = 0.;
              (*smoothprol)(i,i) += edge_weights[i];
              auto verts = e2v[i];
              for(auto face : e2f[i])
                {
                  for(auto e : f2e[face])
                    {
                      auto v2 = e2v[e];
                      auto orientation = verts[0] == v2[0] || verts[1] == v2[1] ? -1. : 1.;
                      if(e == i) orientation = 1.;
                      (*smoothprol)(i,e) += orientation * face_weights[face];
                    }
                }
              auto diag = (*smoothprol)(i,i);
              for(auto j : smoothprol->GetRowIndices(i))
                (*smoothprol)(i,j) *= -alpha/diag;
              (*smoothprol)(i,i) += 1.;
            }});
        prolongation = MatMult(*smoothprol, *prolongation);
        cout << IM(5) << "Smoothed prol nze: " << prolongation->NZE() << endl;
        cout << IM(5) << "Smoothed prol nze per row: " << double(prolongation->NZE())/prolongation->Height() << endl;
      }

    // Node correction only on first level seems enough in most cases
    if(node_on_each_level || level == 0)
      {
        nne = 2;
        for(auto e : Range(ne))
          if(e2f[e].Size() == 0)
            nne[e] = 0;
        gradient = make_shared<SparseMatrix<double, SCAL, SCAL>>(nne, nv);
        for(auto e : Range(ne))
          {
            if(e2f[e].Size() == 0) continue;
            (*gradient)(e, e2v[e][0]) = 1;
            (*gradient)(e, e2v[e][1]) = -1;
          }
      }

    BuildCoarseMat(cinfo, level);
    restriction = dynamic_pointer_cast<SparseMatrixTM<double>>
      (prolongation->CreateTranspose());

    if(node_on_each_level || level == 0)
      {
        trans_gradient = dynamic_pointer_cast<SparseMatrixTM<double>>(gradient->CreateTranspose());
        auto h1mat = mat->Restrict(*gradient);

        Array<double> v_weights(nv);
        v_weights = 0.;

        bool use_h1amg = true;
        // int count = 0;
        auto h1_freedofs = GetH1FreeDofs(e2v, e2f, freedofs);
        if(use_h1amg)
          node_h1 = make_shared<H1AMG_Matrix<SCAL>>
            (dynamic_pointer_cast<SparseMatrixTM<SCAL>>(h1mat), h1_freedofs, e2v,
             edge_weights, v_weights, 0);
        else
          {
            h1mat->SetInverseType(SPARSECHOLESKY);
            node_h1 = h1mat->InverseMatrix(h1_freedofs);
          }
      }
  }

  template<typename SCAL>
  void HCurlAMG_Matrix<SCAL> :: BuildCoarseMat(const AMGInfo& cinfo,
                                               int level)
  {
    auto coarsemat = mat->Restrict(*prolongation);
    cout << IM(5) << "mat nze: " << mat->NZE() << endl;
    cout << IM(5) << "mat nze per row: " << double(mat->NZE())/mat->Height() << endl;
    cout << IM(5) << "coarse mat nze: " << coarsemat->NZE() << endl;
    cout << IM(5) << "coarse mat nze per row: " << double(coarsemat->NZE())/coarsemat->Height() << endl;
    auto nce = cinfo.e2v.Size();
    auto ne = mat->Height();
    if ( (nce < param.max_coarse) || (level >= param.max_level) || (nce == ne))
      {
        if (param.verbose >= 2)
          cout << IM(0) << "coarse direct inverse, size = " << coarsemat->Height() << endl;
        coarsemat->SetInverseType(SPARSECHOLESKY);
        coarse_precond = coarsemat->InverseMatrix(cinfo.freedofs);
      }
    else
      {
        coarse_precond = make_shared<HCurlAMG_Matrix<SCAL>>
          (dynamic_pointer_cast<SparseMatrixTM<SCAL>>(coarsemat),
           cinfo.freedofs, cinfo.f2e, cinfo.e2v, cinfo.edge_weights,
           cinfo.face_weights, level+1, param);
      }
  }

  template<typename SCAL>
  shared_ptr<BitArray> HCurlAMG_Matrix<SCAL> ::
  GetH1FreeDofs(FlatArray<IVec<2>> e2v,
                FlatTable<int> e2f,
                shared_ptr<BitArray> freedofs) const
  {
    int nv = 0;
    for(const auto& verts : e2v)
      nv = max3(nv, verts[0], verts[1]);
    nv++;
    auto ne = e2v.Size();
    auto h1_freedofs = make_shared<BitArray>(nv);
    h1_freedofs->Set();
    for(auto e : Range(ne))
      if(e2f[e].Size() && !freedofs->Test(e))
        {
          h1_freedofs->Clear(e2v[e][0]);
          h1_freedofs->Clear(e2v[e][1]);
        }
    return h1_freedofs;
  }

  template<typename SCAL>
  Array<double> HCurlAMG_Matrix<SCAL> :: CalcEdgeCollapseWeights
  (FlatArray<IVec<3>> f2e, FlatArray<IVec<2>> e2v, FlatArray<double> edge_weights,
   FlatArray<double> face_weights, FlatTable<int> e2f) const
  {
    static Timer timer("HCurlAMG::CalcEdgeCollapseWeights"); RegionTimer rt(timer);
    auto ne = e2v.Size();

    Array<double> coll_eweights(ne);
    LocalHeap clh(20*1000*1000, "HCurlAMG::CalcEdgeCollapse");
    ParallelFor(ne, [&](size_t ei)
    {
      auto lh = clh.Split();
      // flux norm matrix
      auto ncf = e2f[ei].Size();
      if(ncf == 0)
        {
          coll_eweights[ei] = 0;
          return; // unused dof
        }
      
      FlatMatrix<double> top_part(ncf, lh);
      FlatMatrix<double> bot_part(ncf, lh);
      top_part = 0;
      bot_part = 0;
      auto verts = e2v[ei];
      const auto& faces = e2f[ei];
      for(auto i : Range(faces))
        {
          for(auto oe : f2e[faces[i]])
            {
              if(oe == ei) continue;
              auto overts = e2v[oe];
              if(overts[0] == verts[0] || overts[1] == verts[0])
                top_part(i, i) += edge_weights[oe];
              if(overts[0] == verts[1] || overts[1] == verts[1])
                bot_part(i, i) += edge_weights[oe];
            }
        }
      CalcInverse(top_part);
      CalcInverse(bot_part);
      FlatMatrix<double> sum = top_part + bot_part | lh;
      CalcInverse(sum);
      // make m-tilde
      FlatVector<double> ones(ncf, lh);
      ones = 1.;
      FlatVector<double> col_sum = Trans(sum) * ones | lh;
      double denominator = 0.;
      for(auto v : col_sum)
        denominator += v;
      denominator += edge_weights[ei];
      for(auto i : Range(ncf))
        for(auto j : Range(ncf))
          sum(i,j) -= col_sum[i] * col_sum[j] / denominator;
      for(auto i : Range(faces))
        sum(i,i) += face_weights[faces[i]];
      if(sum.Height() != 0)
        {
          top_part = 0.;
          bot_part = 0.;
          for(auto i : Range(faces))
            {
              auto face = faces[i];
              for(auto oe : f2e[face])
                {
                  if(oe == ei) continue;
                  auto overts = e2v[oe];
                  if(verts[0] == overts[0] || verts[0] == overts[1])
                    {
                      for(auto oface : e2f[oe])
                        {
                          if(oface == face) continue;
                          top_part(i,i) += face_weights[oface];
                        }
                      top_part(i,i) += edge_weights[oe];
                    }
                  if(verts[1] == overts[0] || verts[1] == overts[1])
                    {
                      for(auto oface : e2f[oe])
                        {
                          if(oface == face) continue;
                          bot_part(i,i) += face_weights[oface];
                        }
                      bot_part(i,i) += edge_weights[oe];
                    }
                }
            }
          // auto face = faces[0];
          // find center face neighbours - neighbours have 2 common neighbour faces
          Array<IVec<2>> neighbours(faces.Size());
          for(auto i : Range(neighbours))
            neighbours[i] = IVec<2>(-1,-1);
          for(auto fi : Range(faces))
            {
              if(neighbours[fi][0] != -1 && neighbours[fi][1] != -1) continue;
              int te1 = 0;
              for(auto e : f2e[faces[fi]])
                {
                  if(ei == e) continue;
                  if(verts[0] == e2v[e][0] || verts[0] == e2v[e][1])
                    te1 = e;
                }
              for(auto fo : Range(faces))
                {
                  if(neighbours[fo][0] != -1 && neighbours[fo][1] != -1) continue;
                  if(neighbours[fi][0] == fo || fi == fo) continue;
                  int te2 = 0;
                  for(auto e : f2e[faces[fo]])
                    {
                      if(ei == e) continue;
                      if(verts[0] == e2v[e][0] || verts[0] == e2v[e][1])
                        te2 = e;
                    }
                  for(auto f : e2f[te1])
                    {
                      if(f == faces[fi] || f == faces[fo]) continue;
                      auto edges = f2e[f];
                      if(edges[0] == te2 || edges[1] == te2 || edges[2] == te2)
                        {
                          if(neighbours[fi][0] == -1)
                            neighbours[fi][0] = fo;
                          else
                            neighbours[fi][1] = fo;
                          if(neighbours[fo][0] == -1)
                            neighbours[fo][0] = fi;
                          else
                            neighbours[fo][1] = fi;
                        }
                    }
                }
            }
          for(auto i : Range(neighbours))
            neighbours[i].Sort();
          for(auto i : Range(neighbours))
            {
              for(auto j : neighbours[i])
                {
                  if(j == -1) continue;
                  if(i < j)
                    {
                      int top_edge1=-1, bot_edge1=-1, top_edge2=-1, bot_edge2=-1;
                      for(auto e1 : f2e[faces[i]])
                        {
                          if(e1 == ei) continue;
                          if(e2v[e1][0] == verts[0] || e2v[e1][1] == verts[0])
                            top_edge1 = e1;
                          if(e2v[e1][0] == verts[1] || e2v[e1][1] == verts[1])
                            bot_edge1 = e1;
                        }
                      for(auto e2 : f2e[faces[j]])
                        {
                          if(e2 == ei) continue;
                          if(e2v[e2][0] == verts[0] || e2v[e2][1] == verts[0])
                            top_edge2 = e2;
                          if(e2v[e2][0] == verts[1] || e2v[e2][1] == verts[1])
                            bot_edge2 = e2;
                        }
                      for(auto f : e2f[top_edge1])
                        {
                          auto edges = f2e[f];
                          if(edges[0] == top_edge2 || edges[1] == top_edge2 ||
                             edges[2] == top_edge2)
                            {
                              top_part(i,j) -= face_weights[f];
                              top_part(j,i) -= face_weights[f];
                            }
                        }
                      for(auto f : e2f[bot_edge1])
                        {
                          auto edges = f2e[f];
                          if(edges[0] == bot_edge2 || edges[1] == bot_edge2 ||
                             edges[2] == bot_edge2)
                            {
                              bot_part(i,j) -= face_weights[f];
                              bot_part(j,i) -= face_weights[f];
                            }
                        }
                    }
                }
            }
          double eps = 1e-10 * L2Norm(top_part);
          for(auto &d : top_part.Diag())
            d += eps;
          eps = 1e-10 * L2Norm(bot_part);
          for(auto &d : bot_part.Diag())
            d += eps;
          CalcInverse(top_part);
          CalcInverse(bot_part);
          FlatMatrix<double> cor_sum = top_part + bot_part | lh;
          CalcInverse(cor_sum);
          FlatVector<double> col_sum = Trans(cor_sum) * ones | lh;
          double denominator = 0.;
          for(auto v : col_sum)
            denominator += v;
          denominator += edge_weights[ei];
          for(auto i : Range(ncf))
            for(auto j : Range(ncf))
              cor_sum(i,j) -= col_sum[i] * col_sum[j] / denominator;
          for(auto i : Range(faces))
            cor_sum(i,i) += face_weights[faces[i]];
          Vector<double> lami(ncf);
          LapackEigenValuesSymmetric(sum, cor_sum, lami);
          coll_eweights[ei] = lami[0] + double(e2f[ei].Size())/50.;
        }
    }, TasksPerThread(10));
    return coll_eweights;
  }

  template<typename SCAL>
  typename HCurlAMG_Matrix<SCAL>::AMGInfo HCurlAMG_Matrix<SCAL> ::
  CalcCoarsening(FlatArray<double> coll_weights,
                 shared_ptr<BitArray> freedofs,
                 FlatArray<IVec<3>> f2e,
                 FlatArray<IVec<2>> e2v,
                 FlatTable<int> e2f,
                 FlatArray<double> edge_weights,
                 FlatArray<double> face_weights,
                 int nv, int level, int coarsening) const
  {
    AMGInfo info;
    auto ne = edge_weights.Size();
    auto nf = f2e.Size();
    Array<int> indices(ne);
    for(auto i : Range(ne))
      indices[i] = i;
    SampleSortI(coll_weights, indices);

    Array<bool> collapsed_edge(ne);
    Array<bool> collapsed_vertex(nv);
    collapsed_edge = false;
    collapsed_vertex = false;
    Array<bool> unused_vertex(nv);
    unused_vertex = true;

    for(int i = ne-1; i >= 0; i--)
      {
        if(e2f[i].Size() == 0) continue; // unused dof
        auto e = indices[i];
        auto verts = e2v[e];
        unused_vertex[verts[0]] = false;
        unused_vertex[verts[1]] = false;
        if(coll_weights[e] >= 0.1 && !collapsed_vertex[verts[0]] &&
           !collapsed_vertex[verts[1]])
          {
            collapsed_edge[e] = true;
            collapsed_vertex[verts[0]] = true;
            collapsed_vertex[verts[1]] = true;
          }
      }
    collapsed_vertex = false;
    for(auto e : Range(ne))
      if(collapsed_edge[e])
        {
          collapsed_vertex[max2(e2v[e][0], e2v[e][1])] = true;
        }
    size_t ncv = 0; // number of coarse vertices
    Array<int> vert_map(nv);
    for(auto i : Range(nv))
      if(!unused_vertex[i] && !collapsed_vertex[i])
        vert_map[i] = ncv++;
    for(auto e : Range(ne))
      if(collapsed_edge[e])
        vert_map[max2(e2v[e][0], e2v[e][1])] = vert_map[min2(e2v[e][0], e2v[e][1])];
    for(auto v : Range(nv))
      if(unused_vertex[v])
        vert_map[v] = -1;

    size_t nce = 0; // number of coarse edges
    Array<int> edge_map(ne);
    info.e2v.SetAllocSize(ne);
    edge_map = -1;
    HashTable<IVec<2>, int> cv2e(ne);
    Array<bool> inverted_edge(ne);
    inverted_edge = false;
    for(auto e : Range(ne))
      {
        if(e2f[e].Size() == 0) continue; // unused dof
        auto v = e2v[e];
        auto cv = IVec<2>(vert_map[v[0]], vert_map[v[1]]);
        if(cv[0] == cv[1]) continue;
        if(cv[0] > cv[1])
          inverted_edge[e] = true;
        cv.Sort();
        if(cv2e.Used(cv))
          {
            edge_map[e] = cv2e[cv];
          }
        else
          {
            info.e2v.Append(cv);
            cv2e[cv] = nce;
            edge_map[e] = nce++;
          }
      }

    size_t ncf = 0; // number of coarse faces
    Array<int> face_map(nf);
    face_map = -1;
    info.f2e.SetAllocSize(nf);
    HashTable<IVec<3>, int> ce2f(nf);
    for(auto i : Range(nf))
      {
        bool collapsed = collapsed_edge[f2e[i][0]] || collapsed_edge[f2e[i][1]] ||
          collapsed_edge[f2e[i][2]];
        if(!collapsed)
          {
            auto f = IVec<3>(edge_map[f2e[i][0]], edge_map[f2e[i][1]],
                            edge_map[f2e[i][2]]).Sort();
            if(ce2f.Used(f))
              {
                face_map[i] = ce2f[f];
              }
            else
              {
                ce2f[f] = ncf;
                info.f2e.Append(f);
                face_map[i] = ncf++;
              }
          }
      }
    info.edge_weights.SetSize(nce);
    info.edge_weights = 0.;
    for(auto edge : Range(ne))
      if(edge_map[edge] != -1)
        info.edge_weights[edge_map[edge]] += edge_weights[edge];

    info.face_weights.SetSize(ncf);
    info.face_weights = 0.;
    for(auto face : Range(nf))
      if(face_map[face] != -1)
        info.face_weights[face_map[face]] += face_weights[face];

    // build prolongation
    Array<int> nne(ne);
    for(auto i : Range(ne))
      nne[i] = (edge_map[i] != -1) ? 1 : 0;
    info.prolongation = make_shared<SparseMatrix<double, SCAL, SCAL>>(nne, nce);
    for(auto i : Range(ne))
      if(edge_map[i] != -1)
        (*info.prolongation)(i, edge_map[i]) = inverted_edge[i] ? -1. : 1;

    if(need_vertex_prolongation)
      {
        Array<int> nnv(nv);
        nnv = 1;
        for(auto v : Range(nv))
          if(unused_vertex[v])
            nnv[v] = 0;
        info.vert_prolongation = make_shared<SparseMatrix<double, SCAL, SCAL>>(nnv, ncv);
        for(auto i : Range(nv))
          if(!unused_vertex[i])
            (*info.vert_prolongation)(i, vert_map[i]) = 1.;
      }

    info.freedofs = CreateCoarseFreedofs(freedofs, nce, ncv, edge_map, vert_map,
                                         e2v, info.e2v);

    cout << IM(3) << "HCurl Level: " << level << " coarsening " << coarsening << endl;
    cout << IM(5) << "Nr. vertices: " << nv << " Nr. edges: " << ne << endl;
    cout << IM(5) << "Nr. cvertices: " << ncv << " Nr. cedges: " << nce << endl;
    cout << IM(5) << "e/v: " << double(nce)/ncv << endl;
    cout << IM(5) << "coarse/fine edges: " << double(nce)/ne
         << " coarse/fine verts: " << double(ncv)/nv << endl;
    return info;
  }

  template<typename SCAL> shared_ptr<BitArray>
  HCurlAMG_Matrix<SCAL> :: CreateCoarseFreedofs(shared_ptr<BitArray> freedofs,
                                                int nce, int ncv,
                                                FlatArray<int> edge_map,
                                                FlatArray<int> vert_map,
                                                FlatArray<IVec<2>> e2v,
                                                FlatArray<IVec<2>> ce2v) const
  {
    auto cfd = make_shared<BitArray>(nce);
    BitArray cfreeverts(ncv);
    cfreeverts.Set();
    for(auto e : Range(e2v.Size()))
      if(!freedofs->Test(e))
        for(auto v : e2v[e])
          cfreeverts.Clear(vert_map[v]);
    cfd->Set();
    for(auto ce : Range(nce))
      {
        auto verts = ce2v[ce];
        if(!cfreeverts[verts[0]] && !cfreeverts[verts[1]])
          cfd->Clear(ce);
      }
    return cfd;
  }

  template<typename SCAL>
  void HCurlAMG_Matrix<SCAL> :: Mult(const BaseVector& f, BaseVector& u) const
  {
    static Timer timer("HCurlAMG::Mult");
    RegionTimer regt(timer);
    static Timer timer_n("Node correction");
    static Timer timer_c("Coarse correction");

    u = 0.;
    if(smoother)
      for (int k = 0; k < param.smoothing_steps; k++)      
        smoother->GSSmooth(u, f);
    else
      blocksmoother->GSSmooth(u, f, param.smoothing_steps);
    
    auto residuum = f.CreateVector();

    if(gradient)
      {
        residuum = f - (*mat) * u;
        auto op = gradient*node_h1*trans_gradient;
        u += *op * residuum;
      }
    
    {
      RegionTimer rt(timer_c);
      residuum = f - (*mat) * u;

      auto op = prolongation * coarse_precond * restriction;
      u += *op * residuum;      
    }

    if(gradient)
      {
        residuum = f - (*mat) * u;
        auto op = gradient*node_h1*trans_gradient;
        u += *op * residuum;
      }
    
    if(smoother)
      for (int k = 0; k < param.smoothing_steps; k++)
        smoother->GSSmoothBack(u, f);
    else
      blocksmoother->GSSmoothBack(u, f, param.smoothing_steps);
  }


  /* ************************* HCurlAMG preconditioner ************************** */

  
  HCurlAMG :: HCurlAMG(shared_ptr<BilinearForm> bfa, const Flags& flags,
                       const string& name)
    : Preconditioner(bfa, flags, name)
  {
    fes = bfa->GetFESpace();
    node_on_each_level = true; // flags.GetDefineFlag("dirichlet_on_each_level");

    param.verbose = int(flags.GetNumFlag("verbose", 0));
    param.smoothing_steps = int(flags.GetNumFlag("smoothingsteps", 1));    
    param.use_smoothed_prolongation = flags.GetDefineFlagX("smoothedprolongation").IsMaybeTrue();
    param.max_coarse = int(flags.GetNumFlag("maxcoarse", 10));
    param.max_level = int(flags.GetNumFlag("maxlevel", 20));
  }


  DocInfo HCurlAMG :: GetDocu()
  {
    DocInfo docu;
    docu.short_docu = "AMG-predconditioner for edge elements";

    docu.long_docu =
      R"raw_string(using agglomeration-based AMG
)raw_string";      
    
    docu.Arg("smoothingsteps") = "int = 3\n"
      "  number of pre and post-smoothing steps";
    docu.Arg("smoothedprolongation") = "bool = true\n"
      "  use smoothed prolongation";
    docu.Arg("maxcoarse") = "int = 10\n"
      "  maximal dofs on level to switch to direct solver\n";
    docu.Arg("maxlevel") = "int = 20\n"
      "  maximal refinement levels to switch to direct solver\n";
    docu.Arg("verbose") = "int = 3\n"
      "  verbosity level, 0..no output, 5..most output";

    return docu;
  }

  
  template<typename SCAL>
  void HCurlAMG :: AddElementMatrixCommon(FlatArray<int> dnums,
                                          FlatMatrix<SCAL> elmat,
                                          ElementId id, LocalHeap & lh)
  {
    if(L2Norm2(elmat) == 0)
      return;

    auto ndof = dnums.Size();
    auto & ma = fes->GetMeshAccess();
    HeapReset hr(lh);
    BitArray used(ndof, lh);

    FlatMatrix<SCAL> schur_edge(1, lh);
    for(auto i : Range(ndof))
      {
        schur_edge = 0;
        used.Clear();
        used.SetBit(i);
        CalcSchurComplement(elmat, schur_edge, used, lh);
        double weight = fabs(schur_edge(0,0));
        edge_weights_ht.Do(IVec<1>(dnums[i]), [weight] (auto &v) { v+= weight; });
      }

    auto el = ma->GetElement(id);
    ArrayMem<DofId, 3> fdofs;
    ArrayMem<DofId, 3> local_fdofs;
    FlatMatrix<SCAL> schur_face(3, lh);
    for(const auto& fnr : el.Faces())
      {
        auto fdofs = ma->GetFaceEdges(fnr);
        local_fdofs.SetSize(0);
        for(auto d : fdofs)
          local_fdofs.Append(dnums.Pos(d));
        schur_face = 0.;
        used.Clear();
        for(auto ld : local_fdofs)
          used.SetBit(ld);
        CalcSchurComplement(elmat, schur_face, used, lh);
        double face_weight = 0.;
        for(auto i : Range(3))
          face_weight += fabs(schur_face(i,i));
        face_weight /= 3;
        if(!isnan(face_weight))
          face_weights_ht.Do(IVec<3>(fdofs[0], fdofs[1], fdofs[2]).Sort(),
                             [face_weight] (auto& fw) { fw += face_weight; });
      }
  }

  void HCurlAMG :: FinalizeLevel(const BaseMatrix* matrix)
  {
    static Timer timer("HCurlAMG::FinalizeLevel"); RegionTimer rt(timer);

    if (param.verbose >= 1)
      cout << IM(0) << "Create AMG matrix" << endl;
    if (param.verbose >= 2)
      {
        cout << IM(0) << "smoothingsteps = " << param.smoothing_steps << endl;
        cout << IM(0) << "smoothedprolongation = " << int(param.use_smoothed_prolongation) << endl;
        cout << IM(0) << "maxcoarse = " << param.max_coarse << endl;                
        cout << IM(0) << "maxlevel = " << param.max_level << endl;                
      }
    
    auto num_edges = matrix->Height();
    auto num_faces = face_weights_ht.Used();

    Array<double> edge_weights(num_edges);
    Array<double> face_weights(num_faces);
    Array<IVec<3>> f2e(num_faces);

    ParallelForRange (edge_weights.Size(), [&](IntRange r) {
      edge_weights.Range(r) = 0.0;
    });
    
    edge_weights_ht.IterateParallel
      ([&edge_weights](size_t i, IVec<1> key, double weight)
      {
        edge_weights[key[0]] = weight;
      });
    edge_weights_ht = ParallelHashTable<IVec<1>, double>();

    face_weights_ht.IterateParallel
      ([&face_weights, &f2e](size_t i, IVec<3> key, double weight)
      {
        face_weights[i] = weight;
        f2e[i] = key;
      });
    face_weights_ht = ParallelHashTable<IVec<3>, double>();

    Array<IVec<2>> e2v(num_edges);
    for(auto edge : ma->Edges())
      e2v[edge.GetNr()] = ma->GetEdgePNums(edge.GetNr()).Sort();

    if(matrix->IsComplex())
      {
        auto amat = dynamic_pointer_cast<SparseMatrixTM<Complex>>
          (const_cast<BaseMatrix*>(matrix)->shared_from_this());
        mat = make_shared<HCurlAMG_Matrix<Complex>>(amat, freedofs, f2e, e2v,
                                                    edge_weights, face_weights, 0, param);
      }
    else
      {
        auto amat = dynamic_pointer_cast<SparseMatrixTM<double>>
          (const_cast<BaseMatrix*>(matrix)->shared_from_this());
        mat = make_shared<HCurlAMG_Matrix<double>>(amat, freedofs, f2e, e2v,
                                                   edge_weights, face_weights, 0, param);
      }
  }






  /* ***************** A-Phi AMG ***************************** */


  template<typename SCAL>
  class APhiMatrix : public HCurlAMG_Matrix<SCAL>
  {
    using BASE = HCurlAMG_Matrix<SCAL>;
    using BASE::param;
    size_t hcurlsize;
    size_t h1size;
  public:
    APhiMatrix(shared_ptr<SparseMatrixTM<SCAL>> _mat,
               shared_ptr<BitArray> freedofs,
               FlatArray<IVec<3>> f2e,
               FlatArray<IVec<2>> e2v,
               FlatArray<double> edge_weights,
               FlatArray<double> face_weights,
               size_t level, HCurlAMG_Parameters aparam);
  protected:
    void BuildCoarseMat(const typename BASE::AMGInfo& cinfo, int level) override;
    shared_ptr<BitArray>
    GetHCurlFreeDofs(shared_ptr<BitArray> freedofs) const override
    {
      auto fd = make_shared<BitArray>(freedofs->Size());
      fd->Clear();
      for(auto i : Range(hcurlsize))
        if(freedofs->Test(i))
          fd->SetBit(i);
      return fd;
    }
    shared_ptr<BitArray> GetH1FreeDofs(FlatArray<IVec<2>> e2v,
                                       FlatTable<int> e2f,
                                       shared_ptr<BitArray> freedofs) const override
    {
      auto fd = make_shared<BitArray>(h1size);
      fd->Clear();
      for(auto i : Range(h1size))
        if(freedofs->Test(hcurlsize + i))
          fd->SetBit(i);
      return fd;
    }
  shared_ptr<BitArray> CreateCoarseFreedofs(shared_ptr<BitArray> freedofs,
                                            int nce, int ncv,
                                            FlatArray<int> edge_map,
                                            FlatArray<int> vert_map,
                                            FlatArray<IVec<2>> e2v,
                                            FlatArray<IVec<2>> ce2v) const override;
  };
  


  

  APhiHCurlAMG :: APhiHCurlAMG(shared_ptr<BilinearForm> _bfa,
                               const Flags& flags,
                               const string& name)
    : HCurlAMG(_bfa, flags, name), bfa(_bfa)
  {
  }

  void APhiHCurlAMG::AddElementMatrix(FlatArray<int> dnums,
                                      FlatMatrix<double> elmat,
                                      ElementId id, LocalHeap & lh)
  {
    AddElementMatrixCommon(dnums, elmat, id, lh);
  }

  void APhiHCurlAMG::AddElementMatrix(FlatArray<int> dnums,
                                      FlatMatrix<Complex> elmat,
                                      ElementId id, LocalHeap & lh)
  {
    AddElementMatrixCommon(dnums, elmat, id, lh);
  }

  void APhiHCurlAMG :: FinalizeLevel(const BaseMatrix* matrix)
  {
    static Timer timer("APhiHCurlAMG::FinalizeLevel");
    RegionTimer rt(timer);

    auto num_edges = this->ma->GetNEdges();
    auto num_faces = face_weights_ht.Used();

    Array<double> edge_weights(num_edges);
    Array<double> face_weights(num_faces);
    Array<IVec<3>> f2e(num_faces);

    edge_weights_ht.IterateParallel
      ([&edge_weights](size_t i, IVec<1> key, double weight)
      {
        edge_weights[key[0]] = weight;
      });
    edge_weights_ht = ParallelHashTable<IVec<1>, double>();

    face_weights_ht.IterateParallel
      ([&face_weights, &f2e](size_t i, IVec<3> key, double weight)
      {
        face_weights[i] = weight;
        f2e[i] = key;
      });
    face_weights_ht = ParallelHashTable<IVec<3>, double>();

    Array<IVec<2>> e2v(num_edges);
    // for(auto edge : ma->Edges())
    //   e2v[edge.GetNr()] = ma->GetEdgePNums(edge.GetNr()).Sort();
    for(auto i : Range(num_edges))
      e2v[i] = ma->GetEdgePNums(i).Sort();

    int nv = 0;
    for(auto verts : e2v)
      nv = max3(nv, verts[0], verts[1]);

    param.block_smoother = false; // blocks are singular
    param.use_smoothed_prolongation = false;
    
    if(matrix->IsComplex())
      {
        // auto amat = dynamic_pointer_cast<SparseMatrixTM<Complex>>
        //   (const_cast<BaseMatrix*>(matrix)->shared_from_this());
        auto amat = dynamic_pointer_cast<SparseMatrixTM<Complex>>
          (bfa->GetMatrixPtr());
        mat = make_shared<APhiMatrix<Complex>>(amat, freedofs, f2e, e2v,
                                               edge_weights, face_weights, 0, param);
      }
    else
      {
        // auto amat = dynamic_pointer_cast<SparseMatrixTM<double>>
        //   (const_cast<BaseMatrix*>(matrix)->shared_from_this());
        auto amat = dynamic_pointer_cast<SparseMatrixTM<double>>
          (bfa->GetMatrixPtr());
        mat = make_shared<APhiMatrix<double>>(amat, freedofs, f2e, e2v,
                                              edge_weights, face_weights, 0, param);
      }
  }

  template<typename SCAL>
  APhiMatrix<SCAL> :: APhiMatrix(shared_ptr<SparseMatrixTM<SCAL>> _mat,
                                 shared_ptr<BitArray> freedofs,
                                 FlatArray<IVec<3>> f2e,
                                 FlatArray<IVec<2>> e2v,
                                 FlatArray<double> edge_weights,
                                 FlatArray<double> face_weights,
                                 size_t level, HCurlAMG_Parameters param)
    : HCurlAMG_Matrix<SCAL>(param, true /* false, 3, false */)
  {
    this->need_vertex_prolongation = true;
    hcurlsize = e2v.Size();
    h1size = freedofs->Size() - hcurlsize;
    this->nv = h1size;
    this->Init(_mat, freedofs, f2e, e2v, edge_weights, face_weights, level);
  }

  template<typename SCAL>
  void APhiMatrix<SCAL> :: BuildCoarseMat(const typename BASE::AMGInfo& cinfo,
                                          int level)
  {
    auto chcurlsize = this->prolongation->Width();
    auto ch1size = this->vert_prolongation->Width();
    auto csize = chcurlsize + ch1size;

    constexpr bool is_complex = is_same_v<SCAL, Complex>;
    auto e1 = ToSparseMatrix(Embedding(this->size, IntRange(0, hcurlsize), is_complex));
    auto r1 = ToSparseMatrix(EmbeddingTranspose(csize, IntRange(0, chcurlsize), is_complex));
    auto first_prol = MatMult(*this->prolongation, *r1);
    first_prol = MatMult(*e1, *first_prol);

    auto e2 = ToSparseMatrix(Embedding(this->size, IntRange(hcurlsize, this->size), is_complex));
    auto r2 = ToSparseMatrix(EmbeddingTranspose(csize, IntRange(chcurlsize, csize), is_complex));
    auto second_prol = MatMult(*this->vert_prolongation, *r2);
    second_prol = MatMult(*e2, *second_prol);

    this->prolongation = MatAdd<double, SCAL, SCAL>(*first_prol, *second_prol);
    
    this->gradient = e2;

    auto coarsemat = this->mat->Restrict(*this->prolongation);
    cout << IM(5) << "mat nze: " << this->mat->NZE() << endl;
    cout << IM(5) << "mat nze per row: " << double(this->mat->NZE())/this->mat->Height() << endl;
    cout << IM(5) << "coarse mat nze: " << coarsemat->NZE() << endl;
    cout << IM(5) << "coarse mat nze per row: " << double(coarsemat->NZE())/coarsemat->Height() << endl;
    auto nce = chcurlsize;
    auto ne = hcurlsize;
    if((nce < param.max_coarse) || (nce == ne))
      {
        for(auto i : Range(chcurlsize, csize))
          cinfo.freedofs->Clear(i);
        coarsemat->SetInverseType(SPARSECHOLESKY);
        this->coarse_precond = coarsemat->InverseMatrix(cinfo.freedofs);
      }
    else
      {
        this->coarse_precond = make_shared<APhiMatrix<SCAL>>
          (dynamic_pointer_cast<SparseMatrixTM<SCAL>>(coarsemat),
           cinfo.freedofs, cinfo.f2e, cinfo.e2v, cinfo.edge_weights,
           cinfo.face_weights, level+1, this->param);
      }
  }

  template<typename SCAL> shared_ptr<BitArray>
  APhiMatrix<SCAL> :: CreateCoarseFreedofs(shared_ptr<BitArray> freedofs,
                                           int nce, int ncv,
                                           FlatArray<int> edge_map,
                                           FlatArray<int> vert_map,
                                           FlatArray<IVec<2>> e2v,
                                           FlatArray<IVec<2>> ce2v) const
  {
    auto cfd = make_shared<BitArray>(nce + ncv);
    BitArray cfreeverts(ncv);
    cfreeverts.Set();
    for(auto e : Range(e2v.Size()))
      if(!freedofs->Test(e))
        for(auto v : e2v[e])
          cfreeverts.Clear(vert_map[v]);
    cfd->Set();
    for(auto ce : Range(nce))
      {
        auto verts = ce2v[ce];
        if(!cfreeverts[verts[0]] && !cfreeverts[verts[1]])
          cfd->Clear(ce);
      }
    for(auto v : Range(vert_map.Size()))
      if(vert_map[v] != -1 && !freedofs->Test(edge_map.Size() + v))
        cfd->Clear(nce + vert_map[v]);
    return cfd;
  }

  template<typename SCAL>
  void APhiHCurlAMG::AddElementMatrixCommon(FlatArray<int> dnums,
                                            FlatMatrix<SCAL> belmat,
                                            ElementId id, LocalHeap & lh)
  {
    if(dnums.Size() != 10)
      throw Exception("Only tets implemented yet!");
    if(L2Norm(belmat) == 0)
      return;
    auto ndofhc = 6;
    auto ndofh1 = 4;
    auto dnumshc = dnums.Range(0,ndofhc);
    // auto dnumsh1 = dnums.Range(ndofhc, END);
    FlatMatrix<SCAL> elmathc = belmat.Rows(ndofhc).Cols(ndofhc) | lh;
    double reg = 1e-10 * L2Norm(elmathc);
    for(auto i : Range(ndofhc))
      elmathc(i,i) += reg;
    auto elmath1 = belmat.Rows(ndofhc, ndofhc+ndofh1).Cols(ndofhc, ndofhc+ndofh1);
    auto ma = fes->GetMeshAccess();
    HeapReset hr(lh);
    BitArray usedh1(ndofh1, lh);

    FlatMatrix<SCAL> schur_edge(2, lh);
    auto e2v = ElementTopology::GetEdges(ET_TET);
    for(auto edge : Range(ndofhc))
      {
        auto v1v2 = e2v[edge];
        schur_edge = 0;
        usedh1.Clear();
        usedh1.SetBit(v1v2[0]);
        usedh1.SetBit(v1v2[1]);
        CalcSchurComplement(elmath1, schur_edge, usedh1, lh);
        double weight = fabs(schur_edge(0,0)) + fabs(schur_edge(1,1));
        weight /= 2.;
        edge_weights_ht.Do(IVec<1>(dnumshc[edge]), [weight] (auto &v) { v+= weight; });
      }
    FlatMatrix<SCAL> schur_face(3, lh);
    BitArray usedhc(ndofhc, lh);
    static int f2e[4][3] =
      { { 1, 2, 5},
        { 0, 2, 4},
        { 0, 1, 3},
        { 3, 4, 5}
      };
    for(auto face : Range(4))
      {
        schur_face = 0.;
        usedhc.Clear();
        for(auto e : f2e[face])
          usedhc.SetBit(e);
        CalcSchurComplement(elmathc, schur_face, usedhc, lh);
        double face_weight = 0.;
        for(auto i : Range(3))
          face_weight += fabs(schur_face(i,i));
        face_weight /= 3;
        if(!isnan(face_weight))
          face_weights_ht.Do(IVec<3>(dnumshc[f2e[face][0]],
                                    dnumshc[f2e[face][1]],
                                    dnumshc[f2e[face][2]]).Sort(),
                             [face_weight] (auto& fw) { fw += face_weight; });
      }
  }

  static RegisterPreconditioner<HCurlAMG> inithcamg ("hcurlamg");
  static RegisterPreconditioner<APhiHCurlAMG> initaphiamg ("APhiamg");

  
} // namespace hcurlamg
