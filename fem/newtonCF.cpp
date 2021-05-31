/*********************************************************************/
/* File:   newtonCF.cpp                                              */
/* Authors: Joachim Schoeberl, Matthias Rambausek                    */
/* Date:   Feb 2021                                                  */
/*********************************************************************/

#include <cmath>
#include <fem.hpp>
#include <ngstd.hpp>

namespace ngfem {

namespace {
template <typename T> auto LInfNorm(const T &array) -> auto {
  using number_t = remove_reference_t<remove_cv_t<decltype(*(array.Data()))>>;
  number_t s = 0;
  for (const number_t item : array) {
    if (isnan(item))
      return numeric_limits<number_t>::quiet_NaN();
    s = abs(item) > s ? abs(item) : s;
  }
  return s;
}

bool has_vs_embedding(const ProxyFunction *proxy) {
  return proxy->Evaluator()->GetVSEmbedding().has_value();
}

auto get_vs_embedding(const ProxyFunction *proxy) {
  return proxy->Evaluator()->GetVSEmbedding();
}

int proxy_dof_dimension(const ProxyFunction *proxy) {
  if (const auto vsemb = get_vs_embedding(proxy))
    return vsemb->Width();
  else
    return proxy->Dimension();
}

} // namespace

// TODO: (general)
//  * Interpolation into generic compound spaces does not work as expected
//     (only one component is respected)
//  * How to handle consistent linearizations? This is probably a bigger topic
//     because there is also no support for nonlinear equations at the FEM level
//     (only nonlinear energies...). One could mix the linearization from an
//     incremental variational principle with the evolution of via a general
//     nonlinear equation. The effect on the convergence might not be huge in
//     most cases. If it is, a smaller time step might be a good idea anyway.

class NewtonCF : public CoefficientFunction {
  shared_ptr<CoefficientFunction> expression;
  Array<shared_ptr<CoefficientFunction>> startingpoints{};

  Array<ProxyFunction *> proxies{};
  Array<CoefficientFunction *> cachecf{};

  // The total dimension of the linear system (because of VS embeddings, this
  // can be different from this->Dimension())
  int numeric_dim = 0;
  int full_dim = 0;
  Array<int> proxy_dims{};

  // Same parameters as for scipy's newton
  // Alternatively, could one think of ParameterCFs here?
  double tol{1e-8};
  double rtol{0.0};
  int maxiter{10};

public:
  NewtonCF(shared_ptr<CoefficientFunction> aexpression,
           shared_ptr<CoefficientFunction> astartingpoint,
           std::optional<double> atol, std::optional<double> artol,
           std::optional<int> amaxiter)
      : NewtonCF{aexpression,
                 Array<shared_ptr<CoefficientFunction>>{astartingpoint}, atol,
                 artol, amaxiter} {}

  NewtonCF(shared_ptr<CoefficientFunction> aexpression,
           const Array<shared_ptr<CoefficientFunction>> &astartingpoints,
           std::optional<double> atol, std::optional<double> artol,
           std::optional<int> amaxiter)
      : expression(aexpression) {

    // NOTES:

    // All proxies must originate from one FE space. However, there is no way to
    // check this. Available information for proxies:
    //  - block index: via Proxy Evaluator if this can be cast to a
    //  CompoundDifferentialOperator
    //  - block size: via Proxy Evaluator
    //  - dimension: should be the same as block size?
    //  - dims: array with dimensions for each axis

    // For the CF (starting point):
    //  - dim, dims (misleading for GF of CompoundFESpace)
    //  - components (CFs) for GF of CompoundFESpace
    //  --> provide a sequence of stps or use components of a GF

    // Strategy:
    //  1. If more that one proxy has been found, sort proxies by block index;
    //  only allow each block index appear exactly once.
    //     IOW, two different proxies with for some reason same block index are
    //     forbidden.
    //  2. Determine cumulated dimensions of collected proxies and compare
    //    them with dimensions of startingpoints.
    //  4. Call SetDimensions with the appropriate information.

    // NOTE: GFs on generic CompoundSpaces do not provide useful/usable
    // dimension data! They are currently not supported. For that, newtonCF.cpp
    // would have to be part of "comp" instead of "fem"

    expression->TraverseTree([&](CoefficientFunction &nodecf) {
      auto nodeproxy = dynamic_cast<ProxyFunction *>(&nodecf);
      if (nodeproxy) {
        if (!nodeproxy->IsTestFunction()) {
          if (std::find(cbegin(proxies), cend(proxies), nodeproxy) ==
              cend(proxies))
            proxies.Append(nodeproxy);
        }
      } else if (nodecf.StoreUserData() && !cachecf.Contains(&nodecf))
        cachecf.Append(&nodecf);
    });
    if (proxies.Size() == 0)
      throw Exception("NewtonCF: don't have a proxy");
    else if (proxies.Size() > 1) {
      // Check whether all proxies belong to a compound FE space and put them in
      // order
      Array<ProxyFunction *> sorted_proxies(proxies.Size());
      sorted_proxies = nullptr;
      Array<const DifferentialOperator *> diffops(proxies.Size());
      diffops = nullptr;
      //      cout << "\n" << "sorted proxies " << sorted_proxies
      //           << "\n" << "diffops " << diffops << endl;

      for (const auto proxy : proxies) {
        const auto evaluator =
            dynamic_cast<const CompoundDifferentialOperator *>(
                proxy->Evaluator().get());
        if (!evaluator) {
          throw Exception(
              "NewtonCF: More than one proxy has been found but not all proxy "
              "evaluators are of type CompoundDifferentialOperator");
        } else {
          if (sorted_proxies[evaluator->Component()] ||
              std::find(cbegin(diffops), cend(diffops), evaluator) !=
                  cend(diffops))
            throw Exception("NewtonCF: A proxy evaluator (component) has been "
                            "detected twice");
          diffops[evaluator->Component()] = evaluator;
          sorted_proxies[evaluator->Component()] = proxy;
        }
      }
      // Copy over...
      std::copy(begin(sorted_proxies), end(sorted_proxies), begin(proxies));
    }

    // Process proxy dimensions
    for (const auto proxy : proxies) {
      numeric_dim += proxy_dof_dimension(proxy);
      full_dim += proxy->Dimension();
      const auto pdims = proxy->Dimensions();
      if (pdims.Size() > 0)
        for (auto dim : pdims)
          proxy_dims.Append(dim);
      else
        proxy_dims.Append(proxy->Dimension());
    }

    if (expression->Dimension() != full_dim)
      throw Exception(string("NewtonCF: dimension of residual expression (=") +
                      to_string(expression->Dimension()) +
                      ") does not match the accumulated dimension of detected "
                      "trial functions (=" +
                      to_string(full_dim) + ")");

    // Process startingpoints

    /*
      // if we want that we should move it to comp - level ???
    // Handle GF on CompoundFESpace
    if (astartingpoints.Size() == 1 && astartingpoints[0]->Dimension() == 1 &&
        proxies.Size() > 1) {
      const auto startingpoint_gf =
          dynamic_pointer_cast<ngcomp::GridFunction>(astartingpoints[0]);
      if (!startingpoint_gf)
        throw Exception("NewtonCF: number of trial functions greater than one "
                        "requires a GridFunction with corresponding components "
                        "as starting point");

      if (proxies.Size() != startingpoint_gf->GetNComponents())
        throw Exception(string("NewtonCF: number of proxies (=") +
                        to_string(proxies.Size()) +
                        ") does not match the number "
                        "of components of the 'startingpoint' (=" +
                        to_string(startingpoint_gf->GetNComponents()) + ")");

      startingpoints.DeleteAll();
      for (int i : Range(startingpoint_gf->GetNComponents()))
        startingpoints.Append(startingpoint_gf->GetComponent(i));

    } else
    */
      startingpoints = astartingpoints;

    // Check dimensions and/or fill empty startingpoints
    if (startingpoints.Size() == 0) {
      for (const auto proxy : proxies)
        startingpoints.Append(ZeroCF(proxy->Dimensions()));
    } else if (startingpoints.Size() == proxies.Size()) {
      for (int i : Range(proxies)) {
        if (!startingpoints[i])
          startingpoints[i] = ZeroCF(proxies[i]->Dimensions());
        else if (!(proxies[i]->Dimensions() == startingpoints[i]->Dimensions()))
          throw Exception(std::string("NewtonCF: Dimensions of component ") +
                          std::to_string(i) + " do not agree");
      }
    } else if (startingpoints.Size() == 1) {
      if (startingpoints[0]->Dimension() != full_dim)
        throw Exception(
            string("NewtonCF: Total dimension of startingpoints (=") +
            to_string(startingpoints[0]->Dimension()) +
            ") does not match the accumulated dimension of trial "
            "functions (=" +
            to_string(full_dim) + ")");
    } else
      throw Exception(string("NewtonCF: Number of given startingpoints (=") +
                      to_string(startingpoints.Size()) +
                      ") does not match "
                      "number of detected proxies (=" +
                      to_string(proxies.Size()) + ")");

    if (proxies.Size() == 1)
      CoefficientFunction::SetDimensions(
          FlatArray{proxy_dims.Size(), proxy_dims.Data()});
    else {
      Array<int> dims;
      dims.Append(full_dim);
      CoefficientFunction::SetDimensions(dims);
    }

    // Process options
    if (atol)
      tol = *atol;
    if (artol)
      rtol = *artol;
    if (amaxiter)
      maxiter = *amaxiter;
  }

  double Evaluate(const BaseMappedIntegrationPoint &ip) const override {
    cout << "pw eval not overloaded" << endl;
    return 0;
  }

  void Evaluate(const BaseMappedIntegrationRule &mir,
                BareSliceMatrix<double> values) const override {
    // static Timer t("NewtonCF::Eval", 2);
    // static Timer t1("NewtonCF::Eval get Jac", 2);
    // static Timer t2("NewtonCF::Eval solve", 2);
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // RegionTracer regtr(TaskManager::GetThreadId(), t);

    // cout << "eval minimization" << endl;

    LocalHeap lh(1000000);

    // startingpoint -> Evaluate (mir, values);
    // cout << "starting: " << endl << values.AddSize(Dimension(), ir.Size()) <<
    // endl;

    const ElementTransformation &trafo = mir.GetTransformation();
    auto saved_ud = trafo.PushUserData();

    ProxyUserData ud(proxies.Size(), cachecf.Size(), lh);
    for (CoefficientFunction *cf : cachecf)
      ud.AssignMemory(cf, mir.Size(), cf->Dimension(), lh);

    const_cast<ElementTransformation &>(trafo).userdata = &ud;

    for (ProxyFunction *proxy : proxies)
      ud.AssignMemory(proxy, mir.Size(), proxy->Dimension(), lh);

    // Prepare data structures for blocks
    const auto nblocks = proxies.Size();
    FlatArray<FlatMatrix<double>> xk_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> w_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> xold_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> val_blocks(nblocks, lh);
    FlatArray<FlatMatrix<AutoDiff<1, double>>> dval_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> deriv_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> res_blocks(nblocks, lh);
    FlatArray<FlatTensor<3>> lin_blocks(nblocks * nblocks, lh);

    // These are only "independent" for blocks having "vsemb"; otherwise just
    // views
    FlatArray<FlatMatrix<double>> rhs_blocks(nblocks, lh);
    FlatArray<FlatTensor<3>> lhs_blocks(nblocks * nblocks, lh);

    for (int i : Range(nblocks)) {
      const auto proxy = proxies[i];
      xk_blocks[i].Assign(ud.GetMemory(proxy));
      xold_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
      w_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
      val_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
      dval_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
      deriv_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
      res_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);

      if (has_vs_embedding(proxy))
        rhs_blocks[i].AssignMemory(mir.Size(), proxy_dof_dimension(proxy), lh);
      else
        rhs_blocks[i].AssignMemory(mir.Size(), proxy_dof_dimension(proxy),
                                   res_blocks[i].Data());

      for (int j : Range(nblocks)) {
        const auto ij = i * nblocks + j;
        lin_blocks[ij].AssignMemory(lh, mir.Size(), proxies[i]->Dimension(),
                                    proxies[j]->Dimension());
        if (has_vs_embedding(proxies[i]) || has_vs_embedding(proxies[j]))
          lhs_blocks[ij].AssignMemory(lh, mir.Size(),
                                      proxy_dof_dimension(proxies[i]),
                                      proxy_dof_dimension(proxies[j]));
        else
          lhs_blocks[ij].Assign(lin_blocks[ij]);
      }
    }

    // Block-agnostic data structures
    FlatMatrix<> xk(mir.Size(), full_dim, lh);
    FlatMatrix<> xold(mir.Size(), full_dim, lh);
    FlatMatrix<> w(mir.Size(), full_dim, lh);
    FlatMatrix<> val(mir.Size(), full_dim, lh);
    FlatMatrix<AutoDiff<1, double>> dval(mir.Size(), full_dim, lh);

    FlatVector<> rhs(numeric_dim, lh);
    FlatVector<> sol(numeric_dim, lh);
    FlatMatrix<> lhs(numeric_dim, numeric_dim, lh);

    const auto converged = [&](const auto &rhs_vec, double res_0 = 0) {
      const auto res = LInfNorm(rhs_vec);
      return res <= tol || (res_0 > 0 && (res / res_0) <= rtol);
    };

    const auto all_converged = [&](const auto &rhs_blocks, double res_0 = 0) {
      return std::all_of(begin(rhs_blocks), end(rhs_blocks),
                         [=](const auto &block) {
                           return converged(block.AsVector(), res_0);
                         });
    };

    const auto distribute_vec_to_blocks = [&](const auto &src,
                                              auto &dest) -> void {
      for (size_t qi : Range(mir)) {
        auto src_qi = src.Row(qi);
        int offset = 0;
        for (int block : Range(nblocks)) {
          auto dest_qi = dest[block].Row(qi);
          dest_qi = src_qi.Range(offset, offset + dest_qi.Size());
          offset += dest_qi.Size();
        }
      }
    };

    const auto calc_residuals = [&]() -> void {
      // RegionTracer regtr1(TaskManager::GetThreadId(), t1);
      expression->Evaluate(mir, val);
      distribute_vec_to_blocks(val, val_blocks);

      for (int block : Range(nblocks)) {
        auto proxy = proxies[block];
        auto &res = res_blocks[block];
        auto &rhs = rhs_blocks[block];
        auto &valb = val_blocks[block];

        for (int k : Range(proxy->Dimension()))
          for (size_t qi : Range(mir))
            res(qi, k) = valb(qi, k);

        // The actual rhs (respecting VS embeddings)
        if (auto vsemb = get_vs_embedding(proxy); vsemb)
          for (size_t qi : Range(mir))
            rhs.Row(qi) = Trans(vsemb.value()) * res.Row(qi);

        // cout << "res block " << block << " = " << res << endl;
      }
    };

    const auto calc_linearizations = [&]() -> void {
      // RegionTracer regtr2(TaskManager::GetThreadId(), t2);
      for (int block2 : Range(nblocks)) {
        auto proxy2 = proxies[block2];

        for (int l : Range(proxy2->Dimension())) {
          // This means, the derivative is taken wrt proxy 2
          ud.trialfunction = proxy2;
          ud.trial_comp = l;
          expression->Evaluate(mir, dval);
          distribute_vec_to_blocks(dval, dval_blocks);

          for (int block1 : Range(nblocks)) {
            auto proxy1 = proxies[block1];
            auto &dvalb = dval_blocks[block1];
            auto &deriv = deriv_blocks[block1];
            auto &lin = lin_blocks[block1 * nblocks + block2];
            for (int k : Range(proxy1->Dimension())) {
              for (size_t qi : Range(mir)) {
                deriv(qi, k) = dvalb(qi, k).DValue(0);
              }
            }
            lin(STAR, STAR, l) = deriv;
            // cout << "lin block (" << block1 << ", " << block2 << ")[*, *, "
            //      << l << "] = " << lin(STAR, STAR, l) << endl;
          }
          // cout << "lin block (" << block1 << ", " << block2
          //      << ") = " << lin << endl;
        }
      }
    };

    const auto compute_increments = [&]() -> void {
      for (size_t qi : Range(mir)) {

        // NOTE: when to skip something because of convergence?
        //  -> the current approach assumes that evaluation of "expression"
        //  for all qpoints at once is beneficial!

        int offset1 = 0;
        for (int block1 : Range(proxies)) {
          const auto proxy1 = proxies[block1];
          const auto &rhsb = rhs_blocks[block1].Row(qi);

          for (int k : Range(rhsb.Size()))
            rhs[offset1 + k] = rhsb[k];

          int offset2 = 0;
          for (int block2 : Range(proxies)) {
            const auto proxy2 = proxies[block2];
            const auto &linb =
                lin_blocks[block1 * nblocks + block2](qi, STAR, STAR);
            const auto &lhsb =
                lhs_blocks[block1 * nblocks + block2](qi, STAR, STAR);
            if (auto vsemb1 = get_vs_embedding(proxy1); vsemb1)
              if (auto vsemb2 = get_vs_embedding(proxy2); vsemb2) {
                lhsb = Trans(vsemb1.value()) * linb * vsemb2.value();
              } else {
                lhsb = Trans(vsemb1.value()) * linb;
              }
            else if (auto vsemb2 = get_vs_embedding(proxy2); vsemb2) {
              lhsb = linb * vsemb2.value();
            } else {
              // nothing to do
            }

            // cout << "lhs block (" << block1 << ", " << block2 << ") = "
            //      << lhsb << endl;
            for (int k : Range(lhsb.Height()))
              for (int l : Range(lhsb.Width()))
                lhs(offset1 + k, offset2 + l) = lhsb(k, l);

            offset2 += lhsb.Width();
          }
          offset1 += rhsb.Size();
        }

        if (converged(rhs)) {
          w.Row(qi) = 0;
          continue;
        }

        //        cout << "RHS: " << rhs << endl;
        //        cout << "LHS: " << lhs << endl;
        CalcInverse(lhs);
        sol = lhs * rhs;

        // Handle VS-embedding
        auto &wi = w.Row(qi);
        int offset_w = 0;
        int offset_sol = 0;
        for (int block : Range(proxies)) {
          const auto proxy = proxies[block];
          if (const auto vsemb = get_vs_embedding(proxy); vsemb)
            wi.Range(offset_w, offset_w + proxy->Dimension()) =
                vsemb.value() *
                sol.Range(offset_sol, offset_sol + proxy_dof_dimension(proxy));
          else
            wi.Range(offset_w, offset_w + proxy->Dimension()) =
                sol.Range(offset_sol, offset_sol + proxy_dof_dimension(proxy));
          offset_w += proxy->Dimension();
          offset_sol += proxy_dof_dimension(proxy);
        }
      }
    };

    const auto merge_xk_blocks = [&]() -> void {
      for (size_t qi : Range(mir)) {
        auto xk_qi = xk.Row(qi);
        int offset = 0;
        for (int block : Range(this->proxies)) {
          auto xkb_qi = xk_blocks[block].Row(qi);
          xk_qi.Range(offset, offset + xkb_qi.Size()) = xkb_qi;
          offset += xkb_qi.Size();
        }
      }
    };

    // Evaluate starting point
    if (startingpoints.Size() == proxies.Size()) {
      for (int i : Range(startingpoints))
        startingpoints[i]->Evaluate(mir, xk_blocks[i]);
      merge_xk_blocks();
    } else {
      // This is the only case when we should end up here (checked in
      // constructor)
      assert(startingpoints.Size() == 1);

      startingpoints[0]->Evaluate(mir, xk);
      distribute_vec_to_blocks(xk, xk_blocks);
    }

    // cout << "starting value = " << xk << endl;
    // cout << "blocks:" << endl;
    // for (int i : Range(proxies))
    //   cout << i << " -> " << ud.GetMemory(proxies[i]) << endl;
    // cout << endl;

    // The actual algorithm
    //    cout <<
    //    "\n\n------------------------START------------------------------"
    //            "\n";

    calc_residuals();

    //    cout << "(pre) rhs blocks: " << rhs_blocks;

    for (int step : Range(maxiter)) {
      if (all_converged(rhs_blocks))
        break;

      calc_linearizations();
      compute_increments();

      xk -= w;
      // cout << "xk: " << xk << endl;
      distribute_vec_to_blocks(xk, xk_blocks);
      calc_residuals();
      //      cout << "\nstep: " << step << "\n"
      //           << "rhs blocks: " << rhs_blocks;
    }

    //    cout << "rhs blocks (final): " << rhs_blocks;
    // cout << "xk (final): " << xk << endl;

    if (!all_converged(rhs_blocks))
      xk = numeric_limits<double>::quiet_NaN();

    // cout << "result = " << xk << endl;
    values.AddSize(mir.Size(), full_dim) = xk;
    //    cout <<
    //    "\n-------------------------Done------------------------------\n";
  }
};

class MinimizationCF : public CoefficientFunction {

  shared_ptr<CoefficientFunction> expression;
  Array<shared_ptr<CoefficientFunction>> startingpoints{};

  Array<ProxyFunction *> proxies{};
  Array<CoefficientFunction *> cachecf{};

  // The total dimension of the linear system (because of VS embeddings, this
  // can be different from this->Dimension()).
  int numeric_dim = 0;
  int full_dim = 0;
  Array<int> proxy_dims{};

  // Same parameters as for scipy's newton
  // Alternatively, could one think of ParameterCFs here?
  double tol{1e-8};
  double rtol{0.0};
  int maxiter{10};

public:
  MinimizationCF(shared_ptr<CoefficientFunction> aexpression,
                 shared_ptr<CoefficientFunction> astartingpoint,
                 std::optional<double> atol, std::optional<double> artol,
                 std::optional<int> amaxiter)
      : MinimizationCF{aexpression,
                       Array<shared_ptr<CoefficientFunction>>{astartingpoint},
                       atol, artol, amaxiter} {}

  MinimizationCF(shared_ptr<CoefficientFunction> aexpression,
                 const Array<shared_ptr<CoefficientFunction>> &astartingpoints,
                 std::optional<double> atol, std::optional<double> artol,
                 std::optional<int> amaxiter)
      : expression(aexpression) {

    expression->TraverseTree([&](CoefficientFunction &nodecf) {
      auto nodeproxy = dynamic_cast<ProxyFunction *>(&nodecf);
      if (nodeproxy) {
        if (!nodeproxy->IsTestFunction()) {
          if (std::find(cbegin(proxies), cend(proxies), nodeproxy) ==
              cend(proxies))
            proxies.Append(nodeproxy);
        }
      } else if (nodecf.StoreUserData() && !cachecf.Contains(&nodecf))
        cachecf.Append(&nodecf);
    });
    if (proxies.Size() == 0)
      throw Exception("MinimizationCF: don't have a proxy");
    else if (proxies.Size() > 1) {
      // Check whether all proxies belong to a compound FE space and put them in
      // order
      Array<ProxyFunction *> sorted_proxies(proxies.Size());
      sorted_proxies = nullptr;
      Array<const DifferentialOperator *> diffops(proxies.Size());
      diffops = nullptr;
      //      cout << "\n" << "sorted proxies " << sorted_proxies
      //           << "\n" << "diffops " << diffops << endl;

      for (const auto proxy : proxies) {
        const auto evaluator =
            dynamic_cast<const CompoundDifferentialOperator *>(
                proxy->Evaluator().get());
        if (!evaluator) {
          throw Exception(
              "MinimizationCF: More than one proxy has been found but not all "
              "proxy "
              "evaluators are of type CompoundDifferentialOperator");
        } else {
          if (sorted_proxies[evaluator->Component()] ||
              std::find(cbegin(diffops), cend(diffops), evaluator) !=
                  cend(diffops))
            throw Exception(
                "MinimizationCF: A proxy evaluator (component) has been "
                "detected twice");
          diffops[evaluator->Component()] = evaluator;
          sorted_proxies[evaluator->Component()] = proxy;
        }
      }
      // Copy over...
      std::copy(begin(sorted_proxies), end(sorted_proxies), begin(proxies));
    }

    // Process proxy dimensions
    for (const auto proxy : proxies) {
      numeric_dim += proxy_dof_dimension(proxy);
      full_dim += proxy->Dimension();
      const auto pdims = proxy->Dimensions();
      if (pdims.Size() > 0)
        for (auto dim : pdims)
          proxy_dims.Append(dim);
      else
        proxy_dims.Append(proxy->Dimension());
    }

    if (expression->Dimension() != 1)
      throw Exception(string("MinimizationCF: only scalar expressions are allowed"));

    // Process startingpoints

    /*
      // then we shold move to to comp level ???
    // Handle GF on CompoundFESpace
    if (astartingpoints.Size() == 1 && astartingpoints[0]->Dimension() == 1 &&
        proxies.Size() > 1) {
        const auto startingpoint_gf =
            dynamic_pointer_cast<ngcomp::GridFunction>(astartingpoints[0]);
        if (!startingpoint_gf)
          throw Exception(
              "MinimizationCF: number of trial functions greater than one "
              "requires a GridFunction with corresponding components "
              "as starting point");

        if (proxies.Size() != startingpoint_gf->GetNComponents())
          throw Exception(string("MinimizationCF: number of proxies (=") +
                          to_string(proxies.Size()) +
                          ") does not match the number "
                          "of components of the 'startingpoint' (=" +
                          to_string(startingpoint_gf->GetNComponents()) + ")");

        startingpoints.DeleteAll();
        for (int i : Range(startingpoint_gf->GetNComponents()))
          startingpoints.Append(startingpoint_gf->GetComponent(i));

    } else
    */
      startingpoints = astartingpoints;

    // Check dimensions and/or fill empty startingpoints
    if (startingpoints.Size() == 0) {
        for (const auto proxy : proxies)
          startingpoints.Append(ZeroCF(proxy->Dimensions()));
    } else if (startingpoints.Size() == proxies.Size()) {
        for (int i : Range(proxies)) {
          if (!startingpoints[i])
            startingpoints[i] = ZeroCF(proxies[i]->Dimensions());
          else if (!(proxies[i]->Dimensions() ==
                     startingpoints[i]->Dimensions()))
            throw Exception(
                std::string("MinimizationCF: Dimensions of startingpoint "
                            "and proxy component ") +
                std::to_string(i) + " do not agree");
        }
    } else if (startingpoints.Size() == 1) {
        if (startingpoints[0]->Dimension() != full_dim)
          throw Exception(
              string("MinimizationCF: Total dimension of startingpoints (=") +
              to_string(startingpoints[0]->Dimension()) +
              ") does not match the accumulated dimension of trial "
              "functions (=" +
              to_string(full_dim) + ")");
    } else
      throw Exception(string("MinimizationCF: Number of given startingpoints (=") +
                      to_string(startingpoints.Size()) +
                      ") does not match "
                      "number of detected proxies (=" +
                      to_string(proxies.Size()) + ")");

    if (proxies.Size() == 1)
      CoefficientFunction::SetDimensions(
          FlatArray{proxy_dims.Size(), proxy_dims.Data()});
    else {
      Array<int> dims;
      dims.Append(full_dim);
      CoefficientFunction::SetDimensions(dims);
    }

    // Process options
    if (atol)
      tol = *atol;
    if (artol)
      rtol = *artol;
    if (amaxiter)
      maxiter = *amaxiter;
  }

  double Evaluate(const BaseMappedIntegrationPoint &ip) const override {
    cout << "pw eval not overloaded" << endl;
    return 0;
  }

  void Evaluate(const BaseMappedIntegrationRule &mir,
                BareSliceMatrix<double> values) const override {
    // static Timer t("MinimizationCF::Eval", 2);
    // static Timer t1("MinimizationCF::Eval get Jac", 2);
    // static Timer t2("MinimizationCF::Eval solve", 2);
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // RegionTracer regtr(TaskManager::GetThreadId(), t);

    // cout << "eval minimization" << endl;

    LocalHeap lh(1000000);

    const ElementTransformation &trafo = mir.GetTransformation();
    auto saved_ud = trafo.PushUserData();

    ProxyUserData ud(proxies.Size(), cachecf.Size(), lh);
    for (CoefficientFunction *cf : cachecf)
      ud.AssignMemory(cf, mir.Size(), cf->Dimension(), lh);

    const_cast<ElementTransformation &>(trafo).userdata = &ud;

    for (ProxyFunction *proxy : proxies)
      ud.AssignMemory(proxy, mir.Size(), proxy->Dimension(), lh);

    // Prepare data structures for blocks
    const auto nblocks = proxies.Size();
    FlatArray<FlatMatrix<double>> xk_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> w_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> xold_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> diags_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> res_blocks(nblocks, lh);
    FlatArray<FlatTensor<3>> lin_blocks(nblocks * nblocks, lh);

    // These are only "independent" for blocks having "vsemb"; otherwise just
    // views
    FlatArray<FlatMatrix<double>> rhs_blocks(nblocks, lh);
    FlatArray<FlatTensor<3>> lhs_blocks(nblocks * nblocks, lh);

    for (int i : Range(nblocks)) {
      const auto proxy = proxies[i];
      xk_blocks[i].Assign(ud.GetMemory(proxy));
      xold_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
      w_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
      diags_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
      res_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);

      if (has_vs_embedding(proxy))
        rhs_blocks[i].AssignMemory(mir.Size(), proxy_dof_dimension(proxy), lh);
      else
        rhs_blocks[i].AssignMemory(mir.Size(), proxy_dof_dimension(proxy),
                                   res_blocks[i].Data());

      for (int j : Range(nblocks)) {
        const auto ij = i * nblocks + j;
        lin_blocks[ij].AssignMemory(lh, mir.Size(), proxies[i]->Dimension(),
                                    proxies[j]->Dimension());
        if (has_vs_embedding(proxies[i]) || has_vs_embedding(proxies[j]))
          lhs_blocks[ij].AssignMemory(lh, mir.Size(),
                                      proxy_dof_dimension(proxies[i]),
                                      proxy_dof_dimension(proxies[j]));
        else
          lhs_blocks[ij].Assign(lin_blocks[ij]);
      }
    }

    // Block-agnostic data structures
    double energy = 0;
    FlatMatrix<> dderiv(mir.Size(), 1, lh);
    FlatMatrix<AutoDiffDiff<1, double>> ddval(mir.Size(), 1, lh);
    FlatMatrix<> xk(mir.Size(), full_dim, lh);
    FlatMatrix<> xold(mir.Size(), full_dim, lh);
    FlatMatrix<> w(mir.Size(), full_dim, lh);
    FlatVector<> rhs(numeric_dim, lh);
    FlatVector<> sol(numeric_dim, lh);
    FlatMatrix<> lhs(numeric_dim, numeric_dim, lh);

    const auto converged = [&](const auto &rhs_vec, double res_0 = 0) {
      const auto res = LInfNorm(rhs_vec);
      return res <= tol || (res_0 > 0 && (res / res_0) <= rtol);
    };

    const auto all_converged = [&](const auto &rhs_blocks, double res_0 = 0) {
      return std::all_of(begin(rhs_blocks), end(rhs_blocks),
                         [=](const auto &block) {
                           return converged(block.AsVector(), res_0);
                         });
    };

    const auto distribute_vec_to_blocks = [&](const auto &src,
                                              auto &dest) -> void {
      for (size_t qi : Range(mir)) {
        auto src_qi = src.Row(qi);
        int offset = 0;
        for (int block : Range(nblocks)) {
          auto dest_qi = dest[block].Row(qi);
          dest_qi = src_qi.Range(offset, offset + dest_qi.Size());
          offset += dest_qi.Size();
        }
      }
    };


    const auto calc_energy_rhs_and_diags = [&]() -> void {
      energy = 0;
      for (int block : Range(nblocks)) {
        auto proxy = proxies[block];
        auto &res = res_blocks[block];
        auto &rhsb = rhs_blocks[block];
        auto &diags = diags_blocks[block];

        for (int k : Range(proxy->Dimension())) {
          ud.trialfunction = proxy;
          ud.trial_comp = k;
          ud.testfunction = proxy;
          ud.test_comp = k;

          expression->Evaluate(mir, ddval);
          for (size_t qi : Range(mir))
            diags(qi, k) = ddval(qi, 0).DDValue(0);
          for (size_t qi : Range(mir))
            res(qi, k) = ddval(qi, 0).DValue(0);

          if (k == 0 && block == 0)
            for (size_t qi : Range(mir))
              energy += mir[qi].GetWeight() * ddval(qi, 0).Value();
        }

        // The actual rhs (respecting VS embeddings)
        if (auto vsemb = get_vs_embedding(proxy); vsemb)
          for (size_t qi : Range(mir))
            rhsb.Row(qi) = Trans(vsemb.value()) * res.Row(qi);

        // cout << "res block " << block << " = " << res << endl;
      }
    };


    const auto calc_off_diagonals = [&]() -> void {
      // TODO: exploit symmetry
      for (int l1 : Range(nblocks)) {
        for (int k1 : Range(nblocks)) {

          auto proxy1 = proxies[k1];
          auto proxy2 = proxies[l1];
          auto &lin = lin_blocks[l1 * nblocks + k1];

          for (int l : Range(proxy2->Dimension())) {
            for (int k : Range(proxy1->Dimension())) {

              ud.testfunction = proxy2;
              ud.test_comp = l;
              ud.trialfunction = proxy1;
              ud.trial_comp = k;

              expression->Evaluate(mir, ddval);
              for (size_t i : Range(mir)) {
                dderiv(i, 0) = ddval(i, 0).DDValue(0);
              }
              lin(STAR, l, k) = dderiv.Col(0);

              if (proxy1 != proxy2 ||
                  k != l) // computed mixed second derivatives
              {
                lin(STAR, l, k) -= diags_blocks[k1].Col(k);
                lin(STAR, l, k) -= diags_blocks[l1].Col(l);
                lin(STAR, l, k) *= 0.5;
              }

              // cout << "lin block (" << k1 << ", " << l1 << ")[*, *, "
              //      << l << "] = " << lin(STAR, STAR, l) << endl;
            }
            // cout << "lin block (" << k1 << ", " << l1
            //      << ") = " << lin << endl;
          }
        }
      }
    };

    const auto compute_newton_step = [&]() -> void {
      for (size_t qi : Range(mir)) {

        // NOTE: when to skip something because of convergence?
        //  -> the current approach assumes that evaluation of "expression"
        //  for all qpoints at once is beneficial!

        int offset1 = 0;
        for (int block1 : Range(proxies)) {
          const auto proxy1 = proxies[block1];
          const auto &rhsb = rhs_blocks[block1].Row(qi);

          for (int k : Range(rhsb.Size()))
            rhs[offset1 + k] = rhsb[k];

          // TODO: exploit symmetry
          int offset2 = 0;
          for (int block2 : Range(proxies)) {
            const auto proxy2 = proxies[block2];
            const auto &linb =
                lin_blocks[block1 * nblocks + block2](qi, STAR, STAR);
            const auto &lhsb =
                lhs_blocks[block1 * nblocks + block2](qi, STAR, STAR);
            if (auto vsemb1 = get_vs_embedding(proxy1); vsemb1)
              if (auto vsemb2 = get_vs_embedding(proxy2); vsemb2) {
                lhsb = Trans(vsemb1.value()) * linb * vsemb2.value();
              } else {
                lhsb = Trans(vsemb1.value()) * linb;
              }
            else if (auto vsemb2 = get_vs_embedding(proxy2); vsemb2) {
              lhsb = linb * vsemb2.value();
            } else {
              // nothing to do
            }

//            cout << "lhs block (" << block1 << ", " << block2 << ") = "
//                 << lhsb << endl;

            for (int k : Range(lhsb.Height()))
              for (int l : Range(lhsb.Width()))
                lhs(offset1 + k, offset2 + l) = lhsb(k, l);

            offset2 += lhsb.Width();
          }
          offset1 += rhsb.Size();
        }

        if (converged(rhs)) {
          w.Row(qi) = 0;
          continue;
        }

//        cout << "RHS: " << rhs << endl;
//        cout << "LHS: " << lhs << endl;
        CalcInverse(lhs);
        sol = lhs * rhs;

        // Handle VS-embedding
        auto &wi = w.Row(qi);
        int offset_w = 0;
        int offset_sol = 0;
        for (int block : Range(proxies)) {
          const auto proxy = proxies[block];
          if (const auto vsemb = get_vs_embedding(proxy); vsemb)
            wi.Range(offset_w, offset_w + proxy->Dimension()) =
                vsemb.value() *
                sol.Range(offset_sol, offset_sol + proxy_dof_dimension(proxy));
          else
            wi.Range(offset_w, offset_w + proxy->Dimension()) =
                sol.Range(offset_sol, offset_sol + proxy_dof_dimension(proxy));
          offset_w += proxy->Dimension();
          offset_sol += proxy_dof_dimension(proxy);
        }
      }
    };

    const auto linesearch = [&](auto &ud) -> void {
      xold = xk;// linesearch
      double alpha = 1;
      double newenergy = energy + 1;

      auto proxy = proxies[0];
      ud.trialfunction = proxy;
      ud.trial_comp = 0;
      ud.testfunction = proxy;
      ud.test_comp = 0;

      // cout << "w = " << endl << w << endl;
      while (newenergy > energy && alpha > 1e-10) {
        xk = xold - alpha * w;
        distribute_vec_to_blocks(xk, xk_blocks);

        newenergy = 0;
        expression->Evaluate(mir, ddval);
        for (size_t qi : Range(mir)) {
          newenergy += mir[qi].GetWeight() * ddval(qi, 0).Value();
        }

        // cout << "alpha = " << alpha << ", newen = " << newenergy << endl;
        alpha /= 2;
      }
    };

    const auto merge_xk_blocks = [&]() -> void {
      for (size_t qi : Range(mir)) {
        auto xk_qi = xk.Row(qi);
        int offset = 0;
        for (int block : Range(this->proxies)) {
          auto xkb_qi = xk_blocks[block].Row(qi);
          xk_qi.Range(offset, offset + xkb_qi.Size()) = xkb_qi;
          offset += xkb_qi.Size();
        }
      }
    };

    // Evaluate starting point
    if (startingpoints.Size() == proxies.Size()) {
      for (int i : Range(startingpoints))
        startingpoints[i]->Evaluate(mir, xk_blocks[i]);
      merge_xk_blocks();
    } else {
      // This is the only case when we should end up here (checked in
      // constructor)
      assert(startingpoints.Size() == 1);

      startingpoints[0]->Evaluate(mir, xk);
      distribute_vec_to_blocks(xk, xk_blocks);
    }

    // The actual algorithm
//    cout << "\n" << "start newton loop" << "\n";
    calc_energy_rhs_and_diags();
    for (int step = 0; step < maxiter; step++) {
      if (all_converged(rhs_blocks))
        break;

      calc_off_diagonals();
      compute_newton_step();
      linesearch(ud);
      calc_energy_rhs_and_diags();
//      cout << "newton step " << step + 1 << endl;
    }

//    for (int block1 : Range(proxies))
//      cout << "RHS block " << to_string(block1) << ": " << rhs_blocks[block1] << endl;
//
//    cout << "done" << "\n\n";

    if (!all_converged(rhs_blocks))
      xk = numeric_limits<double>::quiet_NaN();

    // cout << "result = " << xk << endl;
    values.AddSize(mir.Size(), Dimension()) = xk;
  }
};

shared_ptr<CoefficientFunction>
CreateMinimizationCF(shared_ptr<CoefficientFunction> expression,
               shared_ptr<CoefficientFunction> startingpoint,
               std::optional<double> atol, std::optional<double> rtol,
               std::optional<int> maxiter) {
  return make_shared<MinimizationCF>(expression, startingpoint, atol, rtol, maxiter);
}

shared_ptr<CoefficientFunction>
CreateMinimizationCF(shared_ptr<CoefficientFunction> expression,
               const Array<shared_ptr<CoefficientFunction>> &startingpoints,
               std::optional<double> tol, std::optional<double> rtol,
               std::optional<int> maxiter) {
  return make_shared<MinimizationCF>(expression, startingpoints, tol, rtol, maxiter);
}

shared_ptr<CoefficientFunction>
CreateNewtonCF(shared_ptr<CoefficientFunction> expression,
               shared_ptr<CoefficientFunction> startingpoint,
               std::optional<double> atol, std::optional<double> rtol,
               std::optional<int> maxiter) {
  return make_shared<NewtonCF>(expression, startingpoint, atol, rtol, maxiter);
}

shared_ptr<CoefficientFunction>
CreateNewtonCF(shared_ptr<CoefficientFunction> expression,
               const Array<shared_ptr<CoefficientFunction>> &startingpoints,
               std::optional<double> tol, std::optional<double> rtol,
               std::optional<int> maxiter) {
  return make_shared<NewtonCF>(expression, startingpoints, tol, rtol, maxiter);
}

} // namespace ngfem
