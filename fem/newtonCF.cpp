/*********************************************************************/
/* File:   newtonCF.cpp                                              */
/* Authors: Joachim Schoeberl, Matthias Rambausek                    */
/* Date:   Feb 2021                                                  */
/*********************************************************************/

#include <cmath>
#include <fem.hpp>
#include <ngstd.hpp>


namespace std {
    
// iterator traits: for compatibility with the std c++ iterator library
// This is only a workaround for presently missing bits in netgen's array.hpp.
template<typename T>
struct iterator_traits<ngcore::ArrayRangeIterator<T>> {
  using difference_type = T;
  using value_type = T;
  using pointer_type = void;
  using reference = T;
  using iterator_category = input_iterator_tag;
};

}

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

template<typename vec_t>
bool converged(const vec_t &vec, double tol, double res_0 = 0, double rtol = 0) {
  const double res = LInfNorm(vec);
  return res <= tol || (res_0 > 0 && (res / res_0) <= rtol);
}

template<typename mat_t, typename vec_t>
bool all_converged_qp(const mat_t &res, double tol, vec_t res_0 = 0, double rtol = 0) {
  auto qp_range = Range(res.Height());
  return std::all_of(begin(qp_range), end(qp_range),
                     [=](const auto &qi) {
                         return converged(res.Row(qi), tol, res_0[qi], rtol);
                     });
}

} // namespace

class NewtonCF : public CoefficientFunction {
  shared_ptr<CoefficientFunction> expression;
  Array<shared_ptr<CoefficientFunction>> startingpoints{};

  Array<ProxyFunction *> proxies{};
  Array<CoefficientFunction *> cachecf{};

  // The total dimension of the linear system (because of VS embeddings, this
  // can be different from this->Dimension())
  size_t numeric_dim = 0;
  size_t full_dim = 0;
  Array<int> proxy_dims{};

  // Dimension of equation
  size_t eq_dim = 0;

  // Same parameters as for scipy's newton
  // Alternatively, could one think of ParameterCFs here?
  double tol{1e-6};
  double rtol{0.0};
  int maxiter{10};

  bool allow_fail{false};

public:
  NewtonCF(shared_ptr<CoefficientFunction> aexpression,
           shared_ptr<CoefficientFunction> astartingpoint,
           std::optional<double> atol, std::optional<double> artol,
           std::optional<int> amaxiter, std::optional<bool> aallow_fail)
      : NewtonCF{aexpression,
                 Array<shared_ptr<CoefficientFunction>>{astartingpoint}, atol,
                 artol, amaxiter, aallow_fail} {}

  NewtonCF(shared_ptr<CoefficientFunction> aexpression,
           const Array<shared_ptr<CoefficientFunction>> &astartingpoints,
           std::optional<double> atol, std::optional<double> artol,
           std::optional<int> amaxiter, std::optional<bool> aallow_fail)
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
    //  1. If more than one proxy has been found, sort proxies by block index;
    //  only allow each block index appear exactly once.
    //     IOW, two different proxies with for some reason same block index are
    //     forbidden.
    //  2. Determine cumulated dimensions of collected proxies and compare
    //    them with dimensions of startingpoints.
    //  4. Call SetDimensions with the appropriate information.

    // NOTE: GFs on generic CompoundSpaces do not provide useful/usable
    // dimension data! They are currently not supported. For that, newtonCF.cpp
    // would have to be part of "comp" instead of "fem"

    eq_dim = expression->Dimension();

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
      std::map<int, ProxyFunction *> sorted_proxies{};
      std::map<int, const DifferentialOperator *> diffops{};
      //      cout << "\n" << "sorted proxies " << sorted_proxies
      //           << "\n" << "diffops " << diffops << endl;

      const auto compound_space = [](const auto& proxy) {
          return static_cast<void*>(proxy->GetFESpace().get());
      };
      const auto space = compound_space(proxies[0]);
      if (!space)
          throw Exception(
                  "NewtonCF: More than one proxy has been found but not all "
                  "belong to a CompoundFESpace.");

      for (const auto proxy : proxies) {
          if (space != compound_space(proxy))
              throw Exception(
                      "NewtonCF: More than one proxy has been found but not all "
                      "belong to the same a CompoundFESpace.");

          const auto evaluator =
                  dynamic_cast<const CompoundDifferentialOperator *>(
                          proxy->Evaluator().get());
          if (!evaluator) {
              throw Exception(
                      "NewtonCF: More than one proxy has been found but not all "
                      "proxy "
                      "evaluators are of type CompoundDifferentialOperator");
          } else {
              if (sorted_proxies[evaluator->Component()] ||
                  std::find_if(cbegin(diffops), cend(diffops), [&](const auto& pair){return pair.second == evaluator;}) !=
                  cend(diffops))
                  throw Exception(
                          "NewtonCF: A proxy evaluator (component) has been "
                          "detected twice");
              diffops[evaluator->Component()] = evaluator;
              sorted_proxies[evaluator->Component()] = proxy;
          }
      }
      // Copy over...
      for (int i1 = 0, i2 = 0; i1 < sorted_proxies.size(); ++i1)
          if  (auto proxy = sorted_proxies[i1])
              proxies[i2++] = proxy;
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

    if (eq_dim < numeric_dim )
      throw Exception(string("NewtonCF: under-determined system of equations detected. "
                             "Dimension of residual expression (=") + to_string(expression->Dimension()) + ")"
                             "is less than the (numeric/dof) dimension of the detected trial functions "
                             "(=" + to_string(full_dim) + ")");
    else if (eq_dim > full_dim)
      cout << IM(3) << "NewtonCF: Over-determined system detected. "
                       "Dimension of residual exceeds (full/symbolic) dimension of trial functions. "
                       "Linear(ized) systems of equations will be solved in a least-square sense." << endl;
    else if (eq_dim > numeric_dim)
        cout << IM(3) << "NewtonCF: Over-determined system detected."
                         "Dimension of residual exceeds (numeric/dof) dimension of trial functions. "
                         "Linear(ized) systems of equations will be solved in a least-square sense." << endl;

    // Process startingpoints

    startingpoints = astartingpoints;

    // Check dimensions and/or fill empty startingpoints
    if (startingpoints.Size() == 0) {
      for (const auto proxy : proxies)
        startingpoints.Append(ZeroCF(proxy->Dimensions()));
    } else if (startingpoints.Size() == proxies.Size()) {
      for (int i : Range(proxies)) {
        if (!startingpoints[i])
          startingpoints[i] = ZeroCF(proxies[i]->Dimensions());
        else if (!(proxies[i]->Dimensions() == startingpoints[i]->Dimensions())) {
            std::stringstream sstr{};
            sstr << "NewtonCF: Dimensions of startingpoint "
                    "and proxy component " << i << " do not agree ("
                 << startingpoints[i]->Dimensions() << " != " << proxies[i]->Dimensions();
            throw Exception(sstr.str());
        }
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
    if (aallow_fail)
      allow_fail = *aallow_fail;
  }

  double Evaluate(const BaseMappedIntegrationPoint &ip) const override {
    cout << "pw eval not overloaded" << endl;
    return 0;
  }

  void Evaluate(const BaseMappedIntegrationRule &mir,
                BareSliceMatrix<double> values) const override {
    static Timer t("NewtonCF::Eval", NoTracing);
    static Timer t1("NewtonCF::Eval get Jac", NoTracing);
    static Timer t2("NewtonCF::Eval solve", NoTracing);
    static Timer t3("NewtonCF::compute increment", NoTracing);    
    RegionTimer reg(t);

    LocalHeap lh(1000000);

    const ElementTransformation &trafo = mir.GetTransformation();
    auto saved_ud = trafo.PushUserData();

    ProxyUserData ud(proxies.Size(), cachecf.Size(), lh);
    for (CoefficientFunction *cf : cachecf)
      ud.AssignMemory(cf, mir.Size(), cf->Dimension(), lh);

    const_cast<ElementTransformation &>(trafo).userdata = &ud;

    for (ProxyFunction *proxy : proxies)
      ud.AssignMemory(proxy, mir.Size(), proxy->Dimension(), lh);

    // needed for evaluation of compiled expressions
    DummyFE<ET_TRIG> dummyfe;
    ud.fel = &dummyfe;

    // Prepare data structures for blocks
    const auto nblocks = proxies.Size();
    FlatArray<FlatMatrix<double>> xk_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> xold_blocks(nblocks, lh);
    FlatArray<FlatTensor<3, double>> lin_blocks(nblocks, lh);

    for (auto i : Range(nblocks)) {
      const auto proxy = proxies[i];
      xk_blocks[i].Assign(ud.GetMemory(proxy));
      xold_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
      lin_blocks[i].AssignMemory(lh, mir.Size(), eq_dim, proxy->Dimension());
    }

    // Block-agnostic data structures
    FlatMatrix<> xk(mir.Size(), full_dim, lh);
    FlatMatrix<> xold(mir.Size(), full_dim, lh);
    FlatMatrix<> w(mir.Size(), full_dim, lh);
    FlatMatrix<> res_all(mir.Size(), eq_dim, lh);
    FlatMatrix<AutoDiff<1, double>> dval(mir.Size(), full_dim, lh);
    FlatArray<double> res_0_qp(mir.Size(), lh);
    res_0_qp = 0;

    // auxiliary data structures for solving
    FlatArray<int> p(numeric_dim, lh);
    FlatMatrix<> Q(eq_dim, numeric_dim, lh);
    FlatMatrix<double> lhs(eq_dim, numeric_dim, lh);
    FlatVector<double> x(numeric_dim, lh);


    // Evaluate starting point
    if (startingpoints.Size() == proxies.Size())
      {
        for (int i : Range(startingpoints))
          startingpoints[i]->Evaluate(mir, xk_blocks[i]);

        // merge blocks
        size_t offset = 0;
        for (auto xkb : xk_blocks)
          {
            auto next = offset + xkb.Width();
            xk.Cols(offset, next) = xkb;
            offset = next;
          }
      }
    else
      {
        // This is the only case when we should end up here (checked in
        // constructor)
        assert(startingpoints.Size() == 1);

        startingpoints[0]->Evaluate(mir, xk);
        distribute_vec_to_blocks(xk, xk_blocks);
      }

//    cout << "starting value = " << xk << endl;
//    cout << "blocks:" << endl;
//    for (int i : Range(proxies))
//      cout << i << " -> " << ud.GetMemory(proxies[i]) << endl;
//    cout << endl;

    // The actual algorithm
    //    cout <<
    //    "\n\n------------------------START------------------------------"
    //            "\n";

    {
      RegionTimer regtr1(t1);
      calc_residuals(res_all, mir);
    }

    for (auto qi : Range(mir))
      res_0_qp[qi] = LInfNorm(res_all.Row(qi));

//    cout << "[0] res: " << res_all << endl;
//    cout << "[0] res_qp: " << res_0_qp << endl;

    bool success = all_converged_qp(res_all, tol, res_0_qp, rtol);
    for ([[maybe_unused]] int step : Range(maxiter))
      {
        if (success)
          break;

        {
          RegionTimer regtr2(t2);
          calc_linearizations(lin_blocks, dval, ud, mir);
          //      cout << "lin: " << lin_blocks << endl;
        }

        {
          // compute increments
          RegionTimer r(t3);
          for (size_t qi : Range(mir))
            {

              // NOTE: when to skip something because of convergence?
              //  -> the current approach assumes that evaluation of "expression"
              //  for all qpoints at once is beneficial!

              auto rhs = res_all.Row(qi);

              // skip if newton for this qp has already converged
              if (converged(rhs, tol, res_0_qp[qi], rtol))
                {
                  w.Row(qi) = 0;
                  continue;
                }

              // handle VS embedding
              form_lhs_from_lin_blocks(lhs, lin_blocks, qi);

              if (lhs.Height() == lhs.Width())
                {
                  p = 0;
                  CalcLU (lhs, p);
                  SolveFromLU (lhs, p, SliceMatrix<double, ColMajor>(rhs.Size(), 1, rhs.Size(), rhs.Data()));
                  expand_increments(rhs, w.Row(qi));
                }
              else
                {
                  assert(lhs.Height() > lhs.Width());
                  QRFactorization(lhs, Q);
                  x = Vector<double>(Trans(Q) * rhs);
                  TriangularSolve<UpperRight, NonNormalized>(lhs.Rows(0, x.Size()), x);
                  expand_increments(x, w.Row(qi));
                }
            }
        }

        xk -= w;
//        cout << "w: " << xk << endl;
//        cout << "xk: " << xk << endl;
        distribute_vec_to_blocks(xk, xk_blocks);

        {
          RegionTimer regtr1(t1);
          calc_residuals(res_all, mir);
//          cout << "res: " << res_all << endl;
        }
        success = all_converged_qp(res_all, tol, res_0_qp, rtol);
      }

//     cout << "xk (final): " << xk << endl;
    if (!success) {
      cout << IM(4) << "The NewtonCF did not converge to tolerance on element " << trafo.GetElementNr() << endl;

      for (auto qi : Range(res_all.Height())) {
        if (!converged(res_all.Row(qi), tol, res_0_qp[qi], rtol)) {
          cout << IM(5) << "Quadrature point index " << qi << ", ||res||_inf=" << LInfNorm(res_all.Row(qi));
          cout << IM(5) << ", ||res_0||_inf=" << res_0_qp[qi] << endl;      
        }
      }

      if (!allow_fail)
        xk = numeric_limits<double>::quiet_NaN();
    }

    // cout << "result = " << xk << endl;
    values.AddSize(mir.Size(), full_dim) = xk;
    //    cout <<
    //    "\n--------------------- NewtonCF done ---------------------------\n";
  }

private:
    template <typename src_t, typename dest_t>
    void expand_increments(const src_t src, dest_t dest) const
    {
      if (src.Size() == dest.Size())
        {
          dest = src;
          return;
        }

      // apply VS embedding to increment
      size_t cola = 0;
      size_t colb = 0;
      size_t colxa = 0;
      size_t colxb = 0;
      for (auto block : Range(proxies))
        {
          colb += proxies[block]->Dimension();
          if (const auto vsemb = get_vs_embedding(proxies[block]); vsemb)
            {
              colxb += vsemb.value().Width();
              dest.Range(cola, colb) = vsemb.value() * src.Range(colxa, colxb);
//              cout << dest.Range(cola, colb) << ", " << vsemb.value() * src.Range(colxa, colxb) << endl;
            }
          else
            {
              colxb += proxies[block]->Dimension();
              dest.Range(cola, colb) = src.Range(colxa, colxb);
//              cout << dest.Range(cola, colb) << ", " << src.Range(colxa, colxb) << endl;
            }
          cola = colb;
          colxa = colxb;
        }
    };


    template <typename src_t, typename dest_t>
    void distribute_vec_to_blocks(const src_t src,  dest_t dest) const {
        size_t offset = 0;
        for (auto destblock : dest)
          {
            auto next = offset + destblock.Width();
            destblock = src.Cols(offset, next);
            offset = next;
          }
    };

    template<typename res_all_t>
    inline void calc_residuals(res_all_t res_all, const BaseMappedIntegrationRule& mir) const
    {
      expression->Evaluate(mir, res_all);
    };

    template<typename lin_blocks_t, typename dval_t>
    void calc_linearizations(
            lin_blocks_t lin_blocks,
            dval_t dval,
            ProxyUserData& ud,
            const BaseMappedIntegrationRule& mir
            ) const
    {
      for (auto block2 : Range(proxies))
        {
          auto proxy2 = proxies[block2];

          for (auto l : Range(proxy2->Dimension()))
            {
              // This means, the derivative is taken wrt proxy 2
              ud.trialfunction = proxy2;
              ud.trial_comp = l;
              expression->Evaluate(mir, dval);

              for (auto qi : Range(mir))
                for (auto k : Range(dval.Width()))
                  lin_blocks[block2](qi, k, l) = dval(qi, k).DValue(0);
            }
        }
    };

    template<typename lhs_t, typename lin_blocks_t>
    void form_lhs_from_lin_blocks(lhs_t lhs, const lin_blocks_t lin_blocks, size_t qi) const
    {
      size_t cola = 0;
      size_t colb = 0;
      for (auto block : Range(proxies))
        {
          if (const auto vsemb = get_vs_embedding(proxies[block]); vsemb)
            {
              colb += vsemb.value().Width();
              lhs.Cols(cola, colb) = lin_blocks[block](qi, STAR, STAR) * vsemb.value();
            }
          else
            {
              colb += proxies[block]->Dimension();
              lhs.Cols(cola, colb) = lin_blocks[block](qi, STAR, STAR);
            }
          cola = colb;
        }
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
  double tol{1e-6};
  double rtol{0.0};
  int maxiter{20};

  bool allow_fail{false};

public:
  MinimizationCF(shared_ptr<CoefficientFunction> aexpression,
                 shared_ptr<CoefficientFunction> astartingpoint,
                 std::optional<double> atol, std::optional<double> artol,
                 std::optional<int> amaxiter, std::optional<bool> aallow_fail)
      : MinimizationCF{aexpression,
                       Array<shared_ptr<CoefficientFunction>>{astartingpoint},
                       atol, artol, amaxiter, aallow_fail} {}

  MinimizationCF(shared_ptr<CoefficientFunction> aexpression,
                 const Array<shared_ptr<CoefficientFunction>> &astartingpoints,
                 std::optional<double> atol, std::optional<double> artol,
                 std::optional<int> amaxiter, std::optional<bool> aallow_fail)
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
      std::map<int, ProxyFunction *> sorted_proxies{};
      std::map<int, const DifferentialOperator *> diffops{};
      //      cout << "\n" << "sorted proxies " << sorted_proxies
      //           << "\n" << "diffops " << diffops << endl;

      const auto compound_space = [](const auto& proxy) {
          return static_cast<void*>(proxy->GetFESpace().get());
      };
      const auto space = compound_space(proxies[0]);
      if (!space)
        throw Exception(
                "MinimizationCF: More than one proxy has been found but not all "
                "belong to a CompoundFESpace.");

      for (const auto proxy : proxies) {
        if (space != compound_space(proxy))
          throw Exception(
                  "MinimizationCF: More than one proxy has been found but not all "
                  "belong to the same a CompoundFESpace.");

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
              std::find_if(cbegin(diffops), cend(diffops), [&](const auto& pair){return pair.second == evaluator;}) !=
                  cend(diffops))
            throw Exception(
                "MinimizationCF: A proxy evaluator (component) has been "
                "detected twice");
          diffops[evaluator->Component()] = evaluator;
          sorted_proxies[evaluator->Component()] = proxy;
        }
      }
      // Copy over...
      for (int i1 = 0, i2 = 0; i1 < sorted_proxies.size(); ++i1)
        if  (auto proxy = sorted_proxies[i1])
            proxies[i2++] = proxy;
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
                     startingpoints[i]->Dimensions())) {
            std::stringstream sstr{};
            sstr << "MinimizationCF: Dimensions of startingpoint "
                    "and proxy component " << i << " do not agree ("
                    << startingpoints[i]->Dimensions() << " != " << proxies[i]->Dimensions();
            throw Exception(sstr.str());
          }
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
    // static Timer t("MinimizationCF::Eval", NoTracing);
    // static Timer t1("MinimizationCF::Eval get Jac", NoTracing);
    // static Timer t2("MinimizationCF::Eval solve", NoTracing);
    // RegionTimer reg(t);
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

    // needed for evaluation of compiled expressions
    DummyFE<ET_TRIG> dummyfe;
    ud.fel = &dummyfe;

    // Prepare data structures for blocks
    const auto nblocks = proxies.Size();
    FlatArray<FlatMatrix<double>> xk_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> xold_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> diags_blocks(nblocks, lh);
    FlatArray<FlatMatrix<double>> res_blocks(nblocks, lh);
    FlatArray<double> res_0_blocks(nblocks, lh);
    FlatArray<double> res_0_qp(mir.Size(), lh);
    FlatArray<FlatTensor<3>> lin_blocks(nblocks * nblocks, lh);

    // These are only "independent" for blocks having "vsemb"; otherwise just
    // views
    FlatArray<FlatMatrix<double>> rhs_blocks(nblocks, lh);
    FlatArray<FlatTensor<3>> lhs_blocks(nblocks * nblocks, lh);

    res_0_blocks = 0;
    res_0_qp = 0;

    for (int i : Range(nblocks)) {
      const auto proxy = proxies[i];
      xk_blocks[i].Assign(ud.GetMemory(proxy));
      xold_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
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
    FlatMatrix<> res_mat(mir.Size(), numeric_dim, lh);
    FlatVector<> rhs(numeric_dim, lh);
    FlatArray<int> p(numeric_dim, lh);
    FlatMatrix<> lhs(numeric_dim, numeric_dim, lh);

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
      // NOTE: this computes only one triangle of the matrix, whereby the diagonal
      // blocks are fully set if there is a VS embedding!
      for (auto l1 : Range(nblocks)) {
        for (auto k1 : Range(l1, nblocks)) {

          auto proxy1 = proxies[k1];
          auto proxy2 = proxies[l1];
          auto &lin = lin_blocks[l1 * nblocks + k1];
//          auto &linT = lin_blocks[k1 * nblocks + l1];

          for (auto l : Range(proxy2->Dimension())) {

            for (auto k : Range(proxy1 == proxy2 ? l : 0, proxy1->Dimension())) {

              if (proxy1 == proxy2 && k == l)
                lin(STAR, k, k) = diags_blocks[k1].Col(k);
              else { // computed mixed second derivatives
                ud.testfunction = proxy2;
                ud.test_comp = l;
                ud.trialfunction = proxy1;
                ud.trial_comp = k;

                expression->Evaluate(mir, ddval);
                for (size_t qi : Range(mir))
                  dderiv(qi, 0) = ddval(qi, 0).DDValue(0);

                lin(STAR, l, k) = dderiv.Col(0);
                lin(STAR, l, k) -= diags_blocks[k1].Col(k);
                lin(STAR, l, k) -= diags_blocks[l1].Col(l);
                lin(STAR, l, k) *= 0.5;

                if (proxy1 == proxy2 && get_vs_embedding(proxy1))
                  lin(STAR, k, l) = lin(STAR, l, k);
              }
              // cout << "lin block (" << k1 << ", " << l1 << ")[*, *, "
              //      << l << "] = " << lin(STAR, STAR, l) << endl;
            }
          }
//          cout << "lin block (" << l1 << ", " << k1
//               << ") = " << lin << endl;
        }
      }
    };

    const auto compute_newton_step = [&]() -> void {
      for (size_t qi : Range(mir)) {

        // NOTE: when to skip something because of convergence?
        //  -> the current approach assumes that evaluation of "expression"
        //  for all qpoints at once is beneficial!

        int offset1 = 0;
        for (size_t block1 : Range(nblocks)) {
          const auto proxy1 = proxies[block1];
          const auto &rhsb = rhs_blocks[block1].Row(qi);

          for (int k : Range(rhsb.Size()))
            rhs[offset1 + k] = rhsb[k];

          int offset2 = offset1;
          for (size_t block2 : Range(block1, nblocks)) {
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
            // Fill q-point LHS with exploitation of symmetry. Unfortunately,
            // CalcInverse does not respect symmetry yet.
            for (auto k : Range(lhsb.Height())) {
              if (offset1 == offset2) {
                lhs(offset1 + k, offset1 + k) = lhsb(k, k);
                for (auto l : Range(k + 1, lhsb.Width())) {
                  lhs(offset1 + k, offset2 + l) = lhsb(k, l);
                  lhs(offset2 + l, offset1 + k) = lhsb(k, l);
                }
              } else {
                for (auto l : Range(lhsb.Width())) {
                  lhs(offset1 + k, offset2 + l) = lhsb(k, l);
                  lhs(offset2 + l, offset1 + k) = lhsb(k, l);
                }
              }
            }

            offset2 += lhsb.Width();
          }
          offset1 += rhsb.Size();
        }

        if (converged(rhs, tol, res_0_qp[qi], rtol)) {
          w.Row(qi) = 0;
          continue;
        }

//        cout << "RHS: " << rhs << endl;
//        cout << "LHS: " << lhs << endl;

        p = 0;
        CalcLU (lhs, p);
        SolveFromLU (lhs, p, SliceMatrix<double, ColMajor>(rhs.Size(), 1, rhs.Size(), rhs.Data()));
        auto &sol = rhs;

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

    const auto linesearch = [&]() -> bool {
      xold = xk;// linesearch
      double alpha = 1;
      double alpha_min = 1e-10;
      double newenergy = energy + 1;
      double energy_eps = 1e-10;

      auto proxy = proxies[0];
      ud.trialfunction = proxy;
      ud.trial_comp = 0;
      ud.testfunction = proxy;
      ud.test_comp = 0;

      // cout << "w = " << endl << w << endl;
      double energy_limit = energy + energy_eps * abs(energy);
      while (newenergy > energy_limit && alpha > alpha_min) {
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

      return !(newenergy > energy_limit);
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

//    for (auto block : Range(nblocks))
//      res_0_blocks[block] = LInfNorm(rhs_blocks[block].AsVector());


    for (auto qi : Range(mir))
      for (auto block : Range(nblocks))
        res_0_qp[qi] = max(LInfNorm(rhs_blocks[block].Row(qi)), res_0_qp[qi]);

    const auto all_converged = [&]() -> bool {
      for (size_t qi : Range(mir)) {
        int offset1 = 0;
        for (size_t block1 : Range(nblocks)) {
          const auto &rhsb = rhs_blocks[block1].Row(qi);
          res_mat.Row(qi).Range(rhsb.Size()) = rhsb;
        }
      }
      return all_converged_qp(res_mat, tol, res_0_qp, rtol);
    };

    bool success = all_converged();

    for ([[maybe_unused]] int step : Range(maxiter)) {
      if (success)
        break;

      calc_off_diagonals();
      compute_newton_step();
      if (!linesearch())
        break;
      calc_energy_rhs_and_diags();
      success = all_converged();
//      cout << "newton step " << step + 1 << endl;
    }

//    for (int block1 : Range(proxies))
//      cout << "RHS block " << to_string(block1) << ": " << rhs_blocks[block1] << endl;
//
//    cout << "MinimizationCF done" << "\n\n";

if (!success) {
  cout << IM(4) << "The MinimizationCF did not converge to tolerance on element " << trafo.GetElementNr() << endl;

  for (auto qi : Range(res_mat.Height())) {
    if (!converged(res_mat.Row(qi), tol, res_0_qp[qi], rtol)) {
      cout << IM(5) << "Quadrature point index " << qi << ", ||res||_inf=" << LInfNorm(res_mat.Row(qi));
      cout << IM(5) << ", ||res_0||_inf=" << res_0_qp[qi] << endl;
    }
  }

  if (!allow_fail)
    xk = numeric_limits<double>::quiet_NaN();
}

    // cout << "result = " << xk << endl;
    values.AddSize(mir.Size(), Dimension()) = xk;
  }
};

shared_ptr<CoefficientFunction>
CreateMinimizationCF(shared_ptr<CoefficientFunction> expression,
               shared_ptr<CoefficientFunction> startingpoint,
               std::optional<double> atol, std::optional<double> rtol,
               std::optional<int> maxiter,
               std::optional<bool> allow_fail) {
  return make_shared<MinimizationCF>(expression, startingpoint, atol, rtol, maxiter, allow_fail);
}

shared_ptr<CoefficientFunction>
CreateMinimizationCF(shared_ptr<CoefficientFunction> expression,
               const Array<shared_ptr<CoefficientFunction>> &startingpoints,
               std::optional<double> tol, std::optional<double> rtol,
               std::optional<int> maxiter,
               std::optional<bool> allow_fail) {
  return make_shared<MinimizationCF>(expression, startingpoints, tol, rtol, maxiter, allow_fail);
}

shared_ptr<CoefficientFunction>
CreateNewtonCF(shared_ptr<CoefficientFunction> expression,
               shared_ptr<CoefficientFunction> startingpoint,
               std::optional<double> atol, std::optional<double> rtol,
               std::optional<int> maxiter,
               std::optional<bool> allow_fail) {
  return make_shared<NewtonCF>(expression, startingpoint, atol, rtol, maxiter, allow_fail);
}

shared_ptr<CoefficientFunction>
CreateNewtonCF(shared_ptr<CoefficientFunction> expression,
               const Array<shared_ptr<CoefficientFunction>> &startingpoints,
               std::optional<double> tol, std::optional<double> rtol,
               std::optional<int> maxiter,
               std::optional<bool> allow_fail) {
  return make_shared<NewtonCF>(expression, startingpoints, tol, rtol, maxiter, allow_fail);
}

} // namespace ngfem
