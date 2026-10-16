#ifndef FILE_MPTOOLS_DEVICE_HPP
#define FILE_MPTOOLS_DEVICE_HPP

/*********************************************************************/
/* File:   mptools_device.hpp                                        */
/* Purpose: the FMM passes on a device of the common GPU layer       */
/*********************************************************************/

/*
  Included at the end of mptools.hpp: with FMM_Parameters::device, CalcMP runs M2M, M2L and L2L on device pools
  (mptools_device.cpp), for an FMM_Operator also the P2M of its device sources and the L2P, M2P and P2P of its device
  targets. A pool holds one tree level; slot k*N+c is component c of node k, so the kernels stay scalar.
*/

#include <functional>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>

namespace ngsbem
{
  // the values must match the variant tests in shiftz_gpu_kernel.h
  enum class FMMTranslation { SS = 0, SR = 1, RR = 2 };

  // the device pools of one FMM tree
  class NGS_DLL_HEADER FMMDeviceTree
  {
  public:
    struct Impl;
    // Cpu: the GPU layer's CPU reference device, Gpu: the registered GPU (throws if there is none); fp64: its pools, P2M and
    // L2P in double (Complex host expansions), else in float (Complex32)
    FMMDeviceTree (FMM_Device device, bool fp64);
    ~FMMDeviceTree ();

    // a pool of `slots` expansions of one order, filled by the plan's zeros or the host's uploads; returns its id
    int AddPool (int order, size_t slots);
    // fill(dst, first, n) produces the entries [first, first+n) of the pool into dst (fill gets pool entry indices), in
    // the tree's real type
    void Upload (int pool, size_t first, size_t n, const std::function<void(Complex32*, size_t, size_t)> & fill);
    void Upload (int pool, size_t first, size_t n, const std::function<void(Complex*, size_t, size_t)> & fill);
    // the entries [0, n) of the pool: consume(src, first, cnt) gets [first, first+cnt)
    void Download (int pool, size_t n, const std::function<void(const Complex32*, size_t, size_t)> & consume);
    void Download (int pool, size_t n, const std::function<void(const Complex*, size_t, size_t)> & consume);

    // host map, filled by the planning code: expansion (&node.mp) -> (pool, node index in the pool)
    std::unordered_map<const void*, std::pair<int,int>> slots;

    std::unique_ptr<Impl> impl;
  };

  // the P2M of singular leaves on the device (p2m_gpu_kernel.h): the FMMDeviceSources kinds into the tree's expansions
  class NGS_DLL_HEADER FMMDeviceP2M
  {
  public:
    struct Impl;
    static constexpr int max_order = 31;     // the kernel's group memory
    // kind, comps as FMMDeviceSources (the values come per matvec by Upload); source j is point order[j], geo[5j..5j+4]
    // its PolarGeometry relative to its leaf center, dir[3j..3j+2] its dipole direction and xyz[3j..3j+2] its position
    FMMDeviceP2M (const FMMDeviceTree & tree, int kind, int comps, std::variant<double,std::complex<double>> kappa,
                  FlatArray<int> order, FlatArray<double> geo, FlatArray<double> dir, FlatArray<double> xyz);
    ~FMMDeviceP2M ();
    // leaf k of the pool: the sources [first[k], first[k]+count[k]) into slot[k]; scale = the leaves' Scale()
    void AddLevel (const FMMDeviceTree & tree, int pool, double scale, FlatArray<int> first, FlatArray<int> count,
                   FlatArray<int> slot);
    // per matvec: the source values (in point order) to the device
    void Upload (const FMMDeviceSources & sources);
    // the sources [first, first+count) of the leaf in slot of pool (AddLevel)
    std::pair<int,int> Sources (int pool, int slot) const;
    std::unique_ptr<Impl> impl;
  };

  // the L2P and M2P at the targets on the device (l2p_gpu_kernel.h): the local expansions of the target leaves and
  // the far singular expansions up to max_order, the values and gradients per target and component
  class NGS_DLL_HEADER FMMDeviceL2P
  {
  public:
    struct Impl;
    static constexpr int max_order = FMMDeviceP2M::max_order;     // the kernels share the P2M tables
    // order, geo, xyz of the targets as FMMDeviceP2M; mode as FMMDeviceTargets
    FMMDeviceL2P (const FMMDeviceTree & tree, int comps, int mode, std::variant<double,std::complex<double>> kappa,
                  FlatArray<int> order, FlatArray<double> geo, FlatArray<double> xyz);
    ~FMMDeviceL2P ();
    // leaf k of the regular pool: the targets [first[k], first[k]+count[k]) of slot[k]; scale = 1/Scale() of the leaves
    void AddLevel (const FMMDeviceTree & tree, int pool, double scale, FlatArray<int> first, FlatArray<int> count,
                   FlatArray<int> slot);
    // from the singular pool: entry e adds the pairs [first[e], first[e]+count[e]) to target[e]; pair q is slot[q] with
    // geo[5q..5q+4] of the target relative to the node center; scale = Scale() of the nodes (the gradient: EvaluateMPDeriv)
    void AddM2P (const FMMDeviceTree & tree, int pool, double scale, FlatArray<int> target, FlatArray<int> first,
                 FlatArray<int> count, FlatArray<int> slot, FlatArray<double> geo);
    // the near field from the sources of p2m: entry e adds to target[e] the source ranges [list[e], list[e]+lcount[e])
    // of ranges (first source, count pairs); no two entries of a call share a target; the gradient only from charges,
    // as EvaluateDeriv
    void AddP2P (const FMMDeviceP2M & p2m, FlatArray<int> target, FlatArray<int> list, FlatArray<int> lcount,
                 FlatArray<int> ranges);
    // per matvec: the values into targets.values (point order)
    void Download (const FMMDeviceTargets & targets);
    std::unique_ptr<Impl> impl;
  };

  // recorded device work (FMMDevicePlan::Record), replayed per matvec by the plan cache
  class NGS_DLL_HEADER FMMDeviceProgram
  {
  public:
    struct Impl;
    FMMDeviceProgram ();
    ~FMMDeviceProgram ();
    // replays the recorded launches on their queue and waits for it
    void Replay ();
    std::unique_ptr<Impl> impl;
  };

  // the translation groups of one CalcMP, run on one in-order queue in the order they were added
  class NGS_DLL_HEADER FMMDevicePlan
  {
  public:
    struct Impl;
    FMMDevicePlan ();
    ~FMMDevicePlan ();
    // the zeros and groups added until the next BeginPass form one pass; a plan starts with one
    void BeginPass ();
    // zero the entries [first, first+count) of a pool before the groups of the current pass
    void AddZero (const FMMDeviceTree & tree, int pool, size_t first, size_t count);
    // one CPU group; lane b: slot tgt_idx[b] of (tgt, tgt_pool) += Transform of slot src_idx[b] of
    // (src, src_pool) with the group's len and theta (after the z-flip) and the lane's phi and flip
    void AddGroup (FMMTranslation variant, std::variant<double,std::complex<double>> kappa,
                   double source_rtyp, double target_rtyp, int os, int ot, double len, double theta,
                   const FMMDeviceTree & src, int src_pool, const FMMDeviceTree & tgt, int tgt_pool,
                   FlatArray<double> phi, FlatArray<int> flip, FlatArray<int> src_idx, FlatArray<int> tgt_idx);
    // launched before the first pass
    void SetP2M (std::shared_ptr<FMMDeviceP2M> p2m);
    // launched after the last pass
    void SetL2P (std::shared_ptr<FMMDeviceL2P> l2p);
    // fuses each pass's groups into units, prepares them, reserves one workspace, enqueues all in order and waits
    void Run ();
    // as Run, but records the launches without running them; one workspace with share (a program replayed
    // alongside), sized for both: share is re-recorded if its own is not that size
    std::shared_ptr<FMMDeviceProgram> Record (FMMDeviceProgram * share = nullptr);
  private:
    std::unique_ptr<Impl> impl;
  };

  // the device plan of one FMM_Operator direction: programs A (P2M, M2M) and B (M2L, L2L, L2P, M2P, P2P) on one
  // workspace, replayed while the host trees' fingerprints match; after its recording a covered plan runs without
  // host trees (FMMDeviceMatvec), relying on the operator's points, normals, kernel and parameters never changing
  class FMMDeviceCache
  {
  public:
    struct Side
    {
      std::vector<double> fingerprint;
      std::shared_ptr<FMMDeviceTree> tree;
      std::shared_ptr<FMMDeviceProgram> program;
      std::vector<int> level_pools;          // pool id per tree level, -1 for a level without expansions
      std::shared_ptr<FMMDeviceP2M> p2m;     // singular: the device P2M at the start of program A
      std::shared_ptr<FMMDeviceL2P> l2p;     // regular: the device L2P + M2P at the end of program B
      bool covered = false;                  // singular: every leaf's P2M on the device; regular: all of every target's
                                             // values and gradients (L2P, M2P and P2P) on the device
    };
    Side singular, regular;
    const void * singular_owner = nullptr;   // the SingularMLExpansion whose values singular.tree holds
    std::vector<char> host_m2p_levels;       // with regular.l2p: the singular levels whose expansions the host evaluates
    std::vector<char> singular_skipped;      // the singular interior levels this matvec did not download
    bool fast = false;                       // the recorded plan is covered
  };


  namespace fmm_device_detail
  {
    template <typename entry_type>
    constexpr bool IsDoubleEntry = std::is_same_v<typename SphericalHarmonics<entry_type>::real_type, double>;

    // the expansions of order >= 0 on one tree level
    template <typename NODE>
    struct Level
    {
      std::vector<NODE*> nodes;   // the nfirst selected ones first, then the others, in traversal order
      size_t nfirst = 0;
      int pool = -1;
      int order = -1;
    };

    template <typename NODE, typename FIRST>
    std::vector<Level<NODE>> CollectLevels (NODE & root, FIRST first)
    {
      std::vector<std::vector<NODE*>> firsts, others;
      root.TraverseTree ([&] (NODE & node)
      {
        if (node.mp.Order() < 0) return;
        if (firsts.size() <= size_t(node.level))
          {
            firsts.resize (node.level+1);
            others.resize (node.level+1);
          }
        (first(node) ? firsts : others)[node.level].push_back (&node);
      });
      std::vector<Level<NODE>> levels (firsts.size());
      for (size_t l = 0; l < levels.size(); l++)
        {
          levels[l].nodes = firsts[l];
          levels[l].nodes.insert (levels[l].nodes.end(), others[l].begin(), others[l].end());
          levels[l].nfirst = firsts[l].size();
        }
      return levels;
    }

    // one pool per nonempty level, and the host map; returns the pool id per level
    template <int N, typename NODE>
    std::vector<int> MakePools (FMMDeviceTree & tree, std::vector<Level<NODE>> & levels)
    {
      std::vector<int> pools (levels.size(), -1);
      for (size_t l = 0; l < levels.size(); l++)
        {
          auto & lev = levels[l];
          if (lev.nodes.empty()) continue;
          int order = lev.nodes[0]->mp.Order();
          double rtyp = lev.nodes[0]->mp.RTyp();
          for (auto * node : lev.nodes)
            if (node->mp.Order() != order || node->mp.RTyp() != rtyp)
              throw Exception ("fmm device: the expansions on tree level " + ToString(l) + " differ in order or radius "
                               + "(order " + ToString(node->mp.Order()) + " vs " + ToString(order) + ", r "
                               + ToString(node->mp.RTyp()) + " vs " + ToString(rtyp)
                               + "); the device pools need one order and radius per level");
          lev.order = order;
          pools[l] = lev.pool = tree.AddPool (order, lev.nodes.size()*N);
          for (size_t k = 0; k < lev.nodes.size(); k++)
            tree.slots[&lev.nodes[k]->mp] = { lev.pool, int(k) };
        }
      return pools;
    }

    // a cached plan's pools for this (structurally identical) tree: the pool ids per level and the host map
    template <typename NODE>
    void AssignSlots (FMMDeviceTree & tree, std::vector<Level<NODE>> & levels, const std::vector<int> & level_pools)
    {
      if (level_pools.size() != levels.size())
        throw Exception ("fmm device cache: the tree levels do not match the cached plan");
      tree.slots.clear();
      for (size_t l = 0; l < levels.size(); l++)
        {
          auto & lev = levels[l];
          lev.pool = level_pools[l];
          if (lev.nodes.empty()) continue;
          lev.order = lev.nodes[0]->mp.Order();
          for (size_t k = 0; k < lev.nodes.size(); k++)
            tree.slots[&lev.nodes[k]->mp] = { lev.pool, int(k) };
        }
    }

    // what a plan depends on, per node in the deterministic TraverseTree order; the node sequence per level
    // also fixes the slot of every node
    template <typename NODE>
    std::vector<double> Fingerprint (NODE & root, int N, std::complex<double> kappa)
    {
      std::vector<double> fp { double(N), kappa.real(), kappa.imag() };
      root.TraverseTree ([&] (NODE & node)
      {
        fp.insert (fp.end(), { double(node.level), node.center(0), node.center(1), node.center(2), node.r,
                               double(node.mp.Order()), node.mp.RTyp(), node.childs[0] ? 1.0 : 0.0 });
      });
      return fp;
    }

    // the complex type of a node's expansion coefficients: Complex32 or Complex
    template <typename NODE>
    using HostComplex = typename decltype(std::declval<NODE&>().mp)::complex_type;

    template <typename NODE>
    HostComplex<NODE> * HostCoefs (const NODE & node)
    {
      return reinterpret_cast<HostComplex<NODE>*> (node.mp.SH().Coefs().Data());
    }

    // the pool entries [first, first+n) as runs within one slot: f(node k, component c, first ii, count, offset)
    template <int N, typename F>
    void ForSlotRuns (size_t rows, size_t first, size_t n, F f)
    {
      for (size_t e = first, end = first+n; e < end; )
        {
          size_t s = e / rows, ii0 = e % rows, cnt = std::min (rows-ii0, end-e);
          f (s / N, s % N, ii0, cnt, e-first);
          e += cnt;
        }
    }

    // the nodes after the first nfirst of every pool (the singular leaves after P2M) from their host expansions,
    // except the levels whose first such node skip accepts
    template <int N, typename NODE, typename SKIP>
    void UploadRest (FMMDeviceTree & tree, const std::vector<Level<NODE>> & levels, SKIP skip)
    {
      for (auto & lev : levels)
        {
          if (lev.pool < 0 || lev.nfirst >= lev.nodes.size() || skip (*lev.nodes[lev.nfirst])) continue;
          size_t rows = sqr (size_t(lev.order)+1);
          size_t begin = lev.nfirst*N*rows, n = (lev.nodes.size()-lev.nfirst)*N*rows;
          tree.Upload (lev.pool, begin, n, [&] (HostComplex<NODE> * dst, size_t first, size_t cnt)
          {
            ForSlotRuns<N> (rows, first, cnt, [&] (size_t k, size_t c, size_t ii0, size_t run, size_t off)
            {
              const auto * h = HostCoefs (*lev.nodes[k]);
              for (size_t i = 0; i < run; i++)
                dst[off+i] = h[(ii0+i)*N+c];
            });
          });
        }
    }

    // the first nfirst nodes of every pool into their host expansions, except the levels l with skip(l)
    template <int N, typename NODE, typename SKIP>
    void DownloadFirst (FMMDeviceTree & tree, const std::vector<Level<NODE>> & levels, SKIP skip)
    {
      for (size_t l = 0; l < levels.size(); l++)
        {
          auto & lev = levels[l];
          if (lev.pool < 0 || lev.nfirst == 0 || skip (l)) continue;
          size_t rows = sqr (size_t(lev.order)+1);
          tree.Download (lev.pool, lev.nfirst*N*rows, [&] (const HostComplex<NODE> * src, size_t first, size_t n)
          {
            ForSlotRuns<N> (rows, first, n, [&] (size_t k, size_t c, size_t ii0, size_t cnt, size_t off)
            {
              auto * h = HostCoefs (*lev.nodes[k]);
              for (size_t i = 0; i < cnt; i++)
                h[(ii0+i)*N+c] = src[off+i];
            });
          });
        }
    }

    inline std::pair<int,int> Slot (const FMMDeviceTree & tree, const void * mp)
    {
      auto it = tree.slots.find (mp);
      if (it == tree.slots.end())
        throw Exception ("fmm device: a recorded expansion has no device slot");
      return it->second;
    }

    // the factor of r in the device geometry: kappa, or |kappa| for a complex kappa (the kernels multiply by kappa/|kappa|)
    inline double GeometryKappa (double kappa) { return kappa; }
    inline double GeometryKappa (std::complex<double> kappa) { return std::abs(kappa); }

    // cos theta, sin theta, cos phi, sin phi and kappa r of x, as SphericalHarmonics::Polar
    inline void PolarGeometry (Vec<3> x, double kappa, double * g)
    {
      double r = L2Norm(x), rxy = sqrt(sqr(x(0))+sqr(x(1)));
      double v[5] = { x(2) > 0 ? 1.0 : -1.0, 0, 1, 0, kappa*r };
      if (rxy > 0)
        { v[0] = x(2)/r; v[1] = rxy/r; v[2] = x(0)/rxy; v[3] = x(1)/rxy; }
      for (int k = 0; k < 5; k++) g[k] = v[k];
    }

    // the points i with leaf_of[i] >= 0 sorted by (leaf, point): order[j] = i, leaf l has [first[l], first[l+1]);
    // throws msg unless matches(l, count of leaf l)
    template <typename MATCHES>
    std::pair<Array<int>, Array<int>> SortByLeaf (FlatArray<int> leaf_of, size_t nleaves, MATCHES matches, const char * msg)
    {
      Array<int> first(nleaves+1), next(nleaves);
      first = 0;
      for (int l : leaf_of)
        if (l >= 0) first[l+1]++;
      for (size_t l = 0; l < nleaves; l++)
        {
          if (!matches (l, size_t(first[l+1])))
            throw Exception (msg);
          first[l+1] += first[l];
          next[l] = first[l];
        }
      Array<int> order(first[nleaves]);
      for (size_t i = 0; i < leaf_of.Size(); i++)
        if (leaf_of[i] >= 0) order[next[leaf_of[i]]++] = int(i);
      return { std::move(order), std::move(first) };
    }

    // one CPU group as lanes (record b, component c), in record order
    template <int N, typename REC, typename SRC, typename TGT>
    void AddGroup (FMMDevicePlan & plan, FMMTranslation variant, FlatArray<REC*> group, double len, double theta,
                   const FMMDeviceTree & src, const FMMDeviceTree & tgt, SRC source_of, TGT target_of)
    {
      if (group.Size() == 0) return;
      auto * s0 = source_of (*group[0]);
      auto * t0 = target_of (*group[0]);
      int spool = Slot (src, s0).first;
      int tpool = Slot (tgt, t0).first;
      size_t lanes = group.Size()*N;
      Array<double> phi(lanes);
      Array<int> flip(lanes), src_idx(lanes), tgt_idx(lanes);
      for (size_t b = 0; b < group.Size(); b++)
        {
          auto [sp, sk] = Slot (src, source_of (*group[b]));
          auto [tp, tk] = Slot (tgt, target_of (*group[b]));
          if (sp != spool || tp != tpool)
            throw Exception ("fmm device: a translation group spans several tree levels");
          for (int c = 0; c < N; c++)
            {
              size_t lane = b*N+c;
              phi[lane] = group[b]->phi;
              flip[lane] = group[b]->flipz ? 1 : 0;
              src_idx[lane] = sk*N+c;
              tgt_idx[lane] = tk*N+c;
            }
        }
      plan.AddGroup (variant, s0->Kappa(), s0->RTyp(), t0->RTyp(), s0->Order(), t0->Order(), len, theta,
                     src, spool, tgt, tpool, phi, flip, src_idx, tgt_idx);
    }
  }


  // M2M: called from SingularMLExpansion::CalcMP after P2M and the CPU grouping
  template <typename entry_type, typename T_Kappa>
  struct FMMDeviceGlue<SingularMLExpansion<entry_type,T_Kappa>>
  {
    using ML = SingularMLExpansion<entry_type,T_Kappa>;
    using Node = typename ML::Node;
    using Rec = typename ML::RecordingSS;

    // zero the interior slots, then the CPU's S->S groups in its order: sorted by len (sqrt(3) r_child), they come deepest
    // level first, so children complete before parents (checked: the pool ids grow with the level, MakePools)
    template <int N>
    static void AddGroups (FMMDevicePlan & plan, const FMMDeviceTree & tree,
                           const std::vector<fmm_device_detail::Level<Node>> & levels,
                           FlatArray<Array<Rec*>> groups, FlatArray<double> lens, FlatArray<double> thetas)
    {
      using namespace fmm_device_detail;
      for (auto & lev : levels)
        if (lev.pool >= 0 && lev.nfirst > 0)
          plan.AddZero (tree, lev.pool, 0, lev.nfirst*N*sqr(size_t(lev.order)+1));
      int prev = std::numeric_limits<int>::max();
      for (size_t g = 0; g < groups.Size(); g++)
        {
          int pool = groups[g].Size() ? Slot (tree, groups[g][0]->mp_source).first : prev;
          if (pool > prev)
            throw Exception ("fmm device: the M2M groups do not come deepest level first");
          prev = pool;
          AddGroup<N> (plan, FMMTranslation::SS, FlatArray<Rec*>(groups[g]), lens[g], thetas[g], tree, tree,
                       [] (Rec & r) { return r.mp_source; }, [] (Rec & r) { return r.mp_target; });
        }
    }

    // the leaves whose P2M runs on the device, so the host P2M skips them: a source type with a device P2M
    // (FMMDeviceSources), and the order within the kernel's bound (the same for a whole level)
    static bool P2MOnDevice (const ML & ml, const Node & node)
    {
      auto * src = ml.fmm_params.device_sources;
      return src && src->kind && ml.fmm_params.device_cache && src->comps == VecLength<entry_type>
        && !node.childs[0] && node.mp.Order() <= FMMDeviceP2M::max_order && node.currents.Size() == 0;
    }

    static bool P2MLevel (const ML & ml, const fmm_device_detail::Level<Node> & lev)
    {
      return lev.pool >= 0 && lev.nfirst < lev.nodes.size() && P2MOnDevice (ml, *lev.nodes[lev.nfirst]);
    }

    // the device P2M of the leaves P2MOnDevice accepts, planned in double: each point's leaf by the descent AddSource
    // takes, the sources sorted by (leaf, point) with their geometry; nullptr if there are none
    template <int N>
    static std::shared_ptr<FMMDeviceP2M> PlanP2M (const ML & ml, const FMMDeviceTree & tree,
                                                  const std::vector<fmm_device_detail::Level<Node>> & levels)
    {
      using namespace fmm_device_detail;
      auto & src = *ml.fmm_params.device_sources;
      std::unordered_map<const Node*, int> index;
      std::vector<const Node*> leaves;
      for (auto & lev : levels)
        if (P2MLevel (ml, lev))
          for (size_t k = lev.nfirst; k < lev.nodes.size(); k++)
            {
              index[lev.nodes[k]] = int(leaves.size());
              leaves.push_back (lev.nodes[k]);
            }
      if (leaves.empty()) return nullptr;

      size_t npts = src.pts.Size();
      Array<int> leaf_of(npts);
      for (size_t i = 0; i < npts; i++)
        {
          const Node * node = &ml.root;
          while (node->childs[0])
            node = node->childs[node->GetChildNum(src.pts[i])].get();
          auto it = index.find (node);
          leaf_of[i] = it == index.end() ? -1 : it->second;
        }
      auto [order, first] = SortByLeaf (leaf_of, leaves.size(), [&] (size_t l, size_t n)
        { return n*src.entries == leaves[l]->charges.Size() + leaves[l]->dipoles.Size() + leaves[l]->chargedipoles.Size(); },
        "fmm device P2M: the source points do not match the tree's leaves");

      size_t nsrc = order.Size();
      Array<double> geo(5*nsrc), dir(3*nsrc), xyz(3*nsrc);
      double gkappa = GeometryKappa (ml.Kappa());
      for (size_t j = 0; j < nsrc; j++)
        {
          int i = order[j];
          PolarGeometry (src.pts[i] - leaves[leaf_of[i]]->center, gkappa, &geo[5*j]);
          for (int k = 0; k < 3; k++)
            {
              xyz[3*j+k] = src.pts[i](k);
              dir[3*j+k] = -src.nv[i](k);
            }
        }

      auto p2m = make_shared<FMMDeviceP2M> (tree, src.kind, N, ml.Kappa(), order, geo, dir, xyz);
      for (auto & lev : levels)
        if (P2MLevel (ml, lev))
          {
            Array<int> lfirst, lcount, lslot;
            for (size_t k = lev.nfirst; k < lev.nodes.size(); k++)
              {
                int l = index[lev.nodes[k]];
                lfirst.Append (first[l]);
                lcount.Append (first[l+1]-first[l]);
                lslot.Append (int(k));
              }
            p2m->AddLevel (tree, lev.pool, lev.nodes[lev.nfirst]->mp.Scale(), lfirst, lcount, lslot);
          }
      return p2m;
    }

    static void Translate (ML & ml, FlatArray<Array<Rec*>> groups, FlatArray<double> lens, FlatArray<double> thetas)
    {
      using namespace fmm_device_detail;
      constexpr int N = VecLength<entry_type>;
      static_assert (sizeof(entry_type) == N*sizeof(HostComplex<Node>));

      // interior nodes first: M2P (Node::EvaluateMP, EvaluateMPDeriv) reads their expansions afterwards
      auto levels = CollectLevels (ml.root, [] (const Node & node) { return bool(node.childs[0]); });

      if (auto cache = ml.fmm_params.device_cache)
        {
          // FMM_Operator: the plan is recorded once (program A) and replayed per matvec
          auto & side = cache->singular;
          auto fp = Fingerprint (ml.root, N, std::complex<double>(ml.Kappa()));
          if (!side.program || !side.tree || side.fingerprint != fp)
            {
              side = {};              // free the old plans first: program B and its L2P read these pools
              cache->regular = {};
              cache->fast = false;
              auto tree = make_shared<FMMDeviceTree> (ml.fmm_params.device, IsDoubleEntry<entry_type>);
              side.level_pools = MakePools<N> (*tree, levels);
              FMMDevicePlan plan;
              if (ml.fmm_params.device_sources && (side.p2m = PlanP2M<N> (ml, *tree, levels)))
                {
                  plan.SetP2M (side.p2m);
                  side.covered = std::all_of (levels.begin(), levels.end(), [&] (auto & lev)
                    { return lev.nfirst == lev.nodes.size() || P2MLevel (ml, lev); });
                }
              AddGroups<N> (plan, *tree, levels, groups, lens, thetas);
              side.program = plan.Record();
              side.tree = tree;
              side.fingerprint = std::move(fp);
            }
          else
            AssignSlots (*side.tree, levels, side.level_pools);   // the M2L plan maps this tree's expansions
          // per matvec: the host P2M leaves or the device P2M's source values up, replay A (P2M, zero the interior,
          // M2M), the interior down
          UploadRest<N> (*side.tree, levels, [&] (const Node & node) { return P2MOnDevice (ml, node); });
          if (side.p2m)
            side.p2m->Upload (*ml.fmm_params.device_sources);
          side.program->Replay();
          // the interior for the host M2P: once the regular side has a device L2P, only the levels it leaves to the host
          // (its plan stays valid: the regular tree of an FMM_Operator changes only with a singular re-plan)
          bool known = cache->regular.program && cache->regular.l2p;
          cache->singular_skipped.assign (levels.size(), 0);
          DownloadFirst<N> (*side.tree, levels, [&] (size_t l)
                            { return cache->singular_skipped[l] =
                                known && (l >= cache->host_m2p_levels.size() || !cache->host_m2p_levels[l]); });
          cache->singular_owner = &ml;
          return;
        }

      auto tree = make_shared<FMMDeviceTree> (ml.fmm_params.device, IsDoubleEntry<entry_type>);
      MakePools<N> (*tree, levels);
      FMMDevicePlan plan;
      AddGroups<N> (plan, *tree, levels, groups, lens, thetas);
      // the P2M results of the leaves; the plan zeros the interior nodes
      UploadRest<N> (*tree, levels, [] (const Node &) { return false; });
      plan.Run();
      DownloadFirst<N> (*tree, levels, [] (size_t) { return false; });
      ml.device_tree = tree;   // the pools are the M2L sources
    }
  };


  // M2L and L2L: called from RegularMLExpansion::CalcMP after the S->R recording and the CPU grouping;
  // replaces ProcessBatchSR and LocalizeExpansionBatched
  template <typename entry_type, typename T_Kappa>
  struct FMMDeviceGlue<RegularMLExpansion<entry_type,T_Kappa>>
  {
    using ML = RegularMLExpansion<entry_type,T_Kappa>;
    using Node = typename ML::Node;
    using RecSR = typename ML::RecordingSR;
    using RecRR = typename ML::RecordingRR;
    using SML = SingularMLExpansion<entry_type,T_Kappa>;

    // Node::RecordLocalizeLevel with a structural nonzero test: a node is nonzero if it received an
    // S->R record or its parent was nonzero (the CPU tests L2Norm(mp) > 0, which differs only by zeros)
    static bool RecordLevel (Node & node, int target_level, bool allow_refine, Array<RecRR> & recording,
                             Array<Node*> & nodes_to_clear, std::unordered_set<const void*> & nonzero)
    {
      if (node.level != target_level)
        {
          bool found = false;
          if (node.childs[0])
            for (auto & child : node.childs)
              found |= RecordLevel (*child, target_level, allow_refine, recording, nodes_to_clear, nonzero);
          return found;
        }

      if (allow_refine && node.mp.Order() > 30 && !node.childs[0])
        node.CreateChilds(allow_refine);

      if (!node.childs[0])
        return true;

      if (node.mp.Order() >= 0 && nonzero.count(&node.mp))
        for (int nr = 0; nr < 8; nr++)
          if (node.childs[nr]->mp.Order() >= 0)
            {
              recording += RecRR(&node.mp, &node.childs[nr]->mp, node.childs[nr]->center-node.center);
              nonzero.insert (&node.childs[nr]->mp);
            }

      nodes_to_clear.Append(&node);
      return true;
    }

    // the R->R records of all levels (with allow_refine this creates nodes), grouped as ProcessLocalizeRecording:
    // by theta, one len per level
    struct LocalizePlan
    {
      std::vector<Array<RecRR>> recordings;
      std::vector<Array<Array<RecRR*>>> groups;
      std::vector<Array<double>> thetas;
      std::vector<double> len;
      Array<Node*> nodes_to_clear;
    };

    static void PlanLocalize (ML & ml, FlatArray<Array<RecSR*>> sr_groups, bool allow_refine, LocalizePlan & lp)
    {
      std::unordered_set<const void*> nonzero;
      for (auto & group : sr_groups)
        for (auto * rec : group)
          nonzero.insert (rec->mpR);
      for (int level = 0; ; level++)
        {
          Array<RecRR> recording;
          if (!RecordLevel (ml.root, level, allow_refine, recording, lp.nodes_to_clear, nonzero))
            break;
          lp.recordings.push_back (std::move(recording));
        }
      lp.groups.resize (lp.recordings.size());
      lp.thetas.resize (lp.recordings.size());
      lp.len.resize (lp.recordings.size());
      for (size_t l = 0; l < lp.recordings.size(); l++)
        {
          auto & recording = lp.recordings[l];
          if (recording.Size() == 0) continue;
          QuickSort (recording, [] (auto & a, auto & b) { return a.theta < b.theta; });
          lp.len[l] = recording[0].len;
          double current_theta = -1e100;
          Array<RecRR*> current_batch;
          for (auto & record : recording)
            {
              bool theta_changed = fabs(record.theta - current_theta) > 1e-8;
              if (theta_changed && current_batch.Size() > 0)
                {
                  lp.groups[l].Append (current_batch);
                  lp.thetas[l].Append (current_theta);
                  current_batch.SetSize(0);
                }
              current_theta = record.theta;
              current_batch.Append (&record);
            }
          if (current_batch.Size() > 0)
            {
              lp.groups[l].Append (current_batch);
              lp.thetas[l].Append (current_theta);
            }
        }
    }

    // zero all regular slots, pass M2L (all S->R groups), then pass L2L (level by level, top-down)
    template <int N>
    static void AddGroups (FMMDevicePlan & plan, const FMMDeviceTree & src, const FMMDeviceTree & tree,
                           const std::vector<fmm_device_detail::Level<Node>> & levels,
                           FlatArray<Array<RecSR*>> sr_groups, FlatArray<double> sr_lens, FlatArray<double> sr_thetas,
                           const LocalizePlan & lp)
    {
      using namespace fmm_device_detail;
      for (auto & lev : levels)
        if (lev.pool >= 0)
          plan.AddZero (tree, lev.pool, 0, lev.nodes.size()*N*sqr(size_t(lev.order)+1));
      for (size_t g = 0; g < sr_groups.Size(); g++)
        AddGroup<N> (plan, FMMTranslation::SR, FlatArray<RecSR*>(sr_groups[g]), sr_lens[g], sr_thetas[g], src, tree,
                     [] (RecSR & r) { return r.mpS; }, [] (RecSR & r) { return r.mpR; });
      plan.BeginPass();
      for (size_t l = 0; l < lp.groups.size(); l++)
        for (size_t g = 0; g < lp.groups[l].Size(); g++)
          AddGroup<N> (plan, FMMTranslation::RR, FlatArray<RecRR*>(lp.groups[l][g]), lp.len[l], lp.thetas[l][g], tree, tree,
                       [] (RecRR & r) { return r.mp_source; }, [] (RecRR & r) { return r.mp_target; });
    }

    // the cache of this call (FMM_Operator) if its singular side holds singmp's pools, else nullptr (the uncached path)
    static FMMDeviceCache * UsableCache (ML & ml)
    {
      auto cache = ml.fmm_params.device_cache.get();
      return cache && cache->singular_owner == ml.singmp.get() && cache->singular.tree ? cache : nullptr;
    }

    // the leaves whose L2P and M2P (up to FMMDeviceL2P::max_order) run on the device, so the host skips them: device
    // targets (FMMDeviceTargets), and the order within the kernel's bound (the same for a whole level); not gradients
    // from charge-dipoles, which the host's EvaluateDeriv rejects
    static bool L2POnDevice (const ML & ml, const Node & node)
    {
      auto * tg = ml.fmm_params.device_targets;
      auto * src = ml.fmm_params.device_sources;
      return tg && ml.fmm_params.device_cache && tg->comps == VecLength<entry_type>
        && !((tg->mode & 2) && src && src->kind == 3)
        && !node.childs[0] && node.mp.Order() >= 0 && node.mp.Order() <= FMMDeviceL2P::max_order;
    }

    static bool L2PLevel (const ML & ml, const fmm_device_detail::Level<Node> & lev)
    {
      return lev.pool >= 0 && lev.nfirst > 0 && L2POnDevice (ml, *lev.nodes[0]);
    }

    // the device L2P and M2P of the targets in the leaves L2POnDevice accepts, planned in double: each point's leaf by the
    // descent AddTarget takes, the targets sorted by (leaf, point) with their geometry, and the (target, singular node)
    // pairs of SingularMLExpansion::Node::EvaluateMP that the device takes; sets cache.host_m2p_levels to the singular
    // levels the host still evaluates; nullptr if no leaf is covered
    template <int N>
    static std::shared_ptr<FMMDeviceL2P> PlanL2P (const ML & ml, const FMMDeviceTree & tree,
                                                  const std::vector<fmm_device_detail::Level<Node>> & levels,
                                                  FMMDeviceCache & cache)
    {
      using namespace fmm_device_detail;
      using SNode = typename SML::Node;
      auto & tg = *ml.fmm_params.device_targets;
      std::unordered_map<const Node*, int> index;
      std::vector<const Node*> leaves;
      for (auto & lev : levels)
        if (L2PLevel (ml, lev))
          for (size_t k = 0; k < lev.nfirst; k++)
            {
              index[lev.nodes[k]] = int(leaves.size());
              leaves.push_back (lev.nodes[k]);
            }
      if (leaves.empty()) return nullptr;

      size_t npts = tg.pts.Size();
      Array<int> leaf_of(npts);
      for (size_t i = 0; i < npts; i++)
        {
          const Node * node = &ml.root;
          while (node->childs[node->GetChildNum(tg.pts[i])])
            node = node->childs[node->GetChildNum(tg.pts[i])].get();
          auto it = index.find (node);
          leaf_of[i] = it == index.end() ? -1 : it->second;
        }
      auto [order, first] = SortByLeaf (leaf_of, leaves.size(), [&] (size_t l, size_t n)
        { return n == leaves[l]->targets.Size() && !leaves[l]->vol_targets.Size(); },
        "fmm device L2P: the target points do not match the tree's leaves");

      size_t ntg = order.Size();
      Array<int> device_index(npts);
      device_index = -1;
      Array<double> geo(5*ntg), xyz(3*ntg);
      double gkappa = GeometryKappa (ml.singmp->Kappa());
      for (size_t j = 0; j < ntg; j++)
        {
          int i = order[j];
          device_index[i] = int(j);
          PolarGeometry (tg.pts[i] - leaves[leaf_of[i]]->center, gkappa, &geo[5*j]);
          for (int k = 0; k < 3; k++)
            xyz[3*j+k] = tg.pts[i](k);
        }
      auto l2p = make_shared<FMMDeviceL2P> (tree, N, tg.mode, ml.singmp->Kappa(), order, geo, xyz);
      for (auto & lev : levels)
        if (L2PLevel (ml, lev))
          {
            Array<int> lfirst, lcount, lslot;
            for (size_t k = 0; k < lev.nfirst; k++)
              {
                int l = index[lev.nodes[k]];
                lfirst.Append (first[l]);
                lcount.Append (first[l+1]-first[l]);
                lslot.Append (int(k));
              }
            l2p->AddLevel (tree, lev.pool, 1/lev.nodes[0]->mp.Scale(), lfirst, lcount, lslot);
          }

      // the singular nodes of every target's singnodes (all nodes on its path), as EvaluateMP recurses: at covered targets
      // the device takes the far ones of order 0 .. max_order and the sources of the device_p2m leaves, the host all others
      auto & stree = *cache.singular.tree;
      auto * p2m = cache.singular.p2m.get();
      auto sources = [&] (const SNode & sn) { return sn.charges.Size() || sn.dipoles.Size() || sn.chargedipoles.Size(); };
      auto range = [&] (const SNode & sn, Array<int> & ranges)
      {
        auto [sp, sk] = Slot (stree, &sn.mp);
        auto [f, c] = p2m->Sources (sp, sk);
        ranges.Append (f);
        ranges.Append (c);
      };
      std::vector<char> host_levels (cache.singular.level_pools.size(), 0);
      std::vector<std::vector<std::pair<int, const SNode*>>> pairs (host_levels.size());   // per level: (target, node)
      std::vector<Array<int>> near (ntg);    // per target: the source ranges below its near interior singnodes
      std::function<void(const SNode&, Vec<3>, int)> walk = [&] (const SNode & sn, Vec<3> p, int j)
      {
        if (sources (sn))
          {
            if (j >= 0 && sn.device_p2m) range (sn, near[j]);
            return;
          }
        if (L2Norm(p-sn.center) > sn.fmm_params.eval_separation*sn.r)
          {
            if (sn.mp.Order() < 0) return;
            if (j >= 0 && sn.mp.Order() <= FMMDeviceL2P::max_order)
              pairs.at(sn.level).push_back ({ j, &sn });
            else
              host_levels.at(sn.level) = 1;
            return;
          }
        if (!sn.childs[0]) return;
        for (auto & child : sn.childs)
          walk (*child, p, j);
      };
      for (size_t i = 0; i < npts; i++)
        {
          Vec<3> p = tg.pts[i];
          for (const Node * node = &ml.root; ; node = node->childs[node->GetChildNum(p)].get())
            {
              for (auto * sn : node->singnodes)
                if (!sources (*sn))
                  walk (*sn, p, device_index[i]);
              if (!node->childs[node->GetChildNum(p)]) break;
            }
        }
      for (auto & list : pairs)
        {
          if (list.empty()) continue;
          std::stable_sort (list.begin(), list.end(), [] (auto & a, auto & b) { return a.first < b.first; });
          Array<int> etarget, efirst, ecount, pslot;
          Array<double> pgeo(5*list.size());
          int pool = -1;
          for (size_t q = 0; q < list.size(); q++)
            {
              auto [j, sn] = list[q];
              auto [sp, sk] = Slot (stree, &sn->mp);
              pool = sp;
              if (etarget.Size() == 0 || etarget.Last() != j)
                {
                  etarget.Append (j);
                  efirst.Append (int(q));
                  ecount.Append (0);
                }
              ecount.Last()++;
              pslot.Append (sk);
              PolarGeometry (tg.pts[order[j]] - sn->center, gkappa, &pgeo[5*q]);
            }
          l2p->AddM2P (stree, pool, list[0].second->mp.Scale(), etarget, efirst, ecount, pslot, pgeo);
        }
      // the near field: the source leaves among the singnodes of a covered leaf's path, the same for all its targets, in one
      // launch; the ones below near interior nodes per target in a second
      if (p2m)
        {
          Array<int> target, list, lcount, ranges;
          for (size_t l = 0; l < leaves.size(); l++)
            {
              size_t r0 = ranges.Size();
              for (const Node * node = &ml.root; ; node = node->childs[node->GetChildNum(leaves[l]->center)].get())
                {
                  for (auto * sn : node->singnodes)
                    if (sources (*sn) && sn->device_p2m)
                      range (*sn, ranges);
                  if (node == leaves[l]) break;
                }
              if (ranges.Size() > r0)
                for (int j = first[l]; j < first[l+1]; j++)
                  {
                    target.Append (j);
                    list.Append (int(r0/2)); lcount.Append (int((ranges.Size()-r0)/2));
                  }
            }
          if (target.Size())
            l2p->AddP2P (*p2m, target, list, lcount, ranges);
          target.SetSize0(); list.SetSize0(); lcount.SetSize0(); ranges.SetSize0();
          for (size_t j = 0; j < ntg; j++)
            if (near[j].Size())
              {
                target.Append (int(j));
                list.Append (int(ranges.Size()/2)); lcount.Append (int(near[j].Size()/2));
                ranges.Append (near[j]);
              }
          if (target.Size())
            l2p->AddP2P (*p2m, target, list, lcount, ranges);
        }
      cache.regular.covered = ntg == npts && cache.singular.covered
        && std::none_of (host_levels.begin(), host_levels.end(), [] (char h) { return h; });
      cache.host_m2p_levels = std::move(host_levels);
      return l2p;
    }

    // after replay B: the regular leaves the host still evaluates down, the device's target values into the targets,
    // the covered leaves marked for the host evaluation
    template <int N>
    static void AfterReplay (const ML & ml, FMMDeviceCache & cache, std::vector<fmm_device_detail::Level<Node>> & levels)
    {
      using namespace fmm_device_detail;
      auto & side = cache.regular;
      DownloadFirst<N> (*side.tree, levels, [&] (size_t l) { return side.l2p && L2PLevel (ml, levels[l]); });
      if (!side.l2p) return;
      side.l2p->Download (*ml.fmm_params.device_targets);
      for (auto & lev : levels)
        if (L2PLevel (ml, lev))
          for (size_t k = 0; k < lev.nfirst; k++)
            lev.nodes[k]->device_order = FMMDeviceL2P::max_order;
    }

    // a hit replays B before CalcMP sorts and groups the S->R records, which only a re-plan needs; a miss does nothing.
    // Only without refinement (FMM_Operator): then no walk changes the tree, so Translate's fingerprint is this one.
    static bool TryReplay (ML & ml)
    {
      using namespace fmm_device_detail;
      constexpr int N = VecLength<entry_type>;
      auto * cache = UsableCache (ml);
      // a hit: program B was recorded for this regular tree (a singular re-plan resets it)
      if (!cache || !cache->regular.program || !cache->regular.tree
          || cache->regular.fingerprint != Fingerprint (ml.root, N, std::complex<double>(ml.singmp->Kappa())))
        return false;
      auto levels = CollectLevels (ml.root, [] (const Node & node) { return !node.childs[0]; });
      AssignSlots (*cache->regular.tree, levels, cache->regular.level_pools);
      cache->regular.program->Replay();
      AfterReplay<N> (ml, *cache, levels);
      ml.root.TraverseTree ([] (Node & node) { if (node.childs[0]) node.ClearLocalExpansion(); });
      return true;
    }

    static void Translate (ML & ml, FlatArray<Array<RecSR*>> sr_groups, FlatArray<double> sr_lens,
                           FlatArray<double> sr_thetas, bool allow_refine)
    {
      using namespace fmm_device_detail;
      constexpr int N = VecLength<entry_type>;
      static_assert (sizeof(entry_type) == N*sizeof(HostComplex<Node>));

      SML & singmp = *ml.singmp;
      LocalizePlan lp;
      PlanLocalize (ml, sr_groups, allow_refine, lp);
      // leaves first: L2P (Node::Evaluate, EvaluateDirectionalDerivative) reads only leaf expansions
      auto levels = CollectLevels (ml.root, [] (const Node & node) { return !node.childs[0]; });

      if (auto * cache = UsableCache (ml))
        {
          // FMM_Operator, a miss (TryReplay replays the hits): record program B, replay it, the leaves down
          auto & side = cache->regular;
          side = {};   // free the old plan first
          cache->fast = false;
          auto tree = make_shared<FMMDeviceTree> (ml.fmm_params.device, IsDoubleEntry<entry_type>);
          side.level_pools = MakePools<N> (*tree, levels);
          FMMDevicePlan plan;
          AddGroups<N> (plan, *cache->singular.tree, *tree, levels, sr_groups, sr_lens, sr_thetas, lp);
          if (ml.fmm_params.device_targets && (side.l2p = PlanL2P<N> (ml, *tree, levels, *cache)))
            plan.SetL2P (side.l2p);
          // a regular-only re-plan (not for an FMM_Operator): the host M2P may need the levels the singular side skipped
          auto & skipped = cache->singular_skipped;
          if (std::find (skipped.begin(), skipped.end(), 1) != skipped.end())
            {
              auto slevels = CollectLevels (singmp.root, [] (const typename SML::Node & node)
                                            { return bool(node.childs[0]); });
              AssignSlots (*cache->singular.tree, slevels, cache->singular.level_pools);
              DownloadFirst<N> (*cache->singular.tree, slevels, [&] (size_t l) { return !skipped[l]; });
            }
          side.program = plan.Record (cache->singular.program.get());   // one workspace with A
          side.tree = tree;
          side.fingerprint = Fingerprint (ml.root, N, std::complex<double>(singmp.Kappa()));
          side.program->Replay();
          AfterReplay<N> (ml, *cache, levels);
          cache->fast = side.covered;
        }
      else
        {
          if (!singmp.device_tree)
            throw Exception ("fmm device: the singular expansion holds no device pools (its CalcMP did not run on the "
                             "device, or another local expansion already used them)");
          FMMDeviceTree tree (ml.fmm_params.device, IsDoubleEntry<entry_type>);
          MakePools<N> (tree, levels);
          FMMDevicePlan plan;
          AddGroups<N> (plan, *singmp.device_tree, tree, levels, sr_groups, sr_lens, sr_thetas, lp);
          plan.Run();
          // the singular pools are only needed for M2L: release them now, not when singmp dies (a potential
          // keeps its singular expansion as long as the local expansion lives)
          singmp.device_tree = nullptr;
          DownloadFirst<N> (tree, levels, [] (size_t) { return false; });
        }
      for (auto node : lp.nodes_to_clear)
        node->ClearLocalExpansion();
    }
  };
}

#endif
