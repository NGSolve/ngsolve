#include "fmmoperator.hpp"
#include <map>
#include <mutex>

namespace ngsbem
{
  namespace
  {
    Timer<> GetFMMTimer(const string & name)
    {
      // Register each runtime kernel name once, not once per matrix instance.
      static std::mutex mutex;
      static std::map<string,Timer<>> timers;
      std::lock_guard<std::mutex> guard(mutex);
      return timers.try_emplace(name, name).first->second;
    }

    // the device P2M input of a matvec (x: one row per point), if the cache is used and the source type has one
    template <typename TSCAL>
    FMMDeviceSources DeviceSources (const FMM_Parameters & params, const BaseFMMInterface & iface, FlatArray<Vec<3>> pts,
                                    FlatArray<Vec<3>> nv, FlatMatrix<TSCAL> x, Array<Complex> & charge,
                                    Array<Complex> & dipole, Array<Complex> & dipvec)
    {
      int kind = params.device_cache ? iface.DeviceSourceKind() : 0;
      if (!kind) return {};
      size_t comps = iface.ExpansionComponents();
      size_t nc = (kind & 1) ? comps : 0, nd = (kind & 2) ? comps : 0, nvec = (kind & 4) ? 3*comps : 0;
      charge.SetSize (pts.Size()*nc);
      dipole.SetSize (pts.Size()*nd);
      dipvec.SetSize (pts.Size()*nvec);
      ParallelFor (pts.Size(), [&](size_t i) {
        iface.DeviceSourceValues (x.Row(i), pts[i], FlatVector<Complex>(nc, charge.Data()+i*nc),
                                  FlatVector<Complex>(nd, dipole.Data()+i*nd),
                                  FlatVector<Complex>(nvec, dipvec.Data()+i*nvec));
      });
      return { kind, int(comps), iface.DeviceSourceEntries(), pts, nv, charge, dipole, dipvec };
    }

    // the device L2P output of a matvec (zeroed), if the cache is used and the target type has a device mode
    FMMDeviceTargets DeviceTargets (const FMM_Parameters & params, const BaseFMMInterface & iface, FlatArray<Vec<3>> pts,
                                    Array<Complex> & values)
    {
      int mode = params.device_cache ? iface.DeviceTargetMode() : 0;
      if (!mode) return {};
      int comps = iface.ExpansionComponents();
      values.SetSize (pts.Size()*comps*FMMDeviceTargetOuts (mode));
      values = Complex(0.0);
      return { comps, mode, pts, values };
    }

    // the points, normals and root box of one side of an FMM_Operator
    struct FMMSide
    {
      const BaseFMMInterface & iface;
      FlatArray<Vec<3>> pts, nv;
      Vec<3> c;
      double r;
    };

    // y = the FMM of x from the sources s to the targets t (Mult and MultTrans); with a device, the direction's plan
    // cache (locked for the whole call), the device P2M input of x and the device L2P output; teval times the
    // evaluation
    template <typename TSCAL>
    void ApplyFMM (const FMM_Parameters & fmm_params, const shared_ptr<FMMDeviceCache> & cache, std::mutex & mutex,
                   FMMSide s, FMMSide t, FlatMatrix<TSCAL> x, FlatMatrix<TSCAL> y, const Timer<> * teval)
    {
      auto lock = cache ? std::unique_lock<std::mutex>(mutex) : std::unique_lock<std::mutex>();
      FMM_Parameters params = fmm_params;
      params.device_cache = cache;
      Array<Complex> charge, dipole, dipvec, values;
      FMMDeviceSources sources = DeviceSources (params, s.iface, s.pts, s.nv, x, charge, dipole, dipvec);
      FMMDeviceTargets targets = DeviceTargets (params, t.iface, t.pts, values);
      if (sources.kind) params.device_sources = &sources;
      if (targets.pts.Size()) params.device_targets = &targets;
      size_t nout = targets.comps*FMMDeviceTargetOuts (targets.mode);

      y = 0;
      // without host trees if the device plan covers everything (FMMDeviceCache)
      shared_ptr<BaseRegularMLExpansion> regmp;
      if (!FMMDeviceMatvec (params))
        {
          auto singmp = s.iface.CreateMultipoleExpansion(s.c, s.r, params);
          ParallelFor (s.pts.Size(), [&](int i){
            s.iface.AddSource(*singmp, s.pts[i], s.nv[i], x.Row(i));
          });
          singmp->CalcMP();
          regmp = t.iface.CreateLocalExpansion(t.c, t.r, params);
          ParallelFor (t.pts.Size(), [&](int i){
            regmp->AddTarget(t.pts[i]);
          });
          regmp->CalcMP(singmp);
        }

      if (teval) teval->Start();
      ParallelFor (t.pts.Size(), [&](int i) {
        if (regmp)
          t.iface.EvaluateMP(*regmp, t.pts[i], t.nv[i], y.Row(i));
        if (params.device_targets)
          t.iface.AddDeviceTargetValue (FlatVector<Complex>(nout, &values[i*nout]), t.pts[i], t.nv[i], y.Row(i));
      }, TasksPerThread(10));
      if (teval) teval->Stop();
    }
  }

  template <typename TSCAL>
  FMM_Operator<TSCAL> ::
  FMM_Operator(string name, IVec<2> shape,
               shared_ptr<const BaseFMMInterface> asource, shared_ptr<const BaseFMMInterface> atarget,
               Array<Vec<3>> xpts, Array<Vec<3>> ypts, Array<Vec<3>> xnv, Array<Vec<3>> ynv,
               const FMM_Parameters & fmm_params)
    : BASE(std::move(xpts), std::move(ypts), std::move(xnv), std::move(ynv), shape, fmm_params),
      kernel_name(std::move(name)), source(std::move(asource)), target(std::move(atarget)),
      tapply(GetFMMTimer("ngbem fmm apply "+kernel_name)), teval(GetFMMTimer("ngbem fmm apply "+kernel_name+" eval")),
      ttrans(GetFMMTimer("ngbem fmm apply Trans "+kernel_name))
  {
    if (fmm_params.device != FMM_Device::None)
      for (auto & cache : device_cache)
        cache = MakeFMMDeviceCache();
  }


  template <typename TSCAL>
  void FMM_Operator<TSCAL> :: Mult(const BaseVector & x, BaseVector & y) const
  {
    RegionTimer reg(tapply);

    auto shape = kernelshape;
    auto matx = x.FV<TSCAL>().AsMatrix(xpts.Size(), shape[1]);
    auto maty = y.FV<TSCAL>().AsMatrix(ypts.Size(), shape[0]);

    ApplyFMM (fmm_params, device_cache[0], device_mutex[0], { *source, xpts, xnv, cx, rx },
              { *target, ypts, ynv, cy, ry }, matx, maty, &teval);
  }


  template <typename TSCAL>
  void FMM_Operator<TSCAL> :: MultTrans(const BaseVector & x, BaseVector & y) const
  {
    RegionTimer reg(ttrans);

    auto shape = kernelshape;
    auto matx = x.FV<TSCAL>().AsMatrix(ypts.Size(), shape[0]);
    auto maty = y.FV<TSCAL>().AsMatrix(xpts.Size(), shape[1]);

    ApplyFMM (fmm_params, device_cache[1], device_mutex[1], { *target, ypts, ynv, cy, ry },
              { *source, xpts, xnv, cx, rx }, matx, maty, nullptr);
  }


  template <typename TSCAL>
  BaseMatrix::OperatorInfo FMM_Operator<TSCAL> :: GetOperatorInfo () const
  {
    return { string("FMM_Operator ")+kernel_name, this->Height(), this->Width() };
  }


  template <typename TSCAL>
  FMMOperatorInfo FMM_Operator<TSCAL> :: GetFMMInfo () const
  {
    static Timer t("ngbem fmm build info trees"); RegionTimer reg(t);

    auto singmp = source->CreateMultipoleExpansion(cx, rx, fmm_params);
    auto shape = kernelshape;
    Vector<TSCAL> zero_vals(shape[1]);
    zero_vals = TSCAL(0);
    for (size_t i = 0; i < xpts.Size(); i++)
      source->AddSource(*singmp, xpts[i], xnv[i], zero_vals);

    auto regmp = target->CreateLocalExpansion(cy, ry, fmm_params);
    for (size_t i = 0; i < ypts.Size(); i++)
      regmp->AddTarget(ypts[i]);

    // S2R walk: allocates target multipoles and counts translations + direct-evaluation pairs.
    auto counts = regmp->CollectM2LStatistics(singmp);

    FMMOperatorInfo info;
    info.kernel_name = kernel_name;
    info.source_size = xpts.Size();
    info.target_size = ypts.Size();
    info.kappa = singmp->GetKappa();
    info.parameters = fmm_params;

    singmp->CollectStatistics(info.source_tree);
    regmp->CollectStatistics(info.target_tree);

    info.total_num_nodes  = info.source_tree.num_nodes  + info.target_tree.num_nodes;
    info.total_num_leaves = info.source_tree.num_leaves + info.target_tree.num_leaves;
    info.total_multipole_coefficients = info.source_tree.total_coefficients + info.target_tree.total_coefficients;
    info.total_memory_bytes = info.source_tree.multipole_bytes + info.target_tree.multipole_bytes;

    info.num_s2r = counts.num_s2r;
    info.num_direct_evaluations = counts.num_direct_evaluations;

    info.source_bbox_radius = rx;
    info.target_bbox_radius = ry;

    double dense = double(info.source_size) * double(info.target_size);
    info.direct_fallback_fraction = (dense > 0) ? double(info.num_direct_evaluations) / dense : 0.0;

    return info;
  }

  template class FMM_Operator<double>;
  template class FMM_Operator<Complex>;
}
