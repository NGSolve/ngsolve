#ifndef FILE_FMMINTERFACE
#define FILE_FMMINTERFACE

#include <ngstd.hpp>
#include <bla.hpp>
#include <variant>

namespace ngsbem
{
  using namespace ngbla;

  // None: native path; Cpu/Gpu: device kernels (mptools_device.cpp) on the CPU device or a registered GPU
  enum class FMM_Device { None, Cpu, Gpu };

  // the device plan of one FMM_Operator direction (mptools_device.hpp)
  class FMMDeviceCache;

  // the sources of one matvec for the device P2M: point i adds the charge charge[i*comps+c], the dipole charge
  // dipole[i*comps+c] along -nv[i] and the dipole vector dipvec[(i*comps+c)*3+j] (BaseFMMInterface::DeviceSourceValues);
  // `entries` host sources per point
  struct FMMDeviceSources
  {
    int kind = 0;              // BaseFMMInterface::DeviceSourceKind
    int comps = 1, entries = 1;
    FlatArray<Vec<3>> pts, nv;
    FlatArray<Complex> charge, dipole, dipvec;
  };

  // the targets of one matvec for the device L2P: mode (BaseFMMInterface::DeviceTargetMode) bit 0 the value, bit 1
  // the gradient; the device writes values[(i*comps+c)*outs+k] (value first, then the derivatives; zero for the
  // targets it skips)
  inline int FMMDeviceTargetOuts (int mode) { return (mode & 1) + 3*(mode >> 1 & 1); }
  struct FMMDeviceTargets
  {
    int comps = 1, mode = 1;
    FlatArray<Vec<3>> pts;
    FlatArray<Complex> values;
  };

  class FMM_Parameters
  {
  public:
    int maxdirect = 100;
    int minorder = 20;
    double order_factor = 2.0;    // order = minorder + order_factor kappa r
    double separation = 2.0;
    double eval_separation = 3.0;
    double split_kr = 5.0;
    int maxlevel = 20;
    FMM_Device device = FMM_Device::None;
    // set by FMM_Operator when device != None
    shared_ptr<FMMDeviceCache> device_cache;
    // set with device_cache, if the types have a device P2M / L2P
    const FMMDeviceSources * device_sources = nullptr;
    const FMMDeviceTargets * device_targets = nullptr;
  };


  struct FMMTreeStats
  {
    int max_level = 0;                           // Deepest tree level reached by any node.
    size_t num_nodes = 0;                        // Total number of nodes in the tree, including internal nodes.
    size_t num_leaves = 0;                       // Number of leaf nodes, including empty leaves.
    size_t active_leaves = 0;                    // Number of leaves containing at least one source or target item.
    Array<size_t> nodes_per_level;               // Histogram of node counts by tree level.
    size_t leaf_size_min = std::numeric_limits<size_t>::max(); // Minimum number of source or target items in a leaf.
    size_t leaf_size_max = 0;                    // Maximum number of source or target items in a leaf.
    double leaf_size_sum = 0;                    // Sum of leaf sizes, used to compute the mean leaf size.
    int order_min = std::numeric_limits<int>::max(); // Minimum allocated spherical expansion order.
    int order_max = -1;                          // Maximum allocated spherical expansion order.
    double order_sum = 0;                        // Sum of allocated orders, used to compute the mean order.
    size_t num_allocated_multipoles = 0;         // Number of nodes whose multipole/local expansion is allocated.
    size_t total_coefficients = 0;               // Total number of stored spherical harmonic coefficients.
    size_t multipole_bytes = 0;                  // Memory used by stored spherical harmonic coefficients.
  };


  // Only tree-level operations cross this interface; coefficients stay typed.
  class NGS_DLL_HEADER BaseSingularMLExpansion
  {
  public:
    virtual ~BaseSingularMLExpansion() = default;
    virtual void CalcMP() = 0;
    virtual void CollectStatistics(FMMTreeStats & stats) const = 0;
    virtual std::variant<double,Complex> GetKappa() const = 0;
  };


  class NGS_DLL_HEADER BaseRegularMLExpansion
  {
  public:
    struct M2LCounts
    {
      size_t num_s2r = 0;                        // Number of recorded source-to-regular box translations.
      size_t num_direct_evaluations = 0;         // Estimated quadrature-point interactions handled by direct fallback.
    };

    virtual ~BaseRegularMLExpansion() = default;
    virtual void AddTarget(Vec<3> p) = 0;
    virtual void AddVolumeTarget(Vec<3> p, double r) = 0;
    virtual void CalcMP(shared_ptr<BaseSingularMLExpansion> singmp, bool onlytargets = true) = 0;
    virtual void CollectStatistics(FMMTreeStats & stats) const = 0;
    virtual M2LCounts CollectM2LStatistics(shared_ptr<BaseSingularMLExpansion> singmp) = 0;
  };


  class NGS_DLL_HEADER BaseFMMInterface
  {
  public:
    virtual ~BaseFMMInterface() = default;
    virtual bool NeedsNormal() const = 0;
    // the components of its expansions
    virtual int ExpansionComponents() const = 0;
    virtual shared_ptr<BaseSingularMLExpansion> CreateMultipoleExpansion(Vec<3> c, double r, FMM_Parameters params) const = 0;
    virtual shared_ptr<BaseRegularMLExpansion> CreateLocalExpansion(Vec<3> c, double r, FMM_Parameters params) const = 0;

    virtual void AddSource(BaseSingularMLExpansion & mp, Vec<3> p, Vec<3> nv, BareSliceVector<double> val) const
    { throw Exception("AddSource not implemented for this FMM type"); }
    virtual void AddSource(BaseSingularMLExpansion & mp, Vec<3> p, Vec<3> nv, BareSliceVector<Complex> val) const
    { throw Exception("AddSource not implemented for this FMM type"); }
    virtual void EvaluateMP(BaseRegularMLExpansion & mp, Vec<3> p, Vec<3> nv, BareSliceVector<double> val) const
    { throw Exception("EvaluateMP not implemented for this FMM type"); }
    virtual void EvaluateMP(BaseRegularMLExpansion & mp, Vec<3> p, Vec<3> nv, BareSliceVector<Complex> val) const
    { throw Exception("EvaluateMP not implemented for this FMM type"); }

    // the device P2M: bits 1 charges, 2 dipoles along -nv (3 as ChargeDipoles), 4 dipole vectors; 0 keeps the host P2M
    virtual int DeviceSourceKind() const { return 0; }
    // the host sources (charges and dipoles) AddSource adds per point
    virtual int DeviceSourceEntries() const { return 1; }
    // the charge, the dipole charge and the dipole vector AddSource adds for val at p, per expansion component,
    // computed in the complex type of its expansions
    virtual void DeviceSourceValues(BareSliceVector<double> val, Vec<3> p, FlatVector<Complex> charge,
                                    FlatVector<Complex> dipole, FlatVector<Complex> dipvec) const
    { throw Exception("DeviceSourceValues not implemented for this FMM type"); }
    virtual void DeviceSourceValues(BareSliceVector<Complex> val, Vec<3> p, FlatVector<Complex> charge,
                                    FlatVector<Complex> dipole, FlatVector<Complex> dipvec) const
    { throw Exception("DeviceSourceValues not implemented for this FMM type"); }
    // the device L2P of what EvaluateMP needs: the mode of FMMDeviceTargets, 0 keeps the host L2P
    virtual int DeviceTargetMode() const { return 0; }
    // val += EvaluateMP of the device's values and gradients dev (FMMDeviceTargets) at the point p with normal nv
    virtual void AddDeviceTargetValue(FlatVector<Complex> dev, Vec<3> p, Vec<3> nv, BareSliceVector<double> val) const
    { throw Exception("AddDeviceTargetValue not implemented for this FMM type"); }
    virtual void AddDeviceTargetValue(FlatVector<Complex> dev, Vec<3> p, Vec<3> nv, BareSliceVector<Complex> val) const
    { throw Exception("AddDeviceTargetValue not implemented for this FMM type"); }
  };
}

#endif
