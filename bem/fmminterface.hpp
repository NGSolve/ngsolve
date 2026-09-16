#ifndef FILE_FMMINTERFACE
#define FILE_FMMINTERFACE

#include <ngstd.hpp>
#include <bla.hpp>
#include <variant>

namespace ngsbem
{
  using namespace ngbla;

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
  };
}

#endif
