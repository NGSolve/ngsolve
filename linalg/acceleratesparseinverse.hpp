#ifndef ACCELERATESPARSEINVERSE_HPP
#define ACCELERATESPARSEINVERSE_HPP

#include "sparsefactorization_interface.hpp"

#ifdef USE_ACCELERATE_SPARSE

namespace ngla
{
  /** Sparse direct solver backed by Apple's Accelerate Sparse Solvers.

      The NGSolve matrix is converted to scalar CSC storage during Analyze().
      Real SPD, real symmetric-indefinite, and real nonsymmetric matrices use
      Cholesky, LDLT, and LU, respectively.  Apple's LU implementation requires
      macOS 15.5 or newer.
   */
  class AccelerateSparseInverse : public SparseFactorizationInterface
  {
    Array<long> column_starts;
    Array<int> row_indices;
    Array<double> csc_values;
    Array<size_t> value_map;

    Array<size_t> pattern_first;
    Array<int> pattern_columns;

    struct AppleState;
    unique_ptr<AppleState> apple;

    void BuildCSC();
    void RefreshValues();
    void CheckPattern() const;
    void Cleanup();

  public:
    explicit AccelerateSparseInverse(shared_ptr<const BaseMatrix> matrix);
    ~AccelerateSparseInverse() override;

    void Analyze() override;
    void Factor() override;
    void Solve(const BaseVector &rhs, BaseVector &solution) const override;
    void SolveTrans(const BaseVector &rhs, BaseVector &solution) const override;

    bool SupportsUpdate() const override { return true; }
    BaseMatrix::OperatorInfo GetOperatorInfo() const override
    {
      return { "Accelerate Sparse", size_t(height), size_t(width) };
    }
  };

  void RegisterAccelerateSparseInverse();
}

#endif // USE_ACCELERATE_SPARSE
#endif // ACCELERATESPARSEINVERSE_HPP
