#include "acceleratesparseinverse.hpp"

#ifdef USE_ACCELERATE_SPARSE

#ifndef __VECLIB__
#define NGSOLVE_DEFINED_VECLIB
#define __VECLIB__
#endif
#ifndef __ACCELERATE__
#define NGSOLVE_DEFINED_ACCELERATE
#endif
#include <vecLib/cblas.h>
#include <vecLib/Sparse/Solve.h>
#ifdef NGSOLVE_DEFINED_VECLIB
#undef __VECLIB__
#undef NGSOLVE_DEFINED_VECLIB
#endif
#ifdef NGSOLVE_DEFINED_ACCELERATE
#undef __ACCELERATE__
#undef NGSOLVE_DEFINED_ACCELERATE
#endif

#include <cstdlib>
#include <limits>

namespace ngla
{
  struct AccelerateSparseInverse::AppleState
  {
    SparseMatrix_Double sparse_matrix = {};
    SparseOpaqueSymbolicFactorization symbolic_factor = {};
    SparseOpaqueFactorization_Double numeric_factor = {};
    SparseFactorization_t factorization_type = SparseFactorizationLDLT;
    bool have_symbolic_factor = false;
    bool have_numeric_factor = false;
  };

  namespace
  {
    thread_local string accelerate_error_message;

    void ReportAccelerateError(const char *message)
    {
      accelerate_error_message = message ? message : "unknown parameter error";
    }

    void ClearAccelerateError()
    {
      accelerate_error_message.clear();
    }

    string StatusText(SparseStatus_t status)
    {
      switch (status)
      {
      case SparseStatusOK: return "success";
      case SparseFactorizationFailed: return "factorization failed";
      case SparseMatrixIsSingular: return "matrix is singular";
      case SparseInternalError: return "internal error";
      case SparseParameterError: return "parameter error";
      case SparseStatusReleased: return "factorization was already released";
      default: return "unknown status " + ToString(int(status));
      }
    }

    void CheckAccelerateStatus(SparseStatus_t status, const char *operation)
    {
      if (status == SparseStatusOK && accelerate_error_message.empty())
        return;

      string message = string("Accelerate Sparse: ") + operation + ": " +
                       StatusText(status);
      if (!accelerate_error_message.empty())
        message += " (" + accelerate_error_message + ")";
      throw Exception(message);
    }

    SparseSymbolicFactorOptions SymbolicOptions()
    {
      return { SparseDefaultControl, SparseOrderDefault, nullptr, nullptr,
               std::malloc, std::free, ReportAccelerateError };
    }
  }

  AccelerateSparseInverse::AccelerateSparseInverse(
      shared_ptr<const BaseMatrix> matrix)
    : SparseFactorizationInterface(std::move(matrix)),
      apple(make_unique<AppleState>())
  { }

  AccelerateSparseInverse::~AccelerateSparseInverse()
  {
    Cleanup();
  }

  void AccelerateSparseInverse::Cleanup()
  {
    if (apple->have_numeric_factor)
    {
      SparseCleanup(apple->numeric_factor);
      apple->have_numeric_factor = false;
    }
    if (apple->have_symbolic_factor)
    {
      SparseCleanup(apple->symbolic_factor);
      apple->have_symbolic_factor = false;
    }
  }

  void AccelerateSparseInverse::BuildCSC()
  {
    auto matrix = dynamic_pointer_cast<const SparseMatrixTM<double>>(inner_mat);
    if (!matrix)
      throw Exception("Accelerate Sparse supports real double matrices only");

    if (matrix->Height() != matrix->Width())
      throw Exception("Accelerate Sparse requires a square matrix");
    if (matrix->Height() > size_t(std::numeric_limits<int>::max()))
      throw Exception("Accelerate Sparse matrix dimension exceeds INT_MAX");

    const bool symmetric = is_symmetric.IsTrue();
    auto first = matrix->GetFirstArray();
    auto columns = matrix->GetColIndices();
    const size_t n = matrix->Height();

    pattern_first.SetSize(first.Size());
    for (auto i : Range(first.Size()))
      pattern_first[i] = first[i];
    pattern_columns.SetSize(columns.Size());
    for (auto i : Range(columns.Size()))
      pattern_columns[i] = columns[i];

    column_starts.SetSize(n+1);
    column_starts = 0;
    for (size_t row = 0; row < n; row++)
      for (size_t position = first[row]; position < first[row+1]; position++)
      {
        int column = columns[position];
        if (column < 0 || size_t(column) >= n)
          throw Exception("Accelerate Sparse encountered an invalid column index");
        if (!symmetric || row >= size_t(column))
          column_starts[column+1]++;
      }

    for (size_t column = 0; column < n; column++)
    {
      if (column_starts[column+1] >
          std::numeric_limits<long>::max() - column_starts[column])
        throw Exception("Accelerate Sparse number of entries exceeds LONG_MAX");
      column_starts[column+1] += column_starts[column];
    }

    const size_t nze = size_t(column_starts[n]);
    row_indices.SetSize(nze);
    csc_values.SetSize(nze);
    value_map.SetSize(nze);

    Array<long> next(n+1);
    next = column_starts;
    for (size_t row = 0; row < n; row++)
      for (size_t position = first[row]; position < first[row+1]; position++)
      {
        int column = columns[position];
        if (symmetric && row < size_t(column))
          continue;
        size_t destination = size_t(next[column]++);
        row_indices[destination] = int(row);
        value_map[destination] = position;
      }

    SparseAttributes_t attributes = {};
    attributes.transpose = false;
    attributes.triangle = SparseLowerTriangle;
    attributes.kind = symmetric ? SparseSymmetric : SparseOrdinary;

    apple->sparse_matrix.structure = {
      int(n), int(n), column_starts.Data(), row_indices.Data(), attributes, 1
    };
    apple->sparse_matrix.data = csc_values.Data();
  }

  void AccelerateSparseInverse::CheckPattern() const
  {
    auto matrix = dynamic_pointer_cast<const SparseMatrixTM<double>>(inner_mat);
    if (!matrix)
      throw Exception("Accelerate Sparse supports real double matrices only");

    auto first = matrix->GetFirstArray();
    auto columns = matrix->GetColIndices();
    if (first.Size() != pattern_first.Size() ||
        columns.Size() != pattern_columns.Size())
      throw Exception("Accelerate Sparse Update requires an unchanged sparsity pattern");

    for (auto i : Range(first.Size()))
      if (first[i] != pattern_first[i])
        throw Exception("Accelerate Sparse Update requires an unchanged sparsity pattern");
    for (auto i : Range(columns.Size()))
      if (columns[i] != pattern_columns[i])
        throw Exception("Accelerate Sparse Update requires an unchanged sparsity pattern");
  }

  void AccelerateSparseInverse::RefreshValues()
  {
    CheckPattern();
    auto matrix = dynamic_pointer_cast<const SparseMatrixTM<double>>(inner_mat);
    auto values = matrix->GetValues();
    for (auto i : Range(value_map))
      csc_values[i] = values[value_map[i]];
  }

  void AccelerateSparseInverse::Analyze()
  {
    static Timer timer("Accelerate Sparse - Analyze");
    RegionTimer region(timer);

    if (is_complex)
      throw Exception("Accelerate Sparse currently supports real matrices only");

    BuildCSC();
    RefreshValues();

    if (IsSPD())
      apple->factorization_type = SparseFactorizationCholesky;
    else if (is_symmetric.IsTrue())
      apple->factorization_type = SparseFactorizationLDLT;
    else
    {
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 150500
      if (__builtin_available(macOS 15.5, *))
        apple->factorization_type = SparseFactorizationLU;
      else
        throw Exception("Accelerate Sparse LU requires macOS 15.5 or newer");
#else
      throw Exception("Accelerate Sparse LU requires an SDK with macOS 15.5 support");
#endif
    }

    ClearAccelerateError();
    apple->symbolic_factor = SparseFactor(apple->factorization_type,
                                          apple->sparse_matrix.structure,
                                          SymbolicOptions());
    apple->have_symbolic_factor = true;
    CheckAccelerateStatus(apple->symbolic_factor.status,
                          "symbolic factorization");
  }

  void AccelerateSparseInverse::Factor()
  {
    static Timer timer("Accelerate Sparse - Factor");
    static Timer refactor_timer("Accelerate Sparse - Refactor");
    RegionTimer region(apple->have_numeric_factor ? refactor_timer : timer);

    const bool refactor = apple->have_numeric_factor;
    RefreshValues();
    ClearAccelerateError();
    if (!apple->have_numeric_factor)
    {
      apple->numeric_factor = SparseFactor(apple->symbolic_factor,
                                           apple->sparse_matrix);
      apple->have_numeric_factor = true;
    }
    else
      SparseRefactor(apple->sparse_matrix, &apple->numeric_factor);

    CheckAccelerateStatus(apple->numeric_factor.status,
                          refactor ? "numeric refactorization"
                                   : "numeric factorization");
  }

  void AccelerateSparseInverse::Solve(const BaseVector &rhs,
                                      BaseVector &solution) const
  {
    static Timer timer("Accelerate Sparse - Solve");
    RegionTimer region(timer);

    if (!apple->have_numeric_factor)
      throw Exception("Accelerate Sparse Solve called before factorization");
    auto input = rhs.FV<double>();
    auto output = solution.FV<double>();
    if (input.Size() != inner_height || output.Size() != inner_width)
      throw Exception("Accelerate Sparse Solve received incompatible vector sizes");

    DenseVector_Double b = { int(input.Size()),
                             const_cast<double *>(input.Data()) };
    DenseVector_Double x = { int(output.Size()), output.Data() };
    SparseSolve(apple->numeric_factor, b, x);
  }

  void AccelerateSparseInverse::SolveTrans(const BaseVector &rhs,
                                           BaseVector &solution) const
  {
    if (is_symmetric.IsTrue())
      return Solve(rhs, solution);

    static Timer timer("Accelerate Sparse - SolveTrans");
    RegionTimer region(timer);
    if (!apple->have_numeric_factor)
      throw Exception("Accelerate Sparse SolveTrans called before factorization");
    auto input = rhs.FV<double>();
    auto output = solution.FV<double>();
    if (input.Size() != inner_width || output.Size() != inner_height)
      throw Exception("Accelerate Sparse SolveTrans received incompatible vector sizes");

    auto transpose_factor = SparseGetTranspose(apple->numeric_factor);
    try
    {
      CheckAccelerateStatus(transpose_factor.status,
                            "transpose factorization view");
      DenseVector_Double b = { int(input.Size()),
                               const_cast<double *>(input.Data()) };
      DenseVector_Double x = { int(output.Size()), output.Data() };
      SparseSolve(transpose_factor, b, x);
    }
    catch (...)
    {
      SparseCleanup(transpose_factor);
      throw;
    }
    SparseCleanup(transpose_factor);
  }

  void RegisterAccelerateSparseInverse()
  {
    static const bool registered = []
    {
      BaseMatrix::RegisterInverseCreator(
        "acceleratesparse",
        [](shared_ptr<BaseMatrix> matrix, shared_ptr<BitArray> subset,
           shared_ptr<const Array<int>> clusters) -> shared_ptr<BaseMatrix>
        {
          auto inverse = make_shared<AccelerateSparseInverse>(matrix);
          inverse->SetSubset(std::move(subset), std::move(clusters));
          inverse->Update();
          return inverse;
        });
      return true;
    }();
    (void)registered;
  }
}

#endif // USE_ACCELERATE_SPARSE
