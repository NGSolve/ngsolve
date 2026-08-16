#ifndef FILE_METAL_BENCHMARK_HPP
#define FILE_METAL_BENCHMARK_HPP


namespace ngsmetal
{
  class Metal_MM_Benchmark
  {
    int n, m, k, lda, ldb;
    MTL::ComputePipelineState* pipelineState;
    
  public:
    Metal_MM_Benchmark (int _n, int _m, int _k, int _lda, int _ldb);
    double Run() const;
  };
}

#endif

