#ifndef FILE_METAL_BENCHMARK_HPP
#define FILE_METAL_BENCHMARK_HPP


namespace ngsmetal
{
  class Metal_MM_Benchmark
  {
    int n, m, k, lda, ldb;
    MTL::ComputePipelineState* pipelineState;

    MTL::Buffer* buffer_A;
    MTL::Buffer* buffer_B;
    MTL::Buffer* buffer_C;
    int blocks, warps;
    
  public:
    Metal_MM_Benchmark (int _n, int _m, int _k, int _lda, int _ldb);
    double Run(bool timing) const;
  };
}

#endif

