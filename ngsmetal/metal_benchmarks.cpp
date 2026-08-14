#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>



#include "metal_btdtb.hpp"
// #include "tinybla.hpp"
extern string code_tinybla;

namespace ngsmetal
{

  static std::string Substitute(std::string src, const std::string &token,
                                const std::string &value)
  {
    for (size_t p = src.find(token); p != std::string::npos; p = src.find(token, p))
      src.replace(p, token.size(), value);
    return src;
  }
  

  double Metal_MM_Benchmark (int n, int m, int k, int lda, int ldb)
  {
    cout << "MM " << n << "x" << m << "x" << k << flush;
    int blocks = 1000;
    int warps = 4;
    int runs = 5e6 / n / m / k + 1;
    
    Matrix<float> a(n,k), b(k,m), c(n,m);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < k; j++)
        a(i,j) = i; // sin(i+j);

    for (int i = 0; i < k; i++)
      for (int j = 0; j < m; j++)
        b(i,j) = j; //  cos(3+i+j);

    c =runs *  a*b;

    MTL::Buffer* buffer_A = GetDevice()->newBuffer(n*k*sizeof(float), MTL::ResourceStorageModeShared);
    FlatMatrix<float> deva(n,k, (float*)buffer_A->contents());
    deva = a;

    MTL::Buffer* buffer_B = GetDevice()->newBuffer(k*m*sizeof(float), MTL::ResourceStorageModeShared);
    FlatMatrix<float> devb(k,m, (float*)buffer_B->contents());
    devb = b;

    MTL::Buffer* buffer_C = GetDevice()->newBuffer(blocks*n*m*sizeof(float), MTL::ResourceStorageModeShared);
    FlatMatrix<float> devc(n,m, (float*)buffer_C->contents());
    
    
    string code = code_tinybla +
      R"(
      #include <metal_stdlib>
      using namespace metal;
      using namespace tinybla;


      kernel void metal_mm_benchmark(
      device const float*  A            [[buffer(0)]],
      device const float*  B            [[buffer(1)]],
      device float* C                   [[buffer(2)]],

      uint   blockIdx           [[threadgroup_position_in_grid]],
      uint   threadIdx          [[thread_position_in_threadgroup]],
      uint   blockDim           [[threads_per_threadgroup]],
      uint   gridDim            [[threadgroups_per_grid]]
      )

      {
      constexpr int n = $N;
      constexpr int m = $M;
      constexpr int k = $K;
      constexpr int lda = $LDA;
      constexpr int ldb = $LDB;
      constexpr int runs = $RUNS;

      alignas(16) threadgroup float sharedA[n][lda];
      alignas(16) threadgroup float sharedB[k][ldb];
      threadgroup float2 sharedB2[k][ldb];
      threadgroup float4 sharedB4[k][ldb];

      for (int i = threadIdx; i < n*lda; i+=blockDim) {
        int r = i/lda;
        int c = i%lda;
        sharedA[r][c] = c < k ? A[r*k+c] : 0;
      }
      for (int i = threadIdx; i < k*ldb; i+=blockDim) {
        int r = i/ldb;
        int c = i%ldb;
        sharedB[r][c] = c < m ? B[r*m+c] : 0;
      }
      for (int i = threadIdx; i < k*ldb; i+=blockDim) {
        int r = i/ldb;
        int c = i%ldb;
        sharedB2[r][c] = c < m ? B[r*m+c] : 0;
      }
      for (int i = threadIdx; i < k*ldb; i+=blockDim) {
        int r = i/ldb;
        int c = i%ldb;
        sharedB4[r][c] = c < m ? B[r*m+c] : 0;
      }

      threadgroup_barrier(mem_flags::mem_threadgroup);

/*
      WarpMatrix<n,m,float> sum = 0.0;

      auto ma = MakeBareMatrix<RowMajor> (&sharedA[0][0], lda);
      auto mb = MakeBareMatrix<RowMajor> (&sharedB[0][0], ldb);
      auto mc = MakeBareMatrix<RowMajor>(C + blockIdx * n * m, m);

      for (int i = 0; i < runs; i++)
        sum.AddMM<k> (ma.Transpose(), mb, threadIdx);

      sum.Store(mc, threadIdx);
*/

      WarpMatrix<n,m/2,float2> sum = 0.0;

      auto ma = MakeBareMatrix<RowMajor> (&sharedA[0][0], lda);
      auto mb = MakeBareMatrix<RowMajor> (&sharedB2[0][0], ldb/2);
      auto mc = MakeBareMatrix<RowMajor> ((device float2*)C, m/2);

      for (int i = 0; i < runs; i++)
        sum.AddMM<k> (ma, mb, threadIdx);

      sum.Store(mc, threadIdx);

/*
      WarpMatrix<n,m/4,float4> sum = 0.0;

      auto ma = MakeBareMatrix<ColMajor> (&sharedA[0][0], lda);
      auto mb = MakeBareMatrix<RowMajor> (&sharedB4[0][0], ldb/4);
      auto mc = MakeBareMatrix<RowMajor> ((device float4*)C, m/4);

      for (int i = 0; i < runs; i++)
        sum.AddMM<k> (ma, mb, threadIdx);

      sum.Store(mc, threadIdx);
*/
      }
       )";
    
    
    code = Substitute(code, "$N", ToString(n)+" /*n*/ ");
    code = Substitute(code, "$M", ToString(m)+" /*m*/ ");
    code = Substitute(code, "$K", ToString(k)+" /*k*/ ");
    code = Substitute(code, "$LDA", ToString(lda)+" /*lda*/ ");
    code = Substitute(code, "$LDB", ToString(ldb)+" /*ldb*/ ");
    code = Substitute(code, "$RUNS", ToString(runs)+" /*runs*/ ");

    
    ofstream codefile("metal_mm_benchmark.metal");
    codefile << code;
    codefile.close();
    
    NS::Error* error = nullptr;
    NS::String* mslCode = NS::String::string(code.c_str(), NS::UTF8StringEncoding);

    MTL::Library* library = GetDevice()->newLibrary(mslCode, nullptr, &error);

    if (!library)
      throw Exception("Shader compile error: " 
                      +string(error->localizedDescription()->utf8String()));

    NS::String* funcName = NS::String::string("metal_mm_benchmark", NS::UTF8StringEncoding);
    auto MM_Func = library->newFunction(funcName);

    auto pipelineState = GetDevice()->newComputePipelineState(MM_Func, &error);
    if (error)
      {
        std::cerr << "Metal Error: " << error->localizedDescription()->utf8String() << std::endl;
        throw Exception (error->localizedDescription()->utf8String());
      }




    
    double durationMs;

    for (int j = 0; j < 10; j++) {

      MTL::CommandBuffer* commandBuffer = GetCommandQueue()->commandBuffer();
      MTL::ComputeCommandEncoder* encoder = commandBuffer->computeCommandEncoder();
      
      encoder->setComputePipelineState(pipelineState);
      encoder->setBuffer(buffer_A, 0, 0);
      encoder->setBuffer(buffer_B, 0, 1);
      encoder->setBuffer(buffer_C, 0, 2);
      
      encoder->dispatchThreadgroups(MTL::Size(blocks,1,1), MTL::Size(warps*32, 1, 1));
      encoder->endEncoding();
      
      commandBuffer->commit();
      commandBuffer->waitUntilCompleted();

      double startTime = commandBuffer->GPUStartTime();
      double endTime = commandBuffer->GPUEndTime();
      durationMs = (endTime - startTime) * 1000.0;
      
      // cout << "time = " << durationMs << endl;
      
    }
    // sleep(5);
    // cout << "c = " << endl << c << endl;
    //    cout << "devc = " << endl << devc << endl;
    // cout << "err = " << endl << L2Norm(devc-c) / L2Norm(c) << endl;
    // std::cout << "Kernel GPU Runtime: " << durationMs << " ms\n";
    double gflops = blocks*warps*2.0*runs*n*m*k / durationMs * 1000 * 1e-9;
    cout << " GFlops = " << gflops << endl;
    return gflops;
  }  

}
