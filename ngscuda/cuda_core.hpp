#ifndef CUDA_CORE_HPP
#define CUDA_CORE_HPP

#include <cuda_runtime.h>

#include <core/array.hpp>
#include <core/exception.hpp>

namespace ngla { void EnsureCuBlasWorkspace(); }

namespace ngs_cuda
{
  // Forward declaration for CudaWhileGraph
  void ConvergenceCheck(double* rz, double tol, cudaGraphConditionalHandle handle, int* iter_count, int maxsteps);


  // from CUDA C++ Programming Guide:
  // https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
#ifdef __CUDA_ARCH__
#if __CUDA_ARCH__ < 600
  inline __device__ double atomicAdd(double* address, double val)
  {
    unsigned long long int* address_as_ull =
      (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed,
                      __double_as_longlong(val +
                                           __longlong_as_double(assumed)));

      // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
  }
#endif
#endif



  extern cudaStream_t ngs_cuda_stream;
  

  // Kernel wrapper only available if we are compiling the current file with the cuda compiler
#ifdef __CUDACC__
   
  template<class F> __global__
  void CUDA_forall(int n, F f)
  {
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
    for (int i = tid; i < n; i += blockDim.x*gridDim.x)
      f(i);
  }

  template<class F> __global__
  void CUDA_forall2(int n, F f)
  {
    int tid = (blockIdx.x*blockDim.y+threadIdx.y)*blockDim.x+threadIdx.x;
    for (int i = tid; i < n; i += blockDim.x*blockDim.y*gridDim.x)
      f(i);
  }
    
#define DEVICE_LAMBDA __device__

  template <class F>
  inline void DeviceParallelFor (int n, F f)
  {
    // CUDA_forall<<<512,256>>> (n, f);
    CUDA_forall<<<n/256+1,256, 0, ngs_cuda_stream>>> (n, f);
    // CUDA_forall<<<4096,32>>> (n, f);           // slower
    // CUDA_forall2<<<512,dim3(16,16)>>> (n, f);  // same performance
  }   

#endif // __CUDACC__

  

  class CudaGraph
  {
    cudaGraph_t graph = nullptr;
    cudaGraphExec_t instance = nullptr;
    cudaStream_t stream;
    cudaStream_t prev_stream;
    bool capture_ok = false;

  public:
    static inline std::function<void(cudaStream_t)> stream_change_callback = nullptr;

    CudaGraph()
    {
      auto err = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
      if (err != cudaSuccess)
        throw ngstd::Exception(std::string("[CudaGraph] cudaStreamCreate FAILED: ")
                               + cudaGetErrorString(err));
    }

    ~CudaGraph()
    {
      if (instance) cudaGraphExecDestroy(instance);
      if (graph)    cudaGraphDestroy(graph);
      cudaStreamDestroy(stream);
    }

    void BeginCapture()
    {
      capture_ok = false;
      ngla::EnsureCuBlasWorkspace();
      prev_stream = ngs_cuda_stream;
      ngs_cuda_stream = stream;
      if (stream_change_callback) {
        stream_change_callback(ngs_cuda_stream);
      }
      auto err = cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
      if (err != cudaSuccess)
        throw ngstd::Exception(std::string("[CudaGraph] cudaStreamBeginCapture FAILED: ")
                               + cudaGetErrorString(err));
    }

    void EndCapture()
    {
      auto err = cudaStreamEndCapture(stream, &graph);
      if (err != cudaSuccess) {
        ngs_cuda_stream = prev_stream;
        throw ngstd::Exception(std::string("[CudaGraph] cudaStreamEndCapture FAILED: ")
                               + cudaGetErrorString(err));
      }

      size_t numnodes = 0;
      cudaGraphGetNodes(graph, nullptr, &numnodes);
      std::cerr << "[CudaGraph] captured nodes: " << numnodes << std::endl;
      if (numnodes == 0)
        std::cerr << "[CudaGraph] WARNING: 0 nodes — ops may not be on capture stream!" << std::endl;

      err = cudaGraphInstantiate(&instance, graph, NULL, NULL, 0);
      if (err != cudaSuccess) {
        ngs_cuda_stream = prev_stream;
        throw ngstd::Exception(std::string("[CudaGraph] cudaGraphInstantiate FAILED: ")
                               + cudaGetErrorString(err));
      }

      capture_ok = true;
      ngs_cuda_stream = prev_stream;
      if (stream_change_callback) {
        stream_change_callback(ngs_cuda_stream);
      }
    }

    void Launch()
    {
      if (!capture_ok || !instance) {
        return;
      }
      auto err = cudaGraphLaunch(instance, ngs_cuda_stream);
      if (err != cudaSuccess)
        std::cerr << "[CudaGraph] cudaGraphLaunch FAILED: "
                  << cudaGetErrorString(err) << std::endl;
    }

    bool IsValid() const { return capture_ok && instance != nullptr; }
    cudaGraph_t GetGraph() const { return graph; }
  };
  





  class CudaWhileGraph
  {
    cudaGraph_t outer_graph = nullptr;
    cudaGraph_t body_graph  = nullptr;
    cudaGraphExec_t instance = nullptr;
    cudaStream_t capture_stream;
    cudaGraphConditionalHandle handle;
    bool capture_ok = false;

  public:

    CudaWhileGraph()
    {
      auto err = cudaStreamCreate(&capture_stream);
      if (err != cudaSuccess)
        throw ngstd::Exception(std::string("[CudaWhileGraph] cudaStreamCreate FAILED: ")
                               + cudaGetErrorString(err));
    }

    ~CudaWhileGraph()
    {
      if (instance)    cudaGraphExecDestroy(instance);
      if (outer_graph) cudaGraphDestroy(outer_graph);
      cudaStreamDestroy(capture_stream);
    }

    // Build WHILE graph using an existing captured graph as the iteration body
    void Build(cudaGraph_t iteration_graph, double* rz_dev_ptr, double tol, int* iter_count, int maxsteps)
    {
      capture_ok = false;
      ngla::EnsureCuBlasWorkspace();

      // 1. Create outer graph
      cudaGraphCreate(&outer_graph, 0);

      // 2. Create conditional handle with default=1 (do-while)
      cudaGraphConditionalHandleCreate(&handle, outer_graph, 1,
                                       cudaGraphCondAssignDefault);

      // // 3. Add WHILE node
      // cudaGraphNode_t while_node;
      // cudaGraphNodeParams cParams = {};
      // cParams.type               = cudaGraphNodeTypeConditional;
      // cParams.conditional.handle = handle;
      // cParams.conditional.type   = cudaGraphCondTypeWhile;
      // cParams.conditional.size   = 1;

      
      // // 6-arg pDependencyData variant exists only in CUDA 12.3–12.8; 12.9+ reverts to 5-arg
      // #if CUDART_VERSION >= 12030 && CUDART_VERSION < 12090
      //   cudaGraphAddNode(&while_node, outer_graph, nullptr, nullptr, 0, &cParams);
      // #else
      //   cudaGraphAddNode(&while_node, outer_graph, nullptr, 0, &cParams);
      // #endif
      // body_graph = cParams.conditional.phGraph_out[0];

      // // 4. Add iteration body as child graph node in body
      // //    Child graphs ARE allowed in conditional bodies
      // cudaGraphNode_t child_node;
      // auto err = cudaGraphAddChildGraphNode(&child_node, body_graph, nullptr, 0, iteration_graph);
      // if (err != cudaSuccess)
      //   throw ngstd::Exception(
      //     std::string("[CudaWhileGraph] cudaGraphAddChildGraphNode FAILED: ")
      //     + cudaGetErrorString(err));

      // 3. Add WHILE node
      cudaGraphNode_t while_node;
      cudaGraphNodeParams cParams = {};
      cParams.type               = cudaGraphNodeTypeConditional;
      cParams.conditional.handle = handle;
      cParams.conditional.type   = cudaGraphCondTypeWhile;
      cParams.conditional.size   = 1;
      #if CUDART_VERSION >= 12030
        cudaGraphAddNode(&while_node, outer_graph, nullptr, nullptr, 0, &cParams);
      #else
        cudaGraphAddNode(&while_node, outer_graph, nullptr, 0, &cParams);
      #endif
      body_graph = cParams.conditional.phGraph_out[0];

      // 4. Add iteration body as child graph node in body
      //    Child graphs ARE allowed in conditional bodies
      cudaGraphNode_t child_node;
      #if CUDART_VERSION >= 12030
        auto err = cudaGraphAddChildGraphNode(&child_node, body_graph, nullptr, 0, iteration_graph);
      #else
        auto err = cudaGraphAddChildGraphNode(&child_node, body_graph, nullptr, nullptr, 0, iteration_graph);
      #endif
      if (err != cudaSuccess)
        throw ngstd::Exception(
          std::string("[CudaWhileGraph] cudaGraphAddChildGraphNode FAILED: ")
          + cudaGetErrorString(err));


      // 5. Capture convergence kernel into body AFTER child node
      //    Custom kernel only - allowed in conditional bodies
      err = cudaStreamBeginCaptureToGraph(capture_stream, body_graph,
                                          &child_node, nullptr, 1,
                                          cudaStreamCaptureModeGlobal);
      if (err != cudaSuccess)
        throw ngstd::Exception(
          std::string("[CudaWhileGraph] cudaStreamBeginCaptureToGraph FAILED: ")
          + cudaGetErrorString(err));

      // Redirect ngs_cuda_stream so ConvergenceCheck goes to capture_stream
      cudaStream_t saved_stream = ngs_cuda_stream;
      ngs_cuda_stream = capture_stream;
      ConvergenceCheck(rz_dev_ptr, tol, handle, iter_count, maxsteps);
      ngs_cuda_stream = saved_stream;

      err = cudaStreamEndCapture(capture_stream, nullptr);
      if (err != cudaSuccess)
        throw ngstd::Exception(
          std::string("[CudaWhileGraph] cudaStreamEndCapture FAILED: ")
          + cudaGetErrorString(err));

      // 6. Debug: count body nodes
      size_t body_nodes = 0;
      cudaGraphGetNodes(body_graph, nullptr, &body_nodes);
      std::cerr << "[CudaWhileGraph] body graph nodes: " << body_nodes << std::endl;

      // 7. Instantiate outer graph
      err = cudaGraphInstantiate(&instance, outer_graph, NULL, NULL, 0);
      if (err != cudaSuccess)
        throw ngstd::Exception(
          std::string("[CudaWhileGraph] cudaGraphInstantiate FAILED: ")
          + cudaGetErrorString(err));

      capture_ok = true;
      std::cerr << "[CudaWhileGraph] Build successful" << std::endl;
    }

    void Launch()
    {
      if (!capture_ok || !instance) {
        return;
      }
      auto err = cudaGraphLaunch(instance, ngs_cuda_stream);
      if (err != cudaSuccess)
        std::cerr << "[CudaWhileGraph] cudaGraphLaunch FAILED: "
                  << cudaGetErrorString(err) << std::endl;
    }
    bool IsValid() const { return capture_ok && instance != nullptr; }
  };

  template <typename T>
  class Dev 
  {
    T data;
    
  public:
    __host__ __device__ Dev() = delete;
    
    static Dev<T> * Malloc(size_t size)
    {
      Dev<T> * ptr;
      if (auto err = cudaMalloc (&ptr, size*sizeof(T)))
        throw ngstd::Exception("cudaMalloc error, ec="+ngcore::ToString(err));
      return ptr;        
    }
    
    static void Free(Dev<T> * data)
    {
      cudaFree (data);   
    }


    
    
    T D2H() const
    {
      T res;
      cudaMemcpy (&res, &data, sizeof(T), cudaMemcpyDeviceToHost);
      return res;
    }

    void H2D (T val)

    {
      cudaMemcpy (this, &val, sizeof(T), cudaMemcpyHostToDevice);
    }

    void D2H (ngcore::FlatArray<T> hosta)
    {
      cudaMemcpy (hosta.Data(), &data, hosta.Size()*sizeof(T), cudaMemcpyDeviceToHost);
    }

    void H2D (ngcore::FlatArray<T> hosta)
    {
      cudaMemcpy (&data, hosta.Data(), hosta.Size()*sizeof(T), cudaMemcpyHostToDevice);
    }
    
    __host__ __device__ Dev<T> & operator= (T d2)
    {
#ifdef __CUDA_ARCH__      
      data = d2;
#else
      H2D(d2);
#endif
      return *this;
    }

    __host__ __device__ operator T() const
    {
#ifdef __CUDA_ARCH__
      return data;
#else
      return D2H();
#endif
    }
    
    __device__ const T& operator*() const { return data; }
    
    template <typename T2>
    __device__ auto & operator+= (T2 other) { data += other; return *this; }
    template <typename T2>
    __device__ auto & operator-= (T2 other) { data -= other; return *this; }
    template <typename T2>
    __device__ auto & operator*= (T2 other) { data *= other; return *this; }
  };


  template <typename T>
  Dev<T> * Host2Device (const T& val)
  {
    auto ptr = Dev<T>::Malloc(1);
    ptr->H2D(val);
    return ptr;
  }

  template <typename T>
  void Free(Dev<T> * dp)
  {
    Dev<T>::Free (dp);
  }  

  

}

#endif

