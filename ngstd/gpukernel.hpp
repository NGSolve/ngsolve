#ifndef FILE_GPUKERNEL
#define FILE_GPUKERNEL

/*********************************************************************/
/* File:   gpukernel.hpp                                             */
/* Date:   29. Aug. 2025                                             */
/*********************************************************************/

/*
  Common kernel syntax for cuda and metal.

  The prelude itself is gpukernel.h, embedded at build time as the
  string code_gpukernel (embed_source.cmake, like tinybla). Prepend it
  to the source given to Device::CompileSource, before code_tinybla. A kernel is then written once:

      KERNEL(saxpy, GLOBAL_IN(float,x), GLOBAL(float,y),
                    VALUE(float,a), VALUE(int,n))
      {
        int i = GLOBAL_ID_X;
        if (i < n) y[i] = a*x[i] + y[i];
      }

  Arguments are positional: the i-th declared argument is the i-th entry
  of the args list passed to Queue::Launch. Metal assigns buffer indices
  in declaration order, cuda takes them as kernel parameters.
*/

#include <string>

namespace ngs_gpu
{
  extern const std::string code_gpukernel;
}

#endif
