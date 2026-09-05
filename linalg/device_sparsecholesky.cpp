/*********************************************************************/
/* File:   device_sparsecholesky.cpp                                 */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   4. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_SPARSECHOLESKY_CPP

#include <la.hpp>
#include <gpukernel.hpp>
#include <climits>
#include <chrono>

namespace ngla
{
  using namespace ngs_gpu;

  namespace
  {
    string SubstituteScal (string src, const string & value)
    {
      const string token = "SCAL";
      for (size_t p = src.find(token); p != string::npos; p = src.find(token, p))
        src.replace (p, token.size(), value);
      return src;
    }

    // MicroTask::BT
    constexpr int L_BLOCK = 0, B_BLOCK = 1, LB_BLOCK = 2;

    /*
      The persistent solves: a group is GROUP_SIZE_X lanes (the simd
      width) times GROUP_SIZE_Y independent rows. Each row takes jobs
      from cnt, waits until incomingdep[job] has been counted down to
      zero by the jobs it depends on, does the work of the micro-task
      and decrements the counters of the dependent jobs.

      hy is passed twice: plainly for the ordinary loads and stores,
      and as an atomic buffer for the accumulations into rows that
      other jobs update at the same time.
    */
    const char * kernel_source =
    R"RAW(

      KERNEL(chol_reorder, GLOBAL_IN(SCAL,src), GLOBAL(SCAL,dst), GLOBAL_IN(int,order), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n && order[i] != -1) dst[order[i]] = src[i];
      }

      KERNEL(chol_reorder_add, GLOBAL_IN(SCAL,src), GLOBAL(SCAL,dst), GLOBAL_IN(int,order),
                               VALUE(SCAL,s), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n && order[i] != -1) dst[i] += s*src[order[i]];
      }

      KERNEL(chol_multdiag, GLOBAL_IN(SCAL,diag), GLOBAL(SCAL,vec), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n) vec[i] *= diag[i];
      }

      // working counters back to their initial values, job counters to zero
      KERNEL(chol_reset, GLOBAL_IN(int,dep0), GLOBAL(int,dep), GLOBAL_IN(int,dep0t), GLOBAL(int,dept),
                         GLOBAL(int,cnt), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n) { dep[i] = dep0[i]; dept[i] = dep0t[i]; }
        if (i < 2) cnt[i] = 0;
      }


      #define CHOL_TILE 32   // register tile of the L solves, at least the lane count

      // per-task descriptor (8 ints, host-side precomputed): type, rfirst, rnext,
      // base into rowindex2, ext_size, mfirst, mnext of the micro-block; loaded
      // before the wait so the dependent latency is off the critical path
      #define CHOL_TASK_LOAD                                                   \
        int type     = taskdesc[8*myjob];                                      \
        int rfirst   = taskdesc[8*myjob+1];                                    \
        int rnext    = taskdesc[8*myjob+2];                                    \
        int base     = taskdesc[8*myjob+3];                                    \
        int ext_size = taskdesc[8*myjob+4];                                    \
        int mfirst   = taskdesc[8*myjob+5];                                    \
        int mnext    = taskdesc[8*myjob+6];                                    \
        int rsize = rnext - rfirst;

    )RAW"
    // msvc: at most 16 KB per string literal, the pieces are concatenated
    R"RAW(
      KERNEL(chol_solveL, GLOBAL_IN(int,depidx), GLOBAL_IN(int,depdata), GLOBAL_ATOMIC(int,incomingdep),
                          GLOBAL(SCAL,hy), GLOBAL_ATOMIC(SCAL,hy_atomic), GLOBAL_ATOMIC(int,cnt),
                          GLOBAL_IN(int,taskdesc), GLOBAL_IN(int,rowindex2),
                          GLOBAL_IN(int,firstinrow), GLOBAL_IN(SCAL,lfact), VALUE(int,njobs))
      {
        SHARED(int, myjobs, 16);
        int lane = int(LOCAL_ID_X), row = int(LOCAL_ID_Y);
        int lanes = int(GROUP_SIZE_X);

        while (true)
          {
            if (lane == 0) myjobs[row] = ATOMIC_ADD(&cnt[0], 1);
            SIMD_BARRIER();
            int myjob = myjobs[row];
            if (myjob >= njobs) break;
            CHOL_TASK_LOAD

            if (lane == 0)
              while (ATOMIC_LOAD(&incomingdep[myjob]) != 0) { }
            SIMD_BARRIER();
            DEVICE_FENCE();

            // L part in tiles of `lanes` columns: lane holds x of its column,
            // the pivot is broadcast, then one trailing update of the rest
            // of the block. U[i,j] = lfact[firstinrow[i]-i-1+j]
            if (type == 0 || type == 2)     // L_BLOCK, LB_BLOCK
              for (int c0 = rfirst; c0 < rnext; c0 += lanes)
                {
                  int col = c0 + lane;
                  bool valid = col < rnext;
                  int nt = (rnext - c0 < lanes) ? rnext - c0 : lanes;
                  SCAL x = valid ? hy[col] : SCAL(0);
                  // the lane's column of the tile into registers, then the pivot chain
                  SCAL u[CHOL_TILE];
                  #pragma unroll
                  for (int i = 0; i < CHOL_TILE; i++)
                    u[i] = (i < nt && lane > i && valid) ? lfact[firstinrow[c0+i] - (c0+i) - 1 + col] : SCAL(0);
                  #pragma unroll
                  for (int i = 0; i < CHOL_TILE-1; i++)
                    if (i < nt-1)
                      {
                        SCAL xi = SIMD_BROADCAST(x, i);
                        x -= u[i] * xi;
                      }
                  if (valid) hy[col] = x;
                  for (int j0 = c0 + nt; j0 < rnext; j0 += lanes)
                    {
                      int j = j0 + lane;
                      bool jv = j < rnext;
                      SCAL v = jv ? hy[j] : SCAL(0);
                      for (int i = 0; i < nt; i++)
                        {
                          SCAL xi = SIMD_BROADCAST(x, i);
                          if (jv) v -= lfact[firstinrow[c0+i] - (c0+i) - 1 + j] * xi;
                        }
                      if (jv) hy[j] = v;
                    }
                  SIMD_BARRIER();   // the next tile reads what other lanes wrote
                }
            SIMD_BARRIER();   // B part reads the hy written by the L part

            if ((type == 1 || type == 2) && ext_size != 0)    // B_BLOCK, LB_BLOCK
              for (int j = lane; j < mnext-mfirst; j += lanes)
                {
                  SCAL temp = 0;
                  for (int i = rfirst; i < rnext; i++)
                    temp += lfact[firstinrow[i] + rnext-i-1 + mfirst + j] * hy[i];
                  ATOMIC_ADD(&hy_atomic[rowindex2[base + mfirst + j]], -temp);
                }

            SIMD_BARRIER();
            DEVICE_FENCE();
            // all lanes release the dependents: the root of the elimination
            // tree has thousands of them in the backward solve
            for (int k = depidx[myjob] + lane; k < depidx[myjob+1]; k += lanes)
              ATOMIC_ADD(&incomingdep[depdata[k]], -1);
          }
      }


    )RAW"
    // msvc: at most 16 KB per string literal, the pieces are concatenated
    R"RAW(
      KERNEL(chol_solveLT, GLOBAL_IN(int,depidx), GLOBAL_IN(int,depdata), GLOBAL_ATOMIC(int,incomingdep),
                           GLOBAL(SCAL,hy), GLOBAL_ATOMIC(SCAL,hy_atomic), GLOBAL_ATOMIC(int,cnt),
                           GLOBAL_IN(int,taskdesc), GLOBAL_IN(int,rowindex2),
                           GLOBAL_IN(int,firstinrow), GLOBAL_IN(SCAL,lfact), VALUE(int,njobs))
      {
        SHARED(int, myjobs, 16);
        int lane = int(LOCAL_ID_X), row = int(LOCAL_ID_Y);
        int lanes = int(GROUP_SIZE_X);

        while (true)
          {
            if (lane == 0) myjobs[row] = ATOMIC_ADD(&cnt[1], 1);
            SIMD_BARRIER();
            int myjob = njobs-1 - myjobs[row];
            if (myjob < 0) break;
            CHOL_TASK_LOAD

            if (lane == 0)
              while (ATOMIC_LOAD(&incomingdep[myjob]) != 0) { }
            SIMD_BARRIER();
            DEVICE_FENCE();

            if ((type == 1 || type == 2) && ext_size != 0)    // B_BLOCK, LB_BLOCK
              for (int i = rfirst + lane; i < rnext; i += lanes)
                {
                  int first = firstinrow[i] + rnext-i-1 + mfirst;
                  SCAL val = 0;
                  for (int j = 0; j < mnext-mfirst; j++)
                    val += lfact[first + j] * hy[rowindex2[base + mfirst + j]];
                  ATOMIC_ADD(&hy_atomic[i], -val);
                }
            SIMD_BARRIER();   // L part reads the hy updated by the B part

            // L part in tiles from the last one down: lane holds x of its row,
            // column j of the tile is broadcast, then the rows before the tile
            // get the tile's contribution
            if (type == 0 || type == 2)     // L_BLOCK, LB_BLOCK
              for (int c0 = rfirst + ((rsize-1)/lanes)*lanes; c0 >= rfirst; c0 -= lanes)
                {
                  int row = c0 + lane;
                  bool valid = row < rnext;
                  int nt = (rnext - c0 < lanes) ? rnext - c0 : lanes;
                  int off = valid ? firstinrow[row] - row - 1 : 0;
                  SCAL x = valid ? hy[row] : SCAL(0);
                  SCAL u[CHOL_TILE];
                  #pragma unroll
                  for (int j = 0; j < CHOL_TILE; j++)
                    u[j] = (j < nt && lane < j && valid) ? lfact[off + c0 + j] : SCAL(0);
                  #pragma unroll
                  for (int j = CHOL_TILE-1; j > 0; j--)
                    if (j < nt)
                      {
                        SCAL xj = SIMD_BROADCAST(x, j);
                        x -= u[j] * xj;
                      }
                  if (valid) hy[row] = x;
                  for (int i0 = rfirst; i0 < c0; i0 += lanes)
                    {
                      int i = i0 + lane;
                      bool iv = i < c0;
                      int offi = iv ? firstinrow[i] - i - 1 : 0;
                      SCAL v = iv ? hy[i] : SCAL(0);
                      for (int j = 0; j < nt; j++)
                        {
                          SCAL xj = SIMD_BROADCAST(x, j);
                          if (iv) v -= lfact[offi + c0 + j] * xj;
                        }
                      if (iv) hy[i] = v;
                    }
                  SIMD_BARRIER();
                }

            SIMD_BARRIER();
            DEVICE_FENCE();
            // all lanes release the dependents: the root of the elimination
            // tree has thousands of them in the backward solve
            for (int k = depidx[myjob] + lane; k < depidx[myjob+1]; k += lanes)
              ATOMIC_ADD(&incomingdep[depdata[k]], -1);
          }
      }


    )RAW"
    // msvc: at most 16 KB per string literal, the pieces are concatenated
    R"RAW(
      /*
        The top of the elimination tree, one threadgroup per block: warp 0
        solves the block's L part, the warps then split the B micro-blocks,
        group barriers in between, one release per block. Block descriptor:
        rfirst, rnext, base, ext_size, nb.
      */
      #define CHOL_BLK_LOAD                                                    \
        int rfirst   = blkdesc[8*myblk];                                       \
        int rnext    = blkdesc[8*myblk+1];                                     \
        int base     = blkdesc[8*myblk+2];                                     \
        int ext_size = blkdesc[8*myblk+3];                                     \
        int nb       = blkdesc[8*myblk+4];                                     \
        int rsize = rnext - rfirst;

      KERNEL(chol_solveL_blk, GLOBAL_IN(int,depidx), GLOBAL_IN(int,depdata), GLOBAL_ATOMIC(int,incomingdep),
                              GLOBAL(SCAL,hy), GLOBAL_ATOMIC(SCAL,hy_atomic), GLOBAL_ATOMIC(int,cnt),
                              GLOBAL_IN(int,blkdesc), GLOBAL_IN(int,rowindex2),
                              GLOBAL_IN(int,firstinrow), GLOBAL_IN(SCAL,lfact), VALUE(int,nblk))
      {
        SHARED(int, myblks, 1);
        int lane = int(LOCAL_ID_X), row = int(LOCAL_ID_Y);
        int lanes = int(GROUP_SIZE_X), rows = int(GROUP_SIZE_Y);

        while (true)
          {
            if (lane == 0 && row == 0) myblks[0] = ATOMIC_ADD(&cnt[0], 1);
            DEVICE_BARRIER();
            int myblk = myblks[0];
            DEVICE_BARRIER();
            if (myblk >= nblk) break;
            CHOL_BLK_LOAD

            if (lane == 0 && row == 0)
              while (ATOMIC_LOAD(&incomingdep[myblk]) != 0) { }
            DEVICE_BARRIER();
            DEVICE_FENCE();

            if (row == 0)
              {
                // L part in tiles of `lanes` columns: lane holds x of its column,
                // the pivot is broadcast, then one trailing update of the rest
                // of the block. U[i,j] = lfact[firstinrow[i]-i-1+j]
                if (true)
                  for (int c0 = rfirst; c0 < rnext; c0 += lanes)
                    {
                      int col = c0 + lane;
                      bool valid = col < rnext;
                      int nt = (rnext - c0 < lanes) ? rnext - c0 : lanes;
                      SCAL x = valid ? hy[col] : SCAL(0);
                      // the lane's column of the tile into registers, then the pivot chain
                      SCAL u[CHOL_TILE];
                      #pragma unroll
                      for (int i = 0; i < CHOL_TILE; i++)
                        u[i] = (i < nt && lane > i && valid) ? lfact[firstinrow[c0+i] - (c0+i) - 1 + col] : SCAL(0);
                      #pragma unroll
                      for (int i = 0; i < CHOL_TILE-1; i++)
                        if (i < nt-1)
                          {
                            SCAL xi = SIMD_BROADCAST(x, i);
                            x -= u[i] * xi;
                          }
                      if (valid) hy[col] = x;
                      for (int j0 = c0 + nt; j0 < rnext; j0 += lanes)
                        {
                          int j = j0 + lane;
                          bool jv = j < rnext;
                          SCAL v = jv ? hy[j] : SCAL(0);
                          for (int i = 0; i < nt; i++)
                            {
                              SCAL xi = SIMD_BROADCAST(x, i);
                              if (jv) v -= lfact[firstinrow[c0+i] - (c0+i) - 1 + j] * xi;
                            }
                          if (jv) hy[j] = v;
                        }
                      SIMD_BARRIER();   // the next tile reads what other lanes wrote
                    }
              }
            DEVICE_BARRIER();   // B parts read the hy written by the L part

            for (int mb = row; mb < nb; mb += rows)
              {
                int mfirst = int((long)ext_size * mb / nb);
                int mnext  = int((long)ext_size * (mb+1) / nb);
                if (ext_size != 0)
                  for (int j = lane; j < mnext-mfirst; j += lanes)
                    {
                      SCAL temp = 0;
                      for (int i = rfirst; i < rnext; i++)
                        temp += lfact[firstinrow[i] + rnext-i-1 + mfirst + j] * hy[i];
                      ATOMIC_ADD(&hy_atomic[rowindex2[base + mfirst + j]], -temp);
                    }
              }

            DEVICE_BARRIER();
            if (row == 0)
              for (int k = depidx[myblk] + lane; k < depidx[myblk+1]; k += lanes)
                ATOMIC_ADD(&incomingdep[depdata[k]], -1);
          }
      }

      KERNEL(chol_solveLT_blk, GLOBAL_IN(int,depidx), GLOBAL_IN(int,depdata), GLOBAL_ATOMIC(int,incomingdep),
                               GLOBAL(SCAL,hy), GLOBAL_ATOMIC(SCAL,hy_atomic), GLOBAL_ATOMIC(int,cnt),
                               GLOBAL_IN(int,blkdesc), GLOBAL_IN(int,rowindex2),
                               GLOBAL_IN(int,firstinrow), GLOBAL_IN(SCAL,lfact), VALUE(int,nblk))
      {
        SHARED(int, myblks, 1);
        int lane = int(LOCAL_ID_X), row = int(LOCAL_ID_Y);
        int lanes = int(GROUP_SIZE_X), rows = int(GROUP_SIZE_Y);

        while (true)
          {
            if (lane == 0 && row == 0) myblks[0] = ATOMIC_ADD(&cnt[1], 1);
            DEVICE_BARRIER();
            int myblk = nblk-1 - myblks[0];
            DEVICE_BARRIER();
            if (myblk < 0) break;
            CHOL_BLK_LOAD

            if (lane == 0 && row == 0)
              while (ATOMIC_LOAD(&incomingdep[myblk]) != 0) { }
            DEVICE_BARRIER();
            DEVICE_FENCE();

            for (int mb = row; mb < nb; mb += rows)
              {
                int mfirst = int((long)ext_size * mb / nb);
                int mnext  = int((long)ext_size * (mb+1) / nb);
                if (ext_size != 0)
                  for (int i = rfirst + lane; i < rnext; i += lanes)
                    {
                      int first = firstinrow[i] + rnext-i-1 + mfirst;
                      SCAL val = 0;
                      for (int j = 0; j < mnext-mfirst; j++)
                        val += lfact[first + j] * hy[rowindex2[base + mfirst + j]];
                      ATOMIC_ADD(&hy_atomic[i], -val);
                    }
              }
            DEVICE_BARRIER();   // L part reads the hy updated by the B parts

            if (row == 0)
              {
                // L part in tiles from the last one down: lane holds x of its row,
                // column j of the tile is broadcast, then the rows before the tile
                // get the tile's contribution
                if (true)
                  for (int c0 = rfirst + ((rsize-1)/lanes)*lanes; c0 >= rfirst; c0 -= lanes)
                    {
                      int row = c0 + lane;
                      bool valid = row < rnext;
                      int nt = (rnext - c0 < lanes) ? rnext - c0 : lanes;
                      int off = valid ? firstinrow[row] - row - 1 : 0;
                      SCAL x = valid ? hy[row] : SCAL(0);
                      SCAL u[CHOL_TILE];
                      #pragma unroll
                      for (int j = 0; j < CHOL_TILE; j++)
                        u[j] = (j < nt && lane < j && valid) ? lfact[off + c0 + j] : SCAL(0);
                      #pragma unroll
                      for (int j = CHOL_TILE-1; j > 0; j--)
                        if (j < nt)
                          {
                            SCAL xj = SIMD_BROADCAST(x, j);
                            x -= u[j] * xj;
                          }
                      if (valid) hy[row] = x;
                      for (int i0 = rfirst; i0 < c0; i0 += lanes)
                        {
                          int i = i0 + lane;
                          bool iv = i < c0;
                          int offi = iv ? firstinrow[i] - i - 1 : 0;
                          SCAL v = iv ? hy[i] : SCAL(0);
                          for (int j = 0; j < nt; j++)
                            {
                              SCAL xj = SIMD_BROADCAST(x, j);
                              if (iv) v -= lfact[offi + c0 + j] * xj;
                            }
                          if (iv) hy[i] = v;
                        }
                      SIMD_BARRIER();
                    }
              }

            DEVICE_BARRIER();
            if (row == 0)
              for (int k = depidx[myblk] + lane; k < depidx[myblk+1]; k += lanes)
                ATOMIC_ADD(&incomingdep[depdata[k]], -1);
          }
      }

    )RAW";


    template <typename T>
    class DeviceCholeskyKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> reorder, reorder_add, multdiag, reset, solveL, solveLT, solveL_blk, solveLT_blk;
      unsigned groupsize;      // elementwise kernels
      unsigned lanes;          // x-size of the solve groups, at most the simd width
      unsigned rows_L = 8, rows_LT = 8;   // y-size: rows of lanes per group (rows in flight = groups*rows)
      unsigned solvegroups;    // persistent groups
      long topblocks = -1;     // blocks at the top of the tree run one per threadgroup; -1: time both (NGS_GPU_CHOL_TOP)
      int maxbs = 32;          // device blocking: columns per block, ext dofs per micro-block
      int maxmubs = 32;        // (NGS_GPU_CHOL_MAXBS, NGS_GPU_CHOL_MAXMUBS)
      unsigned rows_blk = 8;   // warps per group in the block kernels

      DeviceCholeskyKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        if constexpr (is_same_v<T,double>)
          if (!device->HasFloat64())
            throw Exception("DeviceSparseCholesky<double> on "+device->Name()+
                            ", which has no fp64 - use DeviceSparseCholesky<float>");

        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         SubstituteScal (kernel_source, scal));
        reorder     = library->GetKernel ("chol_reorder");
        reorder_add = library->GetKernel ("chol_reorder_add");
        multdiag    = library->GetKernel ("chol_multdiag");
        reset       = library->GetKernel ("chol_reset");
        solveL      = library->GetKernel ("chol_solveL");
        solveLT     = library->GetKernel ("chol_solveLT");
        solveL_blk  = library->GetKernel ("chol_solveL_blk");
        solveLT_blk = library->GetKernel ("chol_solveLT_blk");
        queue = device->DefaultQueue();

        bool gpu = device->SimdWidth() > 1;
        groupsize = gpu ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
        // the host backend runs groups one after another, so all jobs
        // must be taken by the rows of a single group
        lanes = gpu ? device->SimdWidth() : 8;
        if (lanes > 32) throw Exception("DeviceSparseCholesky: simd width exceeds the register tile");
        // few persistent groups per core: the spinning rows compete with the
        // working ones for memory traffic (Metal: 4/core best, 512 up to 10x slower)
        solvegroups = gpu ? 4*device->ComputeUnits() : 1;
        // overrides for experiments: NGS_GPU_CHOL_L=32x8 NGS_GPU_CHOL_LT=32x8 NGS_GPU_CHOL_GROUPS=80
        auto geom = [&] (const char * name, unsigned & rows)
        {
          if (const char * e = getenv(name))
            {
              unsigned x = 0, y = 0;
              if (sscanf (e, "%ux%u", &x, &y) == 2 && x > 0 && y > 0)
                {
                  // a row of lanes must be exactly one simd group: narrower rows
                  // share a barrier with a neighbour that may spin on them (deadlock)
                  if (gpu && x != device->SimdWidth())
                    throw Exception(string(name)+": lanes must equal the simd width");
                  lanes = x; rows = y;
                }
            }
        };
        geom ("NGS_GPU_CHOL_L", rows_L);
        geom ("NGS_GPU_CHOL_LT", rows_LT);
        if (getenv("NGS_GPU_CHOL_GROUPS")) solvegroups = atoi (getenv("NGS_GPU_CHOL_GROUPS"));
        if (getenv("NGS_GPU_CHOL_TOP")) topblocks = atoi (getenv("NGS_GPU_CHOL_TOP"));
        if (getenv("NGS_GPU_CHOL_MAXBS")) maxbs = max (1, atoi (getenv("NGS_GPU_CHOL_MAXBS")));
        if (getenv("NGS_GPU_CHOL_MAXMUBS")) maxmubs = max (1, atoi (getenv("NGS_GPU_CHOL_MAXMUBS")));
        if (getenv("NGS_GPU_CHOL_TOPROWS")) rows_blk = atoi (getenv("NGS_GPU_CHOL_TOPROWS"));
        cout << IM(7) << "chol solve groups = " << solvegroups << " (" << device->ComputeUnits() << " compute units)" << endl
             << IM(7) << "chol_solveL:  " << solveL->Info(lanes*rows_L) << endl
             << IM(7) << "chol_solveLT: " << solveLT->Info(lanes*rows_LT) << endl;
      }

      void LaunchElementwise (const shared_ptr<Kernel> & kernel, size_t n,
                              std::initializer_list<KernelArg> args) const
      {
        if (n == 0) return;
        unsigned groups = (n + groupsize-1) / groupsize;
        queue->Launch (*kernel, Dim3(groups), Dim3(groupsize), args);
      }

      static const DeviceCholeskyKernels & Get()
      {
        static mutex mtx;
        static shared_ptr<DeviceCholeskyKernels> cached;

        auto dev = GetGpuDevice();
        auto lock = lock_guard<mutex>(mtx);
        if (!cached || cached->device != dev)
          cached = make_shared<DeviceCholeskyKernels> (dev);
        return *cached;
      }
    };


    template <typename T>
    TypedBuffer<T> Upload (Device & device, FlatArray<T> data)
    {
      auto buf = device.NewBuffer<T> (max<size_t>(data.Size(),1), MemType::Device);
      if (data.Size()) buf.H2D (data.Data(), data.Size());
      return buf;
    }

    // dependency table as csr
    void TableToCSR (FlatTable<int> table, Array<int> & index, Array<int> & data)
    {
      index.SetSize (table.Size()+1);
      data.SetSize (table.AsArray().Size());
      index[0] = 0;
      for (size_t i = 0; i < table.Size(); i++)
        index[i+1] = index[i] + table[i].Size();
      data = table.AsArray();
    }
  }



  template <typename T>
  template <typename TM>
  DeviceSparseCholesky<T> :: DeviceSparseCholesky (const SparseCholeskyTM<TM> & mat)
  {
    height = mat.Height();
    width = mat.Width();
    nused = mat.GetNUsed();

    const auto & kern = DeviceCholeskyKernels<T>::Get();
    device = kern.device;
    queue = kern.queue;
    memtype = PreferredMemType();

    auto firstinrow = mat.GetFirstInRow();
    auto firstinrow_ri = mat.GetFirstInRowRI();
    auto rowindex2 = mat.GetRowIndex2();

    /*
      The device has its own blocking of the factor: the host blocks are
      supernodes (or pieces of them), cut here to at most `maxbs` columns,
      and the ext dofs of a block are split into micro-blocks of `maxmubs`.
      Block dependencies come from the row indices, the micro-task graph
      is L -> own B micro-blocks -> L of the dependent blocks, as on the
      host, but with sizes tuned for the kernels (a warp per micro-task).
    */
    auto hostblocks = mat.GetBlocks();
    Array<int> blocks;
    for (size_t b = 0; b + 1 < hostblocks.Size(); b++)
      for (int i = hostblocks[b]; i < hostblocks[b+1]; i += kern.maxbs)
        blocks.Append (i);
    if (nused > 0) blocks.Append (nused);
    int nblocks = blocks.Size() > 0 ? blocks.Size()-1 : 0;

    Array<int> block_of_dof(nused);
    for (int b = 0; b < nblocks; b++)
      block_of_dof[Range(blocks[b], blocks[b+1])] = b;

    Table<int> block_dependency;
    {
      TableCreator<int> creator(nblocks);
      for ( ; !creator.Done(); creator++)
        for (int b = 0; b < nblocks; b++)
          {
            Array<int> deps;
            for (int i = blocks[b]; i < blocks[b+1]; i++)
              for (int j : rowindex2.Range(firstinrow_ri[i], firstinrow_ri[i+1]))
                if (block_of_dof[j] != b) deps.Append (block_of_dof[j]);
            QuickSort (deps);
            for (size_t k = 0; k < deps.Size(); k++)
              if (k == 0 || deps[k] != deps[k-1]) creator.Add (b, deps[k]);
          }
      block_dependency = creator.MoveTable();
    }

    struct DevTask { int blocknr, type, bblock, nbblocks; };
    Array<DevTask> tasks;
    Array<int> first_task(nblocks+1);
    for (int b = 0; b < nblocks; b++)
      {
        first_task[b] = tasks.Size();
        int rfirst = blocks[b], rsize = blocks[b+1]-rfirst;
        long ext_size = long(firstinrow[rfirst+1] - firstinrow[rfirst]) - rsize + 1;
        int nb = int((ext_size + kern.maxmubs-1) / kern.maxmubs);
        if (nb <= 1)
          tasks.Append ({b, LB_BLOCK, 0, 1});
        else
          {
            tasks.Append ({b, L_BLOCK, 0, 0});
            for (int j = 0; j < nb; j++) tasks.Append ({b, B_BLOCK, j, nb});
          }
      }
    first_task[nblocks] = tasks.Size();
    njobs = tasks.Size();

    Table<int> dep, dept;
    {
      TableCreator<int> creator(njobs), creator_trans(njobs);
      for ( ; !creator.Done(); creator++, creator_trans++)
        for (int b = 0; b < nblocks; b++)
          {
            int l = first_task[b];
            if (first_task[b+1] == l+1)
              for (int o : block_dependency[b])
                { creator.Add (l, first_task[o]); creator_trans.Add (first_task[o], l); }
            else
              for (int t = l+1; t < first_task[b+1]; t++)
                {
                  creator.Add (l, t); creator_trans.Add (t, l);
                  for (int o : block_dependency[b])
                    { creator.Add (t, first_task[o]); creator_trans.Add (first_task[o], t); }
                }
          }
      dep = creator.MoveTable();
      dept = creator_trans.MoveTable();
    }
    Array<int> depidx, depdata, deptidx, deptdata;
    TableToCSR (dep, depidx, depdata);
    TableToCSR (dept, deptidx, deptdata);

    Array<int> taskdesc(8*njobs);
    for (size_t i = 0; i < njobs; i++)
      {
        const auto & t = tasks[i];
        int rfirst = blocks[t.blocknr], rnext = blocks[t.blocknr+1], rsize = rnext-rfirst;
        long ext_size = long(firstinrow[rfirst+1] - firstinrow[rfirst]) - rsize + 1;
        int * d = &taskdesc[8*i];
        d[0] = t.type;
        d[1] = rfirst;
        d[2] = rnext;
        d[3] = int(firstinrow_ri[rfirst]) + rsize - 1;
        d[4] = int(ext_size);
        d[5] = t.nbblocks ? int(ext_size * t.bblock / t.nbblocks) : 0;
        d[6] = t.nbblocks ? int(ext_size * (t.bblock+1) / t.nbblocks) : 0;
        d[7] = 0;
      }

    // the top `topblocks` blocks run in the block kernels after (forward)
    // or before (backward) the task kernels: cut the job list at the first
    // task of the first top block; about eight blocks per group, so the
    // leaves stay in the task kernels; decided by timing below
    long cand = kern.topblocks >= 0 ? kern.topblocks : 8L * kern.solvegroups;
    int blk0 = max<int> (0, nblocks - int(min<long>(cand, nblocks)));
    ncut = first_task[blk0];
    ntopblocks = nblocks - blk0;

    // the forward solve waits for every predecessor, the backward one for
    // every successor, counted among the jobs below the cut
    Array<int> incoming(njobs), incoming_trans(njobs), incoming_all(njobs), incoming_trans_all(njobs);
    incoming = 0; incoming_trans = 0; incoming_all = 0; incoming_trans_all = 0;
    for (size_t i = 0; i < njobs; i++)
      for (int d : dep[i])
        {
          if (d <= int(i)) throw Exception("DeviceSparseCholesky: dependency graph must be directional");
          incoming_all[d]++; incoming_trans_all[i]++;
          if (d < int(ncut)) { incoming[d]++; incoming_trans[i]++; }
        }

    // block dependencies among the top blocks, from the task edges
    Array<int> bdepidx(ntopblocks+1), bdepdata, bdeptidx(ntopblocks+1), bdeptdata, bincoming(ntopblocks), bincoming_t(ntopblocks);
    Array<int> blkdesc(8*max<size_t>(ntopblocks,1));
    {
      Array<Array<int>> bdep(ntopblocks), bdept(ntopblocks);
      for (size_t i = ncut; i < njobs; i++)
        for (int d : dep[i])
          {
            int bi = tasks[i].blocknr - blk0, bd = tasks[d].blocknr - blk0;
            if (bi != bd) { bdep[bi].Append (bd); bdept[bd].Append (bi); }
          }
      bdepidx[0] = 0; bdeptidx[0] = 0;
      bincoming = 0; bincoming_t = 0;
      for (int b = 0; b < ntopblocks; b++)
        {
          for (auto * r : { &bdep[b], &bdept[b] })
            {
              QuickSort (*r);
              int n = 0;
              for (size_t k = 0; k < r->Size(); k++)
                if (k == 0 || (*r)[k] != (*r)[k-1]) (*r)[n++] = (*r)[k];
              r->SetSize (n);
            }
          for (int c : bdep[b]) { bdepdata.Append (c); bincoming[c]++; }
          for (int c : bdept[b]) { bdeptdata.Append (c); bincoming_t[c]++; }
          bdepidx[b+1] = bdepdata.Size();
          bdeptidx[b+1] = bdeptdata.Size();

          int rfirst = blocks[blk0+b], rnext = blocks[blk0+b+1], rsize = rnext-rfirst;
          long ext_size = long(firstinrow[rfirst+1] - firstinrow[rfirst]) - rsize + 1;
          int nb = 0;
          for (int t = first_task[blk0+b]; t < first_task[blk0+b+1]; t++)
            nb = max (nb, tasks[t].nbblocks);
          int * d = &blkdesc[8*b];
          d[0] = rfirst; d[1] = rnext; d[2] = int(firstinrow_ri[rfirst]) + rsize - 1;
          d[3] = int(ext_size); d[4] = max(nb,1); d[5] = d[6] = d[7] = 0;
        }
    }

    // 64-bit offsets to int32
    if (firstinrow[nused] > size_t(INT_MAX) || firstinrow_ri[nused] > size_t(INT_MAX))
      throw Exception("DeviceSparseCholesky: factor too large for 32-bit offsets");
    Array<int> hfirstinrow(nused+1);
    for (size_t i = 0; i <= nused; i++)
      hfirstinrow[i] = int(firstinrow[i]);

    auto lfact = mat.GetLFact();
    auto diag = mat.GetDiag();
    Array<T> hlfact(lfact.Size()), hdiag(diag.Size());
    for (size_t i = 0; i < lfact.Size(); i++) hlfact[i] = T(lfact[i]);
    for (size_t i = 0; i < diag.Size(); i++) hdiag[i] = T(diag[i]);

    // depth of the dependency graph: the latency floor of the persistent solve;
    // the same path weighted by block columns counts the serialized column
    // steps of the L parts
    size_t levels = 0, maxblock = 0, colsteps = 0;
    {
      Array<int> level(njobs), cols(njobs);
      for (size_t i = 0; i < njobs; i++)
        {
          int bs = blocks[tasks[i].blocknr+1] - blocks[tasks[i].blocknr];
          maxblock = max<size_t> (maxblock, bs);
          cols[i] = (tasks[i].type == B_BLOCK) ? 1 : bs;
        }
      level = 1;
      for (size_t i = 0; i < njobs; i++)
        for (int d : dep[i])
          {
            level[d] = max (level[d], level[i]+1);
            cols[d] = max (cols[d], cols[i] + (tasks[d].type == B_BLOCK ? 1 : blocks[tasks[d].blocknr+1]-blocks[tasks[d].blocknr]));
          }
      for (int l : level) levels = max<size_t> (levels, l);
      for (int c : cols) colsteps = max<size_t> (colsteps, c);
    }
    cout << IM(7) << "DeviceSparseCholesky<" << (is_same_v<T,double> ? "double" : "float")
         << "> height = " << height << ", nused = " << nused << ", lfact = " << lfact.Size()
         << " (" << lfact.Size()*sizeof(T)/1e6 << " MB), blocks = " << nblocks << " (host " << (hostblocks.Size() ? hostblocks.Size()-1 : 0)
         << "), microtasks = " << njobs
         << ", dependency levels = " << levels << ", max block = " << maxblock
         << ", critical path column steps = " << colsteps << endl
         << IM(7) << "top blocks in block kernels: " << ntopblocks << " of " << nblocks
         << ", jobs below the cut: " << ncut << " of " << njobs << endl;

    dev_taskdesc = Upload<int> (*device, taskdesc);
    dev_depidx = Upload<int> (*device, depidx);
    dev_depdata = Upload<int> (*device, depdata);
    dev_depidx_trans = Upload<int> (*device, deptidx);
    dev_depdata_trans = Upload<int> (*device, deptdata);
    dev_incomingdep0 = Upload<int> (*device, incoming);
    dev_incomingdep0_trans = Upload<int> (*device, incoming_trans);
    dev_incomingdep0_all = Upload<int> (*device, incoming_all);
    dev_incomingdep0_trans_all = Upload<int> (*device, incoming_trans_all);
    dev_incomingdep = device->NewBuffer<int> (max<size_t>(njobs,1), MemType::Device);
    dev_incomingdep_trans = device->NewBuffer<int> (max<size_t>(njobs,1), MemType::Device);
    dev_cnt = device->NewBuffer<int> (2, MemType::Device);
    dev_blkdesc = Upload<int> (*device, blkdesc);
    dev_bdepidx = Upload<int> (*device, bdepidx);
    dev_bdepdata = Upload<int> (*device, bdepdata);
    dev_bdepidx_trans = Upload<int> (*device, bdeptidx);
    dev_bdepdata_trans = Upload<int> (*device, bdeptdata);
    dev_bincoming0 = Upload<int> (*device, bincoming);
    dev_bincoming0_trans = Upload<int> (*device, bincoming_t);
    dev_bincoming = device->NewBuffer<int> (max<size_t>(ntopblocks,1), MemType::Device);
    dev_bincoming_trans = device->NewBuffer<int> (max<size_t>(ntopblocks,1), MemType::Device);
    dev_bcnt = device->NewBuffer<int> (2, MemType::Device);

    dev_rowindex2 = Upload<int> (*device, mat.GetRowIndex2());
    dev_firstinrow = Upload<int> (*device, hfirstinrow);
    dev_lfact = Upload<T> (*device, hlfact);
    dev_diag = Upload<T> (*device, hdiag);
    dev_order = Upload<int> (*device, mat.GetOrder());
    dev_hx = device->NewBuffer<T> (max<size_t>(nused,1), MemType::Device);

    // with or without the block kernels: time both on this matrix
    usetop = ntopblocks > 0;
    if (usetop && kern.topblocks < 0 && ncut > 0)
      {
        DeviceVector<T> x(height, memtype), y(height, memtype);
        x = 1.0;
        auto timeit = [&] (bool top)
        {
          usetop = top;
          y = 0.0; MultAdd (1, x, y); queue->Finish();
          auto t0 = std::chrono::steady_clock::now();
          for (int k = 0; k < 3; k++) MultAdd (1, x, y);
          queue->Finish();
          return std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        };
        double t_top = timeit(true), t_tasks = timeit(false);
        usetop = t_top < t_tasks;
        cout << IM(7) << "DeviceSparseCholesky: block kernels " << t_top*1e3/3 << " ms, task kernels "
             << t_tasks*1e3/3 << " ms per solve -> " << (usetop ? "block kernels" : "task kernels") << endl;
      }
  }


  template <typename T>
  void DeviceSparseCholesky<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceSparseCholesky::MultAdd"); RegionTimer reg(t);
    if (x.Size() != width || y.Size() != height)
      throw Exception("DeviceSparseCholesky::MultAdd - size mismatch");
    if (nused == 0) return;

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    const auto & kern = DeviceCholeskyKernels<T>::Get();

    kern.LaunchElementwise (kern.reorder, height,
                            { ux.DevArgRO(), KernelArg(dev_hx), KernelArg(dev_order), KernelArg(int(height)) });
    size_t cut = usetop ? ncut : njobs, ntop = usetop ? ntopblocks : 0;
    kern.LaunchElementwise (kern.reset, max<size_t>(njobs,2),
                            { KernelArg(usetop ? dev_incomingdep0 : dev_incomingdep0_all), KernelArg(dev_incomingdep),
                              KernelArg(usetop ? dev_incomingdep0_trans : dev_incomingdep0_trans_all), KernelArg(dev_incomingdep_trans),
                              KernelArg(dev_cnt), KernelArg(int(njobs)) });
    if (ntop)
      kern.LaunchElementwise (kern.reset, max<size_t>(ntopblocks,2),
                              { KernelArg(dev_bincoming0), KernelArg(dev_bincoming),
                                KernelArg(dev_bincoming0_trans), KernelArg(dev_bincoming_trans),
                                KernelArg(dev_bcnt), KernelArg(int(ntopblocks)) });

    if (cut)
      queue->Launch (*kern.solveL, Dim3(kern.solvegroups), Dim3(kern.lanes, kern.rows_L),
                     { KernelArg(dev_depidx), KernelArg(dev_depdata), KernelArg(dev_incomingdep),
                       KernelArg(dev_hx), KernelArg(dev_hx), KernelArg(dev_cnt),
                       KernelArg(dev_taskdesc), KernelArg(dev_rowindex2),
                       KernelArg(dev_firstinrow), KernelArg(dev_lfact), KernelArg(int(cut)) });
    if (ntop)
      queue->Launch (*kern.solveL_blk, Dim3(kern.solvegroups), Dim3(kern.lanes, kern.rows_blk),
                     { KernelArg(dev_bdepidx), KernelArg(dev_bdepdata), KernelArg(dev_bincoming),
                       KernelArg(dev_hx), KernelArg(dev_hx), KernelArg(dev_bcnt),
                       KernelArg(dev_blkdesc), KernelArg(dev_rowindex2),
                       KernelArg(dev_firstinrow), KernelArg(dev_lfact), KernelArg(int(ntopblocks)) });

    kern.LaunchElementwise (kern.multdiag, nused,
                            { KernelArg(dev_diag), KernelArg(dev_hx), KernelArg(int(nused)) });

    if (ntop)
      queue->Launch (*kern.solveLT_blk, Dim3(kern.solvegroups), Dim3(kern.lanes, kern.rows_blk),
                     { KernelArg(dev_bdepidx_trans), KernelArg(dev_bdepdata_trans), KernelArg(dev_bincoming_trans),
                       KernelArg(dev_hx), KernelArg(dev_hx), KernelArg(dev_bcnt),
                       KernelArg(dev_blkdesc), KernelArg(dev_rowindex2),
                       KernelArg(dev_firstinrow), KernelArg(dev_lfact), KernelArg(int(ntopblocks)) });
    if (cut)
      queue->Launch (*kern.solveLT, Dim3(kern.solvegroups), Dim3(kern.lanes, kern.rows_LT),
                     { KernelArg(dev_depidx_trans), KernelArg(dev_depdata_trans), KernelArg(dev_incomingdep_trans),
                       KernelArg(dev_hx), KernelArg(dev_hx), KernelArg(dev_cnt),
                       KernelArg(dev_taskdesc), KernelArg(dev_rowindex2),
                       KernelArg(dev_firstinrow), KernelArg(dev_lfact), KernelArg(int(cut)) });

    kern.LaunchElementwise (kern.reorder_add, height,
                            { KernelArg(dev_hx), uy.DevArgRW(), KernelArg(dev_order),
                              KernelArg(T(s)), KernelArg(int(height)) });
  }


  template <typename T>
  AutoVector DeviceSparseCholesky<T> :: CreateRowVector () const
  {
    return make_unique<DeviceVector<T>> (width, memtype);
  }

  template <typename T>
  AutoVector DeviceSparseCholesky<T> :: CreateColVector () const
  {
    return make_unique<DeviceVector<T>> (height, memtype);
  }

  template <typename T>
  BaseMatrix::OperatorInfo DeviceSparseCholesky<T> :: GetOperatorInfo () const
  {
    return { string("DeviceSparseCholesky<") + (is_same_v<T,double> ? "double" : "float") + ">",
             height, width };
  }

  template <typename T>
  ostream & DeviceSparseCholesky<T> :: Print (ostream & ost) const
  {
    ost << "DeviceSparseCholesky<" << (is_same_v<T,double> ? "double" : "float")
        << ">, height = " << height << ", nused = " << nused
        << ", microtasks = " << njobs << ", on " << device->Name() << endl;
    return ost;
  }


  template class DeviceSparseCholesky<double>;
  template class DeviceSparseCholesky<float>;
  template DeviceSparseCholesky<double>::DeviceSparseCholesky (const SparseCholeskyTM<double> &);
  template DeviceSparseCholesky<double>::DeviceSparseCholesky (const SparseCholeskyTM<float> &);
  template DeviceSparseCholesky<float>::DeviceSparseCholesky (const SparseCholeskyTM<double> &);
  template DeviceSparseCholesky<float>::DeviceSparseCholesky (const SparseCholeskyTM<float> &);
}
