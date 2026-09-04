#ifndef FILE_GPU_BTDTB_HPP
#define FILE_GPU_BTDTB_HPP

/*********************************************************************/
/* File:   gpu_btdtb.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   3. Aug. 2026                                              */
/*********************************************************************/

#include <comp.hpp>
#include <gpuwrapper.hpp>
#include <gpukernel.hpp>
#include <tinybla.hpp>
#include <devicevector.hpp>
#include <map>
#include <mutex>

namespace ngcomp
{

  inline std::string Substitute(std::string src, const std::string &token,
                                const std::string &value)
  {
    for (size_t p = src.find(token); p != std::string::npos; p = src.find(token, p))
      src.replace(p, token.size(), value);
    return src;
  }


  template <uint N>
  constexpr uint RoundUp (uint i) { return N * ( (i+N-1) / N ); }
  template <uint N>
  constexpr uint RoundDown (uint i) { return N * ( i / N ); }
  
  constexpr uint RoundUp (uint i, uint N) { return N * ( (i+N-1) / N ); }
  constexpr uint RoundDown (uint i, uint N) { return N * ( i / N ); }


  template <typename REAL>
  class GPU_BTDTBMatrix : public BaseMatrix
  {
    int h, w;
    int ne, warps;
    int ngroups;

    shared_ptr<ngs_gpu::Library> library;
    shared_ptr<ngs_gpu::Kernel> kernel;
    shared_ptr<ngs_gpu::Queue> queue;

    ngs_gpu::TypedBuffer<int> buffer_dofx, buffer_dofy, buffer_domain;
    ngs_gpu::TypedBuffer<REAL> buffer_bmatx, buffer_bmaty;
    ngs_gpu::TypedBuffer<REAL> buffer_weights, buffer_geocoefs, buffer_bgeo, buffer_sgeo;
    
  public:
    GPU_BTDTBMatrix (const BaseMatrix& mat);


    
    
    AutoVector CreateRowVector() const override {
      return make_unique<DeviceVector<REAL>>(w, PreferredMemType());
    }
    AutoVector CreateColVector() const override {
      return make_unique<DeviceVector<REAL>>(h, PreferredMemType());
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    
  };


  template <typename REAL>
  GPU_BTDTBMatrix<REAL> :: GPU_BTDTBMatrix (const BaseMatrix& bmat)
  {
    using namespace ngs_gpu;
    h = bmat.Height();
    w = bmat.Width();

    
    auto* pmat = dynamic_cast<const MatrixFreeBTDTB*> (&bmat);
    if (!pmat)
      throw Exception("Expected MatrixFreeBTDTB matrix");

    auto device = ngs_gpu::GetDevice();
    if (!device)
      throw Exception("GPU_BTDTBMatrix: no gpu device");
    queue = device->DefaultQueue();

    auto [locdofsx, dimxref, nip] = pmat->Bx.Shape();
    auto [locdofsy, dimyref, nip_] = pmat->By.Shape();
    auto [numels,dimy,dimx,nipD] = pmat->D.Shape();
    auto [numels2,dimr,dims,nipJ] = pmat->Jacobi.Shape();

    int BS_ipts = pmat->opts.BS_ipts;
    ne = numels;
    if (dimr != dims)
      throw Exception("GPU_BTDTBMatrix: only square element Jacobians supported (volume elements)");

    auto NewSharedBuffer = [&] (size_t nreals)
    { return device->NewBuffer<REAL> (nreals, MemType::Shared); };
    
    /*
      The dofs of an element decompose into intervals of consecutive numbers
      (one per vertex/edge/face/cell), with the same structure on every
      element: one L2 interval, 5 for H(div), at most 15 for H1. We store one
      base per interval and element, the interval number and the offset within
      the interval are global constants. A space violating the structure just
      yields more intervals, never a wrong index.
    */
    auto SplitIntervals = [] (const Table<DofId> & dofs, int locdofs,
                              Array<int> & runof, Array<int> & offof)
    {
      if (locdofs == 0) return 0;
      Array<bool> newrun(locdofs);
      newrun = false;
      newrun[0] = true;
      for (size_t elnr : Range(dofs.Size()))
        for (int k = 1; k < locdofs; k++)
          if (dofs[elnr][k] != dofs[elnr][k-1]+1)
            newrun[k] = true;

      runof.SetSize(locdofs);
      offof.SetSize(locdofs);
      int run = -1, off = 0;
      for (int k = 0; k < locdofs; k++)
        {
          if (newrun[k]) { run++; off = 0; }
          runof[k] = run;
          offof[k] = off++;
        }
      return run+1;
    };

    auto IntervalBases = [&] (const Table<DofId> & dofs, int locdofs, int nruns,
                              FlatArray<int> runof, FlatArray<int> offof)
    {
      if (dofs.Size() != size_t(ne))
        throw Exception("GPU_BTDTBMatrix: dof table does not match number of elements");
      auto buffer = device->NewBuffer<int> (size_t(ne)*nruns, MemType::Shared);
      int * bases = buffer.HostData();
      for (size_t elnr : Range(dofs.Size()))
        for (int k = 0; k < locdofs; k++)
          if (offof[k] == 0)
            bases[elnr*nruns+runof[k]] = dofs[elnr][k];
      return buffer;
    };

    Array<int> runofx, offofx, runofy, offofy;
    int nrunsx = SplitIntervals (pmat->dofx, locdofsx, runofx, offofx);
    int nrunsy = SplitIntervals (pmat->dofy, locdofsy, runofy, offofy);
    if (locdofsx >= 65536 || locdofsy >= 65536 || nip >= 65536)
      throw Exception("GPU_BTDTBMatrix: 16-bit index tables need locdofs, nip < 65536");
    buffer_dofx = IntervalBases (pmat->dofx, locdofsx, nrunsx, runofx, offofx);
    buffer_dofy = IntervalBases (pmat->dofy, locdofsy, nrunsy, runofy, offofy);

    // for the remainder, store the compressed tensor in the last rows of buffer_bmatx
    buffer_bmatx = NewSharedBuffer(dimxref*(RoundUp(nip,BS_ipts)+BS_ipts)*RoundUp<8>(locdofsx));
    FlatTensor<3,REAL> dbmatx(dimxref,RoundUp(nip,BS_ipts),RoundUp<8>(locdofsx), buffer_bmatx.HostData());
    for (int j = 0; j < dimxref; j++)
      for (int i = 0; i < RoundUp(nip,BS_ipts); i++)
        for (int k = 0; k < RoundUp<8>(locdofsx); k++)
          dbmatx(j,i,k) = (i < nip && k < locdofsx) ? pmat->Bx(k,j,i) : 0;

    int bmatx_rem_rows = dimxref*RoundUp(nip,BS_ipts);
    int nip_rem = nip - RoundDown(nip, BS_ipts);
    FlatTensor<3,REAL> dbmatx_rem(dimxref,nip_rem,RoundUp<8>(locdofsx), buffer_bmatx.HostData() + bmatx_rem_rows*RoundUp<8>(locdofsx));
    dbmatx_rem = 0;
    for (int j = 0; j < dimxref; j++)
      for (int i = 0; i < nip_rem; i++)
        for (int k = 0; k < RoundUp<8>(locdofsx); k++)
          dbmatx_rem(j,i,k) = (k < locdofsx) ? pmat->Bx(k,j,i+RoundDown(nip,BS_ipts)) : 0;

    
    buffer_bmaty = NewSharedBuffer(dimyref*(RoundUp(nip,BS_ipts)+BS_ipts)*RoundUp<8>(locdofsy));
    FlatTensor<3,REAL> dbmaty(dimyref,RoundUp(nip,BS_ipts),RoundUp<8>(locdofsy), buffer_bmaty.HostData());
    for (int j = 0; j < dimyref; j++)
      for (int i = 0; i < RoundUp(nip,BS_ipts); i++)
        for (int k = 0; k < RoundUp<8>(locdofsy); k++)
          dbmaty(j,i,k) = (i<nip && k < locdofsy) ? pmat->By(k,j,i) : 0;

    int bmaty_rem_rows = dimyref*RoundUp(nip,BS_ipts);    
    FlatTensor<3,REAL> dbmaty_rem(dimyref,nip_rem,RoundUp<8>(locdofsy), buffer_bmaty.HostData() + bmaty_rem_rows*RoundUp<8>(locdofsy));
    dbmaty_rem = 0;
    for (int j = 0; j < dimyref; j++)
      for (int i = 0; i < nip_rem; i++)
        for (int k = 0; k < RoundUp<8>(locdofsy); k++)
          dbmaty_rem(j,i,k) = (k < locdofsy) ? pmat->By(k,j,i+RoundDown(nip,BS_ipts)) : 0;
    


    buffer_domain = device->template NewBuffer<int> (ne, MemType::Shared);
    for (int e = 0; e < ne; e++)
      buffer_domain.HostData()[e] = (e < pmat->domains.Size()) ? pmat->domains[e] : 0;

    buffer_weights = NewSharedBuffer(RoundUp<8>(nip));
    for (int i = 0; i < RoundUp<8>(nip); i++)
      buffer_weights.HostData()[i] = i < nip ? pmat->weights[i] : 0;

    /*
      geometry: per element the coefficients of the mapping in the L2 basis
      (geocoefs[el][coordinate][node], contiguous), and the gradients of that
      basis at the integration points (bgeo[refdir][ip][node]). The kernel
      evaluates F at every point from these, 9*geo_ndof fma per point.
    */
    auto [numels_g, geo_ndof, dimr_g] = pmat->geocoefs.Shape();
    auto [geo_ndof_b, dims_g, nip_g] = pmat->Bgeo.Shape();
    if (int(numels_g) != ne || int(dimr_g) != dimr || int(nip_g) != nip || geo_ndof_b != geo_ndof)
      throw Exception("GPU_BTDTBMatrix: geometry data does not match");
    /*
      curved classes evaluate F per ip block by the warp matmul into shared
      memory (9*BS_ipts*BS_els reals); that exceeds metal's 32 KB threadgroup
      memory for large blocks, then F is evaluated per point instead - as for
      straight elements, where 4 coefficients make the per-point evaluation
      cheaper than a separate per-element phase.
    */
    bool geo_staged = pmat->geo_order > 1 && BS_ipts <= 16;
    int geo_roundup = RoundUp<8>(geo_ndof);
    int geo_stride = geo_staged ? geo_roundup : int(geo_ndof);   // the matmul needs K padded to 8
    buffer_geocoefs = NewSharedBuffer(size_t(ne)*dimr*geo_stride);
    for (int e = 0; e < ne; e++)
      for (int i = 0; i < dimr; i++)
        for (int n = 0; n < geo_stride; n++)
          buffer_geocoefs.HostData()[(size_t(e)*dimr+i)*geo_stride+n] = (n < int(geo_ndof)) ? pmat->geocoefs(e,n,i) : 0;
    // gradients laid out like bmatx: [refdir][ip (padded to BS_ipts)][node], plus the remainder block
    buffer_bgeo = NewSharedBuffer(dims*(RoundUp(nip,BS_ipts)+BS_ipts)*geo_stride);
    FlatTensor<3,REAL> dbgeo(dims,RoundUp(nip,BS_ipts),geo_stride, buffer_bgeo.HostData());
    for (int j = 0; j < dims; j++)
      for (int i = 0; i < RoundUp(nip,BS_ipts); i++)
        for (int k = 0; k < geo_stride; k++)
          dbgeo(j,i,k) = (i < nip && k < int(geo_ndof)) ? pmat->Bgeo(k,j,i) : 0;
    int bgeo_rem_rows = dims*RoundUp(nip,BS_ipts);
    FlatTensor<3,REAL> dbgeo_rem(dims,nip_rem,geo_stride, buffer_bgeo.HostData() + bgeo_rem_rows*geo_stride);
    dbgeo_rem = 0;
    for (int j = 0; j < dims; j++)
      for (int i = 0; i < nip_rem; i++)
        for (int k = 0; k < geo_stride; k++)
          dbgeo_rem(j,i,k) = (k < int(geo_ndof)) ? pmat->Bgeo(k,j,i+RoundDown(nip,BS_ipts)) : 0;
    // basis values at the points, for the coordinates: sgeo[ip][node]
    buffer_sgeo = NewSharedBuffer(size_t(nip)*geo_ndof);
    for (int i = 0; i < nip; i++)
      for (int n = 0; n < int(geo_ndof); n++)
        buffer_sgeo.HostData()[size_t(i)*geo_ndof+n] = pmat->Sgeo(n,i);


    
    string code = code_gpukernel + code_tinybla +
      R"(
      using namespace tinybla;

      template <uint N>
      constexpr uint RoundUp (uint i) { return N * ( (i+N-1) / N ); }
      template <uint N>
      constexpr uint RoundDown (uint i) { return N * ( i / N ); }

      typedef $REAL Real;

      // dofnr -> (interval, offset within interval), same for all elements;
      // 16 bit: interval < 16 for H1, offset < locdofs
      CONSTANT_ARRAY(unsigned short, runof_x, $LOCDOFSX) = $RUNOFX;
      CONSTANT_ARRAY(unsigned short, offof_x, $LOCDOFSX) = $OFFOFX;
      CONSTANT_ARRAY(unsigned short, runof_y, $LOCDOFSY) = $RUNOFY;
      CONSTANT_ARRAY(unsigned short, offof_y, $LOCDOFSY) = $OFFOFY;
      $NREF_TABLE

      KERNEL(apply_btdtb,
             GLOBAL_IN(Real, x),
             $YARG,
             GLOBAL_IN(int,   basex),
             GLOBAL_IN(int,   basey),
             GLOBAL_IN(Real, bmatx),
             GLOBAL_IN(Real, bmaty),
             GLOBAL_IN(Real, weights),
             GLOBAL_IN(Real, geocoefs),
             GLOBAL_IN(Real, bgeo),
             GLOBAL_IN(Real, sgeo),
             GLOBAL_IN(int,  domain),
             VALUE(Real, s),
             VALUE(int, ne_))
      {
      uint tid = LOCAL_ID_X;
      uint bid  = GROUP_ID_X;
      uint bdim  = GROUP_SIZE_X;

      int warpIdx = tid/32;
      const uint ne = ne_;
      constexpr uint geo_ndof = $GEO_NDOF;
      constexpr uint geo_roundup = $GEO_ROUNDUP;   // node stride of the coefficients and gradients
      constexpr uint dimr = $DIMR;
      constexpr uint dims = $DIMS;
      constexpr uint nip = $NIP;
      constexpr uint bs_els = $BS_ELS;
      constexpr uint bs_ipts = $BS_IPTS;
      constexpr uint nip_padded = RoundUp<bs_ipts> (nip);   // component stride of bmatx/bmaty
      constexpr uint locdofsx = $LOCDOFSX;
      constexpr uint locdofsy = $LOCDOFSY;
      constexpr uint nrunsx = $NRUNSX;
      constexpr uint nrunsy = $NRUNSY;

      constexpr uint locdofsx_roundup = RoundUp<8> (locdofsx);
      constexpr uint locdofsy_roundup = RoundUp<8> (locdofsy);

      SHARED_2D(Real, elvecx, $BS_ELS, locdofsx_roundup);
      SHARED_2D(Real, elvecy, $BS_ELS, locdofsy_roundup);

      constexpr int MAXDIMREF = ($DIMXREF>$DIMYREF) ? $DIMXREF : $DIMYREF;
      SHARED_2D(Real, pointvalsref, MAXDIMREF*$BS_IPTS, $BS_ELS);

      // geometry coefficients of the batch [element][coordinate][node], and for
      // straight elements (constant gradients) F per element [element][coord][refdir]
      SHARED(Real, elgeo, $BS_ELS*dimr*geo_roundup);
      auto mat_bgeo = MakeBareMatrix<RowMajor>(bgeo, geo_roundup);
#if ($GEO_STAGED==1)
      // F at the points of the current ip block [coordinate*dims+refdir][ip][element]
      SHARED_2D(Real, Fvals, dimr*dims*$BS_IPTS, $BS_ELS);
      auto mat_Fvals = MakeBareMatrix<RowMajor>(Fvals);
#endif

      auto mat_elvecx = MakeBareMatrix<RowMajor>(elvecx); 
      auto mat_elvecy = MakeBareMatrix<RowMajor>(elvecy);
      auto mat_pointvalsref = MakeBareMatrix<RowMajor>(pointvalsref);

      auto mat_bmatx = MakeBareMatrix<RowMajor>(bmatx, locdofsx_roundup);
      auto mat_bmaty = MakeBareMatrix<RowMajor>(bmaty, locdofsy_roundup);


      // one group per element batch
      { uint baseelem = bid*$BS_ELS;
      if (baseelem >= ne) return;

      // zero the padding columns of elvecx, the valid ones are loaded below
      if constexpr (locdofsx_roundup > locdofsx)
        {
          constexpr uint npad = locdofsx_roundup - locdofsx;
          for (uint i = tid; i < $BS_ELS*npad; i += bdim)
            elvecx[i/npad][locdofsx + i%npad] = 0;
        }

      // load element vectors
      for (uint i = tid; i < $BS_ELS*locdofsx; i += bdim)
        {
           uint dofnr = i % locdofsx;
           uint locelnr = i / locdofsx;
           uint elnr = baseelem + locelnr;
           elvecx[locelnr][dofnr] = (elnr < ne)
             ? x[basex[elnr*nrunsx+runof_x[dofnr]] + offof_x[dofnr]] : 0;
        }

      // zero elvecy
      for (uint i = tid; i < $BS_ELS*locdofsy_roundup; i += bdim)
        elvecy[i/locdofsy_roundup][i%locdofsy_roundup] = 0;

#if ($ONLY_LOADSTORE==0)
      // geometry coefficients, contiguous per element ([coordinate][node] padded to geo_roundup)
      for (uint i = tid; i < $BS_ELS*dimr*geo_roundup; i += bdim)
        elgeo[i] = (baseelem + i/(dimr*geo_roundup) < ne) ? geocoefs[baseelem*dimr*geo_roundup + i] : 0;
#endif


      BARRIER();

#if ($ONLY_LOADSTORE==1)

      for (uint i = tid; i < $BS_ELS*locdofsy_roundup; i+= bdim)
        {
          uint c = i/$BS_ELS;
          uint r = i%$BS_ELS;
          if (c < locdofsx_roundup)
            elvecy[r][c] = elvecx[r][c];
        }
      BARRIER();
#else


      for (uint baseip = 0; baseip+$BS_IPTS <= $NIP; baseip += $BS_IPTS)
        {
          constexpr uint TS_IPTS = 8;
          constexpr uint TS_ELS = 8;
          constexpr uint BSTS_IPTS = $BS_IPTS/TS_IPTS;
          constexpr uint BSTS_ELS = $BS_ELS/TS_ELS;
          static_assert($BS_IPTS % TS_IPTS==0);

          // multiply with Bx
          for (uint warp = warpIdx; warp < $DIMXREF*BSTS_IPTS*BSTS_ELS; warp += $WARPS)
            {
              uint eltile = warp %  BSTS_ELS;
              uint ipcomp = warp /  BSTS_ELS;

              uint iptile = ipcomp % BSTS_IPTS;
              uint comp = ipcomp / BSTS_IPTS;
              uint ip = baseip + iptile*TS_IPTS;

              WarpMatrix<TS_IPTS,TS_ELS,Real> pointvalsrefxi = 0;
              auto mb = mat_elvecx.SubMatrix(TS_ELS*eltile, 0);
              auto ma = mat_bmatx.SubMatrix(ip+nip_padded*comp, 0);

              pointvalsrefxi.AddMM<8*$DOFX_TILES>(ma, mb.Transpose(), tid);
              pointvalsrefxi.Store(mat_pointvalsref.SubMatrix(TS_IPTS*iptile+bs_ipts*comp, TS_ELS*eltile), tid);
            }

#if ($GEO_STAGED==1)
          // F = geometry coefficients x gradients, one warp tile per (coordinate, refdir, ip tile, el tile)
          for (uint warp = warpIdx; warp < dimr*dims*BSTS_IPTS*BSTS_ELS; warp += $WARPS)
            {
              uint eltile = warp % BSTS_ELS;
              uint rest = warp / BSTS_ELS;
              uint iptile = rest % BSTS_IPTS;
              uint ig = rest / BSTS_IPTS;          // coordinate*dims + refdir
              uint ip = baseip + iptile*TS_IPTS;

              WarpMatrix<TS_IPTS,TS_ELS,Real> f = 0;
              auto mb = MakeBareMatrix<RowMajor>(elgeo + (ig/dims)*geo_roundup, dimr*geo_roundup).SubMatrix(TS_ELS*eltile, 0);
              auto ma = mat_bgeo.SubMatrix(ip+nip_padded*(ig%dims), 0);
              f.AddMM<8*$GEO_TILES>(ma, mb.Transpose(), tid);
              f.Store(mat_Fvals.SubMatrix(TS_IPTS*iptile+bs_ipts*ig, TS_ELS*eltile), tid);
            }
#endif

          BARRIER();
  
          // work on integration points
          for (uint ip = tid; ip < bs_ipts*bs_els ; ip += bdim)
             {
               uint locelnr = ip % $BS_ELS;
               uint locipnr = ip / $BS_ELS;
               uint elnr = baseelem + locelnr;

#if ($ONLY_LOADSTOREB==0)
               Vec<$DIMXREF,Real> xrefvals;
               Vec<$DIMYREF,Real> yrefvals;
               Vec<$DIMX,Real> xvals;
               Vec<$DIMY,Real> yvals;

               Mat<$DIMR,$DIMS,Real> F;
               for (int i = 0; i < $DIMR; i++)
                  for (int j = 0; j < $DIMS; j++)
                    {
#if ($GEO_STAGED==1)
                      F(i,j) = Fvals[(i*dims+j)*bs_ipts + locipnr][locelnr];
#else
                      Real sum = 0;
                      for (uint c = 0; c < geo_ndof; c++)
                        sum += elgeo[(locelnr*dimr+i)*geo_roundup+c] * bgeo[(j*nip_padded + baseip+locipnr)*geo_roundup+c];
                      F(i,j) = sum;
#endif
                    }
               Real J = Det(F);

               for (int comp = 0; comp < $DIMXREF; comp++)
                  xrefvals(comp) = pointvalsref[locipnr+comp*bs_ipts][locelnr];

              const int domain_index = domain[(baseelem + locelnr < ne) ? baseelem + locelnr : 0];

              $MEASURE

              // xrefvals -> xvals
              $TRANSFORMX;   

              $PHYSICS
              yvals = weights[baseip+locipnr] * meas * yvals;

              // yvals -> yrefvals
              $TRANSFORMY;          

               for (int comp = 0; comp < $DIMYREF; comp++)
                  pointvalsref[locipnr+comp*bs_ipts][locelnr] = yrefvals(comp);
#endif
            }


          BARRIER();

          // multiply with By.trans
          for (uint blocknr = warpIdx; blocknr < BSTS_ELS * $DOFY_TILES; blocknr += $WARPS)
            {
              int eltile = blocknr % BSTS_ELS;  // which els
              int ydoftile = blocknr / BSTS_ELS;  // which dofs

              WarpMatrix<TS_ELS,8,Real> sum(mat_elvecy.SubMatrix(TS_ELS*eltile, 8*ydoftile), tid);

              for (int comp = 0; comp < $DIMYREF; comp++)
                {
                   auto ma = mat_pointvalsref.SubMatrix(bs_ipts*comp, TS_ELS*eltile);
                   auto mb = mat_bmaty.SubMatrix(baseip+nip_padded*comp, 8*ydoftile);
                   sum.AddMM<$BS_IPTS> (ma.Transpose(), mb, tid);
                }

              sum.Store(mat_elvecy.SubMatrix(8*eltile, 8*ydoftile), tid);
           }

          BARRIER();

         } // end base_ip



      // ip remainder
     
      constexpr uint baseip = RoundDown<bs_ipts>(nip); 
      constexpr uint numips = nip-baseip;

      if constexpr (numips > 0) {

          // multiply with Bx
        for (uint blocknr = warpIdx; blocknr < RoundUp<8>(numips*$DIMXREF)/8 * $EL_TILES; blocknr += $WARPS)
           {
             uint eltile = blocknr % $EL_TILES;  // which els
             uint iptile = blocknr / $EL_TILES;  // which pts*comp tile

             WarpMatrix<8,8,Real> pointvalsrefxi = 0;

              // slower for convection -> now good!
              auto mb = mat_elvecx.SubMatrix(8*eltile, 0);
              auto ma = mat_bmatx.SubMatrix($BMATX_REM_ROWS+8*iptile, 0);
              pointvalsrefxi.AddMM<RoundUp<8>(locdofsx)>(ma, mb.Transpose(), tid);
              pointvalsrefxi.Store(mat_pointvalsref.SubMatrix(8*iptile,8*eltile), tid);
            }

#if ($GEO_STAGED==1)
          // geometry: per coordinate a padded block of rows (refdir*numips + ip)
          constexpr uint frows = RoundUp<8>(numips*dims);
          for (uint blocknr = warpIdx; blocknr < dimr * (frows/8) * $EL_TILES; blocknr += $WARPS)
            {
              uint eltile = blocknr % $EL_TILES;
              uint rest = blocknr / $EL_TILES;
              uint rowtile = rest % (frows/8);
              uint coord = rest / (frows/8);

              WarpMatrix<8,8,Real> f = 0;
              auto mb = MakeBareMatrix<RowMajor>(elgeo + coord*geo_roundup, dimr*geo_roundup).SubMatrix(8*eltile, 0);
              auto ma = mat_bgeo.SubMatrix($BGEO_REM_ROWS+8*rowtile, 0);
              f.AddMM<8*$GEO_TILES>(ma, mb.Transpose(), tid);
              f.Store(mat_Fvals.SubMatrix(coord*frows + 8*rowtile, 8*eltile), tid);
            }
#endif

          BARRIER();
  
          // work on integration points
          for (uint ip = tid; ip < numips*bs_els; ip += bdim)
             {
               uint locelnr = ip % $BS_ELS;
               uint locipnr = ip / $BS_ELS;
               uint elnr = baseelem + locelnr;

#if ($ONLY_LOADSTOREB==0)
               Vec<$DIMXREF,Real> xrefvals;
               Vec<$DIMYREF,Real> yrefvals;
               Vec<$DIMX,Real> xvals;
               Vec<$DIMY,Real> yvals;

               Mat<$DIMR,$DIMS,Real> F;
               for (uint i = 0; i < $DIMR; i++)
                  for (uint j = 0; j < $DIMS; j++)
                    {
#if ($GEO_STAGED==1)
                      F(i,j) = Fvals[i*frows + j*numips + locipnr][locelnr];
#else
                      Real sum = 0;
                      for (uint c = 0; c < geo_ndof; c++)
                        sum += elgeo[(locelnr*dimr+i)*geo_roundup+c] * bgeo[($BGEO_REM_ROWS + j*numips + locipnr)*geo_roundup+c];
                      F(i,j) = sum;
#endif
                    }
               Real J = Det(F);

               for (uint j = 0; j < $DIMXREF; j++)
                  xrefvals(j) = pointvalsref[locipnr+j*numips][locelnr];

              const int domain_index = domain[(baseelem + locelnr < ne) ? baseelem + locelnr : 0];

              $MEASURE

              // xrefvals -> xvals
              $TRANSFORMX;   

              $PHYSICS
              yvals = weights[baseip+locipnr] * meas * yvals;

              // yvals -> yrefvals
              $TRANSFORMY;          

               for (int j = 0; j < $DIMYREF; j++)
                  pointvalsref[locipnr+j*numips][locelnr] = yrefvals(j);

               for (int j = $DIMYREF*numips; j < $DIMYREF*$BS_IPTS; j++)
                  pointvalsref[j][locelnr] = 0;
#endif
            }


          BARRIER();

          // multiply with By.trans
          for (uint blocknr = warpIdx; blocknr < $EL_TILES * $DOFY_TILES; blocknr += $WARPS)
            {
              uint eltile = blocknr % $EL_TILES;  // which els
              uint ydoftile = blocknr / $EL_TILES;  // which dofs

              WarpMatrix<8,8,Real> sum(mat_elvecy.SubMatrix(8*eltile, 8*ydoftile), tid);

              for (uint ipcomp = 0; ipcomp < numips*$DIMYREF; ipcomp += 8)
                {
                   auto ma = mat_pointvalsref.SubMatrix(ipcomp, 8*eltile);
                   auto mb = mat_bmaty.SubMatrix($BMATY_REM_ROWS+ipcomp, 8*ydoftile);
                   sum.AddMM<8>(ma.Transpose(), mb, tid);
                }

              sum.Store(mat_elvecy.SubMatrix(8*eltile, 8*ydoftile), tid);
           }

          BARRIER();

      }  // if (numips_remainder > 0)

#endif

      // store vector
      for (uint i = tid; i < $BS_ELS*locdofsy; i += bdim)
        {
           uint dofnr = i % locdofsy;
           uint locelnr = i / locdofsy;

           uint elnr = baseelem + locelnr;
           if (elnr < ne)
             {
               uint dof = basey[elnr*nrunsy+runof_y[dofnr]] + offof_y[dofnr];
#if ($ATOMIC==1)
               ATOMIC_ADD(&y[dof], s*elvecy[locelnr][dofnr]);
#else
               y[dof] += s*elvecy[locelnr][dofnr];
#endif
             }
        }

      }
    }

)";



    string phys = "{ // THE PHYSICS (WIP)\n";
    phys += "int i = 0;\n";     // used in generated code, not needed here
    
    // for each test proxy we have one equation
    for (int k : Range(pmat -> test_proxies))
      {
        CoefficientFunction::T_DJC cache;
        auto diffcf = pmat->cf -> DiffJacobi (pmat->test_proxies[k], cache);
        auto compiledcf = Compile (diffcf, false);
        Code code = compiledcf->GenerateProgram(0, false);

        phys += "{  // equation " + ToString(k) +"\n";
        
        for (auto step : Range(compiledcf->Steps()))
          {
            auto stepcf = compiledcf->Steps()[step];
            
            if (auto proxycf = dynamic_cast<ProxyFunction*> (stepcf))
              {
                auto pos = pmat->trial_proxies.Pos(proxycf);
                if (pos != pmat->trial_proxies.ILLEGAL_POSITION)
                  {
                    phys += "auto values_"+ToString(step)+" = [&](int ip, int nr) { return xvals("+ToString(pmat->ranges_x[pos].First())+"+nr); };\n";
                    phys += "bool constexpr has_values_" + ToString(step) + " = true;\n";
                    
                    // Declare dummy comp_ variables to avoid compile errors (won't be used since has_values = true)
                    for(auto i : Range(proxycf->Dimension()))
                      phys += Var("comp", step,i, proxycf->Dimensions()).Declare("double", 0.0);
                  }
              }
          }
    
        phys += code.body;

        int steps = compiledcf->Steps().Size();
        IntRange rangey = pmat -> ranges_y[k];
        if (rangey.Size()==1)
          phys += "yvals("+ToString(rangey[0])+") = var_"+ToString(steps-1)+";\n";          
        else
          {
            auto dims = compiledcf->Dimensions();
            for (int l : Range(rangey))
              phys += "yvals("+ToString(rangey[l])+") = " + Var(steps-1, l, dims).S() + ";\n";
          }
        phys += "}\n";
      }
      
    phys += "} // END PHYSICS \n";

    phys = Substitute(phys, "double", "Real");
    code = Substitute(code, "$PHYSICS", phys);

    

    warps = pmat->opts.warps;

    // y is accumulated atomically only if elements may share dofs
    code = Substitute(code, "$YARG", pmat->opts.atomic
                      ? "GLOBAL_ATOMIC(Real, y)" : "GLOBAL(Real, y)");
    
    code = Substitute(code, "$NIP", ToString(nip)+" /*nip*/ ");
    code = Substitute(code, "$BS_ELS", ToString(pmat->opts.BS_els)+" /*BS_ELS*/ ");
    code = Substitute(code, "$BS_IPTS", ToString(pmat->opts.BS_ipts)+" /*BS_IPTS*/ ");
    code = Substitute(code, "$DIMR", ToString(dimr)+" /*dimr*/ ");
    code = Substitute(code, "$DIMS", ToString(dims)+" /*dims*/ ");
    code = Substitute(code, "$WARPS", ToString(warps)+" /*warps*/ ");

    code = Substitute(code, "$BMATX_REM_ROWS", ToString(bmatx_rem_rows)+" /*bmatx_rem_rows*/ ");      
    code = Substitute(code, "$BMATY_REM_ROWS", ToString(bmaty_rem_rows)+" /*bmaty_rem_rows*/ ");      
    code = Substitute(code, "$GEO_NDOF", ToString(geo_ndof));
    code = Substitute(code, "$GEO_ROUNDUP", ToString(geo_stride));
    code = Substitute(code, "$GEO_TILES", ToString(geo_roundup/8));
    code = Substitute(code, "$GEO_STAGED", ToString(geo_staged ? 1 : 0));
    code = Substitute(code, "$BGEO_REM_ROWS", ToString(bgeo_rem_rows)+" /*bgeo_rem_rows*/ ");
    
    code = Substitute(code, "$EL_TILES", ToString(RoundUp<8>(pmat->opts.BS_els)/8)+" /*EL_TILES*/ ");
    code = Substitute(code, "$DOFX_TILES", ToString(RoundUp<8>(locdofsx)/8) +" /* DOFX_TILES */ ");
    code = Substitute(code, "$DOFY_TILES", ToString(RoundUp<8>(locdofsy)/8) +" /* DOFY_TILES */ ");
    
    auto InitList = [] (FlatArray<int> a)
    {
      string s = "{";
      for (auto i : Range(a))
        { if (i) s += ","; s += ToString(a[i]); }
      return s+"}";
    };
    {
      // element-boundary integrals: weight * |Cof(F) nref_ip| instead of weight * det
      bool eb = pmat->facetnr.Size() > 0;
      string table, measure = "Real meas = J;";
      if (eb)
        {
          // reference normal per facet, facet of each integration point
          int nfacets = pmat->normals_ref.Height();
          Array<double> flat(nfacets*dimr);
          for (int f = 0; f < nfacets; f++)
            for (int d = 0; d < dimr; d++)
              flat[f*dimr+d] = pmat->normals_ref(f,d);
          string init = "{";
          for (auto i : Range(flat)) { if (i) init += ","; init += ToString(flat[i]); }
          table = "CONSTANT_ARRAY(Real, nref, " + ToString(nfacets*dimr) + ") = " + init + "};\n";
          table += "      CONSTANT_ARRAY(unsigned short, facetof, " + ToString(nip) + ") = " + InitList(pmat->facetnr) + ";";
          string nvec = "Vec<" + ToString(dimr) + ",Real>(";
          for (int d = 0; d < dimr; d++)
            nvec += (d ? "," : "") + string("nref[") + ToString(dimr) + "*facetof[baseip+locipnr]+" + ToString(d) + "]";
          nvec += ")";
          // measure |Cof(F) nref| and physical normal; the generated physics
          // accesses the normal as normals(i,k)
          string D = ToString(dimr);
          measure =
            "Vec<"+D+",Real> cofn = Cof(F)*" + nvec + ";\n"
            "Real meas = Norm(cofn);\n"
            "struct { Vec<"+D+",Real> nv; Real operator() (int, int k) const { return nv(k); } }\n"
            "  normals { ((J > 0) ? Real(1) : Real(-1)) / meas * cofn };\n";
        }
      // coordinates x(ip) = sum_node coefs * basis values, accessed as points(i,k)
      if (phys.find("points(") != string::npos)
        {
          string D = ToString(dimr);
          measure +=
            "Vec<"+D+",Real> pnt;\n"
            "for (uint k = 0; k < dimr; k++)\n"
            "  { Real sum = 0;\n"
            "    for (uint c = 0; c < geo_ndof; c++)\n"
            "      sum += elgeo[(locelnr*dimr+k)*geo_roundup+c] * sgeo[(baseip+locipnr)*geo_ndof+c];\n"
            "    pnt(k) = sum; }\n"
            "struct { Vec<"+D+",Real> p; Real operator() (int, int k) const { return p(k); } } points { pnt };\n";
        }
      code = Substitute(code, "$NREF_TABLE", table);
      code = Substitute(code, "$MEASURE", measure);
    }
    code = Substitute(code, "$NRUNSX", ToString(nrunsx));
    code = Substitute(code, "$NRUNSY", ToString(nrunsy));
    code = Substitute(code, "$RUNOFX", InitList(runofx));
    code = Substitute(code, "$OFFOFX", InitList(offofx));
    code = Substitute(code, "$RUNOFY", InitList(runofy));
    code = Substitute(code, "$OFFOFY", InitList(offofy));

    code = Substitute(code, "$LOCDOFSX", ToString(locdofsx));    
    code = Substitute(code, "$LOCDOFSY", ToString(locdofsy));
    code = Substitute(code, "$DIMXREF", ToString(dimxref));
    code = Substitute(code, "$DIMYREF", ToString(dimyref));
    code = Substitute(code, "$DIMX", ToString(dimx));
    code = Substitute(code, "$DIMY", ToString(dimy));

    code = Substitute(code, "$ONLY_LOADSTOREB", ToString(pmat->opts.only_loadstoreB));    
    code = Substitute(code, "$ONLY_LOADSTORE", ToString(pmat->opts.only_loadstore));    
    code = Substitute(code, "$ATOMIC", ToString(pmat->opts.atomic));    

    string transxcode;
    for (uint i = 0 ; i < pmat->diffopsx.Size(); i++)
      {
        IntRange rangexref = pmat->ranges_xref[i];
        IntRange rangex = pmat->ranges_x[i];
        transxcode += "{\n";
        transxcode += "Vec<" + ToString(rangex.Size()) + ",Real> res;\n";
        transxcode += pmat->diffopsx[i] -> GenerateTransformationCode("xrefvals.Range<"+ToString(rangexref.First())+","+ToString(rangexref.Next())+">()",
                                                                      "res", false);
        transxcode += "xvals.SetRange<" + ToString(rangex.First()) +"," + ToString(rangex.Next()) + ">(res);\n";
        transxcode += "}\n";
      }
    code = Substitute(code, "$TRANSFORMX", transxcode);    



    string transycode;
    for (uint i = 0 ; i < pmat->diffopsy.Size(); i++)
      {
        IntRange rangeyref = pmat->ranges_yref[i];
        IntRange rangey = pmat->ranges_y[i];
        transycode += "{\n";
        transycode += "Vec<" + ToString(rangeyref.Size()) + ",Real> res;\n";
        transycode += pmat->diffopsy[i] -> GenerateTransformationCode("yvals.Range<"+ToString(rangey.First())+","+ToString(rangey.Next())+">()",
                                                                      "res", true);
        transycode += "yrefvals.SetRange<" + ToString(rangeyref.First()) +"," + ToString(rangeyref.Next()) + ">(res);\n";
        transycode += "}\n";
      }
    code = Substitute(code, "$TRANSFORMY", transycode);    

    if constexpr (std::is_same_v<REAL,float>)
      code = Substitute(code, "$REAL", "float");
    else
      code = Substitute(code, "$REAL", "double");
    
    if (pmat->opts.write_GPU_kernel)
      ofstream(*pmat->opts.write_GPU_kernel) << code;

    /*
      the generated source depends on space, order and options only (ne is a
      kernel argument), so all element classes and all meshes of one kernel
      shape share one compiled library
    */
    bool cached = false;
    {
      static std::map<std::pair<Device*, string>, shared_ptr<Library>> cache;
      static std::mutex mutex;
      std::lock_guard<std::mutex> lock(mutex);
      auto & entry = cache[{device.get(), code}];
      if (!entry)
        entry = device->CompileSource (code);
      else
        cached = true;
      library = entry;
    }
    kernel = library->GetKernel ("apply_btdtb");

    /*
      one group per element batch: no grid-stride loop in the kernel, so no
      loop-carried state - the register count drops to ~45 for every space
      (measured 64..72 with the loop) and occupancy is no longer an issue.
      The finer granularity also balances the tail; the per-group dispatch
      cost only shows for cheap batches (L2 at BS_els=16, ~5% on cuda).
    */
    ngroups = (ne+pmat->opts.BS_els-1)/pmat->opts.BS_els;
    if (getenv("NGS_BTDTB_INFO"))
      cout << "apply_btdtb: " << (cached ? "cached " : "compiled ") << kernel->Info(warps*32) << " groups=" << ngroups
           << " groupsize=" << warps*32 << " intervals x/y=" << nrunsx << "/" << nrunsy
           << " locdofs x/y=" << locdofsx << "/" << locdofsy << endl;
  }


  template <typename REAL>
  void GPU_BTDTBMatrix<REAL> ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    using namespace ngs_gpu;
    static Timer tfull("GPU_BTDTBMatrix::MultAdd"); RegionTimer rfull(tfull);

    DeviceVectorWrapper<REAL> dvx(x);
    DeviceVectorWrapper<REAL> dvy(y);

    queue->Launch (*kernel, Dim3(ngroups), Dim3(warps*32),
                   { dvx.DevArgRO(), dvy.DevArgRW(),
                     buffer_dofx, buffer_dofy,
                     buffer_bmatx, buffer_bmaty,
                     buffer_weights, buffer_geocoefs, buffer_bgeo, buffer_sgeo,
                     buffer_domain, REAL(s), int(ne) });
  }
}


#endif
