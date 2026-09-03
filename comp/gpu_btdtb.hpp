#ifndef FILE_GPU_BTDTB_HPP
#define FILE_GPU_BTDTB_HPP

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

    shared_ptr<ngs_gpu::Buffer> buffer_dofx, buffer_dofy;
    shared_ptr<ngs_gpu::Buffer> buffer_bmatx, buffer_bmaty;
    shared_ptr<ngs_gpu::Buffer> buffer_weights, buffer_Jacobi, buffer_JacobiDets;
    
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
    { return device->NewBuffer (nreals*sizeof(REAL), MemType::Shared); };
    
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
      auto buffer = device->NewBuffer(size_t(ne)*nruns*sizeof(int), MemType::Shared);
      int * bases = buffer->HostData<int>();
      for (size_t elnr : Range(dofs.Size()))
        for (int k = 0; k < locdofs; k++)
          if (offof[k] == 0)
            bases[elnr*nruns+runof[k]] = dofs[elnr][k];
      return buffer;
    };

    Array<int> runofx, offofx, runofy, offofy;
    int nrunsx = SplitIntervals (pmat->dofx, locdofsx, runofx, offofx);
    int nrunsy = SplitIntervals (pmat->dofy, locdofsy, runofy, offofy);
    buffer_dofx = IntervalBases (pmat->dofx, locdofsx, nrunsx, runofx, offofx);
    buffer_dofy = IntervalBases (pmat->dofy, locdofsy, nrunsy, runofy, offofy);

    // for the remainder, store the compressed tensor in the last rows of buffer_bmatx
    buffer_bmatx = NewSharedBuffer(dimxref*(RoundUp(nip,BS_ipts)+BS_ipts)*RoundUp<8>(locdofsx));
    FlatTensor<3,REAL> dbmatx(dimxref,RoundUp(nip,BS_ipts),RoundUp<8>(locdofsx), buffer_bmatx->HostData<REAL>());
    for (int j = 0; j < dimxref; j++)
      for (int i = 0; i < RoundUp(nip,BS_ipts); i++)
        for (int k = 0; k < RoundUp<8>(locdofsx); k++)
          dbmatx(j,i,k) = (i < nip && k < locdofsx) ? pmat->Bx(k,j,i) : 0;

    int bmatx_rem_rows = dimxref*RoundUp(nip,BS_ipts);
    int nip_rem = nip - RoundDown(nip, BS_ipts);
    FlatTensor<3,REAL> dbmatx_rem(dimxref,nip_rem,RoundUp<8>(locdofsx), buffer_bmatx->HostData<REAL>() + bmatx_rem_rows*RoundUp<8>(locdofsx));
    dbmatx_rem = 0;
    for (int j = 0; j < dimxref; j++)
      for (int i = 0; i < nip_rem; i++)
        for (int k = 0; k < RoundUp<8>(locdofsx); k++)
          dbmatx_rem(j,i,k) = (k < locdofsx) ? pmat->Bx(k,j,i+RoundDown(nip,BS_ipts)) : 0;

    
    buffer_bmaty = NewSharedBuffer(dimyref*(RoundUp(nip,BS_ipts)+BS_ipts)*RoundUp<8>(locdofsy));
    FlatTensor<3,REAL> dbmaty(dimyref,RoundUp(nip,BS_ipts),RoundUp<8>(locdofsy), buffer_bmaty->HostData<REAL>());
    for (int j = 0; j < dimyref; j++)
      for (int i = 0; i < RoundUp(nip,BS_ipts); i++)
        for (int k = 0; k < RoundUp<8>(locdofsy); k++)
          dbmaty(j,i,k) = (i<nip && k < locdofsy) ? pmat->By(k,j,i) : 0;

    int bmaty_rem_rows = dimyref*RoundUp(nip,BS_ipts);    
    FlatTensor<3,REAL> dbmaty_rem(dimyref,nip_rem,RoundUp<8>(locdofsy), buffer_bmaty->HostData<REAL>() + bmaty_rem_rows*RoundUp<8>(locdofsy));
    dbmaty_rem = 0;
    for (int j = 0; j < dimyref; j++)
      for (int i = 0; i < nip_rem; i++)
        for (int k = 0; k < RoundUp<8>(locdofsy); k++)
          dbmaty_rem(j,i,k) = (k < locdofsy) ? pmat->By(k,j,i+RoundDown(nip,BS_ipts)) : 0;
    


    buffer_weights = NewSharedBuffer(RoundUp<8>(nip));
    for (int i = 0; i < RoundUp<8>(nip); i++)
      buffer_weights->HostData<REAL>()[i] = i < nip ? pmat->weights[i] : 0;
    
    buffer_JacobiDets = NewSharedBuffer(ne);
    for (int i = 0; i < ne; i++)
      buffer_JacobiDets->HostData<REAL>()[i] = Det(pmat->Jacobi(i,STAR,STAR,0));


    int ne8 = RoundUp<8>(ne);
    buffer_Jacobi = NewSharedBuffer(dimr*dims*ne8);
    for (int i = 0; i < ne; i++)
      for (int j = 0; j < dimr; j++)
        for (int k = 0; k < dims; k++)
          buffer_Jacobi->HostData<REAL>()[i+k*ne8+j*dims*ne8] = pmat->Jacobi(i,j,k,0);

    
    string code = code_gpukernel + code_tinybla +
      R"(
      using namespace tinybla;

      template <uint N>
      constexpr uint RoundUp (uint i) { return N * ( (i+N-1) / N ); }
      template <uint N>
      constexpr uint RoundDown (uint i) { return N * ( i / N ); }

      typedef $REAL real;

      // dofnr -> (interval, offset within interval), same for all elements
      CONSTANT_ARRAY(int, runof_x, $LOCDOFSX) = $RUNOFX;
      CONSTANT_ARRAY(int, offof_x, $LOCDOFSX) = $OFFOFX;
      CONSTANT_ARRAY(int, runof_y, $LOCDOFSY) = $RUNOFY;
      CONSTANT_ARRAY(int, offof_y, $LOCDOFSY) = $OFFOFY;

      KERNEL(apply_btdtb,
             GLOBAL_IN(real, x),
             $YARG,
             GLOBAL_IN(int,   basex),
             GLOBAL_IN(int,   basey),
             GLOBAL_IN(real, bmatx),
             GLOBAL_IN(real, bmaty),
             GLOBAL_IN(real, weights),
             GLOBAL_IN(real, Jacobi),
             GLOBAL_IN(real, JacobiDets),
             VALUE(real, s),
             VALUE(int, ne_))
      {
      uint tid = LOCAL_ID_X;
      uint bid  = GROUP_ID_X;
      uint bdim  = GROUP_SIZE_X;

      int warpIdx = tid/32;
      const uint ne = ne_;
      const uint ne8 = RoundUp<8>(ne);      // stride of the Jacobi buffer
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

      SHARED_2D(real, elvecx, $BS_ELS, locdofsx_roundup);
      SHARED_2D(real, elvecy, $BS_ELS, locdofsy_roundup);

      constexpr int MAXDIMREF = ($DIMXREF>$DIMYREF) ? $DIMXREF : $DIMYREF;
      SHARED_2D(real, pointvalsref, MAXDIMREF*$BS_IPTS, $BS_ELS);

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

              WarpMatrix<TS_IPTS,TS_ELS,real> pointvalsrefxi = 0;
              auto mb = mat_elvecx.SubMatrix(TS_ELS*eltile, 0);
              auto ma = mat_bmatx.SubMatrix(ip+nip_padded*comp, 0);

              pointvalsrefxi.AddMM<8*$DOFX_TILES>(ma, mb.Transpose(), tid);
              pointvalsrefxi.Store(mat_pointvalsref.SubMatrix(TS_IPTS*iptile+bs_ipts*comp, TS_ELS*eltile), tid);
            }

          BARRIER();
  
          // work on integration points
          for (uint ip = tid; ip < bs_ipts*bs_els ; ip += bdim)
             {
               uint locelnr = ip % $BS_ELS;
               uint locipnr = ip / $BS_ELS;
               uint elnr = baseelem + locelnr;

#if ($ONLY_LOADSTOREB==0)
               Vec<$DIMXREF,real> xrefvals;
               Vec<$DIMYREF,real> yrefvals;
               Vec<$DIMX,real> xvals;
               Vec<$DIMY,real> yvals;

               Mat<$DIMR,$DIMS,real> F;
               for (int i = 0; i < $DIMR; i++)
                  for (int j = 0; j < $DIMS; j++)
                     F(i,j) = Jacobi[elnr + ne8 * (j + $DIMS*i)];
               real J = (elnr < ne) ?  JacobiDets[elnr] : 0;

               for (int comp = 0; comp < $DIMXREF; comp++)
                  xrefvals(comp) = pointvalsref[locipnr+comp*bs_ipts][locelnr];

              // xrefvals -> xvals
              $TRANSFORMX;   

              $PHYSICS
              yvals = weights[baseip+locipnr] * JacobiDets[elnr] * yvals;

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

              WarpMatrix<TS_ELS,8,real> sum(mat_elvecy.SubMatrix(TS_ELS*eltile, 8*ydoftile), tid);

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

             WarpMatrix<8,8,real> pointvalsrefxi = 0;

              // slower for convection -> now good!
              auto mb = mat_elvecx.SubMatrix(8*eltile, 0);
              auto ma = mat_bmatx.SubMatrix($BMATX_REM_ROWS+8*iptile, 0);
              pointvalsrefxi.AddMM<RoundUp<8>(locdofsx)>(ma, mb.Transpose(), tid);
              pointvalsrefxi.Store(mat_pointvalsref.SubMatrix(8*iptile,8*eltile), tid);
            }

          BARRIER();
  
          // work on integration points
          for (uint ip = tid; ip < numips*bs_els; ip += bdim)
             {
               uint locelnr = ip % $BS_ELS;
               uint locipnr = ip / $BS_ELS;
               uint elnr = baseelem + locelnr;

#if ($ONLY_LOADSTOREB==0)
               Vec<$DIMXREF,real> xrefvals;
               Vec<$DIMYREF,real> yrefvals;
               Vec<$DIMX,real> xvals;
               Vec<$DIMY,real> yvals;

               Mat<$DIMR,$DIMS,real> F;
#pragma unroll
               for (uint i = 0; i < $DIMR; i++)
#pragma unroll
                  for (uint j = 0; j < $DIMS; j++)
                     F(i,j) = Jacobi[elnr + ne8 * (j + $DIMS*i)];
               real J = (elnr < ne) ?  JacobiDets[elnr] : 0;

               for (uint j = 0; j < $DIMXREF; j++)
                  xrefvals(j) = pointvalsref[locipnr+j*numips][locelnr];

              // xrefvals -> xvals
              $TRANSFORMX;   

              $PHYSICS
              yvals = weights[baseip+locipnr] * JacobiDets[elnr] * yvals;

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

              WarpMatrix<8,8,real> sum(mat_elvecy.SubMatrix(8*eltile, 8*ydoftile), tid);

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
          for (int l : Range(rangey))
            phys += "yvals("+ToString(rangey[l])+") = var_"+ToString(steps-1)+"_"+ToString(l)+";\n";
        phys += "}\n";
      }
      
    phys += "} // END PHYSICS \n";

    phys = Substitute(phys, "double", "real");
    code = Substitute(code, "$PHYSICS", phys);

    

    warps = pmat->opts.warps;

    // y is accumulated atomically only if elements may share dofs
    code = Substitute(code, "$YARG", pmat->opts.atomic
                      ? "GLOBAL_ATOMIC(real, y)" : "GLOBAL(real, y)");
    
    code = Substitute(code, "$NIP", ToString(nip)+" /*nip*/ ");
    code = Substitute(code, "$BS_ELS", ToString(pmat->opts.BS_els)+" /*BS_ELS*/ ");
    code = Substitute(code, "$BS_IPTS", ToString(pmat->opts.BS_ipts)+" /*BS_IPTS*/ ");
    code = Substitute(code, "$DIMR", ToString(dimr)+" /*dimr*/ ");
    code = Substitute(code, "$DIMS", ToString(dims)+" /*dims*/ ");
    code = Substitute(code, "$WARPS", ToString(warps)+" /*warps*/ ");

    code = Substitute(code, "$BMATX_REM_ROWS", ToString(bmatx_rem_rows)+" /*bmatx_rem_rows*/ ");      
    code = Substitute(code, "$BMATY_REM_ROWS", ToString(bmaty_rem_rows)+" /*bmaty_rem_rows*/ ");      
    
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
        transxcode += "Vec<" + ToString(rangex.Size()) + ",real> res;\n";
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
        transycode += "Vec<" + ToString(rangeyref.Size()) + ",real> res;\n";
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
                     buffer_weights, buffer_Jacobi, buffer_JacobiDets, REAL(s), int(ne) });
  }
}


#endif
