#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>



#include "metal_btdtb.hpp"
#include "tinybla.hpp"

namespace ngsmetal
{

  static std::string Substitute(std::string src, const std::string &token,
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
  
  MetalBTDTBMatrix :: MetalBTDTBMatrix (const BaseMatrix& bmat)
  {
    h = bmat.Height();
    w = bmat.Width();

    
    auto* pmat = dynamic_cast<const MatrixFreeBTDTB*> (&bmat);
    if (!pmat)
      throw Exception("Expected MatrixFreeBTDTB matrix");


    auto [locdofsx, dimxref, nip] = pmat->Bx.Shape();
    auto [locdofsy, dimyref, nip_] = pmat->By.Shape();
    auto [numels,dimy,dimx,nipD] = pmat->D.Shape();
    auto [numels2,dimr,dims,nipJ] = pmat->Jacobi.Shape();

    // cout << "nip = " << nip << endl;
    BS_els = pmat->opts.BS_els;
    BS_ipts = pmat->opts.BS_ipts;  
    ne = numels;
    
    FlatArray<int> hdofx = pmat->dofx.AsArray();
    buffer_dofx = GetDevice()->newBuffer(hdofx.Size()*sizeof(int), MTL::ResourceStorageModeShared);
    FlatArray<int> ddofx(hdofx.Size(), (int*)buffer_dofx->contents());
    ddofx = hdofx;

    FlatArray<int> hdofy = pmat->dofy.AsArray();
    buffer_dofy = GetDevice()->newBuffer(hdofy.Size()*sizeof(int), MTL::ResourceStorageModeShared);
    FlatArray<int> ddofy(hdofy.Size(), (int*)buffer_dofy->contents());
    ddofy = hdofy;

    // for the remainder, store the compressed tensor in the last rows of buffer_bmatx
    buffer_bmatx = GetDevice()->newBuffer(dimxref*(RoundUp(nip,BS_ipts)+BS_ipts)*RoundUp<8>(locdofsx)*sizeof(float), MTL::ResourceStorageModeShared);
    FlatTensor<3,float> dbmatx(dimxref,RoundUp(nip,BS_ipts),RoundUp<8>(locdofsx), (float*)buffer_bmatx->contents());
    for (int j = 0; j < dimxref; j++)
      for (int i = 0; i < RoundUp(nip,BS_ipts); i++)
        for (int k = 0; k < RoundUp<8>(locdofsx); k++)
          dbmatx(j,i,k) = (i < nip && k < locdofsx) ? pmat->Bx(k,j,i) : 0;

    int bmatx_rem_rows = dimxref*RoundUp(nip,BS_ipts);
    int nip_rem = nip - RoundDown(nip, BS_ipts);
    FlatTensor<3,float> dbmatx_rem(dimxref,nip_rem,RoundUp<8>(locdofsx), ((float*)buffer_bmatx->contents()) + bmatx_rem_rows*RoundUp<8>(locdofsx));
    dbmatx_rem = 0;
    for (int j = 0; j < dimxref; j++)
      for (int i = 0; i < nip_rem; i++)
        for (int k = 0; k < RoundUp<8>(locdofsx); k++)
          dbmatx_rem(j,i,k) = (k < locdofsx) ? pmat->Bx(k,j,i+RoundDown(nip,BS_ipts)) : 0;

    
    buffer_bmaty = GetDevice()->newBuffer(dimyref*(RoundUp(nip,BS_ipts)+BS_ipts)*RoundUp<8>(locdofsy)*sizeof(float), MTL::ResourceStorageModeShared);
    FlatTensor<3,float> dbmaty(dimyref,RoundUp(nip,BS_ipts),RoundUp<8>(locdofsy), (float*)buffer_bmaty->contents());
    for (int j = 0; j < dimyref; j++)
      for (int i = 0; i < RoundUp(nip,BS_ipts); i++)
        for (int k = 0; k < RoundUp<8>(locdofsy); k++)
          dbmaty(j,i,k) = (i<nip && k < locdofsy) ? pmat->By(k,j,i) : 0;

    int bmaty_rem_rows = dimyref*RoundUp(nip,BS_ipts);    
    FlatTensor<3,float> dbmaty_rem(dimyref,nip_rem,RoundUp<8>(locdofsy), ((float*)buffer_bmaty->contents()) + bmaty_rem_rows*RoundUp<8>(locdofsy));
    dbmaty_rem = 0;
    for (int j = 0; j < dimyref; j++)
      for (int i = 0; i < nip_rem; i++)
        for (int k = 0; k < RoundUp<8>(locdofsy); k++)
          dbmaty_rem(j,i,k) = (k < locdofsy) ? pmat->By(k,j,i+RoundDown(nip,BS_ipts)) : 0;
    


    buffer_weights = GetDevice()->newBuffer(RoundUp<8>(nip)*sizeof(float), MTL::ResourceStorageModeShared);
    for (int i = 0; i < RoundUp<8>(nip); i++)
      ((float*)buffer_weights->contents())[i] = i < nip ? pmat->weights[i] : 0;
    
    buffer_JacobiDets = GetDevice()->newBuffer(ne*sizeof(float), MTL::ResourceStorageModeShared);    
    for (int i = 0; i < ne; i++)
      ((float*)buffer_JacobiDets->contents())[i] = Det(pmat->Jacobi(i,STAR,STAR,0));


    int ne8 = RoundUp<8>(ne);
    buffer_Jacobi = GetDevice()->newBuffer(dimr*dims*ne8*sizeof(float), MTL::ResourceStorageModeShared);    
    for (int i = 0; i < ne; i++)
      for (int j = 0; j < dimr; j++)
        for (int k = 0; k < dims; k++)
          ((float*)buffer_Jacobi->contents())[i+k*ne8+j*dims*ne8] = pmat->Jacobi(i,j,k,0);

    auto debug = GetDevice()->newBuffer(1000*sizeof(float), MTL::ResourceStorageModeShared);
    FlatVector<float> debugvec(1000, (float*)debug->contents());
    debugvec = 0.;

    
    string code = code_tinybla +
      R"(
      #include <metal_stdlib>
      using namespace metal;
      using namespace tinybla;

      template <uint N>
      constexpr uint RoundUp (uint i) { return N * ( (i+N-1) / N ); }
      template <uint N>
      constexpr uint RoundDown (uint i) { return N * ( i / N ); }

      template <typename T>
      auto MyMin (T a, T b) { return (a < b) ? a : b; }


      kernel void apply_btdtb(
      device const float*  x            [[buffer(0)]],
#if ($ATOMIC==1)
      device atomic_float* y            [[buffer(1)]],
#else
      device float* y            [[buffer(1)]],
#endif
      device const int*    dofx         [[buffer(2)]],
      device const int*    dofy         [[buffer(3)]],
      device const float* bmatx         [[buffer(4)]],
      device const float* bmaty         [[buffer(5)]],
      device const float* weights       [[buffer(6)]],
      device const float* Jacobi        [[buffer(7)]],
      device const float* JacobiDets    [[buffer(8)]],
      device float *    debug           [[buffer(9)]],

      uint   blockIdx           [[threadgroup_position_in_grid]],
      uint   threadIdx          [[thread_position_in_threadgroup]],
      uint   blockDim           [[threads_per_threadgroup]],
      uint   gridDim            [[threadgroups_per_grid]]
      )
      {
      int warpIdx = threadIdx/32;
      constexpr uint ne = $NE;
      constexpr uint ne8 = RoundUp<8>(ne);
      constexpr uint nip = $NIP;
      constexpr uint nip8 = RoundUp<8> (nip); 
      constexpr uint bs_els = $BS_ELS;
      constexpr uint bs_ipts = $BS_IPTS;
      constexpr uint locdofsx = $LOCDOFSX;
      constexpr uint locdofsy = $LOCDOFSY;

      constexpr uint locdofsx_roundup = RoundUp<8> (locdofsx);
      constexpr uint locdofsy_roundup = RoundUp<8> (locdofsy);

      alignas(16) threadgroup float elvecx[$BS_ELS][locdofsx_roundup];
      alignas(16) threadgroup float elvecy[$BS_ELS][locdofsy_roundup];

      constexpr int MAXDIMREF = ($DIMXREF>$DIMYREF) ? $DIMXREF : $DIMYREF;
      alignas(16) threadgroup float pointvalsref[MAXDIMREF*$BS_IPTS][$BS_ELS];

      auto mat_elvecx = MakeBareMatrix<RowMajor>(elvecx); 
      auto mat_elvecy = MakeBareMatrix<RowMajor>(elvecy);
      auto mat_pointvalsref = MakeBareMatrix<RowMajor>(pointvalsref);

      auto mat_bmatx = MakeBareMatrix<RowMajor>(bmatx, locdofsx_roundup);
      auto mat_bmaty = MakeBareMatrix<RowMajor>(bmaty, locdofsy_roundup);


      // zero elvecx (overhead would be enough
      for (uint i = threadIdx; i < $BS_ELS*locdofsx_roundup; i+= blockDim)
        {
          uint c = i%locdofsx_roundup;
          uint r = i/locdofsx_roundup;
          elvecx[r][c] = 0;
        }

      threadgroup_barrier(mem_flags::mem_threadgroup);

      for (uint baseelem = blockIdx*$BS_ELS; baseelem < $NE; baseelem += gridDim*$BS_ELS) { 

      threadgroup_barrier(mem_flags::mem_threadgroup);

      // load element vectors
      for (uint i = threadIdx; i < $BS_ELS*locdofsx; i += blockDim)
        {
           uint dofnr = i % locdofsx;
           uint locelnr = i / locdofsx;
           uint elnr = baseelem + locelnr;
           elvecx[locelnr][dofnr] = (elnr < ne) ? x[dofx[elnr*locdofsx+dofnr]] : 0;
        }

      // zero elvecy
      for (uint i = threadIdx; i < $BS_ELS*locdofsy_roundup; i+= blockDim)
        {
          uint c = i%locdofsy_roundup;
          uint r = i/locdofsy_roundup;
          elvecy[r][c] = 0;
        }


      threadgroup_barrier(mem_flags::mem_threadgroup);

#if ($ONLY_LOADSTORE==1)

      for (uint i = threadIdx; i < $BS_ELS*locdofsy_roundup; i+= blockDim)
        {
          uint c = i/$BS_ELS;
          uint r = i%$BS_ELS;
          if (c < locdofsx_roundup)
            elvecy[r][c] = elvecx[r][c];
        }
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

              WarpMatrix<TS_IPTS,TS_ELS> pointvalsrefxi = 0;
              auto mb = mat_elvecx.SubMatrix(TS_ELS*eltile, 0);
              auto ma = mat_bmatx.SubMatrix(ip+nip8*comp, 0);

              pointvalsrefxi.AddMM<8*$DOFX_TILES>(ma, mb.Transpose(), threadIdx);
              pointvalsrefxi.Store(mat_pointvalsref.SubMatrix(TS_IPTS*iptile+bs_ipts*comp, TS_ELS*eltile), threadIdx);
            }

          threadgroup_barrier(mem_flags::mem_threadgroup);
  
          // work on integration points
          for (uint ip = threadIdx; ip < bs_ipts*bs_els ; ip += blockDim)
             {
               uint locelnr = ip % $BS_ELS;
               uint locipnr = ip / $BS_ELS;
               uint elnr = baseelem + locelnr;

#if ($ONLY_LOADSTOREB==0)
               Vec<$DIMXREF,float> xrefvals;
               Vec<$DIMYREF,float> yrefvals;
               Vec<$DIMX,float> xvals;
               Vec<$DIMY,float> yvals;

               Mat<$DIMR,$DIMS,float> F;
               for (int i = 0; i < $DIMR; i++)
                  for (int j = 0; j < $DIMS; j++)
                     F(i,j) = Jacobi[elnr + ne8 * (j + $DIMS*i)];
               float J = (elnr < ne) ?  JacobiDets[elnr] : 0;

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


          threadgroup_barrier(mem_flags::mem_threadgroup);

          // multiply with By.trans
          for (uint blocknr = warpIdx; blocknr < BSTS_ELS * $DOFY_TILES; blocknr += $WARPS)
            {
              int eltile = blocknr % BSTS_ELS;  // which els
              int ydoftile = blocknr / BSTS_ELS;  // which dofs

              WarpMatrix<TS_ELS,8,float> sum(mat_elvecy.SubMatrix(TS_ELS*eltile, 8*ydoftile), threadIdx);

              for (int comp = 0; comp < $DIMYREF; comp++)
                {
                   auto ma = mat_pointvalsref.SubMatrix(bs_ipts*comp, TS_ELS*eltile);
                   auto mb = mat_bmaty.SubMatrix(baseip+RoundUp<$BS_IPTS>(nip)*comp, 8*ydoftile);
                   sum.AddMM<$BS_IPTS> (ma.Transpose(), mb, threadIdx);
                }

              sum.Store(mat_elvecy.SubMatrix(8*eltile, 8*ydoftile), threadIdx);
           }

          threadgroup_barrier(mem_flags::mem_threadgroup);

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

             WarpMatrix<8,8> pointvalsrefxi = 0;

#ifdef OLD
             for (uint xdoftile = 0; xdoftile < $DOFX_TILES; xdoftile++)
                {
                  auto mb = mat_elvecx.SubMatrix(8*eltile, 8*xdoftile);
                  auto ma = mat_bmatx.SubMatrix($BMATX_REM_ROWS+8*iptile, 8*xdoftile);
                  pointvalsrefxi.AddMM<8>(ma, mb.Transpose(), threadIdx);
                }
#endif

              // slower for convection -> now good!
              auto mb = mat_elvecx.SubMatrix(8*eltile, 0);
              auto ma = mat_bmatx.SubMatrix($BMATX_REM_ROWS+8*iptile, 0);
              pointvalsrefxi.AddMM<RoundUp<8>(locdofsx)>(ma, mb.Transpose(), threadIdx);
              // pointvalsrefxi = AddMM<RoundUp<8>(locdofsx)> (pointvalsrefxi, ma, mb.Transpose(), threadIdx);
              pointvalsrefxi.Store(mat_pointvalsref.SubMatrix(8*iptile,8*eltile), threadIdx);
            }

          threadgroup_barrier(mem_flags::mem_threadgroup);
  
          // work on integration points
          for (uint ip = threadIdx; ip < numips*bs_els; ip += blockDim)
             {
               uint locelnr = ip % $BS_ELS;
               uint locipnr = ip / $BS_ELS;
               uint elnr = baseelem + locelnr;

#if ($ONLY_LOADSTOREB==0)
               Vec<$DIMXREF,float> xrefvals;
               Vec<$DIMYREF,float> yrefvals;
               Vec<$DIMX,float> xvals;
               Vec<$DIMY,float> yvals;

               Mat<$DIMR,$DIMS,float> F;
#pragma unroll
               for (uint i = 0; i < $DIMR; i++)
#pragma unroll
                  for (uint j = 0; j < $DIMS; j++)
                     F(i,j) = Jacobi[elnr + ne8 * (j + $DIMS*i)];
               float J = (elnr < ne) ?  JacobiDets[elnr] : 0;

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


          threadgroup_barrier(mem_flags::mem_threadgroup);

          // multiply with By.trans
          for (uint blocknr = warpIdx; blocknr < $EL_TILES * $DOFY_TILES; blocknr += $WARPS)
            {
              uint eltile = blocknr % $EL_TILES;  // which els
              uint ydoftile = blocknr / $EL_TILES;  // which dofs

              WarpMatrix<8,8,float> sum(mat_elvecy.SubMatrix(8*eltile, 8*ydoftile), threadIdx);

              for (uint ipcomp = 0; ipcomp < numips*$DIMYREF; ipcomp += 8)
                {
                   auto ma = mat_pointvalsref.SubMatrix(ipcomp, 8*eltile);
                   auto mb = mat_bmaty.SubMatrix($BMATY_REM_ROWS+ipcomp, 8*ydoftile);
                   sum.AddMM<8>(ma.Transpose(), mb, threadIdx);
                }

              sum.Store(mat_elvecy.SubMatrix(8*eltile, 8*ydoftile), threadIdx);


#ifdef NEW
              WarpMatrixV2<8,8,float> sum(mat_elvecy.SubMatrix(8*eltile, 8*ydoftile), threadIdx);
              auto ma = mat_pointvalsref.SubMatrix(0, 8*eltile);
              auto mb = mat_bmaty.SubMatrix($BMATY_REM_ROWS+0, 8*ydoftile);
              // sum.AddMM (RoundUp<4>(numips*$DIMYREF), ma.Transpose(), mb, threadIdx);
              sum.AddMM (numips*$DIMYREF, ma.Transpose(), mb, threadIdx);
              sum.Store(mat_elvecy.SubMatrix(8*eltile, 8*ydoftile), threadIdx);
#endif
           }

          threadgroup_barrier(mem_flags::mem_threadgroup);

      }  // if (numips_remainder > 0)

#endif

      threadgroup_barrier(mem_flags::mem_threadgroup);



      // store vector
      for (uint i = threadIdx; i < $BS_ELS*locdofsy; i += blockDim)
        {
           uint dofnr = i % locdofsy;
           uint locelnr = i / locdofsy;

           uint elnr = baseelem + locelnr;
           if (elnr < ne)
#if ($ATOMIC==1)
             atomic_fetch_add_explicit(&y[dofy[elnr*locdofsy+dofnr]], elvecy[locelnr][dofnr], memory_order_relaxed);
#else
             y[dofy[elnr*locdofsy+dofnr]] += elvecy[locelnr][dofnr];
#endif
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

    phys = Substitute(phys, "double", "float");
    code = Substitute(code, "$PHYSICS", phys);

    

    warps = pmat->opts.warps;

    
    code = Substitute(code, "$NE", ToString(ne)+" /*ne*/ ");
    code = Substitute(code, "$NIP", ToString(nip)+" /*nip*/ ");
    code = Substitute(code, "$BS_ELS", ToString(pmat->opts.BS_els)+" /*BS_ELS*/ ");
    code = Substitute(code, "$BS_IPTS", ToString(pmat->opts.BS_ipts)+" /*BS_IPTS*/ ");
    code = Substitute(code, "$DIMR", ToString(3)+" /*dimr*/ ");  
    code = Substitute(code, "$DIMS", ToString(3)+" /*dims*/ ");
    code = Substitute(code, "$WARPS", ToString(warps)+" /*warps*/ ");

    code = Substitute(code, "$BMATX_REM_ROWS", ToString(bmatx_rem_rows)+" /*bmatx_rem_rows*/ ");      
    code = Substitute(code, "$BMATY_REM_ROWS", ToString(bmaty_rem_rows)+" /*bmaty_rem_rows*/ ");      
    
    code = Substitute(code, "$IP_TILES", ToString(RoundDown<8>(nip)/8)+" /* IP_TILES */ ");
    code = Substitute(code, "$EL_TILES", ToString(RoundUp<8>(pmat->opts.BS_els)/8)+" /*EL_TILES*/ ");
    code = Substitute(code, "$DOFX_TILES", ToString(RoundUp<8>(locdofsx)/8) +" /* DOFX_TILES */ ");
    code = Substitute(code, "$DOFY_TILES", ToString(RoundUp<8>(locdofsy)/8) +" /* DOFY_TILES */ ");
    
    code = Substitute(code, "$LOCDOFSX", ToString(locdofsx));    
    code = Substitute(code, "$LOCDOFSY", ToString(locdofsy));
    code = Substitute(code, "$DIMXREF", ToString(dimxref));
    code = Substitute(code, "$DIMYREF", ToString(dimyref));
    code = Substitute(code, "$DIMX", ToString(dimx));
    code = Substitute(code, "$DIMY", ToString(dimy));
    code = Substitute(code, "$BS_IPTS", ToString(pmat->opts.BS_ipts));

    code = Substitute(code, "$ONLY_LOADSTOREB", ToString(pmat->opts.only_loadstoreB));    
    code = Substitute(code, "$ONLY_LOADSTORE", ToString(pmat->opts.only_loadstore));    
    code = Substitute(code, "$ATOMIC", ToString(pmat->opts.atomic));    

    string transxcode;
    for (uint i = 0 ; i < pmat->diffopsx.Size(); i++)
      {
        IntRange rangexref = pmat->ranges_xref[i];
        IntRange rangex = pmat->ranges_x[i];
        transxcode += "{\n";
        transxcode += "Vec<" + ToString(rangex.Size()) + ",float> res;\n";
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
        transycode += "Vec<" + ToString(rangeyref.Size()) + ",float> res;\n";
        transycode += pmat->diffopsy[i] -> GenerateTransformationCode("yvals.Range<"+ToString(rangey.First())+","+ToString(rangey.Next())+">()",
                                                                      "res", true);
        transycode += "yrefvals.SetRange<" + ToString(rangeyref.First()) +"," + ToString(rangeyref.Next()) + ">(res);\n";
        transycode += "}\n";
      }
    code = Substitute(code, "$TRANSFORMY", transycode);    

    
    ofstream codefile("applybtdtb.metal");
    codefile << code;
    codefile.close();
    
    NS::Error* error = nullptr;
    NS::String* mslCode = NS::String::string(code.c_str(), NS::UTF8StringEncoding);

    
    MTL::Library* library = GetDevice()->newLibrary(mslCode, nullptr, &error);

    if (!library)
      throw Exception("Shader compile error: " 
                      +string(error->localizedDescription()->utf8String()));

    NS::String* funcName = NS::String::string("apply_btdtb", NS::UTF8StringEncoding);
    ApplyBTDTB_Func = library->newFunction(funcName);

    pipelineState = GetDevice()->newComputePipelineState(ApplyBTDTB_Func, &error);
    NS::UInteger staticBytes = pipelineState->staticThreadgroupMemoryLength();
    // cout << "shared memory = " << staticBytes << endl;
    // cout << "maxtotalthreads = " << pipelineState->maxTotalThreadsPerThreadgroup() << endl;
    
    if (error)
      {
        std::cerr << "Metal Error: " << error->localizedDescription()->utf8String() << std::endl;
        throw Exception (error->localizedDescription()->utf8String());
      }
  }

  
  void MetalBTDTBMatrix ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer tfull("MetalBTDTBMNatrix::MultAll"); RegionTimer rfull(tfull);
        
    static Timer t("MultBTDTB-GPU");
      
    const MetalVector & mvx = dynamic_cast<const MetalVector&>(x);
    MetalVector & mvy = dynamic_cast<MetalVector&>(y);


    MTL::CommandBuffer* commandBuffer = GetCommandQueue()->commandBuffer();
    MTL::ComputeCommandEncoder* encoder = commandBuffer->computeCommandEncoder();
    
    encoder->setComputePipelineState(pipelineState);
    encoder->setBuffer(mvx.GetBuffer(), 0, 0);
    encoder->setBuffer(mvy.GetBuffer(), 0, 1);
    encoder->setBuffer(buffer_dofx, 0, 2);
    encoder->setBuffer(buffer_dofy, 0, 3);
    encoder->setBuffer(buffer_bmatx, 0, 4);
    encoder->setBuffer(buffer_bmaty, 0, 5);
    encoder->setBuffer(buffer_weights, 0, 6);
    encoder->setBuffer(buffer_Jacobi, 0, 7);
    encoder->setBuffer(buffer_JacobiDets, 0, 8);
    encoder->setBuffer(debug, 0, 9);      
    
    encoder->dispatchThreadgroups(MTL::Size(200,1,1), MTL::Size(warps*32, 1, 1));
    encoder->endEncoding();
    

    mvy.CommitAsync(commandBuffer);
  }
  
}
