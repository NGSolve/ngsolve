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


  int RoundUp8 (int i)
  {
    return 8 * ( (i+7) / 8 );
  }

  
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

    BS_els = pmat->opts.BS_els;
    ne = numels;
    
    FlatArray<int> hdofx = pmat->dofx.AsArray();
    buffer_dofx = GetDevice()->newBuffer(hdofx.Size()*sizeof(int), MTL::ResourceStorageModeShared);
    FlatArray<int> ddofx(hdofx.Size(), (int*)buffer_dofx->contents());
    ddofx = hdofx;

    FlatArray<int> hdofy = pmat->dofy.AsArray();
    buffer_dofy = GetDevice()->newBuffer(hdofy.Size()*sizeof(int), MTL::ResourceStorageModeShared);
    FlatArray<int> ddofy(hdofy.Size(), (int*)buffer_dofy->contents());
    ddofy = hdofy;

    buffer_bmatx = GetDevice()->newBuffer(dimxref*RoundUp8(nip)*RoundUp8(locdofsx)*sizeof(float), MTL::ResourceStorageModeShared);
    FlatTensor<3,float> dbmatx(dimxref,RoundUp8(nip),RoundUp8(locdofsx), (float*)buffer_bmatx->contents());
    for (int j = 0; j < dimxref; j++)
      for (int i = 0; i < RoundUp8(nip); i++)
        for (int k = 0; k < RoundUp8(locdofsx); k++)
          dbmatx(j,i,k) = (i < nip && k < locdofsx) ? pmat->Bx(k,j,i) : 0;

    
    buffer_bmaty = GetDevice()->newBuffer(dimyref*RoundUp8(nip)*dimyref*RoundUp8(locdofsy)*sizeof(float), MTL::ResourceStorageModeShared);
    FlatTensor<3,float> dbmaty(dimyref,RoundUp8(nip),RoundUp8(locdofsy), (float*)buffer_bmaty->contents());
    for (int j = 0; j < dimyref; j++)
      for (int i = 0; i < RoundUp8(nip); i++)
        for (int k = 0; k < RoundUp8(locdofsy); k++)
          dbmaty(j,i,k) = (i<nip && k < locdofsy) ? pmat->By(k,j,i) : 0;
    
    // cout << "dbmatx = " << endl << dbmatx << endl;
    // cout << "dbmaty = " << endl << dbmaty << endl;


    buffer_weights = GetDevice()->newBuffer(RoundUp8(nip)*sizeof(float), MTL::ResourceStorageModeShared);
    for (int i = 0; i < RoundUp8(nip); i++)
      ((float*)buffer_weights->contents())[i] = i < nip ? pmat->weights[i] : 0;
    
    buffer_JacobiDets = GetDevice()->newBuffer(ne*sizeof(float), MTL::ResourceStorageModeShared);    
    for (int i = 0; i < ne; i++)
      ((float*)buffer_JacobiDets->contents())[i] = Det(pmat->Jacobi(i,STAR,STAR,0));


    int ne8 = RoundUp8(ne);
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
      uint   gridDim         [[threadgroups_per_grid]]
      )
      {
      constexpr int ne = $NE;
      constexpr int ne8 = 8 * ((ne+7)/8);
      constexpr int nip = $NIP;
      constexpr int nip8 = $RUNIP;   // roundup
      constexpr int bs_els = $BS_ELS;
      constexpr int bs_ipts = $BS_IPTS;
      constexpr int locdofsx = $LOCDOFSX;
      constexpr int locdofsy = $LOCDOFSY;
      constexpr int locdofsx_lanes = (locdofsx+7)/8;
      constexpr int locdofsy_lanes = (locdofsy+7)/8;
      constexpr int locdofsx_roundup = 8*locdofsx_lanes;
      constexpr int locdofsy_roundup = 8*locdofsy_lanes;

      threadgroup float elvecx[$BS_ELS][locdofsx_roundup];
      threadgroup float elvecy[$BS_ELS][locdofsy_roundup];

      constexpr int MAXDIMREF = ($DIMXREF>$DIMYREF) ? $DIMXREF : $DIMYREF;
      threadgroup float pointvalsref[MAXDIMREF*$BS_IPTS][$BS_ELS];


      // zero elvecx (overhead would be enough
      for (int i = threadIdx; i < $BS_ELS*locdofsx_roundup; i+= blockDim)
        {
          int c = i%locdofsx_roundup;
          int r = i/locdofsx_roundup;
          elvecx[r][c] = 0;
        }

      threadgroup_barrier(mem_flags::mem_threadgroup);

      for (int baseelem = blockIdx*$BS_ELS; baseelem < $NE; baseelem += gridDim*$BS_ELS) { 

      threadgroup_barrier(mem_flags::mem_threadgroup);

      // load element vectors
      for (int i = threadIdx; i < $BS_ELS*locdofsx; i += blockDim)
        {
           int dofnr = i % locdofsx;
           int locelnr = i / locdofsx;
           int elnr = baseelem + locelnr;
           elvecx[locelnr][dofnr] = (elnr < ne) ? x[dofx[elnr*locdofsx+dofnr]] : 0;
        }

      // zero elvecy
      for (int i = threadIdx; i < $BS_ELS*locdofsy_roundup; i+= blockDim)
        {
          int c = i/$BS_ELS;
          int r = i%$BS_ELS;
          elvecy[r][c] = 0;
        }


      threadgroup_barrier(mem_flags::mem_threadgroup);

#if ($ONLY_LOADSTORE==1)

      for (int i = threadIdx; i < $BS_ELS*locdofsy_roundup; i+= blockDim)
        {
          int c = i/$BS_ELS;
          int r = i%$BS_ELS;
          if (c < locdofsx_roundup)
            elvecy[r][c] = elvecx[r][c];
        }
#else

      for (int baseiptile = 0; baseiptile < $IP_TILES; baseiptile += bs_ipts/8)
        {
          int myiptiles = ($IP_TILES-baseiptile < bs_ipts/8) ? $IP_TILES-baseiptile : bs_ipts/8;

          // multiply with Bx
          for (int blocknr = threadIdx/32; blocknr < myiptiles * $EL_TILES; blocknr += blockDim/32)
            {
              int eltile = blocknr % $EL_TILES;  // which els
              int iptile = blocknr / $EL_TILES;  // which pts

              int ip = 8*(baseiptile+iptile);

              simdgroup_float8x8 ma, mb;
              simdgroup_float8x8 pointvalsrefxi[$DIMXREF];
              for (int k = 0; k < $DIMXREF; k++)
                pointvalsrefxi[k] = make_filled_simdgroup_matrix<float, 8, 8>(0.0f);

              for (int xdoftile = 0; xdoftile < $DOFX_TILES; xdoftile++)
                {
                  simdgroup_load(mb, &elvecx[0][0], locdofsx_roundup, ulong2(8*xdoftile, 8*eltile), true);
                  for (int l = 0; l < $DIMXREF; l++)
                     {
                        simdgroup_load(ma, bmatx, locdofsx_roundup, ulong2(8*xdoftile, ip + nip8*l)); 
                        simdgroup_multiply_accumulate(pointvalsrefxi[l], ma, mb, pointvalsrefxi[l]);
                     }
                }

              for (int l = 0; l < $DIMXREF; l++)
                simdgroup_store(pointvalsrefxi[l], &pointvalsref[0][0], bs_els, ulong2(8*eltile, 8*iptile+bs_ipts*l));
            }

          threadgroup_barrier(mem_flags::mem_threadgroup);
  
          // work on integration points
          for (int ip = threadIdx; ip < 8*myiptiles*bs_els ; ip += blockDim)
             {
               int locelnr = ip % $BS_ELS;
               int locipnr = ip / $BS_ELS;
               int elnr = baseelem + locelnr;

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

               for (int j = 0; j < $DIMXREF; j++)
                  xrefvals(j) = pointvalsref[locipnr+j*bs_ipts][locelnr];

              // xrefvals -> xvals
              $TRANSFORMX;   

              $PHYSICS
              yvals = weights[8*baseiptile+locipnr] * JacobiDets[elnr] * yvals;

              // yvals -> yrefvals
              $TRANSFORMY;          

               for (int j = 0; j < $DIMYREF; j++)
                  pointvalsref[locipnr+j*bs_ipts][locelnr] = yrefvals(j);
#endif
            }


          threadgroup_barrier(mem_flags::mem_threadgroup);

          // multiply with By.trans
          for (int blocknr = threadIdx/32; blocknr < $EL_TILES * $DOFY_TILES; blocknr += blockDim/32)
            {
              int eltile = blocknr % $EL_TILES;  // which els
              int ydoftile = blocknr / $EL_TILES;  // which dofs

              simdgroup_float8x8 ma, mb;
              simdgroup_float8x8 sum;
              simdgroup_load(sum, &elvecy[0][0], locdofsy_roundup, ulong2(8*ydoftile, 8*eltile));

              for (int iptile = 0; iptile < myiptiles; iptile++)
                for (int l = 0; l < $DIMYREF; l++)
                  {
                     simdgroup_load(ma, &pointvalsref[0][0], bs_els, ulong2(8*eltile, 8*iptile+bs_ipts*l), true);
                     simdgroup_load(mb, bmaty, locdofsy_roundup, ulong2(8*ydoftile, 8*(baseiptile+iptile+$IP_TILES*l)), false);
                     simdgroup_multiply_accumulate(sum, ma, mb, sum);
                  }

              simdgroup_store(sum, &elvecy[0][0], locdofsy_roundup, ulong2(8*ydoftile, 8*eltile));           
           }

          threadgroup_barrier(mem_flags::mem_threadgroup);

         } // end base_ip
#endif


      threadgroup_barrier(mem_flags::mem_threadgroup);

      // store vector
      for (int i = threadIdx; i < $BS_ELS*locdofsy; i += blockDim)
        {
           int dofnr = i % locdofsy;
           int locelnr = i / locdofsy;

           int elnr = baseelem + locelnr;
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

        
        // auto& code = pmat->physics[k];
        phys += "{  // equation " + ToString(k) +"\n";
        
        // TODO: get all trial proxies and the right indices

        for (auto step : Range(compiledcf->Steps()))
          {
            auto stepcf = compiledcf->Steps()[step];
            
            if (auto proxycf = dynamic_cast<ProxyFunction*> (stepcf))
              {
                auto pos = pmat->trial_proxies.Pos(proxycf);
                if (pos != pmat->trial_proxies.ILLEGAL_POSITION)
                  {
                    phys += "auto values_"+ToString(step)+" = [&](int ip, int nr) { return xvals("+ToString(pmat->ranges_x[pos].First())+"+nr); };\n";
                    // s << "auto values_" << step << " = [dist_input,input](size_t i, int comp)\n"
                    // " { return input[i + (comp+" << proxyoffset[pos] << ")*dist_input]; };\n";
                    phys += "bool constexpr has_values_" + ToString(step) + " = true;\n";
                    // s << "bool constexpr has_values_" << step << " = true;\n" << endl;
                    
                    // Declare dummy com_ variables to avoid compile errors (won't be used since has_values = true)
                    for(auto i : Range(proxycf->Dimension()))
                      phys += Var("comp", step,i, proxycf->Dimensions()).Declare("double", 0.0);
                      // s << Var("comp", step,i, proxycf->Dimensions()).Declare("double", 0.0);
                  }
              }
          }
    
        
        // phys += "auto values_0 = [&](int ip, int nr) { return xvals(nr); };\n";
        // phys += "constexpr bool has_values_0 = true; \n";
        // phys += "float comp_0_0, comp_0_1, comp_0_2;\n";
    
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

    



    
    
    code = Substitute(code, "$NE", ToString(ne)+" /*ne*/ ");
    code = Substitute(code, "$NIP", ToString(nip)+" /*nip*/ ");
    code = Substitute(code, "$RUNIP", ToString(RoundUp8(nip)));
    code = Substitute(code, "$BS_ELS", ToString(pmat->opts.BS_els)+" /*BS_ELS*/ ");
    code = Substitute(code, "$BS_IPTS", ToString(pmat->opts.BS_ipts)+" /*BS_IPTS*/ ");
    code = Substitute(code, "$DIMR", ToString(3)+" /*dimr*/ ");  
    code = Substitute(code, "$DIMS", ToString(3)+" /*sims*/ ");  
    
    code = Substitute(code, "$IP_TILES", ToString(RoundUp8(nip)/8)+" /* IP_TILES */ ");
    code = Substitute(code, "$EL_TILES", ToString(RoundUp8(pmat->opts.BS_els)/8)+" /*EL_TILES*/ ");
    code = Substitute(code, "$DOFX_TILES", ToString(RoundUp8(locdofsx)/8) +" /* DOFX_TILES */ ");
    code = Substitute(code, "$DOFY_TILES", ToString(RoundUp8(locdofsy)/8) +" /* DOFY_TILES */ ");
    
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
    for (int i = 0 ; i < pmat->diffopsx.Size(); i++)
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
    for (int i = 0 ; i < pmat->diffopsy.Size(); i++)
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

    
    ofstream codefile("metalcode.cpp");
    codefile << code;
    codefile.close();
    
    NS::Error* error = nullptr;
    NS::String* mslCode = NS::String::string(code.c_str(), NS::UTF8StringEncoding);
    MTL::Library* library = GetDevice()->newLibrary(mslCode, nullptr, &error);

    if (!library) {
      throw Exception("Shader compile error: " 
                      +string(error->localizedDescription()->utf8String()));
    }

    {
      NS::String* funcName = NS::String::string("apply_btdtb", NS::UTF8StringEncoding);
      ApplyBTDTB_Func = library->newFunction(funcName);
    }

    // NS::Error* error = nullptr;      
    pipelineState = GetDevice()->newComputePipelineState(ApplyBTDTB_Func, &error);
    if (error)
      std::cerr << "Metal Error: " << error->localizedDescription()->utf8String() << std::endl;
    
    
  }

  void MetalBTDTBMatrix ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer tfull("MetalBTDTBMNatrix::MultAll"); RegionTimer rfull(tfull);
        
    static Timer t("MultBTDTB-GPU");
      
    const MetalVector & mvx = dynamic_cast<const MetalVector&>(x);
    MetalVector & mvy = dynamic_cast<MetalVector&>(y);


    // for (int runs = 0; runs < 10; runs++)
      {
      
      
        // 6. Encode and submit compute commands
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

        // encoder->dispatchThreads(MTL::Size(ne/BS_els+1,1,1), MTL::Size(32, BS_els, 1));
        // encoder->dispatchThreadgroups(MTL::Size(ne/BS_els+1,1,1), MTL::Size(16*32, 1, 1));
        //       for (int runs = 0; runs < 10; runs++)
        encoder->dispatchThreadgroups(MTL::Size(200,1,1), MTL::Size(16*32, 1, 1));
        encoder->endEncoding();


      
        // t.Start();
        
        // commandBuffer->commit();
        // commandBuffer->waitUntilCompleted();
        mvy.CommitAsync(commandBuffer);
        // t.Stop();

        /*
        CFTimeInterval gpuStartTime = commandBuffer->GPUStartTime();
        CFTimeInterval gpuEndTime = commandBuffer->GPUEndTime();
        double gpuDurationMS = (gpuEndTime - gpuStartTime) * 1000.0;
        cout << "gpu times measure: " << gpuDurationMS << endl;
        */
        
        // cout << "Applyt BTDTB time = " << t.GetTime() << endl;
      }
    // cout << "debug = " << debugvec.Range(0, 128);
  }

}
