#include <fespace.hpp>
#include <elementbyelement.hpp>


#include "bilinearform.hpp"
#include <diagonalmatrix.hpp>

#include "../fem/h1lofe.hpp"


using namespace ngcomp;

namespace ngcomp
{
  // *******************************************************************************
  // ****************** free-standing function for an operator *********************
  // *******************************************************************************
  
  shared_ptr<BaseMatrix> CreateMatrixFreeOperator (shared_ptr<FESpace> source_fes,
                                                   shared_ptr<FESpace> target_fes,
                                                   std::function<tuple<Matrix<>,Array<int>,Array<int>>(ElementId)> creator,
                                                   LocalHeap & lh)
  {
    // bool mixed = target_fes != nullptr;
    if (!target_fes) target_fes = source_fes;
    auto ma = source_fes->GetMeshAccess();


    const Table<size_t> & table = ma->GetElementsOfClass();

    shared_ptr<BaseMatrix> sum;
  
    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;
      
        ElementId ei(VOL,elclass_inds[0]);
        // auto & trafo = ma->GetTrafo(ei, lh);
        // auto & felx = source_fes->GetFE (ei, lh);
        // auto & fely = GetTestSpace()->GetFE (ei, lh);
        // MixedFiniteElement fel(felx, fely);

        auto [elmat,sind,tind] = creator(ei);
        
        Table<DofId> sdofs(elclass_inds.Size(), sind.Size());
        Table<DofId> tdofs(elclass_inds.Size(), tind.Size());

        Array<DofId> dofs;
        for (auto i : Range(elclass_inds.Size()))
          {
            source_fes->GetDofNrs(ElementId(VOL,elclass_inds[i]), dofs);
            for (int j = 0; j < sind.Size(); j++)
              sdofs[i][j] = dofs[sind[j]];
            for (int j = 0; j < tind.Size(); j++)
              tdofs[i][j] = dofs[tind[j]];
          }

        auto mat = make_shared<ConstantElementByElementMatrix<>>
          (target_fes->GetNDof(), source_fes->GetNDof(),
           elmat, std::move(tdofs), std::move(sdofs));
      
        if (sum)
          sum = make_shared<SumMatrix>(sum, mat);
        else
          sum = mat;
      }

    return sum;
  }






  // ********************************************************************************************
  // **************** Composition of several operators, code gen for integration points 
  // ********************************************************************************************
  


  /*
    Stores product of B2 D B1 matrices

    B matrices ... sum of EBE const 
    D bock-diagonal
   */

  ApplyIntegrationPoints :: 
  ApplyIntegrationPoints (Array<shared_ptr<CoefficientFunction>> acoefs,
                          const Array<ProxyFunction*> & atrialproxies,
                          Matrix<> apoints, Matrix<> anormals,
                          size_t adimx, size_t adimy, size_t anip)
    : coefs(acoefs), trialproxies{atrialproxies}, dimx(adimx), dimy(adimy), nip(anip),
      points(std::move(apoints)), normals(std::move(anormals))
  { 
      // make my own code
    
    Array<int> proxyoffset;
    int starti = 0;
    for (auto proxy : trialproxies)
      {
        proxyoffset.Append (starti);
        starti += proxy->Evaluator()->Dim();
      }
    
    stringstream s;
    s <<
      "#include <cstddef>\n"
      "#include <cmath>\n"
      "using std::sqrt;\n"
      "static inline double L2Norm2 (double v) { return v*v; }\n"
      "extern \"C\" void ApplyIPFunction (size_t nip, double * input, size_t dist_input,\n"
      "                      double * output, size_t dist_output,\n"
      "                      size_t dist, double * pnts, double * nvs,\n"
      "                      double * coef_input, int * elidx) {\n";

    int base_output = 0;
    for (auto cf : coefs)
      {
        auto compiledcf = Compile (cf, false);
        Code code = compiledcf->GenerateProgram(0, false);
        
        s << "{\n";
        // cout << code.header << endl;
        
        for (auto step : Range(compiledcf->Steps()))
          {
            auto stepcf = compiledcf->Steps()[step];

            bool uses_values =
              code.body.find("values_" + ToString(step) + "(") != string::npos;

            if (auto proxycf = dynamic_cast<ProxyFunction*> (stepcf))
              {
                auto pos = trialproxies.Pos(proxycf);
                if (pos != trialproxies.ILLEGAL_POSITION)
                  {
                    s << "auto values_" << step << " = [dist_input,input](size_t i, int comp)\n"
                      " { return input[i + (comp+" << proxyoffset[pos] << ")*dist_input]; };\n";
                    s << "bool constexpr has_values_" << step << " = true;\n" << endl;
                    // Declare dummy com_ variables to avoid compile errors (won't be used since has_values = true)
                    for(auto i : Range(proxycf->Dimension()))
                      s << Var("comp", step,i, proxycf->Dimensions()).Declare("double", 0.0);
                  }
              }
            else if (uses_values)
              {
                // a coefficient function (e.g. a GridFunction) that has to be
                // evaluated at the integration points.
                auto pos = input_coefs.Pos(stepcf);
                int off;
                if (pos == input_coefs.ILLEGAL_POSITION)
                  {
                    off = dim_coef;
                    input_coefs.Append (stepcf);
                    input_coef_offset.Append (off);
                    dim_coef += stepcf->Dimension();
                  }
                else
                  off = input_coef_offset[pos];

                s << "auto values_" << step << " = [dist,coef_input](size_t i, int comp)\n"
                  " { return coef_input[i + (comp+" << off << ")*dist]; };\n" << endl;
              }
          }

        s << "[[maybe_unused]] auto points = [dist,pnts](size_t i, int comp)\n"
          " { return pnts[i+comp*dist]; };\n";
        s << "[[maybe_unused]] auto normals = [dist,nvs](size_t i, int comp)\n"
          " { return nvs[i+comp*dist]; };\n";
        
        s << "for (size_t i = 0; i < nip; i++) {\n";
        if (code.body.find("domain_index") != string::npos)
          {
            needs_element_index = true;
            s << "[[maybe_unused]] int domain_index = elidx[i];\n";
          }
        s << code.body << endl;
        
          // missing: last step nr
        for (int j = 0; j < cf->Dimension(); j++)
          s << "output[i+"<<base_output+j<<"*dist_output] = "
            << Var(compiledcf->Steps().Size()-1, j, cf->Dimensions()).code << ";\n";
        base_output += cf->Dimension();
        
        s << "}\n}";
      }
    s << "}\n";
    
    
    // cout << s.str() << endl;

    try
      {
        library = CompileCode ( { s.str() }, {} );
        compiled_function = library->GetSymbol<lib_function> ("ApplyIPFunction");
      }
    catch (const Exception & e)
      { ; } 
  }

  AutoVector ApplyIntegrationPoints :: CreateColVector() const
  {
    return make_unique<VVector<double>> (nip*dimy);
  }

  AutoVector ApplyIntegrationPoints :: CreateRowVector() const
  {
    return make_unique<VVector<double>> (nip*dimx);
  }

  void ApplyIntegrationPoints :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("ApplyIntegrationPoints"); RegionTimer reg(t);
    static Timer teval("ApplyIntegrationPoints eval");
    static Timer tmir("ApplyIntegrationPoints mir");      
    static Timer ttransx("ApplyIntegrationPoints transx");
    static Timer ttransy("ApplyIntegrationPoints transy");

    if (dim_coef > 0 && coef_values.Height() != size_t(dim_coef))
      throw Exception ("ApplyIntegrationPoints::Mult - coefficient values not set "
                       "(SetCoefValues has to be called for integrands with non-proxy "
                       "coefficient functions)");
    if (needs_element_index && element_index.Size() != nip)
      throw Exception ("ApplyIntegrationPoints::Mult - element indices not set "
                       "(SetElementIndex has to be called for integrands with "
                       "domain-wise coefficient functions)");

    if (compiled_function)
      {
        FlatMatrix<double> mx = x.FV<double>().AsMatrix(dimx, nip);
        FlatMatrix<double> my = y.FV<double>().AsMatrix(dimy, nip);
        ParallelForRange(nip, [this,mx, my, pts=FlatMatrix<>(points), nvs=FlatMatrix<>(normals),
                               cvals=FlatMatrix<>(coef_values)] (IntRange r)
                         {
                           this->compiled_function(r.Size(),
                                                   mx.Cols(r).Data(), mx.Dist(),
                                                   my.Cols(r).Data(), my.Dist(),
                                                   nip, pts.Cols(r).Data(), nvs.Cols(r).Data(),
                                                   cvals.Height() ? cvals.Cols(r).Data() : nullptr,
                                                   element_index.Size() ? element_index.Data()+r.First() : nullptr);
                         });
        return;
      }

    if (needs_element_index)
      throw Exception ("ApplyIntegrationPoints::Mult - integrand contains a "
                       "domain-wise coefficient function, which is only supported "
                       "by the compiled (JIT) path; the JIT compilation apparently "
                       "failed.");

    try
      {
    ParallelForRange
      (nip, [&](IntRange r)
       {
         constexpr size_t BS = 256;
         LocalHeap lh(1000000);
         for (size_t ii = r.First(); ii < r.Next(); ii+=BS)
           {
             IntRange r2(ii, min(ii+BS, r.Next()));               
             HeapReset hr(lh);
             
             SIMD_IntegrationRule simdir(r2.Size(), lh);

             FE_ElementTransformation<2,2> trafo2d(ET_TRIG);
             FE_ElementTransformation<3,3> trafo3d(ET_TET);

             SIMD_MappedIntegrationRule<2,2> simdmir2d(simdir, trafo2d, 0, lh); // don't actually compute
             SIMD_MappedIntegrationRule<3,3> simdmir3d(simdir, trafo3d, 0, lh); // don't actually compute

             const SIMD_BaseMappedIntegrationRule & simdmir =
               (points.Height() == 2) ?
               ((SIMD_BaseMappedIntegrationRule&)simdmir2d) : 
               ((SIMD_BaseMappedIntegrationRule&)simdmir3d);
             
             // tmir.Stop();
             
             ProxyUserData ud(trialproxies.Size(), input_coefs.Size(), lh);
             ScalarFE<ET_TRIG,1> dummyfe2d;
             ScalarFE<ET_TET,1> dummyfe3d;
             if (points.Height()==2)
               {
                 trafo2d.userdata = &ud;
                 ud.fel = &dummyfe2d;
               }
             else
               {
                 trafo3d.userdata = &ud;                 
                 ud.fel = &dummyfe3d;
               }
             
             int starti = 0;
             for (auto proxy : trialproxies)
               {
                 int nexti = starti + proxy->Evaluator()->Dim();
                 ud.AssignMemory (proxy, r2.Size(), proxy->Evaluator()->Dim(), lh);
                 
                 // ttransx.Start();
                 SliceMatrix<double> (nexti-starti, r2.Size(), SIMD<double>::Size()*simdmir.Size(),
                                      (double*)(ud.GetAMemory (proxy)).Data()) = 
                   x.FV<double>().AsMatrix(dimx, nip).Cols(r2).Rows(starti, nexti);
                 // ttransx.Stop();
                 
                 starti = nexti;
               }

             for (auto ci : Range(input_coefs))
               {
                 auto icf = input_coefs[ci];
                 int off = input_coef_offset[ci];
                 int d = icf->Dimension();
                 ud.AssignMemory (icf, r2.Size(), d, lh, false);
                 SliceMatrix<double> (d, r2.Size(), SIMD<double>::Size()*simdmir.Size(),
                                      (double*)(ud.GetAMemory (icf)).Data()) =
                   coef_values.Rows(off, off+d).Cols(r2);
                 ud.SetComputed (icf, true);
               }


             // teval.Start();
             FlatMatrix<SIMD<double>> simdres(dimy, simdmir.Size(), lh);
             starti = 0;
             for (auto coef : coefs)
               {
                 int nexti = starti + coef->Dimension();
                 coef -> Evaluate (simdmir, simdres.Rows(starti, nexti));
                 starti = nexti;
               }
             // teval.Stop();
             
             // ttransy.Start();
             SliceMatrix<> res = y.FV<double>().AsMatrix(dimy, nip).Cols(r2);
             res = SliceMatrix<double> (dimy, r2.Size(), SIMD<double>::Size()*simdmir.Size(),
                                        (double*)simdres.Data());
             // ttransy.Stop();
           }
       }, TasksPerThread(3));
      }
    catch (Exception & e)
      {
        cout << "In ApplyIntegrationRule, e = " << e.What() << endl;
        throw(e);
      }
  }


  MatrixFreeBTDTB ::
  MatrixFreeBTDTB (size_t h, size_t w,
                   Array<size_t> _elnums,
                   Table<DofId> _dofx, Table<DofId> _dofy,
                   Tensor<3> _Bx,  // locdofs, dim, nip
                   Tensor<3> _By,  // locdofs, dim, nip
                   // IntegrationRule _ir,
                   Vector<> _weights,  // ref-element intweights
                   Array<shared_ptr<DifferentialOperator>> _diffopsx,
                   Array<shared_ptr<DifferentialOperator>> _diffopsy,
                   Tensor<4> _D, // element, nip, dimy, dimx;
                   Tensor<4> _Jacobi,
                   MatFreeOptions _opts)
  : height(h), width(w), elnums(std::move(_elnums)), dofx(std::move(_dofx)), dofy(std::move(_dofy)),
    Bx(std::move(_Bx)), By(std::move(_By)),
    weights(std::move(_weights)),
    diffopsx(std::move(_diffopsx)), diffopsy(std::move(_diffopsy)),
    D(std::move(_D)), Jacobi(std::move(_Jacobi)), opts(_opts)
  {

    size_t starti = 0, startiref = 0;
    for (auto &dopx : diffopsx)
      {
        size_t nexti = starti + dopx->Dim();
        size_t nextiref = startiref + dopx->DimRef();
        ranges_x += IntRange(starti, nexti);
        ranges_xref += IntRange(startiref, nextiref);
        starti = nexti;
        startiref = nextiref;
      }

    starti = 0, startiref = 0;
    for (auto &dopy : diffopsy)
      {
        size_t nexti = starti + dopy->Dim();
        size_t nextiref = startiref + dopy->DimRef();
        ranges_y += IntRange(starti, nexti);
        ranges_yref += IntRange(startiref, nextiref);
        starti = nexti;
        startiref = nextiref;
      }

    
    if (opts.generate_code)
      {
        auto [locdofsx, dimxref, nip] = Bx.Shape();
        auto [locdofsy, dimyref, nip_] = By.Shape();
        auto [numels,dimy,dimx,nipD] = D.Shape();
        auto [numels2,dimr,dims,nipJ] = Jacobi.Shape();

        
        stringstream s;
        s <<
          "#include <core/simd.hpp>\n"
          "#include <matrix.hpp>\n"
          "#include <tensor.hpp>\n"
          "using namespace ngbla;\n";


        s << "[[maybe_unused]] static double bmatx[] = { ";
        for (int i = 0; i < nip; i++)
          for (int j = 0; j < dimxref; j++)
            for (int k = 0; k < locdofsx; k++)
              s << setprecision(17) << Bx(k,j,i) << ", ";
        s << "};\n";


        s << "[[maybe_unused]] static double bmaty[] = { ";
        for (int i = 0; i < nip; i++)
          for (int k = 0; k < locdofsy; k++)
            for (int j = 0; j < dimyref; j++)
              s << setprecision(17) << By(k,j,i) << ", ";
        s << "};\n";

        
        
        s << "extern \"C\" void AddBTDTB (double s, FlatVector<> fx, FlatVector<> fy, FlatTable<int> dofxtable, FlatTable<int> dofytable,"
          "FlatTensor<4> Jacobi, FlatVector<double> weights,  size_t numels) {\n";



        if (opts.timers)
          {
            s << 
              "static Timer t(\"MatFreeBTDTB\"); RegionTimer rt(t);\n"
              "static Timer tel(\"MatFreeBTDTB - element\");\n"
              "static Timer tL(\"MatFreeBTDTB-Load\");\n"              
              "static Timer tB(\"MatFreeBTDTB-B"<< opts.BS_els*opts.BS_ipts*dimx*locdofsx <<"\");\n"
              "static Timer tBt(\"MatFreeBTDTB-Bt\");\n"
              "static Timer trem(\"MatFreeBTDTB-remainder\");\n"
              "static Timer tS(\"MatFreeBTDTB-Store\");\n";
          }

        s <<
          "static constexpr int BS_ELS = " << opts.BS_els  << ";\n"          
          "static constexpr int BS_IPTS = " << opts.BS_ipts << ";\n"          
          "ParallelForRange (numels/BS_ELS, [&](IntRange r) {\n"

          "Vec<"+ToString(locdofsx)+",SIMD<double,BS_ELS>> elvecx;\n"
          "Vec<"+ToString(locdofsy)+",SIMD<double,BS_ELS>> elvecy;\n"
          "Mat<"+ToString(dimxref)+",BS_IPTS,SIMD<double,BS_ELS>> pointvalsrefx;\n"
          "Mat<"+ToString(dimyref)+",BS_IPTS,SIMD<double,BS_ELS>> pointvalsrefy;\n"
          "Mat<"+ToString(dimx)+",BS_IPTS,SIMD<double,BS_ELS>> pointvalsx;\n"
          "Mat<"+ToString(dimy)+",BS_IPTS,SIMD<double,BS_ELS>> pointvalsy;\n";

        s << "for (auto i : r) {\n";
        if (opts.timers) s << "RegionTimer rtel(tel);\n";
        if (opts.timers) s << "tL.Start();\n";
        s << 
          "std::array<int,"<<locdofsx<<"> * dofx = (std::array<int,"<<locdofsx<<">*)&dofxtable[i*BS_ELS][0];\n"
          "std::array<int,"<<locdofsy<<"> * dofy = (std::array<int,"<<locdofsy<<">*)&dofytable[i*BS_ELS][0];\n"    
          "for (int k = 0; k < " << locdofsx << "; k++)\n"
          "elvecx(k) = SIMD<double,BS_ELS>([&](int nr) { return fx(dofx[nr][k]); });\n"; 
        if (opts.timers) s << "tL.Stop();\n";

        if (!opts.only_loadstore)
          {
            s << "elvecy = SIMD<double,BS_ELS>(0.0); \n";
            s << "size_t ix = 0;\n";
            s << "size_t iy = 0;\n";
            
            if (nipJ == 1)  // Geometry in one point
              {
                size_t dist0 = Jacobi.GetDist();
                size_t dist1 = Jacobi.GetSubTensor().GetDist();
                size_t dist2 = Jacobi.GetSubTensor().GetSubTensor().GetDist();
                s << "Mat<" << dimr << "," << dims << ", SIMD<double,BS_ELS>> F; \n";
                for (size_t i = 0; i < dimr; i++)
                  for (size_t j = 0; j < dims; j++)
                    // s << "F(" << i << "," << j << ") = SIMD<double,BS_ELS>([&](auto nr) { return Jacobi(i*BS_ELS+nr," << i << "," << j << ",0);});\n";
                    s << "F(" << i << "," << j << ") = SIMD<double,BS_ELS>([&](auto nr) { return Jacobi.Data()[(i*BS_ELS+nr)*"<<dist0<<"+" << i*dist1+j*dist2<< "];});\n";                    
                s << "[[maybe_unused]] SIMD<double,BS_ELS> J = Det(F);\n";
              }
            
            
            s << "for (size_t j1 = 0; j1+BS_IPTS <= " << nip << "; j1+=BS_IPTS) { \n";

            /*
            s << "for (size_t j2 = 0; j2 < BS_IPTS; j2++) {\n";

            s <<
              "for (size_t i = 0; i < " << dimxref << "; i++) { \n"
              "SIMD<double,BS_ELS> sum(0.0);\n"
              "for (size_t k = 0; k < " << locdofsx << "; k++)\n"
              "sum += bmatx[ix++]*elvecx(k); \n"
              "pointvalsrefx(i,j2) = sum;\n"
              "}\n";
            */

            if (opts.timers) s << "tB.Start();\n";
            s << "for (size_t i = 0; i < " << dimxref << "; i++) { \n";
            for (size_t j2 = 0; j2 < opts.BS_ipts; j2++)
              s << "SIMD<double,BS_ELS> sum" << j2 << "(0.0);\n";
            s << "for (size_t k = 0; k < " << locdofsx << "; k++) { \n";
            for (size_t j2 = 0; j2 < opts.BS_ipts; j2++)
              s << "sum" << j2 << " += bmatx[ix+" << dimxref*locdofsx*j2 << "]*elvecx(k); \n";
            // s << "sum" << j2 << " += bmatx[ix+" << j2 << "]*elvecx(k); \n";
            s << "ix++;\n } \n";
            for (size_t j2 = 0; j2 < opts.BS_ipts; j2++)
              s << "pointvalsrefx(i," << j2 << ") = sum" << j2 << ";\n";
            s << "}  // end for i\n"; 
            s << "ix += " << dimxref*(opts.BS_ipts-1)*locdofsx << ";\n";
            if (opts.timers) s << "tB.Stop();\n";
            
            

            if (!opts.only_loadstoreB)
              {
                s << "Vec<"<<dimx<<",SIMD<double,BS_ELS>> hvx;\n";
                s << "Vec<"<<dimy<<",SIMD<double,BS_ELS>> hvy;\n";
                s << "for (size_t j2 = 0; j2 < BS_IPTS; j2++) { \n";
                  {
                    for (size_t i = 0; i < diffopsx.Size(); i++)
                      {
                        IntRange rref = ranges_xref[i];
                        IntRange r = ranges_x[i];
                        s << diffopsx[i]->GenerateTransformationCode ("pointvalsrefx.Col(j2).Range("+ToString(r.First())+","+ToString(r.Next())+")",
                                                                      "hvx.Range("+ToString(rref.First())+","+ToString(rref.Next())+")",
                                                                      false);
                      }
                    // s << "hv *= weights(j1+"<<j2<<") * J;\n";
                    // s << "pointvalsy.Col(" << j2 << ") = weights(j1+" << j2 << ") * J * pointvalsx.Col(" << j2 << ");\n";
                    s << "hvy = weights(j1+j2) * J * hvx;\n";
                    for (size_t i = 0; i < diffopsy.Size(); i++)
                      {
                        // s << diffopsy[i]->GenerateTransformationCode ("hv", "pointvalsrefy.Col("+ToString(j2)+").Range(0,3)", true);
                        IntRange rref = ranges_yref[i];
                        IntRange r = ranges_y[i];
                        s << diffopsy[i]->GenerateTransformationCode ("hvy.Range("+ToString(r.First())+","+ToString(r.Next())+")",
                                                                      "pointvalsrefy.Col(j2).Range("+ToString(rref.First())+","+ToString(rref.Next())+")",
                                                                      true);
                      }
                    s << "}; // for end BS_IPTS\n";
                  }
              }
            else
              {
                s << "pointvalsrefy = pointvalsrefx;\n";
              }
            
            if (opts.timers) s << "tBt.Start();\n";
            s <<
              "for (size_t j2 = 0; j2 < BS_IPTS; j2++)\n"              
              "for (size_t i = 0; i < " << locdofsy << "; i++) { \n"
              "SIMD<double,BS_ELS> sum = elvecy(i);\n"
              "for (size_t k = 0; k < " << dimyref << "; k++)\n"
              "sum += bmaty[iy++]*pointvalsrefy(k,j2); \n"
              "elvecy(i) = sum;\n"
              "}\n";
            if (opts.timers) s << "tBt.Stop();\n";            
            s << "} // end for nip\n"; 
            
            
            // remainder NIP

            s << "size_t rest = " << nip%opts.BS_ipts << ";\n";
            s << "[[maybe_unused]] size_t base = " << nip - (nip%opts.BS_ipts) << ";\n";
            if (opts.timers) s << "trem.Start();\n";            
            s << "for (size_t j2 = 0; j2 < rest; j2++) {\n"; 
            s <<
              "for (size_t i = 0; i < " << dimxref << "; i++) { \n"
              "SIMD<double,BS_ELS> sum(0.0);\n"
              "for (size_t k = 0; k < " << locdofsx << "; k++)\n"
              "sum += bmatx[ix++]*elvecx(k); \n"
              "pointvalsrefx(i,j2) = sum;\n"
              "} }\n";
            
            if  (!opts.only_loadstoreB)
              {
                s << "Vec<3,SIMD<double,BS_ELS>> hv;\n";
                for (size_t j2 = 0; j2 < opts.BS_ipts; j2++)
                  {
                    for (size_t i = 0; i < diffopsx.Size(); i++)
                      s << diffopsx[i]->GenerateTransformationCode ("pointvalsrefx.Col("+ToString(j2)+").Range(0,3)", "hv", false);
                    s << "hv *= weights(base+"<<j2<<") * J;\n";
                    for (size_t i = 0; i < diffopsy.Size(); i++)
                      s << diffopsy[i]->GenerateTransformationCode ("hv", "pointvalsrefy.Col("+ToString(j2)+").Range(0,3)", true);
                    // s << "pointvalsrefy = pointvalsrefx; \n";
                  }
              }
            else
              {
                s << "pointvalsrefy = pointvalsrefx;\n";                
              }

            s <<
              "for (size_t j2 = 0; j2 < rest; j2++)\n"              
              "for (size_t i = 0; i < " << locdofsy << "; i++) { \n"
              "SIMD<double,BS_ELS> sum(0.0);\n"
              "for (size_t k = 0; k < " << dimyref << "; k++)\n"
              "sum += bmaty[iy++]*pointvalsrefy(k,j2); \n"
              "elvecy(i) += sum;\n"
              "}\n";
            if (opts.timers) s << "trem.Stop();\n";            
            // s << "elvecy *= s;\n";
          }

        else  // only_loadstore
          {
            s << "elvecy = elvecx;\n";
          }
        
        if (opts.timers) s << "tS.Start();\n";            
        s << "for (int k = 0; k < " << locdofsy << "; k++) { \n";
        s << "SIMD<double,BS_ELS> selvecy = s * elvecy(k);\n";
        for (int l = 0; l < opts.BS_els; l++)
          if (opts.atomic)
            s << "AtomicAdd (fy(dofy["<<l<<"][k]), selvecy[" << l << "]);\n";
          else
            s << "fy(dofy["<<l<<"][k]) += selvecy[" << l << "];\n";
        s << "}\n";
        if (opts.timers) s << "tS.Stop();\n";                        

        
        // auto compiledcf = Compile (cf, false);
        // Code code = compiledcf->GenerateProgram(0, false);
        
        s << "} } ); }\n";


        // cout << "code = " << endl << s.str() << endl;
        
        try
          {
            library = CompileCode ( { s.str() }, {}, true );
            compiled_function = library->GetSymbol<lib_function> ("AddBTDTB");
          }
        catch (const Exception & e)
          {
            cout << s.str() << endl;
          } 
      }
    ;
  }
  

  AutoVector MatrixFreeBTDTB :: CreateColVector() const
  {
    return make_unique<VVector<double>> (height);
  }

  AutoVector MatrixFreeBTDTB :: CreateRowVector() const
  {
    return make_unique<VVector<double>> (width);
  }


  void MatrixFreeBTDTB :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("MatrixFreeBTDTB"); RegionTimer reg(t);
    
    // LocalHeap lh(10000);
    auto [locdofsx, dimxref, nip] = Bx.Shape();
    auto [locdofsy, dimyref, nip_] = By.Shape();
    auto [numels,dimy,dimx,nipD] = D.Shape();
    auto [numels2,dimr,dims,nipJ] = Jacobi.Shape();
    
    auto fx = x.FV<double>();
    auto fy = y.FV<double>();
    


    // VER1:
    /*
    static IntegrationPoint dummyip;
    static FE_ElementTransformation<2,2> dummytrafo(ET_TRIG);
    [[maybe_unused]] auto CreateMIP = [](SliceMatrix<> jac, LocalHeap &lh) -> const BaseMappedIntegrationPoint&
    {
      if (jac.Height()==2 && jac.Width()==2)
        return *new (lh) MappedIntegrationPoint<2,2> (dummyip, dummytrafo, Vec<2>(0,0), Mat<2,2>(jac));
      if (jac.Height()==3 && jac.Width()==2)
        return *new (lh) MappedIntegrationPoint<2,3> (dummyip, dummytrafo, Vec<3>(0,0,0), Mat<3,2>(jac));
      if (jac.Height()==3 && jac.Width()==3)
        return *new (lh) MappedIntegrationPoint<3,3> (dummyip, dummytrafo, Vec<3>(0,0,0), Mat<3,3>(jac));
      throw Exception("CreateMIP not implemented for that dim");
    };

    ParallelForRange (Range(elnums), [&](IntRange r) {

      auto &lh = TLHeap();
      Vector<> elvecx(locdofsx);
      Vector<> elvecy(locdofsy);
      
      Vector<> pointvalsrefx(dimxref);
      Vector<> pointvalsrefy(dimyref);
      Vector<> pointvalsx(dimx);
      Vector<> pointvalsy(dimy);
      
      Matrix<> tmaty(dimy,dimyref);
      
      for (auto i : r)
        {
          ElementId ei(VOL, elnums[i]);
          
          elvecx = fx(dofx[i]);
          elvecy = 0.0;
          for (size_t j : Range(nip))
            {
              HeapReset hr(lh);
              
              const auto &mip = CreateMIP (Matrix(Jacobi(i, STAR, STAR, 0)), lh);
              
              pointvalsrefx = Trans(Bx(STAR,STAR,j)) * elvecx;
              pointvalsx = 0;
              for (size_t i : Range(diffopsx))
                {
                  auto &dopx = diffopsx[i];
                  FlatMatrix<> tmatx(ranges_x[i].Size(), ranges_xref[i].Size(), lh);
                  dopx->CalcTransformationMatrix(mip, tmatx, lh);
                  pointvalsx.Range(ranges_x[i]) = tmatx * pointvalsrefx.Range(ranges_xref[i]);
                }
              
              if (nipD==1)
                pointvalsy = weights(j) * D(i, STAR, STAR, 0) * pointvalsx;
              else
                pointvalsy = weights(j) * D(i, STAR, STAR, j) * pointvalsx;
              
              for (size_t i : Range(diffopsy))
                {
                  auto &dopy = diffopsy[i];
                  FlatMatrix<> tmaty(ranges_y[i].Size(), ranges_yref[i].Size(),lh);
                  dopy->CalcTransformationMatrix(mip, tmaty, lh);
                  pointvalsrefy.Range(ranges_yref[i]) = Trans(tmaty) * pointvalsy.Range(ranges_y[i]);
                }
              
              elvecy += By(STAR,STAR,j) * pointvalsrefy;
            }
          
          for (size_t j = 0;j < elvecy.Size(); j++)
            AtomicAdd (fy(dofy[i][j]), s*elvecy(j));
        }
    });
    */

    // VER2
    /*
    static IntegrationPoint dummyip;

    ParallelForRange (Range(elnums), [&](IntRange r) {

      auto &lh = TLHeap();
      Vector<> elvecx(locdofsx);
      Vector<> elvecy(locdofsy);
      
      Vector<> pointvalsrefx(dimxref);
      Vector<> pointvalsrefy(dimyref);
      Vector<> pointvalsx(dimx);
      Vector<> pointvalsy(dimy);
      
      Matrix<> tmaty(dimy,dimyref);
      


      Switch<3> (dimr-1, [&] (auto dimrm1) {        
        Switch<dimrm1+1> (dims-1, [&] (auto dimsm1) {

          static constexpr int dims = dimsm1+1;
          static constexpr int dimr = dimrm1+1;
          
          static FE_ElementTransformation<dims,dimr> dummytrafo(dims==2?ET_TRIG:ET_TET);          
          for (auto i : r)
            {
              ElementId ei(VOL, elnums[i]);

              FlatArray<DofId> dofxi = dofx[i];          
              elvecx = fx(dofxi);
              elvecy = 0.0;


              if (!opts.only_loadstore)
                {
                  for (size_t j : Range(nip))
                    {
                      HeapReset hr(lh);
                      
                      MappedIntegrationPoint<dims,dimr> mip(dummyip, dummytrafo, Vec<dimr>(0), Mat<dimr,dims>(Jacobi(i, STAR, STAR, 0)));
                      
                      pointvalsrefx = Trans(Bx(STAR,STAR,j)) * elvecx;
                      // for (size_t k = 0; k < pointvalsrefx.Size(); k++)
                      // pointvalsrefx(k) = InnerProduct(Bx(STAR,k,j), elvecx);
                      
                      if (!opts.only_loadstoreB)
                        {
                          pointvalsx = 0;
                          for (size_t i : Range(diffopsx))
                            {
                              auto &dopx = diffopsx[i];
                              FlatMatrix<> tmatx(ranges_x[i].Size(), ranges_xref[i].Size(), lh);
                              dopx->CalcTransformationMatrix(mip, tmatx, lh);
                              pointvalsx.Range(ranges_x[i]) = tmatx * pointvalsrefx.Range(ranges_xref[i]);
                            }
                          
                          if (nipD==1)
                            pointvalsy = weights(j) * D(i, STAR, STAR, 0) * pointvalsx;
                          else
                            pointvalsy = weights(j) * D(i, STAR, STAR, j) * pointvalsx;
                          
                          for (size_t i : Range(diffopsy))
                            {
                              auto &dopy = diffopsy[i];
                              FlatMatrix<> tmaty(ranges_y[i].Size(), ranges_yref[i].Size(),lh);
                              dopy->CalcTransformationMatrix(mip, tmaty, lh);
                              pointvalsrefy.Range(ranges_yref[i]) = Trans(tmaty) * pointvalsy.Range(ranges_y[i]);
                            }
                          elvecy += By(STAR,STAR,j) * pointvalsrefy;
                        }
                      else
                        {
                          auto mins = min(pointvalsrefy.Size(), pointvalsrefx.Size());                          
                          pointvalsrefy.Range(mins) = pointvalsrefx.Range(mins);
                        }
                    }
                }
              else
                {
                  auto mins = min(elvecx.Size(), elvecy.Size());
                  elvecy.Range(mins) = elvecx.Range(mins);
                }
              
              FlatArray<DofId> dofyi = dofy[i];
              if (opts.atomic)
                for (size_t j = 0;j < elvecy.Size(); j++)
                  AtomicAdd (fy(dofyi[j]), s*elvecy(j));
              else
                for (size_t j = 0;j < elvecy.Size(); j++)
                  fy(dofyi[j]) += s*elvecy(j);
            }
        });
      });
    });
    */


    // VER3: vectorized over elements

    size_t svec = elnums.Size()/SW;
    size_t rest_base = SW*svec;
    if (compiled_function)
      {
        svec = elnums.Size()/opts.BS_els;
        rest_base = svec*opts.BS_els;
        (*compiled_function) (s, fx, fy, dofx, dofy, Jacobi, weights, elnums.Size());
      }
    else
      {
    
        static IntegrationPoint dummyip;
        
        // static constexpr int SW = 4*SIMD<double>::Size();
        
        // size_t rest = elnums.Size()-SW*svec;
        
        ParallelForRange (Range(svec), [&](IntRange r) {

      auto &lh = TLHeap();
      Vector<SIMD<double,SW>> elvecx(locdofsx);
      Vector<SIMD<double,SW>> elvecy(locdofsy);
      
      Vector<SIMD<double,SW>> pointvalsrefx(dimxref);
      Vector<SIMD<double,SW>> pointvalsrefy(dimyref);
      Vector<SIMD<double,SW>> pointvalsx(dimx);
      Vector<SIMD<double,SW>> pointvalsy(dimy);
      FlatMatrixFixWidth<SW> pointvalsrefx_mat(dimxref, (double*)pointvalsrefx.Data());
      FlatMatrixFixWidth<SW> pointvalsrefy_mat(dimyref, (double*)pointvalsrefy.Data());
      FlatMatrixFixWidth<SW> pointvalsx_mat(dimx, (double*)pointvalsx.Data());
      FlatMatrixFixWidth<SW> pointvalsy_mat(dimy, (double*)pointvalsy.Data());

      
      Matrix<> tmaty(dimy,dimyref);

      Switch<3> (dimr-1, [&] (auto dimrm1) {        
        Switch<dimrm1+1> (dims-1, [&] (auto dimsm1) {

          constexpr int dims = dimsm1.value+1;
          constexpr int dimr = dimrm1.value+1;
          static FE_ElementTransformation<dims,dimr> dummytrafo(dims==2?ET_TRIG:ET_TET);
          
          for (auto i : r)
            {
              ElementId ei(VOL, elnums[i]);

              for (int k = 0; k < dofx[i*SW].Size(); k++)
                elvecx(k) = SIMD<double,SW>{ [&](int nr) { return fx(dofx[i*SW+nr][k]); } };
                  
              elvecy = SIMD<double,SW>(0.0);

              if (!opts.only_loadstore)
                {
                  for (size_t j : Range(nip))
                    {
                      HeapReset hr(lh);
                      
                      for (size_t k = 0; k < pointvalsrefx.Size(); k++)
                        pointvalsrefx(k) = InnerProduct(Bx(STAR,k,j), elvecx);
                      
                      if (!opts.only_loadstoreB)
                        {
                          for (size_t k = 0; k < SW; k++)
                            {
                              size_t ii = i*SW+k;
                                                           
                              MappedIntegrationPoint<dims,dimr> mip(dummyip, dummytrafo, Vec<dimr>(0), Mat<dimr,dims>(Jacobi(ii, STAR, STAR, 0)));
                              pointvalsx = 0;
                              for (size_t i : Range(diffopsx))
                                {
                                  auto &dopx = diffopsx[i];
                                  FlatMatrix<> tmatx(ranges_x[i].Size(), ranges_xref[i].Size(), lh);
                                  dopx->CalcTransformationMatrix(mip, tmatx, lh);
                                  pointvalsx_mat.Col(k).Range(ranges_x[i]) = tmatx * pointvalsrefx_mat.Col(k).Range(ranges_xref[i]);
                                }

                              if (nipD==1)
                                pointvalsy_mat.Col(k) = weights(j) * D(ii, STAR, STAR, 0) * pointvalsx_mat.Col(k);
                              else
                                pointvalsy_mat.Col(k) = weights(j) * D(ii, STAR, STAR, j) * pointvalsx_mat.Col(k);
                              
                              for (size_t i : Range(diffopsy))
                                {
                                  auto &dopy = diffopsy[i];
                                  FlatMatrix<> tmaty(ranges_y[i].Size(), ranges_yref[i].Size(),lh);
                                  dopy->CalcTransformationMatrix(mip, tmaty, lh);
                                  pointvalsrefy_mat.Col(k).Range(ranges_yref[i]) = Trans(tmaty) * pointvalsy_mat.Col(k).Range(ranges_y[i]);
                                }
                            }
                        }
                      else
                        {
                          auto mins = min(pointvalsrefy.Size(), pointvalsrefx.Size());                          
                          pointvalsrefy.Range(mins) = pointvalsrefx.Range(mins);
                          // pointvalsrefy = pointvalsrefx;
                        }

                      
                      // elvecy += By(STAR,STAR,j) * pointvalsrefy;
                      for (size_t k = 0; k < elvecy.Size(); k++)
                        elvecy(k) += InnerProduct(Bx(k,STAR,j), pointvalsrefy);
                    }
                }
              else
                {
                  auto mins = min(elvecx.Size(), elvecy.Size());
                  elvecy.Range(mins) = elvecx.Range(mins);
                }
              
              // FlatArray<DofId> dofyi = dofy[i];
              if (opts.atomic)
                for (size_t j = 0;j < elvecy.Size(); j++)
                  for (size_t k = 0; k < SW; k++)
                    AtomicAdd (fy(dofy[i*SW+k][j]), s*elvecy(j)[k]);
              else
                for (size_t j = 0;j < elvecy.Size(); j++)
                  for (size_t k = 0; k < SW; k++)
                    fy(dofy[i*SW+k][j]) += s*elvecy(j)[k];
            }
        });
      });
    });
      }

    // remainder
    ParallelForRange (Range(rest_base,elnums.Size()), [&](IntRange r) {

      auto &lh = TLHeap();
      static IntegrationPoint dummyip;      
      Vector<> elvecx(locdofsx);
      Vector<> elvecy(locdofsy);
      
      Vector<> pointvalsrefx(dimxref);
      Vector<> pointvalsrefy(dimyref);
      Vector<> pointvalsx(dimx);
      Vector<> pointvalsy(dimy);
      
      Matrix<> tmaty(dimy,dimyref);

      Switch<3> (dimr-1, [&] (auto dimrm1) {        
        Switch<dimrm1.value+1> (dims-1, [&] (auto dimsm1) {

          static constexpr int dims = dimsm1.value+1;
          static constexpr int dimr = dimrm1.value+1;
          
          static FE_ElementTransformation<dims,dimr> dummytrafo(dims==2?ET_TRIG:ET_TET);          
          for (auto i : r)
            {
              ElementId ei(VOL, elnums[i]);

              FlatArray<DofId> dofxi = dofy[i];          
              elvecx = fx(dofxi);
              elvecy = 0.0;


              if (!opts.only_loadstore)
                {
                  for (size_t j : Range(nip))
                    {
                      HeapReset hr(lh);
                      
                      MappedIntegrationPoint<dims,dimr> mip(dummyip, dummytrafo, Vec<dimr>(0), Mat<dimr,dims>(Jacobi(i, STAR, STAR, 0)));
                      
                      pointvalsrefx = Trans(Bx(STAR,STAR,j)) * elvecx;
                      // for (size_t k = 0; k < pointvalsrefx.Size(); k++)
                      // pointvalsrefx(k) = InnerProduct(Bx(STAR,k,j), elvecx);
                      
                      if (!opts.only_loadstoreB)
                        {
                          pointvalsx = 0;
                          for (size_t i : Range(diffopsx))
                            {
                              auto &dopx = diffopsx[i];
                              FlatMatrix<> tmatx(ranges_x[i].Size(), ranges_xref[i].Size(), lh);
                              dopx->CalcTransformationMatrix(mip, tmatx, lh);
                              pointvalsx.Range(ranges_x[i]) = tmatx * pointvalsrefx.Range(ranges_xref[i]);
                            }
                          
                          if (nipD==1)
                            pointvalsy = weights(j) * D(i, STAR, STAR, 0) * pointvalsx;
                          else
                            pointvalsy = weights(j) * D(i, STAR, STAR, j) * pointvalsx;
                          
                          for (size_t i : Range(diffopsy))
                            {
                              auto &dopy = diffopsy[i];
                              FlatMatrix<> tmaty(ranges_y[i].Size(), ranges_yref[i].Size(),lh);
                              dopy->CalcTransformationMatrix(mip, tmaty, lh);
                              pointvalsrefy.Range(ranges_yref[i]) = Trans(tmaty) * pointvalsy.Range(ranges_y[i]);
                            }
                          elvecy += By(STAR,STAR,j) * pointvalsrefy;
                        }
                      else
                        {
                          auto mins = min(pointvalsrefy.Size(), pointvalsrefx.Size());                          
                          pointvalsrefy.Range(mins) = pointvalsrefx.Range(mins);
                        }
                    }
                }
              else
                {
                  auto mins = min(elvecx.Size(), elvecy.Size());
                  elvecy.Range(mins) = elvecx.Range(mins);
                }
              
              FlatArray<DofId> dofyi = dofy[i];
              if (opts.atomic)
                for (size_t j = 0;j < elvecy.Size(); j++)
                  AtomicAdd (fy(dofyi[j]), s*elvecy(j));
              else
                for (size_t j = 0;j < elvecy.Size(); j++)
                  fy(dofyi[j]) += s*elvecy(j);
            }
        });
      });
    });

    
  }


  void BilinearForm :: AssembleBDBFused (LocalHeap & lh, bool linear)
  {
    static Timer t("assemble-BDB"); RegionTimer reg(t);
    
    auto fesx = GetTrialSpace();
    auto fesy = GetTestSpace();
    auto ma = GetMeshAccess();

    Array<short> classnr(ma->GetNE(VOL));
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         bool curved = ma->GetElement(el).is_curved;
         classnr[el.Nr()] = 
           2 * SwitchET<ET_SEGM, ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); })
           + (curved ? 1 : 0);
       });
        
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();
    

    shared_ptr<BaseMatrix> sum;

    
    for (auto part : parts)
      {
        auto bfi = dynamic_pointer_cast<SymbolicBilinearFormIntegrator> (part);

        size_t dimR = ma->GetDimension();
        size_t dimS = dimR - int(bfi->VB());

        
        auto & trialproxies = bfi->TrialProxies();
        auto & testproxies = bfi->TestProxies();
        
        int dimx = 0, dimy = 0;
        for (auto proxy : trialproxies)
          dimx += proxy->Evaluator()->Dim();
        for (auto proxy : testproxies)
          dimy += proxy->Evaluator()->Dim();
        
        int dimxref = 0, dimyref = 0;
        for (auto proxy : trialproxies)
          dimxref += proxy->Evaluator()->DimRef();
        for (auto proxy : testproxies)
          dimyref += proxy->Evaluator()->DimRef();
        

        // cout << "dimx = " << dimx << ", dimxref = " << dimxref << endl;
        // cout << "dimy = " << dimy << ", dimyref = " << dimyref << endl;
        
        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
            ElementId ei(VOL,elclass_inds[0]);

            auto & felx = GetTrialSpace()->GetFE (ei, lh);
            auto & fely = GetTestSpace()->GetFE (ei, lh);
            MixedFiniteElement fel(felx, fely);
            int bonus_intorder = bfi->GetBonusIntegrationOrder();

            bool curved = ma->GetElement(ei).is_curved;
            
            IntegrationRule ir;
            if (bfi->ElementVB() == VOL)
              {
                // const IntegrationRule & volir = bfi->GetIntegrationRule(felx.ElementType(), felx.Order()+fely.Order()+bonus_intorder);
                const IntegrationRule & volir = bfi->GetIntegrationRule(fel, lh);
                for (auto ip : volir)
                  ir += ip;
              }
            else
              {
                auto eltype = felx.ElementType();
                
                Facet2ElementTrafo transform(eltype, bfi->ElementVB()); 
                int nfacet = transform.GetNFacets();
                
                for (int k = 0; k < nfacet; k++)
                  {
                    HeapReset hr(lh);
                    ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
                    IntegrationRule ir_facet(etfacet, felx.Order()+fely.Order()+bonus_intorder);
                    IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
                    for (auto ip : ir_facet_vol)
                      ir += ip;
                  }
              }
            
            Tensor<3> bmatx(felx.GetNDof(), dimxref, ir.Size());
            Tensor<3> bmaty(fely.GetNDof(), dimyref, ir.Size());
            
            for (int i : Range(ir.Size()))
              {
                int starti = 0;
                for (auto proxy : trialproxies)
                  {
                    auto diffopx = proxy->Evaluator();
                    int nexti = starti+diffopx->DimRef();
                    Matrix<double,ColMajor> hbmatx(diffopx->DimRef(), felx.GetNDof());
                    diffopx->CalcMatrix(felx, ir[i], hbmatx, lh);
                    bmatx(STAR,STAR,i).Cols(starti, nexti) = Trans(hbmatx);
                    starti = nexti;
                  }
                
                starti = 0;
                for (auto proxy : testproxies)
                  {
                    auto diffopy = proxy->Evaluator();
                    int nexti = starti+diffopy->DimRef();
                    Matrix<double,ColMajor> hbmaty(diffopy->DimRef(), fely.GetNDof());
                    diffopy->CalcMatrix(fely, ir[i], hbmaty, lh);
                    bmaty(STAR,STAR,i).Cols(starti, nexti) = Trans(hbmaty);
                    starti = nexti;
                  }
              }

            // cout << "bmatx = " << bmatx << endl;
            // cout << "bmaty = " << bmaty << endl;

            Table<DofId> dofx(elclass_inds.Size(), felx.GetNDof());
            Table<DofId> dofy(elclass_inds.Size(), fely.GetNDof());
            
            ParallelForRange (elclass_inds.Range(), [&] (IntRange r)
            {
              Array<DofId> dnumsx, dnumsy;
              for (auto i : r)
                {
                  ElementId ei(VOL, elclass_inds[i]);
                  fesx->GetDofNrs(ei, dnumsx);
                  fesy->GetDofNrs(ei, dnumsy);
                  dofx[i] = dnumsx;
                  dofy[i] = dnumsy;
                }
            });

            // cout << "dofx = " << endl << dofx << endl;
            // cout << "dofy = " << endl << dofy << endl;

            shared_ptr<BaseMatrix> mat;
            
            if (linear)
              {
                Tensor<4> diag(elclass_inds.Size(), dimy, dimx, 1 /* ir.Size()*/); 
                Tensor<4> Jacobi(elclass_inds.Size(), dimR, dimS, curved ? ir.Size() : 1);
                
                // for (auto i : Range(elclass_inds))
                ParallelForRange (Range(elclass_inds), [&] (IntRange r)  {
                  auto &lh = TLHeap();
                  for (auto i : r)
                    {
                      HeapReset hr(lh);
                      ElementId ei(VOL, elclass_inds[i]);
                      auto & trafo = ma->GetTrafo(ei, lh);
                      auto & mir = trafo(ir, lh);
                      if (bfi->ElementVB() != VOL) 
                        mir.ComputeNormalsAndMeasure (fel.ElementType());

                      for (auto j : Range(curved?ir.Size():1))
                        Jacobi(i, STAR, STAR, j) = mir[j].GetJacobian();
                      
                      FlatMatrix<> transx(dimx, dimxref, lh);
                      FlatMatrix<> transy(dimy, dimyref, lh);
                      FlatMatrix<> prod(dimyref, dimxref, lh);
                      
                      shared_ptr<CoefficientFunction> cf = bfi -> GetCoefficientFunction();
                      ProxyUserData ud(trialproxies.Size(), bfi->GridFunctionCoefficients().Size(), lh);
                      for (CoefficientFunction * cf : bfi->GridFunctionCoefficients())
                        ud.AssignMemory (cf, ir.GetNIP(), cf->Dimension(), lh,
                                         cf->IsComplex());
                      
                      const_cast<ElementTransformation&>(trafo).userdata = &ud;
                      
                      
                      FlatMatrix<> val(ir.Size(), 1, lh);
                      
                      {
                        int k1 = 0;
                        for (auto proxy1 : trialproxies)
                          {
                            int l1 = 0;
                            for (auto proxy2 : testproxies)
                              {
                                for (int k = 0; k < proxy1->Dimension(); k++)
                                  for (int l = 0; l < proxy2->Dimension(); l++)
                                    {
                                      ud.trialfunction = proxy1;
                                      ud.trial_comp = k;
                                      ud.testfunction = proxy2;
                                      ud.test_comp = l;
                                      
                                      cf -> Evaluate (mir, val);
                                      // proxyvalues(STAR,l1+l,k1+k) = val.Col(0);
                                      diag(i, l1+l, k1+k, 0) = mir[0].GetMeasure()*val.Col(0)(0);
                                    }
                                l1 += proxy2->Dimension();
                              }
                            k1 += proxy1->Dimension();
                          }
                      }
                    }
                });
                
                
                Vector<> weights(ir.Size());
                for (auto i : Range(ir))
                  weights[i] = ir[i].Weight();
                Array<shared_ptr<DifferentialOperator>> diffopsx, diffopsy;
                for (auto proxy : trialproxies)
                  diffopsx.Append (proxy->Evaluator());
                for (auto proxy : testproxies)
                  diffopsy.Append (proxy->Evaluator());

                // cout << "diag = " << endl << diag << endl;
                // cout << "jac = " << Jacobi << endl;
                
                mat = make_shared<MatrixFreeBTDTB> (fesy->GetNDof(), fesx->GetNDof(),
                                                    Array<size_t>(elclass_inds), std::move(dofx), std::move(dofy),
                                                    std::move(bmatx), std::move(bmaty),
                                                    std::move(weights),
                                                    std::move(diffopsx), std::move(diffopsy), std::move(diag), std::move(Jacobi),
                                                    *matfree_opts);
                                                    
              }
            
#ifdef NOTYET       
            else // linear
              {
                Tensor<3> diagx(dimx, dimxref, nip);
                Tensor<3> diagy(dimy, dimyref, nip);
                Matrix<> points(ma->GetDimension(), nip);
                Matrix<> normals(ma->GetDimension(), nip);

                ParallelForRange
                  (elclass_inds.Size(), [&] (IntRange myrange)
                   {
                     // LocalHeap llh(1000000, "assemble-BDB-D");
                     auto &llh = TLHeap();
                    for (auto i : myrange)
                     {
                    HeapReset hr(llh);
                    ElementId ei(VOL, elclass_inds[i]);
                    auto & trafo = ma->GetTrafo(ei, llh);
                    auto & mir = trafo(ir, llh);
                    if (bfi->ElementVB() != VOL)
                      mir.ComputeNormalsAndMeasure (fel.ElementType());

                    FlatMatrix<> transx(dimx, dimxref, llh);
                    FlatMatrix<> transy(dimy, dimyref, llh);

                    transx = 0.0;
                    transy = 0.0;
                    for (int j = 0; j < ir.Size(); j++)
                      {
                        int starti = 0, startiref = 0;
                        for (auto proxy : trialproxies)
                          {
                            auto diffop = proxy->Evaluator();
                            int nexti = starti+diffop->Dim();
                            int nextiref = startiref+diffop->DimRef();
                            diffop->CalcTransformationMatrix(mir[j], transx.Rows(starti,nexti).Cols(startiref,nextiref), llh);
                            starti = nexti;
                            startiref = nextiref;
                          }
                        starti = 0; startiref = 0;
                        for (auto proxy : testproxies)
                          {
                            auto diffop = proxy->Evaluator();
                            int nexti = starti+diffop->Dim();
                            int nextiref = startiref+diffop->DimRef();
                            diffop->CalcTransformationMatrix(mir[j], transy.Rows(starti,nexti).Cols(startiref,nextiref), llh);
                            starti = nexti;
                            startiref = nextiref;
                          }

                        transy *= mir[j].GetWeight();
                        // diagx(STAR,STAR,i*ir.Size()+j) = transx;  // old
                        // diagy(STAR,STAR,i*ir.Size()+j) = transy;  // old
                        diagx(STAR,STAR,i+j*nel) = transx;
                        diagy(STAR,STAR,i+j*nel) = transy; 
                      }

                    /*
                    points.Cols(i*ir.Size(), (i+1)*ir.Size()) = Trans(mir.GetPoints());
                    normals.Cols(i*ir.Size(), (i+1)*ir.Size()) = Trans(mir.GetNormals());
                    */
                    for (int j = 0; j < ir.Size(); j++)
                      {
                        points.Col(i+j*nel) = mir.GetPoints().Row(j);   // untested
                        normals.Col(i+j*nel) = mir.GetNormals().Row(j); // untested
                      }
                     }
                   });
                shared_ptr<CoefficientFunction> coef = bfi -> GetCoefficientFunction();
                Array<shared_ptr<CoefficientFunction>> diffcfs;
                for (auto proxy : testproxies)
                  {
                    CoefficientFunction::T_DJC cache;
                    diffcfs += coef -> DiffJacobi(proxy, cache);                  
                  }

                auto ipop = make_shared<ApplyIntegrationPoints> (std::move(diffcfs), trialproxies, std::move(points), std::move(normals),
                                                                 dimx, dimy, nip);

                auto & input_coefs = ipop->GetInputCoefs();
                if (input_coefs.Size())
                  {
                    Matrix<double> coef_values(ipop->GetDimCoef(), nip);
                    coef_values = 0.0;
                    ParallelForRange
                      (elclass_inds.Size(), [&] (IntRange myrange)
                       {
                         // LocalHeap llh(1000000, "assemble-BDB-coef");
                        auto &llh = TLHeap();                        
                        for (auto i : myrange)
                         {
                        HeapReset hr(llh);
                        ElementId ei(VOL, elclass_inds[i]);
                        auto & trafo = ma->GetTrafo(ei, llh);
                        auto & mir = trafo(ir, llh);
                        if (bfi->ElementVB() != VOL)
                          mir.ComputeNormalsAndMeasure (fel.ElementType());

                        int off = 0;
                        for (auto icf : input_coefs)
                          {
                            int d = icf->Dimension();
                            FlatMatrix<double> vals(ir.Size(), d, llh);
                            icf->Evaluate (mir, vals);
                            for (int j = 0; j < ir.Size(); j++)
                              for (int c = 0; c < d; c++)
                                coef_values(off+c, i+j*nel) = vals(j,c);
                            off += d;
                          }
                         }
                       });
                    ipop->SetCoefValues (std::move(coef_values));
                  }

                if (ipop->NeedsElementIndex())
                  {
                    Array<int> elidx(nip);
                    for (auto i : Range(elclass_inds))
                      {
                        int di = ma->GetElIndex (ElementId(VOL, elclass_inds[i]));
                        for (int j = 0; j < ir.Size(); j++)
                          elidx[i+j*nel] = di;
                      }
                    ipop->SetElementIndex (std::move(elidx));
                  }

                auto diagmatx = make_shared<BlockDiagonalMatrixSoA> (std::move(diagx));
                auto diagmaty = make_shared<BlockDiagonalMatrixSoA> (std::move(diagy));
                
                mat = TransposeOperator(diagmaty * by) * ipop * (diagmatx * bx);
              } // linear
#endif            
            if (sum)
              sum = sum + mat;
            else
              sum = mat;
          }
      }

    mats.SetSize (ma->GetNLevels());
    mats.Last() = sum;
  };



  void BilinearForm :: AssembleBDB (LocalHeap & lh, bool linear)
  {
    if (auto mf = matfree_opts)
      {
        if (mf->fused)
          {
            AssembleBDBFused (lh, linear);
            return;
          }
      }

    static Timer t("assemble-BDB"); RegionTimer reg(t);
    
    auto fesx = GetTrialSpace();
    auto fesy = GetTestSpace();
    auto ma = GetMeshAccess();

    
    static Timer tclass("assmble classify");
    static Timer tgroup("assmble-group");    
    static Timer tgroup2("assmble-group 2");
    static Timer tgroupD("assmble-group D");
    static Timer tgroupT("assmble-group T");
    static Timer tgroupTi("assmble-group Ti");
    static Timer tgroupTi2("assmble-group Ti2");        
    

    Array<short> classnr(ma->GetNE(VOL));
    tclass.Start();
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_SEGM, ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
       });
    tclass.Stop();
        
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();
    

    shared_ptr<BaseMatrix> sum;

    for (auto part : parts)
      {
        auto bfi = dynamic_pointer_cast<SymbolicBilinearFormIntegrator> (part);

        auto & trialproxies = bfi->TrialProxies();
        auto & testproxies = bfi->TestProxies();
        
        int dimx = 0, dimy = 0;
        for (auto proxy : trialproxies)
          dimx += proxy->Evaluator()->Dim();
        for (auto proxy : testproxies)
          dimy += proxy->Evaluator()->Dim();
        
        int dimxref = 0, dimyref = 0;
        for (auto proxy : trialproxies)
          dimxref += proxy->Evaluator()->DimRef();
        for (auto proxy : testproxies)
          dimyref += proxy->Evaluator()->DimRef();
        


        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
            
            RegionTimer rgroup(tgroup);
            
            ElementId ei(VOL,elclass_inds[0]);
            auto & felx = GetTrialSpace()->GetFE (ei, lh);
            auto & fely = GetTestSpace()->GetFE (ei, lh);
        
            MixedFiniteElement fel(felx, fely);
            int bonus_intorder = bfi->GetBonusIntegrationOrder();
            // const IntegrationRule & ir = bfi->GetIntegrationRule(felx.ElementType(), felx.Order()+fely.Order());
            IntegrationRule ir;
            if (bfi->ElementVB() == VOL)
              {
                const IntegrationRule & volir = bfi->GetIntegrationRule(felx.ElementType(), felx.Order()+fely.Order()+bonus_intorder);
                for (auto ip : volir)
                  ir += ip;
              }
            else
              {
                auto eltype = felx.ElementType();
                
                Facet2ElementTrafo transform(eltype, bfi->ElementVB()); 
                int nfacet = transform.GetNFacets();
                
                for (int k = 0; k < nfacet; k++)
                  {
                    HeapReset hr(lh);
                    ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
                    IntegrationRule ir_facet(etfacet, felx.Order()+fely.Order()+bonus_intorder);
                    IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
                    for (auto ip : ir_facet_vol)
                      ir += ip;
                  }
              }
            
            Matrix<double,ColMajor> bmatx_(ir.Size()*dimxref, felx.GetNDof());
            Matrix<double,ColMajor> bmaty_(ir.Size()*dimyref, fely.GetNDof());

            for (int i : Range(ir.Size()))
              {
                int starti = i*dimxref;
                for (auto proxy : trialproxies)
                  {
                    auto diffopx = proxy->Evaluator();
                    int nexti = starti+diffopx->DimRef();
                    diffopx->CalcMatrix(felx, ir[i], bmatx_.Rows(starti, nexti), lh);
                    starti = nexti;
                  }
                
                starti = i*dimyref;
                for (auto proxy : testproxies)
                  {
                    auto diffopy = proxy->Evaluator();
                    int nexti = starti+diffopy->DimRef();
                    diffopy->CalcMatrix(fely, ir[i], bmaty_.Rows(starti, nexti), lh);
                    starti = nexti;
                  }
              }

            Matrix bmatx = bmatx_;
            Matrix bmaty = bmaty_;

            // if (!linear)  // transpose intpnt <--> comp
              {
                for (int i : Range(ir.Size()))
                  for (int j : Range(dimxref))
                    bmatx.Row(j*ir.Size()+i) = bmatx_.Row(i*dimxref+j);
                for (int i : Range(ir.Size()))
                  for (int j : Range(dimyref))
                    bmaty.Row(j*ir.Size()+i) = bmaty_.Row(i*dimyref+j);
              }
            

              static Timer tebe("setup ebe mats");
              static Timer tfillidx("fill idx");

            Table<DofId> xdofsin(elclass_inds.Size(), felx.GetNDof());
            Table<DofId> xdofsout(elclass_inds.Size(), bmatx.Height());

            Table<DofId> ydofsin(elclass_inds.Size(), fely.GetNDof());
            Table<DofId> ydofsout(elclass_inds.Size(), bmaty.Height());

            tebe.Start();            

            ParallelForRange (elclass_inds.Range(), [&] (IntRange r)
            {
              Array<DofId> dnumsx, dnumsy;
              for (auto i : r)
                {
                  ElementId ei(VOL, elclass_inds[i]);
                  fesx->GetDofNrs(ei, dnumsx);
                  fesy->GetDofNrs(ei, dnumsy);
                  xdofsin[i] = dnumsx;
                  ydofsin[i] = dnumsy;
                }
            });

            /*
            Array<DofId> dnumsx, dnumsy;
            for (auto i : Range(elclass_inds))
              {
                ElementId ei(VOL, elclass_inds[i]);
                fesx->GetDofNrs(ei, dnumsx);
                fesy->GetDofNrs(ei, dnumsy);
                xdofsin[i] = dnumsx;
                ydofsin[i] = dnumsy;
              }
            */
            
            tebe.Stop();

            
            int nel = elclass_inds.Size();
            size_t nip = ir.Size()*nel;
            
            tfillidx.Start();
            for (size_t i = 0; i < nel; i++)
              for (size_t k = 0; k < dimxref*ir.Size(); k++)
                xdofsout[i][k] = i+k*nel;
            
            for (size_t i = 0; i < nel; i++)
              for (size_t k = 0; k < dimyref*ir.Size(); k++)
                ydofsout[i][k] = i+k*nel;
            tfillidx.Stop();
            
            auto bx = make_shared<ConstantElementByElementMatrix<>>
              (nip*dimxref, fesx->GetNDof(),
               bmatx, std::move(xdofsout), std::move(xdofsin));
            
            auto by = make_shared<ConstantElementByElementMatrix<>>
              (nip*dimyref, fesy->GetNDof(),
               bmaty, std::move(ydofsout), std::move(ydofsin));

            shared_ptr<BaseMatrix> mat;
            
            if (linear)
              {
                // RegionTimer rgroup2(tgroup2);
                
                Tensor<3> diag(dimyref, dimxref, elclass_inds.Size()*ir.Size());
                
                // for (auto i : Range(elclass_inds))
                ParallelForRange (Range(elclass_inds), [&] (IntRange r)  {
                  auto &lh = TLHeap();
                  for (auto i : r)
                  {
                    HeapReset hr(lh);
                    ElementId ei(VOL, elclass_inds[i]);
                    auto & trafo = ma->GetTrafo(ei, lh);
                    auto & mir = trafo(ir, lh);
                    if (bfi->ElementVB() != VOL) 
                      mir.ComputeNormalsAndMeasure (fel.ElementType());
                    FlatMatrix<> transx(dimx, dimxref, lh);
                    FlatMatrix<> transy(dimy, dimyref, lh);
                    FlatMatrix<> prod(dimyref, dimxref, lh);

                    shared_ptr<CoefficientFunction> cf = bfi -> GetCoefficientFunction();
                    ProxyUserData ud(trialproxies.Size(), bfi->GridFunctionCoefficients().Size(), lh);
                    for (CoefficientFunction * cf : bfi->GridFunctionCoefficients())
                      ud.AssignMemory (cf, ir.GetNIP(), cf->Dimension(), lh,
                                       cf->IsComplex());
                    
                    const_cast<ElementTransformation&>(trafo).userdata = &ud;
                    
                    FlatTensor<3> proxyvalues(lh, ir.Size(), dimy, dimx);
                    FlatMatrix<> val(ir.Size(), 1, lh);

                    {
                      // RegionTimer r(tgroupD);
                    int k1 = 0;
                    for (auto proxy1 : trialproxies)
                      {
                        int l1 = 0;
                        for (auto proxy2 : testproxies)
                          {
                            for (int k = 0; k < proxy1->Dimension(); k++)
                              for (int l = 0; l < proxy2->Dimension(); l++)
                                {
                                  ud.trialfunction = proxy1;
                                  ud.trial_comp = k;
                                  ud.testfunction = proxy2;
                                  ud.test_comp = l;
                                  
                                  cf -> Evaluate (mir, val);
                                  proxyvalues(STAR,l1+l,k1+k) = val.Col(0);
                                }
                            l1 += proxy2->Dimension();
                          }
                        k1 += proxy1->Dimension();
                      }
                    }

                    // RegionTimer rT(tgroupT);
                    FlatMatrix<> prod1(dimy, transx.Width(), lh);
                    
                    for (int j = 0; j < ir.Size(); j++)
                      {
                        /*
                        auto diffopx = trialproxies[0]->Evaluator();
                        auto diffopy = testproxies[0]->Evaluator();
                        diffopx->CalcTransformationMatrix(mir[j], transx, lh);
                        diffopy->CalcTransformationMatrix(mir[j], transy, lh);
                        */
                        transx = 0.0;
                        transy = 0.0;

                        int starti = 0, startiref = 0;
                        {
                          // RegionTimer rT(tgroupTi);                        
                        for (auto proxy : trialproxies)
                          {
                            auto &diffop = proxy->Evaluator();
                            int nexti = starti+diffop->Dim();
                            int nextiref = startiref+diffop->DimRef();
                            diffop->CalcTransformationMatrix(mir[j], transx.Rows(starti,nexti).Cols(startiref,nextiref), lh);
                            starti = nexti;
                            startiref = nextiref;
                          }
                        }

                        starti = 0; startiref = 0;
                        {
                          // RegionTimer rT(tgroupTi2);                        
                        for (auto proxy : testproxies)
                          {
                            auto &diffop = proxy->Evaluator();
                            int nexti = starti+diffop->Dim();
                            int nextiref = startiref+diffop->DimRef();
                            diffop->CalcTransformationMatrix(mir[j], transy.Rows(starti,nexti).Cols(startiref,nextiref), lh);
                            starti = nexti;
                            startiref = nextiref;
                          }
                        }
                        // prod = Trans(transy) * proxyvalues(j,STAR,STAR) * transx;
                        prod1 = proxyvalues(j,STAR,STAR) * transx;
                        prod = Trans(transy)*prod1;
                        
                        prod *= mir[j].GetWeight();
                        diag(STAR,STAR,i+j*nel) = prod;
                      }
                  }
                });
                auto diagmat = make_shared<BlockDiagonalMatrixSoA> (std::move(diag));
                mat = TransposeOperator(by) * diagmat * bx;
              }
            else // linear
              {
                Tensor<3> diagx(dimx, dimxref, nip);
                Tensor<3> diagy(dimy, dimyref, nip);
                Matrix<> points(ma->GetDimension(), nip);
                Matrix<> normals(ma->GetDimension(), nip);

                ParallelForRange
                  (elclass_inds.Size(), [&] (IntRange myrange)
                   {
                     // LocalHeap llh(1000000, "assemble-BDB-D");
                     auto &llh = TLHeap();
                    for (auto i : myrange)
                     {
                    HeapReset hr(llh);
                    ElementId ei(VOL, elclass_inds[i]);
                    auto & trafo = ma->GetTrafo(ei, llh);
                    auto & mir = trafo(ir, llh);
                    if (bfi->ElementVB() != VOL)
                      mir.ComputeNormalsAndMeasure (fel.ElementType());

                    FlatMatrix<> transx(dimx, dimxref, llh);
                    FlatMatrix<> transy(dimy, dimyref, llh);

                    transx = 0.0;
                    transy = 0.0;
                    for (int j = 0; j < ir.Size(); j++)
                      {
                        int starti = 0, startiref = 0;
                        for (auto proxy : trialproxies)
                          {
                            auto diffop = proxy->Evaluator();
                            int nexti = starti+diffop->Dim();
                            int nextiref = startiref+diffop->DimRef();
                            diffop->CalcTransformationMatrix(mir[j], transx.Rows(starti,nexti).Cols(startiref,nextiref), llh);
                            starti = nexti;
                            startiref = nextiref;
                          }
                        starti = 0; startiref = 0;
                        for (auto proxy : testproxies)
                          {
                            auto diffop = proxy->Evaluator();
                            int nexti = starti+diffop->Dim();
                            int nextiref = startiref+diffop->DimRef();
                            diffop->CalcTransformationMatrix(mir[j], transy.Rows(starti,nexti).Cols(startiref,nextiref), llh);
                            starti = nexti;
                            startiref = nextiref;
                          }

                        transy *= mir[j].GetWeight();
                        // diagx(STAR,STAR,i*ir.Size()+j) = transx;  // old
                        // diagy(STAR,STAR,i*ir.Size()+j) = transy;  // old
                        diagx(STAR,STAR,i+j*nel) = transx;
                        diagy(STAR,STAR,i+j*nel) = transy; 
                      }

                    /*
                    points.Cols(i*ir.Size(), (i+1)*ir.Size()) = Trans(mir.GetPoints());
                    normals.Cols(i*ir.Size(), (i+1)*ir.Size()) = Trans(mir.GetNormals());
                    */
                    for (int j = 0; j < ir.Size(); j++)
                      {
                        points.Col(i+j*nel) = mir.GetPoints().Row(j);   // untested
                        normals.Col(i+j*nel) = mir.GetNormals().Row(j); // untested
                      }
                     }
                   });
                shared_ptr<CoefficientFunction> coef = bfi -> GetCoefficientFunction();
                Array<shared_ptr<CoefficientFunction>> diffcfs;
                for (auto proxy : testproxies)
                  {
                    CoefficientFunction::T_DJC cache;
                    diffcfs += coef -> DiffJacobi(proxy, cache);                  
                  }

                auto ipop = make_shared<ApplyIntegrationPoints> (std::move(diffcfs), trialproxies, std::move(points), std::move(normals),
                                                                 dimx, dimy, nip);

                auto & input_coefs = ipop->GetInputCoefs();
                if (input_coefs.Size())
                  {
                    Matrix<double> coef_values(ipop->GetDimCoef(), nip);
                    coef_values = 0.0;
                    ParallelForRange
                      (elclass_inds.Size(), [&] (IntRange myrange)
                       {
                         // LocalHeap llh(1000000, "assemble-BDB-coef");
                        auto &llh = TLHeap();                        
                        for (auto i : myrange)
                         {
                        HeapReset hr(llh);
                        ElementId ei(VOL, elclass_inds[i]);
                        auto & trafo = ma->GetTrafo(ei, llh);
                        auto & mir = trafo(ir, llh);
                        if (bfi->ElementVB() != VOL)
                          mir.ComputeNormalsAndMeasure (fel.ElementType());

                        int off = 0;
                        for (auto icf : input_coefs)
                          {
                            int d = icf->Dimension();
                            FlatMatrix<double> vals(ir.Size(), d, llh);
                            icf->Evaluate (mir, vals);
                            for (int j = 0; j < ir.Size(); j++)
                              for (int c = 0; c < d; c++)
                                coef_values(off+c, i+j*nel) = vals(j,c);
                            off += d;
                          }
                         }
                       });
                    ipop->SetCoefValues (std::move(coef_values));
                  }

                if (ipop->NeedsElementIndex())
                  {
                    Array<int> elidx(nip);
                    for (auto i : Range(elclass_inds))
                      {
                        int di = ma->GetElIndex (ElementId(VOL, elclass_inds[i]));
                        for (int j = 0; j < ir.Size(); j++)
                          elidx[i+j*nel] = di;
                      }
                    ipop->SetElementIndex (std::move(elidx));
                  }

                auto diagmatx = make_shared<BlockDiagonalMatrixSoA> (std::move(diagx));
                auto diagmaty = make_shared<BlockDiagonalMatrixSoA> (std::move(diagy));
                
                mat = TransposeOperator(diagmaty * by) * ipop * (diagmatx * bx);
              } // linear
            
            if (sum)
              sum = sum + mat;
            else
              sum = mat;
          }
      }
    mats.SetSize (ma->GetNLevels());
    mats.Last() = sum;
  }
  
  
}
