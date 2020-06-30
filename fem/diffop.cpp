/*********************************************************************/
/* File:   diffop.cpp                                                */
/* Author: Start                                                     */
/* Date:   24. Nov. 2009                                             */
/*********************************************************************/


 
#include <fem.hpp>
#include "diffop_impl.hpp"



namespace ngfem
{


  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              SliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    throw Exception (string("Error: DifferentialOperator::CalcMatrix called for base class, type = ")
                     + typeid(*this).name());
  }

  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              SliceMatrix<Complex,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    throw Exception (string("Error: DifferentialOperator::CalcMatrix<Complex> called for base class, type = ")
                     + typeid(*this).name());
  }

  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              SliceMatrix<double,ColMajor> mat,   
              LocalHeap & lh) const
  {
    int dim = Dim(); 
    for (int i = 0; i < mir.Size(); i++)
      CalcMatrix (fel, mir[i], mat.Rows(i*dim, (i+1)*dim), lh);
  }

  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              SliceMatrix<Complex,ColMajor> mat,   
              LocalHeap & lh) const
  {
    int dim = Dim(); 
    for (int i = 0; i < mir.Size(); i++)
      CalcMatrix (fel, mir[i], mat.Rows(i*dim, (i+1)*dim), lh);
  }
  
  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const SIMD_BaseMappedIntegrationRule & mir,
              BareSliceMatrix<SIMD<double>> mat) const
  {
    throw ExceptionNOSIMD(string("Error: DifferentialOperator::CalcMatrix does not support SIMD, type = ")
                          + typeid(*this).name());
    
  }
  
  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const SIMD_BaseMappedIntegrationRule & mir,
              BareSliceMatrix<SIMD<Complex>> mat) const
  {
    throw ExceptionNOSIMD(string("Error: DifferentialOperator::CalcMatrix does not support SIMD<Complex>, type = ")
                          + typeid(*this).name());
    
  }

  void DifferentialOperator ::
  CalcLinearizedMatrix (const FiniteElement & fel,
                        const BaseMappedIntegrationRule & mir,
                        BareSliceVector<double> x,
                        SliceMatrix<double,ColMajor> mat,   
                        LocalHeap & lh) const
  {
    CalcMatrix (fel, mir, mat, lh);
  }

  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         BareSliceVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
#ifndef FASTCOMPILE    
    static int cnt = 0;
    if (cnt < 3)
      {
        cnt++;
        cout << "called base class apply, type = " << typeid(*this).name() << endl;
      }
#endif
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), fel.GetNDof(), lh);
    CalcMatrix (fel, mip, mat, lh);
    flux = mat * x;
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         BareSliceVector<Complex> x, 
         FlatVector<Complex> flux,
         LocalHeap & lh) const
  {
#ifndef FASTCOMPILE
    static int cnt = 0;
    if (cnt < 3)
      {
        cnt++;
        cout << "called base class apply, complex, type = " << typeid(*this).name() << endl;
      }
#endif
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), fel.GetNDof(), lh);
    CalcMatrix (fel, mip, mat, lh);
    flux = mat * x;
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationRule & mir,
         BareSliceVector<double> x, 
         BareSliceMatrix<double> flux,
         LocalHeap & lh) const
  {
    for (int i = 0; i < mir.Size(); i++)
      Apply (fel, mir[i], x, flux.Row(i).Range(0,dim), lh);
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationRule & mir,
         BareSliceVector<Complex> x, 
         BareSliceMatrix<Complex> flux,
         LocalHeap & lh) const
  {
    for (int i = 0; i < mir.Size(); i++)
      Apply (fel, mir[i], x, flux.Row(i).Range(0,dim), lh);
  }

  void DifferentialOperator ::
  Apply (const FiniteElement & bfel,
         const SIMD_BaseMappedIntegrationRule & bmir,
         BareSliceVector<double> x, 
         BareSliceMatrix<SIMD<double>> flux) const
  {
    throw ExceptionNOSIMD (string("DifferentialOperator :: Apply ( ... SIMD ... ) not overloaded for class ")
                           + typeid(*this).name());
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & bfel,
         const SIMD_BaseMappedIntegrationRule & bmir,
         BareSliceVector<Complex> x, 
         BareSliceMatrix<SIMD<Complex>> flux) const
  {
    throw ExceptionNOSIMD (string("DifferentialOperator :: Apply ( ... SIMD<Complex> ... ) not overloaded for class ")
                           + typeid(*this).name());
  }
  
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              BareSliceVector<double> x, 
              LocalHeap & lh) const 
  {
#ifndef FASTCOMPILE
    static int cnt = 0;
    if (cnt < 3)
      {
        cnt++;    
        cout << "called base class apply trans, type = " << typeid(*this).name() << endl;
      }
#endif
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), fel.GetNDof(), lh);
    CalcMatrix (fel, mip, mat, lh);
    x.Range(0,fel.GetNDof()) = Trans(mat) * flux;
  }
  
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<Complex> flux,
              BareSliceVector<Complex> x, 
              LocalHeap & lh) const 
  {
#ifndef FASTCOMPILE
    static int cnt = 0;
    if (cnt < 3)
      {
        cnt++;    
        cout << "called base class apply trans complex, type = " << typeid(*this).name() << endl;
      }
#endif
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), fel.GetNDof(), lh);
    CalcMatrix (fel, mip, mat, lh);
    x.Range(0,fel.GetNDof()) = Trans(mat) * flux;
  }
  
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              FlatMatrix<double> flux,
              BareSliceVector<double> x, 
              LocalHeap & lh) const 
  {
    HeapReset hr(lh);
    size_t nd = fel.GetNDof();
    FlatVector<double> hx(nd, lh);
    x.Range(0,nd) = 0.0;
    for (int i = 0; i < mir.Size(); i++)
      {
        ApplyTrans (fel, mir[i], flux.Row(i), hx, lh);
        x.Range(0,nd) += hx;
      }
  }
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              FlatMatrix<Complex> flux,
              BareSliceVector<Complex> x, 
              LocalHeap & lh) const 
  {
    HeapReset hr(lh);
    size_t nd = fel.GetNDof();
    FlatVector<Complex> hx(nd, lh);
    x.Range(0,nd) = 0.0;
    for (int i = 0; i < mir.Size(); i++)
      {
        ApplyTrans (fel, mir[i], flux.Row(i), hx, lh);
        x.Range(0,nd) += hx;
      }
  }

  void DifferentialOperator ::
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & bmir,
            BareSliceMatrix<SIMD<double>> flux,
            BareSliceVector<double> x) const
  {
    throw ExceptionNOSIMD (string("DifferentialOperator :: AddTrans ( ... SIMD ... ) not overloaded for class ") + typeid(*this).name());
  }

  void DifferentialOperator ::
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & bmir,
            BareSliceMatrix<SIMD<Complex>> flux,
            BareSliceVector<Complex> x) const
  {
    throw ExceptionNOSIMD (string("DifferentialOperator :: AddTrans ( ... SIMD<Complex> ... ) not overloaded for class ") + typeid(*this).name());
  }

  list<tuple<string,double>> DifferentialOperator ::
  Timing (const FiniteElement & fel, const BaseMappedIntegrationRule & mir) const
  {
    list<tuple<string,double>> timings;
    LocalHeap lh(100000000);
    Matrix<double,ColMajor> bmat(mir.Size()*Dim(), fel.GetNDof());
    Vector<> coefs(fel.GetNDof());
    Matrix<> values(mir.Size(), Dim());
    
    SIMD_IntegrationRule simd_ir(mir.IR());
    auto & simd_mir = mir.GetTransformation()(simd_ir, lh);
    Matrix<SIMD<double>> simd_bmat(Dim()*fel.GetNDof(), simd_mir.Size());
    Matrix<SIMD<double>> simd_values(Dim(), mir.Size());
    
    double time;
    constexpr size_t steps = 1000;

    time = RunTiming([&]() {
        for (size_t i = 0; i < steps; i++)
          this -> CalcMatrix (fel, mir, bmat, lh);
      });
    timings.push_back(make_tuple("CalcBMatrix", time/steps*1e9/(bmat.Height()*bmat.Width())));
    coefs = 1;
    time = RunTiming([&]() {
        for (size_t i = 0; i < steps; i++)
          this -> Apply (fel, mir, coefs, values, lh);
      });
    timings.push_back(make_tuple("Appy", time/steps*1e9/(bmat.Height()*bmat.Width())));
    values = 1;
    time = RunTiming([&]() {
        for (size_t i = 0; i < steps; i++)
          this -> ApplyTrans (fel, mir, values, coefs, lh);
      });
    timings.push_back(make_tuple("AppyTrans", time/steps*1e9/(bmat.Height()*bmat.Width())));
    try
      {
        time = RunTiming([&]() {
            for (size_t i = 0; i < steps; i++)
              this -> CalcMatrix (fel, simd_mir, simd_bmat);
          });
        timings.push_back(make_tuple("SIMD - CalcBMatrix", time/steps*1e9/(bmat.Height()*bmat.Width())));
      }
    catch (ExceptionNOSIMD e) { ; } 


    coefs = 1;
    try
      {
        time = RunTiming([&]() {
            for (size_t i = 0; i < steps; i++)
              this -> Apply (fel, simd_mir, coefs, simd_values);
          });
        timings.push_back(make_tuple("SIMD - Appy", time/steps*1e9/(bmat.Height()*bmat.Width())));
      }
    catch (ExceptionNOSIMD e) { ; } 

    simd_values = SIMD<double> (1.0);
    coefs = 0.0;
    try
      {
        time = RunTiming([&]() {
            for (size_t i = 0; i < steps; i++)
              this -> AddTrans (fel, simd_mir, simd_values, coefs);
          });
        timings.push_back(make_tuple("SIMD - AppyTrans", time/steps*1e9/(bmat.Height()*bmat.Width())));
      }
    catch (ExceptionNOSIMD e) { ; } 



    
    return timings;
  }


  BlockDifferentialOperator :: ~BlockDifferentialOperator ()  { ; }


  void BlockDifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              SliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat1(diffop->Dim(), fel.GetNDof(), lh);
    diffop->CalcMatrix (fel, mip, mat1, lh);
    mat = 0;
    
    if (comp == -1)
      for (int i = 0; i < mat1.Height(); i++)
        for (int j = 0; j < mat1.Width(); j++)
          for (int k = 0; k < dim; k++)
            mat(dim*i+k, dim*j+k) = mat1(i,j);
    else
      for (int i = 0; i < mat1.Height(); i++)
        for (int j = 0; j < mat1.Width(); j++)
          mat(dim*i+comp, dim*j+comp) = mat1(i,j);
  }
  
  void BlockDifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const SIMD_BaseMappedIntegrationRule & mir,
              BareSliceMatrix<SIMD<double>> mat) const
  {
    diffop->CalcMatrix(fel, mir, mat.RowSlice(0, dim*dim));
    
    size_t hdim = dim;   // how many copies
    size_t dim_diffop = diffop->Dim();
    size_t hdim2 = hdim*hdim;
    size_t dim_dim_do = hdim*dim_diffop;
    size_t dim2_dim_do = hdim2*dim_diffop;

    STACK_ARRAY(SIMD<double>, hval, dim_diffop);
    
    size_t nip = mir.Size();
    if (comp == -1) 
      for (size_t i = 0; i < fel.GetNDof(); i++)
        {
          auto mati = mat.Rows(dim2_dim_do*IntRange(i,i+1));
          for (size_t j = 0; j < nip; j++)
            {
              auto col = mati.Col(j);

              for (size_t l = 0; l < dim_diffop; l++)
                hval[l] = col(l*hdim2);
              
              col.Range(0,dim2_dim_do) = 0;
            
              for (size_t l = 0; l < dim_diffop; l++)
                col.Slice(l*hdim, dim_dim_do+1).Range(0,hdim) = hval[l];
            }
        }
    else
      throw ExceptionNOSIMD("BlockDifferentialOperator::CalcMatrix does not support SIMD");
  }
  

  void BlockDifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         BareSliceVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    HeapReset hr(lh);
    // FlatVector<> hx(fel.GetNDof(), lh);
    FlatVector<> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            // hx = x.Slice(k, dim);
            diffop->Apply(fel, mip, x.Slice(k,dim), hflux, lh);
            flux.Slice(k,dim) = hflux;
          }
      }
    else
      {
        // hx = x.Slice(comp, dim);
        diffop->Apply(fel, mip, x.Slice(comp,dim), hflux, lh);
        flux.Slice(comp,dim) = hflux;
      }
  }


  void BlockDifferentialOperator ::
  Apply (const FiniteElement & fel,
         const SIMD_BaseMappedIntegrationRule & mir,
         BareSliceVector<double> x, 
         BareSliceMatrix<SIMD<double>> flux) const
  // LocalHeap & lh) const
  {
    if (comp == -1)
      for (int k = 0; k < dim; k++)
        diffop->Apply(fel, mir, x.Slice(k, dim), flux.RowSlice(k,dim));
    else
      diffop->Apply(fel, mir, x.Slice(comp, dim), flux.RowSlice(comp,dim));
  }
  
  
  void BlockDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              BareSliceVector<double> x, 
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    // auto x = bx.Range(0,dim*fel.GetNDof());
    // FlatVector<> hx(fel.GetNDof(), lh);
    FlatVector<> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            hflux = flux.Slice(k, dim);
            diffop->ApplyTrans(fel, mip, hflux, x.Slice(k,dim), lh);
            // x.Slice(k,dim) = hx;
          }
      }
    else
      {
        x.Range(0,dim*fel.GetNDof()) = 0.0;
        hflux = flux.Slice(comp, dim);
        diffop->ApplyTrans(fel, mip, hflux, x.Slice(comp,dim), lh);
        // x.Slice(comp,dim) = hx;
      }
  }

  void BlockDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<Complex> flux,
              BareSliceVector<Complex> x, 
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    // auto x = bx;
    // FlatVector<Complex> hx(fel.GetNDof(), lh);
    FlatVector<Complex> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            hflux = flux.Slice(k, dim);
            diffop->ApplyTrans(fel, mip, hflux, x.Slice(k,dim), lh);
            // x.Slice(k,dim) = hx;
          }
      }
    else
      {
        hflux = flux.Slice(comp, dim);
        x.Range(0,dim*fel.GetNDof()) = 0.0;
        diffop->ApplyTrans(fel, mip, hflux, x.Slice(comp,dim), lh);
        // x.Slice(comp,dim) = hx;
      }
  }

  void BlockDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
	      FlatMatrix<double> flux,
	      BareSliceVector<double> x, 
	      LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<double> hflux(flux.Height(), diffop->Dim(), lh);
    for (auto k : (comp == -1) ? Range(0, dim) : Range(comp, comp+1)) {
      for (auto l : Range(diffop->Dim()))
	{ hflux.Col(l) = flux.Col(dim * l + k); }
      diffop->ApplyTrans(fel, mir, hflux, x.Slice(k,dim), lh);
    }
  }

  void BlockDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
	      FlatMatrix<Complex> flux,
	      BareSliceVector<Complex> x, 
	      LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<Complex> hflux(flux.Height(), diffop->Dim(), lh);
    for (auto k : (comp == -1) ? Range(0, dim) : Range(comp, comp+1)) {
      for (auto l : Range(diffop->Dim()))
	{ hflux.Col(l) = flux.Col(dim * l + k); }
      diffop->ApplyTrans(fel, mir, hflux, x.Slice(k,dim), lh);
    }
  }

  void BlockDifferentialOperator ::
  AddTrans (const FiniteElement & fel,
            const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<double>> flux,
            BareSliceVector<double> x) const
  {
    if (comp == -1)
      for (size_t k = 0; k < dim; k++)
        diffop->AddTrans(fel, mir, flux.RowSlice(k,dim), x.Slice(k,dim));
    else
      diffop->AddTrans(fel, mir, flux.RowSlice(comp,dim), x.Slice(comp,dim));
  }

  void BlockDifferentialOperator ::
  AddTrans (const FiniteElement & fel,
            const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<Complex>> flux,
            BareSliceVector<Complex> x) const
  {
    if (comp == -1)
      for (size_t k = 0; k < dim; k++)
        diffop->AddTrans(fel, mir, flux.RowSlice(k,dim), x.Slice(k,dim));
    else
      diffop->AddTrans(fel, mir, flux.RowSlice(comp,dim), x.Slice(comp,dim));
  }


  shared_ptr<CoefficientFunction> BlockDifferentialOperator ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir) const 
  {
    // assume it's a grad ...
    // return (-1)*dir->Operator("grad") * proxy;
    return diffop->DiffShape(proxy, dir);
  }








  BlockDifferentialOperatorTrans :: ~BlockDifferentialOperatorTrans ()  { ; }


  void BlockDifferentialOperatorTrans ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              SliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat1(diffop->Dim(), fel.GetNDof(), lh);
    diffop->CalcMatrix (fel, mip, mat1, lh);
    mat = 0;
    
    if (comp == -1)
      for (int i = 0; i < mat1.Height(); i++)
        for (int j = 0; j < mat1.Width(); j++)
          for (int k = 0; k < dim; k++)
            // mat(dim*i+k, dim*j+k) = mat1(i,j);
            mat(k*mat1.Height()+i, dim*j+k) = mat1(i,j);
    else
      for (int i = 0; i < mat1.Height(); i++)
        for (int j = 0; j < mat1.Width(); j++)
          // mat(dim*i+comp, dim*j+comp) = mat1(i,j);
          mat(comp*mat1.Height()+i, dim*j+comp) = mat1(i,j);
  }
  
  void BlockDifferentialOperatorTrans ::
  CalcMatrix (const FiniteElement & fel,
              const SIMD_BaseMappedIntegrationRule & mir,
              BareSliceMatrix<SIMD<double>> mat) const
  {
    diffop->CalcMatrix(fel, mir, mat.RowSlice(0, dim*dim));
    
    size_t hdim = dim;   // how many copies
    size_t dim_diffop = diffop->Dim();
    size_t hdim2 = hdim*hdim;
    // size_t dim_dim_do = hdim*dim_diffop;
    size_t dim2_dim_do = hdim2*dim_diffop;

    STACK_ARRAY(SIMD<double>, hval, dim_diffop);

    size_t nip = mir.Size();
    if (comp == -1) 
      for (size_t i = 0; i < fel.GetNDof(); i++)
        {
          auto mati = mat.Rows(dim2_dim_do*IntRange(i,i+1));
          for (size_t j = 0; j < nip; j++)
            {
              auto col = mati.Col(j);

              for (size_t l = 0; l < dim_diffop; l++)
                hval[l] = col(l*hdim2);
              
              col.Range(0,dim2_dim_do) = 0;
              for (size_t l = 0; l < dim_diffop; l++)
                col.Slice(l, dim_diffop*(dim+1)).Range(0,hdim) = hval[l];
            }
        }
    else
      throw ExceptionNOSIMD("BlockDifferentialOperatorTrans::CalcMatrix does not support SIMD");
  }
  

  void BlockDifferentialOperatorTrans ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         BareSliceVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    // HeapReset hr(lh);
    // FlatVector<> hx(fel.GetNDof(), lh);
    // FlatVector<> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            // hx = x.Slice(k, dim);
            diffop->Apply(fel, mip, x.Slice(k,dim), flux.Range(diffop->Dim()*IntRange(k,k+1)), lh);
            // flux.Slice(k,dim) = hflux;
          }
      }
    else
      {
        // hx = x.Slice(comp, dim);
        diffop->Apply(fel, mip, x.Slice(comp,dim), flux.Range(diffop->Dim()*IntRange(comp,comp+1)), lh);
        // flux.Slice(comp,dim) = hflux;
      }
  }


  void BlockDifferentialOperatorTrans ::
  Apply (const FiniteElement & fel,
         const SIMD_BaseMappedIntegrationRule & mir,
         BareSliceVector<double> x, 
         BareSliceMatrix<SIMD<double>> flux) const
  {
    if (comp == -1)
      for (int k = 0; k < dim; k++)
        // diffop->Apply(fel, mir, x.Slice(k, dim), flux.RowSlice(k,dim));
        diffop->Apply(fel, mir, x.Slice(k, dim), flux.Rows(diffop->Dim()*IntRange(k,k+1)));
    else
      // diffop->Apply(fel, mir, x.Slice(comp, dim), flux.RowSlice(comp,dim));
      diffop->Apply(fel, mir, x.Slice(comp, dim), flux.Rows(diffop->Dim()*IntRange(comp,comp+1)));
  }
  
  
  void BlockDifferentialOperatorTrans ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              BareSliceVector<double> bx, 
              LocalHeap & lh) const
  {
    // HeapReset hr(lh);
    auto x = bx.Range(0,fel.GetNDof()*dim);
    // FlatVector<> hx(fel.GetNDof(), lh);
    // FlatVector<> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            // hflux = flux.Slice(k, dim);
            diffop->ApplyTrans(fel, mip, flux.Range(diffop->Dim()*IntRange(k,k+1)), x.Slice(k,dim), lh);
            // x.Slice(k,dim) = hx;
          }
      }
    else
      {
        // hflux = flux.Slice(comp, dim);
        x = 0.0;
        diffop->ApplyTrans(fel, mip, flux.Range(diffop->Dim()*IntRange(comp,comp+1)), x.Slice(comp,dim), lh);
        // x.Slice(comp,dim) = hx;
      }
  }

  void BlockDifferentialOperatorTrans ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<Complex> flux,
              BareSliceVector<Complex> bx, 
              LocalHeap & lh) const
  {
    // HeapReset hr(lh);
    auto x = bx.Range(0,dim*fel.GetNDof());
    // FlatVector<Complex> hx(fel.GetNDof(), lh);
    // FlatVector<Complex> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            // hflux = flux.Slice(k, dim);
            diffop->ApplyTrans(fel, mip, flux.Range(diffop->Dim()*IntRange(k,k+1)), x.Slice(k,dim), lh);
            // x.Slice(k,dim) = hx;
          }
      }
    else
      {
        // hflux = flux.Slice(comp, dim);
        x = 0.0;
        diffop->ApplyTrans(fel, mip, flux.Range(diffop->Dim()*IntRange(comp,comp+1)), x.Slice(comp,dim), lh);
        // x.Slice(comp,dim) = hx;
      }
  }

    
  void BlockDifferentialOperatorTrans ::
  AddTrans (const FiniteElement & fel,
            const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<double>> flux,
            BareSliceVector<double> x) const
  {
    if (comp == -1)
      for (size_t k = 0; k < dim; k++)
        diffop->AddTrans(fel, mir, flux.Rows(diffop->Dim()*IntRange(k,k+1)), x.Slice(k,dim));
    else
      diffop->AddTrans(fel, mir, flux.Rows(diffop->Dim()*IntRange(comp,comp+1)), x.Slice(comp,dim));
  }

  void BlockDifferentialOperatorTrans ::
  AddTrans (const FiniteElement & fel,
            const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<Complex>> flux,
            BareSliceVector<Complex> x) const
  {
    if (comp == -1)
      for (size_t k = 0; k < dim; k++)
        diffop->AddTrans(fel, mir, flux.Rows(diffop->Dim()*IntRange(k,k+1)), x.Slice(k,dim));
    else
      diffop->AddTrans(fel, mir, flux.Rows(diffop->Dim()*IntRange(comp,comp+1)), x.Slice(comp,dim));
  }


  shared_ptr<CoefficientFunction> BlockDifferentialOperatorTrans ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir) const 
  {
    return TransposeCF (diffop->DiffShape(TransposeCF(proxy), dir));
  }









  VectorDifferentialOperator :: ~VectorDifferentialOperator ()  { ; }


  void VectorDifferentialOperator ::
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & mip,
              SliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    auto & fel = static_cast<const CompoundFiniteElement&> (bfel)[0];

    size_t ndi = fel.GetNDof();
    size_t dimi = diffop->Dim();

    mat = 0.0;
    diffop->CalcMatrix (fel, mip, mat.Rows(dimi).Cols(ndi), lh);
    for (int i = 1; i < dim; i++)
      mat.Rows(i*dimi, (i+1)*dimi).Cols(i*ndi, (i+1)*ndi) = mat.Rows(dimi).Cols(ndi);
  }
  
  void VectorDifferentialOperator ::
  CalcMatrix (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & mir,
              BareSliceMatrix<SIMD<double>> bmat) const
  {
    auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
    auto & feli = fel[0]; 

    size_t ndi = feli.GetNDof();
    size_t dimi = diffop->Dim();

    auto mat = bmat.AddSize(dimi*dim*bfel.GetNDof(), mir.Size());
    mat = 0.0;

    diffop->CalcMatrix (feli, mir, mat.Rows(dimi*ndi));
    for (int i = 1; i < dim; i++)
      {
        auto mati = mat.Rows(dim*dimi*fel.GetRange(i));
        for (int j = 0; j < feli.GetNDof(); j++)
          mati.Rows(j*dim*dimi+i*dimi, j*dim*dimi+(i+1)*dimi)
            = mat.Rows(j*dimi, (j+1)*dimi);
      }
    for (int j = feli.GetNDof()-1; j >= 0; j--)
      mat.Rows(j*dim*dimi, j*dim*dimi+dimi) = mat.Rows(j*dimi, (j+1)*dimi);
    for (int j = feli.GetNDof()-1; j >= 0; j--)
      mat.Rows(j*dim*dimi+dimi, (j+1)*dim*dimi) = 0.0;

    //    throw ExceptionNOSIMD("VectorDifferentialOperator::CalcMatrix not yet tested for SIMD support");
  }
  

  void VectorDifferentialOperator ::
  Apply (const FiniteElement & bfel,
         const BaseMappedIntegrationPoint & mip,
         BareSliceVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    auto & fel = static_cast<const CompoundFiniteElement&> (bfel)[0];
    size_t ndi = fel.GetNDof();
    size_t dimi = diffop->Dim();    

    for (int k = 0; k < dim; k++)
      diffop->Apply(fel, mip, x.Range(ndi*k, ndi*(k+1)), flux.Range(k*dimi, (k+1)*dimi), lh);
  }


  void VectorDifferentialOperator ::
  Apply (const FiniteElement & bfel,
         const SIMD_BaseMappedIntegrationRule & mir,
         BareSliceVector<double> x, 
         BareSliceMatrix<SIMD<double>> flux) const
  {
    auto & fel = static_cast<const CompoundFiniteElement&> (bfel)[0];
    size_t ndi = fel.GetNDof();
    size_t dimi = diffop->Dim();    

    for (int k = 0; k < dim; k++)
      diffop->Apply(fel, mir, x.Range(k*ndi, (k+1)*ndi), flux.Rows(k*dimi, (k+1)*dimi));
  }
  
  
  void VectorDifferentialOperator ::
  ApplyTrans (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              BareSliceVector<double> x, 
              LocalHeap & lh) const
  {
    auto & fel = static_cast<const CompoundFiniteElement&> (bfel)[0];
    size_t ndi = fel.GetNDof();
    size_t dimi = diffop->Dim();    

    for (int k = 0; k < dim; k++)
      diffop->ApplyTrans(fel, mip, flux.Range(k*dimi, (k+1)*dimi), x.Range(k*ndi, (k+1)*ndi), lh);
  }

  void VectorDifferentialOperator ::
  ApplyTrans (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<Complex> flux,
              BareSliceVector<Complex> x, 
              LocalHeap & lh) const
  {
    auto & fel = static_cast<const CompoundFiniteElement&> (bfel)[0];
    size_t ndi = fel.GetNDof();
    size_t dimi = diffop->Dim();    

    for (int k = 0; k < dim; k++)
      diffop->ApplyTrans(fel, mip, flux.Range(k*dimi, (k+1)*dimi), x.Range(k*ndi, (k+1)*ndi), lh);
  }

    
  void VectorDifferentialOperator ::
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<double>> flux,
            BareSliceVector<double> x) const
  {
    auto & fel = static_cast<const CompoundFiniteElement&> (bfel)[0];
    size_t ndi = fel.GetNDof();
    size_t dimi = diffop->Dim();
    
    for (size_t k = 0; k < dim; k++)
      diffop->AddTrans(fel, mir, flux.Rows(k*dimi, (k+1)*dimi), x.Range(k*ndi, (k+1)*ndi));
  }

  void VectorDifferentialOperator ::
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<Complex>> flux,
            BareSliceVector<Complex> x) const
  {
    auto & fel = static_cast<const CompoundFiniteElement&> (bfel)[0];
    size_t ndi = fel.GetNDof();
    size_t dimi = diffop->Dim();
    
    for (size_t k = 0; k < dim; k++)
      diffop->AddTrans(fel, mir, flux.Rows(k*dimi, (k+1)*dimi), x.Range(k*ndi, (k+1)*ndi));
  }


  shared_ptr<CoefficientFunction> VectorDifferentialOperator ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir) const 
  {
    int ddim = diffop->Dim();
    Array<shared_ptr<CoefficientFunction>> proxys(dim);
    
    for (int i = 0; i < dim; i++)
      {
        Array<shared_ptr<CoefficientFunction>> tmp(ddim);
        for(int j = 0; j < diffop->Dim(); j++)
            tmp[j] = MakeComponentCoefficientFunction(proxy,i*dim+j);
        proxys[i] = MakeVectorialCoefficientFunction(move(tmp));
      }
    
    Array<shared_ptr<CoefficientFunction>> cflist(dim);
    for (int i = 0; i < dim; i++)
        cflist[i] = diffop->DiffShape(proxys[i], dir);
    auto result = MakeVectorialCoefficientFunction(move(cflist));
    result->SetDimensions( Array({dim,diffop->Dim()}) );

    return result;
  }





  





  
  /*
  template class T_DifferentialOperator<DiffOpId<1> >;
  template class T_DifferentialOperator<DiffOpId<2> >;
  template class T_DifferentialOperator<DiffOpId<3> >;

  template class T_DifferentialOperator<DiffOpIdBoundary<1> >;
  template class T_DifferentialOperator<DiffOpIdBoundary<2> >;
  template class T_DifferentialOperator<DiffOpIdBoundary<3> >;


  template class T_DifferentialOperator<DiffOpGradient<1> >;
  template class T_DifferentialOperator<DiffOpGradient<2> >;
  template class T_DifferentialOperator<DiffOpGradient<3> >;

  template class T_DifferentialOperator<DiffOpIdEdge<2> >;
  template class T_DifferentialOperator<DiffOpIdEdge<3> >;

  template class T_DifferentialOperator<DiffOpCurlEdge<2> >;
  template class T_DifferentialOperator<DiffOpCurlEdge<3> >;
  */


}
