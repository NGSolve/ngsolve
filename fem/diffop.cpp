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
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
#ifndef FASTCOMPILE
    cout << "called base class apply, type = " << typeid(*this).name() << endl;
#endif
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    flux = mat * x;
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<Complex> x, 
         FlatVector<Complex> flux,
         LocalHeap & lh) const
  {
#ifndef FASTCOMPILE
    cout << "called base class apply, complex" << endl;
#endif
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    flux = mat * x;
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationRule & mir,
         FlatVector<double> x, 
         FlatMatrix<double> flux,
         LocalHeap & lh) const
  {
    for (int i = 0; i < mir.Size(); i++)
      Apply (fel, mir[i], x, flux.Row(i), lh);
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationRule & mir,
         FlatVector<Complex> x, 
         FlatMatrix<Complex> flux,
         LocalHeap & lh) const
  {
    for (int i = 0; i < mir.Size(); i++)
      Apply (fel, mir[i], x, flux.Row(i), lh);
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
              FlatVector<double> x, 
              LocalHeap & lh) const 
  {
#ifndef FASTCOMPILE
    cout << "called base class apply trans, type = " << typeid(*this).name() << endl;
#endif
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    x = Trans(mat) * flux;
  }
  
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<Complex> flux,
              FlatVector<Complex> x, 
              LocalHeap & lh) const 
  {
#ifndef FASTCOMPILE
    cout << "called base class apply trans complex, type = " << typeid(*this).name() << endl;
#endif
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    x = Trans(mat) * flux;
  }
  
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              FlatMatrix<double> flux,
              FlatVector<double> x, 
              LocalHeap & lh) const 
  {
    HeapReset hr(lh);
    FlatVector<double> hx(x.Size(), lh);
    x = 0.0;
    for (int i = 0; i < mir.Size(); i++)
      {
        ApplyTrans (fel, mir[i], flux.Row(i), hx, lh);
        x += hx;
      }
  }
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              FlatMatrix<Complex> flux,
              FlatVector<Complex> x, 
              LocalHeap & lh) const 
  {
    HeapReset hr(lh);
    FlatVector<Complex> hx(x.Size(), lh);
    x = 0.0;
    for (int i = 0; i < mir.Size(); i++)
      {
        ApplyTrans (fel, mir[i], flux.Row(i), hx, lh);
        x += hx;
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
              
              col.AddSize(dim2_dim_do) = 0;
            
              for (size_t l = 0; l < dim_diffop; l++)
                col.Slice(l*hdim, dim_dim_do+1).AddSize(hdim) = hval[l];
            }
        }
    else
      throw ExceptionNOSIMD("BlockDifferentialOperator::CalcMatrix does not support SIMD");
  }
  

  void BlockDifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<> hx(fel.GetNDof(), lh);
    FlatVector<> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            hx = x.Slice(k, dim);
            diffop->Apply(fel, mip, hx, hflux, lh);
            flux.Slice(k,dim) = hflux;
          }
      }
    else
      {
        hx = x.Slice(comp, dim);
        diffop->Apply(fel, mip, hx, hflux, lh);
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
              FlatVector<double> x, 
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<> hx(fel.GetNDof(), lh);
    FlatVector<> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            hflux = flux.Slice(k, dim);
            diffop->ApplyTrans(fel, mip, hflux, hx, lh);
            x.Slice(k,dim) = hx;
          }
      }
    else
      {
        hflux = flux.Slice(comp, dim);
        diffop->ApplyTrans(fel, mip, hflux, hx,lh);
        x = 0.0;
        x.Slice(comp,dim) = hx;
      }
  }

  void BlockDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<Complex> flux,
              FlatVector<Complex> x, 
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<Complex> hx(fel.GetNDof(), lh);
    FlatVector<Complex> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            hflux = flux.Slice(k, dim);
            diffop->ApplyTrans(fel, mip, hflux, hx, lh);
            x.Slice(k,dim) = hx;
          }
      }
    else
      {
        hflux = flux.Slice(comp, dim);
        diffop->ApplyTrans(fel, mip, hflux, hx,lh);
        x = 0.0;
        x.Slice(comp,dim) = hx;
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
