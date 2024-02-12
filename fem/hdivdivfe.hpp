#ifndef FILE_HDIVDIVFE
#define FILE_HDIVDIVFE

/*********************************************************************/
/* File:   hdivdivfe.hpp                                             */
/* Author: Astrid Pechstein, Joachim Schoeberl                       */
/* Date:   orig 2006, redesign Dec 2016                              */
/*********************************************************************/

#include "finiteelement.hpp"
#include "fe_interfaces.hpp"

#include "hcurlfe.hpp" // for Cross (AD,AD)
#include "recursive_pol.hpp"
#include "recursive_pol_trig.hpp"
#include "recursive_pol_tet.hpp"
#include "shapefunction_utils.hpp"

namespace ngfem
{

  class BaseHDivDivFiniteElement : public FiniteElement
  {
  public:

    using FiniteElement::ndof;
    using FiniteElement::order;
    bool algebraic_mapping = true;
    
    INLINE BaseHDivDivFiniteElement () { ; } 
    INLINE BaseHDivDivFiniteElement (int andof, int aorder)
      : FiniteElement (andof, aorder) { ; }

    void SetAlgebraicMapping (bool am) { algebraic_mapping = am; }
    
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceMatrix<double> shape) const = 0;

    virtual void CalcDivShape (const IntegrationPoint & ip,
			       BareSliceMatrix<double> divshape) const = 0;
  };
  template <int DIM>
  class HDivDivFiniteElement : public BaseHDivDivFiniteElement
  {
  public:
    using BaseHDivDivFiniteElement::BaseHDivDivFiniteElement;
    using BaseHDivDivFiniteElement::ndof;
    using BaseHDivDivFiniteElement::order;

    // old style
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const = 0;

    virtual void CalcDivShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<double> divshape) const = 0;

    // new implementation
    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM> & mip,
      BareSliceMatrix<double> shape) const = 0;

    /*
    virtual void CalcDDMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM> & mip,
                                           BareSliceMatrix<double> shape) const = 0;
    */
    
    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM> & mip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedDivShape (const BaseMappedIntegrationPoint & mip,
      BareSliceMatrix<double> shape) const = 0;


    virtual void CalcMappedShape_Matrix (const SIMD_BaseMappedIntegrationRule & mir, 
                                         BareSliceMatrix<SIMD<double>> shapes) const = 0;
    
    virtual void Evaluate_Matrix (const SIMD_BaseMappedIntegrationRule & ir,
                                  BareSliceVector<> coefs,
                                  BareSliceMatrix<SIMD<double>> values) const = 0;

    virtual void AddTrans_Matrix (const SIMD_BaseMappedIntegrationRule & ir,
                                  BareSliceMatrix<SIMD<double>> values,
                                  BareSliceVector<> coefs) const = 0;

    virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const = 0;
    virtual void CalcDualShape (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> shape) const = 0;
    virtual void EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const = 0;
    virtual void AddDualTrans (const SIMD_BaseMappedIntegrationRule& bmir, BareSliceMatrix<SIMD<double>> values, BareSliceVector<double> coefs) const = 0;


    virtual void CalcMappedDivShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                      BareSliceMatrix<SIMD<double>> divshapes) const = 0;

    virtual void EvaluateDiv (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs,
			      BareSliceMatrix<SIMD<double>> values) const = 0;

    virtual void AddDivTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
			      BareSliceVector<> coefs) const=0;

    virtual void CalcShape_NormalComponent (const SIMD_BaseMappedIntegrationRule & mir, 
                                            BareSliceMatrix<SIMD<double>> shapes) const = 0;
    
    virtual list<tuple<string,double>> Timing () const;
  };

  template <int D,typename VEC,typename MAT>
  void VecToSymMat(const VEC & vec,MAT && mat)
  {
    switch(D)
    {
    case 2:
      mat(0) = vec(0);
      mat(3) = vec(1);
      mat(1) = mat(2) = vec(2);
      break;
    case 3:
      auto v0 = vec(0);
      auto v1 = vec(1);
      auto v2 = vec(2);
      auto v3 = vec(3);
      auto v4 = vec(4);
      auto v5 = vec(5);
      mat(0) = v0;
      mat(1) = v5;
      mat(2) = v4;
      mat(3) = v5;
      mat(4) = v1;
      mat(5) = v3;
      mat(6) = v4;
      mat(7) = v3;
      mat(8) = v2;
      /*
      mat(0) = vec(0);
      mat(4) = vec(1);
      mat(8) = vec(2);
      mat(1) = mat(3) = vec(5);
      mat(2) = mat(6) = vec(4);
      mat(5) = mat(7) = vec(3);
      */
      break;
    }
  }
  
  template <typename T>
  auto SymMatToVecDual (const Mat<2,2,T> & mat) 
  {
    return Vec<3,T> { mat(0,0), mat(1,1), mat(0,1)+mat(1,0) };
  }
  
  template <typename T>
  auto SymMatToVecDual (const Mat<3,3,T> & mat) 
  {
    return Vec<6,T> { mat(0,0), mat(1,1), mat(2,2), mat(1,2)+mat(2,1), mat(0,2)+mat(2,0), mat(0,1)+mat(1,0) };
  }

  
  template <typename T>
  Mat<2,2,T> DyadProd(Vec<2,T> a, Vec<2,T> b)
  {
    // return Matrix<T>({{a(0)*b(0), a(0)*b(1)}, {a(1)*b(0), a(1)*b(1)}} );
    return { a(0)*b(0), a(0)*b(1),  a(1)*b(0), a(1)*b(1) };
  }

  template <typename T>
  Mat<3,3,T> DyadProd(Vec<3,T> a, Vec<3,T> b)
  {
    // return Matrix<T>( {{a(0)*b(0), a(0)*b(1), a(0)*b(2)}, {a(1)*b(0), a(1)*b(1), a(1)*b(2)}, {a(2)*b(0), a(2)*b(1), a(2)*b(2)}} );
    return { a(0)*b(0), a(0)*b(1), a(0)*b(2),   a(1)*b(0), a(1)*b(1), a(1)*b(2),   a(2)*b(0), a(2)*b(1), a(2)*b(2) };
  }

  template <ELEMENT_TYPE ET> class HDivDivFE;

  
  template <ELEMENT_TYPE ET, typename SHAPES = HDivDivFE<ET>>
  class T_HDivDivFE : public HDivDivFiniteElement<ET_trait<ET>::DIM>,
    public VertexOrientedFE<ET>
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
    enum { DIM_STRESS = (DIM*(DIM+1))/2 };
    
    using VertexOrientedFE<ET>::vnums;
    using HDivDivFiniteElement<ET_trait<ET>::DIM>::ndof;
    using HDivDivFiniteElement<ET_trait<ET>::DIM>::order;

    //enum { N_VERTEX = ET_trait<ET>::N_VERTEX };
    //enum { N_FACET   = ET_trait<ET>::N_FACET };    
    //
    //size_t vnums[N_VERTEX];
    IVec<DIM-1> order_facet[ET_trait<ET>::N_FACET];
    IVec<DIM> order_inner;

    // additional div-div free bubbles
    bool plus;

  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    T_HDivDivFE (int aorder, bool _plus = false)
      : plus(_plus)
    {
      order = aorder;
      for (auto & of : order_facet) of = aorder;
      order_inner = aorder;
      //ndof = DIM*(DIM+1)/2 * ET_trait<ET>::PolDimension(aorder);

    }
    
    virtual ELEMENT_TYPE ElementType() const override { return ET; }
    // const HDivDivFE<ET> * Cast() const { return static_cast<const HDivDivFE<ET>*> (this); }
    auto * Cast() const { return static_cast<const SHAPES*> (this); } 
    
    INLINE void SetOrderFacet (int nr, IVec<DIM-1,int> order) { order_facet[nr] = order; }
    INLINE void SetOrderInner (IVec<DIM,int> order) { order_inner = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HDivDivFE<ET>:: ComputeNDof not available, only for ET == TRIG" << endl;
    }

   template <typename T, typename TFA> 
   void T_CalcShape (TIP<DIM,AutoDiff<DIM,T>> tip, TFA & shape) const
    {
      if constexpr (DIM == 2)
                     Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM,T>> (tip), shape);
      else
        Cast() -> T_CalcShape (tip, shape);        
    }

   template <typename T, typename TFA> 
   void T_CalcShape (TIP<DIM,AutoDiffDiff<DIM,T>> tip, TFA & shape) const
    {
      if constexpr (DIM == 2)
                     Cast() -> T_CalcShape (tip, shape);
      else
        throw Exception ("dd shapes are not supported in 3D");
    }

    
    // old style
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const override
    {
      // Vec<DIM, AutoDiff<DIM> > adp = ip;
      /*
      Vec<DIM, AutoDiff<DIM>> adp;
      for ( int i=0; i<DIM; i++)
        adp(i) = AutoDiff<DIM>(ip(i),i);
      auto tip = TIP<DIM, AutoDiff<DIM>> (adp, ip.FacetNr(), ip.VB());
      */
      /* Cast() -> */ T_CalcShape (GetTIPGrad<DIM>(ip), 
                                   SBLambda([&] (int nr, auto val)
                                            {
                                              shape.Row(nr).Range(0,DIM_STRESS) = val.Shape();
                                            }));
    }

    virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const override
    {
      shape.AddSize(ndof, sqr(bmip.DimSpace())) = 0.0;
      Switch<4-DIM>
        (bmip.DimSpace()-DIM,[this, &bmip, shape](auto CODIM)
         {
           auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);

           Cast() -> CalcDualShape2 (mip, SBLambda([&] (size_t nr, auto val)
                                                   {
                                                     shape.Row(nr) = val.AsVector();
                                                   }));
         });
      
    }

    virtual void CalcDualShape (const SIMD_BaseMappedIntegrationRule& bmir, BareSliceMatrix<SIMD<double>> shapes) const override
    {
      Switch<4-DIM>
        (bmir.DimSpace()-DIM,[this, &bmir, shapes](auto CODIM)
         {
           constexpr int DIMSPACE = DIM+CODIM.value;
           auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM+CODIM.value>&> (bmir);

           shapes.AddSize(ndof*sqr(DIMSPACE), mir.Size()) = 0.0;
           for (size_t i = 0; i < mir.Size(); i++)
             {
               Cast() -> CalcDualShape2 (mir[i], SBLambda([shapes,i,DIMSPACE] (size_t j, auto val)
                                                          {
                                                            shapes.Rows(j*sqr(DIMSPACE), (j+1)*sqr(DIMSPACE)).Col(i).Range(0,sqr(DIMSPACE)) = val.AsVector();
                                                          }));
             }
         });
    }
    
    virtual void EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const override
    {
      Switch<4-DIM>
        (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
                             {
                               constexpr int DIMSPACE = DIM+CODIM.value;
                               auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM+CODIM.value>&> (bmir);
                               for (size_t i = 0; i < mir.Size(); i++)
                                 {
                                   Mat<DIMSPACE,DIMSPACE,SIMD<double>> sum (SIMD<double>(0.0));
                                   Cast() -> CalcDualShape2 (mir[i], SBLambda([&sum, coefs] (size_t j, auto val)
                                                                              {
                                                                                sum += coefs(j) * val;
                                                                              }));
                                   for (size_t k = 0; k < sqr(DIMSPACE); k++)
                                     values(k, i) = sum(k);
                                 }});
    }
    
    virtual void AddDualTrans (const SIMD_BaseMappedIntegrationRule& bmir, BareSliceMatrix<SIMD<double>> values, BareSliceVector<double> coefs) const override
    {
      Switch<4-DIM>
        (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
                             {
                               constexpr int DIMSPACE = DIM+CODIM.value;
                               auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM+CODIM.value>&> (bmir);
                               for (size_t i = 0; i < mir.Size(); i++)
                                 {
                                   Mat<DIMSPACE,DIMSPACE,SIMD<double>> value;
                                   for (size_t k = 0; k < sqr(DIMSPACE); k++)
                                     value(k) = values(k, i);
                                   
                                   Cast()-> CalcDualShape2 (mir[i], SBLambda([value, coefs] (size_t j, auto val)
                                                                             {
                                                                               coefs(j) += HSum(InnerProduct(val,value));
                                                                             }));
                                 }});
    }


    virtual void CalcDivShape (const IntegrationPoint & ip,
                               BareSliceMatrix<double> shape) const override
    {
      // MSVC internal compiler error
      // Vec<DIM, AutoDiff<DIM> > adp = ip;
      // TIP<DIM,AutoDiffDiff<DIM>> addp(adp);      
      Vec<DIM, AutoDiff<DIM>> adp;
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM>(ip(i),i);

      /* Cast() -> */ T_CalcShape (TIP<DIM, AutoDiff<DIM>> (adp, ip.FacetNr(), ip.VB()),
                             SBLambda([shape] (int nr, auto val)
                                      {
                                        shape.Row(nr).Range(0,DIM) = val.DivShape();
                                      }));
    }

    // new style
    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM> & mip,
                            BareSliceMatrix<double> shape) const override
    {
      Vec<DIM, AutoDiff<DIM>> adp = mip;
      /*
      Vec<DIM, AutoDiffDiff<DIM>> addp;
      for (int i=0; i<DIM; i++)
      {
        addp[i] = adp[i].Value();
        addp[i].LoadGradient(&adp[i].DValue(0));
      }
      */
      /* Cast() -> */ T_CalcShape (TIP<DIM, AutoDiff<DIM>> (adp, mip.IP().FacetNr(), mip.IP().VB()),
                             SBLambda([&] (int nr, auto val)
                                      {
                                        shape.Row(nr).Range(0,DIM_STRESS) = val.Shape();
                                      }));
    }


    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM> & mip,
                                         BareSliceMatrix<double> shape) const override
    {
      /*
      auto tip = this->algebraic_mapping ? TIP<DIM, AutoDiffDiff<DIM>>(GetTIP(mip)) : GetTIPHesse(mip);
      T_CalcShape (tip, 
                   SBLambda([&](int nr,auto val)
                            {
                              VecToSymMat<DIM> (val.Shape(), shape.Row(nr));
                            }));
      */
      if (this->algebraic_mapping)
        T_CalcShape (GetTIP(mip), 
                     SBLambda([&](int nr,auto val)
                              {
                                VecToSymMat<DIM> (val.Shape(), shape.Row(nr));
                              }));
      else
        T_CalcShape (GetTIPHesse(mip), 
                     SBLambda([&](int nr,auto val)
                              {
                                VecToSymMat<DIM> (val.Shape(), shape.Row(nr));
                              }));
    }

    /*
    virtual void CalcDDMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM> & mip,
                                           BareSliceMatrix<double> shape) const override
    {
      T_CalcShape (GetTIPHesse(mip), 
                   SBLambda([&](int nr,auto val)
                            {
                              VecToSymMat<DIM> (val.Shape(), shape.Row(nr));
                            }));
    }
    */
    

    virtual void CalcMappedDivShape (const BaseMappedIntegrationPoint & bmip,
                                     BareSliceMatrix<double> shape) const override
    {
      auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM>&> (bmip);
      if (!this->algebraic_mapping)
        {
          T_CalcShape (GetTIPHesse(mip), 
                       SBLambda([&](int nr,auto val)
                                {
                                  shape.Row(nr).Range(0,DIM) = val.DivShape();
                                }));
          return;
        }

      if (!mip.GetTransformation().IsCurvedElement()) // non-curved element
        {
          T_CalcShape (GetTIP(mip), 
                       SBLambda([&](int nr,auto val)
                                {
                                  shape.Row(nr).Range(0,DIM) = val.DivShape();
                                }));
          return;
        }
      // curved element
      
      if (false) // eval on physical element
        {
          Mat<DIM> inv = mip.GetJacobianInverse();

          Vec<DIM,Mat<DIM>> hesse, hesse_inv;
          mip.CalcHesse (hesse);

          /*
            
            div ( 1/J F sigma_ref  F^T 1/J ) 
            = deriv (1/J F)  sigma_ref F^T 1/J  +  1/J F div(sigma_ref F^T 1/J) = I + II

            I  ... Hessian : sigma + grad(1/J) sigma 
            II ... by DivShape of SigmaGrad templates, as for non-curved elements

           */ 

          /*
          for (int k = 0; k < DIM; k++)
            hesse_inv(k) = Trans(inv) * hesse(k) * inv;

          Vec<DIM> gradJ_xi;
          for (int k = 0; k < DIM; k++)
            {
              double sum = 0;
              for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++)
                  sum += inv(j,i) * hesse(i)(j,k);
              gradJ_xi(k) = sum;
            }
          Vec<DIM> gradJ = Trans(inv) * gradJ_xi;
          for (int k = 0; k < DIM; k++)
            for (int l = 0; l < DIM; l++)
              hesse_inv(k)(l,k) -= gradJ(l);
          */

          // saving a view operations ...
          Vec<DIM, Mat<DIM>> hesse_inv1;
          
          for (int k = 0; k < DIM; k++)
            hesse_inv1(k) = hesse(k) * inv;          
          
          Vec<DIM> gradJ_xi = 0.0;
          for (int i = 0; i < DIM; i++)
            gradJ_xi += hesse_inv1(i).Col(i);
          
          for (int k = 0; k < DIM; k++)
            hesse_inv1(k).Col(k) -= gradJ_xi;

          for (int k = 0; k < DIM; k++)
            hesse_inv(k) = Trans(inv) * hesse_inv1(k);

          Mat<DIM,DIM_STRESS> hesse_inv_vec;
          for (int k = 0; k < DIM; k++)
            hesse_inv_vec.Row(k) = SymMatToVecDual(hesse_inv(k));

          T_CalcShape (GetTIP(mip),  
                       SBLambda([&] (int nr,auto val)
                                {
                                  shape.Row(nr).Range(0,DIM) = val.DivShape() + hesse_inv_vec * val.Shape();
                                }));
          return;
        }


      if (true) 
        {
          // eval on reference element
          
          Mat<DIM> inv = mip.GetJacobianInverse();
          Mat<DIM> jac = mip.GetJacobian();
          double det = Det(jac);
          
          Vec<DIM,Mat<DIM>> hesse = mip.CalcHesse();
          
          Vec<DIM> gradJ_xi = 0.0;
          for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
              gradJ_xi += inv(j,i) * hesse(i).Col(j);

          for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
              hesse(i).Col(j) -= jac(i,j) * gradJ_xi;
          
          Mat<DIM,DIM_STRESS> hesse_inv_vec;
          for (int k = 0; k < DIM; k++)
            hesse_inv_vec.Row(k) = 1/(det*det) * SymMatToVecDual(hesse(k));            
          Mat<DIM> trans_div = 1/(det*det) * jac;
          
          T_CalcShape (GetTIPGrad<DIM>(mip.IP()),  
                       SBLambda([&] (int nr,auto val) LAMBDA_INLINE
                                {
                                  shape.Row(nr).Range(0,DIM) = trans_div * val.DivShape() + hesse_inv_vec * val.Shape();
                                }));
          return;
        }

      

      
        /*  
        Mat<DIM> jac = mip.GetJacobian();
        Mat<DIM> inv_jac = mip.GetJacobianInverse();
        Mat<DIM> hesse[3],finvT_h_tilde_finv[3];
        mip.CalcHesse (hesse[0],hesse[1],hesse[2]);

        Mat<DIM,DIM,AutoDiff<DIM> > f_tilde;
        for(int i = 0; i < DIM; i++)
        {
          for(int j = 0; j < DIM; j++)
          {
            f_tilde(i,j).Value() = jac(i,j);
            for(int k = 0; k < DIM; k++)
              f_tilde(i,j).DValue(k) = hesse[i](j,k);
          }
        }

        AutoDiff<DIM> ad_det = Det (f_tilde);
        AutoDiff<DIM> iad_det = 1.0 / ad_det;
        f_tilde *= iad_det;

        for(int i=0; i<DIM; i++)
        {
          finvT_h_tilde_finv[i] = 0;
          for(int alpha=0; alpha<DIM; alpha++)
            for(int beta=0; beta<DIM; beta++)
              for(int gamma=0; gamma<DIM; gamma++)
                for(int delta=0; delta<DIM; delta++)
                  finvT_h_tilde_finv[i](alpha,beta) += inv_jac(gamma,alpha)*f_tilde(i,gamma).DValue(delta)*inv_jac(delta,beta);
        }

        T_CalcShape (GetTIP(mip),  // TIP<DIM,AutoDiff<DIM>> (adp),
                               SBLambda([&](int nr,auto val)
                                  {
                                    shape.Row(nr).Range(0,DIM) = val.DivShape();
                                    BareVector<double> divshape = shape.Row(nr);
                                    Vec<DIM_STRESS> vecshape = val.Shape();
                                    Vec<DIM*DIM> matshape;
                                    VecToSymMat<DIM> (vecshape, matshape);

                                    for(int k=0; k<DIM; k++)
                                    {
                                      for(int j=0; j<DIM*DIM; j++)
                                      {
                                        divshape(k) += mip.GetJacobiDet() * finvT_h_tilde_finv[k](j) * matshape(j);
                                      }
                                    }
                                    
                                  }));
      }
        */
      
    }

    template <int DIMSPACE>
    void CalcMappedShape_Matrix2 (const SIMD_MappedIntegrationRule<DIM,DIMSPACE> & mir, 
                                 BareSliceMatrix<SIMD<double>> shapes) const
    {
      // static Timer t("HDivDivFE - Matrix2", NoTracing);
      // RegionTracer regtr(TaskManager::GetThreadId(), t);    

      for (size_t i = 0; i < mir.Size(); i++)
        {
          if (DIM == DIMSPACE)
            {
              const SIMD_BaseMappedIntegrationRule & bmir = mir;
              auto & mir2 = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&>(bmir);
              // Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir2[i];
              auto shapesi = shapes.Col(i);
              /* this->Cast() -> */ T_CalcShape (GetTIP(mir2[i]), // TIP<DIM,AutoDiff<DIM,SIMD<double>>>(adp),
                                           SBLambda ([shapesi] (size_t j, auto val) 
                                                     {
                                                       auto shapeij = shapesi.Range(j*sqr(DIMSPACE),(j+1)*sqr(DIMSPACE));
                                                       // VecToSymMat<DIM> (val.Shape(), shapeij);
                                                       Mat<DIM,DIM,SIMD<double>> shapemat;
                                                       VecToSymMat<DIM> (val.Shape(), shapemat);
                                                       for (size_t i = 0; i < DIMSPACE*DIMSPACE; i++)
                                                         shapeij(i) = shapemat(i);
                                                     }));
            }
          else
            {
              auto jac = mir[i].GetJacobian();
              auto d2 = sqr(mir[i].GetJacobiDet());
              Vec<DIM_STRESS,SIMD<double>> hv;
              Mat<DIM,DIM,SIMD<double>> mat;
              SIMD<double> mem[DIMSPACE*DIMSPACE*DIM_STRESS];
              // FlatMatrix<SIMD<double>> trans(DIMSPACE*DIMSPACE,DIM_STRESS,&mem[0]);
              FlatMatrixFixWidth<DIM_STRESS, SIMD<double>> trans(DIMSPACE*DIMSPACE,&mem[0]);
              for (int k = 0; k < DIM_STRESS; k++)
                {
                  hv = SIMD<double>(0.0);
                  hv(k) = SIMD<double>(1.0);
                  VecToSymMat<DIM> (hv, mat);
                  Mat<DIMSPACE,DIMSPACE,SIMD<double>> physmat = 1/d2*(jac * mat * Trans(jac));
                  trans.Col(k) = physmat.AsVector();
                }
              
              
              // Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir.IR()[i];
              // TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
              
              /* this->Cast() -> */ T_CalcShape (GetTIPGrad<DIM> (mir.IR()[i]), // TIP<DIM,AutoDiff<DIM,SIMD<double>>>(adp),
                                           SBLambda ([i,shapes,trans] (size_t j, auto val) 
                                                     {
                                                       shapes.Rows(j*sqr(DIMSPACE),(j+1)*sqr(DIMSPACE)).Col(i).Range(0,sqr(DIMSPACE)) = trans * val.Shape();
                                                     }));
            }
        }
    }

    virtual void CalcShape_NormalComponent (const SIMD_BaseMappedIntegrationRule & bmir, 
                                            BareSliceMatrix<SIMD<double>> shapes) const override
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);

      for (size_t i = 0; i < mir.Size(); i++)
        {
          // Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir[i];
          Vec<DIM,SIMD<double>> nv = mir[i].GetNV();
          auto shapesi = shapes.Col(i);
          /* this->Cast() -> */ T_CalcShape (GetTIP(mir[i]),   // TIP<DIM,AutoDiff<DIM,SIMD<double>>>(adp),
                                       SBLambda ([shapesi, nv] (size_t j, auto val) 
                                                 {
                                                   auto shapeij = shapesi.Range(j*DIM,(j+1)*DIM);
                                                   Mat<DIM,DIM,SIMD<double>> shapemat;
                                                   VecToSymMat<DIM> (val.Shape(), shapemat);
                                                   Vec<DIM,SIMD<double>> mnv = shapemat * nv;
                                                   for (size_t i = 0; i < DIM; i++)
                                                     shapeij(i) = mnv(i);
                                                 }));
        }
    }
    
    virtual void CalcMappedShape_Matrix (const SIMD_BaseMappedIntegrationRule & bmir, 
                                         BareSliceMatrix<SIMD<double>> shapes) const override
    {
      Iterate<4-DIM>
        ([this, &bmir, shapes](auto CODIM) LAMBDA_INLINE
         {
           constexpr int CD = CODIM.value;
           constexpr int DIMSPACE = DIM+CD;
           if (bmir.DimSpace() == DIMSPACE)
             {
               auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
               this->CalcMappedShape_Matrix2 (mir, shapes);
#ifdef XXX
               
               for (size_t i = 0; i < mir.Size(); i++)
                 {
                   auto jac = mir[i].GetJacobian();
                   auto d2 = sqr(mir[i].GetJacobiDet());

                   Vec<DIM_STRESS,SIMD<double>> hv;
                   Mat<DIM,DIM,SIMD<double>> mat;
                   // Mat<DIMSPACE*DIMSPACE, DIM_STRESS,SIMD<double>> trans;
                   SIMD<double> mem[DIMSPACE*DIMSPACE*DIM_STRESS];
                   FlatMatrix<SIMD<double>> trans(DIMSPACE*DIMSPACE,DIM_STRESS,&mem[0]);
                   for (int k = 0; k < DIM_STRESS; k++)
                     {
                       hv = SIMD<double>(0.0);
                       hv(k) = SIMD<double>(1.0);
                       VecToSymMat<DIM> (hv, mat);
                       Mat<DIMSPACE,DIMSPACE,SIMD<double>> physmat =
                         1/d2 * (jac * mat * Trans(jac));
                       for (int j = 0; j < DIMSPACE*DIMSPACE; j++)
                         trans(j,k) = physmat(j);
                     }
                   
          
                   Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = bmir.IR()[i];
                   // TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);

                   this->Cast() -> T_CalcShape (adp,
                                                SBLambda ([i,shapes,trans] (size_t j, auto val) LAMBDA_INLINE
                                                    {
                                                      /*
                                                      Mat<DIM,DIM,SIMD<double>> mat;
                                                      VecToSymMat<DIM> (val.Shape(), mat);
                                                      Mat<DIMSPACE,DIMSPACE,SIMD<double>> physmat =
                                                        1/d2 * (jac * mat * Trans(jac));
                                                      for (size_t k = 0; k < sqr(DIMSPACE); k++)
                                                        shapes(j*sqr(DIMSPACE)+k,i) = physmat(k);
                                                      */
                                                      Vec<DIMSPACE*DIMSPACE,SIMD<double>> transvec;
                                                      transvec = (trans * val.Shape()).AsVector();
                                                      for (size_t k = 0; k < sqr(DIMSPACE); k++)
                                                        shapes(j*sqr(DIMSPACE)+k,i) = transvec(k);
                                                    }));
                 }
#endif
             }
         });
    }

    virtual void Evaluate_Matrix (const SIMD_BaseMappedIntegrationRule & bmir,
                                  BareSliceVector<> coefs,
                                  BareSliceMatrix<SIMD<double>> values) const override
    {
      if (this->algebraic_mapping == false)
        {
          if (bmir.DimSpace() != DIM)
            throw Exception ("sequential mapping only for volume space");
          
          auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
          for (size_t i = 0; i < bmir.Size(); i++)
            {
              double *pcoefs = &coefs(0);
              const size_t dist = coefs.Dist();
              
              Vec<DIM_STRESS,SIMD<double>> sum(0.0);
              T_CalcShape (GetTIPHesse(mir[i]), 
                           SBLambda ([&sum,&pcoefs,dist] (size_t j, auto val)
                                     {
                                       sum += (*pcoefs)*val.Shape();
                                       pcoefs += dist;
                                     }));
              for (size_t k = 0; k < DIM_STRESS; k++)
                values(k,i) = sum(k);
            }
          return;
        }
      
      for (size_t i = 0; i < bmir.Size(); i++)
        {
          double *pcoefs = &coefs(0);
          const size_t dist = coefs.Dist();
          
          Vec<DIM_STRESS,SIMD<double>> sum(0.0);
          // Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = bmir.IR()[i];
          // TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          
          /* Cast() -> */ T_CalcShape (GetTIPGrad<DIM>(bmir.IR()[i]), // TIP<DIM,AutoDiff<DIM,SIMD<double>>>(adp),
                                 SBLambda ([&sum,&pcoefs,dist] (size_t j, auto val)
                                           {
                                             sum += (*pcoefs)*val.Shape();
                                             pcoefs += dist;
                                           }));

          Mat<DIM,DIM,SIMD<double>> summat;
          VecToSymMat<DIM> (sum, summat);
          
          Iterate<4-DIM>
            ([values,&bmir,i,summat](auto CODIM)
             {
               constexpr auto DIMSPACE = DIM+CODIM.value;
               if (bmir.DimSpace() == DIMSPACE)
                 {
                   auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
                   auto jac = mir[i].GetJacobian();
                   auto d2 = sqr(mir[i].GetJacobiDet());
                   Mat<DIMSPACE,DIMSPACE,SIMD<double>> physmat = 1/d2 * (jac * summat * Trans(jac));
                   for (size_t k = 0; k < sqr(DIMSPACE); k++)
                     values(k,i) = physmat(k);
                 }
             });
        }
    }

    virtual void AddTrans_Matrix (const SIMD_BaseMappedIntegrationRule & bmir,
                                  BareSliceMatrix<SIMD<double>> values,
                                  BareSliceVector<> coefs) const override
    {
       for (size_t i = 0; i < bmir.Size(); i++)
        {
          Mat<DIM,DIM,SIMD<double>> mat;
          
          Iterate<4-DIM>
            ([&bmir,i,&mat,values](auto CODIM)
             {
               constexpr auto DIMSPACE = DIM+CODIM.value;
               if (bmir.DimSpace() == DIMSPACE)
                 {
                   auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);

                   auto jac = mir[i].GetJacobian();
                   auto d2 = sqr(mir[i].GetJacobiDet());

                   Mat<DIMSPACE,DIMSPACE,SIMD<double>> physmat{};
                   // physmat = values.Col(i);
                   physmat.AsVector() = values.Col(i);                   
                   mat = 1/d2 * Trans(jac) * physmat * jac;
                 }
             });
          
          // Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = bmir.IR()[i];
          // TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          double *pcoefs = &coefs(0);
          const size_t dist = coefs.Dist();

          /* Cast() -> */ T_CalcShape (GetTIPGrad<DIM>(bmir.IR()[i]), // TIP<DIM,AutoDiff<DIM,SIMD<double>>>  (adp),
                                 SBLambda ([mat,&pcoefs,dist] (size_t j, auto val)
                                           {
                                             Mat<DIM,DIM,SIMD<double>> mat2;
                                             VecToSymMat<DIM> (val.Shape(), mat2);
                                             
                                             *pcoefs += HSum(InnerProduct(mat,mat2));
                                             pcoefs += dist;
                                           }));
        }
    }


    virtual void CalcMappedDivShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                      BareSliceMatrix<SIMD<double>> divshapes) const override
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);


      if(!mir.GetTransformation().IsCurvedElement()) // non-curved element
      {
	for (size_t i = 0; i < mir.Size(); i++)
	{
          auto jac = mir[i].GetJacobian();
          auto d2 = sqr(mir[i].GetJacobiDet());
          
          Vec<DIM,SIMD<double>> vec;
          SIMD<double> mem[DIM*DIM_STRESS];
          FlatMatrix<SIMD<double>> trans(DIM,DIM,&mem[0]);
          trans = 1/d2 * jac;
            /*
            {
              vec = SIMD<double>(0.0);
              vec(k) = SIMD<double>(1.0);
              Vec<DIM,SIMD<double>> physvec = 1/d2 * (jac * vec);
              trans.Col(k) = physvec;
            }
            */
          // Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir.IR()[i];
          // TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
	  /* Cast() -> */ T_CalcShape
            (GetTIPGrad<DIM>(mir.IR()[i]), // TIP<DIM,AutoDiff<DIM,SIMD<double>>>(adp),
             SBLambda([divshapes,i,trans](int j,auto val)
                      {
                        divshapes.Rows(j*DIM,(j+1)*DIM).Col(i).Range(0,DIM) = trans * val.DivShape();
                      }));
	}
      }
      else
      {
        // throw ExceptionNOSIMD(string("HDivDiv - CalcMappedDivShape SIMD only for noncurved elements"));
        for (size_t i = 0; i < mir.Size(); i++)
          {
            // static Timer t0("HDivDivFE - hesse", NoTracing);
            // static Timer t1("HDivDivFE - prepare div", NoTracing);
            // static Timer t2("HDivDivFE - calc div", NoTracing);
            
            Mat<DIM,DIM,SIMD<double>> jac = mir[i].GetJacobian();
            Mat<DIM,DIM,SIMD<double>> inv_jac = mir[i].GetJacobianInverse();
            Mat<DIM,DIM,SIMD<double>> finvT_h_tilde_finv[DIM];

            // RegionTracer reg0(TaskManager::GetThreadId(), t0);
            
            Vec<DIM, Mat<DIM,DIM,SIMD<double>>> hesse;
            mir.GetTransformation().CalcHesse (mir.IR()[i], hesse);

            // RegionTracer reg1(TaskManager::GetThreadId(), t1);    

            Mat<DIM,DIM,AutoDiff<DIM,SIMD<double>> > f_tilde;
            for(int i = 0; i < DIM; i++)
              for(int j = 0; j < DIM; j++)
                {
                  f_tilde(i,j).Value() = jac(i,j);
                  for(int k = 0; k < DIM; k++)
                    f_tilde(i,j).DValue(k) = hesse[i](j,k);
                }
            
            AutoDiff<DIM,SIMD<double>> ad_det = Det (f_tilde);
            AutoDiff<DIM, SIMD<double>> iad_det = 1.0 / ad_det;
            f_tilde *= iad_det;
            
            for(int i=0; i<DIM; i++)
              {
                finvT_h_tilde_finv[i] = 0;
                for(int alpha=0; alpha<DIM; alpha++)
                  for(int beta=0; beta<DIM; beta++)
                    for(int gamma=0; gamma<DIM; gamma++)
                      for(int delta=0; delta<DIM; delta++)
                        finvT_h_tilde_finv[i](alpha,beta) += inv_jac(gamma,alpha)*f_tilde(i,gamma).DValue(delta)*inv_jac(delta,beta);
              }
            for (int j = 0; j < DIM; j++)
              finvT_h_tilde_finv[j] *= mir[i].GetJacobiDet();
            
            // Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            // TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
            // TIP<DIM,AutoDiff<DIM,SIMD<double>>> addp(adp);
            // RegionTracer reg2(TaskManager::GetThreadId(), t2);    
            /* Cast() -> */ T_CalcShape
              (GetTIP(mir[i]), // TIP<DIM,AutoDiff<DIM,SIMD<double>>>(adp),
               SBLambda([&](int nr,auto val)
                              {
                                BareSliceVector<SIMD<double>> divshape = divshapes.Rows(nr*DIM,(nr+1)*DIM).Col(i);

                                Vec<DIM,SIMD<double>> div1 = val.DivShape();
                                Vec<DIM_STRESS,SIMD<double>> vecshape = val.Shape();
                                Vec<DIM*DIM,SIMD<double>> matshape;
                                VecToSymMat<DIM> (vecshape, matshape);
                                
                                for(size_t k = 0; k < DIM; k++)
                                  {
                                    SIMD<double> sum = div1(k);
                                    for(size_t j = 0; j < DIM*DIM; j++)
                                      sum += finvT_h_tilde_finv[k](j) * matshape(j);
                                    divshape(k) = sum;
                                  }
                              }));
          }
      }
    }

    virtual void EvaluateDiv (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs,
			      BareSliceMatrix<SIMD<double>> values) const override
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
      if(!mir.GetTransformation().IsCurvedElement()) // non-curved element
      {
      for (size_t i = 0; i < bmir.Size(); i++)
        {
          double *pcoefs = &coefs(0);
          const size_t dist = coefs.Dist();
          
          Vec<DIM,SIMD<double>> sum(0.0);
          // Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = bmir.IR()[i];
          // TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          
          /* Cast() -> */ T_CalcShape (GetTIPGrad<DIM>(bmir.IR()[i]), // TIP<DIM,AutoDiff<DIM,SIMD<double>>>(adp),
                                 SBLambda ([&sum,&pcoefs,dist] (size_t j, auto val)
                                           {
                                             sum += (*pcoefs)*val.DivShape();
                                             pcoefs += dist;
                                           }));

          Iterate<4-DIM>
            ([values,&bmir,i,sum](auto CODIM)
             {
               constexpr auto DIMSPACE = DIM+CODIM.value;
               if (bmir.DimSpace() == DIMSPACE)
                 {
                   auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
                   auto jac = mir[i].GetJacobian();
                   auto d2 = sqr(mir[i].GetJacobiDet());
                   Vec<DIMSPACE,SIMD<double>> physvec = 1/d2 * (jac * sum);
                   for (size_t k=0; k < DIMSPACE; k++)
                     values(k,i) = physvec(k);
                 }
             });
        }
      }
      else
      {
        throw ExceptionNOSIMD(string("HDivDiv - EvaluateDiv SIMD only for noncurved elements"));
      }
    }

    virtual void AddDivTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
			      BareSliceVector<> coefs) const override
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
      if(!mir.GetTransformation().IsCurvedElement()) // non-curved element
      {
      for (size_t i = 0; i < bmir.Size(); i++)
        {
          Vec<DIM,SIMD<double>> vec;
          
          Iterate<4-DIM>
            ([&bmir,i,&vec,values](auto CODIM)
             {
               constexpr auto DIMSPACE = DIM+CODIM.value;
               if (bmir.DimSpace() == DIMSPACE)
                 {
                   auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);

                   auto jac = mir[i].GetJacobian();
                   auto d2 = sqr(mir[i].GetJacobiDet());

                   Vec<DIMSPACE,SIMD<double>> physvec{};
                   for (size_t k = 0; k < DIMSPACE; k++)
                     physvec(k) = values(k,i);
                   vec = 1/d2 * Trans(jac) * physvec;
                 }
             });
          
          // Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = bmir.IR()[i];
          // TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          double *pcoefs = &coefs(0);
          const size_t dist = coefs.Dist();

          /* Cast() -> */ T_CalcShape (GetTIPGrad<DIM> (bmir.IR()[i]), // TIP<DIM,AutoDiff<DIM,SIMD<double>>>(adp),
                                 SBLambda ([vec,&pcoefs,dist] (size_t j, auto val)
                                           {
                                             *pcoefs += HSum(InnerProduct(vec,val.DivShape()));
                                             pcoefs += dist;
                                           }));
        }
      }
      else
      {
        throw ExceptionNOSIMD(string("HDivDiv - AddTrans SIMD only for noncurved elements"));

      }
    }

  };



#ifdef FILE_HDIVDIVFE_CPP
#define HDIVDIVFE_EXTERN
#else
#define HDIVDIVFE_EXTERN extern
#endif

  extern template class HDivDivFiniteElement<2>;
  extern template class HDivDivFiniteElement<3>;
  
  
  // ***************** SigmaGrad ****************************** */
  // sigma (nabla u)
  
  template <int D, typename T> class T_SigmaGrad;
  template <typename T> class T_SigmaGrad<2,T>
  {
    AutoDiffDiff<2,T> u;
  public:
    T_SigmaGrad  (AutoDiffDiff<2,T> au) : u(au) { ; }
    Vec<3,T> Shape() { return Vec<3,T> (u.DDValue(1,1), u.DDValue(0,0), -u.DDValue(1,0)); }
    Vec<2,T> DivShape() { return Vec<2,T> (0.0, 0.0); }
  };
  
  template <int D, typename T>
  auto SigmaGrad (AutoDiffDiff<D,T> au) { return T_SigmaGrad<D,T>(au); }
  
  
  // ***************** Sigma_u_Gradv ****************************** */
  // sigma (u nabla v)
  
  template <int D, typename T> class T_Sigma_u_Gradv;
  template <typename T> class T_Sigma_u_Gradv<2,T>
  {
    AutoDiffDiff<2,T> u, v;
  public:
    T_Sigma_u_Gradv  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av) { ; }
    Vec<3,T> Shape() { return Vec<3,T> ((u.Value()*v.DDValue(1,1) + u.DValue(1)*v.DValue(1)),
                                        (u.Value()*v.DDValue(0,0) + u.DValue(0)*v.DValue(0)),
                                        -u.Value()*v.DDValue(1,0) - 0.5 * (u.DValue(0)*v.DValue(1)+u.DValue(1)*v.DValue(0))); }
    Vec<2,T> DivShape()
    {
      T uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      T vx = v.DValue(0), vy = v.DValue(1);
      
      return -0.5 * Vec<2,T> (uyy*vx - uxy*vy + uy*vxy - ux*vyy,
                            -uxy*vx + uxx*vy - uy*vxx + ux*vxy);
    }
  };
  
  template <int D, typename T>
  auto Sigma_u_Gradv (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_Sigma_u_Gradv<D,T>(au, av); }
  
  // ***************** Type2 ****************************** */
  // ????
  
  template <int D, typename T> class T_Type2;
  template <typename T> class T_Type2<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_Type2  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av) { ; }
    Vec<3,T> Shape() { return Vec<3,T> (u.DDValue(1,1)*v.Value() - 2*u.DValue(1)*v.DValue(1) + u.Value()*v.DDValue(1,1),
                                        u.DDValue(0,0)*v.Value() - 2*u.DValue(0)*v.DValue(0) + u.Value()*v.DDValue(0,0),
                                        -(u.DDValue(0,1)*v.Value() - u.DValue(0)*v.DValue(1) -
                                          u.DValue(1)*v.DValue(0) + u.Value()*v.DDValue(1,0))); }

    Vec<2,T> DivShape()
    {
      T uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      T vx = v.DValue(0), vy = v.DValue(1);
      
      return Vec<2,T> (2*uyy*vx + 2*ux*vyy - 2*uxy*vy - 2*uy*vxy, 2*uxx*vy + 2*uy*vxx - 2*uxy*vx - 2*ux*vxy);
    }
  };
  
  template <int D, typename T>
  auto Type2 (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_Type2<D,T>(au, av); }
  
  // ***************** Type3 ****************************** */
  // ????
  
  template <int D, typename T> class T_Type3;
  template <typename T> class T_Type3<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_Type3  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av) { ; }
    Vec<3,T> Shape() { return Vec<3,T> (u.DDValue(1,1)*v.Value() - u.Value()*v.DDValue(1,1),
                                    u.DDValue(0,0)*v.Value() - u.Value()*v.DDValue(0,0),
                                    -(u.DDValue(0,1)*v.Value() - u.Value()*v.DDValue(1,0))); }
    Vec<2,T> DivShape()
    {
      T uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      T vx = v.DValue(0), vy = v.DValue(1);

      return Vec<2,T> (uyy*vx - uxy*vy - ux*vyy + uy*vxy, uxx*vy - uxy*vx - uy*vxx + ux*vxy);
    }
  };
  
  template <int D, typename T>
  auto Type3 (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_Type3<D,T>(au, av); }

 
  template <int D, typename T> class T_vSigmaGradu;
  template <typename T> class T_vSigmaGradu<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_vSigmaGradu  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av) { ; }
    Vec<3,T> Shape() { return Vec<3,T> (u.DDValue(1,1)*v.Value(),
                                      u.DDValue(0,0)*v.Value(),  -(u.DDValue(0,1)*v.Value()));}
    Vec<2,T> DivShape()
    {
      T uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      // T ux = u.DValue(0), uy = u.DValue(1);
      // T vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      T vx = v.DValue(0), vy = v.DValue(1);

      return Vec<2,T> (uyy*vx- uxy*vy, uxx*vy- uxy*vx);
    }
  };
  
  template <int D, typename T>
  auto vSigmaGradu (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_vSigmaGradu<D,T>(au, av); }

  // ***************** Sigma ((vDu - uDv) w) ****************************** */
  // where u, v are NOW POSSIBLY NON-linear hat basis functions (i.e. vDu - uDv is Nedelec0 edge basis function)
  template <int D, typename T> class T_Sigma_Duv_minus_uDv_w;
  template <typename T> class T_Sigma_Duv_minus_uDv_w<2,T>
  {
    AutoDiffDiff<2,T> u,v,w;
  public:
    T_Sigma_Duv_minus_uDv_w  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiffDiff<2,T> aw) : u(au), v(av), w(aw) { ; }
    Vec<3,T> Shape() { return Vec<3,T> (w.DValue(1)*(v.DValue(1)*u.Value()-u.DValue(1)*v.Value()), 
      w.DValue(0)*(v.DValue(0)*u.Value()-u.DValue(0)*v.Value()),
      -0.5*( w.DValue(0)*(v.DValue(1)*u.Value()-u.DValue(1)*v.Value()) +
        w.DValue(1)*(v.DValue(0)*u.Value()-u.DValue(0)*v.Value()) )
      ); }

    Vec<2,T> DivShape()
    {
      T uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      T vx = v.DValue(0), vy = v.DValue(1);
      T wxx = w.DDValue(0,0), wyy = w.DDValue(1,1), wxy = w.DDValue(0,1);
      T wx = w.DValue(0), wy = w.DValue(1);

      return Vec<2,T> (0.5*wxy*(vy*u.Value() - uy*v.Value()) 
        -0.5*wyy*(vx*u.Value() - ux*v.Value())
        +1.5*wy*(vy*ux - uy*vx) + 0.5*wy*(vxy*u.Value()-uxy*v.Value()) 
        -0.5*wx*(vyy*u.Value()-uyy*v.Value()),
        0.5*wxy*(vx*u.Value() - ux*v.Value()) 
        -0.5*wxx*(vy*u.Value() - uy*v.Value())
        +1.5*wx*(vx*uy - ux*vy) + 0.5*wx*(vxy*u.Value()-uxy*v.Value()) 
        -0.5*wy*(vxx*u.Value()-uxx*v.Value())
        ); 
    }

  };
  

  

  class T_SymRotRot_Dl2xDl1_v
  {
    AutoDiff<2> l1,l2,v;
  public:
    T_SymRotRot_Dl2xDl1_v  (AutoDiff<2> lam1, AutoDiff<2> lam2, AutoDiff<2> av) : l1(lam1), l2(lam2), v(av) { ; }
    Vec<3> Shape() { return Vec<3> (v.Value()*(l1.DValue(1)*l2.DValue(1)),
      v.Value()*(l1.DValue(0)*l2.DValue(0)),
      -0.5*v.Value()*(l1.DValue(1)*l2.DValue(0) + l1.DValue(0)*l2.DValue(1))
      ); }

    Vec<2> DivShape()
    {
      // todo
      // double lam1 = l1.Value();
      double lam1x = l1.DValue(0);
      double lam1y = l1.DValue(1);
      // double lam2 = l2.Value();
      double lam2x = l2.DValue(0);
      double lam2y = l2.DValue(1);
      return Vec<2> (
        v.DValue(0)*(lam1y*lam2y) - 0.5*v.DValue(1)*(lam1x*lam2y+lam1y*lam2x),
        -0.5*v.DValue(0)*(lam1x*lam2y+lam1y*lam2x) + v.DValue(1)*(lam1x*lam2x)
        ); 
    }

  };

  template <int D, typename T>
  auto Sigma_Duv_minus_uDv_w (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av, AutoDiffDiff<D,T> aw)
  { return T_Sigma_Duv_minus_uDv_w<D,T>(au, av, aw); }
  
  
  template <ELEMENT_TYPE ET> class HDivDivFE : public T_HDivDivFE<ET> 
  {
  protected:
    using T_HDivDivFE<ET> :: order;
    using T_HDivDivFE<ET> :: ndof;
  public:
    template <typename T, typename TFA> 
    void T_CalcShape (TIP<ET_trait<ET>::DIM,AutoDiffDiff<ET_trait<ET>::DIM,T>> ip, TFA & shape) const
    {
      throw Exception ("Hdivdivfe not implementend for element type");
    }

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hdivdivfe not implementend for element type");
    }
  };

  
  template <> class HDivDivFE<ET_TRIG> : public T_HDivDivFE<ET_TRIG> 
  {
    
  public:
    using T_HDivDivFE<ET_TRIG> :: T_HDivDivFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<3; i++)
      {
        ndof += order_facet[i][0]+1;
        order = max2(order, order_facet[i][0]);
      }
      int ninner = 3*order_inner[0]*(order_inner[0]+1)/2 ;
      order = max2(order, order_inner[0]);
      if (plus)
      { 
        order ++;
        ninner += 2*order_inner[0]; 
      }
      ndof += ninner;

    }
   template <typename T, typename TFA> 
   void T_CalcShape (TIP<2,AutoDiffDiff<2,T>> ip, TFA & shape) const
    {
      // typedef decltype(ip.x.Value()+ip.x.Value()) T;
      typedef AutoDiffDiff<2, T> Tx;
      Tx x = ip.x, y = ip.y;
      Tx ddlami[3] ={ x, y, 1-x-y };
      
      int ii = 0;
      
      int maxorder_facet =
        max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));

      const EDGE * edges = ElementTopology::GetEdges(ET_TRIG);

      ArrayMem<Tx,20> ha(maxorder_facet+1);
      ArrayMem<Tx,20> u(order_inner[0]+2), v(order_inner[0]+2);
      
      for (int i = 0; i < 3; i++)
        {
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es,ee);
          
          Tx ls = ddlami[es], le = ddlami[ee];
          
          // edge functions are all div-free!
          IntegratedLegendreMonomialExt::CalcTrigExt(maxorder_facet+2,
                                                     le-ls, 1-le-ls, ha);

          //ScaledLegendrePolynomial(maxorder_facet,le-ls, 1-le-ls,ha);

          
          for (int l = 0; l <= order_facet[i][0]; l++)
            shape[ii++] = SigmaGrad (ha[l]);
            //shape[ii++] = SymRotRot_Dl2xDl1_v_diffdiff(le, ls, ha[l]);
        }
      
      int es = 0; int ee = 1; int et = 2;
      Tx ls = ddlami[es];
      Tx le = ddlami[ee];
      Tx lt = ddlami[et];
      
      int oi=order_inner[0];
      int oi_plus = oi; //plus ? oi+1 : oi;


      IntegratedLegendreMonomialExt::CalcTrigExt(oi_plus+3,le-ls,1-le-ls,u);
      LegendrePolynomial::EvalMult(oi_plus+1, 2*lt-1, lt, v);
      
      
      for(int i = 0; i <= oi-1; i++)
      {
        for(int j = 0; j+i <= oi-1; j++)
        {
          shape[ii++] = SigmaGrad(u[i]*v[j]);
          shape[ii++] = Type2(u[i],v[j]);
        }
      }
      for(int i = 0; i <= oi_plus-1; i++)
      {
        for(int j = 0; j+i <= oi_plus-1; j++)
        {
          if(j > 0)
            shape[ii++] = Type3(u[i],v[j]);
        }
      }
      
      for (int i = 0; i < oi_plus; i++)
        shape[ii++] = Sigma_Duv_minus_uDv_w (le, -ls, v[i]);
      
      //// element bubbles for Sigma+ space
      if (plus)
        for (int i = 0; i <= oi-1; i++)
          {
            Tx bubble = u[i]*v[oi-1-i];
            shape[ii++] = Sigma_u_Gradv(bubble, x);
            shape[ii++] = Sigma_u_Gradv(bubble, y);
          }
    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      auto & ip = mip.IP();
      typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
      T x = ip(0), y = ip(1);
      T lam[3] = { x, y, 1-x-y };
      Vec<2,T> pnts[3] = { { 1, 0 }, { 0, 1 } , { 0, 0 } };
      int facetnr = ip.FacetNr();

      int ii = 0;

      // const EDGE * edges = ElementTopology::GetEdges(ET_TRIG);

      if (ip.VB() == BND)
        { // facet shapes
          for (int i = 0; i < 3; i++)
            {
              int p = order_facet[i][0];
              
              if (i == facetnr)
                {
                  IVec<2> e = ET_trait<ET_TRIG>::GetEdgeSort (i, vnums);
                  
                  T xi = lam[e[0]]-lam[e[1]];
                  Vec<2,T> tauref = pnts[e[0]] - pnts[e[1]];
                  
                  Vec<2,T> nvref = Vec<2,T>(tauref[1],-tauref[0]);
                  auto nv = Trans(mip.GetJacobianInverse())*nvref;
                  auto nn = DyadProd(nv,nv);
                  
                  LegendrePolynomial::Eval
                    (p, xi,
                     SBLambda([&] (size_t nr, T val)
                              {
                                shape[nr+ii] = mip.GetMeasure()*val*nn;
                              }));
                }
              ii += (p+1);
            }
        }
      else
        {
          for (int i = 0; i < 3; i++)
            ii += order_facet[i][0]+1;
        }
      if (ip.VB() == VOL)
        {
          auto p = order_inner[0]-1;
          if( p >= 0 )
            {
              DubinerBasis::Eval (p, lam[0], lam[1],
                                  SBLambda([&] (size_t nr, T val)
                                           {
                                             shape[ii++] = val*mip.GetMeasure()*mip.GetJacobian()*Mat<2,2>({{1,0},{0,0}})*Trans(mip.GetJacobian());
                                             shape[ii++] = val*mip.GetMeasure()*mip.GetJacobian()*Mat<2,2>({{0,0},{0,1}})*Trans(mip.GetJacobian());
                                             shape[ii++] = val*mip.GetMeasure()*mip.GetJacobian()*Mat<2,2>({{0,1},{1,0}})*Trans(mip.GetJacobian());
                                           }));
            }
        }
    }
    
  };
  
  template <> class HDivDivFE<ET_QUAD> : public T_HDivDivFE<ET_QUAD> 
  {
    
  public:
    using T_HDivDivFE<ET_QUAD> :: T_HDivDivFE;

    enum {incsg = -1};
    enum {incsugv = -1};

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<4; i++)
      {
        ndof += order_facet[i][0]+1;
        order = max2(order, order_facet[i][0]);
      }
      int ninner = (order_inner[0]+1+incsg)*(order_inner[0]+1+incsg) + 
        (order_inner[0]+2)*(order_inner[0]) *2 +
        2*(order_inner[0]+1+incsugv) +1;
      
      /*
      int ninner = (order_inner[0]+1+incsg)*(order_inner[0]+1+incsg) + 
        (order_inner[0]+1)*(order_inner[0]) *2 +
        2*(order_inner[0]+1+incsugv) +2;
      if (plus) ninner += order_inner[0]*2;
      */
      order = max2(order, order_inner[0]);
      order += 1;
      ndof += ninner;

    }
   template <typename T, typename TFA> 
   void T_CalcShape (TIP<2,AutoDiffDiff<2,T>> ip, TFA & shape) const
    {
      // typedef decltype(ip.x.Value()+ip.x.Value()) T;
      typedef AutoDiffDiff<2, T> Tx;
      
      Tx x = ip.x, y = ip.y;
      Tx lx[4] ={1-x, x, x, 1-x};
      Tx ly[4] = {1-y, 1-y, y, y};
      
      int ii = 0;

      const EDGE * edges = ElementTopology::GetEdges(ET_QUAD);

      ArrayMem<Tx,20> u(order+2), v(order+2);
      
      for (int i = 0; i < 4; i++)
        {
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es,ee);
          
          Tx xi = lx[ee]+ly[ee]-lx[es]-ly[es];
          Tx eta = lx[es]*ly[es]+lx[ee]*ly[ee];

	  IntegratedLegendreMonomialExt::Calc(order_facet[i][0]+2,xi,u);

          
          for (int l = 0; l <= order_facet[i][0]; l++)
            // shape[ii++] = SigmaGrad (eta*u[l]);
            shape[ii++] = Sigma_u_Gradv (eta, u[l]);
            // shape[ii++] = vSigmaGradu(u[i],eta);
        }


      int oi=order_inner[0];

      IntegratedLegendreMonomialExt::Calc(oi+3,lx[0]-lx[1],u);
      IntegratedLegendreMonomialExt::Calc(oi+3,ly[0]-ly[2],v);

      // original 
      for(int i = 0; i <= oi+incsg; i++)
      {
        for(int j = 0; j <= oi+incsg; j++)
        {
          shape[ii++] = SigmaGrad(u[i]*v[j]);
        }
      }
      for(int i = 0; i <= oi+1; i++)
      {
        for(int j = 0; j <= oi-1; j++)
        {
          shape[ii++] = vSigmaGradu(u[i],v[j]);
          shape[ii++] = vSigmaGradu(v[i],u[j]);
        }
      }

      shape[ii++] = Sigma_u_Gradv(lx[0], ly[0]);

      for(int i = 0; i <= oi+incsugv; i++)
      {
        shape[ii++] = Sigma_u_Gradv(u[i], ly[0]);
        shape[ii++] = Sigma_u_Gradv(v[i], lx[0]); //
      }
      
      return;


      
      for (int i = 0; i <= oi+incsg+1; i++)
        for (int j = 0; j <= oi+incsg+1; j++)
          shape[ii++] = SigmaGrad(u[i]*v[j]);
      
      /*
      for (int j = 0; j <= oi-1; j++)
        {
          shape[ii++] = SigmaGrad(u[oi+incsg+1]*v[j]);
          shape[ii++] = SigmaGrad(v[oi+incsg+1]*u[j]);
        }
      */

      for (int i = 0; i <= oi-1; i++)
        for (int j = 0; j <= oi-1; j++)
          {
            shape[ii++] = vSigmaGradu(u[i],v[j]);
            shape[ii++] = vSigmaGradu(v[i],u[j]);
          }

      if (plus)
        for (int j = 0; j <= oi-1; j++)
          {
            // shape[ii++] = vSigmaGradu(u[oi+1],v[j]);
            // shape[ii++] = vSigmaGradu(v[oi+1],u[j]);

            shape[ii++] = Sigma_u_Gradv(v[j], u[oi]);
            shape[ii++] = Sigma_u_Gradv(u[j], v[oi]);            
          }

      // shape[ii++] = Sigma_u_Gradv(lx[0], ly[0]);
      shape[ii++] = SigmaGrad((2*x-1)*(2*y-1));
      
      for (int i = 0; i <= oi+incsugv; i++)
        {
          shape[ii++] = Sigma_u_Gradv(u[i], ly[0]);
          shape[ii++] = Sigma_u_Gradv(v[i], lx[0]); //
        }
    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hdivdivfe not implementend for element type");
    }

  };



  class HDivDivFE_QuadFullPol : public T_HDivDivFE<ET_QUAD, HDivDivFE_QuadFullPol> 
  {
    
  public:
    using T_HDivDivFE<ET_QUAD,HDivDivFE_QuadFullPol> :: T_HDivDivFE;
    
    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<4; i++)
        {
          ndof += order_facet[i][0]+1;
          order = max2(order, order_facet[i][0]);
        }
      
      int ninner = 1  //  sigma grad(xy)
        + 2 * (order_inner[0]+1)*(order_inner[0]+2)  // inner nedelec
        + sqr(order_inner[0]+1)
        ;
      if (plus)
        ninner += 4*order_inner[0] + 4;
      order = max2(order, order_inner[0]);
      order += 2;
      if (plus) order++;
      ndof += ninner;
    }
    
    template <typename T, typename TFA> 
    void T_CalcShape (TIP<2,AutoDiffDiff<2,T>> ip, TFA & shape) const
    {
      typedef AutoDiffDiff<2, T> Tx;
      
      Tx x = ip.x, y = ip.y;
      Tx lx[4] ={1-x, x, x, 1-x};
      Tx ly[4] = {1-y, 1-y, y, y};
      
      int ii = 0;

      const EDGE * edges = ElementTopology::GetEdges(ET_QUAD);

      ArrayMem<Tx,20> u(order+2), v(order+2);
      
      for (int i = 0; i < 4; i++)
        {
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es,ee);
          
          Tx xi = lx[ee]+ly[ee]-lx[es]-ly[es];
          Tx eta = lx[es]*ly[es]+lx[ee]*ly[ee];

	  IntegratedLegendreMonomialExt::Calc(order_facet[i][0]+2,xi,u);
          
          for (int l = 0; l <= order_facet[i][0]; l++)
            // shape[ii++] = SigmaGrad (eta*u[l]);
            shape[ii++] = Sigma_u_Gradv (eta, u[l]);
        }

      // shape[ii++] = SigmaGrad((2*x-1)*(2*y-1));
      shape[ii++] = Sigma_u_Gradv((2*x-1),(2*y-1));

      int oi=order_inner[0];

      LegendrePolynomial (oi+1, 2*x-1, u);
      LegendrePolynomial (oi+1, 2*y-1, v);
      auto bubx = x*(1-x);
      auto buby = y*(1-y);
      for (int i = 0; i <= oi; i++)
        for (int j = 0; j <= oi+1; j++)
          {
            shape[ii++] = Sigma_u_Gradv(bubx*u[i]*v[j], 2*y-1);
            shape[ii++] = Sigma_u_Gradv(buby*u[j]*v[i], 2*x-1);
          }

      if (plus)
        {
          for (int i = 0; i <= oi; i++)
            {
              shape[ii++] = Sigma_u_Gradv(bubx*buby*u[i]*v[oi], 2*y-1);            
              shape[ii++] = Sigma_u_Gradv(bubx*buby*u[oi]*v[i], 2*x-1);
              shape[ii++] = Sigma_u_Gradv(bubx*buby*u[i]*v[oi+1], 2*y-1);            
              shape[ii++] = Sigma_u_Gradv(bubx*buby*u[oi+1]*v[i], 2*x-1);
            }
        }
      
      for (int i = 0; i <= oi; i++)
        for (int j = 0; j <= oi; j++)
          shape[ii++] = vSigmaGradu(bubx,u[i]*v[j]*buby);

      
      if (ii != ndof)
        cerr << "Hdivdivfe, full quad, ndof = " << ndof << " != ii = " << ii << endl;
    }


    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hdivdivfe not implementend for quadfullpol");
    }
    
  };










  
  // ***************** S_zz(uvw) ****************************** */
  // write uvw into zz component
  template <int D> class T_S_zz;
  template <> class T_S_zz<3>
  {
    AutoDiff<2> u, v;
    AutoDiff<1> w;
  public:
    T_S_zz ( AutoDiff<2> au, AutoDiff<2> av, AutoDiff<1> aw) : u(au), v(av), w(aw) { ; }
    Vec<6> Shape() 
    { 
      Vec<6> sigma(0.);
      sigma[2] = u.Value()*v.Value()*w.Value();
      return sigma;
    }

    Vec<3> DivShape()
    {
      return Vec<3> (0., 0., u.Value()*v.Value()*w.DValue(0));
    }

  };
  
  template <int D>
  auto S_zz (AutoDiff<D> au, AutoDiff<D> av, AutoDiff<1> aw)
  { return T_S_zz<D+1>(au, av, aw); }
    
  // ***************** S_xz ****************************** */
  template <int D, typename T> class T_S_xz;
  template <typename T> class T_S_xz<3,T>
  {
    AutoDiff<2,T> uv;
    AutoDiff<1,T> w;

    int comp;
  public:
    T_S_xz ( int acomp, AutoDiff<2,T> auv, AutoDiff<1,T> aw) : uv(auv), w(aw), comp(acomp) { ; }
    Vec<6,T> Shape() 
    { 
      Vec<6,T> sigma;
      sigma = 0.;
      if (comp==0)
        sigma[4] = uv.Value()*w.Value();
      else
        sigma[3] = uv.Value()*w.Value();
      return sigma;
    }

    Vec<3,T> DivShape()
    {
      if (comp == 0)
        return Vec<3,T> (uv.Value()*w.DValue(0), 0, uv.DValue(0)*w.Value() );
      else
        return Vec<3,T> (0, uv.Value()*w.DValue(0), uv.DValue(1)*w.Value() );
    }

  };
  
  template <int D, typename T>
  auto S_xz (int comp, AutoDiff<D,T> auv, AutoDiff<1,T> aw)
  { return T_S_xz<D+1,T>(comp,auv, aw); }



  template <typename T>
  class T_Prism_wSigmaGradu
  {
    AutoDiffDiff<2,T> u;
    AutoDiff<1,T> w;
  public:
    T_Prism_wSigmaGradu ( AutoDiffDiff<2,T> au, AutoDiff<1,T> aw) : u(au), w(aw) { ; }
    Vec<6,T> Shape() 
    { 
      Vec<3,T> sigma2d = T_SigmaGrad<2,T>(u).Shape();
      Vec<6,T> sigma(0.);
      sigma[0] = w.Value()*sigma2d[0];
      sigma[1] = w.Value()*sigma2d[1];
      sigma[5] = w.Value()*sigma2d[2];
      return sigma;
    }

    Vec<3,T> DivShape()
    {
      return Vec<3,T> (0., 0., 0);
    }

  };

  template <typename T>
  auto Prism_wSigmaGradu ( AutoDiffDiff<2,T> au, AutoDiff<1,T> aw)
  { return T_Prism_wSigmaGradu<T>(au, aw); }


  template <typename T>
  class T_Prism_wType2
  {
    AutoDiffDiff<2,T> u, v;
    AutoDiff<1,T> w;
  public:
    T_Prism_wType2 ( AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiff<1,T> aw) : u(au), v(av),  w(aw) { ; }
    Vec<6,T> Shape() 
    { 
      Vec<3,T> sigma2d = T_Type2<2,T>(u,v).Shape();
      Vec<6,T> sigma(0.);
      sigma[0] = w.Value()*sigma2d[0];
      sigma[1] = w.Value()*sigma2d[1];
      sigma[5] = w.Value()*sigma2d[2];
      return sigma;
    }

    Vec<3,T> DivShape()
    {
      Vec<2> divsigma2d = w.Value()*T_Type2<2,T>(u,v).DivShape();
      return Vec<3,T> (divsigma2d[0], divsigma2d[1], 0);
    }

  };

  template <typename T>
  auto Prism_wType2 (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiff<1,T> aw)
  { return T_Prism_wType2<T>(au, av, aw); }


  template <typename T>
  class T_Prism_wType3
  {
    AutoDiffDiff<2,T> u, v;
    AutoDiff<1,T> w;
  public:
    T_Prism_wType3 ( AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiff<1,T> aw) : u(au), v(av),  w(aw) { ; }
    Vec<6> Shape() 
    { 
      Vec<3,T> sigma2d = T_Type3<2,T>(u,v).Shape();
      Vec<6,T> sigma(0.);
      sigma[0] = w.Value()*sigma2d[0];
      sigma[1] = w.Value()*sigma2d[1];
      sigma[5] = w.Value()*sigma2d[2];
      return sigma;
    }

    Vec<3,T> DivShape()
    {
      Vec<2,T> divsigma2d = w.Value()*T_Type3<2,T>(u,v).DivShape();
      return Vec<3,T> (divsigma2d[0], divsigma2d[1], 0);
    }
  };

  template <typename T>
  auto Prism_wType3 (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiff<1,T> aw)
  { return T_Prism_wType3<T>(au, av, aw); }


  template <typename T>
  class T_Prism_wType4
  {
    AutoDiffDiff<2,T> u, v, w;
    AutoDiff<1,T> wz;
  public:
    T_Prism_wType4 ( AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiffDiff<2,T> aw, AutoDiff<1,T> awz) : u(au), v(av),  w(aw), wz(awz) { ; }
    Vec<6,T> Shape() 
    { 
      Vec<3,T> sigma2d = wz.Value()*T_Sigma_Duv_minus_uDv_w<2,T>(u,v,w).Shape();
      Vec<6,T> sigma(0.);
      sigma[0] = sigma2d[0];
      sigma[1] = sigma2d[1];
      sigma[5] = sigma2d[2];
      return sigma;
    }

    Vec<3,T> DivShape()
    {
      Vec<2,T> divsigma2d = wz.Value()*T_Sigma_Duv_minus_uDv_w<2,T>(u,v,w).DivShape();
      return Vec<3,T> (divsigma2d[0], divsigma2d[1], 0);
    }

  };

  template <typename T>
  auto Prism_wType4 (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiffDiff<2,T> aw, AutoDiff<1,T> awz)
  { return T_Prism_wType4<T>(au, av, aw, awz); }

  

  class Prism_SymRotRot_Dl2xDl1_vw
  {
    AutoDiff<2> l1,l2,v;
    AutoDiff<1> wz;
  public:
    Prism_SymRotRot_Dl2xDl1_vw ( AutoDiff<2> lam1, AutoDiff<2> lam2, AutoDiff<2> av, AutoDiff<1> awz) : l1(lam1), l2(lam2), v(av), wz(awz) { ; }
    Vec<6> Shape() 
    { 
      Vec<3> sigma2d = wz.Value()*T_SymRotRot_Dl2xDl1_v(l1,l2,v).Shape();
      Vec<6> sigma(0.);
      sigma[0] = sigma2d[0];
      sigma[1] = sigma2d[1];
      sigma[5] = sigma2d[2];
      return sigma;
    }

    Vec<3> DivShape()
    {
      Vec<2> divsigma2d = wz.Value()*T_SymRotRot_Dl2xDl1_v(l1,l2,v).DivShape();
      return Vec<3> (divsigma2d[0], divsigma2d[1], 0);
    }

  };

  template <typename T>
  class Prism_Dl1xDl3_symtensor_Dl2xDl4_u
  {
    AutoDiff<3,T> l1,l2,l3, l4;
    AutoDiff<3,T> u;
  public:
    Prism_Dl1xDl3_symtensor_Dl2xDl4_u ( AutoDiff<3,T> lam1, AutoDiff<3,T> lam2, AutoDiff<3,T> alz1, AutoDiff<3,T> alz2, AutoDiff<3,T> av) 
      : l1(lam1), l2(lam2), l3(alz1), l4(alz2), u(av) { ; }
    INLINE Vec<6,T> Shape() 
    { 
      auto rotlam1 = Cross(l1, l3);
      auto rotlam2 = Cross(l2, l4);

      Vec<6,T> sigma(0.);
      sigma[0] = u.Value()*rotlam1.DValue(0)*rotlam2.DValue(0);
      sigma[1] = u.Value()*rotlam1.DValue(1)*rotlam2.DValue(1);
      sigma[2] = u.Value()*rotlam1.DValue(2)*rotlam2.DValue(2);
      sigma[3] = 0.5*u.Value()*(rotlam1.DValue(2)*rotlam2.DValue(1) + rotlam2.DValue(2)*rotlam1.DValue(1));
      sigma[4] = 0.5*u.Value()*(rotlam1.DValue(2)*rotlam2.DValue(0) + rotlam2.DValue(2)*rotlam1.DValue(0));
      sigma[5] = 0.5*u.Value()*(rotlam1.DValue(0)*rotlam2.DValue(1) + rotlam2.DValue(0)*rotlam1.DValue(1));
      return sigma;
    }

    INLINE AutoDiff<3,T> Cross (const AutoDiff<3,T> & x,
                                const AutoDiff<3,T> & y)
    {
      T hv[3];
      hv[0] = x.DValue(1)*y.DValue(2)-x.DValue(2)*y.DValue(1);
      hv[1] = x.DValue(2)*y.DValue(0)-x.DValue(0)*y.DValue(2);
      hv[2] = x.DValue(0)*y.DValue(1)-x.DValue(1)*y.DValue(0);
      return AutoDiff<3,T> (0,hv);
    }
    Vec<3,T> DivShape()
    {
      auto lam1 = Cross(l1, l3);
      auto lam2 = Cross(l2, l4);
      return Vec<3,T> (u.DValue(0)*lam1.DValue(0)*lam2.DValue(0) + 
        0.5*(u.DValue(1)*(lam1.DValue(0)*lam2.DValue(1)+lam2.DValue(0)*lam1.DValue(1)) + u.DValue(2)*(lam1.DValue(0)*lam2.DValue(2)+lam2.DValue(0)*lam1.DValue(2))),
        u.DValue(1)*lam1.DValue(1)*lam2.DValue(1) + 
        0.5*(u.DValue(0)*(lam1.DValue(0)*lam2.DValue(1)+lam2.DValue(0)*lam1.DValue(1)) + u.DValue(2)*(lam1.DValue(1)*lam2.DValue(2)+lam2.DValue(1)*lam1.DValue(2))),
        u.DValue(2)*lam1.DValue(2)*lam2.DValue(2) + 
        0.5*(u.DValue(0)*(lam1.DValue(0)*lam2.DValue(2)+lam2.DValue(0)*lam1.DValue(2)) + u.DValue(1)*(lam1.DValue(1)*lam2.DValue(2)+lam2.DValue(1)*lam1.DValue(2)))
        );
    }

  };


  template <typename T>
  class t1_symtensor_t2_u
  {
    AutoDiff<3,T> t1, t2;
    AutoDiff<3,T> u;
  public:
    t1_symtensor_t2_u (AutoDiff<3,T> at1, AutoDiff<3,T> at2, AutoDiff<3,T> au) 
      : t1(at1), t2(at2), u(au) { ; } 
    INLINE Vec<6,T> Shape() 
    { 
      //auto rotlam1 = t1;
      //auto rotlam2 = t2;

      Vec<6,T> sigma(0.);
      sigma[0] = t1.DValue(0)*t2.DValue(0);
      sigma[1] = t1.DValue(1)*t2.DValue(1);
      sigma[2] = t1.DValue(2)*t2.DValue(2);
      sigma[3] = 0.5*(t1.DValue(2)*t2.DValue(1) + t2.DValue(2)*t1.DValue(1));
      sigma[4] = 0.5*(t1.DValue(2)*t2.DValue(0) + t2.DValue(2)*t1.DValue(0));
      sigma[5] = 0.5*(t1.DValue(0)*t2.DValue(1) + t2.DValue(0)*t1.DValue(1));
      return u.Value()*sigma;
    }

    Vec<3,T> DivShape()
    {
      T ut1 = 0.5*(u.DValue(0)*t1.DValue(0) + u.DValue(1)*t1.DValue(1) + u.DValue(2)*t1.DValue(2));
      T ut2 = 0.5*(u.DValue(0)*t2.DValue(0) + u.DValue(1)*t2.DValue(1) + u.DValue(2)*t2.DValue(2));
      return Vec<3,T> (ut1*t2.DValue(0) + ut2*t1.DValue(0),
                       ut1*t2.DValue(1) + ut2*t1.DValue(1),
                       ut1*t2.DValue(2) + ut2*t1.DValue(2));
    }

  };



  
  template <> class HDivDivFE<ET_PRISM> : public T_HDivDivFE<ET_PRISM> 
  {
  public:
    // order k+1 for certain components, for inner and boundary shapes
    // analysis from TDNNS paper for case xx1=0, zz1=xx2=zz2=1 for inner and boundary shapes
    // however, works also when boundary order is not increased.. check
    enum { incrorder_xx1 = 0};
    enum { incrorder_zz1 = 1};
    enum { incrorder_xx2 = 1};
    enum { incrorder_zz2 = 1};
    enum { incrorder_xx1_bd = 0};
    enum { incrorder_zz1_bd = 0};
    enum { incrorder_xx2_bd = 0};
    enum { incrorder_zz2_bd = 0};
    using T_HDivDivFE<ET_PRISM> :: T_HDivDivFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<2; i++)
      {
        ndof += (order_facet[i][0]+1+incrorder_zz1_bd)*(order_facet[i][0]+2+incrorder_zz1_bd)/2;
        order = max2(order, order_facet[i][0]+incrorder_zz1_bd);
      }
      for (int i=2; i<5; i++)
      {
        ndof += (order_facet[i][0]+1+incrorder_xx1_bd)*(order_facet[i][1]+1+incrorder_xx2_bd);
        order = max2(order, order_facet[i][0]+incrorder_xx2_bd);
      }
      int oi0 = order_inner[0];
      int oi2 = order_inner[2];
      int ninner = 3*((oi0+1+incrorder_xx1)*(oi0+incrorder_xx1))/2 *(oi2+1+incrorder_xx2) 
        + (oi0+1)*(oi0+2)*(oi2+1) 
        + (oi0+1+incrorder_zz1)*(oi0+2+incrorder_zz1)*(oi2-1+incrorder_zz2)/2;
      ndof += ninner; 

      order = max3(order, oi0+incrorder_zz1, oi2+incrorder_zz2);

    }

    // works only with old-style Transformation
    // does not work with CalcMappedShape
   template <typename Tx, typename TFA> 
    void T_CalcShape_Complex (TIP<3,Tx> ip, TFA & shape) const
    {
      AutoDiffDiff<2> x(ip.x.Value(),0);
      AutoDiffDiff<2> y(ip.y.Value(),1);
      AutoDiff<2> xd(ip.x.Value(),0);
      AutoDiff<2> yd(ip.y.Value(),1);
      AutoDiff<1> z(ip.z.Value(), 0);
      AutoDiffDiff<2> lami[6] ={ x,y,1-x-y,x,y,1-x-y };
      AutoDiff<2> lamid[6] ={ xd,yd,1-xd-yd,xd,yd,1-xd-yd };
      AutoDiff<1> lamiz[6] ={ 1-z,1-z,1-z,z,z,z };

      int ii = 0;
      
      int maxorder_facet =
        max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));

      const FACE * faces = ElementTopology::GetFaces(ET_PRISM);

      ArrayMem<AutoDiffDiff<2>,20> ha(maxorder_facet+2);
      ArrayMem<AutoDiffDiff<2>,20> u(order+2), v(order+3);
      ArrayMem<AutoDiff<2>,20> leg_u(order+2), leg_v(order+3);
      ArrayMem<AutoDiff<1>,20> leg_w(order+2);

      
      // Trig faces, (p+1)(p+2)/2
      for (int fa=0; fa<2; fa++)
      {
        int fav[3] ={faces[fa][0],faces[fa][1],faces[fa][2]};

        //Sort vertices  first edge op minimal vertex
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        
        leg_u.SetSize(order_facet[fa][0]+incrorder_zz1_bd+1);
        leg_v.SetSize(order_facet[fa][0]+incrorder_zz1_bd+1);
        ScaledLegendrePolynomial(order_facet[fa][0]+incrorder_zz1_bd,lamid[fav[0]]-lamid[fav[1]],1-lamid[fav[0]]-lamid[fav[1]],leg_u);
        LegendrePolynomial::Eval(order_facet[fa][0]+incrorder_zz1_bd,2 * lamid[fav[2]] - 1,leg_v);

        for(int j = 0; j <= order_facet[fa][0]+incrorder_zz1_bd; j++)
          for(int k = 0; k <= order_facet[fa][0]+incrorder_zz1_bd-j; k++)
            shape[ii++] = S_zz(leg_u[j],leg_v[k],lamiz[fav[0]]);
      }
      // quad faces -- use face bubbles of trig multiplied by leg_w
      // (px+1)(pz+1)
      for(int fa = 2; fa < 5; fa++)
      {
        int fmax = 0;
        for(int j = 1; j < 4; j++)
          if(vnums[faces[fa][j]] > vnums[faces[fa][fmax]]) fmax = j;

        int fz,ftrig;

        fz = 3 - fmax;
        ftrig = fmax^1;

        fmax = faces[fa][fmax];
        fz = faces[fa][fz];
        ftrig = faces[fa][ftrig];


        // int orderz = order_facet[fa][1];

        bool rotate = false;
        if(vnums[fz] > vnums[ftrig]) rotate = true;
        leg_w.SetSize(order_facet[fa][1]+incrorder_xx2_bd+1);
        ha.SetSize(order_facet[fa][0]+incrorder_xx1_bd+1);
        LegendrePolynomial::Eval(order_facet[fa][1]+incrorder_xx2_bd,lamiz[fmax]*2-1,leg_w);


        // edge functions are all div-free!
        IntegratedLegendreMonomialExt::CalcTrigExt(order_facet[fa][0]+incrorder_xx1_bd+2,
          lami[fmax]-lami[ftrig],1-lami[fmax]-lami[ftrig],ha);

        if(rotate)
          for(int k = 0; k <= order_facet[fa][1]+incrorder_xx2_bd; k++)
            for(int l = 0; l <= order_facet[fa][0]+incrorder_xx1_bd; l++)
            {
              shape[ii++] = Prism_wSigmaGradu(ha[l],leg_w[k]);
            }

        else
          for(int l = 0; l <= order_facet[fa][0]+incrorder_xx1_bd; l++)
            for(int k = 0; k <= order_facet[fa][1]+incrorder_xx2_bd; k++)
            {
              shape[ii++] = Prism_wSigmaGradu(ha[l],leg_w[k]);
            }


      }



      int oi = order_inner[0];
      leg_u.SetSize(oi+incrorder_zz1+1);
      leg_v.SetSize(oi+incrorder_zz1+1);
      leg_w.SetSize(oi+incrorder_xx2+1);
      u.SetSize(oi-1+incrorder_xx1+1);
      v.SetSize(oi-1+incrorder_xx1+1);

      ScaledLegendrePolynomial(oi+incrorder_zz1, lamid[0]-lamid[1], 1-lamid[0]-lamid[1], leg_u);
      LegendrePolynomial::Eval(oi+incrorder_zz1, 2*lamid[2]-1, leg_v);
      LegendrePolynomial::Eval(oi+incrorder_xx2, 2*lamiz[0]-1, leg_w);

      // ------------------------------------
      // based on elasticity-complex-based triangle shapes
      IntegratedLegendreMonomialExt::CalcTrigExt(oi-1+incrorder_xx1+2,lami[0]-lami[1],1-lami[0]-lami[1],u);
      LegendrePolynomial::EvalMult(oi-1+incrorder_xx1,2*lami[2]-1, lami[2], v);
      for(int k=0; k<=oi+incrorder_xx2; k++)
      {
        for(int i = 0; i <= oi-1+incrorder_xx1; i++)
        {
          for(int j = 0; j+i <= oi-1+incrorder_xx1; j++)
          {
            shape[ii++] = Prism_wSigmaGradu(u[i]*v[j],leg_w[k]);
            shape[ii++] = Prism_wType2(u[i],v[j],leg_w[k]);
          }
        }
        for(int i = 0; i <= oi-1+incrorder_xx1; i++)
        {
          for(int j = 0; j+i <= oi-1+incrorder_xx1; j++)
          {
            if(j > 0)
            {
              shape[ii++] = Prism_wType3(u[i],v[j],leg_w[k]);
            }
          }
        }
        for (int i = 0; i < oi+incrorder_xx1; i++)
        {
          shape[ii++] = Prism_wType4 (lami[0], -lami[1], v[i],leg_w[k]);
        }

      }

      // S_xz
      for (int i=0; i<=oi; i++)
      {
        for (int j=0; j+i<=oi; j++)
        {
          AutoDiff<2> uv = leg_u[i]*leg_v[j];
          for (int k=0; k<=oi; k++)
          {
            shape[ii++] = S_xz(0,uv, leg_w[k]);
            shape[ii++] = S_xz(1,uv, leg_w[k]);
          }
        }
      }

      // S_zz
      for(int k=0; k<=oi-2+incrorder_zz2; k++)
      {
        AutoDiff<1> bubw = leg_w[k]*lamiz[0]*(1-lamiz[0]);
        for(int i=0; i<=oi+incrorder_zz1; i++)
        {
          for(int j=0; j<=oi+incrorder_zz1-i; j++)
          {
            shape[ii++] = S_zz(leg_u[i],leg_v[j],bubw);
          }
        }
      }


    };




    // alternative to T_CalcShape, with "simpler" shape functions,
    // that are described in anisotropic paper
    // works with CalcMappedShape etc. routines, also for curved elements
    template <typename T, typename TFA> 
    void T_CalcShape/*_nocomplex*/ (TIP<3,AutoDiff<3,T>> ip, TFA & shape) const
    {
      // Tx x = ip.x, y = ip.y, z = ip.z;
      AutoDiffDiff<3,T> x = ip.x, y = ip.y, z = ip.z;
      // typedef decltype(x.Value()+x.Value()) T;
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> zz(z.Value(), &z.DValue(0));
      AutoDiff<3,T> lx[6] ={ xx, yy, 1-xx-yy, xx, yy, 1-xx-yy };
      AutoDiff<3,T> lz[6] ={ 1-zz,1-zz,1-zz,zz,zz,zz };
      int ii = 0;
      
      // int maxorder_facet =
      // max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));

      const FACE * faces = ElementTopology::GetFaces(ET_PRISM);

      ArrayMem<AutoDiff<3,T>,20> leg_u(order+2), leg_v(order+3);
      ArrayMem<AutoDiff<3,T>,20> leg_w(order+2);

      
      // Trig faces, (p+1)(p+2)/2
      for (int fa=0; fa<2; fa++)
      {
        int fav[3] ={faces[fa][0],faces[fa][1],faces[fa][2]};

        //Sort vertices  first edge op minimal vertex
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        
        leg_u.SetSize(order_facet[fa][0]+incrorder_zz1_bd+1);
        leg_v.SetSize(order_facet[fa][0]+incrorder_zz1_bd+1);
        ScaledLegendrePolynomial(order_facet[fa][0]+incrorder_zz1_bd,lx[fav[0]]-lx[fav[1]],lx[fav[0]]+lx[fav[1]],leg_u);
        LegendrePolynomial::Eval(order_facet[fa][0]+incrorder_zz1_bd,2 * lx[fav[2]] - 1,leg_v);

        for(int j = 0; j <= order_facet[fa][0]+incrorder_zz1_bd; j++)
          for(int k = 0; k <= order_facet[fa][0]+incrorder_zz1_bd-j; k++)
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[fav[0]], lx[fav[1]], lx[fav[2]], lx[fav[2]], leg_u[j]*leg_v[k]*lz[fav[0]]);
      }
      // quad faces -- use face bubbles of trig multiplied by leg_w
      // (px+1)(pz+1)
      for(int fa = 2; fa < 5; fa++)
      {
        int fmax = 0;
        for(int j = 1; j < 4; j++)
          if(vnums[faces[fa][j]] > vnums[faces[fa][fmax]]) fmax = j;

        int fz,ftrig;

        fz = 3 - fmax;
        ftrig = fmax^1;

        fmax = faces[fa][fmax];
        fz = faces[fa][fz];
        ftrig = faces[fa][ftrig];


        // int orderz = order_facet[fa][1];

        bool rotate = false;
        if(vnums[fz] > vnums[ftrig]) rotate = true;
        leg_w.SetSize(order_facet[fa][1]+incrorder_xx2_bd+1);
        LegendrePolynomial::Eval(order_facet[fa][1]+incrorder_xx2_bd,lz[fmax]*2-1,leg_w);


        ScaledLegendrePolynomial(order_facet[fa][0]+incrorder_xx1_bd, lx[fmax]-lx[ftrig], lx[fmax]+lx[ftrig], leg_u);      

        if(rotate)
          for(int k = 0; k <= order_facet[fa][1]+incrorder_xx2_bd; k++)
            for(int l = 0; l <= order_facet[fa][0]+incrorder_xx1_bd; l++)
            {
              shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lz[fmax], lz[fz], lx[fmax], lx[ftrig], leg_u[l]* leg_w[k]);
            }

        else
          for(int l = 0; l <= order_facet[fa][0]+incrorder_xx1_bd; l++)
            for(int k = 0; k <= order_facet[fa][1]+incrorder_xx2_bd; k++)
            {
              shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[fmax], lx[ftrig], lz[fmax], lz[fz], leg_u[l]* leg_w[k]);
            }


      }



      int oi = order_inner[0];
      leg_u.SetSize(oi+incrorder_zz1+1);
      leg_v.SetSize(oi+incrorder_zz1+1);
      leg_w.SetSize(oi+incrorder_xx2+1);

      ScaledLegendrePolynomial(oi+incrorder_zz1, lx[0]-lx[1], lx[0]+lx[1], leg_u);
      LegendrePolynomial::Eval(oi+incrorder_zz1, 2*lx[2]-1, leg_v);
      LegendrePolynomial::Eval(oi+incrorder_xx2, 2*lz[0]-1, leg_w);

      for(int k=0; k<=oi+incrorder_xx2; k++)
      {
        for(int i = 0; i <= oi-1+incrorder_xx1; i++)
        {
          for(int j = 0; j+i <= oi-1+incrorder_xx1; j++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[0], lx[1], lz[0], lz[0], lx[2]*leg_u[i]*leg_v[j]* leg_w[k]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[2], lx[0], lz[0], lz[0], lx[1]*leg_u[i]*leg_v[j]* leg_w[k]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[1], lx[2], lz[0], lz[0], lx[0]*leg_u[i]*leg_v[j]* leg_w[k]);
          }
        }
      }


      // S_xz
      for (int i=0; i<=oi; i++)
      {
        for (int j=0; j+i<=oi; j++)
        {
          for (int k=0; k<=oi; k++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lz[0], lx[1], lx[0], lx[0], leg_u[i]*leg_v[j]*leg_w[k]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lz[0], lx[0], lx[1], lx[1], leg_u[i]*leg_v[j]*leg_w[k]);
          }
        }
      }

      // S_zz
      for(int k=0; k<=oi-2+incrorder_zz2; k++)
      {
        AutoDiff<3,T> bubw = leg_w[k]*lz[0]*(1-lz[0]);
        for(int i=0; i<=oi+incrorder_zz1; i++)
        {
          for(int j=0; j<=oi+incrorder_zz1-i; j++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[0], lx[2], lx[1], lx[1], leg_u[i]*leg_v[j]*bubw);
          }
        }
      }


    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hdivdivfe not implementend for element type");
    }


  };


  template <> class HDivDivFE<ET_TET> : public T_HDivDivFE<ET_TET> 
  {
  public:
    using T_HDivDivFE<ET_TET> :: T_HDivDivFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<4; i++)
      {
        ndof += (order_facet[i][0]+1)*(order_facet[i][0]+2)/2;
        order = max2(order, order_facet[i][0]);
      }
      int p = order_inner[0];
      int ninner = (p+1)*(p+2)*(p+1);
      ndof += ninner; 

      order = max2(order, p);
      if (plus)
        {
          ndof += 2*(p+1)*(p+2);
          order = max2(order, p+1);
        }
    }


    template <typename T, typename TFA> 
    void T_CalcShape (TIP<3,AutoDiff<3,T>> ip, TFA & shape) const
    {
      AutoDiff<3,T> lam[4] = { ip.x, ip.y, ip.z, 1.0-ip.x-ip.y-ip.z };
      size_t ii = 0;
      
      //const FACE * faces = ElementTopology::GetFaces(ET_TET);

      /*
      ArrayMem<AutoDiff<3,T>,20> leg_u(order+2), leg_v(order+3);
      ArrayMem<AutoDiff<3,T>,20> leg_w(order+2);
      */

      typedef AutoDiff<3,T> ADT;
      STACK_ARRAY(ADT, leg_u, order+2);
      STACK_ARRAY(ADT, leg_v, order+3);
      
      for(int fa = 0; fa < 4; fa++)
      {
        int p = order_facet[fa][0];
        /*
        // int fav[3] = {faces[fa][0], faces[fa][1], faces[fa][2]};
        //Sort vertices  first edge op minimal vertex
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0], fav[1]);
        if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1], fav[2]);
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0], fav[1]);
        */
        
        IVec<4> fav = GetVertexOrientedFace (fa);

        ScaledLegendrePolynomial(p+1, lam[fav[0]]-lam[fav[1]],lam[fav[0]]+lam[fav[1]], &leg_u[0]);
        LegendrePolynomial::Eval(p+1, 2 * lam[fav[2]] - 1, &leg_v[0]);

        auto t1 = Cross(lam[fav[0]], lam[fav[2]]);
        auto t2 = Cross(lam[fav[1]], lam[fav[2]]);
        for(int j = 0; j <= p; j++)
          for(int k = 0; k+j <= p; k++)
          {
            /*
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lam[fav[0]], lam[fav[1]],
              lam[fav[2]], lam[fav[2]],leg_u[j]* leg_v[k]);
            */
            shape[ii++] = t1_symtensor_t2_u<T>(t1, t2, leg_u[j]* leg_v[k]);
          }
      }

      int oi = order_inner[0];
      /*
      leg_u.SetSize(oi+1);
      leg_v.SetSize(oi+1);
      leg_w.SetSize(oi+1);

      ScaledLegendrePolynomial(oi+1,lam[0]-lam[1],lam[0]+lam[1],leg_u);
      ScaledLegendrePolynomial(oi+1,lam[2]-lam[0]-lam[1],lam[0]+lam[1]+lam[2],leg_v);
      LegendrePolynomial::Eval(oi+1,2 * lam[3] - 1,leg_w);


      for(int k=0; k<=oi; k++)
      {
        for(int i = 0; i+k <= oi; i++)
        {
          for(int j = 0; j+i+k <= oi; j++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lam[0], lam[1], lam[2], lam[3], leg_u[i]*leg_v[k]* leg_w[j]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lam[0], lam[2], lam[1], lam[3], leg_u[i]*leg_v[j]* leg_w[k]);
          }
        }
      }

      for(int fa = 0; fa < 4; fa++)
      {
        int fav[3] = {faces[fa][0], faces[fa][1], faces[fa][2]};
        for(int k=0; k<=oi-1; k++)
        {
          for(int i = 0; i+k <= oi-1; i++)
          {
            for(int j = 0; j+i+k <= oi-1; j++)
            {
              shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lam[fav[0]],lam[fav[0]],lam[fav[1]],lam[fav[2]],
                (1-lam[fav[0]]-lam[fav[1]]-lam[fav[2]])*leg_u[i]*leg_v[k]* leg_w[j]);
            }
          }
        }

      }
      */
      auto l0 = lam[0];
      auto l1 = lam[1];
      auto l2 = lam[2];
      auto l3 = lam[3];
      auto t01 = Cross(lam[0], lam[1]);
      auto t02 = Cross(lam[0], lam[2]);
      auto t03 = Cross(lam[0], lam[3]);
      auto t12 = Cross(lam[1], lam[2]);      
      auto t13 = Cross(lam[1], lam[3]);
      auto t23 = Cross(lam[2], lam[3]);
      DubinerBasis3D::Eval
	(oi, lam[0], lam[1], lam[2],
         SBLambda([shape, &ii, t02, t13, t01, t23] (size_t nr, auto val) LAMBDA_INLINE
                  {
                    /*
                    shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lam[0], lam[1], lam[2], lam[3], val);
                    shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lam[0], lam[2], lam[1], lam[3], val);
                    */
                    shape[ii++] = t1_symtensor_t2_u<T>(t02, t13, val);
                    shape[ii++] = t1_symtensor_t2_u<T>(t01, t23, val);                    
                  }));

      DubinerBasis3D::Eval
	(oi-1 + (plus?1:0), lam[0], lam[1], lam[2],
         SBLambda([shape, &ii, t02, t03, t12, t13, t01, t23, l0,l1,l2,l3] (size_t nr, auto val) LAMBDA_INLINE
                  {
                    shape[ii++] = t1_symtensor_t2_u<T>(t01, t02, l3*val);
                    shape[ii++] = t1_symtensor_t2_u<T>(t03, t13, l2*val);
                    shape[ii++] = t1_symtensor_t2_u<T>(t23, t02, l1*val);
                    shape[ii++] = t1_symtensor_t2_u<T>(t12, t13, l0*val);
                    /*
                    shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lam[0], lam[0], lam[1], lam[2], lam[3]*val);
                    shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lam[3], lam[3], lam[0], lam[1], lam[2]*val);
                    shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lam[2], lam[2], lam[3], lam[0], lam[1]*val);
                    shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lam[1], lam[1], lam[2], lam[3], lam[0]*val);
                    */
                  }));

      /*
      if (plus)
        {
          // needs AutoDiffDiff !!!
          auto l0 = lam[0];
          auto l1 = lam[1];
          auto l2 = lam[2];
          auto l3 = lam[3];
          shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T> (l0, l1, l1*l2*l3, l2, AutoDiff<3,T>(1.0));
          shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T> (l0, l1, l1*l2*l3, l3, AutoDiff<3,T>(1.0));

          shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T> (l1, l2, l0*l2*l3, l0, AutoDiff<3,T>(1.0));
          shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T> (l1, l2, l0*l2*l3, l3, AutoDiff<3,T>(1.0));

          shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T> (l2, l3, l1*l0*l3, l0, AutoDiff<3,T>(1.0));
          shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T> (l2, l3, l1*l0*l3, l1, AutoDiff<3,T>(1.0));
          
          shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T> (l3, l0, l1*l2*l0, l1, AutoDiff<3,T>(1.0));
          shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T> (l3, l0, l1*l2*l0, l2, AutoDiff<3,T>(1.0));
        }
      */
    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hdivdivfe not implementend for element type");
    }

  };



  template <> class HDivDivFE<ET_HEX> : public T_HDivDivFE<ET_HEX> 
  {
  public:
    using T_HDivDivFE<ET_HEX> :: T_HDivDivFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<6; i++)
      {
        ndof += (order_facet[i][0]+1)*(order_facet[i][0]+1);
        order = max2(order, order_facet[i][0]+1);
      }
      int p = order_inner[0];
      int ninner = 3*p*(p+2)*(p+2) + 3*(p+2)*(p+1)*(p+1);
      ndof += ninner; 

      order = max2(order, p+1);

    }


    template <typename T, typename TFA> 
    void T_CalcShape (TIP<3,AutoDiff<3,T>> ip, TFA & shape) const
    {
      AutoDiffDiff<3,T> x = ip.x, y = ip.y, z = ip.z;
      // typedef decltype(x.Value()+x.Value()) T;            
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> zz(z.Value(), &z.DValue(0));
      AutoDiff<3,T> lx[2] ={ 1-xx, xx};
      AutoDiff<3,T> ly[2] ={ 1-yy, yy};
      AutoDiff<3,T> lz[2] ={ 1-zz, zz};
      AutoDiff<3,T> sigma[8] = {1-xx + 1-yy + 1-zz,
        xx + 1-yy + 1-zz,
        xx + yy + 1-zz,
        1-xx + yy + 1-zz,
        1-xx + 1-yy + zz,
        xx + 1-yy + zz,
        xx + yy + zz,
        1-xx + yy + zz};
      int ii = 0;
      
      // int maxorder_facet =
      //     max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));

      const FACE * faces = ElementTopology::GetFaces(ET_HEX);

      ArrayMem<AutoDiff<3,T>,20> leg_u(order+2), leg_v(order+3);
      ArrayMem<AutoDiff<3,T>,20> leg_w(order+2);
      AutoDiff<3,T> lam_face;
      
      for(int fa = 0; fa < 6; fa++)
      {
        int fmax = 0;
        lam_face = -1 + 0.25*sigma[faces[fa][0]];
        for(int j = 1; j < 4; j++)
        {
          if(vnums[faces[fa][j]] > vnums[faces[fa][fmax]]) fmax = j;
          lam_face += sigma[faces[fa][j]]*0.25;
        }

        int fz,ftrig;

        fz = 3 - fmax;
        ftrig = fmax^1;

        fmax = faces[fa][fmax];
        fz = faces[fa][fz];
        ftrig = faces[fa][ftrig];


        // int orderz = order_facet[fa][1];

        if(vnums[fz] < vnums[ftrig]) swap(fz, ftrig);
        int p = order_facet[fa][0];
        LegendrePolynomial::Eval(p, sigma[fmax] - sigma[ftrig],leg_u);
        LegendrePolynomial::Eval(p, sigma[fmax] - sigma[fz],leg_v);

        for(int k = 0; k <= p; k++)
          for(int l = 0; l <= p; l++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(0.5*(sigma[fmax] - sigma[fz]),-0.5*(sigma[fmax] - sigma[fz]),
              0.5*(sigma[fmax] - sigma[ftrig]),-0.5*(sigma[fmax] - sigma[ftrig]),leg_u[l]* leg_v[k]*lam_face);
          }



      }



      int oi = order_inner[0];
      leg_u.SetSize(oi+2);
      leg_v.SetSize(oi+2);
      leg_w.SetSize(oi+2);

      LegendrePolynomial::Eval(oi+1,sigma[0] - sigma[1],leg_u);
      LegendrePolynomial::Eval(oi+1,sigma[0] - sigma[3],leg_v);
      LegendrePolynomial::Eval(oi+1,sigma[0] - sigma[4],leg_w);

      for(int k=0; k<=oi-1; k++)
      {
        for(int i = 0; i <= oi+1; i++)
        {
          for(int j = 0; j <= oi+1; j++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[0], lx[1], lz[0], lz[1], ly[0]*ly[1]*leg_u[i]*leg_v[k]* leg_w[j]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[0], lx[1], ly[0], ly[1], lz[0]*lz[1]*leg_u[i]*leg_v[j]* leg_w[k]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(ly[0], ly[1], lz[0], lz[1], lx[0]*lx[1]*leg_u[k]*leg_v[j]* leg_w[i]);
          }
        }
      }


      // S_xz
      for (int i=0; i<=oi; i++)
      {
        for (int j=0; j<=oi; j++)
        {
          for (int k=0; k<=oi+1; k++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[0], lx[0], ly[0], lz[0], leg_u[k]*leg_v[j]*leg_w[i]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(ly[0], ly[0], lx[0], lz[0], leg_u[i]*leg_v[k]*leg_w[j]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lz[0], lz[0], lx[0], ly[0], leg_u[i]*leg_v[j]*leg_w[k]);
          }
        }
      }


    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hdivdivfe not implementend for element type");
    }


  };




  ////////////////////// SURFACE ////////////////////////////
    template <int DIM>
  class HDivDivSurfaceFiniteElement : public FiniteElement
  {
  public:
    using FiniteElement::FiniteElement;
    using FiniteElement::ndof;
    using FiniteElement::order;

    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM+1> & mip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM+1> & mip,
      BareSliceMatrix<double> shape) const = 0;

  };


  template <ELEMENT_TYPE ET> class HDivDivSurfaceFE;

  
  template <ELEMENT_TYPE ET>
  class T_HDivDivSurfaceFE : public HDivDivSurfaceFiniteElement<ET_trait<ET>::DIM>,
    public VertexOrientedFE<ET>
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
    enum { DIM_STRESS = ((DIM+2)*(DIM+1))/2 };
    
    using VertexOrientedFE<ET>::vnums;
    using HDivDivSurfaceFiniteElement<ET_trait<ET>::DIM>::ndof;
    using HDivDivSurfaceFiniteElement<ET_trait<ET>::DIM>::order;

    IVec<DIM> order_inner;


  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    T_HDivDivSurfaceFE (int aorder)
    {
      order = aorder;
      order_inner = aorder;
    }
    
    virtual ELEMENT_TYPE ElementType() const { return ET; }
    const HDivDivSurfaceFE<ET> * Cast() const { return static_cast<const HDivDivSurfaceFE<ET>*> (this); } 
    
    INLINE void SetOrderInner (IVec<DIM,int> order) { order_inner = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HDivDivSurfaceFE<ET>:: ComputeNDof not available for base class" << endl;
    }

    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM+1> & mip,
                            BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiff<DIM+1>> adp = mip;
      TIP<DIM, AutoDiffDiff<DIM+1>> addp(adp, mip.IP().FacetNr(), mip.IP().VB());
   
      Cast() -> T_CalcShape (addp, SBLambda([&] (int nr, auto val)
                                          {
                                            shape.Row(nr).Range(0,DIM_STRESS) = val.Shape();
                                          }));
    }


    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM+1> & mip,
                            BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiff<DIM+1>> adp = mip;
      TIP<DIM, AutoDiffDiff<DIM+1>> addp(adp, mip.IP().FacetNr(), mip.IP().VB());
      
      Cast() -> T_CalcShape (addp, SBLambda([&](int nr,auto val)
      {
        Vec<DIM_STRESS> vecshape = val.Shape();
        BareVector<double> matshape = shape.Row(nr);
        VecToSymMat<DIM+1> (vecshape, matshape);
      }));
    }


  };

  template <> class HDivDivSurfaceFE<ET_SEGM> : public T_HDivDivSurfaceFE<ET_SEGM> 
  {
    
  public:
    using T_HDivDivSurfaceFE<ET_SEGM> :: T_HDivDivSurfaceFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += order_inner[0]+1;
      order = max2(order,order_inner[0]);

    }
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<1,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {
      auto x = ip.x;
      AutoDiffDiff<2> ddlami[2] ={ x, 1-x };
      
      int ii = 0;
      
      ArrayMem<AutoDiffDiff<2>,20> u(order_inner[0]+2);
      
      int es = 0,ee = 1;
      if(vnums[es] > vnums[ee]) swap (es,ee);

      AutoDiffDiff<2> ls = ddlami[es],le = ddlami[ee];

      IntegratedLegendreMonomialExt::Calc(order_inner[0]+2, le-ls,u);

      for(int l = 0; l <= order_inner[0]; l++)
        shape[ii++] = SigmaGrad (u[l]);

      
    };
  };


  template <> class HDivDivSurfaceFE<ET_TRIG> : public T_HDivDivSurfaceFE<ET_TRIG> 
  {
    
  public:
    using T_HDivDivSurfaceFE<ET_TRIG> :: T_HDivDivSurfaceFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += (order_inner[0]+1+HDivDivFE<ET_PRISM>::incrorder_zz1_bd)*(order_inner[0]+2+HDivDivFE<ET_PRISM>::incrorder_zz1_bd)/2;
      order = max2(order, order_inner[0]+HDivDivFE<ET_PRISM>::incrorder_zz1_bd);
    }


    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y;
      typedef decltype(x.Value()+x.Value()) T;                  
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> lx[6] ={ xx, yy, 1-xx-yy};
      int ii = 0;

      ArrayMem<AutoDiff<3,T>,20> leg_u(order_inner[0]+2), leg_v(order_inner[0]+3);

      
        int fav[3] ={0,1,2};

        //Sort vertices  first edge op minimal vertex
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        
        leg_u.SetSize(order_inner[0]+HDivDivFE<ET_PRISM>::incrorder_zz1_bd+1);
        leg_v.SetSize(order_inner[0]+HDivDivFE<ET_PRISM>::incrorder_zz1_bd+1);
        ScaledLegendrePolynomial(order_inner[0]+HDivDivFE<ET_PRISM>::incrorder_zz1_bd,lx[fav[0]]-lx[fav[1]],lx[fav[0]]+lx[fav[1]],leg_u);
        LegendrePolynomial::Eval(order_inner[0]+HDivDivFE<ET_PRISM>::incrorder_zz1_bd,2 * lx[fav[2]] - 1,leg_v);

        for(int j = 0; j <= order_inner[0]+HDivDivFE<ET_PRISM>::incrorder_zz1_bd; j++)
          for(int k = 0; k <= order_inner[0]+HDivDivFE<ET_PRISM>::incrorder_zz1_bd-j; k++)
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[fav[0]], lx[fav[1]], lx[fav[2]], lx[fav[2]], leg_u[j]*leg_v[k]);
      }
  };


    template <> class HDivDivSurfaceFE<ET_QUAD> : public T_HDivDivSurfaceFE<ET_QUAD> 
  {
    
  public:
    using T_HDivDivSurfaceFE<ET_QUAD> :: T_HDivDivSurfaceFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += (order_inner[0]+1+HDivDivFE<ET_PRISM>::incrorder_xx1_bd)*(order_inner[1]+1+HDivDivFE<ET_PRISM>::incrorder_xx2_bd);
      order = max2(order, order_inner[0]+HDivDivFE<ET_PRISM>::incrorder_xx1_bd);
      order = max2(order, order_inner[1]+HDivDivFE<ET_PRISM>::incrorder_xx2_bd);
    }


    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      AutoDiffDiff<3> x = ip.x, z = ip.y;
      typedef decltype(x.Value()+x.Value()) T;                  
      AutoDiff<3> xx(x.Value(), &x.DValue(0));
      AutoDiff<3> zz(z.Value(), &z.DValue(0));
      AutoDiff<3> sigma[4] = {1-xx+1-zz, xx+1-zz, xx+zz, 1-xx+zz};
      int ii = 0;
      

      ArrayMem<AutoDiff<3>,20> leg_u(order_inner[0]+2);
      ArrayMem<AutoDiff<3>,20> leg_w(order_inner[1]+2);

      
      int fmax = 0;
      for(int j = 1; j < 4; j++)
        if(vnums[j] > vnums[fmax]) fmax = j;

      int f1, f2;
      f1 = (fmax+1)%4;
      f2 = (fmax+3)%4;


      if(vnums[f1] > vnums[f2])
      {
        swap(f1,f2);
      }

      LegendrePolynomial::Eval(order_inner[0],sigma[fmax] - sigma[f1],leg_u);
      LegendrePolynomial::Eval(order_inner[0],sigma[fmax] - sigma[f2],leg_w);

      for(int k = 0; k <= order_inner[0]+HDivDivFE<ET_PRISM>::incrorder_xx2_bd; k++)
        for(int l = 0; l <= order_inner[0]+HDivDivFE<ET_PRISM>::incrorder_xx1_bd; l++)
        {
          shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(0.5*(sigma[fmax]-sigma[f2]),-0.5*(sigma[fmax]-sigma[f2]),
            0.5*(sigma[fmax]-sigma[f1]),-0.5*(sigma[fmax]-sigma[f1]),leg_u[l]* leg_w[k]);
        }
                

    }
  };


  HDIVDIVFE_EXTERN template class T_HDivDivFE<ET_TRIG>;
  HDIVDIVFE_EXTERN template class T_HDivDivFE<ET_QUAD>;
  HDIVDIVFE_EXTERN template class T_HDivDivFE<ET_TET>;
  HDIVDIVFE_EXTERN template class T_HDivDivFE<ET_PRISM>;
  HDIVDIVFE_EXTERN template class T_HDivDivFE<ET_HEX>;
  HDIVDIVFE_EXTERN template class T_HDivDivSurfaceFE<ET_SEGM>;
  HDIVDIVFE_EXTERN template class T_HDivDivSurfaceFE<ET_QUAD>;
  HDIVDIVFE_EXTERN template class T_HDivDivSurfaceFE<ET_TRIG>;

}


#endif
  
