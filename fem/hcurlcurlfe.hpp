#ifndef FILE_HCURLCURLFE
#define FILE_HCURLCURLFE

/*********************************************************************/
/* File:   hcurlcurlfe.hpp                                           */
/* Author: Michael Neunteufel                                        */
/* Date:   June 2018                                                 */
/*********************************************************************/

namespace ngfem
{

  template <int DIM>
  class HCurlCurlFiniteElement : public FiniteElement
  {
  public:
    using FiniteElement::FiniteElement;
    using FiniteElement::ndof;
    using FiniteElement::order;

    virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceMatrix<> shape) const = 0;
    
    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedCurlShape (const BaseMappedIntegrationPoint & bmip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                                         BareSliceMatrix<SIMD<double>> shapes) const = 0;
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                                  BareSliceVector<> coefs,
                                  BareSliceMatrix<SIMD<double>> values) const = 0;

    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                  BareSliceMatrix<SIMD<double>> values,
                                  BareSliceVector<> coefs) const = 0;

    virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, SliceMatrix<> shape) const = 0;
    virtual void CalcDualShape (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> shape) const = 0;
    virtual void EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const = 0;
    virtual void AddDualTrans (const SIMD_BaseMappedIntegrationRule& bmir, BareSliceMatrix<SIMD<double>> values, BareSliceVector<double> coefs) const = 0;



    const FlatMatrixFixWidth<DIM*DIM> GetShape (const IntegrationPoint & ip, 
					    LocalHeap & lh) const
    {
      FlatMatrixFixWidth<DIM*DIM> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

  };

  template <int D,typename VEC,typename MAT>
  void VecToSymMat(const VEC & vec,MAT & mat)
  {
    switch(D)
    {
    case 1:
      mat(0) = vec(0);
    case 2:
      mat(0) = vec(0);
      mat(3) = vec(1);
      mat(1) = mat(2) = vec(2);
      break;
    case 3:
      mat(0) = vec(0);
      mat(4) = vec(1);
      mat(8) = vec(2);
      mat(1) = mat(3) = vec(5);
      mat(2) = mat(6) = vec(4);
      mat(5) = mat(7) = vec(3);
      break;
    }
  }

  template <int D,typename T ,typename MAT>
  Vec<D*D-D,T> SymMatToVec(MAT & mat)
  {
    if (D == 2)
      return Vec<3,T>(mat(0),mat(3),mat(2));
    else if (D == 3)
      return Vec<6,T>(mat(0),mat(4),mat(8),mat(5),mat(2),mat(1));
    else
      return T(0.0);
  }

  template <typename T>
  Mat<2,2,T> DyadProd(Vec<2,T> a, Vec<2,T> b)
  {
    return Matrix<T>({{a(0)*b(0), a(0)*b(1)}, {a(1)*b(0), a(1)*b(1)}} );
  }

  template <typename T>
  Mat<3,3,T> DyadProd(Vec<3,T> a, Vec<3,T> b)
  {
    return Matrix<T>( {{a(0)*b(0), a(0)*b(1), a(0)*b(2)}, {a(1)*b(0), a(1)*b(1), a(1)*b(2)}, {a(2)*b(0), a(2)*b(1), a(2)*b(2)}} );
  }

  template <typename T>
  Vec<6, AutoDiff<3,T>> SymDyadProd(AutoDiff<3,T> a, AutoDiff<3,T> b)
  {
    return Vec<6, AutoDiff<3,T>>(2*a.DValue(0)*b.DValue(0),2*a.DValue(1)*b.DValue(1),2*a.DValue(2)*b.DValue(2), a.DValue(1)*b.DValue(2)+a.DValue(2)*b.DValue(1), a.DValue(0)*b.DValue(2)+a.DValue(2)*b.DValue(0),a.DValue(1)*b.DValue(0)+a.DValue(0)*b.DValue(1));
  }

  template <typename T>
  Vec<6, AutoDiff<3,T>> SymDyadProd(Vec<3,T> a, Vec<3,T> b)
  {
    return Vec<6, AutoDiff<3,T>>(2*a(0)*b(0),2*a(1)*b(1),2*a(2)*b(2), a(1)*b(2)+a(2)*b(1), a(0)*b(2)+a(2)*b(0),a(1)*b(0)+a(0)*b(1));
  }

  template <typename T>
  Vec<3,AutoDiff<2,T>> SymDyadProd(AutoDiff<2,T> a, AutoDiff<2,T> b)
  {
    return Vec<3,AutoDiff<2,T>>(2*a.DValue(0)*b.DValue(0),2*a.DValue(1)*b.DValue(1),a.DValue(1)*b.DValue(0)+a.DValue(0)*b.DValue(1));
  }

  template <typename T>
  Vec<3,AutoDiff<2,T>> SymDyadProd(Vec<2,T> a, Vec<2,T> b)
  {
    return Vec<3,AutoDiff<2,T>>(2*a(0)*b(0),2*a(1)*b(1),a(1)*b(0)+a(0)*b(1));
  }

  template <typename T>
  AutoDiff<1,T> SymDyadProd(AutoDiff<1,T> a, AutoDiff<1,T> b)
  {
    return a.DValue(0)*b.DValue(0);
  }


  //------------------REGGE_SHAPE---------------------
  template <int D, typename T> class T_REGGE_Shape;
  template <typename T> class T_REGGE_Shape<1,T>
  {
    AutoDiff<1,T> u;
  public:
    T_REGGE_Shape  (AutoDiff<1,T> au) : u(au) { ; }
    Vec<1,T> Shape() { return u.Value(); }
    /*0 2
      2 1*/
    Vec<1,T> CurlShape() { return 0.0; }
  };
  
  template <typename T> class T_REGGE_Shape<2,T>
  {
    Vec<3,AutoDiff<2,T>> u;
  public:
    T_REGGE_Shape  (Vec<3,AutoDiff<2,T>> au) : u(au) { ; }
    Vec<3,T> Shape() { return Vec<3,T> (u(0).Value(), u(1).Value(), u(2).Value()); }
    /*0 2
      2 1*/
    Vec<2,T> CurlShape() { return Vec<2,T> (u(2).DValue(0)-u(0).DValue(1), u(1).DValue(0)-u(2).DValue(1)); }
  };
  
  template <typename T> class T_REGGE_Shape<3,T>
  {
    Vec<6,AutoDiff<3,T>> u;
  public:
    T_REGGE_Shape  (Vec<6,AutoDiff<3,T>> au) : u(au) { ; }
    Vec<6,T> Shape() { return Vec<6,T> (u(0).Value(), u(1).Value(), u(2).Value(), u(3).Value(), u(4).Value(), u(5).Value()); }
    /*0 5 4
      5 1 3
      4 3 2*/
    Vec<9,T> CurlShape() { return Vec<9,T> (u(4).DValue(1)-u(5).DValue(2), -u(4).DValue(0)+u(0).DValue(2), u(5).DValue(0)-u(0).DValue(1),
					    u(3).DValue(1)-u(1).DValue(2), -u(3).DValue(0)+u(5).DValue(2), u(1).DValue(0)-u(5).DValue(1),
					    u(2).DValue(1)-u(3).DValue(2), -u(2).DValue(0)+u(4).DValue(2), u(3).DValue(0)-u(4).DValue(1)); }
  };
  //---------------------------------------------------


    // ***************** EpsGrad ****************************** */
  // eps (nabla u)
  
  template <int D, typename T> class T_EpsGrad;
  template <typename T> class T_EpsGrad<2,T>
  {
    AutoDiffDiff<2,T> u;
  public:
    T_EpsGrad  (AutoDiffDiff<2,T> au) : u(au) { ; }
    Vec<3,T> Shape()
    {
      return Vec<3,T> (u.DDValue(0,0), u.DDValue(1,1), u.DDValue(0,1));
    }
    Vec<2,T> CurlShape() { return Vec<2,T> (0.0, 0.0); }
  };
  
  template <int D, typename T>
  auto EpsGrad (AutoDiffDiff<D,T> au) { return T_EpsGrad<D,T>(au); }

  // ***************** wEpsGrad ****************************** */
  // w*eps (nabla u)
  
  template <int D, typename T> class T_wEpsGrad;
  template <typename T> class T_wEpsGrad<2,T>
  {
    AutoDiffDiff<2,T> u;
    AutoDiff<1,T> w;
  public:
    T_wEpsGrad  (AutoDiffDiff<2,T> au, AutoDiff<1,T> aw) : u(au), w(aw) { ; }
    Vec<6,T> Shape()
    {
      return w.Value()*Vec<6,T> (u.DDValue(0,0), u.DDValue(1,1), u.DDValue(2,2), u.DDValue(1,2), u.DDValue(0,2), u.DDValue(0,1));
    }
    Vec<9,T> CurlShape() { return Vec<9,T> (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); }
  };
  
  template <int D, typename T>
  auto wEpsGrad (AutoDiffDiff<D,T> au, AutoDiff<1,T> aw) { return T_wEpsGrad<D,T>(au, aw); }
  
  
  // ***************** Eps_u_Gradv ****************************** */
  // eps (u nabla v)
  
  template <int D, typename T> class T_Eps_u_Gradv;
  template <typename T> class T_Eps_u_Gradv<2,T>
  {
    AutoDiffDiff<2,T> u, v;
  public:
    T_Eps_u_Gradv  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av) { ; }
    Vec<3,T> Shape() { return Vec<3,T> ((u.Value()*v.DDValue(0,0) + u.DValue(0)*v.DValue(0)),
                                        (u.Value()*v.DDValue(1,1) + u.DValue(1)*v.DValue(1)),
                                        u.Value()*v.DDValue(0,1) + 0.5 * (u.DValue(0)*v.DValue(1)+u.DValue(1)*v.DValue(0))); }
    Vec<2,T> CurlShape()
    {
      T uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      T vx = v.DValue(0), vy = v.DValue(1);
      
      /*return -0.5 * Vec<2,T> (uyy*vx - uxy*vy + uy*vxy - ux*vyy,
      -uxy*vx + uxx*vy - uy*vxx + ux*vxy);*/
      return 0.5 * Vec<2,T>(ux*vxy - uy*vxx - uxy*vx + uxx*vy,
                            ux*vyy + uxy*vy - uyy*vx - uy*vxy);
    }
  };
  
  template <int D, typename T>
  auto Eps_u_Gradv (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_Eps_u_Gradv<D,T>(au, av); }
  
  
  template <int D, typename T> class T_vEpsGradu;
  template <typename T> class T_vEpsGradu<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_vEpsGradu  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av) { ; }
    Vec<3,T> Shape() { return Vec<3,T> (u.DDValue(0,0)*v.Value(),
                                      u.DDValue(1,1)*v.Value(),  (u.DDValue(1,0)*v.Value()));}
    Vec<2,T> CurlShape()
    {
      T uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      T vx = v.DValue(0), vy = v.DValue(1);

      //return Vec<2,T> (uyy*vx- uxy*vy, uxx*vy- uxy*vx);
      return Vec<2,T> (uxy*vx - vy*uxx, uyy*vx - uxy*vy);
    }
  };
  
  template <int D, typename T>
  auto vEpsGradu (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_vEpsGradu<D,T>(au, av); }
    
  template <ELEMENT_TYPE ET> class HCurlCurlFE;

  
  template <ELEMENT_TYPE ET>
  class T_HCurlCurlFE : public HCurlCurlFiniteElement<ET_trait<ET>::DIM>,
    public VertexOrientedFE<ET>
  {
  protected:
    static constexpr int DIM = ET_trait<ET>::DIM;
    enum { DIM_STRESS = (DIM*(DIM+1))/2 };
    enum { DIM_DMAT = 7*DIM-12 };
    
    using VertexOrientedFE<ET>::vnums;
    using HCurlCurlFiniteElement<ET_trait<ET>::DIM>::ndof;
    using HCurlCurlFiniteElement<ET_trait<ET>::DIM>::order;
    

    int order_edge[ET_trait<ET>::N_EDGE];
    INT<DIM-1> order_facet[ET_trait<ET>::N_FACET];
    INT<DIM> order_inner;

    
  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    T_HCurlCurlFE (int aorder)
    {
      order = aorder;
      for (auto & of : order_facet) of = aorder;
      order_inner = aorder;

    }
    
    virtual ELEMENT_TYPE ElementType() const override { return ET; }
    const HCurlCurlFE<ET> * Cast() const { return static_cast<const HCurlCurlFE<ET>*> (this); } 
    
    INLINE void SetOrderFacet (int nr, INT<DIM-1,int> order) { order_facet[nr] = order; }
    INLINE void SetOrderEdge (int nr, int order) { order_edge[nr] = order; }
    INLINE void SetOrderInner (INT<DIM,int> order) { order_inner = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HCurlCurlFE<ET>:: ComputeNDof not available, only for ET == TRIG" << endl;
    }

   virtual void CalcShape (const IntegrationPoint & ip, SliceMatrix<> shape) const override
    {
      Vec<DIM,AutoDiff<DIM>> tip = ip;
      Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM>> (tip), SBLambda ([shape](size_t nr, auto val)      {
                                                                 
          VecToSymMat<DIM> (val.Shape(), shape.Row(nr));
      }));

    }

    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                            BareSliceMatrix<double> shapes) const override
    {

      /*
        Here I would need GetTIP for AutoDiffDiffRec
        Switch<4-DIM>
        (bmip.DimSpace()-DIM,[this, &bmip, shapes](auto CODIM)
         {
           auto & mip = static_cast<const MappedIntegrationPoint<DIM, DIM+CODIM.value>&> (bmip);
           Cast() -> T_CalcShape (GetTIP(mip),SBLambda([shapes](int nr,auto s)
                                                       {
                                                         auto val = s.Value();
                                                         VecToSymMat<val.Size()> (val, shapes.Row(nr));
                                                       }));
         });*/

      auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM>&> (bmip);
      Vec<DIM, AutoDiff<DIM>> adp = mip;
      Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM>> (adp),SBLambda([shapes](int nr,auto val)
                           {
                             VecToSymMat<DIM> (val.Shape(), shapes.Row(nr));
                           }));
      
    }

    virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, SliceMatrix<> shape) const override
    {
      shape = 0.0;
      Switch<4-DIM>
        (bmip.DimSpace()-DIM,[this, &bmip, shape](auto CODIM)
         {
           auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);

           Cast() -> CalcDualShape2 (mip, SBLambda([shape] (size_t nr, auto val)
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
           }
       });
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
           }
       });
  }

    virtual void CalcMappedCurlShape (const BaseMappedIntegrationPoint & bmip,
                                     BareSliceMatrix<double> shape) const override
    {
      auto mip = static_cast<const MappedIntegrationPoint<DIM,DIM> &>(bmip);
      Vec<DIM, AutoDiff<DIM>> adp = mip;
      TIP<DIM, AutoDiffDiff<DIM>> addp(adp);

      if (!mip.GetTransformation().IsCurvedElement()) // non-curved element
      {
        Cast() -> T_CalcShape (addp, SBLambda([&](int nr,auto val)
        {
          shape.Row(nr).AddSize(DIM_DMAT) = val.CurlShape();
        }));
      }
      else // curved element
      {
        throw Exception("CalcMappedCurlShape not implemented for curved elements!");
      }
    }

    template <int DIMSPACE>
    void CalcMappedShape2 (const SIMD_MappedIntegrationRule<DIM,DIMSPACE> & mir, 
                                 BareSliceMatrix<SIMD<double>> shapes) const
    {
      for (size_t i = 0; i < mir.Size(); i++)
        {
          auto jacI = mir[i].GetJacobianInverse();
          
          
          Vec<DIM_STRESS,SIMD<double>> hv;
          Mat<DIM,DIM,SIMD<double>> mat;
          SIMD<double> mem[DIMSPACE*DIMSPACE*DIM_STRESS];
          FlatMatrix<SIMD<double>> trans(DIMSPACE*DIMSPACE,DIM_STRESS,&mem[0]);
          for (int k = 0; k < DIM_STRESS; k++)
            {
              hv = SIMD<double>(0.0);
              hv(k) = SIMD<double>(1.0);
              VecToSymMat<DIM> (hv, mat);
              Mat<DIMSPACE,DIMSPACE,SIMD<double>> physmat = Trans(jacI) * mat * jacI;
              trans.Col(k) = physmat.AsVector();
            }
          
          
          Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir.IR()[i];
          TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          
          this->Cast() -> T_CalcShape (addp,
                                       SBLambda ([i,shapes,trans] (size_t j, auto val) 
                                                 {
                                                   shapes.Rows(j*sqr(DIMSPACE), (j+1)*sqr(DIMSPACE)).Col(i).Range(0,sqr(DIMSPACE)) = trans * val.Shape();
                                                 }));
        }
    }

      
    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                                         BareSliceMatrix<SIMD<double>> shapes) const override
    {

      Switch<4-DIM>
        (bmir.DimSpace()-DIM,[this, &bmir, shapes](auto CODIM)
         {
           auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM+CODIM.value>&> (bmir);
           this->CalcMappedShape2 (mir, shapes);
         });
    }


    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & bmir,
                                  BareSliceVector<> coefs,
                                  BareSliceMatrix<SIMD<double>> values) const override
    {
      for (size_t i = 0; i < bmir.Size(); i++)
        {
          double *pcoefs = &coefs(0);
          const size_t dist = coefs.Dist();
          
          Vec<DIM_STRESS,SIMD<double>> sum(0.0);
          Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = bmir.IR()[i];
          TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          
          Cast() -> T_CalcShape (addp,
                                 SBLambda ([&sum,&pcoefs,dist] (size_t j, auto val)
                                           {
                                             sum += (*pcoefs)*val.Shape();
                                             pcoefs += dist;
                                           }));

          Mat<DIM,DIM,SIMD<double>> summat;
          VecToSymMat<DIM> (sum, summat);
          
          Switch<4-DIM>
            (bmir.DimSpace()-DIM,[values,&bmir,i,summat](auto CODIM)
             {
               constexpr auto DIMSPACE = DIM+CODIM.value;
               auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
               auto jacI = mir[i].GetJacobianInverse();
               Mat<DIMSPACE,DIMSPACE,SIMD<double>> physmat = Trans(jacI) * summat * jacI;
               for (size_t k = 0; k < sqr(DIMSPACE); k++)
                 values(k,i) = physmat(k);
             });
        }
    }

    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & bmir,
                                  BareSliceMatrix<SIMD<double>> values,
                                  BareSliceVector<> coefs) const override
    {
       for (size_t i = 0; i < bmir.Size(); i++)
        {
          Mat<DIM,DIM,SIMD<double>> mat;
          
          Switch<4-DIM>
            (bmir.DimSpace()-DIM,[&bmir,i,&mat,values](auto CODIM)
             {
               constexpr auto DIMSPACE = DIM+CODIM.value;
               auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
               
               auto jacI = mir[i].GetJacobianInverse();
               
               Mat<DIMSPACE,DIMSPACE,SIMD<double>> physmat{};
               physmat.AsVector() = values.Col(i);
               mat = jacI * physmat * Trans(jacI);
             });
          
          Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = bmir.IR()[i];
          TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          double *pcoefs = &coefs(0);
          const size_t dist = coefs.Dist();

          Cast() -> T_CalcShape (addp,
                                 SBLambda ([mat,&pcoefs,dist] (size_t j, auto val)
                                           {
                                             Mat<DIM,DIM,SIMD<double>> mat2;
                                             VecToSymMat<DIM> (val.Shape(), mat2);
                                             
                                             *pcoefs += HSum(InnerProduct(mat,mat2));
                                             pcoefs += dist;
                                           }));
        }
    }
   
  };



  
#ifdef FILE_HCURLCURLFE_CPP
#define HCURLCURLFE_EXTERN
#else
#define HCURLCURLFE_EXTERN extern
#endif
  
  HCURLCURLFE_EXTERN template class HCurlCurlFiniteElement<2>;
  HCURLCURLFE_EXTERN template class HCurlCurlFiniteElement<3>;

  template <> class HCurlCurlFE<ET_SEGM> : public T_HCurlCurlFE<ET_SEGM> 
  {
    
  public:
    using T_HCurlCurlFE<ET_SEGM> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += order_inner[0]+1;
      order = max2(order,order_inner[0]);

    }
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<1,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x;
      typedef decltype(x.Value()+x.Value()) T;                  
      AutoDiff<1,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<1,T> ddlami[2] ={ xx, 1-xx};
      int ii = 0;

      INT<2> e = ET_trait<ET_SEGM>::GetEdgeSort (0, vnums);
      AutoDiff<1,T> ls = ddlami[e[0]],le = ddlami[e[1]];
   
      auto symdyadic = SymDyadProd(ls,le);
          
      LegendrePolynomial::EvalScaled(order_inner[0], ls-le,ls+le,
                                     SBLambda([shape,symdyadic, &ii] (size_t nr, auto val)
                                              {
                                                shape[ii++] = T_REGGE_Shape<1,T>(val*symdyadic);
                                              }));
      
    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hcurlcurlfe calcdualshape2 not implementend for element type ET_SEGM");
    }
  };

  
  template <> class HCurlCurlFE<ET_TRIG> : public T_HCurlCurlFE<ET_TRIG> 
  {
    
  public:
    using T_HCurlCurlFE<ET_TRIG> :: T_HCurlCurlFE;

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
      
      ndof += ninner;

    }
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y;
      // Tx llami[3] ={x, y, 1-x-y};
      typedef decltype(x.Value()+x.Value()) T;      
      AutoDiff<2,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<2,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<2,T> ddlami[3] ={ xx, yy, 1-xx-yy };
      int ii = 0;


      int maxorder_facet =
        max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));
      ArrayMem<Tx,20> ha(maxorder_facet+1);
      ArrayMem<Tx,20> u(order_inner[0]+2), v(order_inner[0]+2);
      
      /*for (int i = 0; i < 3; i++)
        {
          INT<2> e = ET_trait<ET_TRIG>::GetEdgeSort(i,vnums);
	  Tx ls = llami[e[0]], le = llami[e[1]];
          
          // edge functions are all curl-free!
          IntegratedLegendreMonomialExt::CalcTrigExt(maxorder_facet+2,
                                                     le-ls, 1-le-ls, ha);

          for (int l = 0; l <= order_facet[i][0]; l++)
            shape[ii++] = EpsGrad (ha[l]);
            }*/

      for (int i = 0; i < 3; i++)
        {
          INT<2> e = ET_trait<ET_TRIG>::GetEdgeSort (i, vnums);
          AutoDiff<2,T> ls = ddlami[e[0]], le = ddlami[e[1]];

          Vec<3, AutoDiff<2,T>> symdyadic = SymDyadProd(ls,le);

          LegendrePolynomial::EvalScaled(order_facet[i][0], ls-le,ls+le, SBLambda([symdyadic, &ii, shape] (size_t nr, auto val)
                            {
                              shape[ii++] = T_REGGE_Shape<2,T>(val*symdyadic);
                            }));
        }


      if (order_inner[0] > 0)
        {
	  INT<4> f = ET_trait<ET_TRIG>::GetFaceSort(0, vnums); 
	  AutoDiff<2,T> ls = ddlami[f[0]], le = ddlami[f[1]], lt = ddlami[f[2]];
	  
	  Vec<3,AutoDiff<2,T>> symdyadic1 = lt*SymDyadProd(ls,le);
	  Vec<3,AutoDiff<2,T>> symdyadic2 = ls*SymDyadProd(lt,le);
	  Vec<3,AutoDiff<2,T>> symdyadic3 = le*SymDyadProd(ls,lt);
          
	  DubinerBasis::Eval(order_inner[0]-1, ls,le,
			      SBLambda([symdyadic1,symdyadic2,symdyadic3, &ii, shape] (size_t nr, auto val)
				       {
					 shape[ii++] = T_REGGE_Shape<2,T>(2*val*symdyadic1);
					 shape[ii++] = T_REGGE_Shape<2,T>(2*val*symdyadic2);
					 shape[ii++] = T_REGGE_Shape<2,T>(2*val*symdyadic3);
				       }));
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


      if (ip.VB() == BND)
        { // facet shapes
          for (int i = 0; i < 3; i++)
            {
              int p = order_facet[i][0];
              
              if (i == facetnr)
                {             
                  INT<2> e = ET_trait<ET_TRIG>::GetEdgeSort (i, vnums);
                  
                  T xi = lam[e[0]]-lam[e[1]];
                  Vec<2,T> tauref = pnts[e[0]] - pnts[e[1]];
                  
                  
                  auto tv = mip.GetJacobian()*tauref;

                  auto tt = DyadProd(tv,tv);
                  LegendrePolynomial::Eval
                    (p, xi,
                     SBLambda([&] (size_t nr, T val)
                              {
                                shape[nr+ii] = 1/mip.GetMeasure()*val*tt;
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
              INT<4> f =  ET_trait<ET_TRIG>::GetFaceSort(0, vnums);

              DubinerBasis::Eval (p, lam[f[0]], lam[f[1]],
                                   SBLambda([&] (size_t nr, T val)
                                            {
                                              shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<2,2>(Matrix<>({{1,0},{0,0}}))*Trans(mip.GetJacobian());
                                              shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<2,2>(Matrix<>({{0,0},{0,1}}))*Trans(mip.GetJacobian());
                                              shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<2,2>(Matrix<>({{0,1},{1,0}}))*Trans(mip.GetJacobian());
                                            }));
            }
        }
    }

    
  };
  
  template <> class HCurlCurlFE<ET_QUAD> : public T_HCurlCurlFE<ET_QUAD> 
  {
    
  public:
    using T_HCurlCurlFE<ET_QUAD> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<4; i++)
      {
        ndof += order_facet[i][0]+1;
        order = max2(order, order_facet[i][0]);
      }
      int ninner = order_inner[0]*order_inner[0] + (order_inner[0]+2)*order_inner[0]*2  +1;//+ 2*order_inner[0];
      order = max2(order, order_inner[0]);
      order += 1;
      ndof += ninner;

    }
    
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y;
      // Tx llx[4] ={1-x, x, x, 1-x};
      // Tx lly[4] ={1-y, 1-y, y, y};
      typedef decltype(x.Value()+x.Value()) T;
      AutoDiff<2,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<2,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<2,T> lx[4] ={ 1-xx, xx, xx, 1-xx };
      AutoDiff<2,T> ly[4] ={ 1-yy,1-yy, yy,yy };
      AutoDiff<2,T> lami[4] = {(1-xx)*(1-yy),xx*(1-yy),xx*yy,(1-xx)*yy};  
      AutoDiff<2,T> sigma[4] = {(1-xx)+(1-yy),xx+(1-yy),xx+yy,(1-xx)+yy};  
      
      int ii = 0;

      ArrayMem<AutoDiff<2,T>,20> v(order+2), u(order+2);
      
     
      for (int i = 0; i < 4; i++)
        {
          INT<2> e = ET_trait<ET_QUAD>::GetEdgeSort (i, vnums);
          AutoDiff<2,T> xi  = sigma[e[1]]-sigma[e[0]];
          AutoDiff<2,T> lam_e = lami[e[0]]+lami[e[1]];  
          Vec<3, AutoDiff<2,T>> symdyadic = SymDyadProd(xi,xi);


          IntLegNoBubble::
            EvalMult (order_edge[i], 
                      xi, lam_e, SBLambda ([&](int i, auto val)
                                           {
                                             shape[ii++] = T_REGGE_Shape<2,T>(val*symdyadic);
                                           }));
        }



      int oi = order_inner[0];

      Vec<3, AutoDiff<2,T>> symdyadic = SymDyadProd(Vec<2,T>(2,0),Vec<2,T>(0,2)); //(0,0,  0,1) * P(y) * P(x)
      AutoDiff<2,T> eta = ly[2]-ly[1];
      AutoDiff<2,T> xi = lx[1]-lx[0];
      LegendrePolynomial (oi, eta, v);
      LegendrePolynomial (oi, xi, u);

      for (int i = 0; i <= oi; i++)
        for (int j = 0; j <= oi; j++)
          {
            shape[ii++] = T_REGGE_Shape<2,T>(u[i]*v[j]*symdyadic);
          }

      
      auto symdyad = lx[1]*lx[0]*SymDyadProd(Vec<2,T>(0,1),Vec<2,T>(0,1));//x*(1-x)*(0,0,  0,1) * P(y) * P(x)
      for (int i = 0; i < oi; i++)
        for (int j = 0; j <= oi; j++)
          {
            shape[ii++] = T_REGGE_Shape<2,T>(u[i]*v[j]*symdyad);
          }

      symdyad = ly[2]*ly[1]*SymDyadProd(Vec<2,T>(1,0),Vec<2,T>(1,0)); //y*(1-y)*(1,0,  0,0) * P(x) * P(y)
      
      for (int j = 0; j < oi; j++)
        for (int i = 0; i <= oi; i++)
          {
            shape[ii++] = T_REGGE_Shape<2,T>(u[i]*v[j]*symdyad);
          }

      //old version
      //ArrayMem<Tx,20> u(order+2);
      /*for (int i = 0; i < 4; i++)
        {
          INT<2> e = ET_trait<ET_QUAD>::GetEdgeSort (i, vnums);
          Tx xi = llx[e[1]]+lly[e[1]]-llx[e[0]]-lly[e[0]];
          Tx eta = llx[e[0]]*lly[e[0]]+llx[e[1]]*lly[e[1]];

	  IntegratedLegendreMonomialExt::Calc(order_facet[i][0]+2,xi,u);

          
          for (int l = 0; l <= order_facet[i][0]; l++)
            shape[ii++] = Eps_u_Gradv (eta, u[l]);
        }
      
      IntegratedLegendreMonomialExt::Calc(oi+3,llx[0]-llx[1],u);
      IntegratedLegendreMonomialExt::Calc(oi+3,lly[0]-lly[2],v);

      for(int i = 0; i <= oi-1; i++)
        for(int j = 0; j <= oi-1; j++)
          shape[ii++] = EpsGrad(u[i]*v[j]);

      for(int i = 0; i <= oi+1; i++)
        for(int j = 0; j <= oi-1; j++)
        {
          shape[ii++] = vEpsGradu(u[i],v[j]);
          shape[ii++] = vEpsGradu(v[i],u[j]);
        }
      shape[ii++] = Eps_u_Gradv(lx[0], ly[0]);

      for(int i = 0; i <= oi-1; i++)
      {
        shape[ii++] = Eps_u_Gradv(u[i], ly[0]);
        shape[ii++] = Eps_u_Gradv(v[i], lx[0]);
        }*/

      
      
    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      auto & ip = mip.IP();
      typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T; 

      T x = ip(0), y = ip(1);
      T lx[4] = { 1-x, x, x, 1-x };
      T ly[4] = { 1-y, 1-y, y, y };
      // T lam[4] = { 1-x-y+x*y, x*(1-y), x*y, y*(1-x) };
      T sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};

      /*Vec<2,AutoDiff<2,T>> adip = ip;
      auto tip = TIP<2,AutoDiffDiff<2,T>>(adip);
      AutoDiffDiff<2,T> xxx = tip.x, yyy = tip.y;
      AutoDiff<2,T> xx(xxx.Value(), &xxx.DValue(0));
      AutoDiff<2,T> yy(yyy.Value(), &yyy.DValue(0));
      AutoDiff<2,T> lami[4] = {(1-xx)*(1-yy),xx*(1-yy),xx*yy,(1-xx)*yy};  
      AutoDiff<2,T> sigma[4] = {(1-xx)+(1-yy),xx+(1-yy),xx+yy,(1-xx)+yy}; */ 
      
      Vec<2,T> pnts[4] = {  { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
      int facetnr = ip.FacetNr();

      int ii = 0;
      
      ArrayMem<T,20> v(order+2), u(order+2);

      //Mat<2,2,T> tmp(0.0);

      
      if (mip.IP().VB() == BND)
        { // facet shapes
          for (int i = 0; i < 4; i++)
            {
              int p = order_facet[i][0];
              
              if (i == facetnr)
                {             
                  INT<2> e = ET_trait<ET_QUAD>::GetEdgeSort (i, vnums);
                  
                  //T xi = lam[e[0]]-lam[e[1]];
                  T xi  = sigma[e[1]]-sigma[e[0]];
                  Vec<2,T> tauref = pnts[e[0]] - pnts[e[1]];
                  
                  
                  auto tv = mip.GetJacobian()*tauref;

                  auto tt = DyadProd(tv,tv);
                  LegendrePolynomial::Eval
                    (p, xi,
                     SBLambda([&] (size_t nr, T val)
                              {
                                shape[nr+ii] = 1/mip.GetMeasure()*val*tt;
                                }));
                  /*INT<2> e = ET_trait<ET_QUAD>::GetEdgeSort (i, vnums);
                  AutoDiff<2,T> xi  = sigma[e[1]]-sigma[e[0]];
                  AutoDiff<2,T> lam_e = lami[e[0]]+lami[e[1]];  
                  Vec<3, AutoDiff<2,T>> symdyadic = SymDyadProd(xi,xi);


                  IntLegNoBubble::
                    EvalMult (order_edge[i], 
                              xi, lam_e, SBLambda ([&](int nr, auto val)
                                                   {
                                                     VecToSymMat<2>(T_REGGE_Shape<2,T>(val*symdyadic).Shape(),tmp);
                                                     shape[nr + ii] = mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
                                                     }));*/
                  /*AutoDiff<2,T> xi  = sigma[e[1]]-sigma[e[0]];
                  AutoDiff<2,T> lam_e = lami[e[0]]+lami[e[1]];  
                  Vec<3, AutoDiff<2,T>> symdyadic = SymDyadProd(xi,xi);


                  IntLegNoBubble::
                    EvalMult (p,xi, lam_e, SBLambda ([&](int nr, auto val)
                                                   {
                                                     VecToSymMat<2>(T_REGGE_Shape<2,T>(val*symdyadic).Shape(),tmp);
                                                     shape[nr + ii] = 1/mip.GetMeasure()*tmp;
                                                     }));*/
                }
              ii += (p+1);
            }
        }
      else
        {
          for (int i = 0; i < 4; i++)
            ii += order_facet[i][0]+1;
        }
      
      if (mip.IP().VB() == VOL)
        {
          auto p = order_inner[0];
         
          T eta = ly[2]-ly[1];
          T xi = lx[1]-lx[0];
          LegendrePolynomial (p, eta, v);
          LegendrePolynomial (p, xi, u);
          
          for (int i = 0; i <= p; i++)
            for (int j = 0; j <= p; j++)
              {
                shape[ii++] = 1/mip.GetMeasure()*u[i]*v[j]*mip.GetJacobian()*Mat<2,2>(Matrix<>({{0,1},{1,0}}))*Trans(mip.GetJacobian());
              }
          
          
          //auto symdyad = lx[1]*lx[0]*SymDyadProd(Vec<2,T>(0,1),Vec<2,T>(0,1));//x*(1-x)*(0,0,  0,1) * P(y) * P(x)
          for (int i = 0; i < p; i++)
            for (int j = 0; j <= p; j++)
              {
                shape[ii++] = 1/mip.GetMeasure()*u[i]*v[j]*mip.GetJacobian()*Mat<2,2>(Matrix<>({{0,0},{0,1}}))*Trans(mip.GetJacobian());
              }
          
          //symdyad = ly[2]*ly[1]*SymDyadProd(Vec<2,T>(1,0),Vec<2,T>(1,0)); //y*(1-y)*(1,0,  0,0) * P(x) * P(y)
          
          for (int j = 0; j < p; j++)
            for (int i = 0; i <= p; i++)
              {
                shape[ii++] = 1/mip.GetMeasure()*u[i]*v[j]*mip.GetJacobian()*Mat<2,2>(Matrix<>({{1,0},{0,0}}))*Trans(mip.GetJacobian());
              }
      
          //INT<4> f = ET_trait<ET_QUAD>::GetFaceSort(0, vnums);
          
          /*IntegratedLegendreMonomialExt::Calc(p+3,lx[0]-lx[1],u);
          IntegratedLegendreMonomialExt::Calc(p+3,ly[0]-ly[2],v);

          Mat<2,2,T> tmp;
          
          for(int i = 0; i <= p-1; i++)
              for(int j = 0; j <= p-1; j++)
                {
                  VecToSymMat<2>(EpsGrad(u[i]*v[j]).Shape(),tmp);
                  shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
                }
          
          for(int i = 0; i <= p+1; i++)
              for(int j = 0; j <= p-1; j++)
                {
                  VecToSymMat<2>(vEpsGradu(u[i],v[j]).Shape(),tmp);
                  shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
                  VecToSymMat<2>(vEpsGradu(v[i],u[j]).Shape(),tmp);
                  shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
                }

          VecToSymMat<2>(Eps_u_Gradv(lx[0], ly[0]).Shape(),tmp);
          shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
          
          for(int i = 0; i <= p-1; i++)
            {
              VecToSymMat<2>(Eps_u_Gradv(u[i], ly[0]).Shape(),tmp);
              shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
              VecToSymMat<2>(Eps_u_Gradv(v[i], lx[0]).Shape(),tmp);
              shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
              }*/
        }
    }
  };


  template <> class HCurlCurlFE<ET_PRISM> : public T_HCurlCurlFE<ET_PRISM> 
  {
  public:
    enum { incrorder_xx1 = 0};
    enum { incrorder_zz1 = 0};
    enum { incrorder_xx2 = 0};
    enum { incrorder_zz2 = 0};
    enum { incrorder_xx1_bd = 0};
    enum { incrorder_zz1_bd = 0};
    enum { incrorder_xx2_bd = 0};
    enum { incrorder_zz2_bd = 0};
    using T_HCurlCurlFE<ET_PRISM> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;

      for (int i=0; i < 9; i++)
        {
          ndof += order_edge[i]+1;
          order = max2(order,order_edge[i]);
        }

      for (int i=0; i<2; i++)
        {
          ndof += 3*(order_facet[i][0])*(order_facet[i][0]+1)/2;
          order = max2(order, order_facet[i][0]);
        }

      for (int i=2; i<5; i++)
        {
          ndof += order_facet[i][0]*order_facet[i][0] + (order_facet[i][0]+2)*order_facet[i][0]*2 +1;
          order = max2(order, order_facet[i][0]);
        }
      int p = order_inner[0];
      int ninner =  p > 0 ? (3*p*(p)*(p+1)/2 + (p+1)*(2*p*(p-1)/2 + (p+2)*(p-1)/2)) : 0;
              
      ndof += ninner;

      order = 1+max2(order, p);
    }
    
    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y, z = ip.z;
      // Tx llx[6] ={ x, y, 1-x-y, x, y, 1-x-y };
      // Tx llz[6] ={ 1-z,1-z,1-z,z,z,z };
      typedef decltype(x.Value()+x.Value()) T;
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> zz(z.Value(), &z.DValue(0));
      AutoDiff<3,T> lx[6] ={ xx, yy, 1-xx-yy, xx, yy, 1-xx-yy };
      AutoDiff<3,T> lz[6] ={ 1-zz,1-zz,1-zz,zz,zz,zz };

      int ii = 0;
      

      const FACE * faces = ElementTopology::GetFaces(ET_PRISM);

      ArrayMem<AutoDiff<3,T>,20> leg_u(order+2), leg_v(order+3);
      ArrayMem<AutoDiff<3,T>,20> leg_w(order+2);

      //horizontal edge shapes
      for (int i = 0; i < 6; i++)
        {
          INT<2> e = ET_trait<ET_PRISM>::GetEdgeSort (i, vnums);
          AutoDiff<3,T> ls = lx[e[0]], le = lx[e[1]], lm = lz[e[0]];

          Vec<6, AutoDiff<3,T>> symdyadic = lm*SymDyadProd(ls,le);

          LegendrePolynomial::EvalScaled(order_edge[i], ls-le,ls+le, SBLambda([symdyadic, &ii, shape] (size_t nr, auto val)
                            {
                              shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic);
                            }));
        }



      //vertical edge shapes
      for (int i = 6; i < 9; i++)
        {
          INT<2> e = ET_trait<ET_PRISM>::GetEdgeSort (i, vnums);
          AutoDiff<3,T> ls = lx[e[0]], lm1 = lz[e[0]], lm2 = lz[e[1]];
          Vec<6, AutoDiff<3,T>> symdyadic = ls*SymDyadProd(lm1,lm1);
          LegendrePolynomial (order_edge[i],lm1-lm2, leg_v);

          for (int j=0; j <= order_edge[i]; j++)
            shape[ii++] = T_REGGE_Shape<3,T>(leg_v[j]*symdyadic);
        }
      


      //horizontal face shaps
      for(int fa = 0; fa < 2; fa++)
        {
          if (order_facet[fa][0] > 0)
            {
              INT<4> f = ET_trait<ET_PRISM>::GetFaceSort(fa, vnums);
              AutoDiff<3,T> ls = lx[f[0]], le = lx[f[1]], lt = lx[f[2]], lm = lz[f[0]];
              
              Vec<6, AutoDiff<3,T>> symdyadic1 = lm*lt*SymDyadProd(ls,le);
              Vec<6, AutoDiff<3,T>> symdyadic2 = lm*ls*SymDyadProd(lt,le);
              Vec<6, AutoDiff<3,T>> symdyadic3 = lm*le*SymDyadProd(ls,lt);
              
              DubinerBasis::Eval(order_facet[fa][0]-1, ls,le,
                                 SBLambda([symdyadic1,symdyadic2,symdyadic3, &ii, shape] (size_t nr, auto val)
                                          {
                                            shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic1);
                                            shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic2);
                                            shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic3);
                                          }));
            }
        }


      //vertical face shaps
      for(int fa = 2; fa < 5; fa++)
        {
          int of = order_facet[fa][0];
          //INT<4> f = ET_trait<ET_PRISM>::GetFaceSort (fa, vnums);
          
          int fmax = 0;
          for(int j = 1; j < 4; j++)
            if(vnums[faces[fa][j]] > vnums[faces[fa][fmax]]) fmax = j;
          
          int fz,ftrig;
          fz = 3 - fmax;
          ftrig = fmax^1;          
          fmax = faces[fa][fmax];
          fz = faces[fa][fz];
          ftrig = faces[fa][ftrig];
          
          AutoDiff<3,T> eta = lz[fz]-lz[fmax];
          AutoDiff<3,T> xi = lx[ftrig]-lx[fmax];

          LegendrePolynomial (of, eta, leg_v);
          LegendrePolynomial (of, xi, leg_u);

          auto W = uDv_minus_vDu(lx[ftrig],lx[fmax]);
          Vec<6, AutoDiff<3,T>> symdyadic = SymDyadProd(GetGradient(eta),W.Value());   //^= (0,1, 1,0) * P(y)
          for (int j = 0; j <= of; j++)
            shape[ii++] = T_REGGE_Shape<3,T>(leg_v[j]*symdyadic);

          if (of > 0)
            {
              IntLegNoBubble::EvalScaledMult (of-1, xi, lx[ftrig]+lx[fmax], lx[ftrig]*lx[fmax],
                                SBLambda ([&](int nr, auto val)
                                          {
                                            auto tmp = Du (val);
                                            auto symdyadic = SymDyadProd(GetGradient(eta),tmp.Value());
                                            for (int j = 0; j <= of; j++)
                                              shape[ii++] = T_REGGE_Shape<3,T>(leg_v[j]*symdyadic);
                                          }));
            }
          auto symdyad = lx[ftrig]*lx[fmax]*SymDyadProd(eta,eta);  //^= x*(1-x)*(0,0, 0,1) * P(x) * P(y)
          for (int i = 0; i < of; i++)
            for (int j = 0; j <= of; j++)
              {
                shape[ii++] = T_REGGE_Shape<3,T>(leg_u[i]*leg_v[j]*symdyad);
              }

          symdyad = lz[fz]*lz[fmax]*SymDyadProd(lx[ftrig],lx[fmax]);    //^= y*(1-y)*(1,0, 0,0) * P(x)*P(y)
          for (int j = 0; j < of; j++)
            for (int i = 0; i <= of; i++)
              {
                shape[ii++] = T_REGGE_Shape<3,T>(leg_u[i]*leg_v[j]*symdyad);
              }
        }
      
      //inner shapes
      int p = order_inner[0];
      if (p > 0)
        {
          
          INT<4> f = ET_trait<ET_PRISM>::GetFaceSort(0, vnums);

          AutoDiff<3,T> ls = lx[f[0]], le = lx[f[1]], lt = lx[f[2]], lm = lz[0], ln = lz[3];

          Vec<6, AutoDiff<3,T>> symdyadic1 = lm*ln*lt*SymDyadProd(ls,le);
          Vec<6, AutoDiff<3,T>> symdyadic2 = lm*ln*ls*SymDyadProd(lt,le);
          Vec<6, AutoDiff<3,T>> symdyadic3 = lm*ln*le*SymDyadProd(ls,lt);

          AutoDiff<3,T> eta = lz[0]-lz[4];
          LegendrePolynomial (p, eta, leg_w);
                
          DubinerBasis::Eval(p-1, ls,le,
                             SBLambda([symdyadic1,symdyadic2,symdyadic3, &ii, shape,p,leg_w] (size_t nr, auto val)
                                      {
                                        for(int j=0; j < p; j++)
                                          {
                                            shape[ii++] = T_REGGE_Shape<3,T>(leg_w[j]*val*symdyadic1);
                                            shape[ii++] = T_REGGE_Shape<3,T>(leg_w[j]*val*symdyadic2);
                                            shape[ii++] = T_REGGE_Shape<3,T>(leg_w[j]*val*symdyadic3);
                                          }
                                      }));
          

          if(p > 1)
            {
              Vec<6, AutoDiff<3,T>> symdyadic = ls*le*lt*SymDyadProd(eta,eta);
              DubinerBasis::Eval(p-2, ls,le,
                                 SBLambda([symdyadic, &ii, shape,p,leg_w] (size_t nr, auto val)
                                      {
                                        for(int j=0; j <= p; j++)
                                          {
                                            shape[ii++] = T_REGGE_Shape<3,T>(val*leg_w[j]*symdyadic);
                                          }
                                      }));

              DubinerBasis::EvalMult(p-2, ls, le,ls*le*lt, 
                                     SBLambda
                                     ([&](int nr, auto val)
                                      {
                                        auto tmp = Du(val);
                                        Vec<6, AutoDiff<3,T>> symdyadic = SymDyadProd(tmp.Value(),GetGradient(eta));
                                        for(int j=0; j <= p; j++)
                                          {
                                            shape[ii++] = T_REGGE_Shape<3,T>(leg_w[j]*symdyadic);
                                          }
                                      }));

              auto xi  = ls-le;
            	
              TrigShapesInnerLegendre::CalcSplitted(p+1, xi, lt, leg_u,leg_v);
              
              // other combination
              for (int j = 0; j < p-1; j++)
                for (int k = 0; k < p-1-j; k++)
                  {
                    auto tmp = uDv_minus_vDu (leg_v[k], leg_u[j]);
                    auto symdyadic = SymDyadProd(tmp.Value(),GetGradient(eta));
                    for(int l=0; l <= p; l++)
                      {
                        shape[ii++] = T_REGGE_Shape<3,T>(leg_w[l]*symdyadic);
                      }
                  }
              // rec_pol * Nedelec0 
              for (int j = 0; j < p-1; j++)
                {
                  auto tmp = wuDv_minus_wvDu (le, ls, leg_v[j]);
                  auto symdyadic = SymDyadProd(tmp.Value(),GetGradient(eta));
                  for(int l=0; l <= p; l++)
                      {
                        shape[ii++] = T_REGGE_Shape<3,T>(leg_w[l]*symdyadic);
                      }
                }
              

            }
        }

    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hcurlcurlfe calcdualshape2 not implementend for element type ET_PRISM");
    }

  };


  
  template <> class HCurlCurlFE<ET_TET> : public T_HCurlCurlFE<ET_TET> 
  {
  public:
    using T_HCurlCurlFE<ET_TET> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;

      for (int i=0; i<6; i++)
      {
        ndof += order_edge[i]+1;
        order = max2(order, order_edge[i]);
      }
      
      for (int i=0; i<4; i++)
      {
        ndof += 3*(order_facet[i][0])*(order_facet[i][0]+1)/2;
        order = max2(order, order_facet[i][0]);
      }
     
      int p = order_inner[0];
      int ninner = p > 1 ? 6*(p+1)*(p)*(p-1)/6 : 0;
      ndof += ninner; 

      order = max2(order, p);
    }


    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y, z = ip.z;
      typedef decltype(x.Value()+x.Value()) T;      
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> zz(z.Value(), &z.DValue(0));
      AutoDiff<3,T> lam[4] = {xx, yy, zz, 1-xx-yy-zz};
      int ii = 0;

      for (int i = 0; i < 6; i++)
        {
          INT<2> e = ET_trait<ET_TET>::GetEdgeSort (i, vnums);
          AutoDiff<3,T> ls = lam[e[0]], le = lam[e[1]];

          Vec<6, AutoDiff<3,T>> symdyadic = SymDyadProd(ls,le);

          LegendrePolynomial::EvalScaled(order_edge[i], ls-le,ls+le, SBLambda([symdyadic, &ii, shape] (size_t nr, auto val)
                            {
                              shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic);
                            }));
        }

      
      for(int fa = 0; fa < 4; fa++)
        {
          if (order_facet[fa][0] > 0)
            {
              INT<4> f = ET_trait<ET_TET>::GetFaceSort(fa, vnums);
              AutoDiff<3,T> ls = lam[f[0]], le = lam[f[1]], lt = lam[f[2]];
              
              Vec<6, AutoDiff<3,T>> symdyadic1 = lt*SymDyadProd(ls,le);
              Vec<6, AutoDiff<3,T>> symdyadic2 = ls*SymDyadProd(lt,le);
              Vec<6, AutoDiff<3,T>> symdyadic3 = le*SymDyadProd(ls,lt);
              
              DubinerBasis::Eval(order_facet[fa][0]-1, ls,le,
                                  SBLambda([symdyadic1,symdyadic2,symdyadic3, &ii, shape] (size_t nr, auto val)
                                           {
                                             shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic1);
                                             shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic2);
                                             shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic3);
                                           }));
            }
        }

      if (order_inner[0] > 1)
        {
          AutoDiff<3,T> li = lam[0], lj = lam[1], lk = lam[2], ll = lam[3];

          Vec<6, AutoDiff<3,T>> symdyadic1 = li*lj*SymDyadProd(lk,ll);
          Vec<6, AutoDiff<3,T>> symdyadic2 = lj*lk*SymDyadProd(ll,li);
          Vec<6, AutoDiff<3,T>> symdyadic3 = lk*ll*SymDyadProd(li,lj);
          Vec<6, AutoDiff<3,T>> symdyadic4 = ll*li*SymDyadProd(lj,lk);
          Vec<6, AutoDiff<3,T>> symdyadic5 = li*lk*SymDyadProd(lj,ll);
          Vec<6, AutoDiff<3,T>> symdyadic6 = lj*ll*SymDyadProd(li,lk);
          
          DubinerBasis3D::Eval (order_inner[0]-2, lam[0], lam[1], lam[2], SBLambda([&ii, shape, symdyadic1, symdyadic2, symdyadic3, symdyadic4, symdyadic5, symdyadic6](size_t j, auto val)
                                                                                   {
                                                                                     shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic1);
                                                                                     shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic2);
                                                                                     shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic3);
                                                                                     shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic4);
                                                                                     shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic5);
                                                                                     shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic6);                                                                                                                                                   }));
        }
    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      auto & ip = mip.IP();
      typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
      T x = ip(0), y = ip(1), z = ip(2);
      T lam[4] = { x, y, z, 1-x-y-z };
      Vec<3,T> pnts[4] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } , { 0, 0, 0 } };
      int facetnr = ip.FacetNr();

      int ii = 0;

      if (ip.VB() == BBND)
        { // facet shapes
          for (int i = 0; i < 6; i++)
            {
              int p = order_edge[i];
              
              if (i == facetnr)
                {             
                  INT<2> e = ET_trait<ET_TET>::GetEdgeSort (i, vnums);
                  
                  T xi = lam[e[1]]-lam[e[0]];
                  Vec<3,T> tauref = pnts[e[1]] - pnts[e[0]];
                  Vec<3,T> tau = mip.GetJacobian()*tauref;
                  Mat<3,3,T> tt = DyadProd(tau,tau);
                  LegendrePolynomial::Eval
                    (p, xi,
                     SBLambda([&] (size_t nr, T val)
                              {
                                shape[nr+ii] = 1/mip.GetMeasure()*val*tt;
                              }));
                }
              ii += (p+1);
            }
        }
      else
        {
          for (int i = 0; i < 6; i++)
            ii += order_edge[i]+1;
        }
      if (ip.VB() == BND)
        {
          for (int i = 0; i < 4; i++)
            {
              auto p = order_facet[i][0]-1;
              if( p >= 0 && i == facetnr )
                {
                  INT<4> fav = ET_trait<ET_TET>:: GetFaceSort(facetnr, vnums);
                  Vec<3,T> adxi = pnts[fav[0]] - pnts[fav[2]];
                  Vec<3,T> adeta = pnts[fav[1]] - pnts[fav[2]];
                  T xi = lam[fav[0]];
                  T eta = lam[fav[1]];
                  
                  Matrix<T> F(3,2);
                  F.Cols(0,1) = adxi;
                  F.Cols(1,2) = adeta;
		 
                  Matrix<T> Ftmp(2,2);
                  Ftmp = Trans(F)*F;
                  auto det = sqrt(Ftmp(0,0)*Ftmp(1,1)-Ftmp(1,0)*Ftmp(0,1));
                                              
                  DubinerBasis::Eval (p, xi, eta,
                                       SBLambda([&] (size_t nr, T val)
                                                {
                                                  shape[ii++] = 1/(det*mip.GetMeasure())*val*Mat<3,3,T>(mip.GetJacobian()*F*Matrix<>({{1,0},{0,0}})*Trans(mip.GetJacobian()*F));
                                                  shape[ii++] = 1/(det*mip.GetMeasure())*val*Mat<3,3,T>(mip.GetJacobian()*F*Matrix<>({{0,0},{0,1}})*Trans(mip.GetJacobian()*F));
                                                  shape[ii++] = 1/(det*mip.GetMeasure())*val*Mat<3,3,T>(mip.GetJacobian()*F*Matrix<>({{0,1},{1,0}})*Trans(mip.GetJacobian()*F));
                                                }));
                }
              else
                ii += 3*(order_facet[i][0])*(order_facet[i][0]+1)/2;
            }
        }
      else
        {
          for (int i = 0; i < 4; i++)
            ii += 3*(order_facet[i][0])*(order_facet[i][0]+1)/2;
        }
      
      if (ip.VB() == VOL && order_inner[0] >= 2)
        {
          DubinerBasis3D::Eval (order_inner[0]-2, lam[0], lam[1], lam[2], SBLambda([&](size_t j, T val)
                               {
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{1,0,0},{0,0,0},{0,0,0}}))*Trans(mip.GetJacobian());
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{0,0,0},{0,1,0},{0,0,0}}))*Trans(mip.GetJacobian());
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{0,0,0},{0,0,0},{0,0,1}}))*Trans(mip.GetJacobian());
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{0,0,0},{0,0,1},{0,1,0}}))*Trans(mip.GetJacobian());
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{0,0,1},{0,0,0},{1,0,0}}))*Trans(mip.GetJacobian());
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{0,1,0},{1,0,0},{0,0,0}}))*Trans(mip.GetJacobian());
                               }));
          
        }
    }
  };
  


  template <> class HCurlCurlFE<ET_HEX> : public T_HCurlCurlFE<ET_HEX> 
  {
  public:
    using T_HCurlCurlFE<ET_HEX> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i < 12; i++)
        {
          ndof += order_edge[i]+1;
          order = max2(order,order_edge[i]);
        }
      for (int i=0; i<6; i++)
      {
        ndof += order_facet[i][0]*order_facet[i][0] + 2*(order_facet[i][0]+2)*order_facet[i][0]+1;
        order = max2(order, order_facet[i][0]);
      }
      int p = order_inner[0];
      ndof += 3*(p*(p+1)*(p+1) + p*p*(p+1) );

      order = 1 + max2(order, p);
    }


    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y, z = ip.z;
      typedef decltype(x.Value()+x.Value()) T;            
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> zz(z.Value(), &z.DValue(0));
      AutoDiff<3,T> lx[2] ={ 1-xx, xx};
      AutoDiff<3,T> ly[2] ={ 1-yy, yy};
      AutoDiff<3,T> lz[2] ={ 1-zz, zz};
      AutoDiff<3,T> lami[8]={(1-xx)*(1-yy)*(1-zz),xx*(1-yy)*(1-zz),xx*yy*(1-zz),(1-xx)*yy*(1-zz),
                  (1-xx)*(1-yy)*zz,xx*(1-yy)*zz,xx*yy*zz,(1-xx)*yy*zz}; 
      AutoDiff<3,T> sigma[8]={(1-xx)+(1-yy)+(1-zz),xx+(1-yy)+(1-zz),xx+yy+(1-zz),(1-xx)+yy+(1-zz),
                   (1-xx)+(1-yy)+zz,xx+(1-yy)+zz,xx+yy+zz,(1-xx)+yy+zz};
      int ii = 0;
      
      const FACE * faces = ElementTopology::GetFaces(ET_HEX);

      ArrayMem<AutoDiff<3,T>,20> leg_u(order+2), leg_v(order+2), leg_w(order+2);
      
      // edges
      for (int i = 0; i < 12; i++)
        {
          int p = order_edge[i]; 
          INT<2> e = ET_trait<ET_HEX>::GetEdgeSort (i, vnums);
          AutoDiff<3,T> xi  = sigma[e[1]]-sigma[e[0]];
          AutoDiff<3,T> lam_e = lami[e[0]]+lami[e[1]];  
          Vec<6, AutoDiff<3,T>> symdyadic = SymDyadProd(xi,xi);
          
          IntLegNoBubble::
            EvalMult (p, xi, lam_e, SBLambda ([&](int i, auto val)
                                              {
                                                shape[ii++] = T_REGGE_Shape<3,T>(val*symdyadic);
                                              }));
        }
      
      
      for (int i = 0; i<6; i++)
        {
          int p = order_facet[i][0];
          
          AutoDiff<3,T> lam_f(0);
          for (int j = 0; j < 4; j++)
            lam_f += lami[faces[i][j]];
          
          INT<4> f = ET_trait<ET_HEX>::GetFaceSort (i, vnums);	  
          AutoDiff<3,T> xi  = sigma[f[0]] - sigma[f[1]]; 
          AutoDiff<3,T> eta = sigma[f[0]] - sigma[f[3]];
          // auto nv = GetGradient(lam_f);
          
          LegendrePolynomial (p, eta, leg_u);
          LegendrePolynomial (p, xi, leg_v);
          
          Vec<6, AutoDiff<3,T>> symdyadic = lam_f*SymDyadProd(GetGradient(eta),GetGradient(xi));
          for (int j = 0; j <= p; j++)
            for (int k = 0; k <= p; k++)
              shape[ii++] = T_REGGE_Shape<3,T>(leg_u[j]*leg_v[k]*symdyadic);

          symdyadic = lam_f*(1-eta*eta)*SymDyadProd(GetGradient(xi),GetGradient(xi));
          for (int j = 0; j < p; j++)
            for (int k = 0; k <= p; k++)
              shape[ii++] = T_REGGE_Shape<3,T>(leg_u[j]*leg_v[k]*symdyadic);

          symdyadic = lam_f*(1-xi*xi)*SymDyadProd(GetGradient(eta),GetGradient(eta));
          for (int k = 0; k < p; k++)
            for (int j = 0; j <= p; j++)
              shape[ii++] = T_REGGE_Shape<3,T>(leg_u[j]*leg_v[k]*symdyadic);
          
        }

      int p = order_inner[0];
      if (p > 0)
        {
          AutoDiff<3,T> xi  = sigma[0] - sigma[1];
          AutoDiff<3,T> eta = sigma[0] - sigma[3];
          AutoDiff<3,T> nv = sigma[0] - sigma[4];
          
          LegendrePolynomial (p, xi,  leg_u);
          LegendrePolynomial (p, eta, leg_v);
          LegendrePolynomial (p, nv,  leg_w);

          Vec<6, AutoDiff<3,T>> symdyadic1 = lz[0]*lz[1]*SymDyadProd(GetGradient(eta),GetGradient(xi));
          Vec<6, AutoDiff<3,T>> symdyadic2 = lx[0]*lx[1]*SymDyadProd(GetGradient(nv),GetGradient(eta));
          Vec<6, AutoDiff<3,T>> symdyadic3 = ly[0]*ly[1]*SymDyadProd(GetGradient(xi),GetGradient(nv));
          for (int i = 0; i <= p; i++)
            for (int j = 0; j <= p; j++)
              for (int k = 0; k < p; k++)
                {
                  shape[ii++] = T_REGGE_Shape<3,T>(leg_u[i]*leg_v[j]*leg_w[k]*symdyadic1);
                  shape[ii++] = T_REGGE_Shape<3,T>(leg_v[i]*leg_w[j]*leg_u[k]*symdyadic2);
                  shape[ii++] = T_REGGE_Shape<3,T>(leg_w[i]*leg_u[j]*leg_v[k]*symdyadic3);
                }

          symdyadic1 = ly[0]*ly[1]*lz[0]*lz[1]*SymDyadProd(GetGradient(xi),GetGradient(xi));
          symdyadic2 = lz[0]*lz[1]*lx[0]*lx[1]*SymDyadProd(GetGradient(eta),GetGradient(eta));
          symdyadic3 = lx[0]*lx[1]*ly[0]*ly[1]*SymDyadProd(GetGradient(nv),GetGradient(nv));

          for (int i = 0; i <= p; i++)
            for (int j = 0; j < p; j++)
              for (int k = 0; k < p; k++)
                {
                  shape[ii++] = T_REGGE_Shape<3,T>(leg_u[i]*leg_v[j]*leg_w[k]*symdyadic1);
                  shape[ii++] = T_REGGE_Shape<3,T>(leg_v[i]*leg_w[j]*leg_u[k]*symdyadic2);
                  shape[ii++] = T_REGGE_Shape<3,T>(leg_w[i]*leg_u[j]*leg_v[k]*symdyadic3);
                }
        }
    }
    
    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hcurlcurlfe calcdualshape2 not implementend for element type ET_HEX");
    }

  };




  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_SEGM>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_TRIG>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_QUAD>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_TET>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_PRISM>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_HEX>;
}

#endif

  
