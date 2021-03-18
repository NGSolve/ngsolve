/*********************************************************************/
/* File:   coefficient_geo.cpp                                       */
/* Author: Joachim Schoeberl                                         */
/* Date:   Fem 2020                                                  */
/*********************************************************************/

/* 
   Geometric coefficient functions
*/

#include <fem.hpp>

namespace ngfem

{

  template <int D>
  class cl_NormalVectorCF : public CoefficientFunctionNoDerivative
  {
  public:
    cl_NormalVectorCF () : CoefficientFunctionNoDerivative(D,false) { ; }
    // virtual int Dimension() const { return D; }

    virtual string GetDescription() const override
    {
      return "normal vector";
    }
    
      using CoefficientFunctionNoDerivative::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override 
    {
      return 0;
    }
    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> res) const override 
    {
      if (ip.DimSpace() != D)
        throw Exception("illegal dim of normal vector");
      res = static_cast<const DimMappedIntegrationPoint<D>&>(ip).GetNV();
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<> res) const // override 
    {
      const TPMappedIntegrationRule * tpir = dynamic_cast<const TPMappedIntegrationRule *>(&ir);
       if(!tpir)
       {
         if (ir[0].DimSpace() != D)
           throw Exception("illegal dim of normal vector");
         FlatMatrixFixWidth<D> resD(res);
         for (int i = 0; i < ir.Size(); i++)
           resD.Row(i) = static_cast<const DimMappedIntegrationPoint<D>&>(ir[i]).GetNV();
       }
       else
       {
         int facet = tpir->GetFacet();
         auto & mir = *tpir->GetIRs()[facet];
         int dim = mir[0].DimSpace();
         int ii = 0;
         res = 0.0;
         if(facet == 0)
         {
           if(dim == 1)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)
                 res.Row(ii++).Range(0,dim) = static_cast<const DimMappedIntegrationPoint<1>&>(mir[i]).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 2)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(0,dim) = static_cast<const DimMappedIntegrationPoint<2>&>(mir[i]).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 3)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(0,dim) = static_cast<const DimMappedIntegrationPoint<3>&>(mir[i]).GetNV();//res1.Row(i).Range(0,dim);
         }
         else
         {
           if(dim == 1)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)
                 res.Row(ii++).Range(D-dim,D) = static_cast<const DimMappedIntegrationPoint<1>&>(mir[j]).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 2)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(D-dim,D) = static_cast<const DimMappedIntegrationPoint<2>&>(mir[j]).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 3)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(D-dim,D) = static_cast<const DimMappedIntegrationPoint<3>&>(mir[j]).GetNV();//res1.Row(i).Range(0,dim);
         }
      }
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> res) const override 
    {
      if (ir[0].DimSpace() != D)
	throw Exception("illegal dim of normal vector");
      for (int i = 0; i < ir.Size(); i++)
	res.Row(i).AddSize(D) = static_cast<const DimMappedIntegrationPoint<D>&>(ir[i]).GetNV();
    }
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override  {
        string miptype;
        if(code.is_simd)
          miptype = "SIMD<DimMappedIntegrationPoint<"+ToLiteral(D)+">>*";
        else
          miptype = "DimMappedIntegrationPoint<"+ToLiteral(D)+">*";
        auto nv_expr = CodeExpr("static_cast<const "+miptype+">(&ip)->GetNV()");
        auto nv = Var("tmp", index);
        code.body += nv.Assign(nv_expr);
        for( int i : Range(D))
          code.body += Var(index,i).Assign(nv(i));
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override 
    {
      /*
      for (size_t i = 0; i < ir.Size(); i++)
        for (size_t j = 0; j < D; j++)
          values(j,i) = static_cast<const SIMD<DimMappedIntegrationPoint<D>>&>(ir[i]).GetNV()(j).Data();
      */
      values.AddSize(D, ir.Size()) = Trans(ir.GetNormals());
    }

    /*
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const override 
    {
      Evaluate (ir, values);
    }

    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                AFlatMatrix<double> values, AFlatMatrix<double> deriv) const
    {
      Evaluate (ir, values);
      deriv = 0.0;
    }
    
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input,
                                FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    {
      Evaluate (ir, result);
      deriv = 0.0;
    }
    */


    using CoefficientFunction::Operator;
    shared_ptr<CoefficientFunction> Operator (const string & name) const override
    {
      if (name == "grad" || name == "Grad")
        return WeingartenCF (D);
      throw Exception("Normalvector cannot build operator " + name);
    }
    

    virtual shared_ptr<CoefficientFunction>
    Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == shape.get())
        return -TransposeCF(dir->Operator("Gradboundary")) * const_cast<cl_NormalVectorCF*>(this)->shared_from_this();
      return CoefficientFunctionNoDerivative::Diff(var, dir);
    }
    
  };

  shared_ptr<CoefficientFunction> NormalVectorCF (int dim)
  { 
    switch(dim)
      { 
      case 1:
        return make_shared<cl_NormalVectorCF<1>>();
      case 2:
        return make_shared<cl_NormalVectorCF<2>>();
      case 3:
        return make_shared<cl_NormalVectorCF<3>>();
      case 4:
        return make_shared<cl_NormalVectorCF<4>>();
      case 5:
        return make_shared<cl_NormalVectorCF<5>>();
      case 6:
        return make_shared<cl_NormalVectorCF<6>>();
      default:
        throw Exception (string("Normal-vector not implemented for dimension")
                         +ToString(dim));
      }
  }


  template <int D>
  class cl_TangentialVectorCF : public CoefficientFunctionNoDerivative
  {
    bool consistent;
  public:
    cl_TangentialVectorCF (bool aconsistent)
      : CoefficientFunctionNoDerivative(D,false), consistent(aconsistent) { ; }
    // virtual int Dimension() const { return D; }

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      return 0;
    }
    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> res) const override
    {
      if (ip.DimSpace() != D)
        throw Exception("illegal dim of tangential vector");
      
      res = static_cast<const DimMappedIntegrationPoint<D>&>(ip).GetTV();

      if (consistent)
        {
          auto & trafo = ip.GetTransformation();
          int fnr = ip.IP().FacetNr();
          ELEMENT_TYPE et = trafo.GetElementType();
          auto e = ElementTopology::GetEdges(et)[fnr];
          int iavnums[] = { 0, 1, 2, 3 };
          FlatArray<int> vnums(4, &iavnums[0]);
          trafo.GetSort(vnums);
          int invnums[4];
          for (int i = 0; i < 4; i++)
            invnums[iavnums[i]] = i;
          if (invnums[e[0]] > invnums[e[1]])
            res *= -1;
        }
    }
    
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
        string miptype;
        if(code.is_simd)
          miptype = "SIMD<DimMappedIntegrationPoint<"+ToLiteral(D)+">>*";
        else
          miptype = "DimMappedIntegrationPoint<"+ToLiteral(D)+">*";
        auto tv_expr = CodeExpr("static_cast<const "+miptype+">(&ip)->GetTV()");
        auto tv = Var("tmp", index);
        code.body += tv.Assign(tv_expr);
        for( int i : Range(D))
          code.body += Var(index,i).Assign(tv(i));
    }

      using CoefficientFunctionNoDerivative::Evaluate;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override
    {
      if (consistent)
        throw ExceptionNOSIMD("consistent tangent doest not support SIMD");
      
      for (size_t i = 0; i < ir.Size(); i++)
        for (size_t j = 0; j < D; j++)
          values(j,i) = static_cast<const SIMD<DimMappedIntegrationPoint<D>>&>(ir[i]).GetTV()(j).Data();
    }

    virtual shared_ptr<CoefficientFunction>
    Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override
    {
      //d/dt tang|t=0 = dX*tang - ((dX*tang)*tang)*tang
      if (var == shape.get())
        return dir->Operator("Gradboundary") * const_cast<cl_TangentialVectorCF*>(this)->shared_from_this() - InnerProduct(dir->Operator("Gradboundary")*const_cast<cl_TangentialVectorCF*>(this)->shared_from_this(),const_cast<cl_TangentialVectorCF*>(this)->shared_from_this())*const_cast<cl_TangentialVectorCF*>(this)->shared_from_this();    
      return CoefficientFunctionNoDerivative::Diff(var, dir);
    }
  };

  
  shared_ptr<CoefficientFunction> TangentialVectorCF (int dim, bool consistent)
  {
    switch(dim)
      {
      case 1:
        return make_shared<cl_TangentialVectorCF<1>>(consistent);
      case 2:
        return make_shared<cl_TangentialVectorCF<2>>(consistent);
      default:
        return make_shared<cl_TangentialVectorCF<3>>(consistent);
      }
  }


  template <int DIMS, int DIMR>
  class cl_JacobianMatrixCF : public CoefficientFunctionNoDerivative
  {
  public:
    cl_JacobianMatrixCF () : CoefficientFunctionNoDerivative(DIMR*DIMS,false)
    {
      SetDimensions(Array<int>({DIMR,DIMS}));
    }

    using CoefficientFunctionNoDerivative::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override 
    {
      return 0;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> res) const override 
    {
      if (ip.DimSpace() != DIMR)
        throw Exception("illegal dim!");
      res = static_cast<const MappedIntegrationPoint<DIMS,DIMR>&>(ip).GetJacobian().AsVector();
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> res) const override 
    {
      if (ir[0].DimSpace() != DIMR)
      	throw Exception("illegal dim!");
      for (int i = 0; i < ir.Size(); i++)
      	res.Row(i).AddSize(DIMS*DIMR) = static_cast<const MappedIntegrationPoint<DIMS,DIMR>&>(ir[i]).GetJacobian().AsVector();
    }
    
    /*virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override 
    {
      values.AddSize(D*D, ir.Size()) = Trans(ir.GetJacobian());
      }*/

    virtual shared_ptr<CoefficientFunction>
    Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == shape.get())
        throw Exception("Shape derivative not implemented yet for JacobianMatrixCF");
      return CoefficientFunctionNoDerivative::Diff(var, dir);
    }
  };

  shared_ptr<CoefficientFunction> JacobianMatrixCF (int dims,int dimr)
  {
    switch(dimr)
      {
      case 1:
        return make_shared<cl_JacobianMatrixCF<1,1>>();
      case 2:
        switch(dims)
          {
          case 1:
            return make_shared<cl_JacobianMatrixCF<1,2>>();
          default:
            return make_shared<cl_JacobianMatrixCF<2,2>>();
          }
      default:
        switch(dims)
          {
          case 1:
            return make_shared<cl_JacobianMatrixCF<1,3>>();
          case 2:
            return make_shared<cl_JacobianMatrixCF<2,3>>();
          default:
            return make_shared<cl_JacobianMatrixCF<3,3>>();
          }
      }
  }




  
  template <int D>
  class cl_WeingartenCF : public CoefficientFunctionNoDerivative
  {
  public:
    cl_WeingartenCF () : CoefficientFunctionNoDerivative(D*D,false)
    {
      SetDimensions(Array<int> ( { D, D } ));
    }

    void DoArchive(Archive& ar) override
    {
      CoefficientFunctionNoDerivative::DoArchive(ar);
    }
  
    using CoefficientFunctionNoDerivative::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      return 0;
    }

    void Evaluate (const BaseMappedIntegrationPoint & bmip, FlatVector<> res) const override
    {
      if ( (bmip.DimSpace() != D) || D == 1)
        throw Exception("illegal dim of Weingarten tensor");

      const IntegrationPoint& ip = bmip.IP();
      const ElementTransformation & eltrans = bmip.GetTransformation();

      double eps = 1e-4;

      Mat<D,D-1> dshape;
      
      for (int j = 0; j < D-1; j++)   // d / dxj
        {
          IntegrationPoint ipl(ip);
          ipl(j) -= eps;
          IntegrationPoint ipr(ip);
          ipr(j) += eps;
          IntegrationPoint ipll(ip);
          ipll(j) -= 2*eps;
          IntegrationPoint iprr(ip);
          iprr(j) += 2*eps;
        
          MappedIntegrationPoint<D-1,D> sipl(ipl, eltrans);
          MappedIntegrationPoint<D-1,D> sipr(ipr, eltrans);
          MappedIntegrationPoint<D-1,D> sipll(ipll, eltrans);
          MappedIntegrationPoint<D-1,D> siprr(iprr, eltrans);

          dshape.Col(j) = (1.0/(12.0*eps)) * (8.0*sipr.GetNV()-8.0*sipl.GetNV()-siprr.GetNV()+sipll.GetNV());
        }
      
      res = (dshape*static_cast<const MappedIntegrationPoint<D-1,D>&>(bmip).GetJacobianInverse()).AsVector();
    }

    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values) const override 
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (bmir);
      LocalHeapMem<10000> lh("Weingarten-lh");

      auto & ir = mir.IR();
      double eps = 1e-4;

      Mat<D,D-1,SIMD<double>> dshape;
      Mat<D,D,SIMD<double>> phys_dshape;
      
      for (size_t i = 0; i < mir.Size(); i++)
        {
          const SIMD<IntegrationPoint> & ip = ir[i];
          const ElementTransformation & eltrans = mir[i].GetTransformation();

          for (int j = 0; j < D-1; j++)   // d / dxj
            {
              HeapReset hr(lh);
              SIMD<IntegrationPoint> ipts[4];
              ipts[0] = ip;
              ipts[0](j) -= eps;
              ipts[1] = ip;
              ipts[1](j) += eps;              
              ipts[2] = ip;
              ipts[2](j) -= 2*eps;
              ipts[3] = ip;
              ipts[3](j) += 2*eps;

              SIMD_IntegrationRule ir(4, ipts);
              SIMD_MappedIntegrationRule<D-1,D> mirl(ir, eltrans, lh);

              dshape.Col(j) = (1.0/(12.0*eps)) * ( mirl.GetNormals().Row(2) - mirl.GetNormals().Row(3) - 8.0*mirl.GetNormals().Row(0) + 8.0*mirl.GetNormals().Row(1) );
            }

          phys_dshape = dshape*mir[i].GetJacobianInverse();
              
          for (size_t l = 0; l < D*D; l++)
            values(l, i) = phys_dshape(l);
        }
      
    }

    virtual shared_ptr<CoefficientFunction>
    Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == shape.get())
        {
          int dim = dir->Dimension();
          auto n = NormalVectorCF(dim);
          n -> SetDimensions( Array<int> ( { dim, 1 } ) );
          auto Pn = n * TransposeCF(n);

          auto WG = const_cast<cl_WeingartenCF*>(this)->shared_from_this();
          auto dX = dir->Operator("Gradboundary");
          Array<shared_ptr<CoefficientFunction>> cflist(1);
          cflist[0] = TransposeCF(dir->Operator("hesseboundary"))*n;
          auto Hn = MakeVectorialCoefficientFunction(move(cflist));
          Hn->SetDimensions( Array({dim,dim}) );
          
          return -Hn - TransposeCF(dX)*WG + WG*(2*SymmetricCF(Pn*dX)-dX);
        }
      return CoefficientFunctionNoDerivative::Diff(var, dir);
    }

    
  };


  shared_ptr<CoefficientFunction> WeingartenCF (int dim)
  {
    switch(dim)
      {
      case 1:
        throw Exception ("no WeingartenCF in 1D");
        // return make_shared<WeingartenCF<1>>();
      case 2:
        return make_shared<cl_WeingartenCF<2>>();
      default:
        return make_shared<cl_WeingartenCF<3>>();
      }
  }


  template <int D>
  class cl_VertexTangentialVectorsCF : public CoefficientFunctionNoDerivative
  {
  public:
    cl_VertexTangentialVectorsCF () : CoefficientFunctionNoDerivative(2*D,false)
    {
      SetDimensions(Array<int> ( { D, 2 } ));
    }

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      return 0;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> res) const override
    {
      if (ip.DimSpace() != D)
        throw Exception("illegal dim of VertexTangentialVector");

      //cout << "in VertexTangentialVectorsCF, VB = " << ip.IP().VB() << endl;

      if (ip.IP().VB() == BBND)
        {
          auto F = ip.GetJacobian();
          auto & trafo = ip.GetTransformation();
          int vnr = ip.IP().FacetNr();
          auto pnt = ip.IP().Point();
          //ELEMENT_TYPE et = trafo.GetElementType();
          //int iavnums[] = { 0, 1, 2, 3 };
          //FlatArray<int> vnums(4, &iavnums[0]);
          //trafo.GetSort(vnums);

	  //cout << "vnr = " << vnr << endl << "pnt = " << pnt << endl << endl;
	  //cout << "F = " << F << endl;
	  
          Vec<2> tv_ref[] = { Vec<2>(1,0), Vec<2>(0,1), Vec<2>(-1,1) };
          //Vec<2> tv_ref_v0[] = { -tv_ref[0], tv_ref[2] };
          //Vec<2> tv_ref_v1[] = { -tv_ref[2], -tv_ref[1] };
          //Vec<2> tv_ref_v2[] = { tv_ref[1], tv_ref[0] };
          Vec<2> tv_ref_v [] = { -tv_ref[0], tv_ref[2], -tv_ref[2], -tv_ref[1], tv_ref[1], tv_ref[0] };

	  FlatMatrix<double> phys_tv(D,2,&res[0]);
	  Vec<3> tmp = F*tv_ref_v[2*vnr+0];
          phys_tv.Col(0) = 1/L2Norm(tmp)*tmp;
	  tmp = F*tv_ref_v[2*vnr+1];
          phys_tv.Col(1) = 1/L2Norm(tmp)*tmp;

	  //cout << "reft = " << tv_ref_v[2*vnr+0] << " | " << tv_ref_v[2*vnr+1] << endl;
	  //cout << "physt = " << phys_tv << endl;
        }
      else //throw Exception();
        res = 0.0;
      
    }
    
  };


  shared_ptr<CoefficientFunction> VertexTangentialVectorsCF (int dim)
  {
    switch(dim)
      {
      case 1:
        throw Exception ("no VertexTangentialVectors in 1D");
      case 2:
        return make_shared<cl_VertexTangentialVectorsCF<2>>();
      default:
        return make_shared<cl_VertexTangentialVectorsCF<3>>();
      }
  }


  
  template <int D>
  class cl_EdgeCurvatureCF : public CoefficientFunctionNoDerivative
  {
  public:
    cl_EdgeCurvatureCF () : CoefficientFunctionNoDerivative(D,false) { ; }

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      return 0;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationPoint & bmip, FlatVector<> res) const override
    {
      if (bmip.DimSpace() != D)
        throw Exception("illegal dim of EdgeCurvatureCF");


      // (nabla_t t)circ\phi = Ptau\circ\phi*nabla(t\circ\phi) F^-1 t\circ\phi
      //                     = Ptau\circ\phi*nabla(t\circ\phi) F^-1 1/|F t_ref| F t_ref
      //                     = 1/|F t_ref| * Ptau\circ\phi*nabla(t\circ\phi) t_ref
      //                     = 1/|F t_ref| * F(F^TF)^-1F^T * nabla(t\circ\phi) t_ref


      if (bmip.IP().VB() == BND)
        {
	  const IntegrationPoint& ip = bmip.IP();
	  const ElementTransformation & eltrans = bmip.GetTransformation();

	  double eps = 1e-4;
	  
	  Vec<D> dshape;
          auto mip = static_cast<const MappedIntegrationPoint<D-1,D>&>(bmip);
          
          auto F = bmip.GetJacobian();
          auto & trafo = bmip.GetTransformation();
	  ELEMENT_TYPE et = trafo.GetElementType();
          auto e = ElementTopology::GetEdges(et)[bmip.IP().FacetNr()];
	  

	  Vec<2> pnts[3] = { { 1, 0 }, { 0, 1 } , { 0, 0 } };
          Vec<2> tv_ref = pnts[e[1]] - pnts[e[0]];
	  tv_ref /= L2Norm(tv_ref);

          // compute |F t_ref|
          Vec<3> tv_phys = F*tv_ref;
          double measure = L2Norm(tv_phys);

          // compute nabla(t\circ\phi) t_ref numerically
          // This is the directional derivative in direction t_ref

          IntegrationPoint ipl(ip);
          ipl(0) -= tv_ref[0]*eps;
          ipl(1) -= tv_ref[1]*eps;
          IntegrationPoint ipr(ip);
          ipr(0) += tv_ref[0]*eps;
          ipr(1) += tv_ref[1]*eps;
          IntegrationPoint ipll(ip);
          ipll(0) -= 2*tv_ref[0]*eps;
          ipll(1) -= 2*tv_ref[1]*eps;
          IntegrationPoint iprr(ip);
          iprr(0) += 2*tv_ref[0]*eps;
          iprr(1) += 2*tv_ref[1]*eps;

          MappedIntegrationPoint<D-1,D> sipl(ipl, eltrans);
          MappedIntegrationPoint<D-1,D> sipr(ipr, eltrans);
          MappedIntegrationPoint<D-1,D> sipll(ipll, eltrans);
          MappedIntegrationPoint<D-1,D> siprr(iprr, eltrans);

          // Need unit tangent vectors at the stencil points, not directly computed in MappedIntegrationPoint
	  Mat<3,2> Ft = sipr.GetJacobian();
	  Vec<D> tangr = Ft*tv_ref;
	  tangr /= L2Norm(tangr);
	  Ft = sipl.GetJacobian();
	  Vec<D> tangl = Ft*tv_ref;
	  tangl /= L2Norm(tangl);
	  Ft = siprr.GetJacobian();
	  Vec<D> tangrr = Ft*tv_ref;
	  tangrr /= L2Norm(tangrr);
	  Ft = sipll.GetJacobian();
	  Vec<D> tangll = Ft*tv_ref;
	  tangll /= L2Norm(tangll);

          // nabla(t\circ\phi) t_ref
          dshape = (1.0/(12.0*eps)) * (8.0*tangr-8.0*tangl-tangrr+tangll);

          // Ptau\circ\phi =  F(F^TF)^-1F^T = (I - nsurf\circ\phi \otimes nsurf\circ\phi)
          Mat<D> Ptau=0;
	  for (int i=0; i < D; i++)
	    {
	      Ptau(i,i) = 1.0;
	      Ptau.Col(i) -= mip.GetNV()[i]*mip.GetNV();
	    }

          res = 1/measure*Ptau*dshape;
        }
      else //throw Exception();
        res = 0.0;
      
    }
    
  };


  shared_ptr<CoefficientFunction> EdgeCurvatureCF (int dim)
  {
    switch(dim)
      {
      case 1:
        throw Exception ("no EdgeCurvature in 1D");
      case 2:
        return make_shared<cl_EdgeCurvatureCF<2>>();
      default:
        return make_shared<cl_EdgeCurvatureCF<3>>();
      }
  }
}


