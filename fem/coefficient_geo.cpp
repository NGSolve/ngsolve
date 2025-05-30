/*********************************************************************/
/* File:   coefficient_geo.cpp                                       */
/* Author: Joachim Schoeberl                                         */
/* Date:   Fem 2020                                                  */
/*********************************************************************/

/* 
   Geometric coefficient functions
*/

// #include <fem.hpp>
#include <coefficient.hpp>
#include "scalarfe.hpp"
#include "tpintrule.hpp"
#include <core/register_archive.hpp>

namespace ngfem

{

  
  
  // ////////////////////////// Coordinate CF ////////////////////////

  class CoordCoefficientFunction
    : public T_CoefficientFunction<CoordCoefficientFunction, CoefficientFunctionNoDerivative>
  {
    int dir;
    typedef T_CoefficientFunction<CoordCoefficientFunction, CoefficientFunctionNoDerivative> BASE;
  public:
    CoordCoefficientFunction() = default;
    CoordCoefficientFunction (int adir) : BASE(1, false), dir(adir) { SetVariable(true); }

    void DoArchive(Archive& ar) override
    {
      BASE::DoArchive(ar);
      ar & dir;
    }

    virtual string GetDescription () const override
    {
      string dirname;
      switch (dir)
        {
        case 0: dirname = "x"; break;
        case 1: dirname = "y"; break;
        case 2: dirname = "z"; break;
        default: dirname = ToLiteral(dir);
        }
      return string("coordinate ")+dirname;
    }

    using BASE::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      if (!ip.IsComplex())
        return ip.GetPoint()(dir);
      else
        return ip.GetPointComplex()(dir).real();
    }
    /*
      virtual void Evaluate(const BaseMappedIntegrationRule & ir,
      FlatMatrix<> result) const
      {
      if (!ir.IsComplex())
      result.Col(0) = ir.GetPoints().Col(dir);
      else
      {
      auto pnts = ir.GetPointsComplex().Col(dir);
      for (auto i : Range(ir.Size()))
      result(i,0) = pnts(i).real();
      }
      }
      virtual void Evaluate(const BaseMappedIntegrationRule & ir,
      FlatMatrix<Complex> result) const override
      {
      result.Col(0) = ir.GetPoints().Col(dir);
      }
    */

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      auto v = Var(index);
      // code.body += v.Assign(CodeExpr(string("mir.GetPoints()(i,")+ToLiteral(dir)+")"));

      // code.Declare(code.res_type, index, Dimensions());
      code.Declare(index, Dimensions(), IsComplex());
      code.body += v.Assign(CodeExpr(string("points(i,")+ToLiteral(dir)+")"), false);
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
    {
      size_t nv = ir.Size();
      __assume (nv > 0);

      if(dir>=ir.DimSpace())
        {
          for (size_t i = 0; i < nv; i++)
            values(0,i) = 0.0;
          return;
        }

      if(!ir.IsComplex())
        {
          auto points = ir.GetPoints();
          for (size_t i = 0; i < nv; i++)
            values(0,i) = points(i, dir);
        }
      else
        {
          auto cpoints = ir.GetPointsComplex();
          __assume (nv > 0);
          for(auto i : Range(nv))
            values(0,i) = cpoints(i,dir).real();
        }
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    { T_Evaluate (ir, values); }

    /*
      virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
      AFlatMatrix<double> values) const
      {
      Evaluate (ir, values);
      }
    */

    using CoefficientFunction::Operator;
    shared_ptr<CoefficientFunction> Operator (const string & name) const override
    {
      if (spacedim == -1)
        throw Exception("cannot differentiate coordinate since we don't know the space dimension, use 'coef.spacedim=dim'");
      if (name != "grad")
        throw Exception ("cannot apply operator "+name+" for coordinate");
      
      Array<shared_ptr<CoefficientFunction>> funcs(spacedim);
      funcs = ZeroCF(Array<int>());
      funcs[dir] = make_shared<ConstantCoefficientFunction> (1);
      return MakeVectorialCoefficientFunction (std::move(funcs));
    }
    
    
    
    shared_ptr<CoefficientFunction>
    Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dirdiff) const override
    {
      // if (var == shape.get())
      if (dynamic_cast<const DiffShapeCF*>(var))
        return MakeComponentCoefficientFunction (dirdiff, dir);
      // return BASE::Diff (var, dirdiff);
      
      if (auto coordcf = dynamic_cast<const CoordCoefficientFunction*>(var))
        if (coordcf->dir == this->dir)
          return dirdiff;
      
      return ZeroCF(Dimensions());
    }
    
    
  };


  shared_ptr<CoefficientFunction> MakeCoordinateCoefficientFunction (int comp)
  {
    return make_shared<CoordCoefficientFunction> (comp);
  }



  

  template <int D>
  class cl_NormalVectorCF : public CoefficientFunctionNoDerivative
  {
    BitArray inverted_faces;
  public:
    cl_NormalVectorCF (optional<BitArray> ainverted_faces = nullopt)
      : CoefficientFunctionNoDerivative(D,false)
    {
      SetDimensions(Array<int>( { D } ));
      if(ainverted_faces.has_value())
        inverted_faces = ainverted_faces.value();
    }
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
      double fac = 1.;
      if (inverted_faces.Size())
        {
          auto ei = ip.GetTransformation().GetElementIndex();
          if(inverted_faces.Test(ei))
            fac = -1.;
        }
      res = fac * static_cast<const DimMappedIntegrationPoint<D>&>(ip).GetNV();
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<> res) const // override 
    {
      const TPMappedIntegrationRule * tpir = dynamic_cast<const TPMappedIntegrationRule *>(&ir);
      if(!tpir)
        {
          double fac = 1.;
          if (inverted_faces.Size())
            {
              auto ei = ir.GetTransformation().GetElementIndex();
              if(inverted_faces.Test(ei))
                fac = -1.;
            }
          if (ir[0].DimSpace() != D)
            throw Exception("illegal dim of normal vector");
          // FlatMatrixFixWidth<D> resD(res);
          for (int i = 0; i < ir.Size(); i++)
            res.Row(i).Range(D) = fac * static_cast<const DimMappedIntegrationPoint<D>&>(ir[i]).GetNV();
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

      double fac = 1.;
      if (inverted_faces.Size())
        {
          auto ei = ir.GetTransformation().GetElementIndex();
          if(inverted_faces.Test(ei))
            fac = -1.;
        }
      for (int i = 0; i < ir.Size(); i++)
	res.Row(i).Range(D) = fac * static_cast<const DimMappedIntegrationPoint<D>&>(ir[i]).GetNV();
    }
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override  {
      if (inverted_faces.Size())
        throw Exception("Not implemented");
      /*
      string miptype;
      if(code.is_simd)
        miptype = "SIMD<DimMappedIntegrationPoint<"+ToLiteral(D)+">>*";
      else
        miptype = "DimMappedIntegrationPoint<"+ToLiteral(D)+">*";
      auto nv_expr = CodeExpr("static_cast<const "+miptype+">(&ip)->GetNV()");
      auto nv = Var("tmp", index);
      code.body += nv.Assign(nv_expr);
        
      code.Declare (index, this->Dimensions(), this->IsComplex());  
      for( int i : Range(D))
        code.body += Var(index,i).Assign(nv(i), false);
      */

      code.Declare (index, this->Dimensions(), this->IsComplex());  
      for( int i : Range(D))
        code.body += Var(index,i).Assign(CodeExpr("normals(i,"+ToLiteral(i)+")"), false);
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override 
    {
      /*
        for (size_t i = 0; i < ir.Size(); i++)
        for (size_t j = 0; j < D; j++)
        values(j,i) = static_cast<const SIMD<DimMappedIntegrationPoint<D>>&>(ir[i]).GetNV()(j).Data();
      */
      double fac = 1.;
      if (inverted_faces.Size())
        {
          auto ei = ir.GetTransformation().GetElementIndex();
          if(inverted_faces.Test(ei))
            fac = -1.;
        }
      values.AddSize(D, ir.Size()) = fac * Trans(ir.GetNormals());
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
      // if (var == shape.get())
      if (dynamic_cast<const DiffShapeCF*>(var))
        return -TransposeCF(dir->Operator("Gradboundary")) * const_cast<cl_NormalVectorCF*>(this)->shared_from_this();
      return CoefficientFunctionNoDerivative::Diff(var, dir);
    }
    
  };

  shared_ptr<CoefficientFunction> NormalVectorCF (int dim,
                                                  optional<BitArray> inverted_faces)
  { 
    switch(dim)
      { 
      case 1:
        return make_shared<cl_NormalVectorCF<1>>(inverted_faces);
      case 2:
        return make_shared<cl_NormalVectorCF<2>>(inverted_faces);
      case 3:
        return make_shared<cl_NormalVectorCF<3>>(inverted_faces);
      case 4:
        return make_shared<cl_NormalVectorCF<4>>(inverted_faces);
      case 5:
        return make_shared<cl_NormalVectorCF<5>>(inverted_faces);
      case 6:
        return make_shared<cl_NormalVectorCF<6>>(inverted_faces);
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

      if (consistent)
        throw Exception("consistent tangent does not support Compile(True) yet");
        
      string miptype;
      if(code.is_simd)
        miptype = "SIMD<DimMappedIntegrationPoint<"+ToLiteral(D)+">>*";
      else
        miptype = "DimMappedIntegrationPoint<"+ToLiteral(D)+">*";
      auto tv_expr = CodeExpr("static_cast<const "+miptype+">(&ip)->GetTV()");
      auto tv = Var("tmp", index);
      code.body += tv.Assign(tv_expr);

      code.Declare (index, this->Dimensions(), this->IsComplex()); 
      for( int i : Range(D))
        code.body += Var(index,i).Assign(tv(i), false);
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
      // if (var == shape.get())
      if (dynamic_cast<const DiffShapeCF*>(var))                
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
      	res.Row(i).Range(DIMS*DIMR) = static_cast<const MappedIntegrationPoint<DIMS,DIMR>&>(ir[i]).GetJacobian().AsVector();
    }
    
    /*virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override 
      {
      values.AddSize(D*D, ir.Size()) = Trans(ir.GetJacobian());
      }*/

    virtual shared_ptr<CoefficientFunction>
    Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override
    {
      // if (var == shape.get())
      if (dynamic_cast<const DiffShapeCF*>(var))        
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

      if (!eltrans.IsCurvedElement())
	{
	  res = 0;
	  return;
	}
      
      double eps = 1e-4;

      Mat<D,D-1> dshape;
      
      if (bmip.DimSpace() != bmip.DimElement())
	{
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
      else
	{
          LocalHeapMem<10000> lh("Weingarten-lh-nosimd");

	  ELEMENT_TYPE et = eltrans.GetElementType();
          int fnr = bmip.IP().FacetNr();
          auto f2eltrafo = Facet2ElementTrafo(et);

	  for (int j = 0; j < D-1; j++)   // d / dxj
	    {
	      HeapReset hr(lh);
	      IntegrationPoint ipts[4];
              ipts[0] = ip;
	      ipts[0](j) -= eps;
	      ipts[1] = ip;
              ipts[1](j) += eps;
	      ipts[2] = ip;
              ipts[2](j) -= 2*eps;
	      ipts[3] = ip;
	      ipts[3](j) += 2*eps;

              IntegrationPoint ipts_vol[4];
              f2eltrafo(fnr,ipts[0],ipts_vol[0]);
              f2eltrafo(fnr,ipts[1],ipts_vol[1]);
              f2eltrafo(fnr,ipts[2],ipts_vol[2]);
              f2eltrafo(fnr,ipts[3],ipts_vol[3]);

              IntegrationRule ir(4,ipts_vol);
              MappedIntegrationRule<D,D> mir(ir, eltrans, lh);

              dshape.Col(j) = (1.0/(12.0*eps)) * (8.0*mir[1].GetNV()-8.0*mir[0].GetNV()-mir[3].GetNV()+mir[2].GetNV());
	    }

          Mat<D,D-1> F = bmip.GetJacobian()*f2eltrafo.GetJacobian(fnr,lh);
          Mat<D-1,D-1> FtF = Trans(F)*F;
          Mat<D-1,D> Finv = Inv(FtF)*Trans(F);
          res = (dshape*Finv).AsVector();
        }
    }

    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values) const override 
    {
      LocalHeapMem<10000> lh("Weingarten-lh");

      double eps = 1e-4;

      Mat<D,D-1,SIMD<double>> dshape;
      Mat<D,D,SIMD<double>> phys_dshape;

      if (bmir[0].DimSpace() != bmir[0].DimElement())
	{
          auto & mir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (bmir);
          auto & ir = mir.IR();
      
	  for (size_t i = 0; i < mir.Size(); i++)
	    {
	      const SIMD<IntegrationPoint> & ip = ir[i];
	      const ElementTransformation & eltrans = mir[i].GetTransformation();
	      
	      if (!eltrans.IsCurvedElement())
		{
		  values.Col(i).Range(D*D) = 0;
		  continue;
		}
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
      else
	{
          auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D>&> (bmir);
          auto & ir = mir.IR();

	  for (size_t i = 0; i < mir.Size(); i++)
	    {
	      const SIMD<IntegrationPoint> & ip = ir[i];
	      const ElementTransformation & eltrans = mir[i].GetTransformation();
	      
	      if (!eltrans.IsCurvedElement())
		{
		  values.Col(i).Range(D*D) = 0;
		  continue;
		}

	      ELEMENT_TYPE et = eltrans.GetElementType();
              int fnr = ip.FacetNr();
              auto f2eltrafo = Facet2ElementTrafo(et);

              for (int j = 0; j < D-1; j++)   // d / d t_j
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
                  
		  SIMD_IntegrationRule ir(4,ipts);
                  const SIMD_IntegrationRule & f2elir =f2eltrafo(fnr,ir,lh);
                  SIMD_MappedIntegrationRule<D,D> mir(f2elir, eltrans, lh);
                  
                  dshape.Col(j) = (1.0/(12.0*eps)) * (mir.GetNormals().Row(2) - mir.GetNormals().Row(3) - 8.0*mir.GetNormals().Row(0) + 8.0*mir.GetNormals().Row(1));
		}
              
              Mat<D,D-1,SIMD<double>> F = mir[i].GetJacobian()*f2eltrafo.GetJacobian(fnr,lh);
              Mat<D-1,D-1,SIMD<double>> FtF = Trans(F)*F;
              Mat<D-1,D,SIMD<double>> Finv = Inv(FtF)*Trans(F);
              phys_dshape = dshape*Finv;
              
	      for (size_t l = 0; l < D*D; l++)
		values(l, i) = phys_dshape(l);
	    }
	}
      
    }

    virtual shared_ptr<CoefficientFunction>
    Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override
    {
      // if (var == shape.get())
      if (dynamic_cast<const DiffShapeCF*>(var))                
        {
          int dim = dir->Dimension();
          auto n = NormalVectorCF(dim) -> Reshape( Array<int> ( { dim, 1 } ) );
          auto Pn = n * TransposeCF(n);

          auto WG = const_cast<cl_WeingartenCF*>(this)->shared_from_this();
          auto dX = dir->Operator("Gradboundary");
          Array<shared_ptr<CoefficientFunction>> cflist(1);
          cflist[0] = TransposeCF(dir->Operator("hesseboundary"))*n;
          auto Hn = MakeVectorialCoefficientFunction(std::move(cflist))->Reshape(dim, dim);
          
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

    using CoefficientFunctionNoDerivative::Evaluate;
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
          // auto pnt = ip.IP().Point();
          ELEMENT_TYPE et = trafo.GetElementType();
          //int iavnums[] = { 0, 1, 2, 3 };
          //FlatArray<int> vnums(4, &iavnums[0]);
          //trafo.GetSort(vnums);

	  //cout << "vnr = " << vnr << endl << "pnt = " << pnt << endl << endl;
	  //cout << "F = " << F << endl;
	  
          Vec<2> tv_ref[3] = { Vec<2>(1,0), Vec<2>(0,1), Vec<2>(-1,1) };
          Vec<2> tv_ref_v[8];

          switch(et)
            {
            case ET_TRIG:
              tv_ref_v[0] = -tv_ref[0];
              tv_ref_v[1] = tv_ref[2];
              tv_ref_v[2] = -tv_ref[2];
              tv_ref_v[3] = -tv_ref[1];
              tv_ref_v[4] = tv_ref[1];
              tv_ref_v[5] = tv_ref[0];
              break;
            case ET_QUAD:
              tv_ref_v[0] = tv_ref[1];
              tv_ref_v[1] = tv_ref[0];
              tv_ref_v[2] = -tv_ref[0];
              tv_ref_v[3] = tv_ref[1];
              tv_ref_v[4] = -tv_ref[1];
              tv_ref_v[5] = -tv_ref[0];
              tv_ref_v[6] = tv_ref[0];
              tv_ref_v[7] = -tv_ref[1];
              break;
            default:
              throw Exception("VertexTangentialVectorsCF does not support"+ToString(trafo.GetElementType()));
              break;
            }

	  FlatMatrix<double> phys_tv(D,2,&res[0]);
	  Vec<D> tmp = F*tv_ref_v[2*vnr+0];
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
  class cl_EdgeFaceTangentialVectorsCF : public CoefficientFunctionNoDerivative
  {
  public:
    cl_EdgeFaceTangentialVectorsCF () : CoefficientFunctionNoDerivative(2*D,false)
    {
      SetDimensions(Array<int> ( { D, 2 } ));
    }

    using CoefficientFunctionNoDerivative::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      return 0;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> res) const override
    {
      if (ip.DimSpace() != D)
        throw Exception("illegal dim of EdgeFaceTangentialVector");

      //cout << "in EdgeFaceTangentialVectorsCF, VB = " << ip.IP().VB() << endl;
      // assume tets !!!
      
      if (ip.IP().VB() == BBND)
        {
          Mat<3,3> F = ip.GetJacobian();
          int edgenr = ip.IP().FacetNr();

          /*
            edge -> vertex is 
            static const int tet_edges[6][2] =
            { { 3, 0 },
            { 3, 1 },
            { 3, 2 }, 
            { 0, 1 }, 
            { 0, 2 },
            { 1, 2 }};

            face -> vertex is 
            { { 3, 1, 2, -1 },
            { 3, 2, 0, -1 },
            { 3, 0, 1, -1 },
            { 0, 2, 1, -1 } }; // all faces point into interior!

            from this we get edge -> face
          */
          static int edge2face[6][2] =
            { { 1, 2 },
              { 0, 2 },
              { 0, 1 },
              { 2, 3 },
              { 1, 3 },
              { 0, 3 }
            };

          FlatVector<Vec<3>> normals = ElementTopology::GetNormals<3>(ET_TET);

          auto & trafo = ip.GetTransformation();
          ELEMENT_TYPE et = trafo.GetElementType();
          auto edge = ElementTopology::GetEdges(et)[edgenr];
                    
          auto [v1x, v1y, v1z] = ElementTopology::GetVertices(ET_TET)[edge[0]];
          Vec<3> v1(v1x, v1y, v1z);
          auto [v2x, v2y, v2z] = ElementTopology::GetVertices(ET_TET)[edge[1]];
          Vec<3> v2(v2x, v2y, v2z);
          
          Vec<3> nref1 = normals[edge2face[edgenr][0]];
          Vec<3> nref2 = normals[edge2face[edgenr][1]];
          Vec<3> tref = v2-v1;

          // fix orientation of tangential vector (compare TangentialVector with consistent = true)
          int iavnums[] = { 0, 1, 2, 3 };
          FlatArray<int> vnums(4, &iavnums[0]);
          trafo.GetSort(vnums);
          int invnums[4];
          for (int i = 0; i < 4; i++)
            invnums[iavnums[i]] = i;
          if (invnums[edge[0]] > invnums[edge[1]])
            tref *= -1;
          
          Vec<3> tref1 = Cross(tref, nref1);
          Vec<3> tref2 = -Cross(tref, nref2);

          Vec<3> t1 = F * tref1;
          Vec<3> t2 = F * tref2;
          Vec<3> t  = F * tref;

          t /= L2Norm(t);
          t1 -= InnerProduct(t1,t)*t;
          t2 -= InnerProduct(t2,t)*t;

          t1 /= L2Norm(t1);
          t2 /= L2Norm(t2);
          
	  FlatMatrix<double> phys_tv(D,2,&res[0]);

          // after orientation of t has been fixed the orientation of the t1,t2 are fixed by forcing det(t,t1,t2)>0
          if (InnerProduct(Cross(t,t1),t2) > 0)
            {
              phys_tv.Col(0) = t1;
              phys_tv.Col(1) = t2;
            }
          else
            {
              phys_tv.Col(1) = t1;
              phys_tv.Col(0) = t2;
            }
          
        }
      else
        throw Exception("EdgeFaceTangentialVector only makes sense on edges");
    }
    
  };


  shared_ptr<CoefficientFunction> EdgeFaceTangentialVectorsCF (int dim)
  {
    if (dim == 3)
      return make_shared<cl_EdgeFaceTangentialVectorsCF<3>>();
    throw Exception ("EdgeFaceTangentialVectors available only in 3D");
  }




  
  
  template <int D>
  class cl_EdgeCurvatureCF : public CoefficientFunctionNoDerivative
  {
  public:
    cl_EdgeCurvatureCF () : CoefficientFunctionNoDerivative(D,false) { ; }

    using CoefficientFunctionNoDerivative :: Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      return 0;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationPoint & bmip, FlatVector<> res) const override
    {
      if (bmip.DimSpace() != D)
        throw Exception("illegal dim of EdgeCurvatureCF");


      // (nabla_t t)circ\phi = nabla(t\circ\phi) F^-1 t\circ\phi
      //                     = nabla(t\circ\phi) F^-1 1/|F t_ref| F t_ref
      //                     = 1/|F t_ref|*nabla(t\circ\phi) t_ref
      //                     = 1/|F t_ref|*nabla(t\circ\phi) t_ref


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
          Vec<D> tv_phys = F*tv_ref;
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

          MappedIntegrationPoint<2,D> sipl(ipl, eltrans);
          MappedIntegrationPoint<2,D> sipr(ipr, eltrans);
          MappedIntegrationPoint<2,D> sipll(ipll, eltrans);
          MappedIntegrationPoint<2,D> siprr(iprr, eltrans);

          // Need unit tangent vectors at the stencil points, not directly computed in MappedIntegrationPoint
	  Mat<D,2> Ft = sipr.GetJacobian();
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
       

          res = 1/measure*dshape;
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





  
  static RegisterClassForArchive<CoordCoefficientFunction, CoefficientFunction> regcoocf;






  
}


