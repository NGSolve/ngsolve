

#include <comp.hpp>

#include "Hddf.hpp"
#include "JKMspace.hpp"

using namespace ngcomp;


namespace ngcomp
{

  // some utilities ...

  inline FlatVector<Vec<3>> GetTangents()
  {
    static Vec<3> tet_tangents[] =
      {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
        {-1, 1, 0},
        {-1, 0, 1},
        {0, -1, 1},
      };
    
    return FlatVector<Vec<3>>(6, tet_tangents);
  };

  inline FlatVector<Vec<2>> GetTangentsTrig()
  {
    static Vec<2> trig_tangents[] =
      {
        {1, 0},
        {0, -1},
          {-1, 1}
      };
    
    return FlatVector<Vec<2>>(3, trig_tangents);
  };

  
  inline FlatVector<IVec<2>>
  e2f_oriented()
  {
    static IVec<2> faces[] =
      {
        {1, 2},
        {2, 0},
        {0, 1},
        {3, 2},
        {1, 3},
        {3, 0},
      };
    
    return FlatVector<IVec<2>>(6, faces);
  };
  
  
  inline FlatVector<IVec<2>>
  v2e_oriented()
  {
    static IVec<2> edges[] =
      {
        { 2, 0 },
        { 1, 2 },
        { 0, 1 }
      };
    
    return FlatVector<IVec<2>>(3, edges);
  };
  



  // a x b = Skw(a) b  
  inline Mat<3, 3> Skw(Vec<3> v)
  {
    Mat<3, 3> m  { {  0, -v[2], v[1] },
                   { v[2], 0, -v[0] },
                   { -v[1], v[0], 0} };
    return m;
  }
  
  template <int D>
  inline Mat<D, D> OuterProduct(Vec<D> a, Vec<D> b)
  {
    Mat<D, D> result;
    for (int i = 0; i < D; i++)
      for (int j = 0; j < D; j++)
        result(i, j) = a(i) * b(j);
    return result;
  }

  inline Mat<3, 3> SymOuterProduct(Vec<3> a, Vec<3> b)
  {
    Mat<3, 3> result;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        result(i, j) = 0.5 * (a(i) * b(j) + a(j) * b(i));
    return result;
  }






  
  
  class HDivDivFacetElement : public FiniteElement
  {
  public:
    HDivDivFacetElement(int ndof, int order)
        : FiniteElement(ndof, order) {}
    ~HDivDivFacetElement() = default;

    virtual void CalcShape(const IntegrationPoint &ip, SliceMatrix<double> shapes) const = 0;

    // div within same edge/face entity
    virtual void CalcDivShape(const IntegrationPoint &ip, SliceMatrix<double> curl_f_shapes) const = 0;

    // codim-2 divshape: on edge + face jumps
    virtual void CalcDivShapeJump(const IntegrationPoint &ip, SliceMatrix<double> div_e_shapes) const
    {
      if (ip.VB() != BBND) // VB() is the enum of the integration point = {VOL, BND, BBND, BBBND}
        throw Exception("only BBoundary integration points allowed");

      if (Dim() == 3)
        {
          int enr = ip.FacetNr();
          IVec<2> nfaces = e2f_oriented()[enr];
          
          div_e_shapes = 0;
          
          for (int j = 0; j < 2; j++)
            {
              IntegrationPoint ip_face = ip;
              ip_face.SetFacetNr(nfaces[j], BND);
              
              Matrix<double> shape_f(ndof, 9);
              CalcShape(ip_face, shape_f);
              
              Vec<3> t = GetTangents()[enr];
              Vec<3> n = ElementTopology::GetNormals<3>(ET_TET)[nfaces[j]];
              
              Vec<3> nu = Cross(n, t);
              
              double sign = (j == 0) ? -1 : 1;
              
              double factor= 1.0/sqr(L2Norm(n));
              factor *= sign;  
              for (int k = 0; k < ndof; k++)
                {
                  FlatMatrix<> shapei = shape_f.Row(k).AsMatrix(3,3);
                  Vec<3> sigma_nu = shapei * nu;
                  div_e_shapes.Row(k) += factor * sigma_nu;  
                }
            }
        }
      else
        {
          int vnr = ip.FacetNr();
          IVec<2> nedges = v2e_oriented()[vnr];
          
          div_e_shapes = 0;
          
          for (int j = 0; j < 2; j++)
            {
              IntegrationPoint ip_face = ip;
              ip_face.SetFacetNr(nedges[j], BND);
              
              Matrix<double> shape_f(ndof, 4);
              CalcShape(ip_face, shape_f);
              
              Vec<2> n = ElementTopology::GetNormals<2>(ET_TRIG)[nedges[j]];
              double sign = (j == 0) ? -1 : 1;
              double factor = 1/sqr(L2Norm(n));
              Vec<2> tau(-n(1),n(0));

              for (int k = 0; k < ndof; k++)
                {
                  FlatMatrix<> shapei = shape_f.Row(k).AsMatrix(2,2);
                  Vec<2> shapei_tau = shapei * tau;
                  div_e_shapes.Row(k) += -sign * factor * shapei_tau;
                }
            }
        }
    }
  };

  template <ELEMENT_TYPE ET>
  class HDivDivFacetHOElement : public HDivDivFacetElement,
                                public ET_trait<ET>, public VertexOrientedFE<ET>
  {
  protected:
    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_EDGE;
    using ET_trait<ET>::N_FACE;
    using ET_trait<ET>::N_CELL;
    using ET_trait<ET>::FaceType;
    using ET_trait<ET>::DIM;

    using VertexOrientedFE<ET>::GetVertexOrientedFace;
    using VertexOrientedFE<ET>::GetVertexOrientedEdge;
    
  public:
    HDivDivFacetHOElement(int order) :
      HDivDivFacetElement(ET==ET_TRIG ? 3*(order+1) :
                          (ET==ET_TET ? 4 * 9 : 0) , order) {}

    ~HDivDivFacetHOElement() {}

    ELEMENT_TYPE ElementType() const override
    {
      return ET;
    }

    void CalcShape(const IntegrationPoint &ip, SliceMatrix<double> shapes) const override
    {
      throw Exception("HCCFacet, element type " + ToString(ET) + " not implemented");
    }

    virtual void CalcDivShape(const IntegrationPoint &ip, SliceMatrix<double> div_shapes) const override
    {
      throw Exception("HDDFacet, element type " + ToString(ET) + " FaceDivShape not implemented");
    }

  };


  
/*
template <>
void HDivDivFacetHOElement<ET_TET> :: CalcShape(const IntegrationPoint &ip, SliceMatrix<double> shapes) const 
{
  if (ip.VB() != BND)
    throw Exception("only boundary integration points allowed");
  
  shapes = 0;
  
  int facenr = ip.FacetNr();
  auto vertices = GetVertexOrientedFace(facenr);
  
  double lam[4] = {ip(0), ip(1), ip(2), 1 - ip(0) - ip(1) - ip(2)};
  ArrayMem<double, 50> values(ndof);
  
  int vs = vertices[1];
  int ve = vertices[2];
  
  Vec<3> n = ElementTopology::GetNormals<3>(ET_TET)[facenr];
  
  DubinerBasis::Eval(order, lam[vs], lam[ve], values);
  
  int dofs_per_face = (order + 1) * (order + 2) / 2;
  
  Mat<3, 3> nn = OuterProduct(n, n);
  
  for (int i = 0; i < dofs_per_face; i++)
    shapes.Row(dofs_per_face * facenr + i) = values[i] * nn.AsVector();
}

  template <>
  void HDivDivFacetHOElement<ET_TET> ::
  CalcFaceCurlShapes(const IntegrationPoint &ip, SliceMatrix<double> curl_shapes) const
  {
    if (ip.VB() != BND)
      throw Exception("only boundary integration points allowed");
    
    int fnr = ip.FacetNr();
    
    curl_shapes = 0;
    
    auto vertices = GetVertexOrientedFace(fnr);
    
    AutoDiff<3> adx(ip(0), 0), ady(ip(1), 1), adz(ip(2), 2);
    AutoDiff<3> adlam[4] = {adx, ady, adz, 1 - adx - ady - adz};
    
    Vec<3> n = ElementTopology::GetNormals<3>(ET_TET)[fnr];
    
    int vs = vertices[1];
    int ve = vertices[2];
    
    ArrayMem<AutoDiff<3>, 50> values(ndof);
    DubinerBasis::Eval(order, adlam[vs], adlam[ve], values);
    
    int dofs_per_face = (order + 1) * (order + 2) / 2;
    
    Mat<3, 3> nn = OuterProduct(n, n);
    
    for (int i = 0; i < (order + 1) * (order + 2) / 2; i++)
      {
        curl_shapes.Row(dofs_per_face * fnr + i) =
          -(nn * Skw(GetGradient(values[i]))).AsVector();
      }
  }
*/


  template <>
  void HDivDivFacetHOElement<ET_TRIG> :: CalcShape(const IntegrationPoint &ip, SliceMatrix<double> shapes) const 
  {
    if (ip.VB() != BND)
      throw Exception("only boundary integration points allowed");
    
    shapes = 0;
    
    int enr = ip.FacetNr();
    auto vertices = GetVertexOrientedEdge(enr);
    
    double lam[3] = {ip(0), ip(1), 1 - ip(0) - ip(1)};
    ArrayMem<double, 50> values(order+1);
    
    int vs = vertices[0];
    int ve = vertices[1];
    
    auto n = ElementTopology::GetNormals<2>(ET_TRIG)[enr];
    Vec<2> tau(-n(1),n(0));
    LegendrePolynomial::Eval(order, lam[vs]-lam[ve], values);
    
    Mat<2,2> tautau = OuterProduct(tau, tau);
    
    for (int i = 0; i <= order; i++)
      shapes.Row((order+1) * enr + i) = values[i] * tautau.AsVector();
  }



  template <>
  void HDivDivFacetHOElement<ET_TRIG> ::
  CalcDivShape(const IntegrationPoint &ip, SliceMatrix<double> div_shapes) const
  {
    if (ip.VB() != BND)
      throw Exception("only boundary integration points allowed");

    int enr = ip.FacetNr();
    div_shapes = 0;
    
    auto vertices = GetVertexOrientedEdge(enr);
    
    AutoDiff<2> adx(ip(0), 0), ady(ip(1), 1);
    AutoDiff<2> adlam[3] = {adx, ady, 1 - adx - ady};
    
    Vec<2> n = ElementTopology::GetNormals<2>(ET_TRIG)[enr];
    Vec<2> tau(-n(1),n(0));
    
    int vs = vertices[0];
    int ve = vertices[1];
    
    ArrayMem<AutoDiff<2>, 50> values(order+1);
    LegendrePolynomial::Eval(order, adlam[vs]-adlam[ve], values);
    
    int dofs_per_face = order+1;
    
    for (int i = 0; i <= order; i++)
      {
        Vec<2> grad = GetGradient(values[i]);
        div_shapes.Row(dofs_per_face * enr + i) = InnerProduct(tau,grad)*tau;
      }
  }


  // ******************************** TET element *****************************
  
  int TetNDofFace (int order, bool bubble)
  {
    int ndof = 3 * (order+1)*(order+2)/2;
    if (bubble) ndof -= 3*(order+1);
    return ndof;
  }

  int TetNDofEdge (int order, bool bubble)
  {
    int ndof = order+1;
    if (bubble) ndof -= 2;
    return ndof;
  }
  
  template <>
  HDivDivFacetHOElement<ET_TET> ::  HDivDivFacetHOElement(int aorder)
    : HDivDivFacetElement(4*TetNDofFace(aorder, true) + 6*TetNDofEdge (aorder+1, true), aorder) { }
  
  template <>
  void HDivDivFacetHOElement<ET_TET> :: CalcShape(const IntegrationPoint &ip, SliceMatrix<double> shapes) const 
  {
    shapes = 0;
    LocalHeapMem<1000> lh("hddlh");
    
    if (ip.VB() == BND)
      {
        int facenr = ip.FacetNr();
        
        Facet2ElementTrafo f2e(ET_TET,GetVertexNumbers());
        Mat<3,2> jac = f2e.GetJacobian(facenr, lh);
        Mat<2,3> invjac = Inv(Trans(jac)*jac) * Trans(jac);
    
        IntegrationPoint ip0 = f2e(facenr, IntegrationPoint(Vec<2>(0,0)));
        
        Vec<3> vec3d(ip(0)-ip0(0), ip(1)-ip0(1), ip(2)-ip0(2));
        IntegrationPoint ip2d = invjac * vec3d;
    
        HDivDivFE<ET_TRIG> hddtrig(order);
        hddtrig.SetOrderFacet(0, -1);   // make it a bubble
        hddtrig.SetOrderFacet(1, -1);
        hddtrig.SetOrderFacet(2, -1);
        hddtrig.ComputeNDof();

        int ndoftrig = hddtrig.GetNDof();
        Matrix shapes2d(ndoftrig,4);
          
        FE_ElementTransformation<2,2> trafo2d(ET_TRIG);
        MappedIntegrationPoint<2,2> mip2d(ip2d,trafo2d);
        hddtrig.CalcMappedShape_Matrix(mip2d, shapes2d);
        
        int base = ip.FacetNr() * ndoftrig;
        Mat<3,2> trans = jac;  // 1/Jacobidet ?
        
        for (int i = 0; i < ndoftrig; i++)
          {
            Mat<2,2> shape2d = shapes2d.Row(i).AsMatrix(2,2);
            Mat<3,3> shape3d = trans * shape2d * Trans(trans);
            shapes.Row(base+i) = shape3d.AsVector();
          }
      }
    if (ip.VB() == BBND)
      {
        shapes = 0;
        
        int enr = ip.FacetNr();
        auto vertices = GetVertexOrientedEdge(enr);
        
        double lam[4] = {ip(0), ip(1), ip(3), 1-ip(0)-ip(1)-ip(2) };
        ArrayMem<double, 50> values(order);
    
        int vs = vertices[0];
        int ve = vertices[1];
    
        // auto n = ElementTopology::GetNormals<2>(ET_TRIG)[enr];
        // Vec<2> tau(-n(1),n(0));
        auto p1 = ElementTopology::GetVertices(ET_TET)[vertices[1]];
        auto p0 = ElementTopology::GetVertices(ET_TET)[vertices[0]];
        Vec<3> tau;
        for (int j = 0; j < 3; j++) tau(j) = p1[j]-p0[j];
          
        LegendrePolynomial::EvalMult(order-1, lam[vs]-lam[ve], lam[vs]*lam[ve], values);
    
        Mat<3,3> tautau = OuterProduct(tau, tau);

        int base = 4*TetNDofFace(order, true) + enr*order;
        for (int i = 0; i < order; i++)
          shapes.Row(base+i) = values[i] * tautau.AsVector();
      }
  }



  template <>
  void HDivDivFacetHOElement<ET_TET> :: CalcDivShape(const IntegrationPoint &ip, SliceMatrix<double> divshapes) const 
  {
    divshapes = 0;
    LocalHeapMem<1000> lh("hddlh");
    
    if (ip.VB() == BND)
      {
        int facenr = ip.FacetNr();
        
        Facet2ElementTrafo f2e(ET_TET,GetVertexNumbers());
        Mat<3,2> jac = f2e.GetJacobian(facenr, lh);
        Mat<2,3> invjac = Inv(Trans(jac)*jac) * Trans(jac);
    
        IntegrationPoint ip0 = f2e(facenr, IntegrationPoint(Vec<2>(0,0)));
        
        Vec<3> vec3d(ip(0)-ip0(0), ip(1)-ip0(1), ip(2)-ip0(2));
        IntegrationPoint ip2d = invjac * vec3d;
    
        HDivDivFE<ET_TRIG> hddtrig(order);
        hddtrig.SetOrderFacet(0, -1);   // make it a bubble
        hddtrig.SetOrderFacet(1, -1);
        hddtrig.SetOrderFacet(2, -1);
        hddtrig.ComputeNDof();

        int ndoftrig = hddtrig.GetNDof();
        Matrix divshapes2d(ndoftrig,2);
          
        FE_ElementTransformation<2,2> trafo2d(ET_TRIG);
        MappedIntegrationPoint<2,2> mip2d(ip2d,trafo2d);
        hddtrig.CalcMappedDivShape(mip2d, divshapes2d);
        
        int base = ip.FacetNr() * ndoftrig;
        Mat<3,2> trans = jac;  // 1/Jacobidet ?
        
        for (int i = 0; i < ndoftrig; i++)
          divshapes.Row(base+i) = trans*divshapes2d.Row(i);
      }

    if (ip.VB() == BBND)
      {
        divshapes = 0;
        /*
        int enr = ip.FacetNr();
        auto vertices = GetVertexOrientedEdge(enr);
        
        double lam[4] = {ip(0), ip(1), ip(3), 1-ip(0)-ip(1)-ip(2) };
        ArrayMem<double, 50> values(order);
    
        int vs = vertices[0];
        int ve = vertices[1];
    
        // auto n = ElementTopology::GetNormals<2>(ET_TRIG)[enr];
        // Vec<2> tau(-n(1),n(0));
        auto p1 = ElementTopology::GetVertices(ET_TET)[vertices[1]];
        auto p0 = ElementTopology::GetVertices(ET_TET)[vertices[0]];
        Vec<3> tau;
        for (int j = 0; j < 3; j++) tau(j) = p1[j]-p0[j];
          
        LegendrePolynomial::EvalMult(order-1, lam[vs]-lam[ve], lam[vs]*lam[ve], values);
    
        Mat<3,3> tautau = OuterProduct(tau, tau);

        int base = 4*TetNDofFace(order, true) + enr*order;
        for (int i = 0; i < order; i++)
          shapes.Row((order+1) * enr + i) = values[i] * tautau.AsVector();
        */
      }
  }




  
  // ************************** JKM Element ******************************
  

  class HDivDivFacetJKMElement : public HDivDivFacetElement,
                                 public ET_trait<ET_TET>, public VertexOrientedFE<ET_TET>
  {
  protected:
    using ET_trait<ET_TET>::N_VERTEX;
    using ET_trait<ET_TET>::N_EDGE;
    using ET_trait<ET_TET>::N_FACE;
    using ET_trait<ET_TET>::N_CELL;
    using ET_trait<ET_TET>::FaceType;
    using ET_trait<ET_TET>::DIM;

    using VertexOrientedFE<ET_TET>::GetVertexOrientedFace;
    using VertexOrientedFE<ET_TET>::GetVertexOrientedEdge;
    
  public:
    HDivDivFacetJKMElement() : HDivDivFacetElement(4*9, 1) { }
    // ~HDivDivFacetHOElement() {}

    ELEMENT_TYPE ElementType() const override { return ET_TET; }

    void CalcShape(const IntegrationPoint &ip, SliceMatrix<double> shapes) const override
    {
      LocalHeapMem<1000> lh("hddf");
      if (ip.VB() != BND)
        throw Exception("only boundary integration points allowed");
      
      shapes = 0;
      int facenr = ip.FacetNr();
      
      Facet2ElementTrafo f2e(ET_TET,GetVertexNumbers());
      Mat<3,2> jac = f2e.GetJacobian(facenr, lh);
      // Mat<2,3> invjac = Inv(jac);
      Mat<2,3> invjac = Inv(Trans(jac)*jac) * Trans(jac);
      
      IntegrationPoint ip0 = f2e(facenr, IntegrationPoint(Vec<2>(0,0)));
      
      Vec<3> vec3d(ip(0)-ip0(0), ip(1)-ip0(1), ip(2)-ip0(2));
      IntegrationPoint ip2d = invjac * vec3d;
      
      // hardcoded: jkmtrig inner bubbles
      JKMFE_Triangle jkmtrig(1);
      
      Matrix shapes2d(15,4);
      
      FE_ElementTransformation<2,2> trafo2d(ET_TRIG);
      MappedIntegrationPoint<2,2> mip2d(ip2d,trafo2d);
      jkmtrig.CalcMappedShape_Matrix(mip2d, shapes2d);
      
      
      int base = ip.FacetNr() * 9;
      Mat<3,2> trans = jac;  // 1/Jacobidet ?
      
      for (int i = 0; i < 9; i++)
        {
          Mat<2,2> shape2d = shapes2d.Row(6+i).AsMatrix(2,2);
          Mat<3,3> shape3d = trans * shape2d * Trans(trans);
          shapes.Row(base+i) = shape3d.AsVector();
        }
    }

    virtual void CalcDivShape(const IntegrationPoint &ip, SliceMatrix<double> div_shapes) const override
    {
      throw Exception("HDDFacet, element type  JKMTEt FaceDivShape not implemented");
    }

  };




     
  /*

  template <>
  HDivDivFacetHOElement<ET_TET> ::  HDivDivFacetHOElement(int order) :
    HDivDivFacetElement(order >= 1 ? 4*9 : 0, order) { }

  template <>
  void HDivDivFacetHOElement<ET_TET> :: CalcShape(const IntegrationPoint &ip, SliceMatrix<double> shapes) const 
  {
    if (order < 1) return;

    
    LocalHeapMem<1000> lh("hddf");
    if (ip.VB() != BND)
      throw Exception("only boundary integration points allowed");
    
    shapes = 0;
    int facenr = ip.FacetNr();
    
    Facet2ElementTrafo f2e(ET_TET,GetVertexNumbers());
    Mat<3,2> jac = f2e.GetJacobian(facenr, lh);
    // Mat<2,3> invjac = Inv(jac);
    Mat<2,3> invjac = Inv(Trans(jac)*jac) * Trans(jac);
    
    IntegrationPoint ip0 = f2e(facenr, IntegrationPoint(Vec<2>(0,0)));

    Vec<3> vec3d(ip(0)-ip0(0), ip(1)-ip0(1), ip(2)-ip0(2));
    IntegrationPoint ip2d = invjac * vec3d;
    
    // hardcoded: jkmtrig inner bubbles
    JKMFE_Triangle jkmtrig(1);

    Matrix shapes2d(15,4);

    FE_ElementTransformation<2,2> trafo2d(ET_TRIG);
    MappedIntegrationPoint<2,2> mip2d(ip2d,trafo2d);
    jkmtrig.CalcMappedShape_Matrix(mip2d, shapes2d);

    
    int base = ip.FacetNr() * 9;
    Mat<3,2> trans = jac;  // 1/Jacobidet ?
    
    for (int i = 0; i < 9; i++)
      {
        Mat<2,2> shape2d = shapes2d.Row(6+i).AsMatrix(2,2);
        Mat<3,3> shape3d = trans * shape2d * Trans(trans);
        shapes.Row(base+i) = shape3d.AsVector();
      }
  }
  */



  // *************************** DiffOps *********************************
  
  

  /*
  Implementation of the Differential Operators
  */
  template <int D>
  class DiffOpIdHDDFacet : public DiffOp<DiffOpIdHDDFacet<D>>
  {
  public:
    enum
    {
      DIM = 1,
      DIM_SPACE = D,
      DIM_ELEMENT = D,
      DIM_DMAT = D*D,
      DIFFORDER = 0
    };

    static auto &Cast(const FiniteElement &fel)
    {
      return static_cast<const HDivDivFacetElement &>(fel);
    }
    static Array<int> GetDimensions() { return Array<int>({D, D}); }

    template <typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement &fel,
                               const MIP &mip,
                               MAT &&mat, LocalHeap &lh)
    {
      Mat<D, D> jacobian = mip.GetJacobian();
      Cast(fel).CalcShape(mip.IP(), Trans(mat));

      // Vec<D> n = ElementTopology::GetNormals<D>(mip.GetTransformation().GetElementType())[mip.IP().FacetNr()];
      // double factor = 1 / L2Norm(Cof(jacobian) * n);
      double factor = 1.0 / sqr (Det(jacobian));

      for (int i = 0; i < fel.GetNDof(); i++)
        {
          FlatMatrix<> shapei = mat.Col(i).AsMatrix(D, D);
          Mat<D,D> trafo_shapei = factor * jacobian * Mat<D,D>(shapei) * Trans(jacobian);
          shapei = trafo_shapei;
        }
    }
  };

  class DiffOpIdHDDFacetBND : public DiffOp<DiffOpIdHDDFacetBND>
  {
  public:
    enum
    {
      DIM = 1,
      DIM_SPACE = 3,
      DIM_ELEMENT = 2,
      DIM_DMAT = 9,
      DIFFORDER = 0,
    };

    static auto &Cast(const FiniteElement &fel)
    {
      return static_cast<const HDivDivFacetElement &>(fel);
    }
    static Array<int> GetDimensions() { return Array<int>({3, 3}); }

    template <typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement &fel,
                               const MIP &mip,
                               MAT &&mat, LocalHeap &lh)
    {
      Mat<3, 3> jacobian = mip.GetJacobian();
      Cast(fel).CalcShape(mip.IP(), Trans(mat) /* , jacobian*/);
      for (int i = 0; i < fel.GetNDof(); i++)
      {
        FlatMatrix<> shapei(3, 3, &mat(0, i));                                                          // &mat(i, 0)
        Mat<3, 3> trafo_shapei = Cof(mip.GetJacobian()) * Mat<3, 3>(shapei) * Trans(mip.GetJacobian()); // Mat<3, 3> trafo_shapei = (1 / mip.GetJacobiDet()) * mip.GetJacobianInverse() * shapei * mip.GetJacobian(); // check if this is correct
        shapei = trafo_shapei;
      }
    }
  };


  template <int DIM> class DiffOpDivHDDFacet;

  
  template<>
  class DiffOpDivHDDFacet<2> : public DiffOp<DiffOpDivHDDFacet<2>>
  {
  public:
    enum
    {
      DIM = 1,
      DIM_SPACE = 2,
      DIM_ELEMENT = 2,
      DIM_DMAT = 2,
      DIFFORDER = 1,
    };

    static auto &Cast(const FiniteElement &fel)
    {
      return static_cast<const HDivDivFacetElement &>(fel);
    }
    static Array<int> GetDimensions() { return Array<int>({2,}); }

    template <typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement &fel,
                               const MIP &mip,
                               MAT &&mat, LocalHeap &lh)
    {

      Mat<2,2> jacobian = mip.GetJacobian();

      double factor = 1;

      switch (mip.IP().VB())
      {
      case BND:
      {
        Cast(fel).CalcDivShape(mip.IP(), Trans(mat));

        Vec<2> t = GetTangentsTrig()[mip.IP().FacetNr()];
        factor = 1 / L2Norm(jacobian * t);
        break;
      }
      case BBND:
      {
        Cast(fel).CalcDivShapeJump(mip.IP(), Trans(mat));
        break;
      }
      default:
      {
        throw Exception("only boundary integration points allowed");
      }
      }

      for (int i = 0; i < fel.GetNDof(); i++)
      {
        // FlatMatrix<> shapei(3, 3, &mat(0, i));
        FlatVector<> shapei = mat.Col(i);
        // Vec<2> trafo_shapei = factor*Cof(jacobian)*Vec<2>(shapei);
        Vec<2> trafo_shapei = factor*jacobian*Vec<2>(shapei);
        shapei = trafo_shapei;
      }
    }
  };



  template<>
  class DiffOpDivHDDFacet<3> : public DiffOp<DiffOpDivHDDFacet<3>>
  {
  public:
    enum
    {
      DIM = 1,
      DIM_SPACE = 3,
      DIM_ELEMENT = 3,
      DIM_DMAT = 3,
      DIFFORDER = 1,
    };

    static auto &Cast(const FiniteElement &fel)
    {
      return static_cast<const HDivDivFacetElement &>(fel);
    }
    static Array<int> GetDimensions() { return Array<int>({3,}); }

    template <typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement &fel,
                               const MIP &mip,
                               MAT &&mat, LocalHeap &lh)
    {
      Mat<3,3> jacobian = mip.GetJacobian();

      double factor = 1.0/sqr(Det(jacobian));
      
      switch (mip.IP().VB())
        {
        case BND:
          {
            Cast(fel).CalcDivShape(mip.IP(), Trans(mat));
            break;
          }
        case BBND:
          {
            Cast(fel).CalcDivShapeJump(mip.IP(), Trans(mat));
            break;
          }
        default:
          {
            throw Exception("only boundary integration points allowed");
          }
        }

      for (int i = 0; i < fel.GetNDof(); i++)
        {
          FlatVector<> shapei = mat.Col(i);
          Vec<3> trafo_shapei = factor*jacobian*Vec<3>(shapei);
          shapei = trafo_shapei;
        }
    }
  };





  
  
  




  

  
  /*
  Implementation of the FESpace
  */
  HDivDivFacetSpace::HDivDivFacetSpace(shared_ptr<MeshAccess> ama, const Flags &flags)
      : FESpace(ama, flags)
  {
    order = int(flags.GetNumFlag("order", 0));
    JKM = flags.GetDefineFlag("JKM");

    cout << "create HDivDivFacetspace, order=" << order
         << ", JKM = " << (JKM ? "true" : "false") << endl;
    
    switch (ma->GetDimension())
      {
      case 2:
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDDFacet<2>>>();
        additional_evaluators.Set("div", make_shared<T_DifferentialOperator<DiffOpDivHDDFacet<2>>>());        
        break;
      case 3:
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDDFacet<3>>>();
        additional_evaluators.Set("div", make_shared<T_DifferentialOperator<DiffOpDivHDDFacet<3>>>());
        break;
      default:
        throw Exception ("HDDFacet, unsupported dim");
      }
        
    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdHDDFacetBND>>();

    switch(ma->GetDimension())
      {
      case 2:
        {
          dofs_per_face = (order + 1);
          break;
        }
      case 3:
        {
          if (JKM)
            {
              dofs_per_face = (order>=1) ? 9 : 0;
              dofs_per_edge = 0;
            }
          else
            {
              dofs_per_face = TetNDofFace(order, true);
              dofs_per_edge = TetNDofEdge(order+1, true);
              // dofs_per_edge = 0;
            }              
          break;
        }
      default:
        throw Exception("HDviDivFacetspace not supported in dim="+ToString(ma->GetDimension()));
      }
    cout << "dofs_per_edge = " << dofs_per_edge << ", dofs_per_face = " << dofs_per_face << endl;
  }


  DocInfo HDivDivFacetSpace :: GetDocu ()
  {
    DocInfo docu = FESpace::GetDocu();
    docu.short_docu = "H(divdiv)-space in faces and edges.";

    docu.long_docu =
      R"raw_string(XXX.
)raw_string";      

    docu.Arg("JKM") = "bool = false\n"
      "  use JKM elements in faces";
    docu.Arg("bubble") = "bool = true\n"
      "  use bubbles in faces and edges";
    
    return docu;
  }
  

  

  void HDivDivFacetSpace::Update()
  {
    if (ma->GetDimension()==2)
      {
        int nfaces = ma->GetNFacets();
        SetNDof(dofs_per_face * nfaces);
        first_facet_dof = 0;
      }
    else
      {
        int nedges = ma->GetNEdges();
        int nfaces = ma->GetNFaces();
        SetNDof(dofs_per_edge*nedges + dofs_per_face * nfaces);
        first_facet_dof = dofs_per_edge*nedges;
      }
  }

  void HDivDivFacetSpace::GetDofNrs(ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize(0);
    if (ei.VB() != VOL && ei.VB() != BND)
      return;
    
    if (ma->GetDimension() == 2)
      {
        for (auto f : ma->GetElement(ei).Facets())
          for (int i = 0; i < dofs_per_face; i++)
            dnums.Append(dofs_per_face * f + i);
      }
    else
      {
        for (auto e : ma->GetElement(ei).Edges())
          for (int i = 0; i < dofs_per_edge; i++)
            dnums.Append(dofs_per_edge * e + i);
        for (auto f : ma->GetElement(ei).Faces())
          for (int i = 0; i < dofs_per_face; i++)
            dnums.Append(dofs_per_face * f + i);
      }
  }

  
  FiniteElement &HDivDivFacetSpace::GetFE(ElementId ei, Allocator &alloc) const
  {
    if (ei.IsVolume())
      {
        Ngs_Element ngel = ma->GetElement(ei);
        switch (ngel.GetType())
          {
          case ET_TET:
            {
              if (JKM)
                {
                  auto *fe = new (alloc) HDivDivFacetJKMElement();
                  fe->VertexOrientedFE<ET_TET>::SetVertexNumbers(ngel.vertices);
                  return *fe;
                }
              else
                {
                  auto *fe = new (alloc) HDivDivFacetHOElement<ET_TET>(order);
                  fe->VertexOrientedFE<ET_TET>::SetVertexNumbers(ngel.vertices);
                  return *fe;
                }
            }
          case ET_TRIG:
            {
              auto *fe = new (alloc) HDivDivFacetHOElement<ET_TRIG>(order);
              fe->VertexOrientedFE<ET_TRIG>::SetVertexNumbers(ngel.vertices);
              return *fe;          
            }
          default:
            throw Exception("element type "+ToString(ngel.GetType())+" not implemented");
          }
      }

    return SwitchET (ma->GetElement(ei).GetType(), [&] (auto et) -> FiniteElement&
    {
      return *new (alloc) DummyFE<et.ElementType()> ();
    });
  }





  // constraint space for highest-order div-constraints

  
  class DivConstraintTetElement : public HCurlFiniteElement<3>, public VertexOrientedFE<ET_TET>
  {
    
  public:
    DivConstraintTetElement (int aorder)
      : HCurlFiniteElement<3> (4*2*(aorder+1) + 6, aorder)
    {
      ;
    }

    ELEMENT_TYPE ElementType() const override { return ET_TET; }

    void CalcShape (const IntegrationPoint & ip, 
                    BareSliceMatrix<> shape) const override
    {
      shape.Rows(ndof).Cols(3) = 0;
      if (ip.VB() == BND)
        {
          int fanr = ip.FacetNr();
          double lami[4] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
          IVec<4> fav = ET_trait<ET_TET>::GetFaceSort (fanr, vnums);

          Vector dubiner((order+1)*(order+2)/2);
          auto p0 = ElementTopology::GetVertices(ET_TET)[fav[0]];
          auto p1 = ElementTopology::GetVertices(ET_TET)[fav[1]];
          auto p2 = ElementTopology::GetVertices(ET_TET)[fav[2]];

          Vec<3> t1, t2;
          for (int j = 0; j < 3; j++)
            {
              t1(j) = p1[j]-p0[j];
              t2(j) = p2[j]-p0[j];
            }
          
          DubinerBasis::Eval(order, lami[fav[0]], lami[fav[1]], dubiner);
          int ii = fanr*2*(order+1);
          int nr = -1;
          for (int j = 0; j <= order; j++)
            {
              nr += (order+1-j);
              shape.Row(ii++) = dubiner(nr)*t1;
              shape.Row(ii++) = dubiner(nr)*t2;
            }
        }
      else if (ip.VB() == BBND)
        {
          int ednr = ip.FacetNr();
          double lami[4] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
          IVec<2> edv = ET_trait<ET_TET>::GetEdgeSort (ednr, vnums);

          auto p0 = ElementTopology::GetVertices(ET_TET)[edv[0]];
          auto p1 = ElementTopology::GetVertices(ET_TET)[edv[1]];
          Vec<3> t;
          for (int j = 0; j < 3; j++) t(j) = p1[j]-p0[j];
          
          Vector leg(order+2);
          LegendrePolynomial::Eval(order+1, lami[edv[1]]-lami[edv[0]], leg);
          shape.Row(4*2*(order+1) + ednr) = leg(order+1)*t;
          
          // cout << "shape = " << shape.Rows(ndof).Cols(3) << endl;          
        }
      else
        throw Exception("DivConstraint only defined for BND and BBND");

      // cout << "shape = " << shape.Rows(ndof).Cols(3) << endl;
    }

    
  };

  class DivConstraintSpace : public FESpace
  {
  private:
    int nfacedofs;
    int nedgedofs;
  public:
    DivConstraintSpace(shared_ptr<MeshAccess> ama, const Flags &flags)
      : FESpace(ama, flags)
    {
      nfacedofs = 2*(order+1);
      nedgedofs = 1;

      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<3>>>();
    }
    
    virtual ~DivConstraintSpace() { }
    string GetClassName() const override { return "DivConstraintSpace"; }
    void Update() override
    {
      SetNDof(nfacedofs*ma->GetNFaces()+nedgedofs*ma->GetNEdges());
    }
    
    void GetDofNrs(ElementId ei, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
      if (ei.IsVolume())
        {
          Ngs_Element ngel = ma->GetElement(ei);
          for (auto f : ngel.Faces())
            dnums += IntRange(f*nfacedofs, (f+1)*nfacedofs);
          size_t base = nfacedofs*ma->GetNFaces();
          for (auto e : ngel.Edges())
            dnums += IntRange(base+e*nedgedofs, base+(e+1)*nedgedofs);
        }
    }
    
    FiniteElement &GetFE(ElementId ei, Allocator &alloc) const override
    {
      if (ei.IsVolume())
        {
          Ngs_Element ngel = ma->GetElement(ei);
          switch (ngel.GetType())
            {
            case ET_TET:
              {
                auto *fe = new (alloc) DivConstraintTetElement(order);
                return *fe;
              }
            default:
              ;
            }
        }

      if (ei.IsBoundary())    return SwitchET (ma->GetElement(ei).GetType(), [&] (auto et) -> FiniteElement&
      {
        return *new (alloc) DummyFE<et.ElementType()> ();
      });
      
      throw Exception("DivConstraintSpace: no element");
    }
  };
  
  shared_ptr<FESpace> HDivDivFacetSpace::GetDivConstraintSpace() const
  {
    Flags flags;
    flags.SetFlag("order", order);
    return make_shared<DivConstraintSpace>(ma, flags);
  }
};

