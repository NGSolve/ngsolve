

#include <comp.hpp>

#include "Hddf.hpp"

using namespace ngcomp;

/*
Implementation of the Differential Operators
*/

namespace ngcomp
{

  class HDivDivFacetSpace;
  class HDivDivFacetElement;

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

  
  // a x b = Skw(a) b  .... changed sign !!!
  inline Mat<3, 3> Skw(Vec<3> v)
  {
    Mat<3, 3> m  { {  0, -v[2], v[1] },
                   { v[2], 0, -v[0] },
                   { -v[1], v[0], 0} };
    return m;
  }

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


  template <int D>
  inline Mat<D, D> OuterProduct(const Vec<D> &a, const Vec<D> &b)
  {
    Mat<D, D> result;
    for (int i = 0; i < D; i++)
      for (int j = 0; j < D; j++)
        result(i, j) = a(i) * b(j);
    return result;
  }

  inline Mat<3, 3> SymOuterProduct(const Vec<3> &a, const Vec<3> &b)
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

    virtual void CalcShapes(const IntegrationPoint &ip, SliceMatrix<double> shapes) const = 0;
    virtual void CalcEdgeDivShapes(const IntegrationPoint &ip, SliceMatrix<double> curl_e_shapes) const
    {
      if (ip.VB() != BBND) // VB() is the enum of the integration point = {VOL, BND, BBND, BBBND}
        throw Exception("only BBoundary integration points allowed");

      if (Dim() == 3)
        {
          throw Exception("HDivDivFacet 3D not ready");
          int enr = ip.FacetNr();
          IVec<2> nfaces = e2f_oriented()[enr];
          
          curl_e_shapes = 0;
          
          for (int j = 0; j < 2; j++)
            {
              IntegrationPoint ip_face = ip;
              ip_face.SetFacetNr(nfaces[j], BND);
              
              Matrix<double> shape_f(ndof, 9);
              CalcShapes(ip_face, shape_f /* , jacobian*/);
              
              Vec<3> t = GetTangents()[enr];
              Vec<3> n = ElementTopology::GetNormals<3>(ET_TET)[nfaces[j]];
              
              Vec<3> nu = Cross(n, t);
              
              double sign = (j == 0) ? -1 : 1;
              Mat<3,3> Cnu_ = Skw(nu);
              
              double factor = 1; //  / L2Norm((jacobian * t));
              factor /= sqr(L2Norm(n));
              
              for (int k = 0; k < ndof; k++)
                {
                  FlatMatrix<> shapei = shape_f.Row(k).AsMatrix(3,3);
                  // Mat<3, 3> shapeCnu = -Cnu_ * shapei; // foes not need the transpose (look at the theory)
                  Mat<3, 3> shapeCnu = shapei * Cnu_; // foes not need the transpose (look at the theory)
                  shapei = shapeCnu;
                  curl_e_shapes.Row(k) += sign * factor * shape_f.Row(k);
                }
            }
        }
      else
        {
          int vnr = ip.FacetNr();
          IVec<2> nedges = v2e_oriented()[vnr];
          
          curl_e_shapes = 0;
          
          for (int j = 0; j < 2; j++)
            {
              IntegrationPoint ip_face = ip;
              ip_face.SetFacetNr(nedges[j], BND);
              
              Matrix<double> shape_f(ndof, 4);
              CalcShapes(ip_face, shape_f);
              
              Vec<2> n = ElementTopology::GetNormals<2>(ET_TRIG)[nedges[j]];
              double sign = (j == 0) ? -1 : 1;
              double factor = 1/sqr(L2Norm(n));
              Vec<2> tau(-n(1),n(0));

              for (int k = 0; k < ndof; k++)
                {
                  FlatMatrix<> shapei = shape_f.Row(k).AsMatrix(2,2);
                  Vec<2> shapei_tau = shapei * tau;
                  curl_e_shapes.Row(k) += -sign * factor * shapei_tau;
                }
            }
          // *testout << "edge curl shape, ip = " << ip << endl
          // << "shape = " << curl_e_shapes << endl;
        }
    }

    virtual void CalcFaceDivShapes(const IntegrationPoint &ip, SliceMatrix<double> curl_f_shapes) const = 0;
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
                            ET==ET_TET ? 2 * (order + 1) * (order + 2) : 0,
                            order) {}

    ~HDivDivFacetHOElement() {}

    ELEMENT_TYPE ElementType() const override
    {
      return ET;
    }

    void CalcShapes(const IntegrationPoint &ip, SliceMatrix<double> shapes) const override
    {
      throw Exception("HCCFacet, element type " + ToString(ET) + " not implemented");
    }

    virtual void CalcFaceDivShapes(const IntegrationPoint &ip, SliceMatrix<double> div_shapes) const override
    {
      throw Exception("HDDFacet, element type " + ToString(ET) + " FaceDivShape not implemented");
    }

  };

/*
template <>
void HDivDivFacetHOElement<ET_TET> :: CalcShapes(const IntegrationPoint &ip, SliceMatrix<double> shapes) const 
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
void HDivDivFacetHOElement<ET_TRIG> :: CalcShapes(const IntegrationPoint &ip, SliceMatrix<double> shapes) const 
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
  CalcFaceDivShapes(const IntegrationPoint &ip, SliceMatrix<double> div_shapes) const
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
    
    Mat<2,2> tautau = OuterProduct(tau, tau);
    
    for (int i = 0; i <= order; i++)
      {
        Vec<2> grad = GetGradient(values[i]);
        div_shapes.Row(dofs_per_face * enr + i) = InnerProduct(tau,grad)*tau;
      }
  }


  




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
      Cast(fel).CalcShapes(mip.IP(), Trans(mat));

      Vec<D> n = ElementTopology::GetNormals<D>(mip.GetTransformation().GetElementType())[mip.IP().FacetNr()];
      Vec<2> tau(-n(1),n(0));
      double factor = 1 / L2Norm(Cof(jacobian) * n);

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
      Cast(fel).CalcShapes(mip.IP(), Trans(mat) /* , jacobian*/);
      for (int i = 0; i < fel.GetNDof(); i++)
      {
        FlatMatrix<> shapei(3, 3, &mat(0, i));                                                          // &mat(i, 0)
        Mat<3, 3> trafo_shapei = Cof(mip.GetJacobian()) * Mat<3, 3>(shapei) * Trans(mip.GetJacobian()); // Mat<3, 3> trafo_shapei = (1 / mip.GetJacobiDet()) * mip.GetJacobianInverse() * shapei * mip.GetJacobian(); // check if this is correct
        shapei = trafo_shapei;
      }
    }
  };
  // define the square root

  template <int DIM> class DiffOpDivHDDFacet;


  /*
  template<>
  class DiffOpCurlHCCFacet<3> : public DiffOp<DiffOpCurlHCCFacet<3>>
  {
  public:
    enum
    {
      DIM = 1,
      DIM_SPACE = 3,
      DIM_ELEMENT = 3,
      DIM_DMAT = 9,
      DIFFORDER = 1,
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

      double factor = 1;

      switch (mip.IP().VB())
      {
      case BND:
      {
        Vec<3> n = ElementTopology::GetNormals<3>(ET_TET)[mip.IP().FacetNr()];
        factor = 1 / L2Norm(Cof(mip.GetJacobian()) * n);
        // factor = 1;
        Cast(fel).CalcFaceCurlShapes(mip.IP(), Trans(mat));
        break;
      }
      case BBND:
      {
        Vec<3> t = GetTangents()[mip.IP().FacetNr()];
        factor = 1 / L2Norm(jacobian * t);
        Cast(fel).CalcEdgeCurlShapes(mip.IP(), Trans(mat));
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
        FlatMatrix<> shapei = mat.Col(i).AsMatrix(3,3);
        Mat<3, 3> trafo_shapei = factor * Cof(jacobian)*Mat<3, 3>(shapei) * Trans(jacobian);
        shapei = trafo_shapei;
      }
    }
  };
  */
  
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
        Cast(fel).CalcFaceDivShapes(mip.IP(), Trans(mat));

        Vec<2> t = GetTangentsTrig()[mip.IP().FacetNr()];
        factor = 1 / L2Norm(jacobian * t);
        break;
      }
      case BBND:
      {
        Cast(fel).CalcEdgeDivShapes(mip.IP(), Trans(mat));
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






  
  /*
  Implementation of the FESpace
  */
  HDivDivFacetSpace::HDivDivFacetSpace(shared_ptr<MeshAccess> ama, const Flags &flags)
      : FESpace(ama, flags)
  {
    order = int(flags.GetNumFlag("order", 0));

    switch (ma->GetDimension())
      {
      case 2:
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDDFacet<2>>>();
        additional_evaluators.Set("div", make_shared<T_DifferentialOperator<DiffOpDivHDDFacet<2>>>());        
        break;
        /*
      case 3:
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHCCFacet<3>>>();
        additional_evaluators.Set("curl", make_shared<T_DifferentialOperator<DiffOpCurlHCCFacet<3>>>());
        break;
        */
      default:
        throw Exception ("HDDFacet, unsupported dim");
      }
        
    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdHDDFacetBND>>();


    if (ma->GetDimension()==2)
      dofs_per_face = (order + 1);
    else
      dofs_per_face = (order + 1) * (order + 2) / 2;
  }

  void HDivDivFacetSpace::Update()
  {
    int nfaces = ma->GetNFacets();
    SetNDof(dofs_per_face * nfaces);
  }

  void HDivDivFacetSpace::GetDofNrs(ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize(0);
    if (ei.VB() != VOL && ei.VB() != BND)
      return;

    for (auto f : ma->GetElement(ei).Facets())
      for (int i = 0; i < dofs_per_face; i++)
        dnums.Append(dofs_per_face * f + i);
  }
  
  FiniteElement &HDivDivFacetSpace::GetFE(ElementId ei, Allocator &alloc) const
  {
    if (ei.IsVolume())
      {
        Ngs_Element ngel = ma->GetElement(ei);
        switch (ngel.GetType())
          {
            /*
          case ET_TET:
            {
              auto *fe = new (alloc) HDivDivFacetHOElement<ET_TET>(order);
              fe->VertexOrientedFE<ET_TET>::SetVertexNumbers(ngel.vertices);
              return *fe;          
            }
            */
          case ET_TRIG:
            {
              auto *fe = new (alloc) HDivDivFacetHOElement<ET_TRIG>(order);
              fe->VertexOrientedFE<ET_TRIG>::SetVertexNumbers(ngel.vertices);
              return *fe;          
            }
          default:
            throw Exception("element type"+ToString(ngel.GetType())+" not implemented");
          }
      }

    return SwitchET (ma->GetElement(ei).GetType(), [&] (auto et) -> FiniteElement&
    {
      return *new (alloc) DummyFE<et.ElementType()> ();
    });
  }
};

