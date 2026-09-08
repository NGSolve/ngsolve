#ifndef DIFFOPS_HPP
#define DIFFOPS_HPP

#include <bla.hpp>
#include <scalarfe.hpp>
#include <hcurlfe.hpp>
#include <hdivfe.hpp>
#include <diffop_impl.hpp>


namespace ngsbem
{
  using namespace ngbla;
  using namespace ngfem;


  class DiffOpBoundaryRot : public DiffOp<DiffOpBoundaryRot>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 3 };
    enum { DIFFORDER = 1 };

    static bool SupportsVB (VorB checkvb) { return checkvb==BND; }
    
    static string Name() { return "boundaryrot"; }
    
    static const ScalarFiniteElement<2> & Cast (const FiniteElement & fel) 
    { return static_cast<const ScalarFiniteElement<2>&> (fel); }

    ///
    // mat is 3 x ndof
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcMappedDShape (mip, Trans(mat));
      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Vec<3> grad = mat.Col(i);
          mat.Col(i) = Cross(mip.GetNV(), grad);
        }
    }

    static int DimRef() { return 2; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcDShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & bmip,
                                          MAT & mat, LocalHeap & lh)
    {
      auto & mip = static_cast<const MappedIntegrationPoint<2,3>&>(bmip); 
      Vec<3> nv = mip.GetNV();
      mat = Trans(mip.GetJacobianInverse());

      for (int j = 0; j < 2; j++)
        mat.Col(j) = Cross(nv, Vec<3> (mat.Col(j)));
    }
    

    /// mat is (ndof*3) x mip.Size()
    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedDShape (mir, mat);

      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<3,3>>&>(mir[j]).GetNV();
          for (int i = 0; i < fel.GetNDof(); i++)
            {
              Vec<3,SIMD<double>> grad = mat.Col(j).Range(3*i, 3*i+3);
              mat.Col(j).Range(3*i,3*i+3) = Cross(nv, grad);
            }
        }
    }
  };




  class DiffOpRotatedTrace : public DiffOp<DiffOpRotatedTrace>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 3 };
    enum { DIFFORDER = 1 };

    static string Name() { return "rotatedtrace"; }
    static int DimRef() { return 2; }

    static const HCurlFiniteElement<2> & Cast (const FiniteElement & fel) 
    { return static_cast<const HCurlFiniteElement<2>&> (fel); }

    // mat is 2 x ndof
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                    MAT && mat, LocalHeap & lh)
    {
        Cast(fel).CalcShape (ip, Trans(mat));
    }

    ///
    // mat is 3 x ndof
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcMappedShape (mip, Trans(mat));
      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Vec<3> shape = mat.Col(i);
          mat.Col(i) = Cross(mip.GetNV(), shape);
        }
    }

    // mat is 3 x 2
    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & bmip,
                                          MAT & mat, LocalHeap & lh)
    {
        auto & mip = static_cast<const MappedIntegrationPoint<2,3>&>(bmip);
        Vec<3> nv = mip.GetNV();
        mat = Trans(mip.GetJacobianInverse());

        for (int j = 0; j < 2; j++)
          mat.Col(j) = Cross(nv, Vec<3> (mat.Col(j)));
    }

    /// mat is (ndof*3) x mip.Size()
    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedShape (mir, mat);

      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<2,3>>&>(mir[j]).GetNV();
          for (int i = 0; i < fel.GetNDof(); i++)
            {
              Vec<3,SIMD<double>> shape = mat.Col(j).Range(3*i, 3*i+3);
              mat.Col(j).Range(3*i,3*i+3) = Cross(nv, shape);
            }
        }
    }

    using DiffOp<DiffOpRotatedTrace>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(fel).Evaluate (mir, x, y);
      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<2,3>>&>(mir[j]).GetNV();
          Vec<3,SIMD<double>> val = y.Col(j).Range(0,3);
          y.Col(j).Range(0,3) = Cross(nv, val);
        }
    }

    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      Cast(fel).Evaluate (mir, x, y);
      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<2,3>>&>(mir[j]).GetNV();
          Vec<3,SIMD<Complex>> val = y.Col(j).Range(0,3);
          y.Col(j).Range(0,3) = Cross(nv, val);
        }
    }
  };

}


#endif
