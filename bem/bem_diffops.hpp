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
  };

  // rotated + scalar
  class DiffOpHelmholtz : public DiffOp<DiffOpHelmholtz>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 4 };
    enum { DIFFORDER = 1 };

    static string Name() { return "Helmholtz"; }
    static int DimRef() { return 3; }

    static const ScalarFiniteElement<2> & Cast (const FiniteElement & fel) 
    { return static_cast<const ScalarFiniteElement<2>&> (fel); }

    // mat is 3xndof
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                    MAT && mat, LocalHeap & lh)
    {
        auto matvec = mat.Rows(0,2);
        Cast(fel).CalcDShape (ip, Trans(matvec));
        Cast(fel).CalcShape(ip, mat.Row(2));
    }

    ///
    // mat is 4 x ndof
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      // mat.AddSize(4, fel.GetNDof()) = 0.0;
      auto matvec = mat.Rows(0,3);
      Cast(fel).CalcMappedDShape (mip, Trans(matvec));
      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Vec<3> grad = matvec.Col(i);
          matvec.Col(i) = Cross(mip.GetNV(), grad);
        }
      // *testout << "mat1 = " << mat << endl;
      // mat.AddSize(4, fel.GetNDof()) = 0.0;
      Cast(fel).CalcShape(mip.IP(), mat.Row(3));
      // *testout << "scalar mat = " << endl << mat << mat;
      // *testout << "mat2 = " << mat << endl;      
    }

    // mat is 4x3
    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & bmip,
                                          MAT & mat, LocalHeap & lh)
    {
      auto & mip = static_cast<const MappedIntegrationPoint<2,3>&>(bmip); 
      Vec<3> nv = mip.GetNV();
      mat = 0.0;
      mat.Rows(0,3).Cols(0,2) = Trans(mip.GetJacobianInverse());
      auto gradmat = mat.Rows(0,3).Cols(0,2);
      for (int j = 0; j < 2; j++)
        gradmat.Col(j) = Cross(nv, Vec<3> (gradmat.Col(j)));

      mat(3,2) = 1.0;
    }

    /// mat is (ndof*4) x mip.Size()
    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedDShape (mir, mat.Rows(0, 3*fel.GetNDof()));

      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<3,3>>&>(mir[j]).GetNV();
	  for (int i = fel.GetNDof()-1; i >= 0; i--)
            {
	      Vec<3,SIMD<double>> grad = mat.Col(j).Range(3*i, 3*i+3);
              mat.Col(j).Range(4*i,4*i+3) = Cross(nv, grad);
            }
        }
      // *testout << "simd mat1 = " << endl << mat.AddSize(4*fel.GetNDof(), mir.Size()) << endl;            
      // mat.AddSize(4*fel.GetNDof(), mir.Size()) = SIMD<double>(0.0);
      Cast(fel).CalcShape (mir.IR(), mat.RowSlice(3,4));
    }
  };



  // rotated 
  class DiffOpMaxwell : public DiffOp<DiffOpMaxwell>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 4 };
    enum { DIFFORDER = 1 };

    static string Name() { return "Maxwell"; }
    static int DimRef() { return 3; }

    static const HCurlFiniteElement<2> & Cast (const FiniteElement & fel) 
    { return static_cast<const HCurlFiniteElement<2>&> (fel); }

    // mat is 3 x ndof
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                    MAT && mat, LocalHeap & lh)
    {
      // auto matvec = mat.Rows(0,2);
        Cast(fel).CalcShape (ip, Trans(mat));
        Cast(fel).CalcCurlShape(ip, Trans(mat.Rows(2,3)));
    }

    ///
    // mat is 4 x ndof
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      // mat.AddSize(4, fel.GetNDof()) = 0.0;
      auto matvec = mat.Rows(0,3);
      Cast(fel).CalcMappedShape (mip, Trans(matvec));
      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Vec<3> shape = matvec.Col(i);
          matvec.Col(i) = Cross(mip.GetNV(), shape);
        }
      // *testout << "mat1 = " << mat << endl;
      // mat.AddSize(4, fel.GetNDof()) = 0.0;      
      mat.Row(3) =
        1.0/mip.GetJacobiDet() * 
	Cast(fel).GetCurlShape(mip.IP(),lh).Col(0);
      // *testout << "scalar mat = " << endl << mat << mat;
      // *testout << "mat2 = " << mat << endl;      
    }

    // mat is 4x3
    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & bmip,
                                          MAT & mat, LocalHeap & lh)
    {
      auto & mip = static_cast<const MappedIntegrationPoint<2,3>&>(bmip);
      Vec<3> nv = mip.GetNV();
      mat = 0.0;
      mat.Rows(0,3).Cols(0,2) = Trans(mip.GetJacobianInverse());
      auto matvec = mat.Rows(0,3).Cols(0,2);
      for (int j = 0; j < 2; j++)
        matvec.Col(j) = Cross(nv, Vec<3> (matvec.Col(j)));

      mat(3,2) = 1.0 / mip.GetJacobiDet();
    }

    /// mat is (ndof*4) x mip.Size()
    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedShape (mir, mat.Rows(0, 3*fel.GetNDof()));
      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<2,3>>&>(mir[j]).GetNV();
          for (int i = fel.GetNDof()-1; i >= 0; i--)
            {
              Vec<3,SIMD<double>> shape = mat.Col(j).Range(3*i, 3*i+3);
              mat.Col(j).Range(4*i,4*i+3) = Cross(nv, shape);
              // mat.Col(j).Range(4*i+3,4*i+4) = SIMD<double>(0.0);
            }
        }
      // *testout << "simd mat1 = " << endl << mat.AddSize(4*fel.GetNDof(), mir.Size()) << endl;            
      // mat.AddSize(4*fel.GetNDof(), mir.Size()) = SIMD<double>(0.0);
      
      constexpr size_t BS=16;
      LocalHeapMem<BS*SIMD<double>::Size()*sizeof(SIMD<MappedIntegrationPoint<2,2>>)+64> lh("genmatlh");
      FE_ElementTransformation<2,2> trafo2d(fel.ElementType());
      for (size_t first = 0; first < mir.Size(); first += BS)
        {
          HeapReset hr(lh);
          size_t next = std::min(first+BS, mir.Size());
          SIMD_MappedIntegrationRule<2,2> mir2d(mir.IR().Range(first, next), trafo2d, lh);
          Cast(fel).CalcMappedCurlShape (mir2d, mat.RowSlice(3,4).Cols(first, next));
        }
      for (size_t i = 0; i < mir.Size(); i++)
        mat.Col(i).Slice(3,4).Range(fel.GetNDof()) *= 1.0 / mir[i].GetJacobiDet();
      // *testout << "simd mat2 = " << endl << mat.AddSize(4*fel.GetNDof(), mir.Size()) << endl;      
    }

    using DiffOp<DiffOpMaxwell>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      Cast(fel).Evaluate (mir, x, y.Rows(3));
      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<2,3>>&>(mir[j]).GetNV();
          Vec<3,SIMD<Complex>> shape = y.Col(j).Range(0,3);
          y.Col(j).Range(0,3) = Cross(nv, shape);
        }
      y.Row(3).Range(mir.Size()) = SIMD<Complex>(0.0);
      // TODO: curl part
    }
  };






  // copied from ngsolve/fem/hdiv_equations.hpp, DiffOpIdHDivSurface
  class DiffOpMaxwellNew : public DiffOp<DiffOpMaxwellNew>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 4 };
    enum { DIFFORDER = 1 };

    static string Name() { return "Maxwell"; }
    static int DimRef() { return 3; }
    
    static const HDivFiniteElement<2> & Cast (const FiniteElement & fel) 
    { return static_cast<const HDivFiniteElement<2>&> (fel); }

    // mat is 3xndof
    template<typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      auto matvec = mat.Rows(0,2); // 2D vector field
      Cast(fel).CalcShape (ip, Trans(matvec));
      Cast(fel).CalcDivShape (ip, mat.Row(2)); //divergence
    }

    ///
    // mat is 4 x ndof
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      auto matvec = mat.Rows(0,3);
      matvec = (1.0 / mip.GetJacobiDet()) *mip.GetJacobian () * 
        Trans (Cast(fel).GetShape(mip.IP(),lh));

      mat.Row(3) =
        1.0/mip.GetJacobiDet() * 
	Cast(fel).GetDivShape(mip.IP(),lh).Col(0);
    }

    // mat is 4x3
    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & bmip,
                                          MAT & mat, LocalHeap & lh)
    {
        auto & mip = static_cast<const MappedIntegrationPoint<2,3>&>(bmip);
        mat = 0.0;
        mat.Rows(0,3).Cols(0,2) = (1.0 / mip.GetJacobiDet()) * mip.GetJacobian();
        mat(3,2) = 1.0 / mip.GetJacobiDet();
    }

    /// mat is (ndof*4) x mip.Size()
    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedShape (mir, mat.Rows(0, 3*fel.GetNDof()));
      for (int j = 0; j < mir.Size(); j++)
        {
          for (int i = fel.GetNDof()-1; i >= 0; i--)
            {
              Vec<3,SIMD<double>> shape = mat.Col(j).Range(3*i, 3*i+3);
              mat.Col(j).Range(4*i,4*i+3) = shape;
            }
        }

      Cast(fel).CalcMappedDivShape (mir, mat.RowSlice(3, 4));
      /*
      constexpr size_t BS=16;
      LocalHeapMem<BS*SIMD<double>::Size()*sizeof(SIMD<MappedIntegrationPoint<2,2>>)+64> lh("genmatlh");
      FE_ElementTransformation<2,2> trafo2d(fel.ElementType());
      for (size_t first = 0; first < mir.Size(); first += BS)
        {
          HeapReset hr(lh);
          size_t next = std::min(first+BS, mir.Size());
          SIMD_MappedIntegrationRule<2,2> mir2d(mir.IR().Range(first, next), trafo2d, lh);
          Cast(fel).CalcMappedDivShape (mir2d, mat.RowSlice(3,4).Cols(first, next));
        }
      for (size_t i = 0; i < mir.Size(); i++)
        mat.Col(i).Slice(3,4).Range(fel.GetNDof()) *= 1.0 / mir[i].GetJacobiDet();
      */
    }

    using DiffOp<DiffOpMaxwellNew>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      Vector<> xr(fel.GetNDof()), xi(fel.GetNDof());
      Matrix<SIMD<double>> valr(3,mir.Size()), vali(3, mir.Size());

      xr = Real(x);
      xi = Imag(x);
      Cast(fel).Evaluate (mir, xr, valr);
      Cast(fel).Evaluate (mir, xi, vali);
      
      for (int j = 0; j < mir.Size(); j++)
        for (int k = 0; k < 3; k ++)
          y(k,j) = SIMD<Complex> (valr(k,j), vali(k,j));
      y.Row(3).Range(mir.Size()) = SIMD<Complex>(0.0);
    }
  };


  
}


#endif
