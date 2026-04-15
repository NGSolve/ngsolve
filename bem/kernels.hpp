#ifndef KERNELS_hpp
#define KERNELS_hpp

#include "mptools.hpp"
#include <type_traits>
#include <variant>


//  ****************************  The kernels **********************************

namespace ngsbem
{

  struct KernelTerm
  {
    double fac;
    size_t kernel_comp;
    size_t trial_comp;
    size_t test_comp;
  };
  

  class BaseKernel
  {
  public:
    void GetDifferentiatedKernel(const string &name) const {
      throw Exception("GetDifferentiatedKernel not overloaded");
    }
  };

  // *********** FMM Source/Target Interface **********************

  template <typename mp_type, typename value_type, typename T_Kappa = double>
  class FMMInterface
  {
  public:
    T_Kappa kappa;

    FMMInterface (T_Kappa _kappa = 1e-16) : kappa(_kappa) {}

    virtual shared_ptr<SingularMLExpansion<mp_type, T_Kappa>> CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<SingularMLExpansion<mp_type, T_Kappa>> (c, r, kappa, fmm_params);
    }

    virtual shared_ptr<RegularMLExpansion<mp_type, T_Kappa>> CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<RegularMLExpansion<mp_type, T_Kappa>> (c, r, kappa, fmm_params);
    }

    virtual void AddSource (SingularMLExpansion<mp_type, T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const
    {
      throw Exception("AddSource not implemented for this FMM type");
    }

    virtual void EvaluateMP (RegularMLExpansion<mp_type, T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const
    {
      throw Exception("EvaluateMP not implemented for this FMM type");
    }

    virtual ~FMMInterface() = default;
  };

  template <int COMPS, typename value_type, typename T_Kappa = double>
  class Charges : public FMMInterface<Vec<COMPS, Complex>, value_type, T_Kappa>
  {
  public:
    using mp_type = Vec<COMPS, Complex>;
    using Base = FMMInterface<mp_type, value_type, T_Kappa>;
    static constexpr bool needs_normal = false;
    Charges (T_Kappa kappa = 1e-16) : Base(kappa) {}

    void AddSource (SingularMLExpansion<mp_type, T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      mp.AddCharge (pnt, val);
    }

    void EvaluateMP (RegularMLExpansion<mp_type, T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      if constexpr (std::is_same_v<value_type, double>)
        val = Real(mp.Evaluate (pnt));
      else
        val = mp.Evaluate (pnt);
    }
  };

  template <int COMPS, typename value_type, typename T_Kappa = double>
  class Dipoles : public FMMInterface<Vec<COMPS, Complex>, value_type, T_Kappa>
  {
  public:
    using mp_type = Vec<COMPS, Complex>;
    using Base = FMMInterface<mp_type, value_type, T_Kappa>;
    static constexpr bool needs_normal = true;
    Dipoles (T_Kappa kappa = 1e-16) : Base(kappa) {}

    void AddSource (SingularMLExpansion<mp_type, T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      mp.AddDipole (pnt, -nv, val);
    }

    void EvaluateMP (RegularMLExpansion<mp_type, T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      if constexpr (std::is_same_v<value_type, double>)
        val = Real(mp.EvaluateDirectionalDerivative (pnt, nv));
      else
        val = mp.EvaluateDirectionalDerivative (pnt, nv);
    }
  };

  // Evaluates gradient: 3 directional derivatives along unit vectors
  template <int COMPS, typename value_type, typename T_Kappa = double>
  class GradientEval : public FMMInterface<Vec<COMPS, Complex>, value_type, T_Kappa>
  {
  public:
    using mp_type = Vec<COMPS, Complex>;
    using Base = FMMInterface<mp_type, value_type, T_Kappa>;
    static constexpr bool needs_normal = false;
    GradientEval (T_Kappa kappa = 1e-16) : Base(kappa) {}

    void EvaluateMP (RegularMLExpansion<mp_type, T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      for (int i = 0; i < 3; i++)
        {
          Vec<3> ei = 0;
          ei(i) = 1;
          auto deri = mp.EvaluateDirectionalDerivative (pnt, ei);
          for (int c = 0; c < COMPS; c++)
          {
            if constexpr (std::is_same_v<value_type, double>)
              val(3*c+i) = Real(deri(c));
            else
              val(3*c+i) = deri(c);
          }
        }
    }
  };

  // Direct evaluation: val = mp.Evaluate(pnt)
  template <typename mp_type, typename value_type, typename T_Kappa = double>
  class DirectEval : public FMMInterface<mp_type, value_type, T_Kappa>
  {
    using Base = FMMInterface<mp_type, value_type, T_Kappa>;
  public:
    static constexpr bool needs_normal = false;
    DirectEval (T_Kappa kappa = 1e-16) : Base(kappa) {}

    void EvaluateMP (RegularMLExpansion<mp_type, T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      if constexpr (std::is_same_v<mp_type, Complex>)
        val(0) = mp.Evaluate (pnt);
      else
        val = mp.Evaluate (pnt);
    }
  };

  // Combined charge + dipole source (CombinedField)
  template <int COMPS, typename T_Kappa = double>
  class ChargeDipoles : public FMMInterface<Vec<COMPS, Complex>, Complex, T_Kappa>
  {
  public:
    using mp_type = Vec<COMPS, Complex>;
    using Base = FMMInterface<mp_type, Complex, T_Kappa>;
    static constexpr bool needs_normal = true;
    ChargeDipoles (T_Kappa kappa) : Base(kappa) {}

    void AddSource (SingularMLExpansion<mp_type, T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const override
    {
      mp.AddChargeDipole (pnt, -this->kappa * Complex(0, 1)*val, -nv, val);
    }
  };

  // HelmholtzHS source: packs charges + kappa²-scaled normal dipoles
  class HelmholtzHSSource : public FMMInterface<Vec<6,Complex>, Complex>
  {
    double kappa;
    using Base = FMMInterface<Vec<6,Complex>, Complex>;
  public:
    static constexpr bool needs_normal = true;
    HelmholtzHSSource (double _kappa) : Base(_kappa), kappa(_kappa) {}

    void AddSource (SingularMLExpansion<Vec<6,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const override
    {
      Vec<6,Complex> charge;
      charge.Range(0,3) = val.Range(0,3);
      charge.Range(3,6) = -kappa * kappa * val(3) * nv;
      mp.AddCharge(pnt, charge);
    }
  };

  // HelmholtzHS target: unpacks eval + dot product with normal
  class HelmholtzHSTarget : public FMMInterface<Vec<6,Complex>, Complex>
  {
    using Base = FMMInterface<Vec<6,Complex>, Complex>;
  public:
    static constexpr bool needs_normal = true;
    HelmholtzHSTarget (double kappa) : Base(kappa) {}

    void EvaluateMP (RegularMLExpansion<Vec<6,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const override
    {
      Vec<6,Complex> eval = mp.Evaluate (pnt);
      val.Range(0,3) = eval.Range(0,3);
      val(3) = InnerProduct(eval.Range(3,6), nv);
    }
  };

  // MaxwellSL source: kappa-scaled 4-component charges
  class MaxwellSLSource : public FMMInterface<Vec<4,Complex>, Complex>
  {
    double kappa;
    using Base = FMMInterface<Vec<4,Complex>, Complex>;
  public:
    static constexpr bool needs_normal = false;
    MaxwellSLSource (double _kappa) : Base(_kappa), kappa(_kappa) {}

    void AddSource (SingularMLExpansion<Vec<4,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const override
    {
      Vec<4,Complex> charge;
      charge.Range(0,3) = kappa * val.Range(0, 3);
      charge(3) = -1.0/kappa * val(3);
      mp.AddCharge(pnt, charge);
    }
  };

  // MaxwellDL source/target: cross-product dipoles + direct eval (self-symmetric)
  template <typename T_Kappa = double>
  class MaxwellCurlDipoles : public FMMInterface<Vec<3,Complex>, Complex, T_Kappa>
  {
    using Base = FMMInterface<Vec<3,Complex>, Complex, T_Kappa>;
  public:
    static constexpr bool needs_normal = false;
    MaxwellCurlDipoles (T_Kappa kappa) : Base(kappa) {}

    void AddSource (SingularMLExpansion<Vec<3,Complex>,T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const override
    {
      Vec<3,Complex> n_cross_m = val.Range(0, 3);
      for (int k = 0; k < 3; k++)
        {
          Vec<3> ek{0.0}; ek(k) = 1;
          Vec<3> n_cross_m_real = Real(n_cross_m);
          Vec<3> n_cross_m_imag = Imag(n_cross_m);
          mp.AddDipole(pnt, Cross(n_cross_m_real, ek), ek);
          mp.AddDipole(pnt, Cross(n_cross_m_imag, ek), Complex(0,1)*ek);
        }
    }

    void EvaluateMP (RegularMLExpansion<Vec<3,Complex>,T_Kappa> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const override
    {
      val = mp.Evaluate (pnt);
    }
  };

  // Lame source: charges + Jacobi matrix dipoles
  class LameSource : public FMMInterface<Vec<6,Complex>, double>
  {
    using Base = FMMInterface<Vec<6,Complex>, double>;
  public:
    static constexpr bool needs_normal = false;
    LameSource () : Base(1e-16) {}

    void AddSource (SingularMLExpansion<Vec<6,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const override
    {
      Vec<6> charge = 0.0;
      charge.Range(0,3) = val;
      mp.AddCharge(pnt, charge);

      Mat<3,3> jacobi = OuterProduct(pnt, val) + InnerProduct(pnt, val) * Id<3>();

      for (int k = 0; k < 3; k++)
        {
          Vec<6> dipole_charge = 0.0;
          dipole_charge.Range(3,6) = jacobi.Col(k);
          auto ek = UnitVec<3>(k);
          mp.AddDipole(pnt, -ek, dipole_charge);
        }
    }
  };

  // Lame target: multi-step eval with coefficients
  class LameTarget : public FMMInterface<Vec<6,Complex>, double>
  {
    double nu, alpha;
    using Base = FMMInterface<Vec<6,Complex>, double>;
  public:
    static constexpr bool needs_normal = false;
    LameTarget (double _nu, double _alpha) : Base(1e-16), nu(_nu), alpha(_alpha) {}

    void EvaluateMP (RegularMLExpansion<Vec<6,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const override
    {
      Vec<6> mpval = Real(mp.Evaluate (pnt));
      val.Range(0,3) = 0;
      val += (3-4*nu)*alpha * mpval.Range(0,3);

      val -= alpha/2 * mpval.Range(3,6);

      Mat<3,3> jacobi = 0.0;
      for (int k = 0; k < 3; k++)
        {
          auto ek = UnitVec<3>(k);
          jacobi.Col(k) = Real(mp.EvaluateDirectionalDerivative(pnt, ek).Range(0,3));
        }

      val -= alpha/2 * ( Trans(jacobi) * pnt + Trace(jacobi) * pnt);
    }
  };

  // *********** STANDARD KERNELS DEFINITIONS **********************
  /** LaplaceSLkernel is the kernel for the single layer potential of
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM, int COMPS=1> class LaplaceSLKernel;

  /** LaplaceDLkernel is the kernel for the double layer potential of
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM, int COMPS=1> class LaplaceDLKernel;

  /** HelmholtzSLkernel is the kernel for the double layer potential of the
      Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,. $ */
  template <int DIM, int COMPS=1, typename T_Kappa=double> class HelmholtzSLKernel;

  /** HelmholtzDLkernel is the kernel for the double layer potential of
      the Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,.$ */
  template <int DIM, int COMPS=1,typename T_Kappa=double> class HelmholtzDLKernel;

  /** CombinedFieldKernel is a kernel for the combined field integral equation
      is considered for the Helmholtz equation. */
  template <int DIM, int COMPS=1,typename T_Kappa = double> class CombinedFieldKernel;

  template <int D> class MaxwellSLKernel;

  template <int D, typename T_Kappa=double> class MaxwellDLKernel;

  /*
    Dissertation Guenther Of
    "BETI–Gebietszerlegungsmethoden
    mit schnellen Randelementverfahren
    und Anwendungen"
    page 85
   */
  template <int D> class LameSLKernel;

  // *********** DIFF KERNELS **********************
  

  template <int DIM, int COMPS=1> class DiffLaplaceSLKernel;
  template<int COMPS>
  class DiffLaplaceSLKernel<3, COMPS> : public BaseKernel
  {
  public:
    using source_type = Charges<COMPS, double>;
    using target_type = GradientEval<COMPS, double>;

    source_type source;
    target_type target;

    Array<KernelTerm> terms;
    DiffLaplaceSLKernel()
    {
      for (size_t c = 0; c < COMPS; c++)
        for (size_t i = 0; i < 3; i++)
          terms += KernelTerm{1.0, i, c, 3*c+i};
    };
    typedef double value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "DiffLaplaceSL"; }
    static auto Shape() { return IVec<2>(3*COMPS,COMPS); }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto xy = x-y;
      auto kern = 1.0 *  (4 * M_PI * norm*norm*norm) * xy;
      return kern;
    }
  };


  // grad_x G = exp(i*kappa*r) * (i*kappa*r - 1) * (x - y) / (4*pi*r^3)
  template <int DIM, int COMPS=1, typename T_Kappa=double> class DiffHelmholtzSLKernel;
  template<int COMPS, typename T_Kappa>
  class DiffHelmholtzSLKernel<3, COMPS, T_Kappa> : public BaseKernel
  {
      T_Kappa kappa;
  public:
    using source_type = Charges<COMPS, Complex, T_Kappa>;
    using target_type = GradientEval<COMPS, Complex, T_Kappa>;

    source_type source;
    target_type target;

    Array<KernelTerm> terms;
    DiffHelmholtzSLKernel(T_Kappa _kappa) : kappa(_kappa), source(_kappa), target(_kappa)
    {
      for (size_t c = 0; c < COMPS; c++)
        for (size_t i = 0; i < 3; i++)
          terms += KernelTerm{1.0, i, c, 3*c+i};
    };
    typedef Complex value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "DiffHelmholtzSL"; }
    static auto Shape() { return IVec<2>(3*COMPS,COMPS); }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto xy = x-y;
      auto kern = exp(kappa*Complex(0,1)*norm) / (4 * M_PI * norm*norm*norm)
        * (kappa*Complex(0,1)*norm - Complex(1,0)*T(1.)) * xy;
      return kern;
    }
  };

  // *********** STANDARD KERNELS **********************

  /** LaplaceSLkernel in 3D reads 
      $$ G(x-y) = \frac{1}{4\,\pi \, | x-y| }, \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS>
  class LaplaceSLKernel<3, COMPS> : public BaseKernel
  {
  public:
    using source_type = Charges<COMPS, double>;
    using target_type = Charges<COMPS, double>;

    source_type source;
    target_type target;

    LaplaceSLKernel()
    {
      for (size_t i = 0; i < COMPS; i++)
        terms += {1.0, 0, i, i};
    };
    typedef double value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "LaplaceSL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      return Vec<1,T> (1.0 / (4 * M_PI * norm));
    }

    Array<KernelTerm> terms;

    auto GetDifferentiatedKernel(const string &name) const {
      if (name == "grad")
        return DiffLaplaceSLKernel<3,COMPS>();
      else
        throw Exception("don't know how to apply diffop "+name);
    }
  };

  /** LaplaceDLkernel in 3D reads 
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, 
          \frac{ \langle n(y), x-y\rangle }{ | x-y|^3 }, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS>
  class LaplaceDLKernel<3, COMPS> : public BaseKernel
  {
  public:
    using source_type = Dipoles<COMPS, double>;
    using target_type = Charges<COMPS, double>;

    source_type source;
    target_type target;

    LaplaceDLKernel()
    {
      for (size_t i = 0; i < COMPS; i++)
        terms += {1.0, 0, i, i};
    };
    typedef double value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "LaplaceDL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      return Vec<1,T> (nxy / (4 * M_PI * norm*norm*norm));
    }

    Array<KernelTerm> terms;
  };



  /** HelmholtzSLkernel in 3D reads 
      $$ G(x-y) = \frac{1 }{4\,\pi} \,\frac{e^{i\,\kappa \, |x-y| }{|x-y|} \, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS, typename T_Kappa>
  class HelmholtzSLKernel<3, COMPS,T_Kappa> : public BaseKernel
  {
    T_Kappa kappa;
  public:
    using source_type = Charges<COMPS, Complex, T_Kappa>;
    using target_type = Charges<COMPS, Complex, T_Kappa>;

    source_type source;
    target_type target;

    typedef Complex value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "HelmholtzSL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }

    /** Construction of the kernel specifies the wavenumber $\kappa$. */
    HelmholtzSLKernel (T_Kappa _kappa) : kappa(_kappa), source(_kappa), target(_kappa) {
      for (size_t i = 0; i < COMPS; i++)
        terms += {1.0, 0, i, i};
    }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(kappa*Complex(0,1)*norm) / (4 * M_PI * norm);
      return Vec<1,decltype(kern)> (kern);
    }
    T_Kappa GetKappa() const { return kappa; }
    Array<KernelTerm> terms;

    auto GetDifferentiatedKernel(const string &name) const {
      if constexpr (COMPS == 3)
        {
          using diff_t =
            std::variant<DiffHelmholtzSLKernel<3,3,T_Kappa>,
                         MaxwellDLKernel<3,T_Kappa>>;

          if (name == "grad")
            return diff_t(DiffHelmholtzSLKernel<3,3,T_Kappa>(kappa));
          if (name == "curl")
            return diff_t(MaxwellDLKernel<3,T_Kappa>(kappa));
        }
      else
        {
          if (name == "grad")
            return DiffHelmholtzSLKernel<3,COMPS,T_Kappa>(kappa);
        }

      throw Exception("don't know how to apply diffop "+name);
    }
  };


  /** HelmholtzDLkernel in 3D reads
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \langle n(y), x-y\rangle \cdot \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS, typename T_Kappa>
  class HelmholtzDLKernel<3,COMPS,T_Kappa> : public BaseKernel
  {
    T_Kappa kappa;
  public:
    using source_type = Dipoles<COMPS, Complex, T_Kappa>;
    using target_type = Charges<COMPS, Complex, T_Kappa>;

    source_type source;
    target_type target;

    typedef Complex value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "HelmholtzDL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }

    HelmholtzDLKernel (T_Kappa _kappa) : kappa(_kappa), source(_kappa), target(_kappa) {
      for (size_t i = 0; i < COMPS; i++)
        terms += {1.0, 0, i, i};
    }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      auto kern = exp(kappa*Complex(0,1)*norm) / (4 * M_PI * norm*norm*norm)
        * nxy * (Complex(1,0)*T(1.) - kappa*Complex(0,1)*norm);
      return Vec<1,decltype(kern)> (kern);
    }
    T_Kappa GetKappa() const { return kappa; }
    Array<KernelTerm> terms;
  };


  template <int DIM> class HelmholtzHSKernel;
  
  template<>
  class HelmholtzHSKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    using source_type = HelmholtzHSSource;
    using target_type = HelmholtzHSTarget;

    source_type source;
    target_type target;

    typedef Complex value_type;
    static string Name() { return "HelmholtzHS"; }
    static auto Shape() { return IVec<2>(4,4); }

    HelmholtzHSKernel (double _kappa) : kappa(_kappa), source(_kappa), target(_kappa) { }
    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxny = InnerProduct(nx, ny);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      auto kernnxny = -kappa * kappa * kern * nxny;
      // return kern;
      return Vec<2,decltype(kern)> ({kern, kernnxny});
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms =
      {
        KernelTerm{1.0, 0, 0, 0},
        KernelTerm{1.0, 0, 1, 1},
        KernelTerm{1.0, 0, 2, 2},
    	KernelTerm{1.0, 1, 3, 3},
      };
  };


  /** CombinedFieldKernel in 3D reads
      $$ G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \left( \langle n_y, x-y\rangle (1- i\,\kappa\, | x-y|) - i\,\kappa\,|x-y|^2 \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS, typename T_Kappa>
  class CombinedFieldKernel<3,COMPS,T_Kappa> : public BaseKernel
  {
    T_Kappa kappa;
  public:
    using source_type = ChargeDipoles<COMPS, T_Kappa>;
    using target_type = DirectEval<typename source_type::mp_type, Complex, T_Kappa>;

    source_type source;
    target_type target;

    typedef Complex value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "Helmholtz Combined Field"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }

    CombinedFieldKernel (T_Kappa _kappa) : kappa(_kappa), source(_kappa), target(_kappa) {
      for (size_t i = 0; i < COMPS; i++)
        terms += {1.0, 0, i, i};
    }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      auto kern = exp(kappa*Complex(0,1)*norm) / (4 * M_PI * norm*norm*norm)
        * ( nxy * (Complex(1,0)*T(1.) - kappa*Complex(0,1)*norm)  - kappa*Complex(0,1)*norm*norm);
      return Vec<1,decltype(kern)> (kern);
    }
    T_Kappa GetKappa() const { return kappa; }
    Array<KernelTerm> terms;
  };



  template<>
  class MaxwellSLKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    using source_type = MaxwellSLSource;
    using target_type = DirectEval<Vec<4,Complex>, Complex>;

    source_type source;
    target_type target;

    typedef Complex value_type;
    static string Name() { return "MaxwellSL"; }
    static auto Shape() { return IVec<2>(4,4); }

    MaxwellSLKernel (const MaxwellSLKernel&) = default;
    MaxwellSLKernel (MaxwellSLKernel&&) = default;
    MaxwellSLKernel (double _kappa) : kappa(_kappa), source(_kappa), target(_kappa)
    {
      for (size_t i = 0; i < 3; i++)
        terms += KernelTerm { kappa, 0, i, i };
      terms += KernelTerm { -1.0/kappa, 0, 3, 3 };
    }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms;
  };


  // https://weggler.github.io/ngbem/short_and_sweet/Maxwell_Formulations.html
  /** MaxwellDLkernel for 3D in matrix representation reads
      $$  \left( \begin{array}{ccc} 0 & -\frac{\partial G_\kappa(x-y)}{\partial x_3} & \frac{\partial G_\kappa(x-y)}{\partial x_2} \\ \frac{\partial G_\kappa(x-y)}{\partial x_3} & 0 & -\frac{\partial G_\kappa(x-y)}{\partial x_1} \\ -\frac{\partial G_\kappa(x-y)}{\partial x_2} & \frac{\partial G_\kappa(x-y)}{\partial x_1} & 0 \end{array}\right)\,,$$ with 
   $$ G_\kappa(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \langle n(y), x-y\rangle \cdot \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<typename T_Kappa>
  class MaxwellDLKernel<3,T_Kappa> : public BaseKernel
  {
    T_Kappa kappa;
  public:
    using source_type = MaxwellCurlDipoles<T_Kappa>;
    using target_type = MaxwellCurlDipoles<T_Kappa>;

    source_type source;
    target_type target;

    typedef Complex value_type;
    static string Name() { return "MaxwellDL"; }
    static auto Shape() { return IVec<2>(3,3); }

    MaxwellDLKernel (T_Kappa _kappa) : kappa(_kappa), source(_kappa), target(_kappa) { }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(kappa*Complex(0,1)*norm) / (4 * M_PI * norm*norm*norm)
        * (kappa*Complex(0,1)*norm - Complex(1,0)*T(1.)) * (x-y);
      return kern;
    }
    T_Kappa GetKappa() const { return kappa; }
    Array<KernelTerm> terms =
      {
        KernelTerm{ 1.0, 0, 1, 2},  // factor, comp, trial, test
        KernelTerm{-1.0, 0, 2, 1},
        KernelTerm{ 1.0, 1, 2, 0},
        KernelTerm{-1.0, 1, 0, 2},
        KernelTerm{ 1.0, 2, 0, 1},
        KernelTerm{-1.0, 2, 1, 0},
      };
  };


  template<>
  class LameSLKernel<3> : public BaseKernel
  {
    double E, nu;
    double alpha;
  public:
    using source_type = LameSource;
    using target_type = LameTarget;

    source_type source;
    target_type target;

    typedef double value_type;
    static string Name() { return "LameSL"; }
    static auto Shape() { return IVec<2>(3,3); }

    LameSLKernel (const LameSLKernel&) = default;
    LameSLKernel (LameSLKernel&&) = default;
    LameSLKernel (double _E, double _nu) : E(_E), nu(_nu), target(_nu, (1+_nu)/((1-_nu)*2*_E))
    {
      alpha = (1+nu)/((1-nu)*2*E);

      terms += { 3-4*nu, 0, 0, 0 };
      terms += { 3-4*nu, 0, 1, 1 };
      terms += { 3-4*nu, 0, 2, 2 };

      terms += { 1, 1, 0, 0 };
      terms += { 1, 2, 0, 1 };
      terms += { 1, 2, 1, 0 };

      terms += { 1, 3, 0, 2 };
      terms += { 1, 3, 2, 0 };
      terms += { 1, 4, 1, 1 };

      terms += { 1, 5, 1, 2 };
      terms += { 1, 5, 2, 1 };
      terms += { 1, 6, 2, 2 };
    }

    Array<KernelTerm> terms;

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto lapkern = alpha / (4 * M_PI * norm);

      return Vec<7,T> { lapkern,
                        (x(0)-y(0))*(x(0)-y(0))/sqr(norm) * lapkern,
                        (x(0)-y(0))*(x(1)-y(1))/sqr(norm) * lapkern,
                        (x(0)-y(0))*(x(2)-y(2))/sqr(norm) * lapkern,
                        (x(1)-y(1))*(x(1)-y(1))/sqr(norm) * lapkern,
                        (x(1)-y(1))*(x(2)-y(2))/sqr(norm) * lapkern,
                        (x(2)-y(2))*(x(2)-y(2))/sqr(norm) * lapkern
      };
    }
  };

}


#endif
