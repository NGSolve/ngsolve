#ifndef KERNELS_hpp
#define KERNELS_hpp

#include "mptools.hpp"
#include "integralkernel.hpp"
#include <cassert>
#include <type_traits>
#include <variant>


//  ****************************  The kernels **********************************

namespace ngsbem
{


  class BaseKernel
  {
  public:
    static constexpr AnalyticTriangleFormula analytic_triangle_formula =
      AnalyticTriangleFormula::none;

    void GetDifferentiatedKernel(const string &name) const {
      throw Exception("GetDifferentiatedKernel not overloaded");
    }
  };

  template <typename T>
  using fmm_complex_type = std::conditional_t<std::is_same_v<T,float> || std::is_same_v<T,Complex32>, Complex32, Complex>;


  // *********** FMM Source/Target Interface **********************

  template <typename T_MP, typename T_VALUE, typename T_Kappa = double>
  class NGS_DLL_HEADER FMMInterface : public BaseFMMInterface
  {
  public:
    using mp_type = T_MP;
    using value_type = T_VALUE;
    T_Kappa kappa;

    FMMInterface (T_Kappa _kappa = 1e-16) : kappa(_kappa) {}

    int ExpansionComponents () const override { return VecLength<T_MP>; }
    shared_ptr<BaseSingularMLExpansion> CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const override;
    shared_ptr<BaseRegularMLExpansion> CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const override;

  protected:
    // These interfaces create their own expansion type; check it in debug builds.
    static SingularMLExpansion<mp_type,T_Kappa> & Cast(BaseSingularMLExpansion & mp)
    {
      assert((dynamic_cast<SingularMLExpansion<mp_type,T_Kappa>*>(&mp)));
      return static_cast<SingularMLExpansion<mp_type,T_Kappa>&>(mp);
    }

    static RegularMLExpansion<mp_type,T_Kappa> & Cast(BaseRegularMLExpansion & mp)
    {
      assert((dynamic_cast<RegularMLExpansion<mp_type,T_Kappa>*>(&mp)));
      return static_cast<RegularMLExpansion<mp_type,T_Kappa>&>(mp);
    }

    // the device's values as they are: val(i) += dev(i), its real part for a real val
    void AddDeviceTargetValue (FlatVector<Complex> dev, Vec<3> p, Vec<3> nv, BareSliceVector<T_VALUE> val) const override
    {
      for (size_t i = 0; i < dev.Size(); i++)
        if constexpr (!IsComplex<T_VALUE>())
          val(i) += dev(i).real();
        else
          val(i) += dev(i);
    }
  };

  template <int COMPS, typename T_VALUE, typename T_FMM = T_VALUE, typename T_Kappa = double>
  class Charges : public FMMInterface<Vec<COMPS, fmm_complex_type<T_FMM>>, T_VALUE, T_Kappa>
  {
  public:
    using value_type = T_VALUE;
    using coefficient_type = fmm_complex_type<T_FMM>;
    using mp_type = Vec<COMPS, coefficient_type>;
    using Base = FMMInterface<mp_type, value_type, T_Kappa>;
    static constexpr bool needs_normal = false;
    bool NeedsNormal() const override { return needs_normal; }
    Charges (T_Kappa kappa = 1e-16) : Base(kappa) {}

    void AddSource (BaseSingularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      auto & mp = this->Cast(amp);
      mp.AddCharge (pnt, val);
    }

    int DeviceSourceKind () const override { return 1; }

    void DeviceSourceValues (BareSliceVector<value_type> val, Vec<3> p, FlatVector<Complex> charge,
                             FlatVector<Complex> dipole, FlatVector<Complex> dipvec) const override
    {
      charge = mp_type(val);
    }

    int DeviceTargetMode () const override { return 1; }

    void EvaluateMP (BaseRegularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      auto & mp = this->Cast(amp);
      if constexpr (!IsComplex<value_type>())
        val = Real(mp.Evaluate (pnt));
      else
        val = mp.Evaluate (pnt);
    }
  };

  template <int COMPS, typename T_VALUE, typename T_FMM = T_VALUE, typename T_Kappa = double>
  class Dipoles : public FMMInterface<Vec<COMPS, fmm_complex_type<T_FMM>>, T_VALUE, T_Kappa>
  {
  public:
    using value_type = T_VALUE;
    using coefficient_type = fmm_complex_type<T_FMM>;
    using mp_type = Vec<COMPS, coefficient_type>;
    using Base = FMMInterface<mp_type, value_type, T_Kappa>;
    static constexpr bool needs_normal = true;
    bool NeedsNormal() const override { return needs_normal; }
    Dipoles (T_Kappa kappa = 1e-16) : Base(kappa) {}

    void AddSource (BaseSingularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      auto & mp = this->Cast(amp);
      mp.AddDipole (pnt, -nv, val);
    }

    int DeviceSourceKind () const override { return 2; }

    void DeviceSourceValues (BareSliceVector<value_type> val, Vec<3> p, FlatVector<Complex> charge,
                             FlatVector<Complex> dipole, FlatVector<Complex> dipvec) const override
    {
      dipole = mp_type(val);
    }

    void EvaluateMP (BaseRegularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      auto & mp = this->Cast(amp);
      if constexpr (!IsComplex<value_type>())
        val = Real(mp.EvaluateDirectionalDerivative (pnt, nv));
      else
        val = mp.EvaluateDirectionalDerivative (pnt, nv);
    }

    int DeviceTargetMode () const override { return 2; }

    // the derivative along nv from the gradient
    void AddDeviceTargetValue (FlatVector<Complex> dev, Vec<3> p, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      for (int i = 0; i < COMPS; i++)
        {
          Complex d = nv(0)*dev(3*i) + nv(1)*dev(3*i+1) + nv(2)*dev(3*i+2);
          if constexpr (!IsComplex<value_type>())
            val(i) += d.real();
          else
            val(i) += d;
        }
    }
  };

  // Evaluates gradient: 3 directional derivatives along unit vectors
  template <int COMPS, typename T_VALUE, typename T_FMM = T_VALUE, typename T_Kappa = double>
  class GradientEval : public FMMInterface<Vec<COMPS, fmm_complex_type<T_FMM>>, T_VALUE, T_Kappa>
  {
  public:
    using value_type = T_VALUE;
    using coefficient_type = fmm_complex_type<T_FMM>;
    using mp_type = Vec<COMPS, coefficient_type>;
    using Base = FMMInterface<mp_type, value_type, T_Kappa>;
    static constexpr bool needs_normal = false;
    bool NeedsNormal() const override { return needs_normal; }
    GradientEval (T_Kappa kappa = 1e-16) : Base(kappa) {}

    void EvaluateMP (BaseRegularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      auto & mp = this->Cast(amp);
      for (int i = 0; i < 3; i++)
        {
          Vec<3> ei = 0;
          ei(i) = 1;
          auto deri = mp.EvaluateDirectionalDerivative (pnt, ei);
          if constexpr (!IsComplex<value_type>())
            val.Slice(i,3) = Real(deri);
          else
            val.Slice(i,3) = deri;
        }
    }

    int DeviceTargetMode () const override { return 2; }
  };

  // Direct evaluation: val = mp.Evaluate(pnt)
  template <typename mp_type, typename value_type, typename T_Kappa = double>
  class DirectEval : public FMMInterface<mp_type, value_type, T_Kappa>
  {
    using Base = FMMInterface<mp_type, value_type, T_Kappa>;
  public:
    static constexpr bool needs_normal = false;
    bool NeedsNormal() const override { return needs_normal; }
    DirectEval (T_Kappa kappa = 1e-16) : Base(kappa) {}

    void EvaluateMP (BaseRegularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<value_type> val) const override
    {
      auto & mp = this->Cast(amp);
      if constexpr (std::is_same_v<mp_type, Complex> || std::is_same_v<mp_type, Complex32>)
        val(0) = mp.Evaluate (pnt);
      else
        val = mp.Evaluate (pnt);
    }

    int DeviceTargetMode () const override { return 1; }
  };

  // Combined charge + dipole source (CombinedField)
  template <int COMPS, typename T_Kappa = double, typename T_FMM = Complex>
  class ChargeDipoles : public FMMInterface<Vec<COMPS, fmm_complex_type<T_FMM>>, Complex, T_Kappa>
  {
  public:
    using coefficient_type = fmm_complex_type<T_FMM>;
    using mp_type = Vec<COMPS, coefficient_type>;
    using Base = FMMInterface<mp_type, Complex, T_Kappa>;
    static constexpr bool needs_normal = true;
    bool NeedsNormal() const override { return needs_normal; }
    ChargeDipoles (T_Kappa kappa) : Base(kappa) {}

    void AddSource (BaseSingularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const override
    {
      auto & mp = this->Cast(amp);
      Vec<COMPS,coefficient_type> fmm_val = val;
      coefficient_type imag(0,1);
      mp.AddChargeDipole (pnt, -coefficient_type(this->kappa) * imag*fmm_val, -nv, fmm_val);
    }

    int DeviceSourceKind () const override { return 3; }

    void DeviceSourceValues (BareSliceVector<Complex> val, Vec<3> p, FlatVector<Complex> charge,
                             FlatVector<Complex> dipole, FlatVector<Complex> dipvec) const override
    {
      Vec<COMPS,coefficient_type> fmm_val = val;
      coefficient_type imag(0,1);
      Vec<COMPS,coefficient_type> c = -coefficient_type(this->kappa) * imag*fmm_val;
      for (int i = 0; i < COMPS; i++) { charge(i) = c(i); dipole(i) = fmm_val(i); }
    }
  };

  // MaxwellDL source/target: cross-product dipoles + direct eval (self-symmetric)
  template <typename T_Kappa = double, typename T_FMM = Complex>
  class MaxwellCurlDipoles : public FMMInterface<Vec<3,fmm_complex_type<T_FMM>>, Complex, T_Kappa>
  {
    using coefficient_type = fmm_complex_type<T_FMM>;
    using mp_type = Vec<3,fmm_complex_type<T_FMM>>;
    using Base = FMMInterface<mp_type, Complex, T_Kappa>;
  public:
    static constexpr bool needs_normal = false;
    bool NeedsNormal() const override { return needs_normal; }
    MaxwellCurlDipoles (T_Kappa kappa) : Base(kappa) {}

    void AddSource (BaseSingularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const override
    {
      auto & mp = this->Cast(amp);
      Vec<3,coefficient_type> n_cross_m = val.Range(0, 3);
      for (int k = 0; k < 3; k++)
        {
          using real_type = typename scal_traits<coefficient_type>::TSCAL_REAL;
          Vec<3> ek{0.0}; ek(k) = 1;
          Vec<3,real_type> n_cross_m_real = Real(n_cross_m);
          Vec<3,real_type> n_cross_m_imag = Imag(n_cross_m);
          mp.AddDipole(pnt, Cross(n_cross_m_real, ek), ek);
          mp.AddDipole(pnt, Cross(n_cross_m_imag, ek), Complex(0,1)*ek);
        }
    }

    int DeviceSourceKind () const override { return 4; }

    int DeviceSourceEntries () const override { return 6; }

    // component k: the dipole vector v x e_k of the six dipoles AddSource adds
    void DeviceSourceValues (BareSliceVector<Complex> val, Vec<3> p, FlatVector<Complex> charge,
                             FlatVector<Complex> dipole, FlatVector<Complex> dipvec) const override
    {
      Vec<3,coefficient_type> v = val.Range(0, 3);
      coefficient_type z(0.0);
      Vec<3,coefficient_type> d[3] = { { z, v(2), -v(1) }, { -v(2), z, v(0) }, { v(1), -v(0), z } };
      for (int k = 0; k < 3; k++)
        for (int j = 0; j < 3; j++)
          dipvec(3*k+j) = d[k](j);
    }

    void EvaluateMP (BaseRegularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const override
    {
      auto & mp = this->Cast(amp);
      val = mp.Evaluate (pnt);
    }

    int DeviceTargetMode () const override { return 1; }
  };

  // Lame source: charges + Jacobi matrix dipoles
  template <typename T_FMM = double>
  class LameSource : public FMMInterface<Vec<6,fmm_complex_type<T_FMM>>, double>
  {
    using coefficient_type = fmm_complex_type<T_FMM>;
    using mp_type = Vec<6,fmm_complex_type<T_FMM>>;
    using Base = FMMInterface<mp_type, double>;
  public:
    static constexpr bool needs_normal = false;
    bool NeedsNormal() const override { return needs_normal; }
    LameSource () : Base(1e-16) {}

    void AddSource (BaseSingularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const override
    {
      auto & mp = this->Cast(amp);
      Vec<6,coefficient_type> charge = coefficient_type(0.0);
      charge.Range(0,3) = val;
      mp.AddCharge(pnt, charge);

      Mat<3,3> jacobi = OuterProduct(pnt, val) + InnerProduct(pnt, val) * Id<3>();

      for (int k = 0; k < 3; k++)
        {
          Vec<6,coefficient_type> dipole_charge = coefficient_type(0.0);
          dipole_charge.Range(3,6) = jacobi.Col(k);
          auto ek = UnitVec<3>(k);
          mp.AddDipole(pnt, -ek, dipole_charge);
        }
    }

    int DeviceSourceKind () const override { return 5; }

    int DeviceSourceEntries () const override { return 4; }

    // the charges on comps 0..2, and the dipole vectors of the three dipoles on comps 3..5: -(p_i val + (p.val) e_i)
    void DeviceSourceValues (BareSliceVector<double> val, Vec<3> p, FlatVector<Complex> charge,
                             FlatVector<Complex> dipole, FlatVector<Complex> dipvec) const override
    {
      Mat<3,3> jacobi = OuterProduct(p, val) + InnerProduct(p, val) * Id<3>();
      for (int i = 0; i < 3; i++)
        {
          charge(i) = coefficient_type(val(i));
          charge(3+i) = coefficient_type(0.0);
          for (int k = 0; k < 3; k++)
            {
              dipvec(3*i+k) = coefficient_type(0.0);
              dipvec(3*(3+i)+k) = -coefficient_type(jacobi(i,k));
            }
        }
    }
  };

  // Lame target: multi-step eval with coefficients
  template <typename T_FMM = double>
  class LameTarget : public FMMInterface<Vec<6,fmm_complex_type<T_FMM>>, double>
  {
    double nu, alpha;
    using mp_type = Vec<6,fmm_complex_type<T_FMM>>;
    using Base = FMMInterface<mp_type, double>;
  public:
    static constexpr bool needs_normal = false;
    bool NeedsNormal() const override { return needs_normal; }
    LameTarget (double _nu, double _alpha) : Base(1e-16), nu(_nu), alpha(_alpha) {}

    // val += the displacement at pnt from the expansion's values mpval and gradients jacobi(i,k) = d_k mpval(i)
    void AddDisplacement (Vec<6> mpval, const Mat<3,3> & jacobi, Vec<3> pnt, BareSliceVector<double> val) const
    {
      val += (3-4*nu)*alpha * mpval.Range(0,3);
      val -= alpha/2 * mpval.Range(3,6);
      val -= alpha/2 * ( Trans(jacobi) * pnt + Trace(jacobi) * pnt);
    }

    void EvaluateMP (BaseRegularMLExpansion & amp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const override
    {
      auto & mp = this->Cast(amp);
      Vec<6> mpval = Real(mp.Evaluate (pnt));

      Mat<3,3> jacobi = 0.0;
      for (int k = 0; k < 3; k++)
        {
          auto ek = UnitVec<3>(k);
          jacobi.Col(k) = Real(mp.EvaluateDirectionalDerivative(pnt, ek).Range(0,3));
        }

      val.Range(0,3) = 0;
      AddDisplacement (mpval, jacobi, pnt, val);
    }

    int DeviceTargetMode () const override { return 3; }

    // EvaluateMP of the device's values and gradients dev[4c], dev[4c+1+k]
    void AddDeviceTargetValue (FlatVector<Complex> dev, Vec<3> p, Vec<3> nv, BareSliceVector<double> val) const override
    {
      Vec<6> mpval;
      Mat<3,3> jacobi;
      for (int c = 0; c < 6; c++)
        mpval(c) = dev(4*c).real();
      for (int i = 0; i < 3; i++)
        for (int k = 0; k < 3; k++)
          jacobi(i,k) = dev(4*i+1+k).real();
      AddDisplacement (mpval, jacobi, p, val);
    }
  };

  // *********** STANDARD KERNELS DEFINITIONS **********************
  /** LaplaceSLkernel is the kernel for the single layer potential of
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM, int COMPS=1, typename T_VAL=double, typename T_FMM=T_VAL> class LaplaceSLKernel;

  /** LaplaceDLkernel is the kernel for the double layer potential of
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM, int COMPS=1, typename T_VAL=double, typename T_FMM=T_VAL> class LaplaceDLKernel;

  /** HelmholtzSLkernel is the kernel for the double layer potential of the
      Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,. $ */
  template <int DIM, int COMPS=1, typename T_Kappa=double, typename T_FMM=Complex> class HelmholtzSLKernel;

  /** HelmholtzDLkernel is the kernel for the double layer potential of
      the Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,.$ */
  template <int DIM, int COMPS=1,typename T_Kappa=double, typename T_FMM=Complex> class HelmholtzDLKernel;

  /** CombinedFieldKernel is a kernel for the combined field integral equation
      is considered for the Helmholtz equation. */
  template <int DIM, int COMPS=1,typename T_Kappa = double, typename T_FMM = Complex> class CombinedFieldKernel;

  template <int D, typename T_Kappa=double, typename T_FMM=Complex> class MaxwellDLKernel;

  /*
    Dissertation Guenther Of
    "BETI–Gebietszerlegungsmethoden
    mit schnellen Randelementverfahren
    und Anwendungen"
    page 85
   */
  template <int D, typename T_FMM=double> class LameSLKernel;

  // *********** DIFF KERNELS **********************
  

  template <int DIM, int COMPS=1, typename T_VAL=double, typename T_FMM=T_VAL> class DiffLaplaceSLKernel;
  template<int COMPS, typename T_VAL, typename T_FMM>
  class DiffLaplaceSLKernel<3, COMPS, T_VAL, T_FMM> : public BaseKernel
  {
  public:
    using source_type = Charges<COMPS, T_VAL, T_FMM>;
    using target_type = GradientEval<COMPS, T_VAL, T_FMM>;

    source_type source;
    target_type target;

    Array<KernelTerm> terms;
    DiffLaplaceSLKernel()
    {
      for (size_t c = 0; c < COMPS; c++)
        for (size_t i = 0; i < 3; i++)
          terms += KernelTerm{1.0, i, c, 3*c+i};
    };
    typedef T_VAL value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "DiffLaplaceSL"; }
    static auto Shape() { return IVec<2>(3*COMPS,COMPS); }
    static constexpr AnalyticTriangleFormula analytic_triangle_formula =
      AnalyticTriangleFormula::laplace_grad_sl;

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto xy = x-y;
      auto kern = -1.0 / (4 * M_PI * norm*norm*norm) * xy;
      return kern;
    }
  };


  // grad_x G = exp(i*kappa*r) * (i*kappa*r - 1) * (x - y) / (4*pi*r^3)
  template <int DIM, int COMPS=1, typename T_Kappa=double, typename T_FMM=Complex> class DiffHelmholtzSLKernel;
  template<int COMPS, typename T_Kappa, typename T_FMM>
  class DiffHelmholtzSLKernel<3, COMPS, T_Kappa, T_FMM> : public BaseKernel
  {
      T_Kappa kappa;
  public:
    using source_type = Charges<COMPS, Complex, T_FMM, T_Kappa>;
    using target_type = GradientEval<COMPS, Complex, T_FMM, T_Kappa>;

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
    static constexpr AnalyticTriangleFormula analytic_triangle_formula =
      AnalyticTriangleFormula::laplace_grad_sl;

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

  /* The gradient of Helmholtz double layer kernel in 3D reads
      $$ \nabla_x \frac{\partial G_\kappa(x-y)}{\partial n_y}
      = \frac{e^{i\kappa r}}{4\pi r^5} \left[ (1-i\kappa r)r^2 n(y) + \langle n(y),x-y\rangle \left(\kappa^2r^2-3(1-i\kappa r)\right)(x-y) \right],
      \qquad r=|x-y|,\quad x\neq y. $$
  */
  template <int DIM, int COMPS=1, typename T_Kappa=double, typename T_FMM=Complex> class DiffHelmholtzDLKernel;
  template<int COMPS, typename T_Kappa, typename T_FMM>
  class DiffHelmholtzDLKernel<3, COMPS, T_Kappa, T_FMM> : public BaseKernel
  {
    T_Kappa kappa;
  public:
    using source_type = Dipoles<COMPS, Complex, T_FMM, T_Kappa>;
    using target_type = GradientEval<COMPS, Complex, T_FMM, T_Kappa>;

    source_type source;
    target_type target;

    Array<KernelTerm> terms;
    DiffHelmholtzDLKernel(T_Kappa _kappa)
      : kappa(_kappa), source(_kappa), target(_kappa)
    {
      for (size_t c = 0; c < COMPS; c++)
        for (size_t i = 0; i < 3; i++)
          terms += KernelTerm{1.0, i, c, 3*c+i};
    }
    typedef Complex value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "DiffHelmholtzDL"; }
    static auto Shape() { return IVec<2>(3*COMPS,COMPS); }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto xy = x-y;
      T nxy = InnerProduct(ny, xy);
      auto expikr = exp(kappa*Complex(0,1)*norm);
      auto one_minus_ikr = Complex(1,0)*T(1.) - kappa*Complex(0,1)*norm;
      auto first_factor = expikr * one_minus_ikr / (4 * M_PI * norm*norm*norm);
      auto second_factor = expikr *
        (Complex(1,0)*kappa*kappa*T(1.)*norm*norm - 3.0*one_minus_ikr) /
        (4 * M_PI * norm*norm*norm*norm*norm);
      return first_factor*ny + nxy*second_factor*xy;
    }
  };

  // *********** STANDARD KERNELS **********************

  /** LaplaceSLkernel in 3D reads 
      $$ G(x-y) = \frac{1}{4\,\pi \, | x-y| }, \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS, typename T_VAL, typename T_FMM>
  class LaplaceSLKernel<3, COMPS, T_VAL, T_FMM> : public BaseKernel
  {
  public:
    using source_type = Charges<COMPS, T_VAL, T_FMM>;
    using target_type = Charges<COMPS, T_VAL, T_FMM>;

    source_type source;
    target_type target;

    LaplaceSLKernel()
    {
      for (size_t i = 0; i < COMPS; i++)
        terms += {1.0, 0, i, i};
    };
    typedef T_VAL value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "LaplaceSL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }
    static constexpr AnalyticTriangleFormula analytic_triangle_formula =
      AnalyticTriangleFormula::laplace_sl;

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      return Vec<1,T> (1.0 / (4 * M_PI * norm));
    }

    Array<KernelTerm> terms;

    auto GetDifferentiatedKernel(const string &name) const {
      if (name == "grad")
        return DiffLaplaceSLKernel<3,COMPS,T_VAL,T_FMM>();
      else
        throw Exception("don't know how to apply diffop "+name);
    }
  };

  /** LaplaceDLkernel in 3D reads 
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, 
          \frac{ \langle n(y), x-y\rangle }{ | x-y|^3 }, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS, typename T_VAL, typename T_FMM>
  class LaplaceDLKernel<3, COMPS, T_VAL, T_FMM> : public BaseKernel
  {
  public:
    using source_type = Dipoles<COMPS, T_VAL, T_FMM>;
    using target_type = Charges<COMPS, T_VAL, T_FMM>;

    source_type source;
    target_type target;

    LaplaceDLKernel()
    {
      for (size_t i = 0; i < COMPS; i++)
        terms += {1.0, 0, i, i};
    };
    typedef T_VAL value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "LaplaceDL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }
    static constexpr AnalyticTriangleFormula analytic_triangle_formula =
      AnalyticTriangleFormula::laplace_dl;

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
  template<int COMPS, typename T_Kappa, typename T_FMM>
  class HelmholtzSLKernel<3, COMPS,T_Kappa,T_FMM> : public BaseKernel
  {
    T_Kappa kappa;
  public:
    using source_type = Charges<COMPS, Complex, T_FMM, T_Kappa>;
    using target_type = Charges<COMPS, Complex, T_FMM, T_Kappa>;

    source_type source;
    target_type target;

    typedef Complex value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "HelmholtzSL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }
    static constexpr AnalyticTriangleFormula analytic_triangle_formula =
      AnalyticTriangleFormula::laplace_sl;

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
            std::variant<DiffHelmholtzSLKernel<3,3,T_Kappa,T_FMM>,
                         MaxwellDLKernel<3,T_Kappa,T_FMM>>;

          if (name == "grad")
            return diff_t(DiffHelmholtzSLKernel<3,3,T_Kappa,T_FMM>(kappa));
          if (name == "curl")
            return diff_t(MaxwellDLKernel<3,T_Kappa,T_FMM>(kappa));
        }
      else
        {
          if (name == "grad")
            return DiffHelmholtzSLKernel<3,COMPS,T_Kappa,T_FMM>(kappa);
        }

      throw Exception("don't know how to apply diffop "+name);
    }
  };


  /** HelmholtzDLkernel in 3D reads
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \langle n(y), x-y\rangle \cdot \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS, typename T_Kappa, typename T_FMM>
  class HelmholtzDLKernel<3,COMPS,T_Kappa,T_FMM> : public BaseKernel
  {
    T_Kappa kappa;
  public:
    using source_type = Dipoles<COMPS, Complex, T_FMM, T_Kappa>;
    using target_type = Charges<COMPS, Complex, T_FMM, T_Kappa>;

    source_type source;
    target_type target;

    typedef Complex value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "HelmholtzDL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }
    static constexpr AnalyticTriangleFormula analytic_triangle_formula =
      AnalyticTriangleFormula::laplace_dl;

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

    auto GetDifferentiatedKernel(const string &name) const {
      if (name == "grad")
        return DiffHelmholtzDLKernel<3,COMPS,T_Kappa,T_FMM>(kappa);
      throw Exception("don't know how to apply diffop "+name);
    }
  };

  /** CombinedFieldKernel in 3D reads
      $$ G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \left( \langle n_y, x-y\rangle (1- i\,\kappa\, | x-y|) - i\,\kappa\,|x-y|^2 \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS, typename T_Kappa, typename T_FMM>
  class CombinedFieldKernel<3,COMPS,T_Kappa,T_FMM> : public BaseKernel
  {
    T_Kappa kappa;
  public:
    using source_type = ChargeDipoles<COMPS, T_Kappa, T_FMM>;
    using target_type = DirectEval<typename source_type::mp_type, Complex, T_Kappa>;

    source_type source;
    target_type target;

    typedef Complex value_type;
    using mp_type = typename source_type::mp_type;

    static string Name() { return "Helmholtz Combined Field"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }
    static constexpr AnalyticTriangleFormula analytic_triangle_formula =
      AnalyticTriangleFormula::helmholtz_cf;

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

  // https://weggler.github.io/ngbem/short_and_sweet/Maxwell_Formulations.html
  /** MaxwellDLkernel for 3D in matrix representation reads
      $$  \left( \begin{array}{ccc} 0 & -\frac{\partial G_\kappa(x-y)}{\partial x_3} & \frac{\partial G_\kappa(x-y)}{\partial x_2} \\ \frac{\partial G_\kappa(x-y)}{\partial x_3} & 0 & -\frac{\partial G_\kappa(x-y)}{\partial x_1} \\ -\frac{\partial G_\kappa(x-y)}{\partial x_2} & \frac{\partial G_\kappa(x-y)}{\partial x_1} & 0 \end{array}\right)\,,$$ with 
   $$ G_\kappa(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \langle n(y), x-y\rangle \cdot \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<typename T_Kappa, typename T_FMM>
  class MaxwellDLKernel<3,T_Kappa,T_FMM> : public BaseKernel
  {
    T_Kappa kappa;
  public:
    using source_type = MaxwellCurlDipoles<T_Kappa,T_FMM>;
    using target_type = MaxwellCurlDipoles<T_Kappa,T_FMM>;

    source_type source;
    target_type target;

    typedef Complex value_type;
    static string Name() { return "MaxwellDL"; }
    static auto Shape() { return IVec<2>(3,3); }
    static constexpr AnalyticTriangleFormula analytic_triangle_formula =
      AnalyticTriangleFormula::laplace_grad_sl;

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


  template<typename T_FMM>
  class LameSLKernel<3,T_FMM> : public BaseKernel
  {
    double E, nu;
    double alpha;
  public:
    using source_type = LameSource<T_FMM>;
    using target_type = LameTarget<T_FMM>;

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
