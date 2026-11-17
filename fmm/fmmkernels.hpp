#ifndef FMMKERNELS_hpp
#define FMMKERNELS_hpp

#include "mptools.hpp"
#include <cassert>
#include <type_traits>

namespace ngsfmm
{
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
      // grad G x v = sum_k (e_k . grad G) (e_k x v)
      for (int k = 0; k < 3; k++)
        {
          Vec<3> ek{0.0}; ek(k) = 1;
          mp.AddDipole(pnt, ek, Cross(Vec<3,coefficient_type>(ek), n_cross_m));
        }
    }

    int DeviceSourceKind () const override { return 4; }

    int DeviceSourceEntries () const override { return 3; }

    // component k: the dipole vector v x e_k
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
          // the dipoles (LameSource) feed only comps 3..5
          jacobi.Col(k) = Real(mp.EvaluateDirectionalDerivative(pnt, ek, false).Range(0,3));
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
}

#endif
