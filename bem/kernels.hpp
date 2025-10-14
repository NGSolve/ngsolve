#ifndef KERNELS_hpp
#define KERNELS_hpp

#include "mptools.hpp"
#include <type_traits>


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
    shared_ptr<SingularMLExpansion<Complex>> CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      throw Exception("Create Multipole Expansion not implemented");
    }

    shared_ptr<RegularMLExpansion<Complex>> CreateLocalExpansion (Vec<3> c, double r) const
    {
      throw Exception("Create Local Expansion not implemented");      
    }

    template <typename TV>
    void AddSource (SingularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, const TV & val) const
    {
      throw Exception("Addsource not implemented");            
    }

    template <typename entry, typename TV>
    void AddSourceTrans(SingularMLExpansion<entry> & mp, Vec<3> pnt, Vec<3> nv, const TV & val) const
    {
      throw Exception("AddSourceTrans not implemented");
    }

    template <typename TV>
    void EvaluateMP (RegularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, const TV & val) const
    {
      throw Exception("Evaluate not implemented");            
    }

    template <typename entry, typename TV>
    void EvaluateMPTrans(RegularMLExpansion<entry> & mp, Vec<3> pnt, Vec<3> nv, const TV & val) const
    {
      throw Exception("EvaluateMPTrans not implemented");
    }
  };



  
  
  /** LaplaceSLkernel is the kernel for the single layer potential of 
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM, int COMPS=1> class LaplaceSLKernel;

  /** LaplaceSLkernel in 3D reads 
      $$ G(x-y) = \frac{1}{4\,\pi \, | x-y| }, \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS>
  class LaplaceSLKernel<3, COMPS> : public BaseKernel
  {
  public:
    LaplaceSLKernel<3,COMPS>()
    {
      for (size_t i = 0; i < COMPS; i++)
        terms += {1.0, 0, i, i};
    };
    typedef double value_type;
    using mp_type = typename std::conditional<COMPS == 1,
                                              Complex,
                                              Vec<COMPS, Complex>>::type;

    static string Name() { return "LaplaceSL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }
    
    template <typename T>        
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      // return 1.0 / (4 * M_PI * norm);   
      return Vec<1,T> (1.0 / (4 * M_PI * norm));
    }
    
    Array<KernelTerm> terms;
    
    auto CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<SingularMLExpansion<mp_type>> (c, r, 1e-16, fmm_params);
    }

    auto CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<RegularMLExpansion<mp_type>> (c, r, 1e-16, fmm_params);
    }

    void AddSource (SingularMLExpansion<mp_type> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      if constexpr (COMPS == 1)
        mp.AddCharge (pnt, val(0));
      else
        mp.AddCharge (pnt, val);
    }

    void EvaluateMP (RegularMLExpansion<mp_type> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      if constexpr (COMPS == 1)
        val(0) = Real(mp.Evaluate (pnt));
      else
        val = Real(mp.Evaluate (pnt));
    }
  };


  /** LaplaceDLkernel is the kernel for the double layer potential of 
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM, int COMPS=1> class LaplaceDLKernel;

  /** LaplaceDLkernel in 3D reads 
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, 
          \frac{ \langle n(y), x-y\rangle }{ | x-y|^3 }, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS>
  class LaplaceDLKernel<3, COMPS> : public BaseKernel
  {
  public:
    LaplaceDLKernel<3,COMPS>()
    {
      for (size_t i = 0; i < COMPS; i++)
        terms += {1.0, 0, i, i};
    };
    typedef double value_type;
    using mp_type = typename std::conditional<COMPS == 1,
                                                  Complex,
                                                  Vec<COMPS, Complex>>::type;
    static string Name() { return "LaplaceDL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }
    
    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      // return nxy / (4 * M_PI * norm*norm*norm);
      return Vec<1,T> (nxy / (4 * M_PI * norm*norm*norm));
    }

    Array<KernelTerm> terms;

    auto CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<SingularMLExpansion<mp_type>> (c, r, 1e-16, fmm_params);
    }

    auto CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<RegularMLExpansion<mp_type>> (c, r, 1e-16, fmm_params);
    }

    void AddSource (SingularMLExpansion<mp_type> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      if constexpr (COMPS == 1)
        mp.AddDipole(pnt, -nv, val(0));
      else
        mp.AddDipole(pnt, -nv, val);
    }

    void AddSourceTrans(SingularMLExpansion<mp_type> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      if constexpr (COMPS == 1)
        mp.AddCharge (pnt, val(0));
      else
        mp.AddCharge (pnt, val);
    }

    void EvaluateMP (RegularMLExpansion<mp_type> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      if constexpr (COMPS == 1)
        val(0) = Real(mp.Evaluate (pnt));
      else
        val = Real(mp.Evaluate (pnt));
    }

    void EvaluateMPTrans(RegularMLExpansion<mp_type> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      if constexpr (COMPS == 1)
        val(0) = Real(mp.EvaluateDirectionalDerivative(pnt, nv));
      else
        val = Real(mp.EvaluateDirectionalDerivative(pnt, nv));
    }
  };


  /** HelmholtzSLkernel is the kernel for the double layer potential of the 
      Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,. $ */
  template <int DIM, int COMPS=1> class HelmholtzSLKernel;

  /** HelmholtzSLkernel in 3D reads 
      $$ G(x-y) = \frac{1 }{4\,\pi} \,\frac{e^{i\,\kappa \, |x-y| }{|x-y|} \, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<int COMPS>
  class HelmholtzSLKernel<3, COMPS> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    using mp_type = typename std::conditional<COMPS == 1,
                                              Complex,
                                              Vec<COMPS, Complex>>::type;
    static string Name() { return "HelmholtzSL"; }
    static auto Shape() { return IVec<2>(COMPS,COMPS); }
    
    /** Construction of the kernel specifies the wavenumber $\kappa$. */
    HelmholtzSLKernel (double _kappa) : kappa(_kappa) {
      for (size_t i = 0; i < COMPS; i++)
        terms += {1.0, 0, i, i};
    }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      // return kern;
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms;

    auto CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<SingularMLExpansion<mp_type>> (c, r, kappa, fmm_params);
    }

    auto CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<RegularMLExpansion<mp_type>> (c, r, kappa, fmm_params);
    }

    void AddSource (SingularMLExpansion<mp_type> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      if constexpr (COMPS == 1)
        mp.AddCharge (pnt, val(0));
      else
        mp.AddCharge (pnt, val);
    }

    void EvaluateMP (RegularMLExpansion<mp_type> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      if constexpr (COMPS == 1)
        val(0) = mp.Evaluate (pnt);
      else
        val = mp.Evaluate (pnt);
    }
  };


  /** HelmholtzDLkernel is the kernel for the double layer potential of 
      the Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,.$ */
  template <int DIM> class HelmholtzDLKernel;

  /** HelmholtzDLkernel in 3D reads
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \langle n(y), x-y\rangle \cdot \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class HelmholtzDLKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzDL"; }
    static auto Shape() { return IVec<2>(1,1); }
    
    HelmholtzDLKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>    
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
        * nxy * (Complex(1,0)*T(1.) - Complex(0,kappa)*norm);
      // return kern;
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms = { KernelTerm{1.0, 0, 0, 0}, };    

    auto CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<SingularMLExpansion<Complex>> (c, r, kappa, fmm_params);
    }

    auto CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<RegularMLExpansion<Complex>> (c, r, kappa, fmm_params);
    }

    void AddSource (SingularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      mp.AddDipole(pnt, -nv, val(0));
    }

    void AddSourceTrans(SingularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      mp.AddCharge(pnt, val(0));
    }

    void EvaluateMP (RegularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val(0) = mp.Evaluate (pnt);
    }

    void EvaluateMPTrans(RegularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val(0) = mp.EvaluateDirectionalDerivative(pnt, nv);
    }
  };


  template <int DIM> class HelmholtzHSKernel;
  
  template<>
  class HelmholtzHSKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzHS"; }
    static auto Shape() { return IVec<2>(4,4); }
    
    HelmholtzHSKernel (double _kappa) : kappa(_kappa) { }
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

    auto CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<SingularMLExpansion<Vec<6,Complex>>> (c, r, kappa, fmm_params);
    }
 
    auto CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<RegularMLExpansion<Vec<6,Complex>>> (c, r, kappa, fmm_params);
    }

    void AddSource (SingularMLExpansion<Vec<6,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      Vec<6,Complex> charge;
      charge.Range(0,3) = val.Range(0,3);
      charge.Range(3,6) = -kappa * kappa * val(3) * nv;
      mp.AddCharge(pnt, charge);
    }

    void EvaluateMP (RegularMLExpansion<Vec<6,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      Vec<6,Complex> eval = mp.Evaluate (pnt);
      val.Range(0,3) = eval.Range(0,3);
      val(3) = InnerProduct(eval.Range(3,6), nv);
    }
  };




  

  /** CombinedFieldKernel is a kernel for the combined field integral equation 
      is considered for the Helmholtz equation. */
  template <int DIM> class CombinedFieldKernel;

  /** CombinedFieldKernel in 3D reads
      $$ G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \left( \langle n_y, x-y\rangle (1- i\,\kappa\, | x-y|) - i\,\kappa\,|x-y|^2 \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class CombinedFieldKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "Helmholtz Combined Field"; }
    static auto Shape() { return IVec<2>(1,1); }
    
    CombinedFieldKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>    
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
        * ( nxy * (Complex(1,0)*T(1.) - Complex(0,kappa)*norm)  - Complex(0,kappa)*norm*norm);
      // return kern;
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms = { KernelTerm{1.0, 0, 0, 0}, };    

    auto CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<SingularMLExpansion<Complex>> (c, r, kappa, fmm_params);
    }

    auto CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<RegularMLExpansion<Complex>> (c, r, kappa, fmm_params);
    }

    void AddSource (SingularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      // mp.AddCharge(pnt, Complex(0, -kappa)*val(0));
      // mp.AddDipole(pnt, -nv, val(0));

      mp.AddChargeDipole (pnt, Complex(0, -kappa)*val(0), -nv, val(0));
    }

    void EvaluateMP (RegularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val(0) = mp.Evaluate (pnt);
    }
  };

  

  

  template <int D> class MaxwellSLKernel;

  template<>
  class MaxwellSLKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "MaxwellSL"; }
    static auto Shape() { return IVec<2>(4,4); }
    
    MaxwellSLKernel (const MaxwellSLKernel&) = default;
    MaxwellSLKernel (MaxwellSLKernel&&) = default;
    MaxwellSLKernel (double _kappa) : kappa(_kappa)
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

    auto CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<SingularMLExpansion<Vec<4,Complex>>> (c, r, kappa, fmm_params);
    }

    auto CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<RegularMLExpansion<Vec<4,Complex>>> (c, r, kappa, fmm_params);
    }

    void AddSource (SingularMLExpansion<Vec<4,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      Vec<4,Complex> charge;
      charge.Range(0,3) = kappa * val.Range(0, 3);
      charge(3) = -1.0/kappa * val(3);
      mp.AddCharge(pnt, charge);
    }

    void EvaluateMP (RegularMLExpansion<Vec<4,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val = mp.Evaluate (pnt);
    }
  };




  

  template <int D> class MaxwellDLKernel;

  // https://weggler.github.io/ngbem/short_and_sweet/Maxwell_Formulations.html
  /** MaxwellDLkernel for 3D in matrix representation reads
      $$  \left( \begin{array}{ccc} 0 & -\frac{\partial G_\kappa(x-y)}{\partial x_3} & \frac{\partial G_\kappa(x-y)}{\partial x_2} \\ \frac{\partial G_\kappa(x-y)}{\partial x_3} & 0 & -\frac{\partial G_\kappa(x-y)}{\partial x_1} \\ -\frac{\partial G_\kappa(x-y)}{\partial x_2} & \frac{\partial G_\kappa(x-y)}{\partial x_1} & 0 \end{array}\right)\,,$$ with 
   $$ G_\kappa(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \langle n(y), x-y\rangle \cdot \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class MaxwellDLKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "MaxwellDL"; }
    static auto Shape() { return IVec<2>(3,3); }
    
    MaxwellDLKernel (double _kappa) : kappa(_kappa) { }
    
    template <typename T>    
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
        * (Complex(0,kappa)*norm - Complex(1,0)*T(1.)) * (x-y);
      return kern;
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms =
      {
        KernelTerm{ 1.0, 0, 1, 2},  // factor, comp, trial, test
        KernelTerm{-1.0, 0, 2, 1},
        KernelTerm{ 1.0, 1, 2, 0},
        KernelTerm{-1.0, 1, 0, 2},
        KernelTerm{ 1.0, 2, 0, 1},
        KernelTerm{-1.0, 2, 1, 0},
      };    

    auto CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<SingularMLExpansion<Vec<3,Complex>>> (c, r, kappa, fmm_params);
    }

    auto CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<RegularMLExpansion<Vec<3,Complex>>> (c, r, kappa, fmm_params);
    }

    void AddSource (SingularMLExpansion<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
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

    void AddSourceTrans(SingularMLExpansion<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
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

    void EvaluateMP (RegularMLExpansion<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val = mp.Evaluate (pnt);
    }

    void EvaluateMPTrans(RegularMLExpansion<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val = mp.Evaluate (pnt);
    }
  };





  /*
    Dissertation Guenther Of
    "BETI–Gebietszerlegungsmethoden
    mit schnellen Randelementverfahren
    und Anwendungen"
    page 85
   */

  template <int D> class LameSLKernel;

  template<>
  class LameSLKernel<3> : public BaseKernel
  {
    double E, nu;
    double alpha;
  public:
    typedef double value_type;
    
    static string Name() { return "LameSL"; }
    static auto Shape() { return IVec<2>(3,3); }
    
    LameSLKernel (const LameSLKernel&) = default;
    LameSLKernel (LameSLKernel&&) = default;
    LameSLKernel (double _E, double _nu) : E(_E), nu(_nu)
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
      auto lapkern = alpha / (4 * M_PI * norm);  // lapkern times factor

      return Vec<7,T> { lapkern,
                        (x(0)-y(0))*(x(0)-y(0))/sqr(norm) * lapkern,
                        (x(0)-y(0))*(x(1)-y(1))/sqr(norm) * lapkern,
                        (x(0)-y(0))*(x(2)-y(2))/sqr(norm) * lapkern,
                        (x(1)-y(1))*(x(1)-y(1))/sqr(norm) * lapkern,
                        (x(1)-y(1))*(x(2)-y(2))/sqr(norm) * lapkern,
                        (x(2)-y(2))*(x(2)-y(2))/sqr(norm) * lapkern
      };
    }

    
    auto CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<SingularMLExpansion<Vec<6,Complex>>> (c, r, 1e-16, fmm_params);
    }

    auto CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
    {
      return make_shared<RegularMLExpansion<Vec<6,Complex>>> (c, r, 1e-16, fmm_params);
    }

    
    void AddSource (SingularMLExpansion<Vec<6,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      Vec<6> charge = 0.0;
      charge.Range(0,3) = val;
      mp.AddCharge(pnt, charge);   // Row 1+2

      Mat<3,3> jacobi = OuterProduct(pnt, val) + InnerProduct(pnt, val) * Id<3>();

      for (int k = 0; k < 3; k++)
        {
          Vec<6> dipole_charge = 0.0; 
          dipole_charge.Range(3,6) = jacobi.Col(k);
          
          auto ek = UnitVec<3>(k);          
          mp.AddDipole(pnt, -ek, dipole_charge);          
        }
    }

    void EvaluateMP (RegularMLExpansion<Vec<6,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      Vec<6> mpval = Real(mp.Evaluate (pnt));
      val.Range(0,3) = 0;
      val += (3-4*nu)*alpha * mpval.Range(0,3);  // Row 1
      
      val -= alpha/2 * mpval.Range(3,6);        // Row 3

      // Row 2
      Mat<3,3> jacobi = 0.0; 
      for (int k = 0; k < 3; k++)    
        {
          auto ek = UnitVec<3>(k);
          jacobi.Col(k) = Real(mp.EvaluateDirectionalDerivative(pnt, ek).Range(0,3));
        }

      val -= alpha/2 * ( Trans(jacobi) * pnt + Trace(jacobi) * pnt);
    }
  };





  

}


#endif

