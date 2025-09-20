#ifndef KERNELS_hpp
#define KERNELS_hpp

#include "mptools.hpp"


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
  template <int DIM> class LaplaceSLKernel;

  /** LaplaceSLkernel in 3D reads 
      $$ G(x-y) = \frac{1}{4\,\pi \, | x-y| }, \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class LaplaceSLKernel<3> : public BaseKernel
  {
  public:
    LaplaceSLKernel<3> () = default;
    typedef double value_type;
    static string Name() { return "LaplaceSL"; }
    static auto Shape() { return IVec<2>(1,1); }
    
    template <typename T>        
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      // return 1.0 / (4 * M_PI * norm);   
      return Vec<1,T> (1.0 / (4 * M_PI * norm));
    }
    
    Array<KernelTerm> terms = { KernelTerm{1.0, 0, 0, 0}, };
    
    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLExpansion<Complex>> (c, r, 1e-16);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLExpansion<Complex>> (c, r, 1e-16);
    }

    void AddSource (SingularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      mp.AddCharge (pnt, val(0));
    }

    void EvaluateMP (RegularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      val(0) = Real(mp.Evaluate (pnt));
    }
  };


  /** LaplaceDLkernel is the kernel for the double layer potential of 
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM> class LaplaceDLKernel;

  /** LaplaceDLkernel in 3D reads 
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, 
          \frac{ \langle n(y), x-y\rangle }{ | x-y|^3 }, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class LaplaceDLKernel<3> : public BaseKernel
  {
  public:
    typedef double value_type;
    static string Name() { return "LaplaceDL"; }
    static auto Shape() { return IVec<2>(1,1); }
    
    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      // return nxy / (4 * M_PI * norm*norm*norm);
      return Vec<1,T> (nxy / (4 * M_PI * norm*norm*norm));
    }

    Array<KernelTerm> terms = { KernelTerm{1.0, 0, 0, 0}, };    

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLExpansion<Complex>> (c, r, 1e-16);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLExpansion<Complex>> (c, r, 1e-16);
    }

    void AddSource (SingularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      mp.AddDipole(pnt, -nv, val(0));
    }

    void AddSourceTrans(SingularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      mp.AddCharge(pnt, val(0));
    }

    void EvaluateMP (RegularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      val(0) = Real(mp.Evaluate (pnt));
    }

    void EvaluateMPTrans(RegularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      val(0) = Real(mp.EvaluateDirectionalDerivative(pnt, nv));
    }
  };


  template <int DIM> class LaplaceHSKernel;
  
  template<>
  class LaplaceHSKernel<3> : public BaseKernel 
  {
  public:
    typedef double value_type;
    static string Name() { return "LaplaceHL"; }
    static auto Shape() { return IVec<2>(3,3); }
    
    template <typename T>        
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      // return 1.0 / (4 * M_PI * norm);   
      return Vec<1,T> (1.0 / (4 * M_PI * norm));
    }
    
    Array<KernelTerm> terms =
      {
        KernelTerm{1.0, 0, 0, 0},
        KernelTerm{1.0, 0, 1, 1},
        KernelTerm{1.0, 0, 2, 2},
      };

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLExpansion<Vec<3,Complex>>> (c, r, 1e-16);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLExpansion<Vec<3,Complex>>> (c, r, 1e-16);
    }

    void AddSource (SingularMLExpansion<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      mp.AddCharge(pnt, val);
    }

    void EvaluateMP (RegularMLExpansion<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      val = Real(mp.Evaluate (pnt));
    }
  };


  /** HelmholtzSLkernel is the kernel for the double layer potential of the 
      Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,. $ */
  template <int DIM> class HelmholtzSLKernel;

  /** HelmholtzSLkernel in 3D reads 
      $$ G(x-y) = \frac{1 }{4\,\pi} \,\frac{e^{i\,\kappa \, |x-y| }{|x-y|} \, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class HelmholtzSLKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzSL"; }
    static auto Shape() { return IVec<2>(1,1); }
    
    /** Construction of the kernel specifies the wavenumber $\kappa$. */
    HelmholtzSLKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      // return kern;
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms = { KernelTerm{1.0, 0, 0, 0}, };    

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLExpansion<Complex>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLExpansion<Complex>> (c, r, kappa);
    }

    void AddSource (SingularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      mp.AddCharge(pnt, val(0));
    }

    void EvaluateMP (RegularMLExpansion<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val(0) = mp.Evaluate (pnt);
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

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLExpansion<Complex>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLExpansion<Complex>> (c, r, kappa);
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


  template <int DIM> class HelmholtzSLVecKernel;

  template<>
  class HelmholtzSLVecKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzSLVec"; }
    static auto Shape() { return IVec<2>(3,3); }

    HelmholtzSLVecKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms =
      {
        KernelTerm{1.0, 0, 0, 0},
        KernelTerm{1.0, 0, 1, 1},
        KernelTerm{1.0, 0, 2, 2},
      };

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLExpansion<Vec<3,Complex>>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLExpansion<Vec<3,Complex>>> (c, r, kappa);
    }

    void AddSource (SingularMLExpansion<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      mp.AddCharge(pnt, val);
    }

    void EvaluateMP (RegularMLExpansion<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val = mp.Evaluate (pnt);
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

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLExpansion<Vec<6,Complex>>> (c, r, kappa);
    }
 
    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLExpansion<Vec<6,Complex>>> (c, r, kappa);
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

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLExpansion<Complex>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLExpansion<Complex>> (c, r, kappa);
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

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLExpansion<Vec<4,Complex>>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLExpansion<Vec<4,Complex>>> (c, r, kappa);
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

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLExpansion<Vec<3,Complex>>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLExpansion<Vec<3,Complex>>> (c, r, kappa);
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


}


#endif

