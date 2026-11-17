#ifndef KERNELS_hpp
#define KERNELS_hpp

#include "mptools.hpp"
#include "integralkernel.hpp"
#include "fmmkernels.hpp"
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
