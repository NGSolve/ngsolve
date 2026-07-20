#ifndef NGBEM_INTOP_PARAMETERS_HPP
#define NGBEM_INTOP_PARAMETERS_HPP

#include <ngstd.hpp>

#include "mptools.hpp"


namespace ngsbem
{
  class IntOp_Parameters
  {
    bool use_fmm = true;
    int fmm_maxdirect = 100;
    int fmm_minorder = 20;
    double fmm_order_factor = 2.0;
    double fmm_separation = 2.0;
    double fmm_eval_separation = 3.0;
    double fmm_split_kr = 5.0;
    int fmm_maxlevel = 20;
  public:
    IntOp_Parameters () = default;
    IntOp_Parameters (const ngcore::Flags & flags);

    bool UseFMM() const { return use_fmm; }
    int FMMMaxDirect() const { return fmm_maxdirect; }
    int FMMMinOrder() const { return fmm_minorder; }
    double FMMOrderFactor() const { return fmm_order_factor; }
    double FMMSeparation() const { return fmm_separation; }
    double FMMEvalSeparation() const { return fmm_eval_separation; }
    double FMMSplitKR() const { return fmm_split_kr; }
    int FMMMaxLevel() const { return fmm_maxlevel; }

    operator FMM_Parameters() const
    {
      FMM_Parameters fmm_params;
      fmm_params.maxdirect = fmm_maxdirect;
      fmm_params.minorder = fmm_minorder;
      fmm_params.order_factor = fmm_order_factor;
      fmm_params.separation = fmm_separation;
      fmm_params.eval_separation = fmm_eval_separation;
      fmm_params.split_kr = fmm_split_kr;
      fmm_params.maxlevel = fmm_maxlevel;
      return fmm_params;
    }
  };

  inline std::ostream & operator<< (std::ostream & ost, const IntOp_Parameters & ioflags)
  {
    ost << "use_fmm = " << ioflags.UseFMM() << std::endl;
    ost << "fmm_maxdirect = " << ioflags.FMMMaxDirect() << std::endl;
    ost << "fmm_minorder = " << ioflags.FMMMinOrder() << std::endl;
    ost << "fmm_order_factor = " << ioflags.FMMOrderFactor() << std::endl;
    ost << "fmm_separation = " << ioflags.FMMSeparation() << std::endl;
    ost << "fmm_eval_separation = " << ioflags.FMMEvalSeparation() << std::endl;
    ost << "fmm_split_kr = " << ioflags.FMMSplitKR() << std::endl;
    ost << "fmm_maxlevel = " << ioflags.FMMMaxLevel() << std::endl;
    return ost;
  }
}

#endif
