#include <core/register_archive.hpp>
#include "comp.hpp"

using namespace ngfem;

static ngcore::RegisterClassForArchive<DifferentialOperator> reg_diffop;

template<typename DiffOp>
inline void RegDiffop() {
  static ngcore::RegisterClassForArchive<ngfem::T_DifferentialOperator<DiffOp>, DifferentialOperator> reg;
};

template<typename DIFFOP>
using TReg = ngcore::RegisterClassForArchive<ngfem::T_DifferentialOperator<DIFFOP>, DifferentialOperator>;

TReg<DiffOpId<1>> reg_0;
TReg<DiffOpGradient<1>> reg_1;
TReg<DiffOpIdBoundary<1>> reg_2;

TReg<DiffOpIdH1<2,2>> reg_3;
TReg<DiffOpGradient<2>> reg_4;
TReg<DiffOpIdH1<2,1>> reg_5;
TReg<DiffOpGradientBoundary<2>> reg_6;
TReg<DiffOpIdH1<2,0>> reg_7;

TReg<DiffOpIdH1<3,3>> reg_8;
TReg<DiffOpGradient<3>> reg_9;
TReg<DiffOpIdH1<3,2>> reg_10;
TReg<DiffOpGradientBoundary<3>> reg_11;
TReg<DiffOpIdH1<3,1>> reg_12;
TReg<DiffOpGradientBBoundary<3>> reg_13;
TReg<DiffOpIdH1<3,0>> reg_14;
