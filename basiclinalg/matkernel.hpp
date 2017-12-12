template <size_t H, size_t W>
static void MatKernelMultAB
(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
template <size_t H, size_t W>
static void MatKernelAlignedMultAB
(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, SIMD<double> * pc, size_t dc);
template <> void MatKernelMultAB<1, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
}
sum00.Store(pc+SW*0);
pc += dc;
}
template <> void MatKernelMultAB<2, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum10(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
}
template <> void MatKernelMultAB<3, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
sum20.Store(pc+SW*0);
pc += dc;
}
template <> void MatKernelMultAB<4, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
sum20.Store(pc+SW*0);
pc += dc;
sum30.Store(pc+SW*0);
pc += dc;
}
template <> void MatKernelMultAB<6, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b0,sum40);
SIMD<double> a5(pa[5*da]);
FMAasm(a5,b0,sum50);
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
sum20.Store(pc+SW*0);
pc += dc;
sum30.Store(pc+SW*0);
pc += dc;
sum40.Store(pc+SW*0);
pc += dc;
sum50.Store(pc+SW*0);
pc += dc;
}
template <> void MatKernelAlignedMultAB<1, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
}
pc[0]= sum00;
pc += dc;
}
template <> void MatKernelAlignedMultAB<2, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum10(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
}
pc[0]= sum00;
pc += dc;
pc[0]= sum10;
pc += dc;
}
template <> void MatKernelAlignedMultAB<3, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
}
pc[0]= sum00;
pc += dc;
pc[0]= sum10;
pc += dc;
pc[0]= sum20;
pc += dc;
}
template <> void MatKernelAlignedMultAB<4, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
}
pc[0]= sum00;
pc += dc;
pc[0]= sum10;
pc += dc;
pc[0]= sum20;
pc += dc;
pc[0]= sum30;
pc += dc;
}
template <> void MatKernelAlignedMultAB<6, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b0,sum40);
SIMD<double> a5(pa[5*da]);
FMAasm(a5,b0,sum50);
}
pc[0]= sum00;
pc += dc;
pc[0]= sum10;
pc += dc;
pc[0]= sum20;
pc += dc;
pc[0]= sum30;
pc += dc;
pc[0]= sum40;
pc += dc;
pc[0]= sum50;
pc += dc;
}
template <> void MatKernelMultAB<1, 2>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
}
template <> void MatKernelMultAB<2, 2>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
pc += dc;
}
template <> void MatKernelMultAB<3, 2>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
pc += dc;
sum20.Store(pc+SW*0);
sum21.Store(pc+SW*1);
pc += dc;
}
template <> void MatKernelMultAB<4, 2>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum30(0);
SIMD<double> sum31(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
FMAasm(a3,b1,sum31);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
pc += dc;
sum20.Store(pc+SW*0);
sum21.Store(pc+SW*1);
pc += dc;
sum30.Store(pc+SW*0);
sum31.Store(pc+SW*1);
pc += dc;
}
template <> void MatKernelMultAB<6, 2>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum30(0);
SIMD<double> sum31(0);
SIMD<double> sum40(0);
SIMD<double> sum41(0);
SIMD<double> sum50(0);
SIMD<double> sum51(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
FMAasm(a3,b1,sum31);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b0,sum40);
FMAasm(a4,b1,sum41);
SIMD<double> a5(pa[5*da]);
FMAasm(a5,b0,sum50);
FMAasm(a5,b1,sum51);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
pc += dc;
sum20.Store(pc+SW*0);
sum21.Store(pc+SW*1);
pc += dc;
sum30.Store(pc+SW*0);
sum31.Store(pc+SW*1);
pc += dc;
sum40.Store(pc+SW*0);
sum41.Store(pc+SW*1);
pc += dc;
sum50.Store(pc+SW*0);
sum51.Store(pc+SW*1);
pc += dc;
}
template <> void MatKernelAlignedMultAB<1, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> b1(pb[1]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
}
pc[0]= sum00;
pc[1]= sum01;
pc += dc;
}
template <> void MatKernelAlignedMultAB<2, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> b1(pb[1]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
}
pc[0]= sum00;
pc[1]= sum01;
pc += dc;
pc[0]= sum10;
pc[1]= sum11;
pc += dc;
}
template <> void MatKernelAlignedMultAB<3, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> b1(pb[1]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
}
pc[0]= sum00;
pc[1]= sum01;
pc += dc;
pc[0]= sum10;
pc[1]= sum11;
pc += dc;
pc[0]= sum20;
pc[1]= sum21;
pc += dc;
}
template <> void MatKernelAlignedMultAB<4, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum30(0);
SIMD<double> sum31(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> b1(pb[1]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
FMAasm(a3,b1,sum31);
}
pc[0]= sum00;
pc[1]= sum01;
pc += dc;
pc[0]= sum10;
pc[1]= sum11;
pc += dc;
pc[0]= sum20;
pc[1]= sum21;
pc += dc;
pc[0]= sum30;
pc[1]= sum31;
pc += dc;
}
template <> void MatKernelAlignedMultAB<6, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum30(0);
SIMD<double> sum31(0);
SIMD<double> sum40(0);
SIMD<double> sum41(0);
SIMD<double> sum50(0);
SIMD<double> sum51(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> b1(pb[1]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
FMAasm(a3,b1,sum31);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b0,sum40);
FMAasm(a4,b1,sum41);
SIMD<double> a5(pa[5*da]);
FMAasm(a5,b0,sum50);
FMAasm(a5,b1,sum51);
}
pc[0]= sum00;
pc[1]= sum01;
pc += dc;
pc[0]= sum10;
pc[1]= sum11;
pc += dc;
pc[0]= sum20;
pc[1]= sum21;
pc += dc;
pc[0]= sum30;
pc[1]= sum31;
pc += dc;
pc[0]= sum40;
pc[1]= sum41;
pc += dc;
pc[0]= sum50;
pc[1]= sum51;
pc += dc;
}
template <> void MatKernelMultAB<1, 3>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
pc += dc;
}
template <> void MatKernelMultAB<2, 3>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
FMAasm(a1,b2,sum12);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
sum12.Store(pc+SW*2);
pc += dc;
}
template <> void MatKernelMultAB<3, 3>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum22(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
FMAasm(a1,b2,sum12);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
FMAasm(a2,b2,sum22);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
sum12.Store(pc+SW*2);
pc += dc;
sum20.Store(pc+SW*0);
sum21.Store(pc+SW*1);
sum22.Store(pc+SW*2);
pc += dc;
}
template <> void MatKernelMultAB<4, 3>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum22(0);
SIMD<double> sum30(0);
SIMD<double> sum31(0);
SIMD<double> sum32(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
FMAasm(a1,b2,sum12);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
FMAasm(a2,b2,sum22);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
FMAasm(a3,b1,sum31);
FMAasm(a3,b2,sum32);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
sum12.Store(pc+SW*2);
pc += dc;
sum20.Store(pc+SW*0);
sum21.Store(pc+SW*1);
sum22.Store(pc+SW*2);
pc += dc;
sum30.Store(pc+SW*0);
sum31.Store(pc+SW*1);
sum32.Store(pc+SW*2);
pc += dc;
}
template <> void MatKernelMultAB<6, 3>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum22(0);
SIMD<double> sum30(0);
SIMD<double> sum31(0);
SIMD<double> sum32(0);
SIMD<double> sum40(0);
SIMD<double> sum41(0);
SIMD<double> sum42(0);
SIMD<double> sum50(0);
SIMD<double> sum51(0);
SIMD<double> sum52(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
FMAasm(a1,b2,sum12);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
FMAasm(a2,b2,sum22);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
FMAasm(a3,b1,sum31);
FMAasm(a3,b2,sum32);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b0,sum40);
FMAasm(a4,b1,sum41);
FMAasm(a4,b2,sum42);
SIMD<double> a5(pa[5*da]);
FMAasm(a5,b0,sum50);
FMAasm(a5,b1,sum51);
FMAasm(a5,b2,sum52);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
sum12.Store(pc+SW*2);
pc += dc;
sum20.Store(pc+SW*0);
sum21.Store(pc+SW*1);
sum22.Store(pc+SW*2);
pc += dc;
sum30.Store(pc+SW*0);
sum31.Store(pc+SW*1);
sum32.Store(pc+SW*2);
pc += dc;
sum40.Store(pc+SW*0);
sum41.Store(pc+SW*1);
sum42.Store(pc+SW*2);
pc += dc;
sum50.Store(pc+SW*0);
sum51.Store(pc+SW*1);
sum52.Store(pc+SW*2);
pc += dc;
}
template <> void MatKernelAlignedMultAB<1, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> b1(pb[1]);
SIMD<double> b2(pb[2]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
}
pc[0]= sum00;
pc[1]= sum01;
pc[2]= sum02;
pc += dc;
}
template <> void MatKernelAlignedMultAB<2, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> b1(pb[1]);
SIMD<double> b2(pb[2]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
FMAasm(a1,b2,sum12);
}
pc[0]= sum00;
pc[1]= sum01;
pc[2]= sum02;
pc += dc;
pc[0]= sum10;
pc[1]= sum11;
pc[2]= sum12;
pc += dc;
}
template <> void MatKernelAlignedMultAB<3, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum22(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> b1(pb[1]);
SIMD<double> b2(pb[2]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
FMAasm(a1,b2,sum12);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
FMAasm(a2,b2,sum22);
}
pc[0]= sum00;
pc[1]= sum01;
pc[2]= sum02;
pc += dc;
pc[0]= sum10;
pc[1]= sum11;
pc[2]= sum12;
pc += dc;
pc[0]= sum20;
pc[1]= sum21;
pc[2]= sum22;
pc += dc;
}
template <> void MatKernelAlignedMultAB<4, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum22(0);
SIMD<double> sum30(0);
SIMD<double> sum31(0);
SIMD<double> sum32(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> b1(pb[1]);
SIMD<double> b2(pb[2]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
FMAasm(a1,b2,sum12);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
FMAasm(a2,b2,sum22);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
FMAasm(a3,b1,sum31);
FMAasm(a3,b2,sum32);
}
pc[0]= sum00;
pc[1]= sum01;
pc[2]= sum02;
pc += dc;
pc[0]= sum10;
pc[1]= sum11;
pc[2]= sum12;
pc += dc;
pc[0]= sum20;
pc[1]= sum21;
pc[2]= sum22;
pc += dc;
pc[0]= sum30;
pc[1]= sum31;
pc[2]= sum32;
pc += dc;
}
template <> void MatKernelAlignedMultAB<6, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> * hpc = pc;
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum22(0);
SIMD<double> sum30(0);
SIMD<double> sum31(0);
SIMD<double> sum32(0);
SIMD<double> sum40(0);
SIMD<double> sum41(0);
SIMD<double> sum42(0);
SIMD<double> sum50(0);
SIMD<double> sum51(0);
SIMD<double> sum52(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> b1(pb[1]);
SIMD<double> b2(pb[2]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
FMAasm(a1,b2,sum12);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
FMAasm(a2,b2,sum22);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b0,sum30);
FMAasm(a3,b1,sum31);
FMAasm(a3,b2,sum32);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b0,sum40);
FMAasm(a4,b1,sum41);
FMAasm(a4,b2,sum42);
SIMD<double> a5(pa[5*da]);
FMAasm(a5,b0,sum50);
FMAasm(a5,b1,sum51);
FMAasm(a5,b2,sum52);
}
pc[0]= sum00;
pc[1]= sum01;
pc[2]= sum02;
pc += dc;
pc[0]= sum10;
pc[1]= sum11;
pc[2]= sum12;
pc += dc;
pc[0]= sum20;
pc[1]= sum21;
pc[2]= sum22;
pc += dc;
pc[0]= sum30;
pc[1]= sum31;
pc[2]= sum32;
pc += dc;
pc[0]= sum40;
pc[1]= sum41;
pc[2]= sum42;
pc += dc;
pc[0]= sum50;
pc[1]= sum51;
pc[2]= sum52;
pc += dc;
}
template <size_t H>
static void MatKernelMultABMask
(size_t n, SIMD<mask64> mask, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
template <> void MatKernelMultABMask<1>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum0(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b(pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
}
sum0.Store(pc,mask);
pc += dc;
}
template <> void MatKernelMultABMask<2>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum0(0);
SIMD<double> sum1(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b(pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
}
template <> void MatKernelMultABMask<3>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b(pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b,sum2);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
}
template <> void MatKernelMultABMask<4>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b(pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b,sum2);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b,sum3);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
}
template <size_t H, size_t W> static auto MatKernelScalAB
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db);
template <> INLINE auto MatKernelScalAB<3, 4>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum03(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
SIMD<double> sum13(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum22(0);
SIMD<double> sum23(0);
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> a1(pa+1*da+i);
SIMD<double> a2(pa+2*da+i);
SIMD<double> b0(pb+0*db+i);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
SIMD<double> b1(pb+1*db+i);
FMAasm(a0,b1,sum01);
FMAasm(a1,b1,sum11);
FMAasm(a2,b1,sum21);
SIMD<double> b2(pb+2*db+i);
FMAasm(a0,b2,sum02);
FMAasm(a1,b2,sum12);
FMAasm(a2,b2,sum22);
SIMD<double> b3(pb+3*db+i);
FMAasm(a0,b3,sum03);
FMAasm(a1,b3,sum13);
FMAasm(a2,b3,sum23);
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> a1(pa+1*da+i, mask);
SIMD<double> a2(pa+2*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
SIMD<double> b1(pb+1*db+i, mask);
SIMD<double> b2(pb+2*db+i, mask);
SIMD<double> b3(pb+3*db+i, mask);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
FMAasm(a0,b3,sum03);
FMAasm(a1,b0,sum10);
FMAasm(a1,b1,sum11);
FMAasm(a1,b2,sum12);
FMAasm(a1,b3,sum13);
FMAasm(a2,b0,sum20);
FMAasm(a2,b1,sum21);
FMAasm(a2,b2,sum22);
FMAasm(a2,b3,sum23);
}
return make_tuple(HSum(sum00,sum01,sum02,sum03),HSum(sum10,sum11,sum12,sum13),HSum(sum20,sum21,sum22,sum23));
}
template <> INLINE auto MatKernelScalAB<1, 4>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum03(0);
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> b0(pb+0*db+i);
FMAasm(a0,b0,sum00);
SIMD<double> b1(pb+1*db+i);
FMAasm(a0,b1,sum01);
SIMD<double> b2(pb+2*db+i);
FMAasm(a0,b2,sum02);
SIMD<double> b3(pb+3*db+i);
FMAasm(a0,b3,sum03);
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
SIMD<double> b1(pb+1*db+i, mask);
SIMD<double> b2(pb+2*db+i, mask);
SIMD<double> b3(pb+3*db+i, mask);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
FMAasm(a0,b2,sum02);
FMAasm(a0,b3,sum03);
}
return make_tuple(HSum(sum00,sum01,sum02,sum03));
}
template <> INLINE auto MatKernelScalAB<3, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> a1(pa+1*da+i);
SIMD<double> a2(pa+2*da+i);
SIMD<double> b0(pb+0*db+i);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> a1(pa+1*da+i, mask);
SIMD<double> a2(pa+2*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
}
return make_tuple(HSum(sum00),HSum(sum10),HSum(sum20));
}
template <> INLINE auto MatKernelScalAB<1, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> b0(pb+0*db+i);
FMAasm(a0,b0,sum00);
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
FMAasm(a0,b0,sum00);
}
return make_tuple(HSum(sum00));
}
template <size_t H, size_t W>
static void MyScalTrans
(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
template <> void MyScalTrans<1, 4>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
SIMD<double> sum03(pc+SW*3);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa += da, pb += db) {
SIMD<double> a0(pa[0]);
SIMD<double> b0(pb+0*SW);
FMAasm(b0,a0,sum00);
SIMD<double> b1(pb+1*SW);
FMAasm(b1,a0,sum01);
SIMD<double> b2(pb+2*SW);
FMAasm(b2,a0,sum02);
SIMD<double> b3(pb+3*SW);
FMAasm(b3,a0,sum03);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
sum03.Store(pc+SW*3);
pc += dc;
}
template <> void MyScalTrans<2, 4>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
SIMD<double> sum03(pc+SW*3);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
SIMD<double> sum13(pc+SW*3);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa += da, pb += db) {
SIMD<double> a0(pa[0]);
SIMD<double> a1(pa[1]);
SIMD<double> b0(pb+0*SW);
FMAasm(b0,a0,sum00);
FMAasm(b0,a1,sum10);
SIMD<double> b1(pb+1*SW);
FMAasm(b1,a0,sum01);
FMAasm(b1,a1,sum11);
SIMD<double> b2(pb+2*SW);
FMAasm(b2,a0,sum02);
FMAasm(b2,a1,sum12);
SIMD<double> b3(pb+3*SW);
FMAasm(b3,a0,sum03);
FMAasm(b3,a1,sum13);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
sum03.Store(pc+SW*3);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
sum12.Store(pc+SW*2);
sum13.Store(pc+SW*3);
pc += dc;
}
template <> void MyScalTrans<3, 4>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
SIMD<double> sum03(pc+SW*3);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
SIMD<double> sum13(pc+SW*3);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
SIMD<double> sum23(pc+SW*3);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa += da, pb += db) {
SIMD<double> a0(pa[0]);
SIMD<double> a1(pa[1]);
SIMD<double> a2(pa[2]);
SIMD<double> b0(pb+0*SW);
FMAasm(b0,a0,sum00);
FMAasm(b0,a1,sum10);
FMAasm(b0,a2,sum20);
SIMD<double> b1(pb+1*SW);
FMAasm(b1,a0,sum01);
FMAasm(b1,a1,sum11);
FMAasm(b1,a2,sum21);
SIMD<double> b2(pb+2*SW);
FMAasm(b2,a0,sum02);
FMAasm(b2,a1,sum12);
FMAasm(b2,a2,sum22);
SIMD<double> b3(pb+3*SW);
FMAasm(b3,a0,sum03);
FMAasm(b3,a1,sum13);
FMAasm(b3,a2,sum23);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
sum03.Store(pc+SW*3);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
sum12.Store(pc+SW*2);
sum13.Store(pc+SW*3);
pc += dc;
sum20.Store(pc+SW*0);
sum21.Store(pc+SW*1);
sum22.Store(pc+SW*2);
sum23.Store(pc+SW*3);
pc += dc;
}
template <> void MyScalTrans<4, 4>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
SIMD<double> sum03(pc+SW*3);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
SIMD<double> sum13(pc+SW*3);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
SIMD<double> sum23(pc+SW*3);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
SIMD<double> sum33(pc+SW*3);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa += da, pb += db) {
SIMD<double> a0(pa[0]);
SIMD<double> a1(pa[1]);
SIMD<double> a2(pa[2]);
SIMD<double> a3(pa[3]);
SIMD<double> b0(pb+0*SW);
FMAasm(b0,a0,sum00);
FMAasm(b0,a1,sum10);
FMAasm(b0,a2,sum20);
FMAasm(b0,a3,sum30);
SIMD<double> b1(pb+1*SW);
FMAasm(b1,a0,sum01);
FMAasm(b1,a1,sum11);
FMAasm(b1,a2,sum21);
FMAasm(b1,a3,sum31);
SIMD<double> b2(pb+2*SW);
FMAasm(b2,a0,sum02);
FMAasm(b2,a1,sum12);
FMAasm(b2,a2,sum22);
FMAasm(b2,a3,sum32);
SIMD<double> b3(pb+3*SW);
FMAasm(b3,a0,sum03);
FMAasm(b3,a1,sum13);
FMAasm(b3,a2,sum23);
FMAasm(b3,a3,sum33);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
sum03.Store(pc+SW*3);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
sum12.Store(pc+SW*2);
sum13.Store(pc+SW*3);
pc += dc;
sum20.Store(pc+SW*0);
sum21.Store(pc+SW*1);
sum22.Store(pc+SW*2);
sum23.Store(pc+SW*3);
pc += dc;
sum30.Store(pc+SW*0);
sum31.Store(pc+SW*1);
sum32.Store(pc+SW*2);
sum33.Store(pc+SW*3);
pc += dc;
}
template <> void MyScalTrans<5, 4>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
SIMD<double> sum03(pc+SW*3);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
SIMD<double> sum13(pc+SW*3);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
SIMD<double> sum23(pc+SW*3);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
SIMD<double> sum33(pc+SW*3);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
SIMD<double> sum42(pc+SW*2);
SIMD<double> sum43(pc+SW*3);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa += da, pb += db) {
SIMD<double> a0(pa[0]);
SIMD<double> a1(pa[1]);
SIMD<double> a2(pa[2]);
SIMD<double> a3(pa[3]);
SIMD<double> a4(pa[4]);
SIMD<double> b0(pb+0*SW);
FMAasm(b0,a0,sum00);
FMAasm(b0,a1,sum10);
FMAasm(b0,a2,sum20);
FMAasm(b0,a3,sum30);
FMAasm(b0,a4,sum40);
SIMD<double> b1(pb+1*SW);
FMAasm(b1,a0,sum01);
FMAasm(b1,a1,sum11);
FMAasm(b1,a2,sum21);
FMAasm(b1,a3,sum31);
FMAasm(b1,a4,sum41);
SIMD<double> b2(pb+2*SW);
FMAasm(b2,a0,sum02);
FMAasm(b2,a1,sum12);
FMAasm(b2,a2,sum22);
FMAasm(b2,a3,sum32);
FMAasm(b2,a4,sum42);
SIMD<double> b3(pb+3*SW);
FMAasm(b3,a0,sum03);
FMAasm(b3,a1,sum13);
FMAasm(b3,a2,sum23);
FMAasm(b3,a3,sum33);
FMAasm(b3,a4,sum43);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
sum03.Store(pc+SW*3);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
sum12.Store(pc+SW*2);
sum13.Store(pc+SW*3);
pc += dc;
sum20.Store(pc+SW*0);
sum21.Store(pc+SW*1);
sum22.Store(pc+SW*2);
sum23.Store(pc+SW*3);
pc += dc;
sum30.Store(pc+SW*0);
sum31.Store(pc+SW*1);
sum32.Store(pc+SW*2);
sum33.Store(pc+SW*3);
pc += dc;
sum40.Store(pc+SW*0);
sum41.Store(pc+SW*1);
sum42.Store(pc+SW*2);
sum43.Store(pc+SW*3);
pc += dc;
}
template <> void MyScalTrans<6, 4>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
SIMD<double> sum03(pc+SW*3);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
SIMD<double> sum13(pc+SW*3);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
SIMD<double> sum23(pc+SW*3);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
SIMD<double> sum33(pc+SW*3);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
SIMD<double> sum42(pc+SW*2);
SIMD<double> sum43(pc+SW*3);
pc += dc;
SIMD<double> sum50(pc+SW*0);
SIMD<double> sum51(pc+SW*1);
SIMD<double> sum52(pc+SW*2);
SIMD<double> sum53(pc+SW*3);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa += da, pb += db) {
SIMD<double> a0(pa[0]);
SIMD<double> a1(pa[1]);
SIMD<double> a2(pa[2]);
SIMD<double> a3(pa[3]);
SIMD<double> a4(pa[4]);
SIMD<double> a5(pa[5]);
SIMD<double> b0(pb+0*SW);
FMAasm(b0,a0,sum00);
FMAasm(b0,a1,sum10);
FMAasm(b0,a2,sum20);
FMAasm(b0,a3,sum30);
FMAasm(b0,a4,sum40);
FMAasm(b0,a5,sum50);
SIMD<double> b1(pb+1*SW);
FMAasm(b1,a0,sum01);
FMAasm(b1,a1,sum11);
FMAasm(b1,a2,sum21);
FMAasm(b1,a3,sum31);
FMAasm(b1,a4,sum41);
FMAasm(b1,a5,sum51);
SIMD<double> b2(pb+2*SW);
FMAasm(b2,a0,sum02);
FMAasm(b2,a1,sum12);
FMAasm(b2,a2,sum22);
FMAasm(b2,a3,sum32);
FMAasm(b2,a4,sum42);
FMAasm(b2,a5,sum52);
SIMD<double> b3(pb+3*SW);
FMAasm(b3,a0,sum03);
FMAasm(b3,a1,sum13);
FMAasm(b3,a2,sum23);
FMAasm(b3,a3,sum33);
FMAasm(b3,a4,sum43);
FMAasm(b3,a5,sum53);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
sum03.Store(pc+SW*3);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
sum12.Store(pc+SW*2);
sum13.Store(pc+SW*3);
pc += dc;
sum20.Store(pc+SW*0);
sum21.Store(pc+SW*1);
sum22.Store(pc+SW*2);
sum23.Store(pc+SW*3);
pc += dc;
sum30.Store(pc+SW*0);
sum31.Store(pc+SW*1);
sum32.Store(pc+SW*2);
sum33.Store(pc+SW*3);
pc += dc;
sum40.Store(pc+SW*0);
sum41.Store(pc+SW*1);
sum42.Store(pc+SW*2);
sum43.Store(pc+SW*3);
pc += dc;
sum50.Store(pc+SW*0);
sum51.Store(pc+SW*1);
sum52.Store(pc+SW*2);
sum53.Store(pc+SW*3);
pc += dc;
}
