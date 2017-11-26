template <size_t H, size_t W>
void MatKernelMultAB
(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
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
template <size_t H>
void MatKernelMultABMask
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
template <size_t H, size_t W>
void MyScalTrans
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
