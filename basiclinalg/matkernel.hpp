enum OPERATION { ADD, SUB, SET, SETNEG };
template <size_t H, size_t W, OPERATION OP>
inline void MatKernelMultAB
(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
template <size_t H, size_t W, OPERATION OP>
inline void MatKernelMultAB
(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, double * pc, size_t dc);
template <size_t H, size_t W>
inline void MatKernelAlignedMultAB
(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, SIMD<double> * pc, size_t dc);
template <> INLINE void MatKernelMultAB<1, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
}
sum00.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
}
sum00.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 1, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
}
sum00.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
}
sum00.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<2, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<2, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<2, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<2, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<2, 1, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
template <> INLINE void MatKernelMultAB<2, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<2, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
template <> INLINE void MatKernelMultAB<2, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<3, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<3, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
sum20.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<3, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<3, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
sum20.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<3, 1, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
template <> INLINE void MatKernelMultAB<3, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
sum20.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<3, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
template <> INLINE void MatKernelMultAB<3, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
}
sum00.Store(pc+SW*0);
pc += dc;
sum10.Store(pc+SW*0);
pc += dc;
sum20.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<4, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<4, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
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
template <> INLINE void MatKernelMultAB<4, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<4, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
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
template <> INLINE void MatKernelMultAB<4, 1, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
template <> INLINE void MatKernelMultAB<4, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
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
template <> INLINE void MatKernelMultAB<4, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
template <> INLINE void MatKernelMultAB<4, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
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
template <> INLINE void MatKernelMultAB<5, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
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
}
template <> INLINE void MatKernelMultAB<5, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
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
}
template <> INLINE void MatKernelMultAB<5, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
pc = hpc;
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
}
template <> INLINE void MatKernelMultAB<5, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
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
}
template <> INLINE void MatKernelMultAB<5, 1, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
}
template <> INLINE void MatKernelMultAB<5, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
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
}
template <> INLINE void MatKernelMultAB<5, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
}
template <> INLINE void MatKernelMultAB<5, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
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
}
template <> INLINE void MatKernelMultAB<6, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<6, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
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
template <> INLINE void MatKernelMultAB<6, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<6, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
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
template <> INLINE void MatKernelMultAB<6, 1, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
template <> INLINE void MatKernelMultAB<6, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
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
template <> INLINE void MatKernelMultAB<6, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
template <> INLINE void MatKernelMultAB<6, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
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
template <> inline void MatKernelAlignedMultAB<1, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> sum00(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb[0]);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
}
pc[0]= sum00;
pc += dc;
}
template <> inline void MatKernelAlignedMultAB<2, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<3, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<4, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<5, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
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
}
template <> inline void MatKernelAlignedMultAB<6, 1>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> INLINE void MatKernelMultAB<1, 2, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<1, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<1, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 2, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b0,sum00);
FMAasm(a0,b1,sum01);
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
}
template <> INLINE void MatKernelMultAB<2, 2, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<2, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
pc += dc;
}
template <> INLINE void MatKernelMultAB<2, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<2, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
pc += dc;
}
template <> INLINE void MatKernelMultAB<2, 2, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
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
template <> INLINE void MatKernelMultAB<2, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
pc += dc;
}
template <> INLINE void MatKernelMultAB<2, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
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
template <> INLINE void MatKernelMultAB<2, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
pc += dc;
sum10.Store(pc+SW*0);
sum11.Store(pc+SW*1);
pc += dc;
}
template <> INLINE void MatKernelMultAB<3, 2, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<3, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
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
template <> INLINE void MatKernelMultAB<3, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<3, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
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
template <> INLINE void MatKernelMultAB<3, 2, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
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
template <> INLINE void MatKernelMultAB<3, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
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
template <> INLINE void MatKernelMultAB<3, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
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
template <> INLINE void MatKernelMultAB<3, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
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
template <> INLINE void MatKernelMultAB<4, 2, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<4, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
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
template <> INLINE void MatKernelMultAB<4, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<4, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
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
template <> INLINE void MatKernelMultAB<4, 2, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum30(0);
SIMD<double> sum31(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
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
template <> INLINE void MatKernelMultAB<4, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
SIMD<double> sum30(0);
SIMD<double> sum31(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
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
template <> INLINE void MatKernelMultAB<4, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
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
template <> INLINE void MatKernelMultAB<4, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
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
template <> INLINE void MatKernelMultAB<5, 2, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
}
template <> INLINE void MatKernelMultAB<5, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
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
}
template <> INLINE void MatKernelMultAB<5, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
pc += dc;
pc = hpc;
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
}
template <> INLINE void MatKernelMultAB<5, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
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
}
template <> INLINE void MatKernelMultAB<5, 2, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
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
}
template <> INLINE void MatKernelMultAB<5, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
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
}
template <> INLINE void MatKernelMultAB<5, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
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
}
template <> INLINE void MatKernelMultAB<5, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
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
}
template <> INLINE void MatKernelMultAB<6, 2, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<6, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
sum51 -= a5 * b1;
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
template <> INLINE void MatKernelMultAB<6, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
pc += dc;
SIMD<double> sum50(pc+SW*0);
SIMD<double> sum51(pc+SW*1);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<6, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
pc += dc;
SIMD<double> sum50(pc+SW*0);
SIMD<double> sum51(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
sum51 -= a5 * b1;
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
template <> INLINE void MatKernelMultAB<6, 2, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
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
template <> INLINE void MatKernelMultAB<6, 2, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
sum51 -= a5 * b1;
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
template <> INLINE void MatKernelMultAB<6, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
pc += dc;
SIMD<double> sum50(pc+SW*0);
SIMD<double> sum51(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
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
template <> INLINE void MatKernelMultAB<6, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
pc += dc;
SIMD<double> sum50(pc+SW*0);
SIMD<double> sum51(pc+SW*1);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
sum51 -= a5 * b1;
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
template <> inline void MatKernelAlignedMultAB<1, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<2, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<3, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<4, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<5, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
}
template <> inline void MatKernelAlignedMultAB<6, 2>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> INLINE void MatKernelMultAB<1, 3, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<1, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 3, ADD>
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
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<1, 3, SUB>
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
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 3, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
template <> INLINE void MatKernelMultAB<1, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
pc += dc;
}
template <> INLINE void MatKernelMultAB<1, 3, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
template <> INLINE void MatKernelMultAB<1, 3, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
}
sum00.Store(pc+SW*0);
sum01.Store(pc+SW*1);
sum02.Store(pc+SW*2);
pc += dc;
}
template <> INLINE void MatKernelMultAB<2, 3, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<2, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
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
template <> INLINE void MatKernelMultAB<2, 3, ADD>
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
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<2, 3, SUB>
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
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
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
template <> INLINE void MatKernelMultAB<2, 3, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
template <> INLINE void MatKernelMultAB<2, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum12(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
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
template <> INLINE void MatKernelMultAB<2, 3, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
template <> INLINE void MatKernelMultAB<2, 3, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
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
template <> INLINE void MatKernelMultAB<3, 3, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<3, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
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
template <> INLINE void MatKernelMultAB<3, 3, ADD>
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
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<3, 3, SUB>
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
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
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
template <> INLINE void MatKernelMultAB<3, 3, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
template <> INLINE void MatKernelMultAB<3, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
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
template <> INLINE void MatKernelMultAB<3, 3, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
template <> INLINE void MatKernelMultAB<3, 3, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
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
template <> INLINE void MatKernelMultAB<4, 3, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<4, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
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
template <> INLINE void MatKernelMultAB<4, 3, ADD>
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
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<4, 3, SUB>
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
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
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
template <> INLINE void MatKernelMultAB<4, 3, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
template <> INLINE void MatKernelMultAB<4, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
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
template <> INLINE void MatKernelMultAB<4, 3, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
template <> INLINE void MatKernelMultAB<4, 3, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
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
template <> INLINE void MatKernelMultAB<5, 3, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
}
template <> INLINE void MatKernelMultAB<5, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
sum42 -= a4 * b2;
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
}
template <> INLINE void MatKernelMultAB<5, 3, ADD>
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
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
SIMD<double> sum42(pc+SW*2);
pc += dc;
pc = hpc;
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
}
template <> INLINE void MatKernelMultAB<5, 3, SUB>
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
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
SIMD<double> sum42(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
sum42 -= a4 * b2;
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
}
template <> INLINE void MatKernelMultAB<5, 3, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
}
template <> INLINE void MatKernelMultAB<5, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
sum42 -= a4 * b2;
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
}
template <> INLINE void MatKernelMultAB<5, 3, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
SIMD<double> sum42(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
}
template <> INLINE void MatKernelMultAB<5, 3, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
SIMD<double> sum42(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
sum42 -= a4 * b2;
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
}
template <> INLINE void MatKernelMultAB<6, 3, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
template <> INLINE void MatKernelMultAB<6, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
sum42 -= a4 * b2;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
sum51 -= a5 * b1;
sum52 -= a5 * b2;
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
template <> INLINE void MatKernelMultAB<6, 3, ADD>
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
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
SIMD<double> sum42(pc+SW*2);
pc += dc;
SIMD<double> sum50(pc+SW*0);
SIMD<double> sum51(pc+SW*1);
SIMD<double> sum52(pc+SW*2);
pc += dc;
pc = hpc;
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
template <> INLINE void MatKernelMultAB<6, 3, SUB>
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
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
SIMD<double> sum42(pc+SW*2);
pc += dc;
SIMD<double> sum50(pc+SW*0);
SIMD<double> sum51(pc+SW*1);
SIMD<double> sum52(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> b1(pb+1*SW);
SIMD<double> b2(pb+2*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
sum42 -= a4 * b2;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
sum51 -= a5 * b1;
sum52 -= a5 * b2;
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
template <> INLINE void MatKernelMultAB<6, 3, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
template <> INLINE void MatKernelMultAB<6, 3, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
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
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
sum42 -= a4 * b2;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
sum51 -= a5 * b1;
sum52 -= a5 * b2;
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
template <> INLINE void MatKernelMultAB<6, 3, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
SIMD<double> sum42(pc+SW*2);
pc += dc;
SIMD<double> sum50(pc+SW*0);
SIMD<double> sum51(pc+SW*1);
SIMD<double> sum52(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
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
template <> INLINE void MatKernelMultAB<6, 3, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
SIMD<double> sum01(pc+SW*1);
SIMD<double> sum02(pc+SW*2);
pc += dc;
SIMD<double> sum10(pc+SW*0);
SIMD<double> sum11(pc+SW*1);
SIMD<double> sum12(pc+SW*2);
pc += dc;
SIMD<double> sum20(pc+SW*0);
SIMD<double> sum21(pc+SW*1);
SIMD<double> sum22(pc+SW*2);
pc += dc;
SIMD<double> sum30(pc+SW*0);
SIMD<double> sum31(pc+SW*1);
SIMD<double> sum32(pc+SW*2);
pc += dc;
SIMD<double> sum40(pc+SW*0);
SIMD<double> sum41(pc+SW*1);
SIMD<double> sum42(pc+SW*2);
pc += dc;
SIMD<double> sum50(pc+SW*0);
SIMD<double> sum51(pc+SW*1);
SIMD<double> sum52(pc+SW*2);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> b1 = pb[1];
SIMD<double> b2 = pb[2];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
sum01 -= a0 * b1;
sum02 -= a0 * b2;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
sum11 -= a1 * b1;
sum12 -= a1 * b2;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
sum21 -= a2 * b1;
sum22 -= a2 * b2;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
sum31 -= a3 * b1;
sum32 -= a3 * b2;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
sum41 -= a4 * b1;
sum42 -= a4 * b2;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
sum51 -= a5 * b1;
sum52 -= a5 * b2;
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
template <> inline void MatKernelAlignedMultAB<1, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<2, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<3, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<4, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> inline void MatKernelAlignedMultAB<5, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
}
template <> inline void MatKernelAlignedMultAB<6, 3>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     SIMD<double> * pc, size_t dc)
{
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
template <> INLINE void MatKernelMultAB<8, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
SIMD<double> sum60(0);
SIMD<double> sum70(0);
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
SIMD<double> a6(pa[6*da]);
FMAasm(a6,b0,sum60);
SIMD<double> a7(pa[7*da]);
FMAasm(a7,b0,sum70);
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<8, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
SIMD<double> sum60(0);
SIMD<double> sum70(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
SIMD<double> a6(pa[6*da]);
sum60 -= a6 * b0;
SIMD<double> a7(pa[7*da]);
sum70 -= a7 * b0;
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<8, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
SIMD<double> sum60(pc+SW*0);
pc += dc;
SIMD<double> sum70(pc+SW*0);
pc += dc;
pc = hpc;
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
SIMD<double> a6(pa[6*da]);
FMAasm(a6,b0,sum60);
SIMD<double> a7(pa[7*da]);
FMAasm(a7,b0,sum70);
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<8, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
SIMD<double> sum60(pc+SW*0);
pc += dc;
SIMD<double> sum70(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
SIMD<double> a6(pa[6*da]);
sum60 -= a6 * b0;
SIMD<double> a7(pa[7*da]);
sum70 -= a7 * b0;
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<8, 1, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
SIMD<double> sum60(0);
SIMD<double> sum70(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
SIMD<double> a6(pa[6*da]);
FMAasm(a6,b0,sum60);
SIMD<double> a7(pa[7*da]);
FMAasm(a7,b0,sum70);
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<8, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
SIMD<double> sum60(0);
SIMD<double> sum70(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
SIMD<double> a6(pa[6*da]);
sum60 -= a6 * b0;
SIMD<double> a7(pa[7*da]);
sum70 -= a7 * b0;
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<8, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
SIMD<double> sum60(pc+SW*0);
pc += dc;
SIMD<double> sum70(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
SIMD<double> a6(pa[6*da]);
FMAasm(a6,b0,sum60);
SIMD<double> a7(pa[7*da]);
FMAasm(a7,b0,sum70);
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<8, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
SIMD<double> sum60(pc+SW*0);
pc += dc;
SIMD<double> sum70(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
SIMD<double> a6(pa[6*da]);
sum60 -= a6 * b0;
SIMD<double> a7(pa[7*da]);
sum70 -= a7 * b0;
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<12, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
SIMD<double> sum60(0);
SIMD<double> sum70(0);
SIMD<double> sum80(0);
SIMD<double> sum90(0);
SIMD<double> sum100(0);
SIMD<double> sum110(0);
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
SIMD<double> a6(pa[6*da]);
FMAasm(a6,b0,sum60);
SIMD<double> a7(pa[7*da]);
FMAasm(a7,b0,sum70);
SIMD<double> a8(pa[8*da]);
FMAasm(a8,b0,sum80);
SIMD<double> a9(pa[9*da]);
FMAasm(a9,b0,sum90);
SIMD<double> a10(pa[10*da]);
FMAasm(a10,b0,sum100);
SIMD<double> a11(pa[11*da]);
FMAasm(a11,b0,sum110);
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
sum80.Store(pc+SW*0);
pc += dc;
sum90.Store(pc+SW*0);
pc += dc;
sum100.Store(pc+SW*0);
pc += dc;
sum110.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<12, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
SIMD<double> sum60(0);
SIMD<double> sum70(0);
SIMD<double> sum80(0);
SIMD<double> sum90(0);
SIMD<double> sum100(0);
SIMD<double> sum110(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
SIMD<double> a6(pa[6*da]);
sum60 -= a6 * b0;
SIMD<double> a7(pa[7*da]);
sum70 -= a7 * b0;
SIMD<double> a8(pa[8*da]);
sum80 -= a8 * b0;
SIMD<double> a9(pa[9*da]);
sum90 -= a9 * b0;
SIMD<double> a10(pa[10*da]);
sum100 -= a10 * b0;
SIMD<double> a11(pa[11*da]);
sum110 -= a11 * b0;
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
sum80.Store(pc+SW*0);
pc += dc;
sum90.Store(pc+SW*0);
pc += dc;
sum100.Store(pc+SW*0);
pc += dc;
sum110.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<12, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
SIMD<double> sum60(pc+SW*0);
pc += dc;
SIMD<double> sum70(pc+SW*0);
pc += dc;
SIMD<double> sum80(pc+SW*0);
pc += dc;
SIMD<double> sum90(pc+SW*0);
pc += dc;
SIMD<double> sum100(pc+SW*0);
pc += dc;
SIMD<double> sum110(pc+SW*0);
pc += dc;
pc = hpc;
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
SIMD<double> a6(pa[6*da]);
FMAasm(a6,b0,sum60);
SIMD<double> a7(pa[7*da]);
FMAasm(a7,b0,sum70);
SIMD<double> a8(pa[8*da]);
FMAasm(a8,b0,sum80);
SIMD<double> a9(pa[9*da]);
FMAasm(a9,b0,sum90);
SIMD<double> a10(pa[10*da]);
FMAasm(a10,b0,sum100);
SIMD<double> a11(pa[11*da]);
FMAasm(a11,b0,sum110);
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
sum80.Store(pc+SW*0);
pc += dc;
sum90.Store(pc+SW*0);
pc += dc;
sum100.Store(pc+SW*0);
pc += dc;
sum110.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<12, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
SIMD<double> sum60(pc+SW*0);
pc += dc;
SIMD<double> sum70(pc+SW*0);
pc += dc;
SIMD<double> sum80(pc+SW*0);
pc += dc;
SIMD<double> sum90(pc+SW*0);
pc += dc;
SIMD<double> sum100(pc+SW*0);
pc += dc;
SIMD<double> sum110(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0(pb+0*SW);
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
SIMD<double> a6(pa[6*da]);
sum60 -= a6 * b0;
SIMD<double> a7(pa[7*da]);
sum70 -= a7 * b0;
SIMD<double> a8(pa[8*da]);
sum80 -= a8 * b0;
SIMD<double> a9(pa[9*da]);
sum90 -= a9 * b0;
SIMD<double> a10(pa[10*da]);
sum100 -= a10 * b0;
SIMD<double> a11(pa[11*da]);
sum110 -= a11 * b0;
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
sum80.Store(pc+SW*0);
pc += dc;
sum90.Store(pc+SW*0);
pc += dc;
sum100.Store(pc+SW*0);
pc += dc;
sum110.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<12, 1, SET>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
SIMD<double> sum60(0);
SIMD<double> sum70(0);
SIMD<double> sum80(0);
SIMD<double> sum90(0);
SIMD<double> sum100(0);
SIMD<double> sum110(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
SIMD<double> a6(pa[6*da]);
FMAasm(a6,b0,sum60);
SIMD<double> a7(pa[7*da]);
FMAasm(a7,b0,sum70);
SIMD<double> a8(pa[8*da]);
FMAasm(a8,b0,sum80);
SIMD<double> a9(pa[9*da]);
FMAasm(a9,b0,sum90);
SIMD<double> a10(pa[10*da]);
FMAasm(a10,b0,sum100);
SIMD<double> a11(pa[11*da]);
FMAasm(a11,b0,sum110);
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
sum80.Store(pc+SW*0);
pc += dc;
sum90.Store(pc+SW*0);
pc += dc;
sum100.Store(pc+SW*0);
pc += dc;
sum110.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<12, 1, SETNEG>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
SIMD<double> sum60(0);
SIMD<double> sum70(0);
SIMD<double> sum80(0);
SIMD<double> sum90(0);
SIMD<double> sum100(0);
SIMD<double> sum110(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
SIMD<double> a6(pa[6*da]);
sum60 -= a6 * b0;
SIMD<double> a7(pa[7*da]);
sum70 -= a7 * b0;
SIMD<double> a8(pa[8*da]);
sum80 -= a8 * b0;
SIMD<double> a9(pa[9*da]);
sum90 -= a9 * b0;
SIMD<double> a10(pa[10*da]);
sum100 -= a10 * b0;
SIMD<double> a11(pa[11*da]);
sum110 -= a11 * b0;
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
sum80.Store(pc+SW*0);
pc += dc;
sum90.Store(pc+SW*0);
pc += dc;
sum100.Store(pc+SW*0);
pc += dc;
sum110.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<12, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
SIMD<double> sum60(pc+SW*0);
pc += dc;
SIMD<double> sum70(pc+SW*0);
pc += dc;
SIMD<double> sum80(pc+SW*0);
pc += dc;
SIMD<double> sum90(pc+SW*0);
pc += dc;
SIMD<double> sum100(pc+SW*0);
pc += dc;
SIMD<double> sum110(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
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
SIMD<double> a6(pa[6*da]);
FMAasm(a6,b0,sum60);
SIMD<double> a7(pa[7*da]);
FMAasm(a7,b0,sum70);
SIMD<double> a8(pa[8*da]);
FMAasm(a8,b0,sum80);
SIMD<double> a9(pa[9*da]);
FMAasm(a9,b0,sum90);
SIMD<double> a10(pa[10*da]);
FMAasm(a10,b0,sum100);
SIMD<double> a11(pa[11*da]);
FMAasm(a11,b0,sum110);
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
sum80.Store(pc+SW*0);
pc += dc;
sum90.Store(pc+SW*0);
pc += dc;
sum100.Store(pc+SW*0);
pc += dc;
sum110.Store(pc+SW*0);
pc += dc;
}
template <> INLINE void MatKernelMultAB<12, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * hpc = pc;
SIMD<double> sum00(pc+SW*0);
pc += dc;
SIMD<double> sum10(pc+SW*0);
pc += dc;
SIMD<double> sum20(pc+SW*0);
pc += dc;
SIMD<double> sum30(pc+SW*0);
pc += dc;
SIMD<double> sum40(pc+SW*0);
pc += dc;
SIMD<double> sum50(pc+SW*0);
pc += dc;
SIMD<double> sum60(pc+SW*0);
pc += dc;
SIMD<double> sum70(pc+SW*0);
pc += dc;
SIMD<double> sum80(pc+SW*0);
pc += dc;
SIMD<double> sum90(pc+SW*0);
pc += dc;
SIMD<double> sum100(pc+SW*0);
pc += dc;
SIMD<double> sum110(pc+SW*0);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b0 = pb[0];
SIMD<double> a0(pa[0*da]);
sum00 -= a0 * b0;
SIMD<double> a1(pa[1*da]);
sum10 -= a1 * b0;
SIMD<double> a2(pa[2*da]);
sum20 -= a2 * b0;
SIMD<double> a3(pa[3*da]);
sum30 -= a3 * b0;
SIMD<double> a4(pa[4*da]);
sum40 -= a4 * b0;
SIMD<double> a5(pa[5*da]);
sum50 -= a5 * b0;
SIMD<double> a6(pa[6*da]);
sum60 -= a6 * b0;
SIMD<double> a7(pa[7*da]);
sum70 -= a7 * b0;
SIMD<double> a8(pa[8*da]);
sum80 -= a8 * b0;
SIMD<double> a9(pa[9*da]);
sum90 -= a9 * b0;
SIMD<double> a10(pa[10*da]);
sum100 -= a10 * b0;
SIMD<double> a11(pa[11*da]);
sum110 -= a11 * b0;
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
sum60.Store(pc+SW*0);
pc += dc;
sum70.Store(pc+SW*0);
pc += dc;
sum80.Store(pc+SW*0);
pc += dc;
sum90.Store(pc+SW*0);
pc += dc;
sum100.Store(pc+SW*0);
pc += dc;
sum110.Store(pc+SW*0);
pc += dc;
}
template <size_t H, OPERATION OP>
inline void MatKernelMultABMask
(size_t n, SIMD<mask64> mask, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
template <size_t H, OPERATION OP>
inline void MatKernelMultABMask
(size_t n, SIMD<mask64> mask, double * pa, size_t da, SIMD<double> * pb, size_t db, double * pc, size_t dc);
template <> inline void MatKernelMultABMask<1, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
}
sum0.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<1, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
}
sum0.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<1, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
}
sum0.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<1, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
}
sum0.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<1, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
}
sum0.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<1, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
}
sum0.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<1, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
}
sum0.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<1, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
}
sum0.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<2, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<2, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<2, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<2, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<2, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<2, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<2, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<2, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<3, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<3, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<3, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<3, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<3, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<3, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<3, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<3, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<4, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<4, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
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
template <> inline void MatKernelMultABMask<4, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<4, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
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
template <> inline void MatKernelMultABMask<4, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<4, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
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
template <> inline void MatKernelMultABMask<4, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
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
template <> inline void MatKernelMultABMask<4, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
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
template <> inline void MatKernelMultABMask<5, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
SIMD<double> sum4(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b,sum2);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b,sum3);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b,sum4);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<5, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
SIMD<double> sum4(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
SIMD<double> a4(pa[4*da]);
sum4 -= a4*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<5, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
SIMD<double> sum4(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b,sum2);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b,sum3);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b,sum4);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<5, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
SIMD<double> sum4(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
SIMD<double> a4(pa[4*da]);
sum4 -= a4*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<5, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
SIMD<double> sum4(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b,sum2);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b,sum3);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b,sum4);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<5, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
SIMD<double> sum4(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
SIMD<double> a4(pa[4*da]);
sum4 -= a4*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<5, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
SIMD<double> sum4(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b,sum2);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b,sum3);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b,sum4);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<5, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
SIMD<double> sum4(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
SIMD<double> a4(pa[4*da]);
sum4 -= a4*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<6, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
SIMD<double> sum4(0);
SIMD<double> sum5(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b,sum2);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b,sum3);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b,sum4);
SIMD<double> a5(pa[5*da]);
FMAasm(a5,b,sum5);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
sum5.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<6, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
SIMD<double> sum4(0);
SIMD<double> sum5(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
SIMD<double> a4(pa[4*da]);
sum4 -= a4*b;
SIMD<double> a5(pa[5*da]);
sum5 -= a5*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
sum5.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<6, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
SIMD<double> sum4(pc, mask);
pc += dc;
SIMD<double> sum5(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b,sum2);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b,sum3);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b,sum4);
SIMD<double> a5(pa[5*da]);
FMAasm(a5,b,sum5);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
sum5.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<6, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
SIMD<double> sum4(pc, mask);
pc += dc;
SIMD<double> sum5(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
SIMD<double> a4(pa[4*da]);
sum4 -= a4*b;
SIMD<double> a5(pa[5*da]);
sum5 -= a5*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
sum5.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<6, SET>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
SIMD<double> sum4(0);
SIMD<double> sum5(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b,sum2);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b,sum3);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b,sum4);
SIMD<double> a5(pa[5*da]);
FMAasm(a5,b,sum5);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
sum5.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<6, SETNEG>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
SIMD<double> sum0(0);
SIMD<double> sum1(0);
SIMD<double> sum2(0);
SIMD<double> sum3(0);
SIMD<double> sum4(0);
SIMD<double> sum5(0);
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
SIMD<double> a4(pa[4*da]);
sum4 -= a4*b;
SIMD<double> a5(pa[5*da]);
sum5 -= a5*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
sum5.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<6, ADD>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
SIMD<double> sum4(pc, mask);
pc += dc;
SIMD<double> sum5(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
FMAasm(a0,b,sum0);
SIMD<double> a1(pa[1*da]);
FMAasm(a1,b,sum1);
SIMD<double> a2(pa[2*da]);
FMAasm(a2,b,sum2);
SIMD<double> a3(pa[3*da]);
FMAasm(a3,b,sum3);
SIMD<double> a4(pa[4*da]);
FMAasm(a4,b,sum4);
SIMD<double> a5(pa[5*da]);
FMAasm(a5,b,sum5);
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
sum5.Store(pc,mask);
pc += dc;
}
template <> inline void MatKernelMultABMask<6, SUB>
    (size_t n, SIMD<mask64> mask,
     double * pa, size_t da,
     SIMD<double> * pb, size_t db,
     double * pc, size_t dc)
{
double * hpc = pc;
SIMD<double> sum0(pc, mask);
pc += dc;
SIMD<double> sum1(pc, mask);
pc += dc;
SIMD<double> sum2(pc, mask);
pc += dc;
SIMD<double> sum3(pc, mask);
pc += dc;
SIMD<double> sum4(pc, mask);
pc += dc;
SIMD<double> sum5(pc, mask);
pc += dc;
pc = hpc;
for (size_t i = 0; i < n; i++, pa++, pb += db) {
SIMD<double> b((double*)pb,mask);
SIMD<double> a0(pa[0*da]);
sum0 -= a0*b;
SIMD<double> a1(pa[1*da]);
sum1 -= a1*b;
SIMD<double> a2(pa[2*da]);
sum2 -= a2*b;
SIMD<double> a3(pa[3*da]);
sum3 -= a3*b;
SIMD<double> a4(pa[4*da]);
sum4 -= a4*b;
SIMD<double> a5(pa[5*da]);
sum5 -= a5*b;
}
sum0.Store(pc,mask);
pc += dc;
sum1.Store(pc,mask);
pc += dc;
sum2.Store(pc,mask);
pc += dc;
sum3.Store(pc,mask);
pc += dc;
sum4.Store(pc,mask);
pc += dc;
sum5.Store(pc,mask);
pc += dc;
}
template <size_t H, size_t W> inline auto MatKernelScalAB
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db);
template <size_t H, size_t W> inline auto MatKernelScalAB
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db);
template <> INLINE auto MatKernelScalAB<6, 4>
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
SIMD<double> sum30(0);
SIMD<double> sum31(0);
SIMD<double> sum32(0);
SIMD<double> sum33(0);
SIMD<double> sum40(0);
SIMD<double> sum41(0);
SIMD<double> sum42(0);
SIMD<double> sum43(0);
SIMD<double> sum50(0);
SIMD<double> sum51(0);
SIMD<double> sum52(0);
SIMD<double> sum53(0);
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> a1(pa+1*da+i);
SIMD<double> a2(pa+2*da+i);
SIMD<double> a3(pa+3*da+i);
SIMD<double> a4(pa+4*da+i);
SIMD<double> a5(pa+5*da+i);
SIMD<double> b0(pb+0*db+i);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
FMAasm(a3,b0,sum30);
FMAasm(a4,b0,sum40);
FMAasm(a5,b0,sum50);
SIMD<double> b1(pb+1*db+i);
FMAasm(a0,b1,sum01);
FMAasm(a1,b1,sum11);
FMAasm(a2,b1,sum21);
FMAasm(a3,b1,sum31);
FMAasm(a4,b1,sum41);
FMAasm(a5,b1,sum51);
SIMD<double> b2(pb+2*db+i);
FMAasm(a0,b2,sum02);
FMAasm(a1,b2,sum12);
FMAasm(a2,b2,sum22);
FMAasm(a3,b2,sum32);
FMAasm(a4,b2,sum42);
FMAasm(a5,b2,sum52);
SIMD<double> b3(pb+3*db+i);
FMAasm(a0,b3,sum03);
FMAasm(a1,b3,sum13);
FMAasm(a2,b3,sum23);
FMAasm(a3,b3,sum33);
FMAasm(a4,b3,sum43);
FMAasm(a5,b3,sum53);
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> a1(pa+1*da+i, mask);
SIMD<double> a2(pa+2*da+i, mask);
SIMD<double> a3(pa+3*da+i, mask);
SIMD<double> a4(pa+4*da+i, mask);
SIMD<double> a5(pa+5*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
FMAasm(a3,b0,sum30);
FMAasm(a4,b0,sum40);
FMAasm(a5,b0,sum50);
SIMD<double> b1(pb+1*db+i, mask);
FMAasm(a0,b1,sum01);
FMAasm(a1,b1,sum11);
FMAasm(a2,b1,sum21);
FMAasm(a3,b1,sum31);
FMAasm(a4,b1,sum41);
FMAasm(a5,b1,sum51);
SIMD<double> b2(pb+2*db+i, mask);
FMAasm(a0,b2,sum02);
FMAasm(a1,b2,sum12);
FMAasm(a2,b2,sum22);
FMAasm(a3,b2,sum32);
FMAasm(a4,b2,sum42);
FMAasm(a5,b2,sum52);
SIMD<double> b3(pb+3*db+i, mask);
FMAasm(a0,b3,sum03);
FMAasm(a1,b3,sum13);
FMAasm(a2,b3,sum23);
FMAasm(a3,b3,sum33);
FMAasm(a4,b3,sum43);
FMAasm(a5,b3,sum53);
}
return make_tuple(HSum(sum00,sum01,sum02,sum03),HSum(sum10,sum11,sum12,sum13),HSum(sum20,sum21,sum22,sum23),HSum(sum30,sum31,sum32,sum33),HSum(sum40,sum41,sum42,sum43),HSum(sum50,sum51,sum52,sum53));
}
template <> INLINE auto MatKernelScalAB<6, 4>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
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
SIMD<double> sum30(0);
SIMD<double> sum31(0);
SIMD<double> sum32(0);
SIMD<double> sum33(0);
SIMD<double> sum40(0);
SIMD<double> sum41(0);
SIMD<double> sum42(0);
SIMD<double> sum43(0);
SIMD<double> sum50(0);
SIMD<double> sum51(0);
SIMD<double> sum52(0);
SIMD<double> sum53(0);
size_t i = 0;
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> a1(pa[1*da+i]);
SIMD<double> a2(pa[2*da+i]);
SIMD<double> a3(pa[3*da+i]);
SIMD<double> a4(pa[4*da+i]);
SIMD<double> a5(pa[5*da+i]);
SIMD<double> b0(pb[0*db+i]);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
FMAasm(a3,b0,sum30);
FMAasm(a4,b0,sum40);
FMAasm(a5,b0,sum50);
SIMD<double> b1(pb[1*db+i]);
FMAasm(a0,b1,sum01);
FMAasm(a1,b1,sum11);
FMAasm(a2,b1,sum21);
FMAasm(a3,b1,sum31);
FMAasm(a4,b1,sum41);
FMAasm(a5,b1,sum51);
SIMD<double> b2(pb[2*db+i]);
FMAasm(a0,b2,sum02);
FMAasm(a1,b2,sum12);
FMAasm(a2,b2,sum22);
FMAasm(a3,b2,sum32);
FMAasm(a4,b2,sum42);
FMAasm(a5,b2,sum52);
SIMD<double> b3(pb[3*db+i]);
FMAasm(a0,b3,sum03);
FMAasm(a1,b3,sum13);
FMAasm(a2,b3,sum23);
FMAasm(a3,b3,sum33);
FMAasm(a4,b3,sum43);
FMAasm(a5,b3,sum53);
}
return make_tuple(HSum(sum00,sum01,sum02,sum03),HSum(sum10,sum11,sum12,sum13),HSum(sum20,sum21,sum22,sum23),HSum(sum30,sum31,sum32,sum33),HSum(sum40,sum41,sum42,sum43),HSum(sum50,sum51,sum52,sum53));
}
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
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
SIMD<double> b1(pb+1*db+i, mask);
FMAasm(a0,b1,sum01);
FMAasm(a1,b1,sum11);
FMAasm(a2,b1,sum21);
SIMD<double> b2(pb+2*db+i, mask);
FMAasm(a0,b2,sum02);
FMAasm(a1,b2,sum12);
FMAasm(a2,b2,sum22);
SIMD<double> b3(pb+3*db+i, mask);
FMAasm(a0,b3,sum03);
FMAasm(a1,b3,sum13);
FMAasm(a2,b3,sum23);
}
return make_tuple(HSum(sum00,sum01,sum02,sum03),HSum(sum10,sum11,sum12,sum13),HSum(sum20,sum21,sum22,sum23));
}
template <> INLINE auto MatKernelScalAB<3, 4>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
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
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> a1(pa[1*da+i]);
SIMD<double> a2(pa[2*da+i]);
SIMD<double> b0(pb[0*db+i]);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
SIMD<double> b1(pb[1*db+i]);
FMAasm(a0,b1,sum01);
FMAasm(a1,b1,sum11);
FMAasm(a2,b1,sum21);
SIMD<double> b2(pb[2*db+i]);
FMAasm(a0,b2,sum02);
FMAasm(a1,b2,sum12);
FMAasm(a2,b2,sum22);
SIMD<double> b3(pb[3*db+i]);
FMAasm(a0,b3,sum03);
FMAasm(a1,b3,sum13);
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
sum00 += a0 * b0;
SIMD<double> b1(pb+1*db+i);
sum01 += a0 * b1;
SIMD<double> b2(pb+2*db+i);
sum02 += a0 * b2;
SIMD<double> b3(pb+3*db+i);
sum03 += a0 * b3;
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
FMAasm(a0,b0,sum00);
SIMD<double> b1(pb+1*db+i, mask);
FMAasm(a0,b1,sum01);
SIMD<double> b2(pb+2*db+i, mask);
FMAasm(a0,b2,sum02);
SIMD<double> b3(pb+3*db+i, mask);
FMAasm(a0,b3,sum03);
}
return make_tuple(HSum(sum00,sum01,sum02,sum03));
}
template <> INLINE auto MatKernelScalAB<1, 4>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum02(0);
SIMD<double> sum03(0);
size_t i = 0;
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> b0(pb[0*db+i]);
sum00 += a0 * b0;
SIMD<double> b1(pb[1*db+i]);
sum01 += a0 * b1;
SIMD<double> b2(pb[2*db+i]);
sum02 += a0 * b2;
SIMD<double> b3(pb[3*db+i]);
sum03 += a0 * b3;
}
return make_tuple(HSum(sum00,sum01,sum02,sum03));
}
template <> INLINE auto MatKernelScalAB<6, 2>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db)
{
constexpr int SW = SIMD<double>::Size();
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
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> a1(pa+1*da+i);
SIMD<double> a2(pa+2*da+i);
SIMD<double> a3(pa+3*da+i);
SIMD<double> a4(pa+4*da+i);
SIMD<double> a5(pa+5*da+i);
SIMD<double> b0(pb+0*db+i);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
FMAasm(a3,b0,sum30);
FMAasm(a4,b0,sum40);
FMAasm(a5,b0,sum50);
SIMD<double> b1(pb+1*db+i);
FMAasm(a0,b1,sum01);
FMAasm(a1,b1,sum11);
FMAasm(a2,b1,sum21);
FMAasm(a3,b1,sum31);
FMAasm(a4,b1,sum41);
FMAasm(a5,b1,sum51);
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> a1(pa+1*da+i, mask);
SIMD<double> a2(pa+2*da+i, mask);
SIMD<double> a3(pa+3*da+i, mask);
SIMD<double> a4(pa+4*da+i, mask);
SIMD<double> a5(pa+5*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
FMAasm(a3,b0,sum30);
FMAasm(a4,b0,sum40);
FMAasm(a5,b0,sum50);
SIMD<double> b1(pb+1*db+i, mask);
FMAasm(a0,b1,sum01);
FMAasm(a1,b1,sum11);
FMAasm(a2,b1,sum21);
FMAasm(a3,b1,sum31);
FMAasm(a4,b1,sum41);
FMAasm(a5,b1,sum51);
}
return make_tuple(HSum(sum00,sum01),HSum(sum10,sum11),HSum(sum20,sum21),HSum(sum30,sum31),HSum(sum40,sum41),HSum(sum50,sum51));
}
template <> INLINE auto MatKernelScalAB<6, 2>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
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
size_t i = 0;
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> a1(pa[1*da+i]);
SIMD<double> a2(pa[2*da+i]);
SIMD<double> a3(pa[3*da+i]);
SIMD<double> a4(pa[4*da+i]);
SIMD<double> a5(pa[5*da+i]);
SIMD<double> b0(pb[0*db+i]);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
FMAasm(a3,b0,sum30);
FMAasm(a4,b0,sum40);
FMAasm(a5,b0,sum50);
SIMD<double> b1(pb[1*db+i]);
FMAasm(a0,b1,sum01);
FMAasm(a1,b1,sum11);
FMAasm(a2,b1,sum21);
FMAasm(a3,b1,sum31);
FMAasm(a4,b1,sum41);
FMAasm(a5,b1,sum51);
}
return make_tuple(HSum(sum00,sum01),HSum(sum10,sum11),HSum(sum20,sum21),HSum(sum30,sum31),HSum(sum40,sum41),HSum(sum50,sum51));
}
template <> INLINE auto MatKernelScalAB<3, 2>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> a1(pa+1*da+i);
SIMD<double> a2(pa+2*da+i);
SIMD<double> b0(pb+0*db+i);
sum00 += a0 * b0;
sum10 += a1 * b0;
sum20 += a2 * b0;
SIMD<double> b1(pb+1*db+i);
sum01 += a0 * b1;
sum11 += a1 * b1;
sum21 += a2 * b1;
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
SIMD<double> b1(pb+1*db+i, mask);
FMAasm(a0,b1,sum01);
FMAasm(a1,b1,sum11);
FMAasm(a2,b1,sum21);
}
return make_tuple(HSum(sum00,sum01),HSum(sum10,sum11),HSum(sum20,sum21));
}
template <> INLINE auto MatKernelScalAB<3, 2>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
SIMD<double> sum00(0);
SIMD<double> sum01(0);
SIMD<double> sum10(0);
SIMD<double> sum11(0);
SIMD<double> sum20(0);
SIMD<double> sum21(0);
size_t i = 0;
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> a1(pa[1*da+i]);
SIMD<double> a2(pa[2*da+i]);
SIMD<double> b0(pb[0*db+i]);
sum00 += a0 * b0;
sum10 += a1 * b0;
sum20 += a2 * b0;
SIMD<double> b1(pb[1*db+i]);
sum01 += a0 * b1;
sum11 += a1 * b1;
sum21 += a2 * b1;
}
return make_tuple(HSum(sum00,sum01),HSum(sum10,sum11),HSum(sum20,sum21));
}
template <> INLINE auto MatKernelScalAB<8, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
SIMD<double> sum60(0);
SIMD<double> sum70(0);
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> a1(pa+1*da+i);
SIMD<double> a2(pa+2*da+i);
SIMD<double> a3(pa+3*da+i);
SIMD<double> a4(pa+4*da+i);
SIMD<double> a5(pa+5*da+i);
SIMD<double> a6(pa+6*da+i);
SIMD<double> a7(pa+7*da+i);
SIMD<double> b0(pb+0*db+i);
sum00 += a0 * b0;
sum10 += a1 * b0;
sum20 += a2 * b0;
sum30 += a3 * b0;
sum40 += a4 * b0;
sum50 += a5 * b0;
sum60 += a6 * b0;
sum70 += a7 * b0;
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> a1(pa+1*da+i, mask);
SIMD<double> a2(pa+2*da+i, mask);
SIMD<double> a3(pa+3*da+i, mask);
SIMD<double> a4(pa+4*da+i, mask);
SIMD<double> a5(pa+5*da+i, mask);
SIMD<double> a6(pa+6*da+i, mask);
SIMD<double> a7(pa+7*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
FMAasm(a3,b0,sum30);
FMAasm(a4,b0,sum40);
FMAasm(a5,b0,sum50);
FMAasm(a6,b0,sum60);
FMAasm(a7,b0,sum70);
}
return make_tuple(HSum(sum00, sum10, sum20, sum30),HSum(sum40, sum50, sum60, sum70));
}
template <> INLINE auto MatKernelScalAB<8, 1>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
SIMD<double> sum60(0);
SIMD<double> sum70(0);
size_t i = 0;
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> a1(pa[1*da+i]);
SIMD<double> a2(pa[2*da+i]);
SIMD<double> a3(pa[3*da+i]);
SIMD<double> a4(pa[4*da+i]);
SIMD<double> a5(pa[5*da+i]);
SIMD<double> a6(pa[6*da+i]);
SIMD<double> a7(pa[7*da+i]);
SIMD<double> b0(pb[0*db+i]);
sum00 += a0 * b0;
sum10 += a1 * b0;
sum20 += a2 * b0;
sum30 += a3 * b0;
sum40 += a4 * b0;
sum50 += a5 * b0;
sum60 += a6 * b0;
sum70 += a7 * b0;
}
return make_tuple(HSum(sum00, sum10, sum20, sum30),HSum(sum40, sum50, sum60, sum70));
}
template <> INLINE auto MatKernelScalAB<6, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> a1(pa+1*da+i);
SIMD<double> a2(pa+2*da+i);
SIMD<double> a3(pa+3*da+i);
SIMD<double> a4(pa+4*da+i);
SIMD<double> a5(pa+5*da+i);
SIMD<double> b0(pb+0*db+i);
sum00 += a0 * b0;
sum10 += a1 * b0;
sum20 += a2 * b0;
sum30 += a3 * b0;
sum40 += a4 * b0;
sum50 += a5 * b0;
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> a1(pa+1*da+i, mask);
SIMD<double> a2(pa+2*da+i, mask);
SIMD<double> a3(pa+3*da+i, mask);
SIMD<double> a4(pa+4*da+i, mask);
SIMD<double> a5(pa+5*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
FMAasm(a3,b0,sum30);
FMAasm(a4,b0,sum40);
FMAasm(a5,b0,sum50);
}
return make_tuple(HSum(sum00),HSum(sum10),HSum(sum20),HSum(sum30),HSum(sum40),HSum(sum50));
}
template <> INLINE auto MatKernelScalAB<6, 1>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
SIMD<double> sum40(0);
SIMD<double> sum50(0);
size_t i = 0;
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> a1(pa[1*da+i]);
SIMD<double> a2(pa[2*da+i]);
SIMD<double> a3(pa[3*da+i]);
SIMD<double> a4(pa[4*da+i]);
SIMD<double> a5(pa[5*da+i]);
SIMD<double> b0(pb[0*db+i]);
sum00 += a0 * b0;
sum10 += a1 * b0;
sum20 += a2 * b0;
sum30 += a3 * b0;
sum40 += a4 * b0;
sum50 += a5 * b0;
}
return make_tuple(HSum(sum00),HSum(sum10),HSum(sum20),HSum(sum30),HSum(sum40),HSum(sum50));
}
template <> INLINE auto MatKernelScalAB<4, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> a1(pa+1*da+i);
SIMD<double> a2(pa+2*da+i);
SIMD<double> a3(pa+3*da+i);
SIMD<double> b0(pb+0*db+i);
sum00 += a0 * b0;
sum10 += a1 * b0;
sum20 += a2 * b0;
sum30 += a3 * b0;
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> a1(pa+1*da+i, mask);
SIMD<double> a2(pa+2*da+i, mask);
SIMD<double> a3(pa+3*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
FMAasm(a2,b0,sum20);
FMAasm(a3,b0,sum30);
}
return make_tuple(HSum(sum00, sum10, sum20, sum30));
}
template <> INLINE auto MatKernelScalAB<4, 1>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
SIMD<double> sum30(0);
size_t i = 0;
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> a1(pa[1*da+i]);
SIMD<double> a2(pa[2*da+i]);
SIMD<double> a3(pa[3*da+i]);
SIMD<double> b0(pb[0*db+i]);
sum00 += a0 * b0;
sum10 += a1 * b0;
sum20 += a2 * b0;
sum30 += a3 * b0;
}
return make_tuple(HSum(sum00, sum10, sum20, sum30));
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
sum00 += a0 * b0;
sum10 += a1 * b0;
sum20 += a2 * b0;
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
template <> INLINE auto MatKernelScalAB<3, 1>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
SIMD<double> sum00(0);
SIMD<double> sum10(0);
SIMD<double> sum20(0);
size_t i = 0;
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> a1(pa[1*da+i]);
SIMD<double> a2(pa[2*da+i]);
SIMD<double> b0(pb[0*db+i]);
sum00 += a0 * b0;
sum10 += a1 * b0;
sum20 += a2 * b0;
}
return make_tuple(HSum(sum00),HSum(sum10),HSum(sum20));
}
template <> INLINE auto MatKernelScalAB<2, 1>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> sum00(0);
SIMD<double> sum10(0);
size_t i = 0;
for ( ; i+SW <= n; i+=SW) {
SIMD<double> a0(pa+0*da+i);
SIMD<double> a1(pa+1*da+i);
SIMD<double> b0(pb+0*db+i);
sum00 += a0 * b0;
sum10 += a1 * b0;
}
size_t r = n % SW;
if (r) {
SIMD<mask64> mask(r);
SIMD<double> a0(pa+0*da+i, mask);
SIMD<double> a1(pa+1*da+i, mask);
SIMD<double> b0(pb+0*db+i, mask);
FMAasm(a0,b0,sum00);
FMAasm(a1,b0,sum10);
}
return make_tuple(HSum(sum00),HSum(sum10));
}
template <> INLINE auto MatKernelScalAB<2, 1>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
SIMD<double> sum00(0);
SIMD<double> sum10(0);
size_t i = 0;
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> a1(pa[1*da+i]);
SIMD<double> b0(pb[0*db+i]);
sum00 += a0 * b0;
sum10 += a1 * b0;
}
return make_tuple(HSum(sum00),HSum(sum10));
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
sum00 += a0 * b0;
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
template <> INLINE auto MatKernelScalAB<1, 1>
    (size_t n,
     SIMD<double> * pa, size_t da,
     SIMD<double> * pb, size_t db)
{
SIMD<double> sum00(0);
size_t i = 0;
for ( ; i < n; i++) {
SIMD<double> a0(pa[0*da+i]);
SIMD<double> b0(pb[0*db+i]);
sum00 += a0 * b0;
}
return make_tuple(HSum(sum00));
}
template <size_t H, size_t W>
inline void MyScalTrans
(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
template <> inline void MyScalTrans<1, 4>
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
template <> inline void MyScalTrans<2, 4>
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
template <> inline void MyScalTrans<3, 4>
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
template <> inline void MyScalTrans<4, 4>
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
template <> inline void MyScalTrans<5, 4>
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
template <> inline void MyScalTrans<6, 4>
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
template <size_t H, size_t W, OPERATION OP>
inline void MatKernelDaxpy
(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
template <size_t H, size_t W, OPERATION OP>
inline void MatKernelDaxpy
(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, SIMD<double> * pc, size_t dc);
template <> INLINE void MatKernelDaxpy<1, 0, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * pc0 = pc+0*dc;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 0, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * pc0 = pc+0*dc;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 0, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
double * pc0 = pc+0*dc;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 2, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 3, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 3, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 3, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 4, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 4, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 4, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 5, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 5, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 5, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i);
c0 -= a04 * b4;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 -= a04 * b4;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 6, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 6, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 6, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i);
c0 -= a05 * b5;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 -= a05 * b5;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 7, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 7, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 7, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i);
c0 -= a06 * b6;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 -= a06 * b6;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 8, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
SIMD<double> b7(pb7+i);
c0 += a07 * b7;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 += a07 * b7;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 8, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
SIMD<double> b7(pb7+i);
c0 += a07 * b7;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 += a07 * b7;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 8, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i);
c0 -= a06 * b6;
SIMD<double> b7(pb7+i);
c0 -= a07 * b7;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 -= a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 -= a07 * b7;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 9, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
SIMD<double> b7(pb7+i);
c0 += a07 * b7;
SIMD<double> b8(pb8+i);
c0 += a08 * b8;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 += a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 += a08 * b8;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 9, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
SIMD<double> b7(pb7+i);
c0 += a07 * b7;
SIMD<double> b8(pb8+i);
c0 += a08 * b8;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 += a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 += a08 * b8;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 9, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i);
c0 -= a06 * b6;
SIMD<double> b7(pb7+i);
c0 -= a07 * b7;
SIMD<double> b8(pb8+i);
c0 -= a08 * b8;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 -= a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 -= a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 -= a08 * b8;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 10, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
SIMD<double> a09(pa[0*da+9]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
double * pb9 = pb+9*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
SIMD<double> b7(pb7+i);
c0 += a07 * b7;
SIMD<double> b8(pb8+i);
c0 += a08 * b8;
SIMD<double> b9(pb9+i);
c0 += a09 * b9;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 += a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 += a08 * b8;
SIMD<double> b9(pb9+i, mask);
c0 += a09 * b9;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 10, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
SIMD<double> a09(pa[0*da+9]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
double * pb9 = pb+9*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
SIMD<double> b7(pb7+i);
c0 += a07 * b7;
SIMD<double> b8(pb8+i);
c0 += a08 * b8;
SIMD<double> b9(pb9+i);
c0 += a09 * b9;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 += a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 += a08 * b8;
SIMD<double> b9(pb9+i, mask);
c0 += a09 * b9;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 10, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
SIMD<double> a09(pa[0*da+9]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
double * pb9 = pb+9*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i);
c0 -= a06 * b6;
SIMD<double> b7(pb7+i);
c0 -= a07 * b7;
SIMD<double> b8(pb8+i);
c0 -= a08 * b8;
SIMD<double> b9(pb9+i);
c0 -= a09 * b9;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 -= a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 -= a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 -= a08 * b8;
SIMD<double> b9(pb9+i, mask);
c0 -= a09 * b9;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 11, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
SIMD<double> a09(pa[0*da+9]);
SIMD<double> a010(pa[0*da+10]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
double * pb9 = pb+9*db;
double * pb10 = pb+10*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
SIMD<double> b7(pb7+i);
c0 += a07 * b7;
SIMD<double> b8(pb8+i);
c0 += a08 * b8;
SIMD<double> b9(pb9+i);
c0 += a09 * b9;
SIMD<double> b10(pb10+i);
c0 += a010 * b10;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 += a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 += a08 * b8;
SIMD<double> b9(pb9+i, mask);
c0 += a09 * b9;
SIMD<double> b10(pb10+i, mask);
c0 += a010 * b10;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 11, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
SIMD<double> a09(pa[0*da+9]);
SIMD<double> a010(pa[0*da+10]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
double * pb9 = pb+9*db;
double * pb10 = pb+10*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
SIMD<double> b7(pb7+i);
c0 += a07 * b7;
SIMD<double> b8(pb8+i);
c0 += a08 * b8;
SIMD<double> b9(pb9+i);
c0 += a09 * b9;
SIMD<double> b10(pb10+i);
c0 += a010 * b10;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 += a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 += a08 * b8;
SIMD<double> b9(pb9+i, mask);
c0 += a09 * b9;
SIMD<double> b10(pb10+i, mask);
c0 += a010 * b10;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 11, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
SIMD<double> a09(pa[0*da+9]);
SIMD<double> a010(pa[0*da+10]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
double * pb9 = pb+9*db;
double * pb10 = pb+10*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i);
c0 -= a06 * b6;
SIMD<double> b7(pb7+i);
c0 -= a07 * b7;
SIMD<double> b8(pb8+i);
c0 -= a08 * b8;
SIMD<double> b9(pb9+i);
c0 -= a09 * b9;
SIMD<double> b10(pb10+i);
c0 -= a010 * b10;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 -= a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 -= a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 -= a08 * b8;
SIMD<double> b9(pb9+i, mask);
c0 -= a09 * b9;
SIMD<double> b10(pb10+i, mask);
c0 -= a010 * b10;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 12, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
SIMD<double> a09(pa[0*da+9]);
SIMD<double> a010(pa[0*da+10]);
SIMD<double> a011(pa[0*da+11]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
double * pb9 = pb+9*db;
double * pb10 = pb+10*db;
double * pb11 = pb+11*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
SIMD<double> b7(pb7+i);
c0 += a07 * b7;
SIMD<double> b8(pb8+i);
c0 += a08 * b8;
SIMD<double> b9(pb9+i);
c0 += a09 * b9;
SIMD<double> b10(pb10+i);
c0 += a010 * b10;
SIMD<double> b11(pb11+i);
c0 += a011 * b11;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 += a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 += a08 * b8;
SIMD<double> b9(pb9+i, mask);
c0 += a09 * b9;
SIMD<double> b10(pb10+i, mask);
c0 += a010 * b10;
SIMD<double> b11(pb11+i, mask);
c0 += a011 * b11;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 12, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
SIMD<double> a09(pa[0*da+9]);
SIMD<double> a010(pa[0*da+10]);
SIMD<double> a011(pa[0*da+11]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
double * pb9 = pb+9*db;
double * pb10 = pb+10*db;
double * pb11 = pb+11*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
SIMD<double> b4(pb4+i);
c0 += a04 * b4;
SIMD<double> b5(pb5+i);
c0 += a05 * b5;
SIMD<double> b6(pb6+i);
c0 += a06 * b6;
SIMD<double> b7(pb7+i);
c0 += a07 * b7;
SIMD<double> b8(pb8+i);
c0 += a08 * b8;
SIMD<double> b9(pb9+i);
c0 += a09 * b9;
SIMD<double> b10(pb10+i);
c0 += a010 * b10;
SIMD<double> b11(pb11+i);
c0 += a011 * b11;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 += a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 += a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 += a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 += a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 += a08 * b8;
SIMD<double> b9(pb9+i, mask);
c0 += a09 * b9;
SIMD<double> b10(pb10+i, mask);
c0 += a010 * b10;
SIMD<double> b11(pb11+i, mask);
c0 += a011 * b11;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<1, 12, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a04(pa[0*da+4]);
SIMD<double> a05(pa[0*da+5]);
SIMD<double> a06(pa[0*da+6]);
SIMD<double> a07(pa[0*da+7]);
SIMD<double> a08(pa[0*da+8]);
SIMD<double> a09(pa[0*da+9]);
SIMD<double> a010(pa[0*da+10]);
SIMD<double> a011(pa[0*da+11]);
double * pc0 = pc+0*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
double * pb4 = pb+4*db;
double * pb5 = pb+5*db;
double * pb6 = pb+6*db;
double * pb7 = pb+7*db;
double * pb8 = pb+8*db;
double * pb9 = pb+9*db;
double * pb10 = pb+10*db;
double * pb11 = pb+11*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i);
c0 -= a06 * b6;
SIMD<double> b7(pb7+i);
c0 -= a07 * b7;
SIMD<double> b8(pb8+i);
c0 -= a08 * b8;
SIMD<double> b9(pb9+i);
c0 -= a09 * b9;
SIMD<double> b10(pb10+i);
c0 -= a010 * b10;
SIMD<double> b11(pb11+i);
c0 -= a011 * b11;
c0.Store(pc0+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
SIMD<double> b4(pb4+i, mask);
c0 -= a04 * b4;
SIMD<double> b5(pb5+i, mask);
c0 -= a05 * b5;
SIMD<double> b6(pb6+i, mask);
c0 -= a06 * b6;
SIMD<double> b7(pb7+i, mask);
c0 -= a07 * b7;
SIMD<double> b8(pb8+i, mask);
c0 -= a08 * b8;
SIMD<double> b9(pb9+i, mask);
c0 -= a09 * b9;
SIMD<double> b10(pb10+i, mask);
c0 -= a010 * b10;
SIMD<double> b11(pb11+i, mask);
c0 -= a011 * b11;
c0.Store(pc0+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a10(pa[1*da+0]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a10(pa[1*da+0]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a10(pa[1*da+0]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
c1 -= a10 * b0;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
c1 -= a10 * b0;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 2, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
c1 -= a10 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
c1 -= a11 * b1;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
c1 -= a10 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
c1 -= a11 * b1;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 3, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
c1 += a12 * b2;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
c1 += a12 * b2;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 3, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
c1 += a12 * b2;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
c1 += a12 * b2;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 3, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
c1 -= a10 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
c1 -= a11 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
c1 -= a12 * b2;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
c1 -= a10 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
c1 -= a11 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
c1 -= a12 * b2;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 4, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
SIMD<double> a13(pa[1*da+3]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
c1 += a12 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
c1 += a13 * b3;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
c1 += a12 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
c1 += a13 * b3;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 4, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
SIMD<double> a13(pa[1*da+3]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
c1 += a12 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
c1 += a13 * b3;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
c1 += a12 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
c1 += a13 * b3;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<2, 4, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
SIMD<double> a13(pa[1*da+3]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
c1 -= a10 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
c1 -= a11 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
c1 -= a12 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
c1 -= a13 * b3;
c0.Store(pc0+i);
c1.Store(pc1+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
c1 -= a10 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
c1 -= a11 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
c1 -= a12 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
c1 -= a13 * b3;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 1, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a20(pa[2*da+0]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> c2(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> c2(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 1, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a20(pa[2*da+0]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> c2(pc2+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> c2(pc2+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 1, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a20(pa[2*da+0]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> c2(pc2+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
c1 -= a10 * b0;
c2 -= a20 * b0;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> c2(pc2+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
c1 -= a10 * b0;
c2 -= a20 * b0;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 2, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a20(pa[2*da+0]);
SIMD<double> a21(pa[2*da+1]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> c2(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> c2(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 2, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a20(pa[2*da+0]);
SIMD<double> a21(pa[2*da+1]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> c2(pc2+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> c2(pc2+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 2, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a20(pa[2*da+0]);
SIMD<double> a21(pa[2*da+1]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> c2(pc2+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
c1 -= a10 * b0;
c2 -= a20 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
c1 -= a11 * b1;
c2 -= a21 * b1;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> c2(pc2+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
c1 -= a10 * b0;
c2 -= a20 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
c1 -= a11 * b1;
c2 -= a21 * b1;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 3, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
SIMD<double> a20(pa[2*da+0]);
SIMD<double> a21(pa[2*da+1]);
SIMD<double> a22(pa[2*da+2]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> c2(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
c1 += a12 * b2;
c2 += a22 * b2;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> c2(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
c1 += a12 * b2;
c2 += a22 * b2;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 3, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
SIMD<double> a20(pa[2*da+0]);
SIMD<double> a21(pa[2*da+1]);
SIMD<double> a22(pa[2*da+2]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> c2(pc2+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
c1 += a12 * b2;
c2 += a22 * b2;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> c2(pc2+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
c1 += a12 * b2;
c2 += a22 * b2;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 3, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
SIMD<double> a20(pa[2*da+0]);
SIMD<double> a21(pa[2*da+1]);
SIMD<double> a22(pa[2*da+2]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> c2(pc2+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
c1 -= a10 * b0;
c2 -= a20 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
c1 -= a11 * b1;
c2 -= a21 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
c1 -= a12 * b2;
c2 -= a22 * b2;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> c2(pc2+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
c1 -= a10 * b0;
c2 -= a20 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
c1 -= a11 * b1;
c2 -= a21 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
c1 -= a12 * b2;
c2 -= a22 * b2;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 4, SET>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
SIMD<double> a13(pa[1*da+3]);
SIMD<double> a20(pa[2*da+0]);
SIMD<double> a21(pa[2*da+1]);
SIMD<double> a22(pa[2*da+2]);
SIMD<double> a23(pa[2*da+3]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> c2(0);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
c1 += a12 * b2;
c2 += a22 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
c1 += a13 * b3;
c2 += a23 * b3;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(0);
SIMD<double> c1(0);
SIMD<double> c2(0);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
c1 += a12 * b2;
c2 += a22 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
c1 += a13 * b3;
c2 += a23 * b3;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 4, ADD>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
SIMD<double> a13(pa[1*da+3]);
SIMD<double> a20(pa[2*da+0]);
SIMD<double> a21(pa[2*da+1]);
SIMD<double> a22(pa[2*da+2]);
SIMD<double> a23(pa[2*da+3]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> c2(pc2+i);
SIMD<double> b0(pb0+i);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
SIMD<double> b2(pb2+i);
c0 += a02 * b2;
c1 += a12 * b2;
c2 += a22 * b2;
SIMD<double> b3(pb3+i);
c0 += a03 * b3;
c1 += a13 * b3;
c2 += a23 * b3;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> c2(pc2+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 += a00 * b0;
c1 += a10 * b0;
c2 += a20 * b0;
SIMD<double> b1(pb1+i, mask);
c0 += a01 * b1;
c1 += a11 * b1;
c2 += a21 * b1;
SIMD<double> b2(pb2+i, mask);
c0 += a02 * b2;
c1 += a12 * b2;
c2 += a22 * b2;
SIMD<double> b3(pb3+i, mask);
c0 += a03 * b3;
c1 += a13 * b3;
c2 += a23 * b3;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
template <> INLINE void MatKernelDaxpy<3, 4, SUB>
    (size_t n,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
SIMD<double> a00(pa[0*da+0]);
SIMD<double> a01(pa[0*da+1]);
SIMD<double> a02(pa[0*da+2]);
SIMD<double> a03(pa[0*da+3]);
SIMD<double> a10(pa[1*da+0]);
SIMD<double> a11(pa[1*da+1]);
SIMD<double> a12(pa[1*da+2]);
SIMD<double> a13(pa[1*da+3]);
SIMD<double> a20(pa[2*da+0]);
SIMD<double> a21(pa[2*da+1]);
SIMD<double> a22(pa[2*da+2]);
SIMD<double> a23(pa[2*da+3]);
double * pc0 = pc+0*dc;
double * pc1 = pc+1*dc;
double * pc2 = pc+2*dc;
double * pb0 = pb+0*db;
double * pb1 = pb+1*db;
double * pb2 = pb+2*db;
double * pb3 = pb+3*db;
size_t i = 0;
for ( ; i+SW < n; i+=SW) {
SIMD<double> c0(pc0+i);
SIMD<double> c1(pc1+i);
SIMD<double> c2(pc2+i);
SIMD<double> b0(pb0+i);
c0 -= a00 * b0;
c1 -= a10 * b0;
c2 -= a20 * b0;
SIMD<double> b1(pb1+i);
c0 -= a01 * b1;
c1 -= a11 * b1;
c2 -= a21 * b1;
SIMD<double> b2(pb2+i);
c0 -= a02 * b2;
c1 -= a12 * b2;
c2 -= a22 * b2;
SIMD<double> b3(pb3+i);
c0 -= a03 * b3;
c1 -= a13 * b3;
c2 -= a23 * b3;
c0.Store(pc0+i);
c1.Store(pc1+i);
c2.Store(pc2+i);
}
SIMD<mask64> mask(n%SW);
SIMD<double> c0(pc0+i, mask);
SIMD<double> c1(pc1+i, mask);
SIMD<double> c2(pc2+i, mask);
SIMD<double> b0(pb0+i, mask);
c0 -= a00 * b0;
c1 -= a10 * b0;
c2 -= a20 * b0;
SIMD<double> b1(pb1+i, mask);
c0 -= a01 * b1;
c1 -= a11 * b1;
c2 -= a21 * b1;
SIMD<double> b2(pb2+i, mask);
c0 -= a02 * b2;
c1 -= a12 * b2;
c2 -= a22 * b2;
SIMD<double> b3(pb3+i, mask);
c0 -= a03 * b3;
c1 -= a13 * b3;
c2 -= a23 * b3;
c0.Store(pc0+i, mask);
c1.Store(pc1+i, mask);
c2.Store(pc2+i, mask);
}
// C = A * B,  with short inner loop
template <size_t WA, OPERATION OP>
inline void MatKernelShortSum
(size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
// C = A * B,  with short inner loop, unroll width B
template <size_t WA, OPERATION OP>
inline void MatKernelShortSum2
(size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
template <> INLINE void MatKernelShortSum<0, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<0, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<0, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<0, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<1, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<1, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<1, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<1, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<2, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<2, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<2, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<2, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<3, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<3, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<3, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<3, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<4, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<4, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<4, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<4, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<5, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<5, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<5, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<5, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<6, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<6, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<6, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<6, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<7, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<7, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<7, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<7, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<8, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<8, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<8, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<8, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<9, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<9, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW,mask); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<9, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<9, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW,mask); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<10, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<10, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW,mask); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW,mask); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<10, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<10, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW,mask); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW,mask); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<11, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
SIMD<double> b10(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
SIMD<double> b10(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<11, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW); pb2 += db;
SIMD<double> b100(pb2);
SIMD<double> b101(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0 += SIMD<double>(pa2[10]) * b100;
sum1 += SIMD<double>(pa2[10]) * b101;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW,mask); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW,mask); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW,mask); pb2 += db;
SIMD<double> b100(pb2);
SIMD<double> b101(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0 += SIMD<double>(pa2[10]) * b100;
sum1 += SIMD<double>(pa2[10]) * b101;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
SIMD<double> b10(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
SIMD<double> b10(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<11, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
SIMD<double> b10(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
SIMD<double> b10(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<11, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW); pb2 += db;
SIMD<double> b100(pb2);
SIMD<double> b101(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0 += SIMD<double>(pa2[10]) * b100;
sum1 += SIMD<double>(pa2[10]) * b101;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW,mask); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW,mask); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW,mask); pb2 += db;
SIMD<double> b100(pb2);
SIMD<double> b101(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0 += SIMD<double>(pa2[10]) * b100;
sum1 += SIMD<double>(pa2[10]) * b101;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
SIMD<double> b10(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
SIMD<double> b10(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<12, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
SIMD<double> b10(pb2); pb2 += db;
SIMD<double> b11(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum += SIMD<double>(pa2[11]) * b11;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
SIMD<double> b10(pb2, mask); pb2 += db;
SIMD<double> b11(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum += SIMD<double>(pa2[11]) * b11;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<12, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW); pb2 += db;
SIMD<double> b100(pb2);
SIMD<double> b101(pb2+SW); pb2 += db;
SIMD<double> b110(pb2);
SIMD<double> b111(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0 += SIMD<double>(pa2[10]) * b100;
sum1 += SIMD<double>(pa2[10]) * b101;
sum0 += SIMD<double>(pa2[11]) * b110;
sum1 += SIMD<double>(pa2[11]) * b111;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW,mask); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW,mask); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW,mask); pb2 += db;
SIMD<double> b100(pb2);
SIMD<double> b101(pb2+SW,mask); pb2 += db;
SIMD<double> b110(pb2);
SIMD<double> b111(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0 += SIMD<double>(pa2[10]) * b100;
sum1 += SIMD<double>(pa2[10]) * b101;
sum0 += SIMD<double>(pa2[11]) * b110;
sum1 += SIMD<double>(pa2[11]) * b111;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
SIMD<double> b10(pb2); pb2 += db;
SIMD<double> b11(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum += SIMD<double>(pa2[11]) * b11;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
SIMD<double> b10(pb2, mask); pb2 += db;
SIMD<double> b11(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum += SIMD<double>(pa2[11]) * b11;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum<12, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
SIMD<double> b10(pb2); pb2 += db;
SIMD<double> b11(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum += SIMD<double>(pa2[11]) * b11;
sum.Store(pc2);
} }
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
SIMD<double> b10(pb2, mask); pb2 += db;
SIMD<double> b11(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum(pc2, mask);
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum += SIMD<double>(pa2[11]) * b11;
sum.Store(pc2, mask);
} }
template <> INLINE void MatKernelShortSum2<12, ADD>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
{
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW); pb2 += db;
SIMD<double> b100(pb2);
SIMD<double> b101(pb2+SW); pb2 += db;
SIMD<double> b110(pb2);
SIMD<double> b111(pb2+SW); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0 += SIMD<double>(pa2[10]) * b100;
sum1 += SIMD<double>(pa2[10]) * b101;
sum0 += SIMD<double>(pa2[11]) * b110;
sum1 += SIMD<double>(pa2[11]) * b111;
sum0.Store(pc2);
sum1.Store(pc2+SW);
} }
size_t rest = wb % (2*SW); 
if (rest == 0) return; 
if (rest >= SW) 
{
if (rest > SW)
{
SIMD<mask64> mask(rest-SW); 
double * pb2 = pb;
SIMD<double> b00(pb2);
SIMD<double> b01(pb2+SW,mask); pb2 += db;
SIMD<double> b10(pb2);
SIMD<double> b11(pb2+SW,mask); pb2 += db;
SIMD<double> b20(pb2);
SIMD<double> b21(pb2+SW,mask); pb2 += db;
SIMD<double> b30(pb2);
SIMD<double> b31(pb2+SW,mask); pb2 += db;
SIMD<double> b40(pb2);
SIMD<double> b41(pb2+SW,mask); pb2 += db;
SIMD<double> b50(pb2);
SIMD<double> b51(pb2+SW,mask); pb2 += db;
SIMD<double> b60(pb2);
SIMD<double> b61(pb2+SW,mask); pb2 += db;
SIMD<double> b70(pb2);
SIMD<double> b71(pb2+SW,mask); pb2 += db;
SIMD<double> b80(pb2);
SIMD<double> b81(pb2+SW,mask); pb2 += db;
SIMD<double> b90(pb2);
SIMD<double> b91(pb2+SW,mask); pb2 += db;
SIMD<double> b100(pb2);
SIMD<double> b101(pb2+SW,mask); pb2 += db;
SIMD<double> b110(pb2);
SIMD<double> b111(pb2+SW,mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum0 = 0.0;
SIMD<double> sum1 = 0.0;
sum0 += SIMD<double>(pa2[0]) * b00;
sum1 += SIMD<double>(pa2[0]) * b01;
sum0 += SIMD<double>(pa2[1]) * b10;
sum1 += SIMD<double>(pa2[1]) * b11;
sum0 += SIMD<double>(pa2[2]) * b20;
sum1 += SIMD<double>(pa2[2]) * b21;
sum0 += SIMD<double>(pa2[3]) * b30;
sum1 += SIMD<double>(pa2[3]) * b31;
sum0 += SIMD<double>(pa2[4]) * b40;
sum1 += SIMD<double>(pa2[4]) * b41;
sum0 += SIMD<double>(pa2[5]) * b50;
sum1 += SIMD<double>(pa2[5]) * b51;
sum0 += SIMD<double>(pa2[6]) * b60;
sum1 += SIMD<double>(pa2[6]) * b61;
sum0 += SIMD<double>(pa2[7]) * b70;
sum1 += SIMD<double>(pa2[7]) * b71;
sum0 += SIMD<double>(pa2[8]) * b80;
sum1 += SIMD<double>(pa2[8]) * b81;
sum0 += SIMD<double>(pa2[9]) * b90;
sum1 += SIMD<double>(pa2[9]) * b91;
sum0 += SIMD<double>(pa2[10]) * b100;
sum1 += SIMD<double>(pa2[10]) * b101;
sum0 += SIMD<double>(pa2[11]) * b110;
sum1 += SIMD<double>(pa2[11]) * b111;
sum0.Store(pc2);
sum1.Store(pc2+SW,mask);
}
return;
}
double * pb2 = pb;
SIMD<double> b0(pb2); pb2 += db;
SIMD<double> b1(pb2); pb2 += db;
SIMD<double> b2(pb2); pb2 += db;
SIMD<double> b3(pb2); pb2 += db;
SIMD<double> b4(pb2); pb2 += db;
SIMD<double> b5(pb2); pb2 += db;
SIMD<double> b6(pb2); pb2 += db;
SIMD<double> b7(pb2); pb2 += db;
SIMD<double> b8(pb2); pb2 += db;
SIMD<double> b9(pb2); pb2 += db;
SIMD<double> b10(pb2); pb2 += db;
SIMD<double> b11(pb2); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum += SIMD<double>(pa2[11]) * b11;
sum.Store(pc2);
}
return;
}
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> b0(pb2, mask); pb2 += db;
SIMD<double> b1(pb2, mask); pb2 += db;
SIMD<double> b2(pb2, mask); pb2 += db;
SIMD<double> b3(pb2, mask); pb2 += db;
SIMD<double> b4(pb2, mask); pb2 += db;
SIMD<double> b5(pb2, mask); pb2 += db;
SIMD<double> b6(pb2, mask); pb2 += db;
SIMD<double> b7(pb2, mask); pb2 += db;
SIMD<double> b8(pb2, mask); pb2 += db;
SIMD<double> b9(pb2, mask); pb2 += db;
SIMD<double> b10(pb2, mask); pb2 += db;
SIMD<double> b11(pb2, mask); pb2 += db;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
{
SIMD<double> sum = 0.0;
sum += SIMD<double>(pa2[0]) * b0;
sum += SIMD<double>(pa2[1]) * b1;
sum += SIMD<double>(pa2[2]) * b2;
sum += SIMD<double>(pa2[3]) * b3;
sum += SIMD<double>(pa2[4]) * b4;
sum += SIMD<double>(pa2[5]) * b5;
sum += SIMD<double>(pa2[6]) * b6;
sum += SIMD<double>(pa2[7]) * b7;
sum += SIMD<double>(pa2[8]) * b8;
sum += SIMD<double>(pa2[9]) * b9;
sum += SIMD<double>(pa2[10]) * b10;
sum += SIMD<double>(pa2[11]) * b11;
sum.Store(pc2, mask);
} }
// C = A^t * B,  with short inner loop
template <size_t WA, OPERATION OP>
inline void MatKernelAtB_SmallWA
(size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
// C = A^t * B,  with short inner loop, unroll width B
template <size_t WA, OPERATION OP>
inline void MatKernelAtB_SmallWA
(size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);
template <> INLINE void MatKernelAtB_SmallWA<0, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
}
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
}
}
template <> INLINE void MatKernelAtB_SmallWA<1, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
}
sum0.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
}
sum0.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<2, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<3, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
sum2.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
sum2.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<4, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
sum2.Store(pc2); pc2 += dc;
sum3.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
sum2.Store(pc2, mask); pc2 += dc;
sum3.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<5, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
sum2.Store(pc2); pc2 += dc;
sum3.Store(pc2); pc2 += dc;
sum4.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
sum2.Store(pc2, mask); pc2 += dc;
sum3.Store(pc2, mask); pc2 += dc;
sum4.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<6, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
sum2.Store(pc2); pc2 += dc;
sum3.Store(pc2); pc2 += dc;
sum4.Store(pc2); pc2 += dc;
sum5.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
sum2.Store(pc2, mask); pc2 += dc;
sum3.Store(pc2, mask); pc2 += dc;
sum4.Store(pc2, mask); pc2 += dc;
sum5.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<7, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
sum2.Store(pc2); pc2 += dc;
sum3.Store(pc2); pc2 += dc;
sum4.Store(pc2); pc2 += dc;
sum5.Store(pc2); pc2 += dc;
sum6.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
sum2.Store(pc2, mask); pc2 += dc;
sum3.Store(pc2, mask); pc2 += dc;
sum4.Store(pc2, mask); pc2 += dc;
sum5.Store(pc2, mask); pc2 += dc;
sum6.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<8, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
SIMD<double> sum7(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
FMAasm (bjk,SIMD<double>(pa2[7]), sum7);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
sum2.Store(pc2); pc2 += dc;
sum3.Store(pc2); pc2 += dc;
sum4.Store(pc2); pc2 += dc;
sum5.Store(pc2); pc2 += dc;
sum6.Store(pc2); pc2 += dc;
sum7.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
SIMD<double> sum7(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
FMAasm (bjk,SIMD<double>(pa2[7]), sum7);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
sum2.Store(pc2, mask); pc2 += dc;
sum3.Store(pc2, mask); pc2 += dc;
sum4.Store(pc2, mask); pc2 += dc;
sum5.Store(pc2, mask); pc2 += dc;
sum6.Store(pc2, mask); pc2 += dc;
sum7.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<9, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
SIMD<double> sum7(0.0);
SIMD<double> sum8(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
FMAasm (bjk,SIMD<double>(pa2[7]), sum7);
FMAasm (bjk,SIMD<double>(pa2[8]), sum8);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
sum2.Store(pc2); pc2 += dc;
sum3.Store(pc2); pc2 += dc;
sum4.Store(pc2); pc2 += dc;
sum5.Store(pc2); pc2 += dc;
sum6.Store(pc2); pc2 += dc;
sum7.Store(pc2); pc2 += dc;
sum8.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
SIMD<double> sum7(0.0);
SIMD<double> sum8(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
FMAasm (bjk,SIMD<double>(pa2[7]), sum7);
FMAasm (bjk,SIMD<double>(pa2[8]), sum8);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
sum2.Store(pc2, mask); pc2 += dc;
sum3.Store(pc2, mask); pc2 += dc;
sum4.Store(pc2, mask); pc2 += dc;
sum5.Store(pc2, mask); pc2 += dc;
sum6.Store(pc2, mask); pc2 += dc;
sum7.Store(pc2, mask); pc2 += dc;
sum8.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<10, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
SIMD<double> sum7(0.0);
SIMD<double> sum8(0.0);
SIMD<double> sum9(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
FMAasm (bjk,SIMD<double>(pa2[7]), sum7);
FMAasm (bjk,SIMD<double>(pa2[8]), sum8);
FMAasm (bjk,SIMD<double>(pa2[9]), sum9);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
sum2.Store(pc2); pc2 += dc;
sum3.Store(pc2); pc2 += dc;
sum4.Store(pc2); pc2 += dc;
sum5.Store(pc2); pc2 += dc;
sum6.Store(pc2); pc2 += dc;
sum7.Store(pc2); pc2 += dc;
sum8.Store(pc2); pc2 += dc;
sum9.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
SIMD<double> sum7(0.0);
SIMD<double> sum8(0.0);
SIMD<double> sum9(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
FMAasm (bjk,SIMD<double>(pa2[7]), sum7);
FMAasm (bjk,SIMD<double>(pa2[8]), sum8);
FMAasm (bjk,SIMD<double>(pa2[9]), sum9);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
sum2.Store(pc2, mask); pc2 += dc;
sum3.Store(pc2, mask); pc2 += dc;
sum4.Store(pc2, mask); pc2 += dc;
sum5.Store(pc2, mask); pc2 += dc;
sum6.Store(pc2, mask); pc2 += dc;
sum7.Store(pc2, mask); pc2 += dc;
sum8.Store(pc2, mask); pc2 += dc;
sum9.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<11, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
SIMD<double> sum7(0.0);
SIMD<double> sum8(0.0);
SIMD<double> sum9(0.0);
SIMD<double> sum10(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
FMAasm (bjk,SIMD<double>(pa2[7]), sum7);
FMAasm (bjk,SIMD<double>(pa2[8]), sum8);
FMAasm (bjk,SIMD<double>(pa2[9]), sum9);
FMAasm (bjk,SIMD<double>(pa2[10]), sum10);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
sum2.Store(pc2); pc2 += dc;
sum3.Store(pc2); pc2 += dc;
sum4.Store(pc2); pc2 += dc;
sum5.Store(pc2); pc2 += dc;
sum6.Store(pc2); pc2 += dc;
sum7.Store(pc2); pc2 += dc;
sum8.Store(pc2); pc2 += dc;
sum9.Store(pc2); pc2 += dc;
sum10.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
SIMD<double> sum7(0.0);
SIMD<double> sum8(0.0);
SIMD<double> sum9(0.0);
SIMD<double> sum10(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
FMAasm (bjk,SIMD<double>(pa2[7]), sum7);
FMAasm (bjk,SIMD<double>(pa2[8]), sum8);
FMAasm (bjk,SIMD<double>(pa2[9]), sum9);
FMAasm (bjk,SIMD<double>(pa2[10]), sum10);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
sum2.Store(pc2, mask); pc2 += dc;
sum3.Store(pc2, mask); pc2 += dc;
sum4.Store(pc2, mask); pc2 += dc;
sum5.Store(pc2, mask); pc2 += dc;
sum6.Store(pc2, mask); pc2 += dc;
sum7.Store(pc2, mask); pc2 += dc;
sum8.Store(pc2, mask); pc2 += dc;
sum9.Store(pc2, mask); pc2 += dc;
sum10.Store(pc2, mask); pc2 += dc;
}
template <> INLINE void MatKernelAtB_SmallWA<12, SET>
    (size_t ha, size_t wb,
     double * pa, size_t da,
     double * pb, size_t db,
     double * pc, size_t dc)
{
constexpr int SW = SIMD<double>::Size();
for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
{
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
SIMD<double> sum7(0.0);
SIMD<double> sum8(0.0);
SIMD<double> sum9(0.0);
SIMD<double> sum10(0.0);
SIMD<double> sum11(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
FMAasm (bjk,SIMD<double>(pa2[7]), sum7);
FMAasm (bjk,SIMD<double>(pa2[8]), sum8);
FMAasm (bjk,SIMD<double>(pa2[9]), sum9);
FMAasm (bjk,SIMD<double>(pa2[10]), sum10);
FMAasm (bjk,SIMD<double>(pa2[11]), sum11);
}
sum0.Store(pc2); pc2 += dc;
sum1.Store(pc2); pc2 += dc;
sum2.Store(pc2); pc2 += dc;
sum3.Store(pc2); pc2 += dc;
sum4.Store(pc2); pc2 += dc;
sum5.Store(pc2); pc2 += dc;
sum6.Store(pc2); pc2 += dc;
sum7.Store(pc2); pc2 += dc;
sum8.Store(pc2); pc2 += dc;
sum9.Store(pc2); pc2 += dc;
sum10.Store(pc2); pc2 += dc;
sum11.Store(pc2); pc2 += dc;
}
size_t rest = wb % SW; 
if (rest == 0) return; 
SIMD<mask64> mask(rest); 
double * pb2 = pb;
SIMD<double> sum0(0.0);
SIMD<double> sum1(0.0);
SIMD<double> sum2(0.0);
SIMD<double> sum3(0.0);
SIMD<double> sum4(0.0);
SIMD<double> sum5(0.0);
SIMD<double> sum6(0.0);
SIMD<double> sum7(0.0);
SIMD<double> sum8(0.0);
SIMD<double> sum9(0.0);
SIMD<double> sum10(0.0);
SIMD<double> sum11(0.0);
double * pa2 = pa;
double * pc2 = pc;
__assume(ha>0);
#pragma unroll 1
for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
{
SIMD<double> bjk(pb2, mask);
FMAasm (bjk,SIMD<double>(pa2[0]), sum0);
FMAasm (bjk,SIMD<double>(pa2[1]), sum1);
FMAasm (bjk,SIMD<double>(pa2[2]), sum2);
FMAasm (bjk,SIMD<double>(pa2[3]), sum3);
FMAasm (bjk,SIMD<double>(pa2[4]), sum4);
FMAasm (bjk,SIMD<double>(pa2[5]), sum5);
FMAasm (bjk,SIMD<double>(pa2[6]), sum6);
FMAasm (bjk,SIMD<double>(pa2[7]), sum7);
FMAasm (bjk,SIMD<double>(pa2[8]), sum8);
FMAasm (bjk,SIMD<double>(pa2[9]), sum9);
FMAasm (bjk,SIMD<double>(pa2[10]), sum10);
FMAasm (bjk,SIMD<double>(pa2[11]), sum11);
}
sum0.Store(pc2, mask); pc2 += dc;
sum1.Store(pc2, mask); pc2 += dc;
sum2.Store(pc2, mask); pc2 += dc;
sum3.Store(pc2, mask); pc2 += dc;
sum4.Store(pc2, mask); pc2 += dc;
sum5.Store(pc2, mask); pc2 += dc;
sum6.Store(pc2, mask); pc2 += dc;
sum7.Store(pc2, mask); pc2 += dc;
sum8.Store(pc2, mask); pc2 += dc;
sum9.Store(pc2, mask); pc2 += dc;
sum10.Store(pc2, mask); pc2 += dc;
sum11.Store(pc2, mask); pc2 += dc;
}
// y = A * x,  with fix width
template <size_t WA, OPERATION OP>
inline void KernelMatVec
(size_t ha, double * pa, size_t da, double * x, double * y);
template <> INLINE void KernelMatVec<0, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<1, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<mask64> mask(1);
SIMD<double> x0(x+0, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0, mask) * x0;
sum1 += SIMD<double>(pa+da+0, mask) * x0;
sum2 += SIMD<double>(pa+2*da+0, mask) * x0;
sum3 += SIMD<double>(pa+3*da+0, mask) * x0;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0, mask) * x0;
sum1 += SIMD<double>(pa+da+0, mask) * x0;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0, mask) * x0;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<2, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<mask64> mask(2);
SIMD<double> x0(x+0, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0, mask) * x0;
sum1 += SIMD<double>(pa+da+0, mask) * x0;
sum2 += SIMD<double>(pa+2*da+0, mask) * x0;
sum3 += SIMD<double>(pa+3*da+0, mask) * x0;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0, mask) * x0;
sum1 += SIMD<double>(pa+da+0, mask) * x0;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0, mask) * x0;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<3, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<mask64> mask(3);
SIMD<double> x0(x+0, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0, mask) * x0;
sum1 += SIMD<double>(pa+da+0, mask) * x0;
sum2 += SIMD<double>(pa+2*da+0, mask) * x0;
sum3 += SIMD<double>(pa+3*da+0, mask) * x0;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0, mask) * x0;
sum1 += SIMD<double>(pa+da+0, mask) * x0;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0, mask) * x0;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<4, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<5, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<mask64> mask(1);
SIMD<double> x1(x+4, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4, mask) * x1;
sum1 += SIMD<double>(pa+da+4, mask) * x1;
sum2 += SIMD<double>(pa+2*da+4, mask) * x1;
sum3 += SIMD<double>(pa+3*da+4, mask) * x1;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4, mask) * x1;
sum1 += SIMD<double>(pa+da+4, mask) * x1;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4, mask) * x1;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<6, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<mask64> mask(2);
SIMD<double> x1(x+4, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4, mask) * x1;
sum1 += SIMD<double>(pa+da+4, mask) * x1;
sum2 += SIMD<double>(pa+2*da+4, mask) * x1;
sum3 += SIMD<double>(pa+3*da+4, mask) * x1;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4, mask) * x1;
sum1 += SIMD<double>(pa+da+4, mask) * x1;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4, mask) * x1;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<7, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<mask64> mask(3);
SIMD<double> x1(x+4, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4, mask) * x1;
sum1 += SIMD<double>(pa+da+4, mask) * x1;
sum2 += SIMD<double>(pa+2*da+4, mask) * x1;
sum3 += SIMD<double>(pa+3*da+4, mask) * x1;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4, mask) * x1;
sum1 += SIMD<double>(pa+da+4, mask) * x1;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4, mask) * x1;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<8, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<9, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<mask64> mask(1);
SIMD<double> x2(x+8, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8, mask) * x2;
sum1 += SIMD<double>(pa+da+8, mask) * x2;
sum2 += SIMD<double>(pa+2*da+8, mask) * x2;
sum3 += SIMD<double>(pa+3*da+8, mask) * x2;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8, mask) * x2;
sum1 += SIMD<double>(pa+da+8, mask) * x2;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8, mask) * x2;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<10, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<mask64> mask(2);
SIMD<double> x2(x+8, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8, mask) * x2;
sum1 += SIMD<double>(pa+da+8, mask) * x2;
sum2 += SIMD<double>(pa+2*da+8, mask) * x2;
sum3 += SIMD<double>(pa+3*da+8, mask) * x2;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8, mask) * x2;
sum1 += SIMD<double>(pa+da+8, mask) * x2;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8, mask) * x2;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<11, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<mask64> mask(3);
SIMD<double> x2(x+8, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8, mask) * x2;
sum1 += SIMD<double>(pa+da+8, mask) * x2;
sum2 += SIMD<double>(pa+2*da+8, mask) * x2;
sum3 += SIMD<double>(pa+3*da+8, mask) * x2;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8, mask) * x2;
sum1 += SIMD<double>(pa+da+8, mask) * x2;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8, mask) * x2;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<12, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<13, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<mask64> mask(1);
SIMD<double> x3(x+12, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12, mask) * x3;
sum1 += SIMD<double>(pa+da+12, mask) * x3;
sum2 += SIMD<double>(pa+2*da+12, mask) * x3;
sum3 += SIMD<double>(pa+3*da+12, mask) * x3;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12, mask) * x3;
sum1 += SIMD<double>(pa+da+12, mask) * x3;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12, mask) * x3;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<14, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<mask64> mask(2);
SIMD<double> x3(x+12, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12, mask) * x3;
sum1 += SIMD<double>(pa+da+12, mask) * x3;
sum2 += SIMD<double>(pa+2*da+12, mask) * x3;
sum3 += SIMD<double>(pa+3*da+12, mask) * x3;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12, mask) * x3;
sum1 += SIMD<double>(pa+da+12, mask) * x3;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12, mask) * x3;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<15, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<mask64> mask(3);
SIMD<double> x3(x+12, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12, mask) * x3;
sum1 += SIMD<double>(pa+da+12, mask) * x3;
sum2 += SIMD<double>(pa+2*da+12, mask) * x3;
sum3 += SIMD<double>(pa+3*da+12, mask) * x3;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12, mask) * x3;
sum1 += SIMD<double>(pa+da+12, mask) * x3;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12, mask) * x3;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<16, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<double> x3(x+12);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum2 += SIMD<double>(pa+2*da+12) * x3;
sum3 += SIMD<double>(pa+3*da+12) * x3;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12) * x3;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<17, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<double> x3(x+12);
SIMD<mask64> mask(1);
SIMD<double> x4(x+16, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum2 += SIMD<double>(pa+2*da+12) * x3;
sum3 += SIMD<double>(pa+3*da+12) * x3;
sum0 += SIMD<double>(pa+16, mask) * x4;
sum1 += SIMD<double>(pa+da+16, mask) * x4;
sum2 += SIMD<double>(pa+2*da+16, mask) * x4;
sum3 += SIMD<double>(pa+3*da+16, mask) * x4;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum0 += SIMD<double>(pa+16, mask) * x4;
sum1 += SIMD<double>(pa+da+16, mask) * x4;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12) * x3;
sum += SIMD<double>(pa+16, mask) * x4;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<18, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<double> x3(x+12);
SIMD<mask64> mask(2);
SIMD<double> x4(x+16, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum2 += SIMD<double>(pa+2*da+12) * x3;
sum3 += SIMD<double>(pa+3*da+12) * x3;
sum0 += SIMD<double>(pa+16, mask) * x4;
sum1 += SIMD<double>(pa+da+16, mask) * x4;
sum2 += SIMD<double>(pa+2*da+16, mask) * x4;
sum3 += SIMD<double>(pa+3*da+16, mask) * x4;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum0 += SIMD<double>(pa+16, mask) * x4;
sum1 += SIMD<double>(pa+da+16, mask) * x4;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12) * x3;
sum += SIMD<double>(pa+16, mask) * x4;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<19, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<double> x3(x+12);
SIMD<mask64> mask(3);
SIMD<double> x4(x+16, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum2 += SIMD<double>(pa+2*da+12) * x3;
sum3 += SIMD<double>(pa+3*da+12) * x3;
sum0 += SIMD<double>(pa+16, mask) * x4;
sum1 += SIMD<double>(pa+da+16, mask) * x4;
sum2 += SIMD<double>(pa+2*da+16, mask) * x4;
sum3 += SIMD<double>(pa+3*da+16, mask) * x4;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum0 += SIMD<double>(pa+16, mask) * x4;
sum1 += SIMD<double>(pa+da+16, mask) * x4;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12) * x3;
sum += SIMD<double>(pa+16, mask) * x4;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<20, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<double> x3(x+12);
SIMD<double> x4(x+16);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum2 += SIMD<double>(pa+2*da+12) * x3;
sum3 += SIMD<double>(pa+3*da+12) * x3;
sum0 += SIMD<double>(pa+16) * x4;
sum1 += SIMD<double>(pa+da+16) * x4;
sum2 += SIMD<double>(pa+2*da+16) * x4;
sum3 += SIMD<double>(pa+3*da+16) * x4;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum0 += SIMD<double>(pa+16) * x4;
sum1 += SIMD<double>(pa+da+16) * x4;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12) * x3;
sum += SIMD<double>(pa+16) * x4;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<21, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<double> x3(x+12);
SIMD<double> x4(x+16);
SIMD<mask64> mask(1);
SIMD<double> x5(x+20, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum2 += SIMD<double>(pa+2*da+12) * x3;
sum3 += SIMD<double>(pa+3*da+12) * x3;
sum0 += SIMD<double>(pa+16) * x4;
sum1 += SIMD<double>(pa+da+16) * x4;
sum2 += SIMD<double>(pa+2*da+16) * x4;
sum3 += SIMD<double>(pa+3*da+16) * x4;
sum0 += SIMD<double>(pa+20, mask) * x5;
sum1 += SIMD<double>(pa+da+20, mask) * x5;
sum2 += SIMD<double>(pa+2*da+20, mask) * x5;
sum3 += SIMD<double>(pa+3*da+20, mask) * x5;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum0 += SIMD<double>(pa+16) * x4;
sum1 += SIMD<double>(pa+da+16) * x4;
sum0 += SIMD<double>(pa+20, mask) * x5;
sum1 += SIMD<double>(pa+da+20, mask) * x5;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12) * x3;
sum += SIMD<double>(pa+16) * x4;
sum += SIMD<double>(pa+20, mask) * x5;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<22, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<double> x3(x+12);
SIMD<double> x4(x+16);
SIMD<mask64> mask(2);
SIMD<double> x5(x+20, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum2 += SIMD<double>(pa+2*da+12) * x3;
sum3 += SIMD<double>(pa+3*da+12) * x3;
sum0 += SIMD<double>(pa+16) * x4;
sum1 += SIMD<double>(pa+da+16) * x4;
sum2 += SIMD<double>(pa+2*da+16) * x4;
sum3 += SIMD<double>(pa+3*da+16) * x4;
sum0 += SIMD<double>(pa+20, mask) * x5;
sum1 += SIMD<double>(pa+da+20, mask) * x5;
sum2 += SIMD<double>(pa+2*da+20, mask) * x5;
sum3 += SIMD<double>(pa+3*da+20, mask) * x5;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum0 += SIMD<double>(pa+16) * x4;
sum1 += SIMD<double>(pa+da+16) * x4;
sum0 += SIMD<double>(pa+20, mask) * x5;
sum1 += SIMD<double>(pa+da+20, mask) * x5;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12) * x3;
sum += SIMD<double>(pa+16) * x4;
sum += SIMD<double>(pa+20, mask) * x5;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<23, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<double> x3(x+12);
SIMD<double> x4(x+16);
SIMD<mask64> mask(3);
SIMD<double> x5(x+20, mask);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum2 += SIMD<double>(pa+2*da+12) * x3;
sum3 += SIMD<double>(pa+3*da+12) * x3;
sum0 += SIMD<double>(pa+16) * x4;
sum1 += SIMD<double>(pa+da+16) * x4;
sum2 += SIMD<double>(pa+2*da+16) * x4;
sum3 += SIMD<double>(pa+3*da+16) * x4;
sum0 += SIMD<double>(pa+20, mask) * x5;
sum1 += SIMD<double>(pa+da+20, mask) * x5;
sum2 += SIMD<double>(pa+2*da+20, mask) * x5;
sum3 += SIMD<double>(pa+3*da+20, mask) * x5;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum0 += SIMD<double>(pa+16) * x4;
sum1 += SIMD<double>(pa+da+16) * x4;
sum0 += SIMD<double>(pa+20, mask) * x5;
sum1 += SIMD<double>(pa+da+20, mask) * x5;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12) * x3;
sum += SIMD<double>(pa+16) * x4;
sum += SIMD<double>(pa+20, mask) * x5;
y[i] = HSum(sum);
} }
template <> INLINE void KernelMatVec<24, SET>
(size_t ha, double * pa, size_t da, double * x, double * y) {
constexpr int SW = SIMD<double>::Size();
SIMD<double> x0(x+0);
SIMD<double> x1(x+4);
SIMD<double> x2(x+8);
SIMD<double> x3(x+12);
SIMD<double> x4(x+16);
SIMD<double> x5(x+20);
size_t i = 0;
for ( ; i+4 <= ha; i+=4, pa += 4*da) {
SIMD<double> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum2 += SIMD<double>(pa+2*da+0) * x0;
sum3 += SIMD<double>(pa+3*da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum2 += SIMD<double>(pa+2*da+4) * x1;
sum3 += SIMD<double>(pa+3*da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum2 += SIMD<double>(pa+2*da+8) * x2;
sum3 += SIMD<double>(pa+3*da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum2 += SIMD<double>(pa+2*da+12) * x3;
sum3 += SIMD<double>(pa+3*da+12) * x3;
sum0 += SIMD<double>(pa+16) * x4;
sum1 += SIMD<double>(pa+da+16) * x4;
sum2 += SIMD<double>(pa+2*da+16) * x4;
sum3 += SIMD<double>(pa+3*da+16) * x4;
sum0 += SIMD<double>(pa+20) * x5;
sum1 += SIMD<double>(pa+da+20) * x5;
sum2 += SIMD<double>(pa+2*da+20) * x5;
sum3 += SIMD<double>(pa+3*da+20) * x5;
SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);
vsum.Store(y+i);
}
if (ha & 2) {
SIMD<double> sum0(0.0), sum1(0.0);
sum0 += SIMD<double>(pa+0) * x0;
sum1 += SIMD<double>(pa+da+0) * x0;
sum0 += SIMD<double>(pa+4) * x1;
sum1 += SIMD<double>(pa+da+4) * x1;
sum0 += SIMD<double>(pa+8) * x2;
sum1 += SIMD<double>(pa+da+8) * x2;
sum0 += SIMD<double>(pa+12) * x3;
sum1 += SIMD<double>(pa+da+12) * x3;
sum0 += SIMD<double>(pa+16) * x4;
sum1 += SIMD<double>(pa+da+16) * x4;
sum0 += SIMD<double>(pa+20) * x5;
sum1 += SIMD<double>(pa+da+20) * x5;
SIMD<double,2> vsum = HSum(sum0,sum1);
vsum.Store(y+i);
i += 2; pa += 2*da;
}
if (ha & 1) {
SIMD<double> sum(0.0);
sum += SIMD<double>(pa+0) * x0;
sum += SIMD<double>(pa+4) * x1;
sum += SIMD<double>(pa+8) * x2;
sum += SIMD<double>(pa+12) * x3;
sum += SIMD<double>(pa+16) * x4;
sum += SIMD<double>(pa+20) * x5;
y[i] = HSum(sum);
} }
