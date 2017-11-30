// ngscxx generate_mat_kernels.cpp ; a.out

#include <ngstd.hpp>
using namespace ngstd;

enum OP { ADD, SUB, SET };

string ToString (OP op)
{
  switch (op)
    {
    case SET: return "SET";
    case ADD: return "ADD";
    case SUB: return "SUB";
    }
}

/*
  C = A * B
  C += A * B
  C -= A * B

  A ... h x n
  B ... n x w*SIMD.Size
 */
void GenerateMultAB (ostream & out, int h, int w, OP op, bool aligned_b)
{
  out << "template <> INLINE void MatKernelMultAB<" << h << ", " << w << ", " << ToString(op) << ">" << endl;
  out << "    (size_t n," << endl
      << "     double * pa, size_t da," << endl
      << "     " << (aligned_b ? "SIMD<double>" : "double") << " * pb, size_t db," << endl
      << "     double * pc, size_t dc)" << endl
      << "{" << endl;
  out << "constexpr int SW = SIMD<double>::Size();" << endl
      << "double * hpc = pc;" << endl;

  if (op == SET)
    {
      for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++)
          out << "SIMD<double> sum" << i << j << "(0);" << endl;
    }
  else
    {
      for (int i = 0; i < h; i++)
        {
          for (int j = 0; j < w; j++)
            out << "SIMD<double> sum" << i << j << "(pc+SW*" << j << ");" << endl;
          out << "pc += dc;" << endl;
        }
      out << "pc = hpc;" << endl;
    }
  
  out << "for (size_t i = 0; i < n; i++, pa++, pb += db) {" << endl;
  if (aligned_b)
    for (int i = 0; i < w; i++)
      out << "SIMD<double> b" << i << " = pb[" << i << "];" << endl;
  else
    for (int i = 0; i < w; i++)
      out << "SIMD<double> b" << i << "(pb+" << i << "*SW);" << endl;

  for (int i = 0; i < h; i++)
    {
      out << "SIMD<double> a" << i << "(pa["<< i << "*da]);" << endl;
      for (int j = 0; j < w; j++)
        // out << "sum" << j << i << " += a" << j << " * b" << i << ";" << endl;
        out << "FMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;
    }
  out << "}" << endl;

  for (int i = 0; i < h; i++)
    {
      for (int j = 0; j < w; j++)
        out << "sum"<< i << j << ".Store(pc+SW*" << j << ");" << endl;
      out << "pc += dc;" << endl;
    }
  
  out << "}" << endl;
}


void GenerateMultAB (ostream & out, int h, int w)
{
  GenerateMultAB (out, h, w, SET, false);
  GenerateMultAB (out, h, w, ADD, false);
  GenerateMultAB (out, h, w, SET, true);
  GenerateMultAB (out, h, w, ADD, true);
}





/*
  C = A * B
  C += A * B
  C -= A * B

  A ... h x n
  B ... n x w*SIMD.Size
 */
void AlignedGenerateMultAB (ostream & out, int h, int w, OP op)
{
  
  out << "template <> void MatKernelAlignedMultAB<" << h << ", " << w << ">" << endl
      << "    (size_t n," << endl
      << "     double * pa, size_t da," << endl
      << "     SIMD<double> * pb, size_t db," << endl
      << "     SIMD<double> * pc, size_t dc)" << endl
      << "{" << endl;

  out << "SIMD<double> * hpc = pc;" << endl;

  if (op == SET)
    {
      for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++)
          out << "SIMD<double> sum" << i << j << "(0);" << endl;
    }
  else
    {
      for (int i = 0; i < h; i++)
        {
          for (int j = 0; j < w; j++)
            out << "SIMD<double> sum" << i << j << "(pc+" << j << ");" << endl;
          out << "pc += dc;" << endl;
        }
      out << "pc = hpc;" << endl;
    }
  
  out << "for (size_t i = 0; i < n; i++, pa++, pb += db) {" << endl;
  for (int i = 0; i < w; i++)
    out << "SIMD<double> b" << i << "(pb[" << i << "]);" << endl;

  for (int i = 0; i < h; i++)
    {
      out << "SIMD<double> a" << i << "(pa["<< i << "*da]);" << endl;
      for (int j = 0; j < w; j++)
        // out << "sum" << j << i << " += a" << j << " * b" << i << ";" << endl;
        out << "FMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;
    }
  out << "}" << endl;

  for (int i = 0; i < h; i++)
    {
      for (int j = 0; j < w; j++)
        // out << "sum"<< i << j << ".Store(pc+" << j << ");" << endl;
        out << "pc[" << j << "]= sum"  << i << j << ";" << endl;
      out << "pc += dc;" << endl;
    }
  
  out << "}" << endl;
}




void GenerateMultABMask (ostream & out, int h, OP op, bool aligned_b)
{
  out << "template <> void MatKernelMultABMask<" << h << ", " << ToString(op) << ">" << endl;
    
  out << "    (size_t n, SIMD<mask64> mask," << endl
      << "     double * pa, size_t da," << endl
      << "     " << (aligned_b ? "SIMD<double>" : "double") << " * pb, size_t db," << endl    
      << "     double * pc, size_t dc)" << endl
      << "{" << endl;
  out << "constexpr int SW = SIMD<double>::Size();" << endl
      << "double * hpc = pc;" << endl;

  if (op == SET)
    {
      for (int i = 0; i < h; i++)
        out << "SIMD<double> sum" << i << "(0);" << endl;
    }
  else
    {
      for (int i = 0; i < h; i++)
        {
          out << "SIMD<double> sum" << i << "(pc, mask);" << endl;
          out << "pc += dc;" << endl;
        }
      out << "pc = hpc;" << endl;
    }
  
  out << "for (size_t i = 0; i < n; i++, pa++, pb += db) {" << endl;
  out << "SIMD<double> b((double*)pb,mask);" << endl;

  for (int i = 0; i < h; i++)
    {
      out << "SIMD<double> a" << i << "(pa["<< i << "*da]);" << endl;
      out << "FMAasm(a"<<i<<",b,sum" << i << ");" << endl;
    }
  out << "}" << endl;

  for (int i = 0; i < h; i++)
    {
      out << "sum"<< i << ".Store(pc,mask);" << endl;
      out << "pc += dc;" << endl;
    }
  
  out << "}" << endl;
}

void GenerateMultABMask (ostream & out, int h)
{
  GenerateMultABMask (out, h, SET, false);
  GenerateMultABMask (out, h, ADD, false);
  GenerateMultABMask (out, h, SET, true);
  GenerateMultABMask (out, h, ADD, true);
}


/*
  C = A * B^t
  A ... h x n
  B ... w * n
 */
void GenerateScalAB (ostream & out, int h, int w)
{
  out << "template <> INLINE auto MatKernelScalAB<" << h << ", " << w << ">" << endl
      << "    (size_t n," << endl
      << "     double * pa, size_t da," << endl
      << "     double * pb, size_t db)" << endl
      << "{" << endl;
  out << "constexpr int SW = SIMD<double>::Size();" << endl;

  for (int i = 0; i < h; i++)
    for (int j = 0; j < w; j++)
      out << "SIMD<double> sum" << i << j << "(0);" << endl;

  out << "size_t i = 0;" << endl;
  out << "for ( ; i+SW <= n; i+=SW) {" << endl;
  for (int i = 0; i < h; i++)
    out << "SIMD<double> a" << i << "(pa+" << i << "*da+i);" << endl;
  // for (int i = 0; i < w; i++)
  // out << "SIMD<double> b" << i << "(pb+" << i << "*db+i);" << endl;

  for (int j = 0; j < w; j++)
    {
      out << "SIMD<double> b" << j << "(pb+" << j << "*db+i);" << endl;    
      for (int i = 0; i < h; i++)
        // out << "sum" << j << i << " += a" << j << " * b" << i << ";" << endl;
        out << "FMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;
    }
  out << "}" << endl;

  out << "size_t r = n % SW;" << endl;
  out << "if (r) {" << endl;
  out << "SIMD<mask64> mask(r);" << endl;
  for (int i = 0; i < h; i++)
    out << "SIMD<double> a" << i << "(pa+" << i << "*da+i, mask);" << endl;
  for (int i = 0; i < w; i++)
    out << "SIMD<double> b" << i << "(pb+" << i << "*db+i, mask);" << endl;

  for (int i = 0; i < h; i++)
    {
      for (int j = 0; j < w; j++)
        // out << "sum" << j << i << " += a" << j << " * b" << i << ";" << endl;
        out << "FMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;
    }
  out << "}" << endl;
    
  
  out << "return make_tuple(";
  for (int i = 0; i < h; i++)
    {
      out << "HSum(";
      for (int j = 0; j < w; j++)
        {
          out << "sum"<< i << j;
          if (j < w-1)
            out << ",";
          else
            out << ")";
        }
      if (i < h-1)
        out << ",";
      else
        out << ");" << endl;
    }
  
  out << "}" << endl;
}





void GenKernel (ofstream & out, int h, int w)
{
  out << "template <> void MyScalTrans<" << h << ", " << w << ">" << endl
      << "    (size_t n," << endl
      << "     double * pa, size_t da," << endl
      << "     double * pb, size_t db," << endl
      << "     double * pc, size_t dc)" << endl
      << "{" << endl;
  out << "constexpr int SW = SIMD<double>::Size();" << endl
      << "double * hpc = pc;" << endl;
  for (int i = 0; i < h; i++)
    {
      for (int j = 0; j < w; j++)
        out << "SIMD<double> sum" << i << j << "(pc+SW*" << j << ");" << endl;
      out << "pc += dc;" << endl;
    }
  out << "pc = hpc;" << endl;

  out << "for (size_t i = 0; i < n; i++, pa += da, pb += db) {" << endl;
  for (int i = 0; i < h; i++)
    out << "SIMD<double> a" << i << "(pa[" << i << "]);" << endl;

  for (int i = 0; i < w; i++)
    {
      out << "SIMD<double> b" << i << "(pb+" << i << "*SW);" << endl;
      for (int j = 0; j < h; j++)
        // out << "sum" << j << i << " += a" << j << " * b" << i << ";" << endl;
        out << "FMAasm(b"<<i<<",a" << j << ",sum" << j << i << ");" << endl;
    }
  out << "}" << endl;

  for (int i = 0; i < h; i++)
    {
      for (int j = 0; j < w; j++)
        out << "sum"<< i << j << ".Store(pc+SW*" << j << ");" << endl;
      out << "pc += dc;" << endl;
    }
  
  out << "}" << endl;
}


int main ()
{
  ofstream out("matkernel.hpp");

  out << "enum OPERATION { ADD, SUB, SET };" << endl;

  out << "template <size_t H, size_t W, OPERATION OP>" << endl
      << "static void MatKernelMultAB" << endl
      << "(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);" << endl;
  out << "template <size_t H, size_t W, OPERATION OP>" << endl
      << "static void MatKernelMultAB" << endl
      << "(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, double * pc, size_t dc);" << endl;
  out << "template <size_t H, size_t W>" << endl
      << "static void MatKernelAlignedMultAB" << endl
      << "(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, SIMD<double> * pc, size_t dc);" << endl;

  for (int i = 1; i <= 3; i++)
    {
      GenerateMultAB (out, 1, i);  
      GenerateMultAB (out, 2, i);
      GenerateMultAB (out, 3, i);
      GenerateMultAB (out, 4, i);
      GenerateMultAB (out, 6, i);
      
      AlignedGenerateMultAB (out, 1, i, SET);  
      AlignedGenerateMultAB (out, 2, i, SET);
      AlignedGenerateMultAB (out, 3, i, SET);
      AlignedGenerateMultAB (out, 4, i, SET);
      AlignedGenerateMultAB (out, 6, i, SET);
    }


  out << "template <size_t H, OPERATION OP>" << endl
      << "static void MatKernelMultABMask" << endl
      << "(size_t n, SIMD<mask64> mask, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);" << endl;
  out << "template <size_t H, OPERATION OP>" << endl
      << "static void MatKernelMultABMask" << endl
      << "(size_t n, SIMD<mask64> mask, double * pa, size_t da, SIMD<double> * pb, size_t db, double * pc, size_t dc);" << endl;

  GenerateMultABMask (out, 1);  
  GenerateMultABMask (out, 2);
  GenerateMultABMask (out, 3);
  GenerateMultABMask (out, 4);

  
  // Scal AB
  
  out << "template <size_t H, size_t W> static auto MatKernelScalAB" << endl
      << "    (size_t n," << endl
      << "     double * pa, size_t da," << endl
      << "     double * pb, size_t db);" << endl;
  
  GenerateScalAB (out, 3, 4);  
  GenerateScalAB (out, 1, 4);  
  GenerateScalAB (out, 3, 1);  
  GenerateScalAB (out, 1, 1);  


  
  out << "template <size_t H, size_t W>" << endl
      << "static void MyScalTrans" << endl
      << "(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);" << endl;
  
  GenKernel (out, 1, 4);
  GenKernel (out, 2, 4);
  GenKernel (out, 3, 4);
  GenKernel (out, 4, 4);
  GenKernel (out, 5, 4);
  GenKernel (out, 6, 4);
}
