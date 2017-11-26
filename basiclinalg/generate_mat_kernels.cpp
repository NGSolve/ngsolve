// ngscxx generate_mat_kernels.cpp ; a.out

#include <ngstd.hpp>
using namespace ngstd;

enum OP { ADD, SUB, SET };
/*
  C = A * B
  C += A * B
  C -= A * B

  A ... h x n
  B ... n x w*SIMD.Size
 */
void GenerateMultAB (ostream & out, int h, int w, OP op)
{
  
  out << "template <> void MatKernelMultAB<" << h << ", " << w << ">" << endl
      << "    (size_t n," << endl
      << "     double * pa, size_t da," << endl
      << "     double * pb, size_t db," << endl
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


void GenerateMultABMask (ostream & out, int h, OP op)
{
  
  out << "template <> void MatKernelMultABMask<" << h << ">" << endl
      << "    (size_t n, SIMD<mask64> mask," << endl
      << "     double * pa, size_t da," << endl
      << "     double * pb, size_t db," << endl
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
  out << "SIMD<double> b(pb,mask);" << endl;

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

  out << "template <size_t H, size_t W>" << endl
      << "void MatKernelMultAB" << endl
      << "(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);" << endl;

  GenerateMultAB (out, 1, 1, SET);  
  GenerateMultAB (out, 2, 1, SET);
  GenerateMultAB (out, 3, 1, SET);
  GenerateMultAB (out, 4, 1, SET);
  GenerateMultAB (out, 1, 2, SET);  
  GenerateMultAB (out, 2, 2, SET);
  GenerateMultAB (out, 3, 2, SET);
  GenerateMultAB (out, 4, 2, SET);
  GenerateMultAB (out, 1, 3, SET);  
  GenerateMultAB (out, 2, 3, SET);
  GenerateMultAB (out, 3, 3, SET);
  GenerateMultAB (out, 4, 3, SET);



  out << "template <size_t H>" << endl
      << "void MatKernelMultABMask" << endl
      << "(size_t n, SIMD<mask64> mask, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);" << endl;

  GenerateMultABMask (out, 1, SET);  
  GenerateMultABMask (out, 2, SET);
  GenerateMultABMask (out, 3, SET);
  GenerateMultABMask (out, 4, SET);

  
  
  out << "template <size_t H, size_t W>" << endl
      << "void MyScalTrans" << endl
      << "(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);" << endl;
  
  GenKernel (out, 1, 4);
  GenKernel (out, 2, 4);
  GenKernel (out, 3, 4);
  GenKernel (out, 4, 4);
  GenKernel (out, 5, 4);
  GenKernel (out, 6, 4);
}
