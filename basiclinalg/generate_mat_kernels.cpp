// ngscxx generate_mat_kernels.cpp ; a.out

#include "../include/ngs_stdcpp_include.hpp"
#define NGS_DLL_HEADER

using namespace std;

 

#ifndef SIMD_SIZE
  #include "../ngstd/simd_complex.hpp"
  #define SIMD_SIZE ngcore::SIMD<double>::Size()
#endif // SIMD_SIZE

enum OP { ADD, SUB, SET, SETNEG };
enum ORDERING { ColMajor, RowMajor };

string ToString (OP op)
{
  switch (op)
    {
    case SET: return "SET";
    case SETNEG: return "SETNEG";
    case ADD: return "ADD";
    case SUB: return "SUB";
    }
  return "none";  // make the compile happy
}

/*
  callasm: 
  enforce the update fma instruction, important for 12 register version
 */
string FMAOp (OP op, bool callasm = true)
{
  if (callasm)
    return (op == SET || op == ADD) ? "FMAasm" : "FNMAasm";
  else
    return (op == SET || op == ADD) ? "FMAnonasm" : "FNMAnonasm";
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
  out << "constexpr int SW = SIMD<double>::Size();" << endl;

  if (op == SET || op == SETNEG)
    {
      for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++)
          out << "SIMD<double> sum" << i << j << "(0);" << endl;
    }
  else
    {
      out << "double * hpc = pc;" << endl;
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
        out << FMAOp(op, true) << "(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;        
        /*
        if (op == ADD || op == SET)
          out << "FMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;
        else
          out << "FNMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;          
        */
      // out << "sum" << i << j << " -= a" << i << " * b" << j << ";" << endl;
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
  GenerateMultAB (out, h, w, SETNEG, false);
  GenerateMultAB (out, h, w, ADD, false);
  GenerateMultAB (out, h, w, SUB, false);
  GenerateMultAB (out, h, w, SET, true);
  GenerateMultAB (out, h, w, SETNEG, true);
  GenerateMultAB (out, h, w, ADD, true);
  GenerateMultAB (out, h, w, SUB, true);
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
  
  out << "template <> inline void MatKernelAlignedMultAB<" << h << ", " << w << ">" << endl
      << "    (size_t n," << endl
      << "     double * pa, size_t da," << endl
      << "     SIMD<double> * pb, size_t db," << endl
      << "     SIMD<double> * pc, size_t dc)" << endl
      << "{" << endl;

  if (op == SET || op == SETNEG)
    {
      for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++)
          out << "SIMD<double> sum" << i << j << "(0);" << endl;
    }
  else
    {
      out << "SIMD<double> * hpc = pc;" << endl;
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
        out << FMAOp(op, true) << "(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;        
        /*
        if (op == ADD || op == SET)
          out << "FMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;
        else
          out << "FNMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;          
        */
      // out << "sum" << i << j << " -= a" << i << " * b" << j << ";" << endl;          
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
  out << "template <> inline void MatKernelMultABMask<" << h << ", " << ToString(op) << ">" << endl;
    
  out << "    (size_t n, SIMD<mask64> mask," << endl
      << "     double * pa, size_t da," << endl
      << "     " << (aligned_b ? "SIMD<double>" : "double") << " * pb, size_t db," << endl    
      << "     double * pc, size_t dc)" << endl
      << "{" << endl;
  // out << "constexpr int SW = SIMD<double>::Size();" << endl;

  if (op == SET || op == SETNEG)
    {
      for (int i = 0; i < h; i++)
        out << "SIMD<double> sum" << i << "(0);" << endl;
    }
  else
    {
      out << "double * hpc = pc;" << endl;
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
      if (op == SET || op == ADD)
        out << "FMAasm(a"<<i<<",b,sum" << i << ");" << endl;
      else
        out << "FNMAasm(a"<<i<<",b,sum" << i << ");" << endl;        
      // out << "sum" << i << " -= a" << i << "*b;" << endl;
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
  GenerateMultABMask (out, h, SETNEG, false);
  GenerateMultABMask (out, h, ADD, false);
  GenerateMultABMask (out, h, SUB, false);
  GenerateMultABMask (out, h, SET, true);
  GenerateMultABMask (out, h, SETNEG, true);
  GenerateMultABMask (out, h, ADD, true);
  GenerateMultABMask (out, h, SUB, true);
}


/*
  C = A * B^t
  A ... h x n
  B ... w * n
 */
void GenerateScalAB (ostream & out, int h, int w, bool simded)
{
  out << "template <> INLINE auto MatKernelScalAB<" << h << ", " << w << ">" << endl
      << "    (size_t n," << endl
      << "     " << (simded ? "SIMD<double>" : "double") << " * pa, size_t da," << endl
      << "     " << (simded ? "SIMD<double>" : "double") << " * pb, size_t db)" << endl
      << "{" << endl;
  if (!simded)
    out << "constexpr int SW = SIMD<double>::Size();" << endl;

  for (int i = 0; i < h; i++)
    for (int j = 0; j < w; j++)
      out << "SIMD<double> sum" << i << j << "(0);" << endl;

  out << "size_t i = 0;" << endl;
  if (!simded)
    out << "for ( ; i+SW <= n; i+=SW) {" << endl;
  else
    out << "for ( ; i < n; i++) {" << endl;
  
  for (int i = 0; i < h; i++)
    if (simded)
      out << "SIMD<double> a" << i << "(pa[" << i << "*da+i]);" << endl;
    else
      out << "SIMD<double> a" << i << "(pa+" << i << "*da+i);" << endl;
  // for (int i = 0; i < w; i++)
  // out << "SIMD<double> b" << i << "(pb+" << i << "*db+i);" << endl;

  for (int j = 0; j < w; j++)
    {
      if (simded)
        out << "SIMD<double> b" << j << "(pb[" << j << "*db+i]);" << endl;
      else
        out << "SIMD<double> b" << j << "(pb+" << j << "*db+i);" << endl;    
      for (int i = 0; i < h; i++)
        {
          if (h*w < 12)  // with 12 we are on the limit of registers -> fmaasm better
            out << "sum" << i << j << " += a" << i << " * b" << j << ";" << endl;
          else
            out << "FMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;
        }
    }
  out << "}" << endl;

  if (!simded)
    {
      out << "size_t r = n % SW;" << endl;
      out << "if (r) {" << endl;
      out << "SIMD<mask64> mask(r);" << endl;
      for (int i = 0; i < h; i++)
        out << "SIMD<double> a" << i << "(pa+" << i << "*da+i, mask);" << endl;
      
      for (int j = 0; j < w; j++)
        {
          out << "SIMD<double> b" << j << "(pb+" << j << "*db+i, mask);" << endl;
          for (int i = 0; i < h; i++)
            out << "FMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;
        }
      out << "}" << endl;
    }

  if (w == 1 && (h % 4 == 0))
    {
      // out << "return SIMD<double," << h << "> ("; // make_tuple(";
      out << "return Concat(make_tuple(";
      for (int i = 0; i < h; i+=4)
        {
          out << "HSum(sum" << i << "0, sum" << i+1 << "0, sum" << i+2 << "0, sum" << i+3 <<"0)";
          if (i+4 < h) out << ",";
        }
      out << "));"  << endl;
    }

  else

    {
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
    }
  out << "}" << endl;
}


void GenerateScalAB (ostream & out, int h, int w)
{
  GenerateScalAB(out, h, w, false);
  GenerateScalAB(out, h, w, true);
}


/*
  C = A * B^t
  A ... h x n
  B ... w * n

  similar to ScalAB but with nonconstant distances between vectors
*/
void GenerateMultiVecScalAB (ostream & out, int h, int w)
{
  out << "template <> INLINE auto MultiVecScalAB<" << h << ", " << w << ">" << endl
      << "    (size_t n," << endl
      << "     " << "double** ppa," << endl
      << "     " << "double** ppb)" << endl
      << "{" << endl;

    out << "constexpr int SW = SIMD<double>::Size();" << endl;

  for (int i = 0; i < h; i++)
    for (int j = 0; j < w; j++)
      out << "SIMD<double> sum" << i << j << "(0);" << endl;



  for (int i = 0; i < h; i++) {
    out << "double* pa" << i << " = ppa[" << i << "];" << endl;
  }
  for (int i = 0; i < w; i++) {
    out << "double* pb" << i << " = ppb[" << i << "];" << endl;
  }


  out << "size_t i = 0;" << endl;
  out << "for ( ; i+SW <= n; i+=SW) {" << endl;

  for (int i = 0; i < h; i++)
    out << "SIMD<double> a" << i << "(pa" << i << "+i);" << endl;

  for (int j = 0; j < w; j++)
    {
	    out << "SIMD<double> b" << j << "(pb" << j << "+i);" << endl;
      for (int i = 0; i < h; i++)
        {
          if (h*w < 12)  // with 12 we are on the limit of registers -> fmaasm better
            out << "sum" << i << j << " += a" << i << " * b" << j << ";" << endl;
          else
            out << "FMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;
        }
    }
  out << "}" << endl;

  out << "size_t r = n % SW;" << endl;
  out << "if (r) {" << endl;
  out << "SIMD<mask64> mask(r);" << endl;
  for (int i = 0; i < h; i++)
	out << "SIMD<double> a" << i << "(pa" << i << "+i, mask);" << endl;


      for (int j = 0; j < w; j++)
        {
          out << "SIMD<double> b" << j << "(pb" << j << "+i, mask);" << endl;
          for (int i = 0; i < h; i++)
            out << "FMAasm(a"<<i<<",b" << j << ",sum" << i << j << ");" << endl;
        }
      out << "}" << endl;


  if (w == 1 && (h % 4 == 0))
    {
      out << "return make_tuple(";
      for (int i = 0; i < h; i+=4)
        {
          out << "HSum(sum" << i << "0, sum" << i+1 << "0, sum" << i+2 << "0, sum" << i+3 <<"0)";
          if (i+4 < h) out << ",";
        }
      out << ");"  << endl;
    }

  else

    {
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
    }
  out << "}" << endl;
}



/*
  used in GenerateMultiVecScalC if __FMA__ is not defined
*/
void GenerateMultiVecScalC_nofma (ostream & out, int h, int w, bool c)
{
  out << "constexpr int SW = SIMD<double>::Size();" << endl;

  for (int i = 0; i < h; i++)
    for (int j = 0; j < w; j++)
      out << "SIMD<Complex> sum" << i << "_" << j << "(0);" << endl;

  for (int i = 0; i < h; i++) {
    out << "Complex* pa" << i << " = ppa[" << i << "];" << endl;
  }
  for (int i = 0; i < w; i++) {
    out << "Complex* pb" << i << " = ppb[" << i << "];" << endl;
  }


  out << "size_t i = 0;" << endl;
  out << "for ( ; i+SW <= n; i+=SW) {" << endl;

  for (int i = 0; i < h; i++)
  {
    out << "SIMD<Complex> a" << i << ";" << endl;
    out << "a" << i << ".LoadFast(pa" << i << "+i);" << endl;
  }

  for (int j = 0; j < w; j++)
    {
	    out << "SIMD<Complex> b" << j << ";" << endl;
      if (c) {
        out << "b" << j << "=" << "Conj(b" << j << ");" << endl;
      }
      out << "b" << j << ".LoadFast(pb" << j << "+i);" << endl;
      for (int i = 0; i < h; i++)
      {
        if (c) {
          out << "sum" << i << "_" << j << " += a" << i << " * Conj(b" << j << ");" << endl;
        }
        else {
          out << "sum" << i << "_" << j << " += a" << i << " * b" << j << ";" << endl;
        }
      }
    }
  out << "}" << endl;

  out << "int r = n % SW;" << endl;
  out << "if (r) {" << endl;
  for (int i = 0; i < h; i++) {
	   out << "SIMD<Complex> a" << i << ";" << endl;
     out << "a" << i << ".LoadFast(pa" << i << "+i, r);" << endl;
  }

      for (int j = 0; j < w; j++)
        {
          out << "SIMD<Complex> b" << j << ";" << endl;
          if (c) {
            out << "b" << j << "=" << "Conj(b" << j << ");" << endl;
          }
          out << "b" << j << ".LoadFast(pb" << j << "+i, r);" << endl;
          for (int i = 0; i < h; i++) {
            if (c) {
              out << "sum" << i << "_" << j << " += a" << i << " * Conj(b" << j << ");" << endl;
            }
            else {
              out << "sum" << i << "_" << j << " += a" << i << " * b" << j << ";" << endl;
            }
          }
        }
      out << "}" << endl;


  for (int i = 0; i < h; i++) {
    for(int j = 0; j < w; j++) {
      out << "*(pc+" << i << "*dc+" << j << ") = HSum(sum" << i << "_" << j << ");" << endl;
    }
  }

  out << "}" << endl;
}

/*
  used in GenerateMultiVecScalC if __FMA__ is defined
*/
#ifdef __FMA__
void GenerateMultiVecScalC_fma (ostream & out, int h, int w, bool c)
{
  #if defined __AVX512F__
    const string SIMD_TYPE = "__m512d";
    const string SIMD_SHUFFLE = "_mm512_shuffle_pd";
    const string SIMD_MUL = "_mm512_mul_pd";
    const string SIMD_FMAADDSUB = "_mm512_fmaddsub_pd";

    constexpr int shuffle1 = 0b11111111;
    constexpr int shuffle2 = 0b01010101;

    if (c) out << "SIMD<double> conj(_mm512_set_pd(-1,1,-1,1,-1,1,-1,1));" << endl;

  #elif defined __AVX__
    const string SIMD_TYPE = "__m256d";
    const string SIMD_SHUFFLE = "_mm256_shuffle_pd";
    const string SIMD_MUL = "_mm256_mul_pd";
    const string SIMD_FMAADDSUB = "_mm256_fmaddsub_pd";

    constexpr int shuffle1 = 0b1111;
    constexpr int shuffle2 = 0b0101;

    if (c) out << "SIMD<double> conj(1,-1,1,-1);" << endl;

  #elif defined __SSE__
    const string SIMD_TYPE = "__m128d";
    const string SIMD_SHUFFLE = "_mm_shuffle_pd";
    const string SIMD_MUL = "_mm_mul_pd";
    const string SIMD_FMAADDSUB = "_mm_fmaddsub_pd";

    constexpr int shuffle1 = 0b11;
    constexpr int shuffle2 = 0b01;

    if (c) out << "SIMD<double> conj(1,-1);" << endl;
  #endif

  out << "constexpr int SW = SIMD<double>::Size();" << endl;

  for (int i = 0; i < h; i++)
    for (int j = 0; j < w; j++)
      out << "SIMD<double> sum" << i << "_" << j << "(0);" << endl;

  for (int i = 0; i < h; i++)
    out << "double* pa" << i << " = (double*) ppa[" << i << "];" << endl;
  for (int j = 0; j < w; j++)
    out << "double* pb" << j << " = (double*) ppb[" << j << "];" << endl;

  // calculate inner-product
  out << "size_t i = 0;" << endl;
  out << "for (; i+SW <= 2*n; i+=SW) {" << endl; // 2*n since the vectors are complex

  for (int i = 0; i < h; i++) {
    out << "SIMD<double> a" << i << "(pa" << i << " + i);" << endl;
    out << SIMD_TYPE << " a" << i << "Im = " << SIMD_SHUFFLE << "(a" << i << ".Data(), a" << i << ".Data(), " << shuffle1 << ");" << endl;
    out << SIMD_TYPE << " a" << i << "Re = " << SIMD_SHUFFLE << "(a" << i << ".Data(), a" << i << ".Data(), 0);" << endl;
  }
  for (int j = 0; j < w; j++) {
    out << "SIMD<double> b" << j << " (pb" << j << " + i);" << endl;
    if (c) {
      out << "b" << j << " *= conj;" << endl;
    }
    out << SIMD_TYPE << " b" << j << "Swap = " << SIMD_SHUFFLE <<"(b" << j << ".Data(), b" << j << ".Data(), " << shuffle2 << ");" << endl;
    for (int i = 0; i < h; i++) {
      out << SIMD_TYPE << " a" << i << "Im_b" << j << "Swap = "<<SIMD_MUL<<"(a" << i << "Im, b" << j << "Swap);" << endl;
    }
    for (int i = 0; i < h; i++) {
      out << "sum" << i << "_" << j << " += SIMD<double> (" << SIMD_FMAADDSUB << "(a" << i << "Re, b" << j << ".Data(), a" << i << "Im_b" << j << "Swap));" << endl;
    }
  }
  out << "}" << endl;

  // remaining coefficients
  out << "int r = (2*n) % SW;" << endl;
  out << "if (r) {" << endl;
  out << "SIMD<mask64> mask(r);" << endl;

  for (int i = 0; i < h; i++) {
    out << "SIMD<double> a" << i << "(pa" << i << " + i, mask);" << endl;
    out << SIMD_TYPE << " a" << i << "Im = " << SIMD_SHUFFLE << "(a" << i << ".Data(), a" << i << ".Data(), " << shuffle1 << ");" << endl;
    out << SIMD_TYPE << " a" << i << "Re = " << SIMD_SHUFFLE << "(a" << i << ".Data(), a" << i << ".Data(), 0);" << endl;

  }
  for (int j = 0; j < w; j++) {
    out << "SIMD<double> b" << j << " (pb" << j << " + i, mask);" << endl;
    if (c) {
      out << "b" << j << " *= conj;" << endl;
    }
    out << SIMD_TYPE << " b" << j << "Swap = " << SIMD_SHUFFLE << "(b" << j << ".Data(), b" << j << ".Data(), " << shuffle2 << ");" << endl;
    for (int i = 0; i < h; i++) {
      out << SIMD_TYPE << " a" << i << "Im_b" << j << "Swap = " << SIMD_MUL << "(a" << i << "Im, b" << j << "Swap);" << endl;
    }
    for (int i = 0; i < h; i++) {
      out << "sum" << i << "_" << j << " += SIMD<double> (" << SIMD_FMAADDSUB << "(a" << i << "Re, b" << j << ".Data(), a" << i << "Im_b" << j << "Swap));" << endl;
    }
  }
  out << "}" << endl;

  // store results
  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      #if defined __AVX512F__
        out << "pc[" << j << "+" << i << "*dc].real((sum" << i << "_" << j << "[0] + sum" << i << "_" << j << "[2]) + (sum" << i << "_" << j << "[4] + sum" << i << "_" << j << "[6]));" << endl;
        out << "pc[" << j << "+" << i << "*dc].imag((sum" << i << "_" << j << "[1] + sum" << i << "_" << j << "[3]) + (sum" << i << "_" << j << "[5] + sum" << i << "_" << j << "[7]));" << endl;
      #elif defined __AVX__
        out << "pc[" << j << "+" << i << "*dc].real(sum" << i << "_" << j << "[0] + sum" << i << "_" << j << "[2]);" << endl;
        out << "pc[" << j << "+" << i << "*dc].imag(sum" << i << "_" << j << "[1] + sum" << i << "_" << j << "[3]);" << endl;
      #elif defined __SSE__
        out << "pc[" << j << "+" << i << "*dc].real(sum" << i << "_" << j << "[0]);" << endl;
        out << "pc[" << j << "+" << i << "*dc].imag(sum" << i << "_" << j << "[1]);" << endl;
      #endif
    }
  }

out << "}" << endl;
}
#endif

/*
  C = A * B^t
  A ... h x n
  B ... w * n

  bool c for conjugate
*/
void GenerateMultiVecScalC (ostream & out, int h, int w, bool c)
{
  out << "template <> INLINE void MultiVecScalC<" << h << ", " << w << ", " << c << ">" << endl
      << "    (size_t n," << endl
      << "     Complex** ppa," << endl
      << "     Complex** ppb," << endl
      << "     Complex* pc, size_t dc)" << endl
      << "{" << endl;

  // The alternative version using fmaaddsub turned out to be faster
  // If FMA is not available we need
  // SIMD<Complex>, LoadFast and StoreFast
  #ifdef __FMA__
    GenerateMultiVecScalC_fma (out, h, w, c);
  #else
    GenerateMultiVecScalC_nofma (out, h, w, c);
  #endif

}




/*
  A[i] += sum_j c(j,i) * y[j]
  A ... h x n
  B ... w x n
*/
void GenerateMultiScaleAdd (ostream & out, int h, int w)
{
  out << "template <> INLINE void MultiScaleAdd<" << h << ", " << w << ">" << endl
      << "    (size_t n," << endl
      << "     double ** ppa, " << endl
      << "     double ** ppb, " << endl
      << "     double * pc, size_t dc)" << endl
      << "{" << endl;
  out << "constexpr int SW = SIMD<double>::Size();" << endl;

  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      out << "double c" << i << "_" << j << " = pc[" << i << "+" << j << "*dc];" << endl;
    }
  }

  for (int i = 0; i < h; i++) {
    out << "double* pa" << i << " = ppa[" << i << "];" << endl;
  }
  for (int j = 0; j < w; j++) {
    out << "double* pb" << j << " = ppb[" << j << "];" << endl;
  }

  out << "size_t i = 0;" << endl;
  out << "for ( ; i+SW <= n; i+=SW) {" << endl;

  for (int i = 0; i < h; i++)
    out << "SIMD<double> a" << i << "(pa" << i << "+i);" << endl;

  for (int j = 0; j < w; j++)
    {
      out << "SIMD<double> b" << j << "(pb" << j << "+i);" << endl;
      for (int i = 0; i < h; i++)
        {
          out << "a" << i << " += c" << i << "_" << j << " * b" << j << ";" << endl;
        }
    }

  for (int i = 0; i < h; i++)
    out << "a" << i << ".Store(pa" << i << "+i);" << endl;

  out << "}" << endl;

  out << "size_t r = n % SW;" << endl;
  out << "if (r) {" << endl;
  out << "SIMD<mask64> mask(r);" << endl;
  for (int i = 0; i < h; i++)
    out << "SIMD<double> a" << i << "(pa" << i << "+i, mask);" << endl;

  for (int j = 0; j < w; j++)
    {
      out << "SIMD<double> b" << j << "(pb" << j << "+i, mask);" << endl;
      for (int i = 0; i < h; i++)
          out << "a" << i << " += c" << i << "_" << j << " * b" << j << ";" << endl;
    }
  for (int i = 0; i < h; i++)
    out << "a" << i << ".Store(pa" << i << "+i, mask);" << endl;

  out << "}" << endl;

  out << "}" << endl;
}


/*
  used in GenerateMultiScaleAddC if __FMA__ is not defined
*/
void GenerateMultiScaleAddC_nofma (ostream & out, int h, int w)
{
  out << "constexpr int SW = SIMD<double>::Size();" << endl;

  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      out << "Complex c" << i << "_" << j << " = pc[" << i << "+" << j << "*dc];" << endl;
    }
  }

  for (int i = 0; i < h; i++) {
    out << "Complex* pa" << i << " = ppa[" << i << "];" << endl;
  }
  for (int j = 0; j < w; j++) {
    out << "Complex* pb" << j << " = ppb[" << j << "];" << endl;
  }

  out << "size_t i = 0;" << endl;
  out << "for ( ; i+SW <= n; i+=SW) {" << endl;

  for (int i = 0; i < h; i++) {
    out << "SIMD<Complex> a" << i << ";" << endl;
    out << "a" << i << ".LoadFast(pa" << i << "+i);" << endl;
  }

  for (int j = 0; j < w; j++)
    {
      out << "SIMD<Complex> b" << j << ";" << endl;
      out << "b" << j << ".LoadFast(pb" << j << "+i);" << endl;
      for (int i = 0; i < h; i++)
        {
          out << "a" << i << " += c" << i << "_" << j << " * b" << j << ";" << endl;
        }
    }

  for (int i = 0; i < h; i++)
    out << "a" << i << ".StoreFast(pa" << i << "+i);" << endl;

  out << "}" << endl;

  out << "int r = n % SW;" << endl;
  out << "if (r) {" << endl;
  for (int i = 0; i < h; i++) {
    out << "SIMD<Complex> a" << i << ";" << endl;
    out << "a" << i << ".LoadFast(pa" << i << "+i, r);" << endl;
  }

  for (int j = 0; j < w; j++)
    {
      out << "SIMD<Complex> b" << j << ";" << endl;
      out << "b" << j << ".LoadFast(pb" << j << "+i, r);" << endl;
      for (int i = 0; i < h; i++)
          out << "a" << i << " += c" << i << "_" << j << " * b" << j << ";" << endl;
    }
  for (int i = 0; i < h; i++)
    out << "a" << i << ".StoreFast(pa" << i << "+i, r);" << endl;

  out << "}" << endl;

  out << "}" << endl;

}

/*
  used in GenerateMultiScaleAddC if __FMA__ is defined
*/
#ifdef __FMA__
void GenerateMultiScaleAddC_fma (ostream & out, int h, int w)
{
  #if defined __AVX512F__
    const string SIMD_TYPE = "__m512d";
    const string SIMD_SET = "_mm512_set1_pd";
    const string SIMD_SHUFFLE = "_mm512_shuffle_pd";
    const string SIMD_MUL = "_mm512_mul_pd";
    const string SIMD_FMAADDSUB = "_mm512_fmaddsub_pd";

    constexpr int swap_pairs = 0b01010101;

  #elif defined __AVX__
    const string SIMD_TYPE = "__m256d";
    const string SIMD_SET = "_mm256_set1_pd";
    const string SIMD_SHUFFLE = "_mm256_shuffle_pd";
    const string SIMD_MUL = "_mm256_mul_pd";
    const string SIMD_FMAADDSUB = "_mm256_fmaddsub_pd";

    constexpr int swap_pairs = 0b0101;

  #elif defined __SSE__
    const string SIMD_TYPE = "__m128d";
    const string SIMD_SET = "_mm_set1_pd";
    const string SIMD_SHUFFLE = "_mm_shuffle_pd";
    const string SIMD_MUL = "_mm_mul_pd";
    const string SIMD_FMAADDSUB = "_mm_fmaddsub_pd";

    constexpr int swap_pairs = 0b01;
  #endif

  out << "constexpr int SW = SIMD<double>::Size();" << endl;

  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      out << SIMD_TYPE << " c" << i << "_" << j << "Re = " << SIMD_SET << "(pc[" << i << "+" << j << "*dc].real());" << endl;
      out << SIMD_TYPE << " c" << i << "_" << j << "Im = " << SIMD_SET << "(pc[" << i << "+" << j << "*dc].imag());" << endl;
    }
  }

  for (int i = 0; i < h; i++) {
    out << "double* pa" << i << " = (double*) ppa[" << i << "];" << endl;
  }
  for (int j = 0; j < w; j++) {
    out << "double* pb" << j << " = (double*) ppb[" << j << "];" << endl;
  }

  out << "size_t i = 0;" << endl;
  out << "for ( ; i+SW <= 2*n; i+=SW) {" << endl;

  for (int i = 0; i < h; i++)
    out << "SIMD<double> a" << i << "(pa" << i << "+i);" << endl;

  for (int j = 0; j < w; j++)
    {
      out << "SIMD<double> b" << j << "(pb" << j << "+i);" << endl;
      out << SIMD_TYPE << " b" << j << "Swap = " << SIMD_SHUFFLE << "(b" << j << ".Data(), b" << j << ".Data(), " << swap_pairs << ");" << endl;
      for (int i = 0; i < h; i++)
        {
          out << SIMD_TYPE << " c" << i << "_" << j << "Im_b" << j << "Swap = " << SIMD_MUL << "(c" << i << "_" << j << "Im, b" << j << "Swap);" << endl;
          out << "a" << i << " += SIMD<double> (" << SIMD_FMAADDSUB << "(c" << i << "_" << j << "Re, b" << j << ".Data(), c" << i << "_" << j << "Im_b" << j << "Swap));" << endl;
        }
    }

  for (int i = 0; i < h; i++)
    out << "a" << i << ".Store(pa" << i << "+i);" << endl;

  out << "}" << endl;

  out << "size_t r = (2 * n) % SW;" << endl;
  out << "if (r) {" << endl;
  out << "SIMD<mask64> mask(r);" << endl;
  for (int i = 0; i < h; i++)
    out << "SIMD<double> a" << i << "(pa" << i << "+i, mask);" << endl;

  for (int j = 0; j < w; j++)
    {
      out << "SIMD<double> b" << j << "(pb" << j << "+i, mask);" << endl;
      out << SIMD_TYPE << " b" << j << "Swap = " << SIMD_SHUFFLE << "(b" << j << ".Data(), b" << j << ".Data(), " << swap_pairs << ");" << endl;
      for (int i = 0; i < h; i++) {
        out << SIMD_TYPE << " c" << i << "_" << j << "Im_b" << j << "Swap = " << SIMD_MUL << "(c" << i << "_" << j << "Im, b" << j << "Swap);" << endl;
        out << "a" << i << " += SIMD<double> (" << SIMD_FMAADDSUB << "(c" << i << "_" << j << "Re, b" << j << ".Data(), c" << i << "_" << j << "Im_b" << j << "Swap));" << endl;
      }
    }
  for (int i = 0; i < h; i++)
    out << "a" << i << ".Store(pa" << i << "+i, mask);" << endl;

  out << "}" << endl;

  out << "}" << endl;

}
#endif

/*
  A[i] += sum_j c(j,i) * y[j]
  A ... h x n
  B ... w x n
*/
void GenerateMultiScaleAddC (ostream & out, int h, int w)
{
  out << "template <> INLINE void MultiScaleAddC<" << h << ", " << w << ">" << endl
      << "    (size_t n," << endl
      << "     Complex ** ppa, " << endl
      << "     Complex ** ppb, " << endl
      << "     Complex * pc, size_t dc)" << endl
      << "{" << endl;

  // The alternative version using fmaaddsub turned out to be faster.
  // If FMA is not available we need
  // SIMD<Complex>, LoadFast and StoreFast
  #ifdef __FMA__
    GenerateMultiScaleAddC_fma(out, h, w);
  #else
    GenerateMultiScaleAddC_nofma(out, h, w);
  #endif
}





void GenKernel (ofstream & out, int h, int w)
{
  out << "template <> inline void MyScalTrans<" << h << ", " << w << ">" << endl
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







/*
  C = A * B
  C += A * B
  C -= A * B

  A ... h x w
  B ... w x n
 */
void GenerateDaxpy (ostream & out, int h, int w, OP op, bool aligned_b)
{
  out << "template <> INLINE void MatKernelDaxpy<" << h << ", " << w << ", " << ToString(op) << ">" << endl;
  out << "    (size_t n," << endl
      << "     double * pa, size_t da," << endl
      << "     " << (aligned_b ? "SIMD<double>" : "double") << " * pb, size_t db," << endl
      << "     " << (aligned_b ? "SIMD<double>" : "double") << " * pc, size_t dc)" << endl
      << "{" << endl;
  out << "constexpr int SW = SIMD<double>::Size();" << endl;

  for (int i = 0; i < h; i++)
    for (int j = 0; j < w; j++)
      out << "SIMD<double> a" << i << j << "(pa[" << i << "*da+"<< j << "]);" << endl;


  for (int i = 0; i < h; i++)
    out << "double * pc" << i << " = pc+" << i << "*dc;" << endl;
  for (int i = 0; i < w; i++)
    out << "double * pb" << i << " = pb+" << i << "*db;" << endl;

  out << "size_t i = 0;" << endl;
  out << "for ( ; i+SW <= n; i+=SW) {" << endl;
  
  
  if (op == SET || op == SETNEG)
    {
      for (int i = 0; i < h; i++)
        out << "SIMD<double> c" << i << "(0);" << endl;
    }
  else
    {
      for (int i = 0; i < h; i++)
        out << "SIMD<double> c" << i << "(pc" << i << "+i);" << endl;
    }
  
  /*
    if (aligned_b)
    for (int i = 0; i < w; i++)
    out << "SIMD<double> b" << i << " = pb[" << i << "];" << endl;
    else
    for (int i = 0; i < w; i++)
    out << "SIMD<double> b" << i << "(pb+" << i << "*SW);" << endl;
  */
  
  for (int j = 0; j < w; j++)
    {
      out << "SIMD<double> b" << j << "(pb" << j << "+i);" << endl;
      for (int i = 0; i < h; i++)
        out << FMAOp(op, false) << "(a" << i << j << ", b" << j << ", c" << i << ");\n";
        /*
        if (op == ADD || op == SET)
          out << "c" << i << " += a" << i  << j << " * b" << j << ";" << endl;
        else
          out << "c" << i << " -= a" << i  << j << " * b" << j << ";" << endl;
        */
    }

  for (int i = 0; i < h; i++)
    out << "c" << i << ".Store(pc" << i << "+i);" << endl;
  
  out << "}" << endl;



  out << "SIMD<mask64> mask(n%SW);" << endl;
  if (op == SET || op == SETNEG)
    {
      for (int i = 0; i < h; i++)
        out << "SIMD<double> c" << i << "(0);" << endl;
    }
  else
    {
      for (int i = 0; i < h; i++)
        out << "SIMD<double> c" << i << "(pc" << i << "+i, mask);" << endl;
    }
  
  /*
    if (aligned_b)
    for (int i = 0; i < w; i++)
    out << "SIMD<double> b" << i << " = pb[" << i << "];" << endl;
    else
    for (int i = 0; i < w; i++)
    out << "SIMD<double> b" << i << "(pb+" << i << "*SW);" << endl;
  */
  
  for (int j = 0; j < w; j++)
    {
      out << "SIMD<double> b" << j << "(pb" << j << "+i, mask);" << endl;
      for (int i = 0; i < h; i++)
        if (op == ADD || op == SET)
          out << "c" << i << " += a" << i  << j << " * b" << j << ";" << endl;
        else
          out << "c" << i << " -= a" << i  << j << " * b" << j << ";" << endl;
    }

  for (int i = 0; i < h; i++)
    out << "c" << i << ".Store(pc" << i << "+i, mask);" << endl;

  

  
  out << "}" << endl;
}

void GenerateDaxpy (ostream & out, int h, int w)
{
  GenerateDaxpy (out, h, w, SET, false);
  GenerateDaxpy (out, h, w, ADD, false);
  GenerateDaxpy (out, h, w, SUB, false);
  /*
  GenerateDaxpy (out, h, w, SET, true);
  GenerateDaxpy (out, h, w, ADD, true);
  GenerateDaxpy (out, h, w, SUB, true);
  */
}







void GenerateShortSum (ostream & out, int wa, OP op)
{
  out << "template <> INLINE void MatKernelShortSum<" << wa << ", " << ToString(op) << ">" << endl;
  out << "    (size_t ha, size_t wb," << endl
      << "     double * pa, size_t da," << endl
      << "     double * pb, size_t db," << endl
      << "     double * pc, size_t dc)" << endl
      << "{" << endl;
  out << "constexpr int SW = SIMD<double>::Size();\n" 
      << "for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)\n"
      << "{\n";
  if (wa > 0)
    out << "double * pb2 = pb;\n";
  for (int k = 0; k < wa; k++)
    out << "SIMD<double> b" << k << "(pb2); pb2 += db;\n";
  out << "double * pa2 = pa;\n"
      << "double * pc2 = pc;\n"
      << "__assume(ha>0);\n";
  out << "#pragma unroll 1\n";
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)\n"
      << "{\n";
  if (op == SET || op == SETNEG)
    out << "SIMD<double> sum = 0.0;\n";
  else
    out << "SIMD<double> sum(pc2);\n";
    
  for (int k = 0; k < wa; k++)
    if (op == SET || op == ADD)
      // out << "sum += SIMD<double>(pa2[" << k << "]) * b"<< k << ";\n";
      out << "FMAasm (b" << k << ",SIMD<double>(pa2[" << k << "]), sum" <<");\n";      
    else
      // out << "sum -= SIMD<double>(pa2[" << k << "]) * b"<< k << ";\n";
      out << "FNMAasm (b" << k << ",SIMD<double>(pa2[" << k << "]), sum" <<");\n";        
  out << "sum.Store(pc2);\n"
      << "} }\n";

  out << "size_t rest = wb % SW; \n"
      << "if (rest == 0) return; \n"
      << "SIMD<mask64> mask(rest); \n";

  if (wa > 0)
    out << "double * pb2 = pb;\n";
  for (int k = 0; k < wa; k++)
    out << "SIMD<double> b" << k << "(pb2, mask); pb2 += db;\n";
  out << "double * pa2 = pa;\n"
      << "double * pc2 = pc;\n"
      << "__assume(ha>0);\n";

  out << "#pragma unroll 1\n";  
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)\n"
      << "{\n";
  if (op == SET || op == SETNEG)
    out << "SIMD<double> sum = 0.0;\n";
  else
    out << "SIMD<double> sum(pc2, mask);\n";
  for (int k = 0; k < wa; k++)
    if (op == SET || op == ADD)
      out << "sum += SIMD<double>(pa2[" << k << "]) * b"<< k << ";\n";
    else
      out << "sum -= SIMD<double>(pa2[" << k << "]) * b"<< k << ";\n";      
  out << "sum.Store(pc2, mask);\n"
      << "} }\n";





  // unroll B width 2

  out << "template <> INLINE void MatKernelShortSum2<" << wa << ", " << ToString(op) << ">" << endl;
  out << "    (size_t ha, size_t wb," << endl
      << "     double * pa, size_t da," << endl
      << "     double * pb, size_t db," << endl
      << "     double * pc, size_t dc)" << endl
      << "{" << endl;
  out << "constexpr int SW = SIMD<double>::Size();\n" 
      << "for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)\n"
      << "{\n";
  if (wa > 0)
    out << "double * pb2 = pb;\n";
  for (int k = 0; k < wa; k++)
    out << "SIMD<double> b" << k << "0(pb2);\n"
        << "SIMD<double> b" << k << "1(pb2+SW); pb2 += db;\n";
  out << "double * pa2 = pa;\n"
      << "double * pc2 = pc;\n"
      << "__assume(ha>0);\n";
  
  // out << "#pragma unroll 2\n";  
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)\n"
      << "{\n";
  if (op == SET || op == SETNEG)
    {
      out << "SIMD<double> sum0 = 0.0;\n"
          << "SIMD<double> sum1 = 0.0;\n";
    }
  else
    {
      out << "SIMD<double> sum0(pc2);\n";      
      out << "SIMD<double> sum1(pc2+SW);\n";      
    }
  for (int k = 0; k < wa; k++)
      out << FMAOp(op) << "(SIMD<double>(pa2[" << k << "])," << "b" << k << "0, sum0);\n"
          << FMAOp(op) << "(SIMD<double>(pa2[" << k << "])," << "b" << k << "1, sum1);\n";

      /*
      out << "sum0 += SIMD<double>(pa2[" << k << "]) * b"<< k << "0;\n"
          << "sum1 += SIMD<double>(pa2[" << k << "]) * b"<< k << "1;\n";
    else
      out << "sum0 -= SIMD<double>(pa2[" << k << "]) * b"<< k << "0;\n"
          << "sum1 -= SIMD<double>(pa2[" << k << "]) * b"<< k << "1;\n";
      */
  out << "sum0.Store(pc2);\n"
      << "sum1.Store(pc2+SW);\n"
      << "} }\n";
  
  out << "size_t rest = wb % (2*SW); \n"
      << "if (rest == 0) return; \n";

  
  for (int r : { 8, 4, 2, 1})
    {
      if (r > SIMD_SIZE) continue;
      
      out << "if (rest & " << r << ") {  \n";
      if (wa > 0)
        out << "double * pb2 = pb;\n";
      for (int k = 0; k < wa; k++)
        out << "SIMD<double,"<<r<<"> b" << k << "(pb2); pb2 += db;\n";
      out << "double * pa2 = pa;\n"
          << "double * pc2 = pc;\n"
          << "__assume(ha>0);\n";
  
      out << "#pragma unroll 1\n";  
      out << "for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)\n"
          << "{\n";
      if (op == SET || op == SETNEG)
        out << "SIMD<double,"<<r<<"> sum = 0.0;\n";
      else
        out << "SIMD<double,"<<r<<"> sum(pc2);\n";
      
      for (int k = 0; k < wa; k++)
        if (op == SET || op == ADD)    
          out << "sum += SIMD<double,"<<r<<">(pa2[" << k << "]) * b"<< k << ";\n";
        else
          out << "sum -= SIMD<double,"<<r<<">(pa2[" << k << "]) * b"<< k << ";\n";
      
      out << "sum.Store(pc2);\n"
          << "}\n";

      out << "pc += " << r << ";\n";
      out << "pb += " << r << ";\n";      
      out << "}\n";
    }
  out << "return; \n";

  


  
  out << "if (rest >= SW) \n"
      << "{\n"
      << "if (rest > SW)\n"
      << "{\n";

  out << "SIMD<mask64> mask(rest-SW); \n";
  if (wa > 0)
    out << "double * pb2 = pb;\n";
  for (int k = 0; k < wa; k++)
    out << "SIMD<double> b" << k << "0(pb2);\n"
        << "SIMD<double> b" << k << "1(pb2+SW,mask); pb2 += db;\n";
  out << "double * pa2 = pa;\n"
      << "double * pc2 = pc;\n"
      << "__assume(ha>0);\n";

  out << "#pragma unroll 1\n";    
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)\n"
      << "{\n";

  if (op == SET || op == SETNEG)
    {
      out << "SIMD<double> sum0 = 0.0;\n"
          << "SIMD<double> sum1 = 0.0;\n";
    }
  else
    {
      out << "SIMD<double> sum0(pc2);\n";      
      out << "SIMD<double> sum1(pc2+SW,mask);\n";      
    }
  for (int k = 0; k < wa; k++)
    if (op == SET || op == ADD)
      out << "sum0 += SIMD<double>(pa2[" << k << "]) * b"<< k << "0;\n"
          << "sum1 += SIMD<double>(pa2[" << k << "]) * b"<< k << "1;\n";
    else
      out << "sum0 -= SIMD<double>(pa2[" << k << "]) * b"<< k << "0;\n"
          << "sum1 -= SIMD<double>(pa2[" << k << "]) * b"<< k << "1;\n";
      
  out << "sum0.Store(pc2);\n"
      << "sum1.Store(pc2+SW,mask);\n"
      << "}\n";
    
  out << "return;\n"
      << "}\n";
    
    // rest == SW
  if (wa > 0)
    out << "double * pb2 = pb;\n";
  for (int k = 0; k < wa; k++)
    out << "SIMD<double> b" << k << "(pb2); pb2 += db;\n";
  out << "double * pa2 = pa;\n"
      << "double * pc2 = pc;\n"
      << "__assume(ha>0);\n";
  
  out << "#pragma unroll 1\n";  
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)\n"
      << "{\n";
  if (op == SET || op == SETNEG)
    out << "SIMD<double> sum = 0.0;\n";
  else
    out << "SIMD<double> sum(pc2);\n";
    
  for (int k = 0; k < wa; k++)
    if (op == SET || op == ADD)    
      out << "sum += SIMD<double>(pa2[" << k << "]) * b"<< k << ";\n";
    else
      out << "sum -= SIMD<double>(pa2[" << k << "]) * b"<< k << ";\n";
      
  out << "sum.Store(pc2);\n"
      << "}\n";
  
  out << "return;\n"
      << "}\n";
  
  
  // rest < SW
  out << "SIMD<mask64> mask(rest); \n";
  if (wa > 0)
    out << "double * pb2 = pb;\n";
  for (int k = 0; k < wa; k++)
    out << "SIMD<double> b" << k << "(pb2, mask); pb2 += db;\n";
  out << "double * pa2 = pa;\n"
      << "double * pc2 = pc;\n"
      << "__assume(ha>0);\n";
  
  out << "#pragma unroll 1\n";  
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)\n"
      << "{\n";
  if (op == SET || op == SETNEG)
    out << "SIMD<double> sum = 0.0;\n";
  else
    out << "SIMD<double> sum(pc2, mask);\n";
  
  for (int k = 0; k < wa; k++)
    if (op == SET || op == ADD)    
      out << "sum += SIMD<double>(pa2[" << k << "]) * b"<< k << ";\n";
    else
      out << "sum -= SIMD<double>(pa2[" << k << "]) * b"<< k << ";\n";
  out << "sum.Store(pc2, mask);\n"
      << "} }\n";




}





void GenerateAtB_SmallWA (ostream & out, int wa, OP op)
{
  out << "template <> INLINE void MatKernelAtB_SmallWA<" << wa << ", " << ToString(op) << ">" << endl;
  out << "    (size_t ha, size_t wb," << endl
      << "     double * pa, size_t da," << endl
      << "     double * pb, size_t db," << endl
      << "     double * pc, size_t dc)" << endl
      << "{" << endl;
  out << "constexpr int SW = SIMD<double>::Size();\n" 
      << "for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)\n"
      << "{\n"
      << "double * pb2 = pb;\n";
        
  out << "double * pa2 = pa;\n"
      << "[[maybe_unused]] double * pc2 = pc;\n"
      << "__assume(ha>0);\n";

  for (int k = 0; k < wa; k++)
    if (op == SET || op == SETNEG)
      out << "SIMD<double> sum" << k << "(0.0);\n";
    else
      out << "SIMD<double> sum" << k << "(pc2); pc2 += dc;\n";
  
  out << "#pragma unroll 1\n";
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)\n"
      << "{\n";
  out << "SIMD<double> bjk(pb2);\n";
  for (int k = 0; k < wa; k++)
    out << FMAOp(op) << "(bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";
    // else
      // out << "FNMAasm (bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";      
  // out << "sum" << k << " -= bjk * SIMD<double>(pa2[" << k << "]);\n";
  out << "}\n";
  out << "pc2 = pc;\n";
  for (int k = 0; k < wa; k++)
    out << "sum" << k << ".Store(pc2); pc2 += dc;\n";
  out << "}\n";

  out << "size_t rest = wb % SW; \n"
      << "if (rest == 0) return; \n"
      << "SIMD<mask64> mask(rest); \n";
  
  out << "double * pb2 = pb;\n";
  // for (int k = 0; k < wa; k++)
  // out << "SIMD<double> sum" << k << "(0.0);\n";    
  out << "double * pa2 = pa;\n"
      << "[[maybe_unused]] double * pc2 = pc;\n"
      << "__assume(ha>0);\n";

  for (int k = 0; k < wa; k++)
    if (op == SET || op == SETNEG)
      out << "SIMD<double> sum" << k << "(0.0);\n";
    else
      out << "SIMD<double> sum" << k << "(pc2,mask); pc2 += dc;\n";
  
  out << "#pragma unroll 1\n";  
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)\n"
      << "{\n"
      << "SIMD<double> bjk(pb2, mask);\n";    
  for (int k = 0; k < wa; k++)
    out << FMAOp(op) << "(bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";    
    // out << "FMAasm (bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";
/*
    if (op == ADD || op == SET)    
      out << "FMAasm (bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";
    else
      out << "FNMAasm (bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";      
  */
  out << "}\n";
  out << "pc2 = pc;\n";
  for (int k = 0; k < wa; k++)
    out << "sum" << k << ".Store(pc2, mask); pc2 += dc;\n";
  
  out << "}\n";
}






void GenerateAtB_SmallWA2 (ostream & out, int wa, OP op)
{
  out << "template <> INLINE void MatKernelAtB_SmallWA2<" << wa << ", " << ToString(op) << ">" << endl;
  out << "    (size_t ha, size_t wb," << endl
      << "     double * pa, size_t da," << endl
      << "     double * pb, size_t db," << endl
      << "     double * pc, size_t dc)" << endl
      << "{" << endl;
  out << "constexpr int SW = SIMD<double>::Size();\n";
  out << "size_t i = 0;\n";
  
  out << "for ( ; i+3*SW <= wb; i += 3*SW, pb += 3*SW, pc += 3*SW)\n"
      << "{\n"
      << "double * pb2 = pb;\n";
        
  out << "double * pa2 = pa;\n"
      << "[[maybe_unused]] double * pc2 = pc;\n"
      << "__assume(ha>0);\n";

  for (int k = 0; k < wa; k++)
    if (op == SET || op == SETNEG)
      {
        out << "SIMD<double> sum" << k << "0(0.0);\n";
        out << "SIMD<double> sum" << k << "1(0.0);\n";
        out << "SIMD<double> sum" << k << "2(0.0);\n";
      }
    else
      {
        out << "SIMD<double> sum" << k << "0(pc2);\n";
        out << "SIMD<double> sum" << k << "1(pc2+SW);\n";
        out << "SIMD<double> sum" << k << "2(pc2+2*SW); pc2 += dc;\n";
      }
  
  out << "#pragma unroll 1\n";
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)\n"
      << "{\n";
  out << "SIMD<double> bjk0(pb2);\n";
  out << "SIMD<double> bjk1(pb2+SW);\n";
  out << "SIMD<double> bjk2(pb2+2*SW);\n";
  for (int k = 0; k < wa; k++)
    {
      out << FMAOp(op, true) << "(bjk0,SIMD<double>(pa2[" << k << "]), sum" << k <<"0);\n";
      out << FMAOp(op, true) << "(bjk1,SIMD<double>(pa2[" << k << "]), sum" << k <<"1);\n";
      out << FMAOp(op, true) << "(bjk2,SIMD<double>(pa2[" << k << "]), sum" << k <<"2);\n";
    }
    /*
    if (op == ADD || op == SET)
      {
        out << "FMAasm (bjk0,SIMD<double>(pa2[" << k << "]), sum" << k <<"0);\n";
        out << "FMAasm (bjk1,SIMD<double>(pa2[" << k << "]), sum" << k <<"1);\n";
        out << "FMAasm (bjk2,SIMD<double>(pa2[" << k << "]), sum" << k <<"2);\n";
      }
    else
      {
        out << "FNMAasm (bjk0,SIMD<double>(pa2[" << k << "]), sum" << k <<"0);\n";
        out << "FNMAasm (bjk1,SIMD<double>(pa2[" << k << "]), sum" << k <<"1);\n";
        out << "FNMAasm (bjk2,SIMD<double>(pa2[" << k << "]), sum" << k <<"2);\n";
      }
    */
  out << "}\n";
  out << "pc2 = pc;\n";
  for (int k = 0; k < wa; k++)
    {
      out << "sum" << k << "0.Store(pc2);\n";
      out << "sum" << k << "1.Store(pc2+SW);\n";
      out << "sum" << k << "2.Store(pc2+2*SW); pc2 += dc;\n";
    }      
  out << "}\n";


  out << "for ( ; i+SW <= wb; i += SW, pb += SW, pc += SW)\n"
      << "{\n"
      << "double * pb2 = pb;\n";
        
  out << "double * pa2 = pa;\n"
      << "[[maybe_unused]] double * pc2 = pc;\n"
      << "__assume(ha>0);\n";

  for (int k = 0; k < wa; k++)
    if (op == SET || op == SETNEG)
      out << "SIMD<double> sum" << k << "(0.0);\n";
    else
      out << "SIMD<double> sum" << k << "(pc2); pc2 += dc;\n";
  
  out << "#pragma unroll 1\n";
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)\n"
      << "{\n";
  out << "SIMD<double> bjk(pb2);\n";
  for (int k = 0; k < wa; k++)
    if (op == ADD || op == SET)    
      out << "FMAasm (bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";
    else
      out << "FNMAasm (bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";
  out << "}\n";
  out << "pc2 = pc;\n";
  for (int k = 0; k < wa; k++)
    out << "sum" << k << ".Store(pc2); pc2 += dc;\n";
  out << "}\n";







  
  out << "size_t rest = wb % SW; \n"
      << "if (rest == 0) return; \n"
      << "SIMD<mask64> mask(rest); \n";
  
  out << "double * pb2 = pb;\n";
  // for (int k = 0; k < wa; k++)
  // out << "SIMD<double> sum" << k << "(0.0);\n";    
  out << "double * pa2 = pa;\n"
      << "[[maybe_unused]] double * pc2 = pc;\n"
      << "__assume(ha>0);\n";

  for (int k = 0; k < wa; k++)
    if (op == SET || op == SETNEG)
      out << "SIMD<double> sum" << k << "(0.0);\n";
    else
      out << "SIMD<double> sum" << k << "(pc2,mask); pc2 += dc;\n";
  
  out << "#pragma unroll 1\n";  
  out << "for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)\n"
      << "{\n"
      << "SIMD<double> bjk(pb2, mask);\n";    
  for (int k = 0; k < wa; k++)
    // out << "FMAasm (bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";
    if (op == ADD || op == SET)    
      out << "FMAasm (bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";
    else
      out << "FNMAasm (bjk,SIMD<double>(pa2[" << k << "]), sum" << k <<");\n";      
  
  out << "}\n";
  out << "pc2 = pc;\n";
  for (int k = 0; k < wa; k++)
    out << "sum" << k << ".Store(pc2, mask); pc2 += dc;\n";
  
  out << "}\n";
}









void  GenerateMatVec (ostream & out, int wa, OP op)
{
  out << "template <> INLINE void KernelMatVec<" << wa << ", " << ToString(op) << ">" << endl
      << "(size_t ha, double * pa, size_t da, double * x, double * y) {" << endl;

  int SW = SIMD_SIZE;  // generate optimal code for my host
  // out << "constexpr int SW = SIMD<double>::Size();" << endl;
  int i = 0;
  for ( ; SW*(i+1) <= wa; i++)
    out << "SIMD<double," << SW << "> x" << i << "(x+" << i*SW << ");" << endl;
  
  if (SW == 4 && (wa % SW == 1))
    {
      out << "double x" << i << " = x[" << i*SW << "];" << endl;      
    }
  else if (SW == 4 && (wa % SW == 2))
    {
      out << "SIMD<double,2> x" << i << "(x+" << i*SW << ");" << endl;      
    }
  else if (wa % SW)  // do the mask load :-(
    {
      out << "SIMD<mask64," << SW << "> mask(size_t(" << wa % SW << "));" << endl;
      out << "SIMD<double," << SW << "> x" << i << "(x+" << i*SW << ", mask);" << endl;
    }
  out << "size_t i = 0;" << endl;
  out << "for ( ; i+4 <= ha; i+=4, pa += 4*da) {" << endl;
  out << "SIMD<double," << SW << "> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);" << endl;
  i = 0;
  for ( ; SW*(i+1) <= wa; i++)
    {
      out << "sum0 += SIMD<double," << SW << ">(pa+" << i*SW << ") * x" << i << ";" << endl;
      out << "sum1 += SIMD<double," << SW << ">(pa+da+" << i*SW << ") * x" << i << ";" << endl;
      out << "sum2 += SIMD<double," << SW << ">(pa+2*da+" << i*SW << ") * x" << i << ";" << endl;
      out << "sum3 += SIMD<double," << SW << ">(pa+3*da+" << i*SW << ") * x" << i << ";" << endl;
    }

  if (SW == 4 && (wa % SW == 1))
    {
      /*
      for (int k = 0; k < 4; k++)
        out << "sum"<<k<< " += SIMD<double,4> (pa[" << k << "*da+" << i*SW << "] * x" << i << ", 0,0,0);" << endl;
      */
      ;
    }
  else if (SW == 4 && (wa % SW == 2))
    {
      out << "SIMD<double,2> zero(0.0);" << endl;
      out << "sum0 += SIMD<double,4> (SIMD<double,2>(pa+" << i*SW << ") * x" << i << ", zero);" << endl;      
      out << "sum1 += SIMD<double,4> (SIMD<double,2>(pa+da+" << i*SW << ") * x" << i << ", zero);" << endl;      
      out << "sum2 += SIMD<double,4> (SIMD<double,2>(pa+2*da+" << i*SW << ") * x" << i << ", zero);" << endl;      
      out << "sum3 += SIMD<double,4> (SIMD<double,2>(pa+3*da+" << i*SW << ") * x" << i << ", zero);" << endl;      
    }
  else if (wa % SW)
    {
      out << "sum0 += SIMD<double," << SW << ">(pa+" << i*SW << ", mask) * x" << i << ";" << endl;
      out << "sum1 += SIMD<double," << SW << ">(pa+da+" << i*SW << ", mask) * x" << i << ";" << endl;
      out << "sum2 += SIMD<double," << SW << ">(pa+2*da+" << i*SW << ", mask) * x" << i << ";" << endl;
      out << "sum3 += SIMD<double," << SW << ">(pa+3*da+" << i*SW << ", mask) * x" << i << ";" << endl;
    }
  out << "SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);" << endl;

  if (SW == 4 && (wa % SW == 1))
    {
      out << "vsum += x" << i << "*SIMD<double,4> ("
          << "pa[0*da+" << i*SW << "], "
          << "pa[1*da+" << i*SW << "], "
          << "pa[2*da+" << i*SW << "], "
          << "pa[3*da+" << i*SW << "]);" << endl;
    }

  if (op == ADD)
    out << "vsum += SIMD<double,4> (y+i);" << endl;
  out << "vsum.Store(y+i);" << endl;
  out << "}" << endl;


  out << "if (ha & 2) {" << endl;
  out << "SIMD<double," << SW << "> sum0(0.0), sum1(0.0);" << endl;
  i = 0;
  for ( ; SW*(i+1) <= wa; i++)
    {
      out << "sum0 += SIMD<double," << SW << ">(pa+" << i*SW << ") * x" << i << ";" << endl;
      out << "sum1 += SIMD<double," << SW << ">(pa+da+" << i*SW << ") * x" << i << ";" << endl;
    }
  
  if (SW == 4 && (wa % SW == 1))
    {
      for (int k = 0; k < 2; k++)
        out << "sum"<<k<< " += SIMD<double,4> (pa[" << k << "*da+" << i*SW << "] * x" << i << ", 0,0,0);" << endl;      
    }
  else if (SW == 4 && (wa % SW == 2))
    {
      out << "SIMD<double,2> zero(0.0);" << endl;
      out << "sum0 += SIMD<double,4> (SIMD<double,2>(pa+" << i*SW << ") * x" << i << ", zero);" << endl;      
      out << "sum1 += SIMD<double,4> (SIMD<double,2>(pa+da+" << i*SW << ") * x" << i << ", zero);" << endl;      
    }
  else if (wa % SW)
    {
      out << "sum0 += SIMD<double," << SW << ">(pa+" << i*SW << ", mask) * x" << i << ";" << endl;
      out << "sum1 += SIMD<double," << SW << ">(pa+da+" << i*SW << ", mask) * x" << i << ";" << endl;
    }
  out << "SIMD<double,2> vsum = HSum(sum0,sum1);" << endl;
  if (op == ADD)
    out << "vsum += SIMD<double,2> (y+i);" << endl;  
  out << "vsum.Store(y+i);" << endl;
  out << "i += 2; pa += 2*da;" << endl;
  out << "}" << endl;
  
  
  out << "if (ha & 1) {" << endl;
  out << "SIMD<double," << SW << "> sum(0.0);" << endl;
  i = 0;
  for ( ; SW*(i+1) <= wa; i++)
    out << "sum += SIMD<double," << SW << ">(pa+" << i*SW << ") * x" << i << ";" << endl;

  
  if (SW == 4 && (wa % SW == 1))
    {
      out << "sum += SIMD<double,4> (pa[" << i*SW << "] * x" << i << ", 0,0,0);" << endl;      
    }
  else if (SW == 4 && (wa % SW == 2))
    {
      out << "SIMD<double,2> zero(0.0);" << endl;
      out << "sum += SIMD<double,4> (SIMD<double,2>(pa+" << i*SW << ") * x" << i << ", zero);" << endl;      
    }
  else if (wa % SW)
    out << "sum += SIMD<double," << SW << ">(pa+" << i*SW << ", mask) * x" << i << ";" << endl;
  if (op == ADD)
    out << "y[i] += HSum(sum);" << endl;
  else
    out << "y[i] = HSum(sum);" << endl;

  out << "} }" << endl;
}



void GenerateAddMatVec (ostream & out, int wa)
{
  out << "template <> INLINE void KernelAddMatVec<" << wa << ">" << endl
      << "(double s, size_t ha, double * pa, size_t da, double * x, double * y) {" << endl;

  int SW = SIMD_SIZE;  // generate optimal code for my host
  int i = 0;
  for ( ; SW*(i+1) <= wa; i++)
    out << "SIMD<double," << SW << "> x" << i << "(x+" << i*SW << ");" << endl;
  
  if (SW == 4 && (wa % SW == 1))
    {
      out << "double x" << i << " = x[" << i*SW << "];" << endl;      
    }
  else if (SW == 4 && (wa % SW == 2))
    {
      out << "SIMD<double,2> x" << i << "(x+" << i*SW << ");" << endl;      
    }
  else if (wa % SW)  // do the mask load :-(
    {
      out << "SIMD<mask64," << SW << "> mask(size_t(" << wa % SW << "));" << endl;
      out << "SIMD<double," << SW << "> x" << i << "(x+" << i*SW << ", mask);" << endl;
    }
  out << "size_t i = 0;" << endl;
  out << "for ( ; i+4 <= ha; i+=4, pa += 4*da) {" << endl;
  out << "SIMD<double," << SW << "> sum0(0.0), sum1(0.0), sum2(0.0), sum3(0.0);" << endl;
  i = 0;
  for ( ; SW*(i+1) <= wa; i++)
    {
      out << "sum0 += SIMD<double," << SW << ">(pa+" << i*SW << ") * x" << i << ";" << endl;
      out << "sum1 += SIMD<double," << SW << ">(pa+da+" << i*SW << ") * x" << i << ";" << endl;
      out << "sum2 += SIMD<double," << SW << ">(pa+2*da+" << i*SW << ") * x" << i << ";" << endl;
      out << "sum3 += SIMD<double," << SW << ">(pa+3*da+" << i*SW << ") * x" << i << ";" << endl;
    }

  if (SW == 4 && (wa % SW == 1))
    {
      /*
      for (int k = 0; k < 4; k++)
        out << "sum"<<k<< " += SIMD<double,4> (pa[" << k << "*da+" << i*SW << "] * x" << i << ", 0,0,0);" << endl;
      */
      ;
    }
  else if (SW == 4 && (wa % SW == 2))
    {
      out << "SIMD<double,2> zero(0.0);" << endl;
      out << "sum0 += SIMD<double,4> (SIMD<double,2>(pa+" << i*SW << ") * x" << i << ", zero);" << endl;      
      out << "sum1 += SIMD<double,4> (SIMD<double,2>(pa+da+" << i*SW << ") * x" << i << ", zero);" << endl;      
      out << "sum2 += SIMD<double,4> (SIMD<double,2>(pa+2*da+" << i*SW << ") * x" << i << ", zero);" << endl;      
      out << "sum3 += SIMD<double,4> (SIMD<double,2>(pa+3*da+" << i*SW << ") * x" << i << ", zero);" << endl;      
    }
  else if (wa % SW)
    {
      out << "sum0 += SIMD<double," << SW << ">(pa+" << i*SW << ", mask) * x" << i << ";" << endl;
      out << "sum1 += SIMD<double," << SW << ">(pa+da+" << i*SW << ", mask) * x" << i << ";" << endl;
      out << "sum2 += SIMD<double," << SW << ">(pa+2*da+" << i*SW << ", mask) * x" << i << ";" << endl;
      out << "sum3 += SIMD<double," << SW << ">(pa+3*da+" << i*SW << ", mask) * x" << i << ";" << endl;
    }
  out << "SIMD<double,4> vsum = HSum(sum0,sum1,sum2,sum3);" << endl;

  if (SW == 4 && (wa % SW == 1))
    {
      out << "vsum += x" << i << "*SIMD<double,4> ("
          << "pa[0*da+" << i*SW << "], "
          << "pa[1*da+" << i*SW << "], "
          << "pa[2*da+" << i*SW << "], "
          << "pa[3*da+" << i*SW << "]);" << endl;
    }

  out << "vsum *= SIMD<double,4>(s);" << endl;
  out << "vsum += SIMD<double,4>(y+i);" << endl;
  out << "vsum.Store(y+i);" << endl;
  out << "}" << endl;


  out << "if (ha & 2) {" << endl;
  out << "SIMD<double," << SW << "> sum0(0.0), sum1(0.0);" << endl;
  i = 0;
  for ( ; SW*(i+1) <= wa; i++)
    {
      out << "sum0 += SIMD<double," << SW << ">(pa+" << i*SW << ") * x" << i << ";" << endl;
      out << "sum1 += SIMD<double," << SW << ">(pa+da+" << i*SW << ") * x" << i << ";" << endl;
    }
  
  if (SW == 4 && (wa % SW == 1))
    {
      for (int k = 0; k < 2; k++)
        out << "sum"<<k<< " += SIMD<double,4> (pa[" << k << "*da+" << i*SW << "] * x" << i << ", 0,0,0);" << endl;      
    }
  else if (SW == 4 && (wa % SW == 2))
    {
      out << "SIMD<double,2> zero(0.0);" << endl;
      out << "sum0 += SIMD<double,4> (SIMD<double,2>(pa+" << i*SW << ") * x" << i << ", zero);" << endl;      
      out << "sum1 += SIMD<double,4> (SIMD<double,2>(pa+da+" << i*SW << ") * x" << i << ", zero);" << endl;      
    }
  else if (wa % SW)
    {
      out << "sum0 += SIMD<double," << SW << ">(pa+" << i*SW << ", mask) * x" << i << ";" << endl;
      out << "sum1 += SIMD<double," << SW << ">(pa+da+" << i*SW << ", mask) * x" << i << ";" << endl;
    }
  out << "SIMD<double,2> vsum = HSum(sum0,sum1);" << endl;
  out << "vsum *= SIMD<double,2>(s);" << endl;
  out << "vsum += SIMD<double,2>(y+i);" << endl;
  out << "vsum.Store(y+i);" << endl;
  out << "i += 2; pa += 2*da;" << endl;
  out << "}" << endl;
  
  
  out << "if (ha & 1) {" << endl;
  out << "SIMD<double," << SW << "> sum(0.0);" << endl;
  i = 0;
  for ( ; SW*(i+1) <= wa; i++)
    out << "sum += SIMD<double," << SW << ">(pa+" << i*SW << ") * x" << i << ";" << endl;

  
  if (SW == 4 && (wa % SW == 1))
    {
      out << "sum += SIMD<double,4> (pa[" << i*SW << "] * x" << i << ", 0,0,0);" << endl;      
    }
  else if (SW == 4 && (wa % SW == 2))
    {
      out << "SIMD<double,2> zero(0.0);" << endl;
      out << "sum += SIMD<double,4> (SIMD<double,2>(pa+" << i*SW << ") * x" << i << ", zero);" << endl;      
    }
  else if (wa % SW)
    out << "sum += SIMD<double," << SW << ">(pa+" << i*SW << ", mask) * x" << i << ";" << endl;
  out << "y[i] += s*HSum(sum);" << endl;

  out << "} }" << endl;
}




void GenerateAddMatTransVecI (ostream & out, int wa)
{
  out << "template <>" << endl
      << "inline void KernelAddMatTransVecI<" << wa << ">" << endl
      << "(double s, size_t ha, double * pa, size_t da, double * x, double * y, int * ind) {" << endl;

  int SW = SIMD_SIZE;  // generate optimal code for my host

  int nfull = wa / SW;
  int rest = wa % SW;
  int unroll = 1;
  if (nfull <= 4) unroll = 2;
  if (nfull <= 1) unroll = 4;

  for (int j = 0; j < unroll; j++)
    {
      for (int i = 0; i < nfull; i++)
        out << "SIMD<double> sy" << i << j << "(0);" << endl;
      if (rest)
        {
          if (rest == 1)
            out << "SIMD<double,1> syrest" << j << "(0);" << endl;
          else if (rest == 2 && SW==4)
            out << "SIMD<double,2> syrest" << j << "(0);" << endl;
          else
            {
              out << "SIMD<double> syrest" << j << "(0);" << endl;
              if (j == 0)
                out << "SIMD<mask64> mask(" << rest << ");" << endl;
            }
        }
    }

  out << "size_t i;" << endl;
  out << "for (i = 0; i+" << unroll << "<= ha; i+=" << unroll << ") {" << endl;
  // out << "for (int j = 0; j < " << unroll << "; j++) {" << endl;
  for (int j = 0; j < unroll; j++)
    {
      out << "SIMD<double> sx" << j << "(x[ind[i+" << j << "]]);" << endl;
      for (int i = 0; i < nfull; i++)
        out << "sy" << i << j << " += SIMD<double>(pa+" << i*SW << ") * sx"<<j<<";" << endl;
      if (rest)
        {
          if (rest == 1)
            out << "syrest" << j << " += SIMD<double,1>(pa+" << nfull*SW << ") * sx"<<j<<"[0];" << endl;
          else if (rest == 2 && SW==4)
            out << "syrest" << j << " += SIMD<double,2>(pa+" << nfull*SW << ") * sx"<<j<<".Lo();" << endl;
          else
            out << "syrest" << j << " += SIMD<double>(pa+" << nfull*SW << ", mask) * sx"<<j<<";" << endl;
        }
      out << "pa += da;" << endl;
    }
  out << "}" << endl;

  /*
  out << "for (  ; i < ha; i++) {" << endl;
  out << "SIMD<double> sx(x[ind[i]]);" << endl;
  for (int i = 0; i < nfull; i++)
    out << "sy" << i << "0 += SIMD<double>(pa+" << i*SW << ") * sx;" << endl;
  if (rest)
    {
      if (rest == 1)
        out << "syrest0 += SIMD<double,1>(pa+" << nfull*SW << ") * sx[0];" << endl;
      else
        out << "syrest0 += SIMD<double>(pa+" << nfull*SW << ", mask) * sx;" << endl;
    }
  out << "pa += da;" << endl;
  out << "}" << endl;
  */

  out << "switch (ha-i) {" << endl;
  for (int j = unroll-1; j >= 1; j--)
    {
      out << "case " << j << ": {" << endl;
      out << "SIMD<double> sx(x[ind[i]]);" << endl;
      for (int i = 0; i < nfull; i++)
        out << "sy" << i << j << " += SIMD<double>(pa+" << i*SW << ") * sx;" << endl;
      if (rest)
        {
          if (rest == 1)
            out << "syrest" << j << " += SIMD<double,1>(pa+" << nfull*SW << ") * sx[0];" << endl;
          else if (rest == 2 && SW==4)
            out << "syrest" << j << " += SIMD<double,2>(pa+" << nfull*SW << ") * sx.Lo();" << endl;
          else
            out << "syrest" << j << " += SIMD<double>(pa+" << nfull*SW << ", mask) * sx;" << endl;
        }
      out << "pa += da; i++; }" << endl;
    }
  out << "default: ;}; " << endl;
  
  if (unroll > 1)
    {
      for (int i = 0; i < nfull; i++)
        {
          out << "sy" << i << "0 += ";
          for (int j = 1; j < unroll-1; j++)
            out << "sy" << i << j << "+";
          out << "sy" << i << unroll-1 << ";" << endl;
        }
      if (rest)
        {
          out << "syrest0 += ";
          for (int j = 1; j < unroll-1; j++)
            out << "syrest" << j << "+";
          out << "syrest" << unroll-1 << ";" << endl;
        }
    }
      
  for (int i = 0; i < nfull; i++)
    out << "(s * sy" << i << "0 + SIMD<double>(y+" << i*SW << ")).Store(y+" << i*SW << ");" << endl;
  if (rest)
    {
      if (rest == 1)
        out << "(s * syrest0 + SIMD<double,1>(y+" << nfull*SW << ")).Store(y+" << nfull*SW << ");" << endl;        
      else if (rest == 2 && SW==4)
        out << "(s * syrest0 + SIMD<double,2>(y+" << nfull*SW << ")).Store(y+" << nfull*SW << ");" << endl;        
      else
        out << "(s * syrest0 + SIMD<double>(y+" << nfull*SW << ",mask)).Store(y+" << nfull*SW << ",mask);" << endl;
    }
  out << "; }" << endl;
}



/* ********************* Triangular kernels ********************************** */

void GenerateTriangular (ofstream & out, bool solve, bool lowerleft, bool normalized, ORDERING order, int dim)
{
  out << "template <> " << endl
      << "inline void " << (solve ? "KernelTriangularSolve" : "KernelTriangularMult")
      << "<" << (lowerleft ? "LowerLeft" : "UpperRight") << ","
      << (normalized ? "Normalized" : "NonNormalized") << ","
      << (order == RowMajor ? "RowMajor" : "ColMajor") << ","
      << dim << ">" 
      << "(size_t wx, double * pt, size_t dt, double * px, size_t dx) { " << endl;

  if (lowerleft)
    {
      for (int i = 0; i < dim; i++)
        for (int j = 0; j < (normalized ? i : i+1); j++)
          if (order == RowMajor)
            out << "SIMD<double> L" << i << j << " = pt[" << i << "*dt+" << j << "];" << endl;
          else
            out << "SIMD<double> L" << i << j << " = pt[" << i << "+dt*" << j << "];" << endl;
    }
  else
    {
      for (int i = 0; i < dim; i++)
        for (int j = (normalized ? i+1 : i); j < dim; j++)
          if (order == RowMajor)          
            out << "SIMD<double> U" << i << j << " = pt[" << i << "*dt+" << j << "];" << endl;
          else
            out << "SIMD<double> U" << i << j << " = pt[" << i << "+dt*" << j << "];" << endl;
    }

  if (solve && !normalized)
    {
      if (lowerleft)
        for (int i = 0; i < dim; i++)
          out << "SIMD<double> Linv" << i << i << " = 1.0 / L" << i << i << ";\n";
      else
        for (int i = 0; i < dim; i++)
          out << "SIMD<double> Uinv" << i << i << " = 1.0 / U" << i << i << ";\n";
    }
  

  stringstream operation;
  if (!solve)
    { // mult
      if (lowerleft)
        { //  lowerleft mult
          for (int i = dim-1; i >= 0; i--)
            {
              if (!normalized)
                operation << "x" << i << " *= L" << i << i << ";\n";
              for (int j = 0; j < i; j++)
                operation << "x" << i << " += L" << i << j << " * x" << j << ";\n";
            }
        }
      else
        {  // upperright mult
          for (int i = 0; i < dim; i++)
            {
              if (!normalized)
                operation << "x" << i << " *= " << "U" << i << i << ";\n";
              for (int j = i+1; j < dim; j++)
                operation << "x" << i << " += U" << i << j << "*" << "x" << j << ";\n";
            }
        }
    }
  else
    { // solve
      if (lowerleft)
        { //  lowerleft solve
          for (int i = 0; i < dim; i++)
            {
              for (int j = 0; j < i; j++)
                operation << "x" << i << " -= L" << i << j << "*x" << j << ";\n";
              if (!normalized)
                operation << "x" << i << " *= " << "Linv" << i << i << ";\n";
            }
        }
      else
        { //  upperright solve
          for (int i = dim-1; i >= 0; i--)
            {
              for (int j = i+1; j < dim; j++)
                operation << "x" << i << " -= U" << i << j << "*x" << j << ";\n";
              if (!normalized)
                operation << "x" << i << " *= " << "Uinv" << i << i << ";\n";
            }
        }
    }

  out << "constexpr size_t SW = SIMD<double>::Size(); \n"
      << "size_t i = 0; \n"
      << "for ( ; i+SW <= wx; i+=SW) { \n";
  
  for (int i = 0; i < dim; i++)
    out << "SIMD<double> x" << i << "(px+" << i << "*dx+i); \n";

  out << operation.str();

  int begin = 0, end = dim;
  if (lowerleft && normalized) begin++;
  if (!lowerleft && normalized) end--;
  for (int i = begin; i < end; i++)
    out << "x" << i << ".Store(px+" << i << "*dx+i); \n";
  
  out << "}\n";

  // remainder

  // mask = ...
  out << "size_t rest = wx % SW; \n"
      << "if (rest == 0) return; \n"
      << "SIMD<mask64> mask(rest);" << endl;
  
  for (int i = 0; i < dim; i++)
    out << "SIMD<double> x" << i << "(px+" << i << "*dx+i, mask); \n";

  out << operation.str();

  for (int i = begin; i < end; i++)    
    out << "x" << i << ".Store(px+" << i << "*dx+i, mask); \n";

  out << "}\n";
}



void GenerateTriangularXY (ofstream & out, bool solve, bool lowerleft, bool normalized, ORDERING order,
                           OP op, int dim)
{
  out << "template <> " << endl
      << "inline void " << (solve ? "KernelTriangularSolveXY" : "KernelTriangularMultXY")
      << "<" << (lowerleft ? "LowerLeft" : "UpperRight") << ","
      << (normalized ? "Normalized" : "NonNormalized") << ","
      << (order == RowMajor ? "RowMajor" : "ColMajor") << ","
      << ToString(op) << ", "
      << dim << ">" 
      << "(size_t wx, double * pt, size_t dt, double * px, size_t dx, double * py, size_t dy) { " << endl;

  if (lowerleft)
    {
      for (int i = 0; i < dim; i++)
        for (int j = 0; j < (normalized ? i : i+1); j++)
          if (order == RowMajor)
            out << "SIMD<double> L" << i << j << " = pt[" << i << "*dt+" << j << "];" << endl;
          else
            out << "SIMD<double> L" << i << j << " = pt[" << i << "+dt*" << j << "];" << endl;
      if (normalized)
        for (int i = 0; i < dim; i++)
          out << "SIMD<double> L" << i << i << "(1.0);\n";
    }
  else
    {
      for (int i = 0; i < dim; i++)
        for (int j = (normalized ? i+1 : i); j < dim; j++)
          if (order == RowMajor)          
            out << "SIMD<double> U" << i << j << " = pt[" << i << "*dt+" << j << "];" << endl;
          else
            out << "SIMD<double> U" << i << j << " = pt[" << i << "+dt*" << j << "];" << endl;
      if (normalized)
        for (int i = 0; i < dim; i++)
          out << "SIMD<double> U" << i << i << "(1.0);\n";
    }

  if (solve && !normalized)
    {
      if (lowerleft)
        for (int i = 0; i < dim; i++)
          out << "SIMD<double> Linv" << i << i << " = 1.0 / L" << i << i << ";\n";
      else
        for (int i = 0; i < dim; i++)
          out << "SIMD<double> Uinv" << i << i << " = 1.0 / U" << i << i << ";\n";
    }
  

  stringstream operation;
  string addop = (op == SET || op == ADD) ? "+=" : "-=";
  if (!solve)
    { // mult
      if (lowerleft)
        { //  lowerleft mult
          for (int i = dim-1; i >= 0; i--)
            {
              for (int j = 0; j <= i; j++)
                operation << FMAOp (op, false) << "(L" << i << j << ", x" << j << ", y" << i << ");\n";
                // operation << "y" << i << addop << "L" << i << j << " * x" << j << ";\n";
            }
        }
      else
        {  // upperright mult
          for (int i = 0; i < dim; i++)
            {
              for (int j = i; j < dim; j++)
                operation << FMAOp (op, false) << "(U" << i << j << ", x" << j << ", y" << i << ");\n";                
              // operation << "y" << i << addop << "U" << i << j << "*" << "x" << j << ";\n";
            }
        }
    }
  else
    {
      throw std::runtime_error ("solvetrig, xy not implemented");
      // solve
      if (lowerleft)
        { //  lowerleft solve
          for (int i = 0; i < dim; i++)
            {
              for (int j = 0; j < i; j++)
                operation << "x" << i << " -= L" << i << j << "*x" << j << ";\n";
              if (!normalized)
                operation << "x" << i << " *= " << "Linv" << i << i << ";\n";
            }
        }
      else
        { //  upperright solve
          for (int i = dim-1; i >= 0; i--)
            {
              for (int j = i+1; j < dim; j++)
                operation << "x" << i << " -= U" << i << j << "*x" << j << ";\n";
              if (!normalized)
                operation << "x" << i << " *= " << "Uinv" << i << i << ";\n";
            }
        }
    }

  out << "constexpr size_t SW = SIMD<double>::Size(); \n"
      << "size_t i = 0; \n"
      << "for ( ; i+SW <= wx; i+=SW) { \n";
  
  for (int i = 0; i < dim; i++)
    {
      out << "SIMD<double> x" << i << "(px+" << i << "*dx+i); \n";
      if (op == SET || op == SETNEG)
        out << "SIMD<double> y" << i << "(0);\n";
      else
        out << "SIMD<double> y" << i << "(py+" << i << "*dy+i); \n";
    }

  out << operation.str();

  for (int i = 0; i < dim; i++)
    out << "y" << i << ".Store(py+" << i << "*dy+i); \n";
  
  out << "}\n";

  // remainder

  // mask = ...
  out << "size_t rest = wx % SW; \n"
      << "if (rest == 0) return; \n"
      << "SIMD<mask64> mask(rest);" << endl;
  
  for (int i = 0; i < dim; i++)
    {
      out << "SIMD<double> x" << i << "(px+" << i << "*dx+i, mask); \n";
      if (op == SET || op == SETNEG)
        out << "SIMD<double> y" << i << "(0);\n";
      else
        out << "SIMD<double> y" << i << "(py+" << i << "*dy+i, mask); \n";
    }
  out << operation.str();

  for (int i = 0; i < dim; i++)    
    out << "y" << i << ".Store(py+" << i << "*dy+i, mask); \n";

  out << "}\n";
}








int main (int argn, char **argv)
{
  ofstream out(argv[1]);


  out <<
R"raw_string(
template <typename Tuple, std::size_t ... Is>
auto pop_front2_impl(const Tuple& tuple, std::index_sequence<Is...>)
{
    return std::make_tuple(std::get<2 + Is>(tuple)...);
}

template <typename Tuple>
auto pop_front2(const Tuple& tuple)
{
    return pop_front2_impl(tuple,
                          std::make_index_sequence<std::tuple_size<Tuple>::value - 2>());
}

template<int N1, int N2>
auto Concat2(SIMD<double,N1> simd1, SIMD<double,N2> simd2)
{
    if constexpr (IsNativeSIMDSize(simd1.Size())) {
        return SIMD<double,simd1.Size()+simd2.Size()>(simd1, simd2);
    }
    else {
      auto lo = simd1.Lo();
      auto hi = simd1.Hi();
  
      SIMD<double,hi.Size()+simd2.Size()> res1(hi, simd2);
      SIMD<double,simd1.Size()+simd2.Size()> res2(lo, res1);
      return res2;
    }
}

template <typename ...Args, int N>
auto Concat (tuple<SIMD<double,N>, Args...> tup)
{
  if constexpr (tuple_size<tuple<SIMD<double,N>, Args...>>() == 1)
                 return get<0>(tup);
  else if constexpr (tuple_size<tuple<SIMD<double,N>, Args...>>() == 2)
      return Concat2(get<0>(tup), get<1>(tup));
  else
    {
      auto front = Concat2(get<0>(tup), get<1>(tup));
      auto rest = pop_front2(tup);
      return Concat(std::tuple_cat(make_tuple(front), rest));
    }
}
)raw_string"
      << endl;


  
  out << "template <int N>\n"
    "void FMAnonasm (SIMD<double,N> a, SIMD<double,N> b, SIMD<double,N> & sum)\n"
    "{ sum = FMA(a,b,sum); } ";
  out << "template <int N>\n"
    "void FNMAnonasm (SIMD<double,N> a, SIMD<double,N> b, SIMD<double,N> & sum)\n"
    "{ sum = FNMA(a,b,sum); }";

  
  out << "static_assert(SIMD<double>::Size() == " << SIMD_SIZE << ", \"inconsistent compile flags for generate_mat_kernels.cpp and matkernel.hpp\");" << endl;
  out << "enum OPERATION { ADD, SUB, SET, SETNEG };" << endl;

  out << " /* *********************** MatKernelMultAB ********************* */" << endl
      << " /* A,B,C ... row major storage                                   */" << endl
      << " /* dim C = H * (SW*W)     SW .. SIMD<double>::Size()             */" << endl
      << " /* OP == SET:    C = A * B                                       */" << endl
      << " /* OP == ADD:    C += A * B                                      */" << endl
      << " /* OP == NEG:    C = -A * B                                      */" << endl
      << " /* OP == SUB:    C -= A * B                                      */" << endl
      << " /* ************************************************************* */" << endl;   
  out << "template <size_t H, size_t W, OPERATION OP>" << endl
      << "inline void MatKernelMultAB" << endl
      << "(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);" << endl;
  out << "template <size_t H, size_t W, OPERATION OP>" << endl
      << "inline void MatKernelMultAB" << endl
      << "(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, double * pc, size_t dc);" << endl;
  out << "template <size_t H, size_t W>" << endl
      << "inline void MatKernelAlignedMultAB" << endl
      << "(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, SIMD<double> * pc, size_t dc);" << endl;

  
  for (int i = 1; i <= 4; i++)
    {
      GenerateMultAB (out, 1, i);  
      GenerateMultAB (out, 2, i);
      GenerateMultAB (out, 3, i);
      GenerateMultAB (out, 4, i);
      GenerateMultAB (out, 5, i);
      GenerateMultAB (out, 6, i);
      
      AlignedGenerateMultAB (out, 1, i, SET);  
      AlignedGenerateMultAB (out, 2, i, SET);
      AlignedGenerateMultAB (out, 3, i, SET);
      AlignedGenerateMultAB (out, 4, i, SET);
      AlignedGenerateMultAB (out, 5, i, SET);
      AlignedGenerateMultAB (out, 6, i, SET);
    }
  
  GenerateMultAB (out, 8, 1);
  GenerateMultAB (out, 12, 1);
  

  out << "template <size_t H, OPERATION OP>" << endl
      << "inline void MatKernelMultABMask" << endl
      << "(size_t n, SIMD<mask64> mask, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);" << endl;
  out << "template <size_t H, OPERATION OP>" << endl
      << "inline void MatKernelMultABMask" << endl
      << "(size_t n, SIMD<mask64> mask, double * pa, size_t da, SIMD<double> * pb, size_t db, double * pc, size_t dc);" << endl;

  GenerateMultABMask (out, 1);  
  GenerateMultABMask (out, 2);
  GenerateMultABMask (out, 3);
  GenerateMultABMask (out, 4);
  GenerateMultABMask (out, 5);
  GenerateMultABMask (out, 6);

  
  // Scal AB
  out << " /* *********************** MatKernelScalAB ********************* */" << endl
      << " /* Inner products of rows of A with rows of B                    */" << endl
      << " /* A,B ... row major storage                                     */" << endl
      << " /* dim A = H * n                                                 */" << endl
      << " /* dim B = W * n                                                 */" << endl
      << " /* ************************************************************* */" << endl;
  
  out << "template <size_t H, size_t W> inline auto MatKernelScalAB" << endl
      << "    (size_t n," << endl
      << "     double * pa, size_t da," << endl
      << "     double * pb, size_t db);" << endl;
  out << "template <size_t H, size_t W> inline auto MatKernelScalAB" << endl
      << "    (size_t n," << endl
      << "     SIMD<double> * pa, size_t da," << endl
      << "     SIMD<double> * pb, size_t db);" << endl;

  GenerateScalAB (out, 6, 4);  
  GenerateScalAB (out, 3, 4);  
  GenerateScalAB (out, 1, 4);
  GenerateScalAB (out, 6, 2);  
  GenerateScalAB (out, 3, 2);  
  GenerateScalAB (out, 8, 1);  
  GenerateScalAB (out, 6, 1);  
  GenerateScalAB (out, 4, 1);  
  GenerateScalAB (out, 3, 1);  
  GenerateScalAB (out, 2, 1);  
  GenerateScalAB (out, 1, 1);  
  
  
    // MultiVecScalAB

  out << "template <size_t H, size_t W> inline auto MultiVecScalAB" << endl
      << "    (size_t n," << endl
      << "     double ** ppa," << endl
      << "     double ** ppb);" << endl;

  GenerateMultiVecScalAB (out, 6, 4);
  GenerateMultiVecScalAB (out, 3, 4);
  GenerateMultiVecScalAB (out, 1, 4);
  GenerateMultiVecScalAB (out, 6, 2);
  GenerateMultiVecScalAB (out, 3, 2);
  GenerateMultiVecScalAB (out, 8, 1);
  GenerateMultiVecScalAB (out, 6, 1);
  GenerateMultiVecScalAB (out, 4, 1);
  GenerateMultiVecScalAB (out, 3, 1);
  GenerateMultiVecScalAB (out, 2, 1);
  GenerateMultiVecScalAB (out, 1, 1);


  // MultiVecScalC
  out << "template <size_t H, size_t W, bool conjugate> inline void MultiVecScalC" << endl
      << "    (size_t n," << endl
      << "     Complex ** ppa," << endl
      << "     Complex ** ppb," << endl
      << "     Complex* pc, size_t dc);" << endl;

  GenerateMultiVecScalC(out, 6, 2, 0);
  GenerateMultiVecScalC(out, 6, 1, 0);
  GenerateMultiVecScalC(out, 3, 2, 0);
  GenerateMultiVecScalC(out, 3, 1, 0);
  GenerateMultiVecScalC(out, 1, 8, 0);
  GenerateMultiVecScalC(out, 1, 4, 0);
  GenerateMultiVecScalC(out, 1, 1, 0);
  GenerateMultiVecScalC(out, 6, 2, 1);
  GenerateMultiVecScalC(out, 6, 1, 1);
  GenerateMultiVecScalC(out, 3, 2, 1);
  GenerateMultiVecScalC(out, 3, 1, 1);
  GenerateMultiVecScalC(out, 1, 8, 1);
  GenerateMultiVecScalC(out, 1, 4, 1);
  GenerateMultiVecScalC(out, 1, 1, 1);




    // MultiScaleAdd
  out << "template <size_t H, size_t W> inline void MultiScaleAdd" << endl
      << "    (size_t n," << endl
      << "     double ** pa, " << endl
      << "     double ** pb, " << endl
      << "     double * pc, size_t dc);" << endl;

  GenerateMultiScaleAdd(out, 6, 6);
  GenerateMultiScaleAdd(out, 2, 6);
  GenerateMultiScaleAdd(out, 1, 6);
  GenerateMultiScaleAdd(out, 6, 2);
  GenerateMultiScaleAdd(out, 2, 2);
  GenerateMultiScaleAdd(out, 1, 2);
  GenerateMultiScaleAdd(out, 6, 1);
  GenerateMultiScaleAdd(out, 2, 1);
  GenerateMultiScaleAdd(out, 1, 1);


    // MultiScaleAddC
  out << "template <size_t H, size_t W> inline void MultiScaleAddC" << endl
      << "    (size_t n," << endl
      << "     Complex ** pa, " << endl
      << "     Complex ** pb, " << endl
      << "     Complex * pc, size_t dc);" << endl;


  GenerateMultiScaleAddC (out, 3, 4);
  GenerateMultiScaleAddC (out, 2, 4);
  GenerateMultiScaleAddC (out, 1, 4);
  GenerateMultiScaleAddC (out, 3, 1);
  GenerateMultiScaleAddC (out, 2, 1);
  GenerateMultiScaleAddC (out, 1, 1);


  
  out << "template <size_t H, size_t W>" << endl
      << "inline void MyScalTrans" << endl
      << "(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);" << endl;
  
  GenKernel (out, 1, 4);
  GenKernel (out, 2, 4);
  GenKernel (out, 3, 4);
  GenKernel (out, 4, 4);
  GenKernel (out, 5, 4);
  GenKernel (out, 6, 4);

  out << "template <size_t H, size_t W, OPERATION OP>" << endl
      << "inline void MatKernelDaxpy" << endl
      << "(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);" << endl;
  out << "template <size_t H, size_t W, OPERATION OP>" << endl
      << "inline void MatKernelDaxpy" << endl
      << "(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, SIMD<double> * pc, size_t dc);" << endl;

  for (int i = 0; i <= 12; i++)
    GenerateDaxpy (out, 1, i);

  GenerateDaxpy (out, 2, 1);  
  GenerateDaxpy (out, 2, 2);  
  GenerateDaxpy (out, 2, 3);  
  GenerateDaxpy (out, 2, 4);  
  GenerateDaxpy (out, 3, 1);  
  GenerateDaxpy (out, 3, 2);  
  GenerateDaxpy (out, 3, 3);  
  GenerateDaxpy (out, 3, 4);

  out << "// C = A * B,  with short inner loop\n"
      << "template <size_t WA, OPERATION OP>\n"
      << "inline void MatKernelShortSum\n"
      << "(size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);\n";

  out << "// C = A * B,  with short inner loop, unroll width B\n"
      << "template <size_t WA, OPERATION OP>\n"
      << "inline void MatKernelShortSum2\n"
      << "(size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);\n";

  for (int i = 0; i <= 12; i++)
    for (auto op : { SET, SETNEG, ADD, SUB })
      GenerateShortSum (out, i, op);  


  out << "// C = A^t * B,  with short inner loop\n"
      << "template <size_t WA, OPERATION OP>\n"
      << "inline void MatKernelAtB_SmallWA\n"
      << "(size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);\n";

  out << "// C = A^t * B,  with short inner loop, unroll width B\n"
      << "template <size_t WA, OPERATION OP>\n"
      << "inline void MatKernelAtB_SmallWA2\n"
      << "(size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc);\n";

  for (int i = 0; i <= 12; i++)
    for (OP op : { SET, ADD, SUB, SETNEG })
      GenerateAtB_SmallWA (out, i, op);
  
  for (int i = 0; i <= 8; i++)    
    for (OP op : { SET, ADD, SUB, SETNEG })
      GenerateAtB_SmallWA2 (out, i, op);
  
  out << "// y = A * x,  with fix width" << endl;
  out << "template <size_t WA, OPERATION OP>" << endl
      << "inline void KernelMatVec" << endl
      << "(size_t ha, double * pa, size_t da, double * x, double * y);" << endl;
  for (int i = 0; i <= 24; i++)
    {
      GenerateMatVec (out, i, SET);
      GenerateMatVec (out, i, ADD);
    }

  out << "// y += s * A * x,  with fix width" << endl;
  out << "template <size_t WA>" << endl
      << "inline void KernelAddMatVec" << endl
      << "(double s, size_t ha, double * pa, size_t da, double * x, double * y);" << endl;
  for (int i = 0; i <= 24; i++)
    GenerateAddMatVec (out, i);


  out << "// y += s * A^t * x(ind),  with fix width" << endl;
  out << "template <size_t WA>" << endl
      << "inline void KernelAddMatTransVecI" << endl
      << "(double s, size_t ha, double * pa, size_t da, double * x, double * y, int * ind);" << endl;
  for (int i = 0; i <= 24; i++)
    GenerateAddMatTransVecI (out, i);




  /* *********************** MatKernelTriangularMult ***************** */
  out << "template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, int DIM>" << endl
      << "inline void KernelTriangularMult (size_t wx, double * pt, size_t dt, double * px, size_t dx);" << endl;
  out << "template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, OPERATION OP, int DIM>" << endl
      << "inline void KernelTriangularMultXY (size_t wx, double * pt, size_t dt, "
      << "double * px, size_t dx, double * py, size_t dy) { " << endl
      << "cerr << \"missing implementation, side = \" << SIDE << \" norm = \" << NORM" << endl
      << "<< \" ORD = \" << ORD << \", OP = \" << OP << \", DIM = \" << DIM << endl; }" << endl;
  out << "template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, int DIM>" << endl
      << "inline void KernelTriangularSolve (size_t wx, double * pt, size_t dt, double * px, size_t dx);" << endl;


  for (int i = 0; i <= 4; i++)
    for (bool lowerleft : { false, true })
      for (bool normalized : { false, true })
        {
          // bool solve, bool lowerleft, bool normalized, bool trans, int dim
          
          GenerateTriangular (out, false, lowerleft, normalized, RowMajor, i);
          GenerateTriangular (out, false, lowerleft, normalized, ColMajor, i);
          GenerateTriangular (out, true, lowerleft, normalized, RowMajor, i);

          GenerateTriangularXY (out, false, lowerleft, normalized, RowMajor, SET, i);
          GenerateTriangularXY (out, false, lowerleft, normalized, ColMajor, SET, i);
          GenerateTriangularXY (out, false, lowerleft, normalized, RowMajor, ADD, i);
          GenerateTriangularXY (out, false, lowerleft, normalized, ColMajor, ADD, i);
          GenerateTriangularXY (out, false, lowerleft, normalized, RowMajor, SUB, i);
          GenerateTriangularXY (out, false, lowerleft, normalized, ColMajor, SUB, i);
          
          /*
          GenerateTriangular (out, false, true, true, RowMajor, i);
          GenerateTriangular (out, false, false, false, RowMajor, i);
          GenerateTriangular (out, false, false, true, RowMajor, i);
          GenerateTriangular (out, false, true, false, ColMajor, i);
          GenerateTriangular (out, false, true, true, ColMajor, i);
          GenerateTriangular (out, false, false, false, ColMajor, i);
          GenerateTriangular (out, false, false, true, ColMajor, i);
          GenerateTriangular (out, true, true, false, RowMajor, i);
          GenerateTriangular (out, true, true, true, RowMajor, i);
          GenerateTriangular (out, true, false, false, RowMajor, i);
          GenerateTriangular (out, true, false, true, RowMajor, i);
          */
        }
}
