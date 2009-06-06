#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include "opti.hpp"

namespace netgen
{

void LinearOptimize (const DenseMatrix & a, const Vector & b, 
    const Vector & c, Vector & x)
    
  {
  int i1, i2, i3, j;
  DenseMatrix m(3), inv(3);
  Vector rs(3), hx(3), res(a.Height()), res2(3);
  double f, fmin;
  int nrest;
    
  if (a.Width() != 3)
    {
    cerr << "LinearOptimize only implemented for 3 unknowns" << endl;
    return;
    }
    
  fmin = 1e10;
  x = 0;
  nrest = a.Height();
  for (i1 = 1; i1 <= nrest; i1++)
    for (i2 = i1 + 1; i2 <= nrest; i2++)
      for (i3 = i2 + 1; i3 <= nrest; i3++)
        {
        for (j = 1; j <= 3; j++)
          {
          m.Elem(1, j) = a.Get(i1, j);
          m.Elem(2, j) = a.Get(i2, j);
          m.Elem(3, j) = a.Get(i3, j);
          }
          
        rs(0) = b(i1-1);
        rs(1) = b(i2-1);
        rs(2) = b(i3-1);
        
        if (fabs (m.Det()) < 1e-12) continue;
        
        CalcInverse (m, inv);
        inv.Mult (rs, hx);
        
        a.Residuum (hx, b, res);
//        m.Residuum (hx, rs, res2);
        f = c * hx;

/*        
        testout -> precision(12);
        (*testout) << "i = (" << i1 << "," << i2 << "," << i3 
           << "), f = " << f << " x = " << x << " res = " << res 
           <<  " resmin = " << res.Min() 
           << " res2 = " << res2 << " prod = " << prod << endl;
*/

	
	double rmin = res(0);
	for (int hi = 1; hi < res.Size(); hi++)
	  if (res(hi) < rmin) rmin = res(hi);
        
        if ( (f < fmin) && rmin >= -1e-8)
          {
          fmin = f;
          x = hx;
          }
        }
  }
}
