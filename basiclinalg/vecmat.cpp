#include <bla.hpp>

namespace ngbla
{
  void test()
  {
    Array<int> ia(3);
    Vector<Vec<2> > hv(10);
    cout << hv(ia) << endl;
  }



  void CheckMatRange(int h, int w, int i)
  {
    if (i < 0 || i > h*w)
      {
        stringstream st;
        st << "FlatMatrix: " << i << "out of range (0, " << h*w << ")" << endl;
        throw Exception (st.str());
      }
  }

  void CheckMatRange(int h, int w, int i, int j)
  {
    if (i < 0 || i > h || j < 0 || j > w)
      {
	stringstream st;
	st << "FlatMatrix: (" << i << ", " << j << ") out of range (0, " << h << ") x (0, " << w << ")" << endl;
	throw Exception (st.str());
      }
  }

  void CheckVecRange(int s, int i) 
  {
    if (i < 0 || i > s)
      {
	stringstream st;
	st << "Vec: " << i << "out of range (0, " << s << ")" << endl;
	throw Exception (st.str());
      }
  }
  
  void CheckVecRange(int s, int i, int j) 
  {
    if (i < 0 || i > s || j != 0)
      {
	stringstream st;
	st << "Vec: (" << i << ", " << j << ") out of range (0, " << s << ") x (0,0)" << endl;
	throw Exception (st.str());
      }
  }
}
