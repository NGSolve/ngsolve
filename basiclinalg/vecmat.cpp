#include <bla.hpp>

namespace ngbla
{

  void CheckMatRange(size_t h, size_t w, size_t i)
  {
    if ( /* i < 0 || */ i > h*w)
      {
        stringstream st;
        st << "Matrix: " << i << "out of range [0, " << h*w << ")" << endl;
        throw Exception (st.str());
      }
  }

  void CheckMatRange(size_t h, size_t w, size_t i, size_t j)
  {
    if ( /* i < 0 || */ i > h ||  /* j < 0 || */j > w)
      {
	stringstream st;
	st << "FlatMatrix: (" << i << ", " << j << ") out of range [0, " << h << ") x [0, " << w << ")" << endl;
	throw Exception (st.str());
      }
  }

  void CheckVecRange(size_t s, size_t i) 
  {
    if (/* i < 0 || */ i > s)
      {
	stringstream st;
	st << "Vec: " << i << "out of range [0, " << s << ")" << endl;
	throw Exception (st.str());
      }
  }
  
  void CheckVecRange(size_t s, size_t i, size_t j) 
  {
    if (/* i < 0 || */ i > s || j != 0)
      {
	stringstream st;
	st << "Vec: (" << i << ", " << j << ") out of range [0, " << s << ") x [0,1)" << endl;
	throw Exception (st.str());
      }
  }
}
