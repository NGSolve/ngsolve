#include <ngstd.hpp>
#include <bla.hpp>




using namespace ngstd;
using namespace ngbla;

class Mesh1D_Base
{
public:
  int ne;
  Array<INT<1,double>> points;
  Array<INT<2>> els;
  Array<INT<1>> bound_els;
};

extern Mesh1D_Base mesh1d;




namespace ngcomp
{
#ifdef FILE_MESHACCESS
  class Mesh1D : public MeshAccess
  {
  public:
    Mesh1D (int ane)
    { 
      mesh1d.ne = ane;

      mesh1d.els.SetSize(ane);
      for (int i = 0; i < ane; i++)
        {
          mesh1d.els[i][0] = i+2;
          mesh1d.els[i][1] = i+1;
        }


      mesh1d.points.SetSize(ane+1);
      for (int i = 0; i <= ane; i++)
        mesh1d.points[i][0] = double(i)/ane;

      mesh1d.bound_els.SetSize(2);
      mesh1d.bound_els[0][0] = 1;
      mesh1d.bound_els[1][0] = ane+1;

      UpdateBuffers();
    }
  };
#endif  

  
}
