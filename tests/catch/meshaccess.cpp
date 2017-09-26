
#include "catch.hpp"
#include <comp.hpp>

using namespace ngcomp;


TEST_CASE ("ElementVolume", "[elementvolume]")
{
  SECTION ("ElementVolume", "[elementvolume]")
    {
      for(int dim : {2,3})
	{
	  SECTION ("ElementVolume, dim = "+to_string(dim), "[elementvolume]")
	    {
	      string meshfile = dim==2 ? "square.vol" : "cube.vol";
	      double elvol = dim==2 ? 1./2 : 1./12;
	      auto ma = MeshAccess(meshfile);
	      for(auto el : ma.Elements(VOL))
		{
		  CHECK(ma.ElementVolume(el.Nr()) == Approx(elvol).epsilon(1e-9));
		  // CHECK(fabs(ma.ElementVolume(el.Nr())-elvol) < 1e-9);
		}
	    }
	}
    }
}

