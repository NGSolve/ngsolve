
#include "catch.hpp"
#include <comp.hpp>

using namespace ngcomp;

TEST_CASE ("ElementVolume", "[elementvolume]")
{
  netgen::printmessage_importance = 0;
  SECTION ("ElementVolume", "[elementvolume]")
    {
      for(int dim : {1,2,3})
	{
	  SECTION ("ElementVolume, dim = "+to_string(dim), "[elementvolume]")
	    {
	      string meshfile = dim==1 ? "line.vol" : (dim==2 ? "square.vol" : "cube.vol");
	      double elvol = dim==1 ? 1.0 : (dim==2 ? 1./2 : 1./12);
	      auto ma = MeshAccess(meshfile);
	      for(auto el : ma.Elements(VOL))
		{
		  CHECK(ma.ElementVolume(el.Nr()) == Approx(elvol).epsilon(1e-9));
		}
	    }
	}
    }
}

TEST_CASE ("Region")
{
  SECTION("Constructors")
    {
      auto ma = make_shared<MeshAccess>("2_doms.vol");
      auto nreg = ma->GetNRegions(VOL);
      CHECK(nreg == 2);
      BitArray compare(nreg);
      compare.Clear();
      auto none = Region(ma, VOL);
      CHECK(none.Mask() == compare);
      auto inner = Region(ma, VOL, "inner");
      compare.SetBit(1);
      CHECK(inner.Mask() == compare);
      auto outer = Region(ma, VOL, "outer");
      compare.Invert();
      CHECK(outer.Mask() == compare);
      auto all = Region(ma, VOL, true);
      compare.Set();
      CHECK(all.Mask() == compare);
    }
}

