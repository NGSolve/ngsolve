from ngsolve import *
from netgen.csg import unit_cube


mesh=Mesh(unit_cube.GenerateMesh(maxh=0.1))
pml_brick=pml.BrickRadial((0.2,0.2,0.2),(0.4,0.4,0.4),1j,origin=(0.3,0.3,0.3))
pml_halfspace=pml.HalfSpace((0.2,0.2,0.2),(2,1,1),1j)
pml_halfspace2=pml.HalfSpace((0.7,0.7,0.7),(-2,1,1),1j)
pml_cartesian=pml.Cartesian((0.2,0.2,0.2),(0.5,0.5,0.5),1j)
pml_radial=pml.Radial(rad=0.3,alpha=1j,origin=(0.3,0.3,0.3))
pml_cylinder=pml.Compound(pml_radial,pml_cartesian,[1,2])


cf_brick=pml_brick.PML_CF(3)
cf_halfspace=pml_halfspace.PML_CF(3)
cf_cartesian=pml_cartesian.PML_CF(3)
cf_radial=pml_radial.PML_CF(3)
cf_cylinder=pml_cylinder.PML_CF(3)
cf_sum=(pml_halfspace+pml_halfspace2).PML_CF(3)

Draw(cf_brick-Conj(cf_brick),mesh,"brick")
Draw(cf_halfspace-Conj(cf_halfspace),mesh,"halfspace")
Draw(cf_cartesian-Conj(cf_cartesian),mesh,"cartesian")
Draw(cf_radial-Conj(cf_radial),mesh,"radial")
Draw(cf_cylinder-Conj(cf_cylinder),mesh,"cylinder")
Draw(cf_sum-Conj(cf_sum),mesh,"sum")
