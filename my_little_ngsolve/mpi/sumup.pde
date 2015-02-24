mesh = square.vol.gz
# mesh = cube.vol.gz

shared = mymip

fespace v -type=nodal -order=1
gridfunction u -fespace=v

numproc SumUp np1
numproc SumUp np2



