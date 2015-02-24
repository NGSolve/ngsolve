objects = all_in_one.o demo_instat.o demo_stokes.o myElement.o \
myHOElement.o myIntegrator.o demo_coupling.o demo_coupling_adv.o demo_nonlinear.o  \
myFESpace.o myHOFESpace.o myPreconditioner.o myAssembling.o linhypDG.o 

# linhypDG_par.o 


%.o : %.cpp
	ngscxx -I. -c $? -o $@


libmyngsolve.so : $(objects)
	ngscxx -shared $(objects) -lngsolve -o $@

clean:
	rm *.o libmyngsolve.so



dist:
	cd ..; tar -czf MyLittleNGSolve-5.1.tar.gz my_little_ngsolve-5.1/Makefile my_little_ngsolve-5.1/*.cpp my_little_ngsolve-5.1/*.hpp my_little_ngsolve-5.1/*.pde my_little_ngsolve-5.1/*.vol my_little_ngsolve-5.1/*.in2d my_little_ngsolve-5.1/windows/my_little_ngsolve* \
my_little_ngsolve-5.1/mixed/*.cpp my_little_ngsolve-5.1/mixed/*.hpp my_little_ngsolve-5.1/mixed*.pde \
my_little_ngsolve-5.1/mixed*.vol my_little_ngsolve-5.1/mixed*.in2d \
my_little_ngsolve-5.1/precond/*.cpp my_little_ngsolve-5.1/precond/*.hpp my_little_ngsolve-5.1/precond*.pde \
my_little_ngsolve-5.1/precond*.vol my_little_ngsolve-5.1/precond*.in2d \
my_little_ngsolve-5.1/tcl/*.cpp my_little_ngsolve-5.1/tcl/*.hpp my_little_ngsolve-5.1/tcl*.pde \
my_little_ngsolve-5.1/tcl*.vol my_little_ngsolve-5.1/tcl*.in2d \
my_little_ngsolve-5.1/windows/*.sln my_little_ngsolve-5.1/windows/*.vcproj
