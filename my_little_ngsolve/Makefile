objects = all_in_one.o demo_instat.o demo_stokes.o myElement.o \
myHOElement.o myIntegrator.o demo_coupling.o demo_coupling_adv.o demo_nonlinear.o  \
myFESpace.o myHOFESpace.o myPreconditioner.o myAssembling.o linhypDG.o 


%.o : %.cpp
	gcc  -O2 -fopenmp -fpic -DUSE_TIMEOFDAY -DLAPACK -I. -I$(NETGENDIR)/../include -c $? -o $@

libmyngsolve.so : $(objects)
	gcc -shared -fopenmp -fpic $(objects) -L$(NETGENDIR) -lngsolve -o $@

clean:
	rm *.o libmyngsolve.so



dist:
	cd ..; tar -czf MyLittleNGSolve-4.9.14.tar.gz my_little_ngsolve/Makefile my_little_ngsolve/*.cpp my_little_ngsolve/*.hpp my_little_ngsolve/*.pde my_little_ngsolve/*.vol my_little_ngsolve/*.in2d my_little_ngsolve/windows/my_little_ngsolve*
