objects = all_in_one.o demo_instat.o demo_stokes.o myElement.o	 \
myHOElement.o myIntegrator.o demo_coupling.o demo_coupling_adv.o \
demo_nonlinear.o myFESpace.o myHOFESpace.o myPreconditioner.o	 \
myAssembling.o linhypDG.o
# linhypDG_par.o 


%.o : %.cpp
	ngscxx -I. -c $? -o $@


libmyngsolve.so : $(objects)
	ngscxx -shared $(objects) -lngfem -lngcomp -lngsolve -o $@

clean:
	rm *.o libmyngsolve.so



dist:
	cd ..; tar -czf MyLittleNGSolve-5.2.tar.gz my_little_ngsolve-5.2/Makefile my_little_ngsolve-5.2/*.cpp my_little_ngsolve-5.2/*.hpp my_little_ngsolve-5.2/*.pde my_little_ngsolve-5.2/*.vol my_little_ngsolve-5.2/*.in2d my_little_ngsolve-5.2/windows/my_little_ngsolve*
