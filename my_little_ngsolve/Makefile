NGSCXX = /opt/netgen-5.3/bin/ngscxx

objects = demo_instat.o myElement.o myHOElement.o myIntegrator.o \
demo_coupling.o myFESpace.o myHOFESpace.o myPreconditioner.o	 \
myAssembling.o linhypDG.o all_in_one.o demo_stokes.o		 \
demo_coupling_adv.o demo_nonlinear.o 



%.o : %.cpp
	$(NGSCXX) -I. -c $? -o $@


libmyngsolve.so : $(objects)
	$(NGSCXX) -shared $(objects) -lngfem -lngcomp -lngsolve -o $@

clean:
	rm *.o libmyngsolve.so



dist:
	tar --transform='s,^,my_little_ngsolve-5.3/,' -czf MyLittleNGSolve-5.3.tar.gz Makefile *.cpp *.hpp *.pde *.vol *.in2d windows/*.vcxproj windows/*.sln windows/*.bat */Makefile */*.cpp  */*.pde  */*.vol*  */*.geo  */*.in2d  */*.tcl


#	cd ..; \
#	tar -czf MyLittleNGSolve-5.3.tar.gz my_little_ngsolve-5.3/Makefile my_little_ngsolve-5.3/*.cpp my_little_ngsolve-5.3/*.hpp my_little_ngsolve-5.3/*.pde my_little_ngsolve-5.3/*.vol my_little_ngsolve-5.3/*.in2d my_little_ngsolve-5.3/windows/my_little_ngsolve* my_little_ngsolve-5.3/*/Makefile my_little_ngsolve-5.3/*/*.cpp  my_little_ngsolve-5.3/*/*.pde  my_little_ngsolve-5.3/*/*.vol*  my_little_ngsolve-5.3/*/*.geo  my_little_ngsolve-5.3/*/*.in2d  my_little_ngsolve-5.3/*/*.tcl

