objects = asm_precond.o 

%.o : %.cpp
	ngscxx -I. -c $? -o $@

libprecond.so : $(objects)
	ngscxx -shared $(objects) -lngsolve -lngcomp -o $@

clean:
	rm *.o libprecond.so

