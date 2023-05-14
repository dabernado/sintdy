SHELL = /bin/sh

FC = gfortran
FFLAGS = -ffixed-form -std=legacy

# Use implicit rules to compile the source .f to object .o

example: example.o
	gfortran -Wl,-rpath,.. -L.. -lodepack -o example example.f

example.o: example.f
	gfortran -c example.f

.PHONY : clean

clean:
	rm -f *.o *.dll *.a *~ core
