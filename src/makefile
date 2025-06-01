FORTRAN=mpif90 # Define the compiler
OPTIONS = -O3 -fdefault-real-8 -lfftw3

# source files
EXE = laps
MAIN = mhd.f90 
SRC4 = mhdrhs.f90 restart.f90
SRC3 = fftw.f90 mhdoutput.f90 rktmod.f90 dealiasing.f90 mhdrms.f90 AEBmod.f90
SRC2 = mhdinit.f90
SRC1 = parallel.f90 

OBJMAIN = ${MAIN:.f90=.o}
OBJ1  = ${SRC1:.f90=.o}
OBJ2  = ${SRC2:.f90=.o}
OBJ3  = ${SRC3:.f90=.o}
OBJ4  = ${SRC4:.f90=.o}
OBJ   = $(OBJMAIN) $(OBJ4) $(OBJ3) $(OBJ2) $(OBJ1)

# compilation rule
$(EXE): $(OBJ)
	$(FORTRAN) -o $(EXE) $(OBJ) $(OPTIONS)
$(OBJMAIN): $(OBJ4) $(OBJ3) $(OBJ2) $(MAIN)
	$(FORTRAN) -c $(MAIN) $(OPTIONS)
$(OBJ4): $(OBJ1) $(OBJ2) $(OBJ3) $(SRC4)
	$(FORTRAN) -c $(SRC4) $(OPTIONS)
$(OBJ3): $(OBJ1) $(OBJ2) $(SRC3)
	$(FORTRAN) -c $(SRC3) $(OPTIONS)
$(OBJ2): $(OBJ1) $(SRC2)
	$(FORTRAN) -c $(SRC2) $(OPTIONS)
$(OBJ1): $(SRC1) 
	$(FORTRAN) -c $(SRC1) $(OPTIONS)

install: $(EXE)
	install -m 755 $(EXE) $(CONDA_PREFIX)/bin/
	rm $(EXE)

clean:
	rm *.o *.mod


