# SET THE USED FORTRAN COMPILER
#FC = gfortran -std=legacy
FC = ifort 

# DEBUGGING OPTIONS
FDEBUG  = -g -traceback -check all -debug all -ftrapuv

# ADD ANY REQUIRED FLAGS HERE (UNCOMMENT FDEBUG FOR DEBUGGING, REMOVES -O2 cCOMPILER FLAG)
FFLAGS  = -O0 $(FDEBUG) 

# ADD ANY REQUIRED LIBARYS HERE
MKLLIB = -mkl=sequential
LIBS = ${MKLLIB}

SRCDIR=src/
SRC=${SRCDIR}align.F90

# OBJECT FILES
MAIN=$(SRC:.F90=.o)

VERSION=$(shell git describe)

# PROGRAMS
INST=bin
EXE=.x
PROGS=align$(EXE)

#
.PHONY:	all
all:	$(PROGS)
	mv *$(EXE) $(INST)
#
.PHONY:	clean
clean:	
	for f in $(PROGS) ; do rm -f $(INST)/$${f} ; done
	rm $(SRCDIR)*.o #$(SRCDIR)*.mod
#
align.x: $(MAIN)
	$(FC) $(FFLAGS) -o $@ ${MAIN} $(LIBS)
#
%.o: %.F90 
	$(FC) -Dversion=$(VERSION) -c $(FFLAGS) -o $@ $<
