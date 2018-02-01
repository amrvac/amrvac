#default fortran compiler
fc = ifort

#executable name
exec = test

#default libraries
lib = -llapack

#guess if you have MKL
GREP_MKL = $(shell export | grep 'mkl')
ifneq ($(GREP_MKL),)
	lib = -mkl
endif

#if dgesv is not present non need for LAPACK
GREP = $(shell grep -i 'dgesv' krome_all.f90)
ifeq ($(GREP),)
	lib =
endif

#flags
switchOPT = -O3 -ipo -ip -unroll -xHost -g
switchDBG = -O0 -check all -warn all -fpe0 -u -traceback -warn nounused
switchPRO = $(switchOPT) -pg -traceback
switchOMP = $(switchOPT) -openmp

#default switch
switch = $(switchOPT)

#no warning switch (for DLSODES)
nowarnIfort = -nowarn
nowarnGNU = -w

#objects
objs = opkda2.o
objs += opkda1.o
objs += opkdmain.o
objs += krome_user_commons.o
objs += krome_all.o
objs += test.o

#default target
all: 	$(objs)
	$(fc) $(objs) -o $(exec) $(switch) $(lib)

#ifort full debug target
debug: switch = $(switchDBG)
debug: nowarn = $(nowarnIfort)
debug: all

#ifort profile target
profile: switch = $(switchPRO)
profile: nowarn = $(nowarnIfort)
profile: all

#gfortran target
gfortran: fc = gfortran
gfortran: switch = -ffree-line-length-none
gfortran: nowarn = $(nowarnGNU)
gfortran: all

#gfortran debug target
gfortran_dbg: fc = gfortran
gfortran_dbg: switch = -fbacktrace -g
gfortran_dbg: switch += -ffpe-trap=zero,overflow,invalid
gfortran_dbg: switch += -fbounds-check -ffree-line-length-none
gfortran_dbg: nowarn = $(nowarnGNU)
gfortran_dbg: all

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exec)

#rule for f90
%.o:%.f90
	$(fc) $(switch) -c $^ -o $@

#rule for f
%.o:%.f
	$(fc) $(switch) $(nowarn) -c $^ -o $@

