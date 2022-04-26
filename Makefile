# this is the modified version of original makefile

# Compiler options inclusion
# Compile and Compiler flags
f90invoke = gfortran 
F90       = -Og -fdefault-real-8 -fdefault-integer-8 -g #-fbacktrace -fcheck=all -fdump-core


# General rules
.SUFFIXES: 
.SUFFIXES: .exe .o .F90 .f90 .mod .f

.F90.o:
	$(f90invoke)  -c $< -o $*.o $(F90)

#.f90.o:
#	$(f90invoke) -c $< -o $*.o $(F90) 

#.f.o:
#	$(f90invoke) -c $< -o $*.o $(F90)

# Object code list
TREE_OBJS =  num_pars.o defined_types.o kind_numbers.o parameter_modules.o   memory_modules.o tree_routines.o modified_merger_tree.o deltcrit.o memory.o sigmacdm_spline.o interp.o locate.o hyperbolic.o split_PCH.o ran3.o spline.o make_tree.o indexxx.o transfer_function.o unresolved_mass.o parameters.o


all:	trees.exe

clean:
	\rm -f ./*.o *.o *.mod trees.exe

#Dependencies
memory.o: memory_modules.o defined_types.o
memory_modules.o: defined_types.o
make_tree.o: defined_types.o
sigmacdm_spline.o: num_pars.o parameters.o
deltcrit.o: num_pars.o
defined_types.o: kind_numbers.o
tree_routines.o: defined_types.o
trees.o: defined_types.o parameter_modules.o   memory_modules.o tree_routines.o modified_merger_tree.o parameters.o

#Rule for making the executable 
trees.exe: $(TREE_OBJS) trees.o 
	$(f90invoke)  $(TREE_OBJS) trees.o -o trees.exe $(F90)

