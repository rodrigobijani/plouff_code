#############################################
# Makefile for Rodrigo's programs.          #
#                                           #
# Targets are:                              #
#   make all   - test_pmoga                 #
#                                           #
#  make clean - rm all .o and .mod files    #
#                all executables            #
#                and .dSYM directories      #
#############################################

#######################
# Define directories: #
#######################

# PROG = programs's directory
# BLD = Build directory, where executable files (.o and .mod) will go
# TEST = Test programs directory

SRC  =  modules/
TEST = tests/
PROG = programs/

# Directories where build files (.o and .mod) and executables will go:
BLD  = build/
EXE  = exe/

###############################
# Define compilers and flags: #
###############################

# MacOSX: gfortran
# If you are not on a 64 bit machine then remove the -m64 flag.
F90 = gfortran
#F90FLAGS = -g -m64 -J$(BLD)

#F90FLAGS = -g -fbounds-check -fbacktrace -Wall -Wextra -Wno-compare-reals -Wno-unused-function -pedantic -fopenmp -m64 -J$(BLD)
F90FLAGS = -g -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -Wall -Wextra -Wno-compare-reals -Wno-unused-function -pedantic -fopenmp -m64 -J$(BLD)

##########################
# List all the .o files: #
##########################

OBJS = $(BLD)geometric_mod.o $(BLD)gravmag_mod.o $(BLD)kernels_mod.o

##########################
# List all the programs: #
##########################

TESTS = test_kernels_mod

#######################
# Define the targets: #
#######################

tests: $(TESTS)

############################
# cleanning folders
############################
clean:
	rm -rf $(BLD)*.o $(BLD)*.mod

clear:
	rm -rf $(BLD)*.dSYM $(TESTS) 

###############################
# List all the test programs: #
###############################				# files and exe. should be in the build directory

test_kernels_mod: $(TEST)test_kernels_mod.f95 $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) -o $(EXE)test_kernels_mod $(TEST)test_kernels_mod.f95

test_geometric_mod: $(TEST)test_geometric_mod.f95 $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) -o $(EXE)test_geometric_mod $(TEST)test_geometric_mod.f95

test_gravmag_mod: $(TEST)test_gravmag_mod.f95 $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) -o $(EXE)test_gravmag_mod $(TEST)test_gravmag_mod.f95

###############################
# List all the main programs: #
###############################	

##########################
# List the dependencies: #
##########################

$(BLD)geometric_mod.o: $(SRC)geometric_mod.f95
	$(F90) $(F90FLAGS) -c $(SRC)geometric_mod.f95 -o $(BLD)geometric_mod.o  

$(BLD)kernels_mod.o: $(SRC)kernels_mod.f95 $(BLD) geometric_mod.o
	$(F90) $(F90FLAGS) -c $(SRC)kernels_mod.f95 -o $(BLD)kernels_mod.o 

$(BLD)gravmag_mod.o: $(SRC)gravmag_mod.f95 $(BLD) geometric_mod.o $(BLD)kernels_mod.o
	$(F90) $(F90FLAGS) -c $(SRC)gravmag_mod.f95 -o $(BLD)gravmag_mod.o
