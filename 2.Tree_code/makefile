## Commands below for pink - must have loaded modules with 
#              module load lampi intel

## Linux cluster (ie, pink) options - select only one FFLAGS definition
#----------------------------------------
  FC = gfortran   ## compiler
  #MPILIB =  -L${MPI_ROOT}/lib64 -lmpi -I${MPI_ROOT}/include/
  #MPILIB =  -I${MPI_ROOT}/include/
  FFLAGS = -O2 -ffixed-line-length-none 
#  FFLAGS = -g -extend_source 
#   FFLAGS = -g -extend_source -warn all -check bounds -check overflow 

## COMPAQ options - select only one FFLAGS definition
#----------------------------------------
# FC = f90    ## compiler
# MPILIB =  -lmpi -lfmpi -lelan
# FMACH = ${MPI_COMPILE_FLAGS} ${MPI_LD_FLAGS} 
# DEBUG_OPTIONS = -check bounds -check overflow
# FFLAGS =   -g  ${FMACH} ${DEBUG_OPTIONS} #-convert big_endian
# FFLAGS = -O3 ${FMACH} # ${DEBUG_OPTIONS} #-convert little_endian


##SGI Origin 2000 options - select only one FFLAGS definition
#----------------------------------------
# MPILIB =  -lmpi
# FMACH =  -64
# DEBUG_OPTIONS = -C -DEBUG:trap_uninitialized=ON -DEBUG:div_check
# FFLAGS =  -g -r8
# FFLAGS = ${FMACH} -g -r8 -extend_source ${DEBUG_OPTIONS}
# FFLAGS = ${FMACH} -O2  -OPT:Olimit=0 -extend_source
# FFLAGS = ${FMACH} -O3 -r8 -OPT:Olimit=0 -extend_source

OBJS = variables.o main.o define_variables.o \
			 io.o metryc.o shapes.o distribution.o \
			 baseline.o treatments.o

trees: ${OBJS}
	${FC} -o $@  ${FFLAGS} ${OBJS} ${MPILIB}

clean: 
	rm trees *.o *.mod *.dat

# optimized suffix rules
.SUFFIXES: .f .f90

.f.o:
	${FC} ${FFLAGS} -c $<
.f90.o:
	${FC} ${FFLAGS} -c $<
