target = poisson3d

.SUFFIXES:
.SUFFIXES: .f90 .o .c

CC = gcc
FC = gfortran
FFLAGS = -O3 -ffast-math -funroll-loops -fopenmp

OBJS = poisson3d.o precision.o m_init.o m_linspace.o  m_alloc.o  m_input_output.o poisson_methods.o m_gldat.o write_vtk.o 

LIBS =

.PHONY: all
all: $(target)

.PHONY: new
new: clean $(target)

.PHONY: clean realclean
clean:
	@/bin/rm -f $(OBJS) *.mod

realclean: clean
	@/bin/rm -f $(target)

# linking: the target depends on the objects
$(target): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(target)

.f90.o:
	$(FC) -c $(FFLAGS) $<

.c.o:
	$(CC) -c $<

# dependencies:
poisson3d.o: precision.o poisson_methods.o write_vtk.o m_alloc.o m_gldat.o m_init.o m_input_output.o
poisson_methods.o: precision.o m_gldat.o m_alloc.o
m_init.o: precision.o m_gldat.o m_linspace.o
m_alloc.o: precision.o
m_linspace.o: precision.o
m_gldat.o: precision.o
m_input_output.o: precision.o m_gldat.o



