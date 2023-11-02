CC	= mpicc
CFLAGS	= -O0 -g -Wall -fopenmp
LDFLAGS = -lm -lpthread -lX11 -fopenmp -g 

CUDA_CC = nvcc
CUDA_CFLAGS = -O0 -g
CUDA_LDFLAGS =  -lm -lpthread -lX11 -g
CUDA_FLAGS = -arch=sm_37 -gencode=arch=compute_37,code=sm_37 -Wno-deprecated-gpu-targets 

VERBOSE	=
TARGET	= mpi_nbody_brute_force omp_nbody_brute_force sequential_nbody_brute_force sequential_nbody_barnes_hut
OBJS	= ui.o xstuff.o nbody_tools.o nbody_alloc.o

DISPLAY = -DDISPLAY
DISPLAY =

DUMP = -DDUMP_RESULT
#DUMP =

all: $(TARGET)

define RULE_TEMPLATE

$1: $1.o $(OBJS)
	$(CC) $(VERBOSE) -o $$@ $$< $(OBJS) $(LDFLAGS)

endef

cuda_nbody_brute_force: cuda_nbody_brute_force.o $(OBJS)
	$(CUDA_CC) $(CUDA_FLAGS) $(VERBOSE) -o $@ $< $(OBJS) $(CUDA_LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c  $< $(VERBOSE) $(DISPLAY) $(DUMP) 

%.o: %.cu
	$(CUDA_CC) $(CUDA_CFLAGS) -c  $< $(VERBOSE) $(DISPLAY) $(DUMP)

clean:
	rm -f *.o $(TARGET)

$(foreach tgt,$(TARGET),$(eval $(call RULE_TEMPLATE,$(tgt))))