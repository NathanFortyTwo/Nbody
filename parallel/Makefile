CC	= mpicc
CFLAGS	= -O0 -g -Wall -fopenmp
LDFLAGS = -lm -lpthread -lX11 -fopenmp -g 
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

%.o: %.c
	$(CC) $(CFLAGS) -c  $< $(VERBOSE) $(DISPLAY) $(DUMP) 

clean:
	rm -f *.o $(TARGET)

$(foreach tgt,$(TARGET),$(eval $(call RULE_TEMPLATE,$(tgt))))