LIBS = -lgomp -lm
ARCH = native
NUM_CORES = `cat /proc/cpuinfo | grep processor | wc -l`
FLAGS = -mtune=$(ARCH) -O3 -fopenmp -DNUM_THREADS=$(NUM_CORES) # -pg

all: capsule.o paralela.o mymath.o
	gcc $? $(FLAGS) $(LIBS) -o ep

%.o: %.c
	gcc $(FLAGS) -c $^ -o $@

clean:
	rm *.o; rm ep


