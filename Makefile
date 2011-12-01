LIBS = -lgomp -lm
ARCH = native
OTIMIZACAO = -O3
DEBUG = #-DDEBUG
PROFILING = #-pg
NUM_CORES = `cat /proc/cpuinfo | grep processor | wc -l`
NUM_THREADS = -DNUM_THREADS=$(NUM_CORES)
MIN_TILES_TO_PARALLEL = 100
FLAGS = -mtune=$(ARCH) $(OTIMIZACAO) -fopenmp $(NUM_THREADS) $(DEBUG) $(PROFILING) -DMIN_TILES_TO_PARALLEL=$(MIN_TILES_TO_PARALLEL)

all: capsule.o paralela.o mymath.o
	gcc $? $(FLAGS) $(LIBS) -o ep

%.o: %.c
	gcc $(FLAGS) -c $^ -o $@

clean:
	rm *.o; rm ep


