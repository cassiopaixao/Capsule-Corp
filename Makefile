LIBS = -lgomp -lm
ARCH = native
OTIMIZACAO = -O3
DEBUG = -DDEBUG
NUM_CORES = `cat /proc/cpuinfo | grep processor | wc -l`
FLAGS = -mtune=$(ARCH) $(OTIMIZACAO) -fopenmp -DNUM_THREADS=$(NUM_CORES) $(DEBUG) # -pg

all: capsule.o paralela.o mymath.o
	gcc $? $(FLAGS) $(LIBS) -o ep

%.o: %.c
	gcc $(FLAGS) -c $^ -o $@

clean:
	rm *.o; rm ep


