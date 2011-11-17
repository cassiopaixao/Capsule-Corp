
LIBS = -lm

all: capsule.o paralela.o mymath.o
	gcc $? $(LIBS) -o ep

%.o: %.c
	gcc -c $^ -o $@
