
all: capsule.o paralela.o mymath.o
	gcc $? -o ep

%.o: %.c
	gcc -c $^ -o $@
