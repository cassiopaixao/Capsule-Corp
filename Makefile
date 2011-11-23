
LIBS = -lm

all: capsule.o paralela.o mymath.o
	gcc -pg $? $(LIBS) -o ep

%.o: %.c
	gcc -pg -c $^ -o $@

clean:
	rm *.o; rm ep
