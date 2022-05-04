CC=g++
CFLAGS=-I. -I../src -g -pthread -fopenmp
DEPS = bitops.hpp sparse_format.h
OBJ = alto_test.o sparse_format.o

%.o: %.c $(DEPS)
	$(CC) -g -c -o $@ $< $(CFLAGS)

alto: $(OBJ)
	$(CC) -g -o $@ $^ $(CFLAGS)

clean:
	rm -rf *.o alto
