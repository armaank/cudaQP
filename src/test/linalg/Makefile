CC = gcc
CFLAGS  = -Wall

test.out: test_linalg.c test_linalg.h
	$(CC) $(CFLAGS) test_linalg.c ldl.o linalg.o csc.o -o tests.out -lm

csc.o: ../../csc.c ../../../include/csc.h
	$(CC) $(CFLAGS) -c ../../csc.c ../../linalg.c

linalg.o: ../../linalg.c ../../../include/linalg.h
	$(CC) $(CFLAGS) -c ../../linalg.c -lm

ldl.o: ../../../ldl/ldl.c ../../../ldl/ldl.h
	$(CC) $(CFLAGS) -c ../../../ldl/ldl.c

all: ldl.o linalg.o csc.o

clean:
	rm *.o *.out
