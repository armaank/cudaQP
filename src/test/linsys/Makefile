CC = gcc
CFLAGS  = -Wall

test.out: test_linsys.c test_linsys.h
	$(CC) $(CFLAGS) test_linsys.c ldl_interface.o kkt.o csc.o linalg.o ldl.o linsys.o -o tests.out -lm -lamd

linsys.o: ../../linsys.c ../../../include/linsys.h
	$(CC) $(CFLAGS) -c ../../linsys.c ../../../ldl/ldl_interface.c

ldl_interface.o: ../../../ldl/ldl_interface.c ../../../ldl/ldl_interface.c
	$(CC) $(CFLAGS) -c ../../../ldl/ldl_interface.c ../../../ldl/ldl_src/ldl.c ../../kkt.c -lamd

kkt.o: ../../kkt.c ../../../include/kkt.h
	$(CC) $(CFLAGS) -c ../../kkt.c ../../csc.c

csc.o: ../../csc.c ../../../include/csc.h
	$(CC) $(CFLAGS) -c ../../csc.c ../../linalg.c

linalg.o: ../../linalg.c ../../../include/linalg.h
	$(CC) $(CFLAGS) -c ../../linalg.c -lm

ldl.o: ../../../ldl/ldl_src/ldl.c ../../../ldl/ldl_src/ldl.h
	$(CC) $(CFLAGS) -c ../../../ldl/ldl_src/ldl.c

all: ldl.o linalg.o csc.o kkt.o ldl_interface.o linsys.o

clean:
	rm *.o *.out
