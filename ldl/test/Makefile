CC = gcc
CFLAGS  = -Wall

make: test.c ldl.o
	$(CC) $(CFLAGS) test.c ldl_unittest.h ldl.o -o tests.out

ldl.o: ../ldl.c ../ldl.h
	$(CC) $(CFLAGS) -c ../ldl.c ../ldl.h

clean:
	rm *.o *.out

