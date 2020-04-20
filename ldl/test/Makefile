CC = gcc
CFLAGS  = -Wall

make: test_ldl.c ldl.o
	$(CC) $(CFLAGS) test_ldl.c ../../include/unittest.h ldl.o -o tests.out

ldl.o: ../ldl.c ../ldl.h
	$(CC) $(CFLAGS) -c ../ldl.c 

clean:
	rm *.o *.out

