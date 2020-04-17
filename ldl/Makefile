CC = gcc
CFLAGS  = -Wall


make: test.c ldl.o 
	$(CC) $(CFLAGS) test.c ../include/ldl_unittest.h ldl.o -o tests.out

ldl.o: ldl.c ../include/ldl.h
	$(CC) $(CFLAGS) -c ldl.c ../include/ldl.h

clean:
	rm *.o *.out

